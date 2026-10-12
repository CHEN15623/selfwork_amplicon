#!/usr/bin/env Rscript
# ------------------------------------------------------------------
# sintax_rank_summary.R
#
# Free replacement for `usearch -sintax_summary` (that subcommand is
# available only in the licensed USEARCH build).
#
# Input : OTU table (OTU x sample, counts)  +  sintax taxonomy
# Output: result/tax/sum_<rank>.txt for each requested rank, in the
#         exact format expected by amplicon::tax_stackplot /
#         tax_circlize, i.e.
#             #OTUID  <sample1>  <sample2> ...   (relative abundance %)
#
# The taxonomy column of a sintax table can be either:
#   * the *tabbed* output of `vsearch --sintax --tabbedout` (col 4), or
#   * the 8-field `taxonomy.txt` style (Kingdom..Species)
# Both are auto-detected.
# ------------------------------------------------------------------

suppressWarnings(suppressMessages(library(optparse)))

option_list <- list(
    make_option(c("-i", "--input"), type="character", default="result/otutab.txt",
                help="OTU/ASV table (counts, OTU x sample) [default %default]"),
    make_option(c("-t", "--taxonomy"), type="character", default="result/otus.sintax",
                help="taxonomy file (sintax tabbed or 8-field taxonomy.txt) [default %default]"),
    make_option(c("-o", "--output"), type="character", default="result/tax/",
                help="output directory [default %default]"),
    make_option(c("-r", "--ranks"), type="character", default="p,c,o,f,g",
                help="ranks to summarise, comma separated [default %default]")
)
opts <- parse_args(OptionParser(option_list=option_list))
dir.create(opts$output, showWarnings=FALSE, recursive=TRUE)

ranks <- strsplit(opts$ranks, ",")[[1]]

## ------------------------------------------------------------------
## 1. read the OTU table (counts) and make it relative (percent)
## ------------------------------------------------------------------
otutab <- read.table(opts$input, header=TRUE, row.names=1, sep="\t",
                     comment.char="", quote="", check.names=FALSE,
                     stringsAsFactors=FALSE)
otutab[is.na(otutab)] <- 0
rel <- sweep(as.matrix(otutab), 2, colSums(as.matrix(otutab)), "/") * 100
rel[is.na(rel)] <- 0
samples <- colnames(rel)
cat("OTUs:", nrow(rel), " Samples:", length(samples), "\n")

## ------------------------------------------------------------------
## 2. read taxonomy and build a OTU -> named vector map
## ------------------------------------------------------------------
tax <- read.table(opts$taxonomy, header=FALSE, row.names=1, sep="\t",
                  comment.char="", quote="", check.names=FALSE,
                  stringsAsFactors=FALSE, fill=TRUE)

# pick the column that actually holds a taxonomy string (contains ':' or ';')
tax_col <- NULL
for (j in seq_len(ncol(tax))) {
    v <- tax[[j]]
    if (any(grepl(":", v), na.rm=TRUE) || any(grepl("__", v), na.rm=TRUE)) { tax_col <- j; break }
}
if (is.null(tax_col)) stop("Cannot locate a taxonomy column in ", opts$taxonomy)
cat("Using taxonomy column", tax_col, "of", opts$taxonomy, "\n")

raw <- as.character(tax[[tax_col]])
names(raw) <- rownames(tax)

## parse a single taxonomy string -> named vector keyed by rank letter
parse_lineage <- function(s) {
    out <- c(k=NA, p=NA, c=NA, o=NA, f=NA, g=NA, s=NA)
    if (is.na(s) || s == "") return(out)
    # normalise separators
    s <- gsub(";", ",", s)
    parts <- strsplit(s, ",", fixed=TRUE)[[1]]
    for (p in parts) {
        p <- trimws(p)
        if (p == "") next
        # sintax style : "p:Proteobacteria"  /  "d:Bacteria"
        m <- regmatches(p, regexec("^([a-zA-Z])[:|](.+)$", p))[[1]]
        if (length(m) == 3) {
            key <- tolower(m[2]); val <- m[3]
            if (key == "d") key <- "k"
            if (key %in% names(out)) out[key] <- val
            next
        }
        # 8-field style : "p__Proteobacteria"
        m <- regmatches(p, regexec("^([a-zA-Z])__(.*)$", p))[[1]]
        if (length(m) == 3) {
            key <- tolower(m[2]); val <- m[3]
            if (key == "d") key <- "k"
            if (key %in% names(out)) out[key] <- val
        }
    }
    out
}

lineages <- t(vapply(raw, parse_lineage, character(7)))
colnames(lineages) <- c("k","p","c","o","f","g","s")

## ------------------------------------------------------------------
## 3. collapse to each rank and write
## ------------------------------------------------------------------
for (rk in ranks) {
    if (!rk %in% colnames(lineages)) { warning("unknown rank: ", rk); next }
    lab <- lineages[, rk]
    lab[is.na(lab) | lab == ""] <- "Unassigned"
    lab <- gsub('"', "", lab)
    lab <- gsub("^[a-z]__", "", lab)
    # drop sintax confidence values, e.g. "Bacillus(0.98)" -> "Bacillus"
    lab <- trimws(gsub("\\([0-9.]+\\)", "", lab))
    # RDP marks unknowns like "unclassified_Bacteria"
    lab[grepl("^unclassified", lab, ignore.case=TRUE) |
        grepl("^incertae", lab, ignore.case=TRUE)] <- "Unassigned"

    agg <- rowsum(rel, group=lab, reorder=TRUE)
    # drop all-zero taxa
    agg <- agg[rowSums(agg) > 0, , drop=FALSE]
    # order by overall abundance (descending) - keeps plots readable
    agg <- agg[order(-rowSums(agg)), , drop=FALSE]

    outfile <- file.path(opts$output, paste0("sum_", rk, ".txt"))
    out <- data.frame("#OTUID"=rownames(agg), agg, check.names=FALSE)
    write.table(out, outfile, sep="\t", quote=FALSE, row.names=FALSE)
    cat(sprintf("  %s: %d taxa -> %s\n", rk, nrow(agg), outfile))
}

cat("Done.\n")
