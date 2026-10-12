#!/usr/bin/env Rscript
# ------------------------------------------------------------------
# tax_maptree_treemap.R
#
# Lightweight treemap of the top-N most abundant ASVs, coloured by
# phylum. Used as a fallback for amplicon::format2maptree(), which
# requires the heavy Bioconductor stack (phyloseq + ggtree + ggraph).
#
# Input : OTU table (counts, OTU x sample) and 8-field taxonomy.txt
# Output: a single PDF treemap
# ------------------------------------------------------------------

suppressWarnings(suppressMessages({
    library(optparse)
    library(treemap)
}))

option_list <- list(
    make_option(c("-i", "--input"),    type="character", default="result/otutab.txt"),
    make_option(c("-t", "--taxonomy"), type="character", default="result/taxonomy.txt"),
    make_option(c("-o", "--output"),   type="character", default="result/tax/tax_maptree.pdf"),
    make_option(c("-N", "--topN"),     type="numeric",   default=100),
    make_option(c("-w", "--width"),    type="numeric",   default=183),
    make_option(c("-e", "--height"),   type="numeric",   default=118)
)
opts <- parse_args(OptionParser(option_list=option_list))
dir.create(dirname(opts$output), showWarnings=FALSE, recursive=TRUE)

otutab <- read.table(opts$input, header=TRUE, row.names=1, sep="\t",
                     comment.char="", quote="", check.names=FALSE)
tax <- read.table(opts$taxonomy, header=TRUE, row.names=1, sep="\t",
                  comment.char="", quote="", check.names=FALSE)

mean_abund <- rowMeans(otutab)
ord <- order(mean_abund, decreasing=TRUE)
ids <- rownames(otutab)[ord]
ids <- head(ids, opts$topN)

tax <- tax[ids, , drop=FALSE]
df <- data.frame(
    ID     = ids,
    Phylum = ifelse(is.na(tax[["Phylum"]]) | tax[["Phylum"]]=="", "Unassigned", tax[["Phylum"]]),
    Genus  = ifelse(is.na(tax[["Genus"]])  | tax[["Genus"]] =="", "Unassigned", tax[["Genus"]]),
    Abundance = as.numeric(mean_abund[ids]),
    stringsAsFactors = FALSE
)
df <- df[df$Abundance > 0, , drop=FALSE]
df$Label <- paste0(df$Genus, " (", df$ID, ")")

pal <- c("#CBD588","#599861","#DA5724","#508578","#CD9BCD","#AD6F3B",
         "#673770","#D14285","#652926","#C84248","#8569D5","#5E738F",
         "#D1A33D","#8A7C64","#B0B0B0","#7EBFE5","#FFB444","#98D66D",
         "#F76FC3","#45B1F9")

pdf(opts$output, width=opts$width/25.4, height=opts$height/25.4)
treemap(df,
        index=c("Phylum","Label"),
        vSize="Abundance",
        type="index",
        palette=pal,
        border.col="#4a4a4a",
        fontsize.labels=c(10, 5),
        fontcolor.labels=c("white","grey10"),
        bg.labels=0,
        align.labels=list(c("centre","centre"), c("centre","bottom")),
        overlap.labels=0.5,
        inflate.labels=FALSE,
        title="Top abundant ASVs (area = mean relative abundance, colour = phylum)")
dev.off()
cat("Treemap written:", opts$output, "\n")
