#!/usr/bin/env Rscript
# Generate rarefaction curve data (replacement for usearch -alpha_div_rare)
# Output format matches the amplicon::alpha_rare_curve expectation:
#   rows = depth points, columns = samples

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("vegan"))

option_list <- list(
    make_option(c("-i", "--input"), type="character", default="result/otutab_rare.txt",
                help="Input rarefied OTU table [default %default]"),
    make_option(c("-o", "--output"), type="character", default="result/alpha/alpha_rare.txt",
                help="Output rarefaction file [default %default]")
)
opts <- parse_args(OptionParser(option_list=option_list))

# Read rarefied OTU table (OTUs as rows, samples as columns)
otutab <- read.table(opts$input, header=TRUE, row.names=1, sep="\t", comment.char="", check.names=FALSE)
# Transpose: samples as rows (vegan format)
otutab_t <- as.data.frame(t(otutab))
samples <- rownames(otutab_t)
cat("Samples:", length(samples), "\n")
cat("OTUs:", ncol(otutab_t), "\n")

# Compute rarefaction curves
depth <- max(rowSums(otutab_t))
step <- max(1, floor(depth / 50))  # finer resolution
out <- rarecurve(otutab_t, step=step, sample=depth, tidy=FALSE)

# Format as matrix: rows = depth, cols = samples
all_x <- sort(unique(unlist(lapply(out, function(x) attr(x, "Subsample")))))
result <- matrix(NA, nrow=length(all_x), ncol=length(samples),
                 dimnames=list(as.character(all_x), samples))

for (i in seq_along(samples)) {
    x_vals <- attr(out[[i]], "Subsample")
    y_vals <- out[[i]]
    idx <- match(x_vals, all_x)
    result[idx, i] <- y_vals
}

# Fill NA: forward fill then backward fill (by column)
for (i in seq_len(ncol(result))) {
    col <- result[, i]
    # forward fill
    last_good <- NA
    for (j in seq_len(nrow(result))) {
        if (!is.na(col[j])) {
            last_good <- col[j]
        } else if (!is.na(last_good)) {
            col[j] <- last_good
        }
    }
    # backward fill
    last_good <- NA
    for (j in rev(seq_len(nrow(result)))) {
        if (!is.na(col[j])) {
            last_good <- col[j]
        } else if (!is.na(last_good)) {
            col[j] <- last_good
        }
    }
    col[is.na(col)] <- 0
    result[, i] <- col
}

# amplicon::alpha_rare_curve reads with:
#   alpha_rare = alpha_rare[, rownames(metadata)] -- columns must match sample IDs
# So we write with samples as columns (no row name column header = "")

write.table(result, file=opts$output, sep="\t", quote=FALSE, col.names=NA)
cat("Rarefaction curve data written to", opts$output, "\n")
cat("Dimensions:", nrow(result), "depth points x", ncol(result), "samples\n")
