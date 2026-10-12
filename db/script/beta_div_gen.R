#!/usr/bin/env Rscript
# Generate beta diversity distance matrices (replacement for usearch -beta_div)
# Using vegan::vegdist and ape::bionj for tree-based metrics

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("vegan"))
suppressPackageStartupMessages(library("ape"))

option_list <- list(
    make_option(c("-i", "--input"), type="character", default="result/otutab_rare.txt",
                help="Rarefied OTU table [default %default]"),
    make_option(c("-t", "--tree"), type="character", default=NULL,
                help="Phylogenetic tree (optional) [default %default]"),
    make_option(c("-p", "--prefix"), type="character", default="result/beta/",
                help="Output prefix of distance matrices [default %default]")
)
opts <- parse_args(OptionParser(option_list=option_list))

suppressWarnings(dir.create(opts$prefix, showWarnings=FALSE))

# Read OTU table (OTUs x samples)
otutab <- read.table(opts$input, header=TRUE, row.names=1, sep="\t", comment.char="", check.names=FALSE)
# Transpose to samples x OTUs
otutab_t <- as.data.frame(t(otutab))
cat("Samples:", nrow(otutab_t), "OTUs:", ncol(otutab_t), "\n")

# Compute Bray-Curtis distance
bray <- vegdist(otutab_t, method="bray")
write.table(as.matrix(bray), file=file.path(opts$prefix, "bray_curtis.txt"),
            sep="\t", quote=FALSE, col.names=NA)
cat("bray_curtis.txt written\n")

# Compute other distances
# Jaccard
jac <- vegdist(otutab_t, method="jaccard", binary=TRUE)
write.table(as.matrix(jac), file=file.path(opts$prefix, "jaccard.txt"),
            sep="\t", quote=FALSE, col.names=NA)
cat("jaccard.txt written\n")

# Euclidean
euc <- vegdist(otutab_t, method="euclidean")
write.table(as.matrix(euc), file=file.path(opts$prefix, "euclidean.txt"),
            sep="\t", quote=FALSE, col.names=NA)
cat("euclidean.txt written\n")

# Phylogenetic distances (require tree)
if (!is.null(opts$tree) && file.exists(opts$tree)) {
    tree <- read.tree(opts$tree)
    cat("Tree loaded:", length(tree$tip.label), "tips\n")
    if (requireNamespace("picante", quietly=TRUE)) {
        library(picante)
        # NOTE: the picante build available here exposes
        #   unifrac(comm, tree)
        # i.e. it computes the UNWEIGHTED UniFrac only; there is no
        # `weighted=` argument (that is the phyloseq/GG2 API). Passing
        # weighted=TRUE aborts the script, so it is not used.
        # FastTree produces an UNROOTED tree; picante requires a rooted one.
        if (!is.rooted(tree)) {
            tree <- root(tree, outgroup=tree$tip.label[1], resolve.root=TRUE)
            cat("Tree rooted for UniFrac\n")
        }
        uuf <- unifrac(otutab_t, tree)
        write.table(as.matrix(uuf), file=file.path(opts$prefix, "unweighted_unifrac.txt"),
                    sep="\t", quote=FALSE, col.names=NA)
        cat("unweighted_unifrac.txt written\n")
    } else {
        cat("picante not installed, skipping unifrac\n")
    }
} else {
    cat("No tree provided, skipping phylogenetic distances\n")
}
