#!/usr/bin/env Rscript
# ------------------------------------------------------------------
# function_faprotax_plot.R
#
# Visualization of FAPROTAX functional predictions for 16S data.
#
# Inputs
#   -i/--input     FAPROTAX collapsed table (function x sample, counts)
#   -d/--design    metadata (first column = sample id, must contain Group)
#   -n/--group     grouping column in metadata              [Group]
#   -c/--compare   comparison pair "A-B" (A = case, B = control)  [auto]
#   -o/--output    output directory                         [result/function]
#
# Outputs (PDF + the numeric tables behind each figure)
#   1. function_stackplot.group.pdf   mean relative abundance per group
#   2. function_stackplot.sample.pdf  relative abundance per sample
#   3. function_heatmap.pdf           top-N functions, z-score by row
#   4. function_boxplot.pdf           per-sample abundance, key functions
#   5. function_diff_bar.pdf          case-vs-control mean difference
#   6. function_volcano.pdf           log2FC vs -log10(FDR)
#   7. function_diff.txt              full differential statistics
#   8. function_pca.pdf               PCA of CLR-transformed profiles
# ------------------------------------------------------------------

suppressWarnings(suppressMessages({
    library(optparse)
    library(ggplot2)
    library(reshape2)
    library(vegan)
}))

option_list <- list(
    make_option(c("-i", "--input"),   type="character", default="result/function/faprotax.txt",
                help="FAPROTAX collapsed table [default %default]"),
    make_option(c("-d", "--design"),  type="character", default="result/metadata.txt",
                help="Metadata file [default %default]"),
    make_option(c("-n", "--group"),   type="character", default="Group",
                help="Group column in metadata [default %default]"),
    make_option(c("-c", "--compare"), type="character", default="",
                help="Comparison pair, e.g. 'M-C' (case-control) [default: first two groups]"),
    make_option(c("-o", "--output"),  type="character", default="result/function/",
                help="Output directory [default %default]"),
    make_option(c("-l", "--legend"),  type="numeric",   default=12,
                help="Number of functions shown separately in stackplot [default %default]"),
    make_option(c("-N", "--topN"),    type="numeric",   default=30,
                help="Number of top functions shown in the heatmap [default %default]"),
    make_option(c("-w", "--width"),   type="numeric",   default=120,
                help="Figure width in mm [default %default]"),
    make_option(c("-e", "--height"),  type="numeric",   default=90,
                help="Figure height in mm [default %default]")
)
opts <- parse_args(OptionParser(option_list=option_list))
outdir <- sub("/$", "", opts$output)
dir.create(outdir, showWarnings=FALSE, recursive=TRUE)
figdir <- file.path(outdir, "figure")
dir.create(figdir, showWarnings=FALSE, recursive=TRUE)

## ------------------------------------------------------------------
## 1. Read data
## ------------------------------------------------------------------
feat <- read.table(opts$input, header=TRUE, row.names=1, sep="\t",
                   comment.char="", quote="", check.names=FALSE,
                   stringsAsFactors=FALSE)
## FAPROTAX writes a "taxonomy" column when --column_names_are_in is used;
## if the table is transposed (samples x functions) fix it automatically.
if (opts$group %in% rownames(feat) && !(opts$group %in% colnames(feat))) {
    feat <- as.data.frame(t(feat), check.names=FALSE)
}
num_cols <- vapply(feat, is.numeric, logical(1))
if (!all(num_cols)) {
    cat("Dropping non-numeric column(s):", paste(names(feat)[!num_cols], collapse=", "), "\n")
    feat <- feat[, num_cols, drop=FALSE]
}

feat <- as.matrix(feat)
feat[is.na(feat)] <- 0
## drop functions with no reads at all
feat <- feat[rowSums(feat) > 0, , drop=FALSE]
rel <- sweep(feat, 2, colSums(feat), "/") * 100     # relative abundance (%)

design <- read.table(opts$design, header=TRUE, sep="\t", comment.char="",
                     stringsAsFactors=FALSE, check.names=FALSE)
rownames(design) <- as.character(design[[1]])
grp_col <- opts$group
if (!grp_col %in% colnames(design)) stop("Column '", grp_col, "' not found in ", opts$design)

common <- intersect(colnames(rel), rownames(design))
if (length(common) < 2) stop("No matching samples between feature table and metadata")
rel <- rel[, common, drop=FALSE]
design <- design[common, , drop=FALSE]
groups <- factor(as.character(design[[grp_col]]))
cat("Samples:", length(common), " Functions:", nrow(rel), "\n")
cat("Groups :", paste(levels(groups), table(groups), sep="=", collapse="  "), "\n")

## comparison pair
if (opts$compare != "" && grepl("-", opts$compare)) {
    cmp <- strsplit(opts$compare, "-")[[1]]
} else {
    cmp <- levels(groups)[1:2]
}
case <- cmp[1]; ctrl <- cmp[2]
cat("Comparison: ", case, " (case) vs ", ctrl, " (control)\n", sep="")

## group means (relative abundance)
grp_mean <- t(apply(rel, 1, function(v) tapply(v, groups, mean)))
grp_mean[is.na(grp_mean)] <- 0
grp_mean <- as.data.frame(grp_mean, check.names=FALSE)
grp_mean <- data.frame(Function=rownames(grp_mean), grp_mean,
                       check.names=FALSE, stringsAsFactors=FALSE)
write.table(grp_mean, file.path(outdir, "function_group_mean.txt"),
            sep="\t", quote=FALSE, row.names=FALSE)

theme_pub <- theme_classic(base_size=7) +
    theme(axis.text=element_text(colour="black"),
          axis.text.x=element_text(angle=45, hjust=1))

## ------------------------------------------------------------------
## 2. Stackplot (top N functions + Other), by group and by sample
## ------------------------------------------------------------------
plot_stack <- function(mat, fill_levels, fname, w, h) {
    ord <- rownames(mat)
    df <- as.data.frame(mat, check.names=FALSE)
    df$Function <- factor(rownames(df), levels=ord)
    m <- melt(df, id.vars="Function", variable.name="Sample", value.name="Abundance")
    m$Function <- factor(as.character(m$Function), levels=fill_levels)
    p <- ggplot(m, aes(x=Sample, y=Abundance, fill=Function)) +
        geom_bar(stat="identity", position="fill", width=0.8) +
        scale_y_continuous(labels=scales::percent) +
        labs(x=NULL, y="Relative abundance (%)", fill=NULL) +
        theme_pub +
        theme(axis.text.x=element_text(angle=45, hjust=1))
    ggsave(fname, p, width=w, height=h, units="mm")
    invisible(p)
}

## functions sorted by overall mean abundance
ord_fun <- rownames(rel)[order(-rowSums(rel))]
topN <- min(opts$legend, length(ord_fun))
keep <- ord_fun[seq_len(topN - 1)]
mat_g <- as.matrix(grp_mean[, -1, drop=FALSE])
rownames(mat_g) <- grp_mean$Function
other <- colSums(mat_g[setdiff(rownames(mat_g), keep), , drop=FALSE])
mat_g2 <- rbind(mat_g[keep, , drop=FALSE], Other=other)
fill_levels <- c(keep, "Other")

## by group
plot_stack(mat_g2, fill_levels, file.path(figdir, "function_stackplot.group.pdf"),
           opts$width, opts$height)
## by sample (all samples, same colour key)
rel2 <- rbind(rel[keep, , drop=FALSE],
              Other=colSums(rel[setdiff(rownames(rel), keep), , drop=FALSE]))
plot_stack(rel2, fill_levels, file.path(figdir, "function_stackplot.sample.pdf"),
           opts$width * 1.4, opts$height)

## ------------------------------------------------------------------
## 3. Heatmap of top N functions (pheatmap, if available)
## ------------------------------------------------------------------
topNfun <- head(ord_fun, opts$topN)
hm_mat <- rel[topNfun, , drop=FALSE]
## annotate columns by group
ann <- data.frame(Group=as.character(groups), row.names=colnames(hm_mat))
ann_colors <- list(Group=setNames(rep(c("#4E79A7","#F28E2B","#59A14F","#E15759",
                                        "#B07AA1","#76B7B2","#EDC948","#FF9DA7"),
                                      length(levels(groups)))[seq_along(levels(groups))],
                                      levels(groups)))
if (requireNamespace("pheatmap", quietly=TRUE)) {
    pdf(file.path(figdir, "function_heatmap.pdf"),
        width=opts$width/25.4, height=max(3, opts$height/25.4))
    pheatmap::pheatmap(log10(hm_mat + 1),
                       scale="row",
                       cluster_rows=TRUE, cluster_cols=TRUE,
                       annotation_col=ann, annotation_colors=ann_colors,
                       fontsize=6, border_color=NA,
                       color=colorRampPalette(c("#2C7BB6","white","#D7191C"))(100),
                       main="Top functions (log10 abundance, row z-score)")
    dev.off()
    cat("Heatmap written\n")
} else {
    cat("pheatmap not available - heatmap skipped\n")
}
write.table(data.frame(Function=rownames(hm_mat), hm_mat, check.names=FALSE),
            file.path(outdir, "function_heatmap_data.txt"),
            sep="\t", quote=FALSE, row.names=FALSE)

## ------------------------------------------------------------------
## 4. Boxplot of the most abundant / most differential functions
## ------------------------------------------------------------------
diff_tab <- NULL
lev <- levels(groups)
if (case %in% lev && ctrl %in% lev && case != ctrl) {
    idx_case <- which(groups == case)
    idx_ctrl <- which(groups == ctrl)
    res <- t(apply(rel, 1, function(v) {
        a <- v[idx_case]; b <- v[idx_ctrl]
        p <- tryCatch(t.test(a, b)$p.value, error=function(e) NA_real_)
        c(mean_case=mean(a), mean_ctrl=mean(b),
          log2FC=log2((mean(a) + 1e-6) / (mean(b) + 1e-6)),
          pvalue=p)
    }))
    res <- as.data.frame(res, check.names=FALSE)
    res$Function <- rownames(res)
    res$FDR <- p.adjust(res$pvalue, method="BH")

    # FAPROTAX groups are NOT mutually exclusive: e.g. human_gut,
    # mammal_gut and animal_parasites_or_symbionts describe the same taxa,
    # so they naturally show identical profiles. Flag them to avoid
    # over-interpreting them as independent evidence.
    prof_key <- apply(rel, 1, function(v) paste(round(v, 8), collapse="|"))
    dup <- duplicated(prof_key) | duplicated(prof_key, fromLast=TRUE)
    if (any(dup)) {
        cat("NOTE: redundant (identical-profile) functional groups:\n")
        cat(paste0("   ", names(prof_key)[dup]), sep="\n")
    }
    res <- res[order(res$FDR, -abs(res$log2FC)), ]
    diff_tab <- res
    write.table(res[, c("Function","mean_case","mean_ctrl","log2FC","pvalue","FDR")],
                file.path(outdir, "function_diff.txt"),
                sep="\t", quote=FALSE, row.names=FALSE)

    ## differential barplot: top 15 by FDR then |log2FC|
    top_d <- head(res[order(res$FDR, -abs(res$log2FC)), ], 15)
    top_d <- top_d[order(top_d$log2FC), ]
    top_d$Function <- factor(top_d$Function, levels=top_d$Function)
    top_d$Direction <- ifelse(top_d$log2FC > 0, paste0("Enriched in ", case),
                                              paste0("Enriched in ", ctrl))
    pd <- melt(top_d[, c("Function","mean_case","mean_ctrl")], id.vars="Function")
    pd$Group <- ifelse(pd$variable == "mean_case", case, ctrl)
    p_bar <- ggplot(pd, aes(x=Function, y=value, fill=Group)) +
        geom_bar(stat="identity", position=position_dodge(width=0.8), width=0.75) +
        coord_flip() +
        labs(x=NULL, y="Mean relative abundance (%)", fill="Group") +
        theme_pub + theme(axis.text.x=element_text(angle=0, hjust=0.5))
    ggsave(file.path(figdir, "function_diff_bar.pdf"), p_bar,
           width=opts$width*1.3, height=max(80, 12*nrow(top_d)), units="mm")

    ## volcano
    vt <- res
    vt$sig <- "ns"
    vt$sig[vt$FDR < 0.05 & vt$log2FC > 0] <- paste0("Up in ", case)
    vt$sig[vt$FDR < 0.05 & vt$log2FC < 0] <- paste0("Up in ", ctrl)
    p_vol <- ggplot(vt, aes(x=log2FC, y=-log10(FDR + 1e-300), colour=sig)) +
        geom_point(size=1.2, alpha=0.85) +
        geom_vline(xintercept=c(-0.5, 0.5), linetype=2, colour="grey50") +
        geom_hline(yintercept=-log10(0.05), linetype=2, colour="grey50") +
        scale_colour_manual(values=setNames(c("#B0B0B0", "#D7191C", "#2C7BB6"),
                                            c("ns", paste0("Up in ", case),
                                                   paste0("Up in ", ctrl)))) +
        labs(x=expression(log[2]~"fold change"), y=expression(-log[10]~"FDR"),
             colour=NULL) +
        theme_pub + theme(axis.text.x=element_text(angle=0, hjust=0.5))
    if (requireNamespace("ggrepel", quietly=TRUE)) {
        lab <- head(vt[vt$sig != "ns", ], 10)
        if (nrow(lab) > 0)
            p_vol <- p_vol + ggrepel::geom_text_repel(data=lab, aes(label=Function),
                                                      size=1.8, max.overlaps=20)
    }
    ggsave(file.path(figdir, "function_volcano.pdf"), p_vol,
           width=opts$width, height=opts$height, units="mm")
}

## ------------------------------------------------------------------
## 5. Boxplot of key functions (top 12 by abundance)
## ------------------------------------------------------------------
box_fun <- head(ord_fun, 12)
bx <- as.data.frame(t(rel[box_fun, , drop=FALSE]), check.names=FALSE)
bx$Sample <- rownames(bx)
bx$Group <- as.character(groups)
bxm <- melt(bx, id.vars=c("Sample","Group"), variable.name="Function", value.name="Abundance")
p_box <- ggplot(bxm, aes(x=Group, y=Abundance, fill=Group)) +
    geom_boxplot(outlier.size=0.5, linewidth=0.3) +
    geom_jitter(width=0.15, size=0.6, alpha=0.7) +
    facet_wrap(~Function, scales="free_y", ncol=4) +
    labs(x=NULL, y="Relative abundance (%)") +
    theme_pub + theme(legend.position="none",
                      axis.text.x=element_text(angle=0, hjust=0.5))
ggsave(file.path(figdir, "function_boxplot.pdf"), p_box,
       width=opts$width*2, height=opts$height*1.5, units="mm")

## ------------------------------------------------------------------
## 6. PCA on centred-log-ratio transformed profiles
## ------------------------------------------------------------------
X <- t(rel)                                  # samples x functions
X <- X[, colSums(X) > 0, drop=FALSE]
Xc <- X + 0.5                                # pseudo-count for CLR
clr <- log(Xc) - rowMeans(log(Xc))
pca <- prcomp(clr, center=TRUE, scale.=FALSE)
ve <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
pd <- data.frame(PC1=pca$x[,1], PC2=pca$x[,2], Group=as.character(groups))
p_pca <- ggplot(pd, aes(x=PC1, y=PC2, colour=Group)) +
    stat_ellipse(level=0.68, linewidth=0.3) +
    geom_point(size=1.8) +
    labs(x=paste0("PC1 (", ve[1], "%)"), y=paste0("PC2 (", ve[2], "%)"),
         colour=NULL) +
    theme_pub + theme(axis.text.x=element_text(angle=0, hjust=0.5))
ggsave(file.path(figdir, "function_pca.pdf"), p_pca,
       width=opts$width, height=opts$height, units="mm")

cat("\nFAPROTAX figures written to: ", figdir, "\n", sep="")
cat("Tables written to: ", outdir, "\n", sep="")
if (!is.null(diff_tab)) {
    cat("\nDifferential functions (FDR<0.05):",
        sum(diff_tab$FDR < 0.05), "of", nrow(diff_tab), "\n")
}
cat("Done.\n")
