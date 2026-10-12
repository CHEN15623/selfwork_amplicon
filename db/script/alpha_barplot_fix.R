#!/usr/bin/env Rscript
# ------------------------------------------------------------------
# alpha_barplot_fix.R
#
# `amplicon::alpha_barplot()` (v1.19.0) still calls the dplyr verb
# `summarise_()`, which was made defunct in dplyr >= 1.1 and therefore
# always errors out. This script re-implements the exact same interface
# and produces an equivalent bar plot (mean +/- sd with the same
# significance letters as the boxplot), using only current APIs.
#
# It is invoked with the same command-line arguments as the original
# script, so it can be dropped into the pipeline unchanged.
# ------------------------------------------------------------------

suppressWarnings(suppressMessages({
    library(optparse)
    library(ggplot2)
    library(dplyr)
    library(rlang)
    library(agricolae)
}))

option_list <- list(
    make_option(c("-i", "--input"),       type="character", default="result/alpha/vegan.txt"),
    make_option(c("-a", "--alpha_index"), type="character", default="richness"),
    make_option(c("-d", "--design"),      type="character", default="result/metadata.txt"),
    make_option(c("-n", "--group"),       type="character", default="Group"),
    make_option(c("-o", "--output"),      type="character", default="result/alpha/"),
    make_option(c("-w", "--width"),       type="numeric",   default=89),
    make_option(c("-e", "--height"),      type="numeric",   default=59)
)
opts <- parse_args(OptionParser(option_list=option_list))
dir.create(opts$output, showWarnings=FALSE, recursive=TRUE)

alpha_div <- read.table(opts$input, header=TRUE, row.names=1, sep="\t",
                        comment.char="", check.names=FALSE)
metadata <- read.table(opts$design, header=TRUE, row.names=1, sep="\t",
                       comment.char="", stringsAsFactors=FALSE, check.names=FALSE)

idx <- rownames(metadata) %in% rownames(alpha_div)
metadata <- metadata[idx, , drop=FALSE]
alpha_div <- alpha_div[rownames(metadata), , drop=FALSE]

df <- data.frame(value=alpha_div[[opts$alpha_index]],
                 group=factor(metadata[[opts$group]]))
colnames(df) <- c(opts$alpha_index, "group")

## significance letters (LSD for 2 groups, Tukey HSD for >2)
model <- aov(df[[opts$alpha_index]] ~ group, data=df)
if (length(unique(df$group)) == 2) {
    out <- agricolae::LSD.test(model, "group", p.adj="none")
    stat <- out$groups
} else {
    th <- TukeyHSD(model, ordered=TRUE, conf.level=0.95)
    lv <- th$group[, 4]
    stat <- data.frame(groups=multcompView::multcompLetters(lv)[["Letters"]])
}
df$stat <- stat[as.character(df$group), "groups"]

## mean +/- sd per group
summ <- df %>%
    group_by(group) %>%
    summarise(mean=mean(.data[[opts$alpha_index]], na.rm=TRUE),
              sd=sd(.data[[opts$alpha_index]], na.rm=TRUE),
              .groups="drop")

p <- ggplot(summ, aes(x=group, y=mean, fill=group)) +
    geom_col(width=0.65, colour="grey20", linewidth=0.3) +
    geom_errorbar(aes(ymin=mean - sd, ymax=mean + sd), width=0.2, linewidth=0.3) +
    geom_text(aes(y=mean + sd, label=df$stat[match(group, df$group)]),
              vjust=-0.5, size=2.5) +
    labs(x="Groups", y=paste(opts$alpha_index, "index"), fill=opts$group) +
    theme_classic(base_size=7) +
    theme(axis.text=element_text(colour="black"),
          legend.position="none")

ggsave(file.path(opts$output, paste0("barplot_", opts$alpha_index, ".pdf")),
       p, width=opts$width, height=opts$height, units="mm")
cat("Barplot written:",
    file.path(opts$output, paste0("barplot_", opts$alpha_index, ".pdf")), "\n")
