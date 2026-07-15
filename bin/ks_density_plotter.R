#!/usr/bin/env Rscript

# Ks density distribution, one curve per genome pair, from the long-format table written by
# ks_summary.py. The dashed vertical line on each curve is that pair's median Ks -- the same
# value that becomes its entry in ks_matrix.tsv and its branch length in the Ks tree.

library(ggplot2)
library(optparse)

option_list <- list(
    make_option(c("-i", "--input"), type = "character", default = "ks_all.tsv",
                help = "Long-format Ks table from ks_summary.py [default %default]"),
    make_option(c("-o", "--output"), type = "character", default = "ks_density.pdf",
                help = "Output PDF [default %default]"),
    make_option(c("--max_ks"), type = "numeric", default = 2.0,
                help = "Saturation cutoff the data was filtered at; caps the x-axis [default %default]"),
    make_option(c("--xmax"), type = "numeric", default = NULL,
                help = "Force the x-axis limit. Default: fit to the data (99.5th percentile), so a
                large -ks_max does not leave most of the panel empty.")
)
# -in / -out / -max_ks, to match the flag style of the other bin/ scripts.
opt <- parse_args(OptionParser(option_list = option_list))

ks <- read.table(opt$input, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
if (nrow(ks) == 0) {
    stop("No Ks values in ", opt$input)
}
ks$pair <- factor(ks$pair, levels = unique(ks$pair))

medians <- aggregate(Ks ~ pair, data = ks, FUN = median)

# Scale the x-axis to where the data actually is. -ks_max is a saturation FILTER (2.0 by
# default), not a statement about where the Ks peaks sit -- for closely related genomes the
# peaks are below 0.2, and plotting out to 2.0 would leave most of the panel empty.
xmax <- if (!is.null(opt$xmax)) {
    opt$xmax
} else {
    min(opt$max_ks, ceiling(quantile(ks$Ks, 0.995) * 20) / 20)  # round up to the next 0.05
}

# Okabe-Ito: colourblind-safe, and distinguishable in greyscale.
palette <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7",
             "#E69F00", "#56B4E9", "#F0E442", "#000000")
palette <- rep_len(palette, nlevels(ks$pair))

p <- ggplot(ks, aes(x = Ks, colour = pair, fill = pair)) +
    geom_density(alpha = 0.15, linewidth = 0.7, na.rm = TRUE) +
    geom_vline(data = medians, aes(xintercept = Ks, colour = pair),
               linetype = "dashed", linewidth = 0.5, show.legend = FALSE) +
    scale_colour_manual(values = palette) +
    scale_fill_manual(values = palette) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    # coord_cartesian, not scale_x_continuous(limits=): a scale limit DROPS the rows outside it
    # before the density is estimated, which would change the curve. This zooms instead, so the
    # density is still computed over every gene pair that passed the -ks_max filter.
    coord_cartesian(xlim = c(0, xmax)) +
    labs(x = expression(K[s]), y = "Density", colour = NULL, fill = NULL) +
    theme_minimal(base_size = 11) +
    theme(
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
        panel.grid.minor = element_blank(),
        axis.text = element_text(colour = "black"),
        legend.position = c(0.98, 0.98),
        legend.justification = c(1, 1),
        legend.background = element_rect(fill = "white", colour = NA)
    )

pdf(opt$output, width = 6, height = 4.5, useDingbats = FALSE)
print(p)
invisible(dev.off())

cat("Saved Ks density plot to", opt$output, "\n")
for (i in seq_len(nrow(medians))) {
    cat(sprintf("  %s: median Ks = %.4f\n", medians$pair[i], medians$Ks[i]))
}
# [END]
