#!/usr/bin/env Rscript

# Load required libraries
library(ape)
library(optparse)

# Builds a UPGMA tree, given a distance matrix
BuildUPGMATree <- function(D) {
    D <- as.dist(D)
    hc <- hclust(D, method = "average")
    phylo_tree <- as.phylo(hc)
    newick <- write.tree(phylo_tree)
    return(newick)
}

# Rescales a genetic distance matrix into years, for UPGMA to turn into a time tree.
#
# The divergence time is T = D / (2 * rate): a distance counts substitutions accumulated down
# BOTH lineages since their ancestor, so the 2 converts "both lineages" to "one lineage".
#
# But UPGMA applies that same halving itself -- it sets a node's height (ancestor to ONE tip)
# to half the pairwise value. So the matrix handed to it must be the full tip-to-tip path,
# 2T = D / rate. Dividing by 2 here as well would halve the same quantity twice and yield
# node heights of D / (4 * rate), i.e. divergence times exactly 2x too young.
compute_divergence_time <- function(D, mutation_rate) {
    time_matrix <- D / mutation_rate
    time_matrix <- round(time_matrix)  # Round to the nearest whole integer
    return(time_matrix)
}

# Load the distance matrix from a file
load_distance_matrix <- function(file) {
    D <- read.table(file, sep="\t", header=TRUE, row.names=1)
    return(as.matrix(D))
}

# Load the name key from a file and return a named vector
load_name_key <- function(file) {
    name_key <- read.table(file, header=FALSE, stringsAsFactors=FALSE, sep="\t")
    return(setNames(name_key$V2, name_key$V1))
}

# Apply name key to a tree
apply_name_key <- function(tree, name_key) {
    tree$tip.label <- sapply(tree$tip.label, function(x) {
        if (x %in% names(name_key)) {
            return(name_key[[x]])
        } else {
            return(x)
        }
    })
    return(tree)
}

# Plot tree and save as PDF using ape
plot_tree_ape <- function(newick_file, output_pdf, name_key=NULL, tree_type="divergence", mutation_rate=NULL) {
    tree <- read.tree(newick_file)
    
    if (length(tree$tip.label) <= 2) {
        cat("Skipping file", newick_file, "- only two inputs\n")
        return()
    }
    
    # Apply name key if provided
    if (!is.null(name_key)) {
        tree <- apply_name_key(tree, name_key)
    }

    # Set plot title and x-axis label
    title <- if (tree_type == "time") {
        paste0("synLTR ", opt$method, " tree - Mutation rate: ", mutation_rate)
    } else {
        paste0("synLTR ", opt$method, " tree")
    }

    x_axis_label <- if (tree_type == "time") "Time (years)" else opt$xlab
    
    # Plot using ape and save as PDF
    pdf(output_pdf)
    plot(tree, main=title)
    axisPhylo()
    title(xlab = x_axis_label)
    dev.off()
    cat("Saved ape plot to", output_pdf, "\n")
}

# Parse command-line arguments
option_list = list(
    make_option(c("--method"), type="character", default="upgma",
                help="method to build tree [upgma|nj]", metavar="character"),
    make_option(c("--mutation_rate"), type="numeric", default=1.3e-8,
                help="substitution rate for divergence time calculation", metavar="numeric"),
    make_option(c("--name_key"), type="character", default=NULL,
                help="optional key file to rename taxa in the plot", metavar="character"),
    make_option(c("--prefix"), type="character", default="k2p_matrix",
                help="process matrices named <prefix>*.tsv [default %default]", metavar="character"),
    make_option(c("--xlab"), type="character", default="Divergence K",
                help="x-axis label for the distance tree [default %default]", metavar="character")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

# Load name key if provided
name_key <- NULL
if (!is.null(opt$name_key)) {
    name_key <- load_name_key(opt$name_key)
}

# Process all <prefix>*.tsv matrices in the directory (k2p_matrix by default; the Ks pass
# re-runs this with --prefix ks_matrix).
files <- list.files(pattern=paste0("^", opt$prefix, ".*\\.tsv$"))
if (length(files) == 0) {
    stop("No matrices matching '", opt$prefix, "*.tsv' in the working directory.")
}

for (file in files) {
    D <- load_distance_matrix(file)
    
    # Build tree using the specified method
    if (opt$method == "upgma") {
        cat("Using UPGMA method for", file, "\n")
        newick_tree <- BuildUPGMATree(D)
    } else if (opt$method == "nj") {
        cat("Using NJ method for", file, "\n")
        nj_tree <- nj(as.dist(D))
        newick_tree <- write.tree(nj_tree)
    } else {
        stop("Invalid method. Choose 'upgma' or 'nj'.")
    }
    
    # Save the raw Newick tree. The patterns are anchored: unanchored '.tsv'/'.nwk' are
    # regexes whose '.' matches any character, which mangles dotted genome IDs.
    output_file <- sub("\\.tsv$", ".nwk", file)
    write(newick_tree, file=output_file)
    cat("Saved Newick tree to", output_file, "\n")

    # Plot raw distance tree and save as ape PDF
    newick_pdf <- sub("\\.nwk$", ".nwk.ape.pdf", output_file)
    plot_tree_ape(output_file, newick_pdf, name_key, tree_type="divergence")

    # Compute divergence times and save the time-based Newick file without decimals
    divergence_time_matrix <- compute_divergence_time(as.dist(D), opt$mutation_rate)
    time_tree <- BuildUPGMATree(divergence_time_matrix)

    # Save divergence time Newick tree with rounded values (no decimals)
    time_output_file <- sub("\\.tsv$", ".time.nwk", file)
    time_tree_rounded <- gsub("\\.[0-9]+", "", time_tree)  # Remove decimals in Newick
    write(time_tree_rounded, file=time_output_file)
    cat("Saved divergence time Newick tree to", time_output_file, "\n")

    # Plot time-based tree and save as ape PDF. '.time' is already carried by
    # time_output_file, so the substitution must not add a second one.
    time_pdf <- sub("\\.nwk$", ".nwk.ape.pdf", time_output_file)
    plot_tree_ape(time_output_file, time_pdf, name_key, tree_type="time", mutation_rate=opt$mutation_rate)
}
# [END]
