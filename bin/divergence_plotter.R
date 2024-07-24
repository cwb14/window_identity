# Chris Benson
# July 11th, 2024.

#!/usr/bin/env Rscript

# Required Libraries.
library(ggplot2)
library(dplyr)
library(tidyr)

# Function to display help message.
show_help <- function() {
  cat("
Description:
  This script reads a TSV file containing chromosome data with associated sequence divergence values.
  It then plots these values across the chromosomes with line breaks at positions without alignment data (denoted by 'NA').
  Optionally, it includes mean lines from .report files if '-include_mean_line yes' is specified.
  Theres a known bug, where '-include_mean_line yes', causes chromosome panels to be presented in non-numeric order. Not sure why. Issue is occurs in ggplot or immediately prior.

Usage:
  Rscript divergence_plotter.R -in <input_file> -include_mean_line <yes|no> -ymax [INT]

Arguments:
  -in <input_file>              A TSV file with columns: chromosome, start_position, end_position, accession_1, accession_2, ..., accession_n.
                                Columns for accession values should be numeric or 'NA' to denote missing data.
  -include_mean_line <yes|no>   Optional argument. If 'yes', includes mean lines from .report files. Default is 'no'.
  -ymax [INT]                   Optional integer argument. If provided, overrides the Y-axis maximum value with the specified integer.

Input Data Example:
  chromosome    start_position  end_position  Y476h2       WWSL
  NIPT2T_chr1   0               1000000       0.007243465  0.007775874
  NIPT2T_chr1   500000          1500000       0.007263594  NA
  NIPT2T_chr1   1000000         2000000       NA           0.009170159
  NIPT2T_chr2   0               1000000       0.005857224  0.009170159
  NIPT2T_chr2   500000          1500000       0.006539817  0.007775875

Output:
  A PDF file named <input_file>_weighted_avg_de_across_chromosomes.pdf containing the plotted data.

Example:
  Rscript divergence_plotter.R -in weighted_average.bed -include_mean_line yes -ymax 0.04

")
}

# Read arguments from the command line.
args <- commandArgs(trailingOnly = TRUE)

# Initialize default values.
input_file <- NULL
include_mean_line <- "no"
ymax <- NULL

# Parse arguments.
for (i in seq(1, length(args), by = 2)) {
  if (args[i] == "-in") {
    input_file <- args[i + 1]
  } else if (args[i] == "-ymax") {
    ymax <- as.numeric(args[i + 1])
  } else if (args[i] == "-include_mean_line") {
    include_mean_line <- args[i + 1]
  }
}

# Check if input file is provided.
if (is.null(input_file) || !(include_mean_line %in% c("yes", "no"))) {
  show_help()
  quit(status = 1)
}

# Read data.
data <- read.csv(input_file, sep = "\t", header = TRUE, na.strings = "NA")

# Get the column names for accession values.
accession_cols <- names(data)[!(names(data) %in% c("chromosome", "start_position", "end_position"))]

# Sort accession columns.
accession_cols <- sort(accession_cols)

# Convert chromosome names to a factor with the correct order.
chrom_levels <- unique(data$chromosome)
chrom_levels_sorted <- chrom_levels[order(as.numeric(gsub("[^0-9]", "", chrom_levels)))]
data$chromosome <- factor(data$chromosome, levels = chrom_levels_sorted)

# Melt data to long format for ggplot.
data_long <- data %>%
  pivot_longer(cols = all_of(accession_cols), names_to = "accession", values_to = "value")

# Define a shuffled color palette using hcl.colors with the 'Dark 3' palette.
num_categories <- length(unique(data_long$accession))
set.seed(37)
color_palette <- sample(hcl.colors(num_categories, "Dark 3"))

# Define custom plot theme.
plot_border_theme <- theme(
  panel.border = element_rect(colour = "grey", fill = NA, linewidth = 1),
  axis.text.x = element_text(),
  axis.ticks.x = element_line(),
  plot.margin = margin(1, 1, 1, 1, "lines")
)

# Create the base plot.
p <- ggplot(data_long, aes(x = start_position / 1e6, y = value, color = accession, group = accession)) +
  geom_line(na.rm = TRUE) + # Remove NA values to break the line.
  geom_point(size = 0, show.legend = TRUE) +
  facet_wrap(~ chromosome, scales = "free_x") +
  labs(title = "Gap Compressed Sequence Identity Extracted from Minimap2 de Scores",
       x = "Position on chromosome (Mb)",
       y = "Sequence Divergence",
       color = "Accession") +
  theme_minimal() +
  plot_border_theme +
  scale_color_manual(values = color_palette, guide = guide_legend(override.aes = list(shape = 15, size = 6))) +
  scale_x_continuous(expand = c(0, 0))

# Apply Y-axis scaling if ymax is provided.
if (!is.null(ymax)) {
  p <- p + scale_y_continuous(limits = c(0, ymax), expand = expansion(mult = c(0.005, 0.01)))
} else {
  p <- p + scale_y_continuous(expand = expansion(mult = c(0.005, 0.01)))
}

# If include_mean_line is "yes", add geom_hline for each accession based on .report files.
if (include_mean_line == "yes") {
  # Read .report files with flexible naming.
  report_files <- list.files(pattern = "\\.report$")
  for (file in report_files) {
    report_data <- read.csv(file, sep = "\t", header = TRUE)
    accession <- sub(".*\\.([^\\.]+)\\.report", "\\1", file)
    
    # Adjust column name if necessary.
    if ("avg_snp_fract" %in% colnames(report_data)) {
      colnames(report_data)[colnames(report_data) == "avg_snp_fract"] <- "avg_weighted_de"
    }
    
    # Filter out WholeGenome rows.
    report_data <- report_data[report_data$chromosome != "WholeGenome",]
    
    # Match the color of the accession line.
    accession_color <- color_palette[which(sort(unique(data_long$accession)) == accession)]
    
    # Add horizontal mean lines for each chromosome.
    for (chrom in unique(report_data$chromosome)) {
      avg_de <- report_data$avg_weighted_de[report_data$chromosome == chrom]
      p <- p + geom_hline(aes(yintercept = avg_de), data = data.frame(chromosome = chrom), color = accession_color, linetype = "dotted")
    }
  }
}

# Save the plot to a PDF file.
output_file <- sub(".txt$", "_weighted_avg_de_across_chromosomes.pdf", input_file)
output_file <- sub(".bed$", "_weighted_avg_de_across_chromosomes.pdf", output_file)
pdf(output_file, width = 16, height = 6)
print(p)
dev.off()
