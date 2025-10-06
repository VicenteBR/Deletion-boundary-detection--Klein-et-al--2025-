#install.packages("BiocManager")
#BiocManager::install("Gviz")
#BiocManager::install("Rsamtools")
#BiocManager::install("GenomicAlignments")
#BiocManager::install("GenomicRanges")
#BiocManager::install("Rsamtools")
#BiocManager::install("ggplot2")
#BiocManager::install("dplyr")
#BiocManager::install("patchwork")
library(Gviz)
library(Rsamtools)
library(GenomicAlignments)
library(GenomicRanges)
library(Rsamtools)
library(ggplot2)
library(dplyr)
library(patchwork)

setwd("/bam_outs/standard_bam/spliced_reads/filtered_bams/") #Set working directory to the one containing the filtered spliced reads (ROI_filtering.sh output)

# Define the genomic region of interest
gr <- GRanges(seqnames = "CP047231.1", ranges = IRanges(start = 300000, end = 355000))

# BAM files to load
HNH <- "filtered_HNH_sorted.bam"
HNH_n <- "filtered_HNHdeltaNTD_sorted.bam"
ifvding <- "filtered_IFvDinG_sorted.bam"
ifvmut <- "filtered_IFvDinGMut_sorted.bam"
wtjnc <- "filtered_WT_sorted.bam"

# Extract reads in the specified region
param <- ScanBamParam(which = gr, tag = "SA")

HNH_alignments <- readGAlignments(HNH, param = param)
HNH_n_alignments <- readGAlignments(HNH_n, param = param)
ifvding_alignments <- readGAlignments(ifvding, param = param)
ifvmut_alignments <- readGAlignments(ifvmut, param = param)
wtjnc_alignments <- readGAlignments(wtjnc, param = param)

# Identify common junctions from split reads
HNH_j <- summarizeJunctions(HNH_alignments)
HNH_n_j <- summarizeJunctions(HNH_n_alignments)
ifvding_j <- summarizeJunctions(ifvding_alignments)
ifvmut_j <- summarizeJunctions(ifvmut_alignments)
wtjnc_j <- summarizeJunctions(wtjnc_alignments)

# Filter junctions based on a threshold (e.g., at least 3 occurrences)
threshold <- 3
HNH_cj <- HNH_j[HNH_j$score >= threshold, ]
HNH_n_cj <- HNH_n_j[HNH_n_j$score >= threshold, ]
ifvding_cj <- ifvding_j[ifvding_j$score >= threshold, ]
ifvmut_cj <- ifvmut_j[ifvmut_j$score >= threshold, ]
wtjnc_cj <- wtjnc_j[wtjnc_j$score >= threshold, ]

# Create a data track for junction frequency of each library
# Convert junctions to a GRanges object with the score as a metadata column
HNH_junction_frequency <- GRanges(
  seqnames = seqnames(HNH_cj),
  ranges = ranges(HNH_cj),
  score = HNH_cj$score
)

HNH_n_junction_frequency <- GRanges(
  seqnames = seqnames(HNH_n_cj),
  ranges = ranges(HNH_n_cj),
  score = HNH_n_cj$score
)

ifvding_junction_frequency <- GRanges(
  seqnames = seqnames(ifvding_cj),
  ranges = ranges(ifvding_cj),
  score = ifvding_cj$score
)

ifvmut_junction_frequency <- GRanges(
  seqnames = seqnames(ifvmut_cj),
  ranges = ranges(ifvmut_cj),
  score = ifvmut_cj$score
)

wt_junction_frequency <- GRanges(
  seqnames = seqnames(wtjnc_cj),
  ranges = ranges(wtjnc_cj),
  score = wtjnc_cj$score
)

HNH_junc_freq = as.data.frame(HNH_junction_frequency)
HNH_n_junc_freq = as.data.frame(HNH_n_junction_frequency)
ifvding_junc_freq = as.data.frame(ifvding_junction_frequency)
ifvmut_junc_freq = as.data.frame(ifvmut_junction_frequency)
wt_junc_freq = as.data.frame(wt_junction_frequency)

# Filter junctions to keep only those within the 300000 to 360000 range
HNH_junc_freq_filter <- HNH_junc_freq[HNH_junc_freq$start >= 300000 & HNH_junc_freq$end <= 360000, ]
HNH_n_junc_freq_filter <- HNH_n_junc_freq[HNH_n_junc_freq$start >= 300000 & HNH_n_junc_freq$end <= 360000, ]
ifvding_junc_freq_filter <- ifvding_junc_freq[ifvding_junc_freq$start >= 300000 & ifvding_junc_freq$end <= 360000, ]
ifvmut_junc_freq_filter <- ifvmut_junc_freq[ifvmut_junc_freq$start >= 300000 & ifvmut_junc_freq$end <= 360000, ]
wt_junc_freq_filter <- wt_junc_freq[wt_junc_freq$start >= 300000 & wt_junc_freq$end <= 360000, ]

# Calculate the total count of junctions
HNH_total_junctions <- sum(HNH_junc_freq_filter$score)
HNH_n_junctions <- sum(HNH_n_junc_freq_filter$score)
ifvding_junctions <- sum(ifvding_junc_freq_filter$score)
ifvmut_junctions <- sum(ifvmut_junc_freq_filter$score)
wt_junctions <- sum(wt_junc_freq_filter$score)

# Calculate the relative percentage of each junction
HNH_junc_freq_filter$percentage <- (HNH_junc_freq_filter$score / HNH_total_junctions) * 100
HNH_n_junc_freq_filter$percentage <- (HNH_n_junc_freq_filter$score / HNH_n_junctions) * 100
ifvding_junc_freq_filter$percentage <- (ifvding_junc_freq_filter$score / ifvding_junctions) * 100
ifvmut_junc_freq_filter$percentage <- (ifvmut_junc_freq_filter$score / ifvmut_junctions) * 100
wt_junc_freq_filter$percentage <- (wt_junc_freq_filter$score / wt_junctions) * 100

# Sort by width (longer bars at the bottom)
HNH_junc_freq_filter <- HNH_junc_freq_filter[order(HNH_junc_freq_filter$width, decreasing = TRUE), ]
HNH_n_junc_freq_filter <- HNH_n_junc_freq_filter[order(HNH_n_junc_freq_filter$width, decreasing = TRUE), ]
ifvding_junc_freq_filter <- ifvding_junc_freq_filter[order(ifvding_junc_freq_filter$width, decreasing = TRUE), ]
ifvmut_junc_freq_filter <- ifvmut_junc_freq_filter[order(ifvmut_junc_freq_filter$width, decreasing = TRUE), ]
wt_junc_freq_filter <- wt_junc_freq_filter[order(wt_junc_freq_filter$width, decreasing = TRUE), ]

# Add a grouping variable to each table
HNH_junc_freq_filter$group <- "HNH"
HNH_n_junc_freq_filter$group <- "HNH NTD"
ifvding_junc_freq_filter$group <- "IFv DinG"
ifvmut_junc_freq_filter$group <- "IFv mutant"
wt_junc_freq_filter$group <- "WT"

# List of data frames and their corresponding group names
data_frames <- list(
  HNH = HNH_junc_freq_filter,
  HNH_NTD = HNH_n_junc_freq_filter,
  IFv_DinG = ifvding_junc_freq_filter,
  IFv_mutant = ifvmut_junc_freq_filter,
  WT = wt_junc_freq_filter
)

# Assign fixed height and stack bars vertically
bar_height <- 1  # Fixed height for all bars

# Initialize a list to store plots
plots <- list()

# Create plots
plots <- lapply(names(data_frames), function(group_name) {
  df <- data_frames[[group_name]]
  df$group <- group_name
  df$ypos <- seq(0, (nrow(df) - 1) * bar_height, by = bar_height)
  
  total_junctions <- sum(df$score)
  
  ggplot(df) +
    geom_rect(aes(xmin = start, xmax = end, ymin = ypos, ymax = ypos + bar_height), fill = "lightblue", color = "black") +
    geom_text(aes(x = (start + end) / 2, y = ypos + bar_height / 2, label = paste0(round(percentage, 1), "%")), size = 4, color = "black") +
    labs(title = paste(group_name, "(Total Junctions:", total_junctions, ")"), x = "Genomic Position", y = "Junctions") +
    theme_minimal() +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), panel.grid.major.x = element_line(color = "black", linewidth = 1), panel.grid.minor.x = element_line(color = "gray80", linetype = "dashed"), panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), axis.line.x = element_line(color = "black", linewidth = 1)) +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 10), limits = c(300000, 360000)) +
    scale_y_continuous(breaks = NULL)
})

# Arrange plots in grid
grid_layout <- plots[[1]] / 
  plots[[2]] / 
  plots[[5]]
# Display the grid
print(grid_layout)

# Arrange plots in grid with the coverage plot at the bottom
common_theme <- theme_minimal() +
  theme(
    axis.text.x = element_text(color = "black"),
    axis.line.x = element_line(color = "black", linewidth = 1),
    panel.grid.major.x = element_line(color = "black", linewidth = 1),
    panel.grid.minor.x = element_line(color = "gray80", linetype = "dashed"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank()
  )

a <- a + common_theme +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10), limits = c(300000, 360000))

plots[[3]] <- plots[[3]] + common_theme +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10), limits = c(300000, 360000))

grid_layout <- a / plots[[3]] + plot_layout(heights = c(1,4))

# Display the updated grid
print(grid_layout)

