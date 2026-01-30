# Load necessary libraries
library(ggplot2)
library(reshape2)
library(viridis)
library(cowplot)
library(dplyr)
library(vegan)
library(tidyr)

# Read command line arguments
args <- commandArgs(trailingOnly = TRUE)
pav_table_path = args[1]
infection_values_path = args[2]
trait_column = args[3]
output_tag = args[4]
subset_num = as.numeric(args[5])
Ytitle = args[6]

# Example usage (commented out)
# result <- perform_pav_analysis(
#   pav_table_path = "path/to/pav_table.txt",
#   infection_values_path = "path/to/infection_values.txt",
#   trait_column = "YourTraitColumn",
#   output_tag = "MyAnalysis",
#   subset_num = 25
# )

# Number of samples to subset
# Improved K-mer Presence-Absence Variation (PAV) Analysis Script

# Load population data
pop <- read.table("/mnt/data/DanaS/pop_code.txt", header = TRUE, sep = "\t")

# Enhanced function to perform PAV analysis

  # Error handling for input files
  if (!file.exists(pav_table_path)) {
    stop(paste("PAV table file not found:", pav_table_path))
  }
  if (!file.exists(infection_values_path)) {
    stop(paste("Infection values file not found:", infection_values_path))
  }

  # Read input files
  pav_df <- read.table(pav_table_path, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
  infection_df <- read.table(infection_values_path, header = TRUE, sep = "\t", check.names = FALSE)
  infection_df <- merge(infection_df, pop, by = "File")
  # Validate trait column
  if (!(trait_column %in% colnames(infection_df))) {
    stop(paste("Trait column", trait_column, "not found in infection values data"))
  }

  # Preprocessing PAV data
  # Remove k-mers that are in all or no samples
  pav_df <- pav_df[, colSums(pav_df) > 0 & colSums(pav_df) < nrow(pav_df)]

  # Transpose PAV table
  pav_df_transposed <- t(pav_df)

  # Ensure trait values are numeric
  infection_df[[trait_column]] <- as.numeric(infection_df[[trait_column]])

  # Match samples between infection values and PAV data
  common_samples <- intersect(rownames(pav_df_transposed), infection_df$File)
  if (length(common_samples) == 0) {
    stop("No common samples found between PAV table and infection values.")
  }

  # Filter data to common samples
  pav_df_filtered <- pav_df_transposed[common_samples, , drop = FALSE]
  pheno_df_filtered <- infection_df[infection_df$File %in% common_samples, ]



  # Custom Linkage Disequilibrium (LD) calculation function
  calculate_r2 <- function(x, y) {
    p_x1 <- mean(x == 1)
    p_y1 <- mean(y == 1)
    p_x1_y1 <- mean(x == 1 & y == 1)
    r2 <- (p_x1_y1 - (p_x1 * p_y1))^2 / (p_x1 * (1 - p_x1) * p_y1 * (1 - p_y1))
    return(ifelse(is.nan(r2), 0, r2))
  }

  # Calculate pairwise LD
  num_columns <- ncol(pav_df_filtered)
  ld_results <- matrix(0, nrow = num_columns, ncol = num_columns)
  colnames(ld_results) <- colnames(pav_df_filtered)
  rownames(ld_results) <- colnames(pav_df_filtered)

  for (i in 1:(num_columns - 1)) {
    for (j in (i + 1):num_columns) {
      r2_value <- calculate_r2(pav_df_filtered[, i], pav_df_filtered[, j])
      ld_results[i, j] <- r2_value
      ld_results[j, i] <- r2_value
    }
  }

  # Order k-mers by average LD
  ld_means <- rowMeans(ld_results, na.rm = TRUE)
  ordered_kmers <- names(sort(ld_means, decreasing = TRUE))
  pav_df_filtered <- pav_df_filtered[, ordered_kmers]

  # Clustering
  #dist_matrix_pav <- vegdist(pav_df_filtered, method = "jaccard")
  #hclust_result_pav <- hclust(dist_matrix_pav, method = "complete")
  #dend_pav <- as.dendrogram(hclust_result_pav)
  #cluster_order_pav <- order.dendrogram(dend_pav)

  # Order data
    # Combine PAV and trait data
  data_combined <- cbind(pav_df_filtered, trait_value = pheno_df_filtered[[trait_column]])
  #data_combined_ordered <- pav_df_filtered[rownames(pav_df_filtered), , drop = FALSE]
  data_combined_ordered <- data_combined[order(data_combined[, "trait_value"], decreasing = FALSE), , drop = FALSE]
#  data_combined_ordered <- data_combined_ordered[order(pav_df_filtered[, "trait_value"], decreasing = FALSE), , drop = FALSE]

  # Create subset if requested
  if (!is.null(subset_num) && subset_num > 0) {
    subset_indices <- c(1:min(subset_num, nrow(data_combined_ordered)), 
                        (nrow(data_combined_ordered) - subset_num + 1):nrow(data_combined_ordered))
    data_combined_subset <- data_combined_ordered[subset_indices, , drop = FALSE]
  } else {
    data_combined_subset <- data_combined_ordered
  }

  # Prepare melted data for plotting
  df_full_melted <- melt(as.matrix(data_combined_ordered[, -ncol(data_combined_ordered)]))
  df_subset_melted <- melt(as.matrix(data_combined_subset[, -ncol(data_combined_subset)]))

  # Visualization functions
  create_pav_plot <- function(melted_data, is_subset = FALSE) {
    melted_data$value <- as.factor(melted_data$value)
    
    plot <- ggplot(melted_data, aes(x = Var1, y = Var2)) +
      geom_tile(aes(fill = factor(value))) +
      scale_fill_manual(values = c("0" = "lightgrey", "1" = "red"), name = "Presence/Absence") +
      ylab('K-mers') + 
      xlab("Accession") +
      theme_bw() +
      theme(
        axis.text.y = element_blank(),
        axis.title.y = element_text(size = 30),
        axis.text.x = element_text(size = 6, angle = 90, hjust = 1, vjust = 0.5),
        axis.title.x = element_text(size = 30),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank(),
        legend.position = "right",
        legend.text = element_text(size = 40),
        legend.title = element_text(size = 20)
      )
    return(plot)
  }

  create_trait_scatter <- function(data, trait_name, is_subset = FALSE) {
    ggplot(data, 
           aes(x = factor(rownames(data), levels = rownames(data)), 
               y = as.numeric(trait_value))) +
      geom_point(shape = 20, size = 5) +
      ylab(Ytitle) +
      ggtitle(paste0(output_tag, " ", trait_name)) + 
      theme(
        axis.title.y = element_text(size = 30), 
        axis.title.x = element_blank(),
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank()
      )
  }

  # Create plots
  full_pav_plot <- create_pav_plot(df_full_melted)
  subset_pav_plot <- create_pav_plot(df_subset_melted, is_subset = TRUE)
  
  full_scatter_plot <- create_trait_scatter(data_combined_ordered, trait_column)
  subset_scatter_plot <- create_trait_scatter(data_combined_subset, trait_column, is_subset = TRUE)

  # Combine plots
  full_combined_plot <- plot_grid(full_scatter_plot, full_pav_plot, 
                                   align = "v", nrow = 2, rel_heights = c(1, 3))
  subset_combined_plot <- plot_grid(subset_scatter_plot, subset_pav_plot, 
                                     align = "v", nrow = 2, rel_heights = c(1, 3))

  # Save plots
  ggsave(
    filename = paste0(output_tag,"_",trait_column,"_PAV_full1.jpeg"), 
    plot = full_combined_plot, 
    height = 10, width = 25, dpi = 300
  )
  ggsave(
    filename = paste0(output_tag,"_",trait_column,"_PAV_subset1.jpeg"), 
    plot = subset_combined_plot, 
    height = 10, width = 25, dpi = 300
  )



