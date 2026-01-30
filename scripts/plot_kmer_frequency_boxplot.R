# ===============================================================================
# K-mer Frequency Boxplot Analysis in R
# X-axis: Resistant/Susceptible accessions (predefined groups)
# Y-axis: K-mer frequency (proportion present)
# ===============================================================================

# library(ggplot2)
# library(dplyr)
# library(readr)
# library(tidyr)
# library(RColorBrewer)

# # ===============================================================================
# # CONFIGURATION - CHOOSE CHROMOSOME
# # ===============================================================================
 setwd("/mnt/data/DanaS/sam_broom_gwas/gwas/haplotyps/")
# Set which chromosome to analyze
CHROMOSOME <- "Chr07"  # Change to "Chr07" for chromosome 7

if (CHROMOSOME == "Chr09") {
  # Chr09 configuration
  KMER_FILE <- "/mnt/data/DanaS/sam_broom_gwas/gwas/haplotyps/pav_table_introg09.txt"
  
  KMER_SEQUENCES <- c(
    "ACTTGACTTTTTAAAAGTCAAAAGTGGACAG",
    "ATCTGTCCACTTTTGACTTTTAAAAAGTCAA", 
    "CGCTGATATGTAGCCCTCGATACTAAACCTC",
    "CTTGACTTTTTAAAAGTCAAAAGTGGACAGA"

  )


  
} else if (CHROMOSOME == "Chr07") {
  # Chr07 configuration
  KMER_FILE <- "/mnt/data/DanaS/sam_broom_gwas/gwas/haplotyps/pav_table_introg07.txt"
  
  KMER_SEQUENCES <- c(
    "TTATCAAAGAGACATGATTTGTGGGTCCAAA",
    "TTTACCAGGTACCATTCTTCTTTACGCAAAA",
    "TCGTTCGGTTCGTTTACAACCCTAATGGTAA",
    "TTAAAAAGTTCATTATCTTTATAGTTTAAAA",
    #"TGATCTTGAATTCGGATTAGGCCCTGATCAA",
    #"TTCCGGCTATTAACCCATATTTCCGGCTAAA",
    "TTTAACATGTTTTACACACTTAACCCCTAAA",
    #"GTAAAATCCATGAAAATTACAAAGTCACACA",
    "GTCGTTCAGAAGAATTCTGATCAATACGAGA"
    #"TCAATAACATCTTCACCTTCTAGAGTTACAA",
    #"ACCTAATCCATAACCAATAACTTAACCATAA"
  )
  
} else {
  stop("CHROMOSOME must be either 'Chr07' or 'Chr09'")
}

# # ===============================================================================
# # PREDEFINED SAMPLE GROUPS
# # ===============================================================================
# # Susceptible and resistant Gadot
# SUSCEPTIBLE_SAMPLES <- c("SAM287", "SAM027",
# "SAM029", "SAM039", "SAM069", "SAM072",
# "SAM076", "SAM089", "SAM104", "SAM112",
# "SAM114", "SAM116", "SAM117", "SAM118",
# "SAM130", "SAM136", "SAM155", "SAM184",
# "SAM189", "SAM200", "SAM226", "SAM265"
# )

# RESISTANT_SAMPLES <- c("SAM002", "SAM005", "SAM012",
# "SAM014", "SAM015", "SAM021", "SAM024",
# "SAM025", "SAM032", "SAM034", "SAM038",
# "SAM040", "SAM042", "SAM048", "SAM049",
# "SAM050", "SAM055", "SAM058", "SAM059",
# "SAM063", "SAM065", "SAM066", "SAM067",
# "SAM079", "SAM081", "SAM091", "SAM092",
# "SAM094", "SAM098", "SAM100", "SAM101",
# "SAM107", "SAM119", "SAM120", "SAM121",
# "SAM122", "SAM123", "SAM124", "SAM125",
# "SAM126", "SAM129", "SAM137", "SAM138",
# "SAM140", "SAM141", "SAM142", "SAM147",
# "SAM152", "SAM153", "SAM158", "SAM161",
# "SAM162", "SAM164", "SAM170", "SAM171",
# "SAM172", "SAM173", "SAM174", "SAM177",
# "SAM178", "SAM180", "SAM188", "SAM191",
# "SAM193", "SAM194", "SAM197", "SAM198",
# "SAM199", "SAM201", "SAM202", "SAM205",
# "SAM206", "SAM207", "SAM209", "SAM210",
# "SAM211", "SAM212", "SAM213", "SAM216",
# "SAM227", "SAM229", "SAM230", "SAM231",
# "SAM232", "SAM234", "SAM238", "SAM241",
# "SAM243", "SAM253", "SAM257", "SAM261",
# "SAM270", "SAM278", "SAM281", "SAM282",
# "SAM283", "SAM284", "SAM287"
# )


 # Resistant and susceptible Yavor
SUSCEPTIBLE_SAMPLES <- c("SAM003", "SAM008", "SAM012",
"SAM020", "SAM023", "SAM024", "SAM027",
"SAM033", "SAM034", "SAM035", "SAM039",
"SAM048", "SAM050", "SAM057", "SAM059",
"SAM062", "SAM063", "SAM065", "SAM072",
"SAM074", "SAM076", "SAM077", "SAM078",
"SAM080", "SAM081", "SAM082", "SAM084",
"SAM087", "SAM089", "SAM091", "SAM099",
"SAM103", "SAM104", "SAM106", "SAM109",
"SAM112", "SAM113", "SAM114", "SAM115",
"SAM118", "SAM122", "SAM136", "SAM139",
"SAM142", "SAM143", "SAM158", "SAM167",
"SAM178", "SAM179", "SAM180", "SAM185",
"SAM189", "SAM192", "SAM200", "SAM212",
"SAM214", "SAM220", "SAM221", "SAM229",
"SAM243", "SAM253", "SAM256", "SAM258",
"SAM283", "SAM284"
)

RESISTANT_SAMPLES <- c("SAM001", "SAM002", "SAM011",
"SAM013", "SAM014", "SAM021", "SAM028",
"SAM030", "SAM032", "SAM038", "SAM049",
"SAM051", "SAM056", "SAM058", "SAM079",
"SAM097", "SAM098", "SAM100", "SAM107",
"SAM111", "SAM127", "SAM129", "SAM133",
"SAM137", "SAM146", "SAM147", "SAM149",
"SAM150", "SAM151", "SAM153", "SAM160",
"SAM162", "SAM165", "SAM166", "SAM168",
"SAM173", "SAM174", "SAM187", "SAM188",
"SAM191", "SAM197", "SAM198", "SAM199",
"SAM203", "SAM204", "SAM205", "SAM206",
"SAM208", "SAM209", "SAM210", "SAM211",
"SAM215", "SAM217", "SAM230", "SAM231",
"SAM238", "SAM240", "SAM241", "SAM259",
"SAM263", "SAM274", "SAM275", "SAM278",
"SAM281", "SAM282", "SAM285", "SAM287"
)

# ===============================================================================
# LOAD DATA
# ===============================================================================
# ===============================================================================
# CORRECTED: K-mer Frequency Boxplot Analysis in R
# Using actual PAV file data and corrected sample groups
# ===============================================================================

library(ggplot2)
library(dplyr)
library(readr)
library(tidyr)

# ===============================================================================
# CONFIGURATION
# ===============================================================================

# Chr09 configuration (using actual PAV file data)
# KMER_FILE <- "/mnt/data/DanaS/sam_broom_gwas/gwas/haplotyps/pav_table_introg09.txt"

# # Actual k-mer sequences from PAV file
# KMER_SEQUENCES <- c(
#   "ACTTGACTTTTTAAAAGTCAAAAGTGGACAG",
#   "ATCTGTCCACTTTTGACTTTTAAAAAGTCAA", 
#   "CGCTGATATGTAGCCCTCGATACTAAACCTC",
#   "CTTGACTTTTTAAAAGTCAAAAGTGGACAGA"
# )

# ===============================================================================
# CORRECTED SAMPLE GROUPS
# ===============================================================================

# Corrected sample groups (22 susceptible + 84 resistant = 106 total)
# # Note: Removed duplicate SAM287 from resistant group
# SUSCEPTIBLE_SAMPLES <- c("SAM027", "SAM029", "SAM039", "SAM069", "SAM072", "SAM076",
#                         "SAM089", "SAM104", "SAM112", "SAM114", "SAM116", "SAM117",
#                         "SAM118", "SAM130", "SAM136", "SAM155", "SAM184", "SAM189",
#                         "SAM200", "SAM226", "SAM265", "SAM287")  # 22 samples

# RESISTANT_SAMPLES <- c("SAM002", "SAM005", "SAM012", "SAM014", "SAM015", "SAM021", "SAM024",
#                       "SAM025", "SAM032", "SAM034", "SAM038", "SAM040", "SAM042", "SAM048", 
#                       "SAM049", "SAM050", "SAM055", "SAM058", "SAM059", "SAM063", "SAM065", 
#                       "SAM066", "SAM067", "SAM079", "SAM081", "SAM091", "SAM092", "SAM094", 
#                       "SAM098", "SAM100", "SAM101", "SAM107", "SAM119", "SAM120", "SAM121",
#                       "SAM122", "SAM123", "SAM124", "SAM125", "SAM126", "SAM129", "SAM137", 
#                       "SAM138", "SAM140", "SAM141", "SAM142", "SAM147", "SAM152", "SAM153", 
#                       "SAM158", "SAM161", "SAM162", "SAM164", "SAM170", "SAM171", "SAM172", 
#                       "SAM173", "SAM174", "SAM177", "SAM178", "SAM180", "SAM188", "SAM191",
#                       "SAM193", "SAM194", "SAM197", "SAM198", "SAM199", "SAM201", "SAM202", 
#                       "SAM205", "SAM206", "SAM207", "SAM209", "SAM210", "SAM211", "SAM212", 
#                       "SAM213", "SAM216", "SAM227", "SAM229", "SAM230", "SAM231", "SAM232", 
#                       "SAM234", "SAM238", "SAM241", "SAM243", "SAM253", "SAM257", "SAM261",
#                       "SAM270", "SAM278", "SAM281", "SAM282", "SAM283", "SAM284")  # 84 samples

# ===============================================================================
# LOAD DATA
# ===============================================================================

cat("===============================================================================\n")
cat("LOADING Chr09 K-MER FREQUENCY DATA (CORRECTED)\n")
cat("===============================================================================\n\n")

# Load K-mer PAV data
cat("Reading k-mer PAV:", KMER_FILE, "\n")
pav_df <- read_tsv(KMER_FILE, show_col_types = FALSE)

# Set k-mer column as row names
kmers <- pav_df[[1]]
pav_matrix <- as.matrix(pav_df[, -1])
rownames(pav_matrix) <- kmers

cat("Found", nrow(pav_matrix), "total k-mers in file\n")
cat("Found", ncol(pav_matrix), "samples in file\n")

# Verify target k-mers
cat("\nVerifying k-mer sequences:\n")
available_kmers <- intersect(KMER_SEQUENCES, rownames(pav_matrix))
for (kmer in KMER_SEQUENCES) {
  if (kmer %in% available_kmers) {
    cat("  ✓ Found k-mer:", substr(kmer, 1, 20), "...\n")
  } else {
    cat("  ✗ K-mer not found:", substr(kmer, 1, 20), "...\n")
  }
}

if (length(available_kmers) == 0) {
  stop("ERROR: No target k-mers found!")
}

# Filter PAV matrix
pav_subset <- pav_matrix[available_kmers, , drop = FALSE]

# Find common samples
common_samples <- intersect(colnames(pav_subset), c(RESISTANT_SAMPLES, SUSCEPTIBLE_SAMPLES))
pav_subset <- pav_subset[, common_samples, drop = FALSE]

cat("\nFiltered to", nrow(pav_subset), "k-mers and", ncol(pav_subset), "samples\n")

# ===============================================================================
# CALCULATE K-MER FREQUENCIES BY GROUP
# ===============================================================================

cat("\n===============================================================================\n")
cat("CALCULATING K-MER FREQUENCIES\n")
cat("===============================================================================\n\n")

# Create sample groups
resistant_common <- intersect(RESISTANT_SAMPLES, common_samples)
susceptible_common <- intersect(SUSCEPTIBLE_SAMPLES, common_samples)

cat("Sample groups (corrected):\n")
cat("Resistant group:", length(resistant_common), "samples\n")
cat("Susceptible group:", length(susceptible_common), "samples\n")
cat("Total:", length(resistant_common) + length(susceptible_common), "samples\n\n")

# Calculate frequencies for each group and k-mer
freq_data <- data.frame()

for (group_name in c("Resistant", "Susceptible")) {
  group_samples <- if (group_name == "Resistant") resistant_common else susceptible_common
  
  if (length(group_samples) == 0) {
    cat("WARNING: No", group_name, "samples found in PAV data!\n")
    next
  }
  
  cat(group_name, "group:", length(group_samples), "samples\n")
  
  for (kmer in available_kmers) {
    kmer_values <- pav_subset[kmer, group_samples]
    frequency <- sum(kmer_values, na.rm = TRUE) / length(kmer_values)
    
    # Create short k-mer label
    kmer_short <- if (nchar(kmer) > 15) paste0(substr(kmer, 1, 15), "...") else kmer
    
    freq_data <- rbind(freq_data, data.frame(
      Group = group_name,
      Kmer = kmer_short,
      Kmer_Full = kmer,
      Frequency = frequency,
      N_Samples = length(group_samples),
      N_Present = sum(kmer_values, na.rm = TRUE),
      N_Absent = length(kmer_values) - sum(kmer_values, na.rm = TRUE),
      stringsAsFactors = FALSE
    ))
    
    cat("  ", substr(kmer, 1, 15), "...:", sum(kmer_values, na.rm = TRUE), "/", 
        length(kmer_values), "=", sprintf("%.3f", frequency), "\n")
  }
}

cat("\nCreated frequency dataframe:", nrow(freq_data), "observations\n\n")

# ===============================================================================
# CREATE BOXPLOT
# ===============================================================================

cat("===============================================================================\n")
cat("CREATING BOXPLOT\n")
cat("===============================================================================\n\n")

# Set factor levels for consistent ordering
freq_data$Group <- factor(freq_data$Group, levels = c("Susceptible", "Resistant"))

# Create main boxplot
p1 <- ggplot(freq_data, aes(x = Group, y = Frequency)) +
  geom_boxplot(aes(fill = Group), linewidth = 1, alpha = 0.7, outlier.shape = NA) +
  geom_jitter(size = 3, alpha = 0.8, width = 0.2, color = "#aaa1a1") +
  scale_fill_manual(values = c("Susceptible" = "#d26e77", "Resistant" = "#6285b3")) +
  labs(
    title = paste0(tolower(CHROMOSOME), " K-mer Frequencies"),
                    
    x = paste("Resistant (n =", length(resistant_common), ") vs Susceptible (n =",length(susceptible_common), ")"),
    y = "K-mer Frequency"
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5),
    axis.title.x = element_text(size = 12),
    axis.text.x = element_text(size = 14, face = "bold"),
    axis.text.y = element_text(size = 14, face = "bold"),
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  ylim(0, max(freq_data$Frequency) + 0.15)

# Add statistical test
resistant_freqs <- freq_data[freq_data$Group == "Resistant", "Frequency"]
susceptible_freqs <- freq_data[freq_data$Group == "Susceptible", "Frequency"]

if (length(resistant_freqs) > 0 & length(susceptible_freqs) > 0) {
  # Mann-Whitney U test
  test_result <- wilcox.test(resistant_freqs, susceptible_freqs)
  p_val <- test_result$p.value
  
  # Add significance annotation
  y_max <- max(freq_data$Frequency)
  y_bar <- y_max + 0.05
  
  sig_text <- ifelse(p_val < 0.0001, "****",
                    ifelse(p_val < 0.001, "***",
                          ifelse(p_val < 0.01, "**",
                                ifelse(p_val < 0.05, "*", "ns"))))
  
  p1 <- p1 + 
    annotate("segment", x = 1, xend = 2, y = y_bar, yend = y_bar, color = "black", linewidth = 1) +
    annotate("segment", x = 1, xend = 1, y = y_bar, yend = y_bar - 0.01, color = "black", linewidth = 1) +
    annotate("segment", x = 2, xend = 2, y = y_bar, yend = y_bar - 0.01, color = "black", linewidth = 1) +
    annotate("text", x = 1.5, y = y_bar + 0.03, label = sig_text, size = 6, fontface = "bold") +
    annotate("text", x = 1.5, y = y_bar + 0.08, label = paste("p =", sprintf("%.4e", p_val)), 
             size = 4, fontface = "italic")
}

# Add mean lines
for (group in c("Susceptible", "Resistant")) {
  group_freqs <- freq_data[freq_data$Group == group, "Frequency"]
  if (length(group_freqs) > 0) {
    mean_freq <- mean(group_freqs)
    
    group_x <- ifelse(group == "Susceptible", 1, 2)
    
    p1 <- p1 + 
      annotate("segment", x = group_x - 0.3, xend = group_x + 0.3, 
               y = mean_freq, yend = mean_freq, color = "red", linewidth = 0.8, alpha = 0.8) +
      annotate("text", x = group_x, y = mean_freq + 0.02, 
               label = paste("μ =", sprintf("%.3f", mean_freq)), 
               size = 3.5, fontface = "bold", color = "#414040")
  }
}

# Save main plot
output_file <- paste0(tolower(CHROMOSOME), "_kmer_frequency_boxplot_corrected.png")
ggsave(output_file, p1, width = 8, height = 8, dpi = 300, bg = "white")
cat("✓ Saved:", output_file, "\n")

# ===============================================================================
# CREATE INDIVIDUAL K-MER COMPARISON
# ===============================================================================

cat("Creating individual k-mer comparison...\n")

# Prepare data for grouped bar plot
freq_wide <- freq_data %>%
  select(Group, Kmer, Frequency) %>%
  pivot_wider(names_from = Group, values_from = Frequency) %>%
  pivot_longer(cols = c("Resistant", "Susceptible"), names_to = "Group", values_to = "Frequency")

freq_wide$Group <- factor(freq_wide$Group, levels = c("Resistant", "Susceptible"))

# Create grouped bar plot
p2 <- ggplot(freq_wide, aes(x = Kmer, y = Frequency, fill = Group)) +
  geom_col(position = position_dodge(width = 0.8), alpha = 0.8, color = "black", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.2f", Frequency)), 
            position = position_dodge(width = 0.8), vjust = -0.3, size = 3, fontface = "bold") +
  scale_fill_manual(values = c("Resistant" = "#2E4057", "Susceptible" = "#E63946")) +
  labs(
    title = paste0(tolower(CHROMOSOME), " Individual K-mer Frequencies"),
    x = "K-mer",
    y = "K-mer Frequency",
    fill = "Group"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "top"
  )

# Save individual comparison plot
output_file2 <- paste0(tolower(CHROMOSOME), "_kmer_individual_comparison_corrected.png")
ggsave(output_file2, p2, width = 14, height = 8, dpi = 300, bg = "white")
cat("✓ Saved:", output_file2, "\n")

# ===============================================================================
# SUMMARY STATISTICS
# ===============================================================================

cat("\n===============================================================================\n")
cat("SUMMARY STATISTICS\n")
cat("===============================================================================\n\n")

if (length(resistant_freqs) > 0 & length(susceptible_freqs) > 0) {
  cat("Overall k-mer frequency comparison:\n")
  cat("  Resistant group mean frequency:", sprintf("%.4f", mean(resistant_freqs)), "\n")
  cat("  Susceptible group mean frequency:", sprintf("%.4f", mean(susceptible_freqs)), "\n")
  cat("  Difference:", sprintf("%.4f", abs(mean(resistant_freqs) - mean(susceptible_freqs))), "\n")
  cat("  Mann-Whitney U test: p =", sprintf("%.4e", p_val), "\n")
  
  # Effect size (Cohen's d)
  if (var(resistant_freqs) > 0 | var(susceptible_freqs) > 0) {
    pooled_sd <- sqrt((var(resistant_freqs) + var(susceptible_freqs)) / 2)
    if (pooled_sd > 0) {
      effect_size <- abs(mean(resistant_freqs) - mean(susceptible_freqs)) / pooled_sd
      cat("  Effect size:", sprintf("%.4f", effect_size), "\n")
    }
  }
}

cat("\nIndividual k-mer frequencies:\n")
cat(paste(rep("-", 60), collapse = ""), "\n")

# Create summary table
summary_table <- freq_data %>%
  select(Kmer_Full, Group, Frequency) %>%
  pivot_wider(names_from = Group, values_from = Frequency, names_prefix = "Freq_") %>%
  mutate(Difference = abs(Freq_Resistant - Freq_Susceptible)) %>%
  arrange(desc(Difference))

for (i in 1:nrow(summary_table)) {
  row <- summary_table[i, ]
  cat(substr(row$Kmer_Full, 1, 20), "...\n")
  cat("  Resistant:", sprintf("%.3f", row$Freq_Resistant), "\n")
  cat("  Susceptible:", sprintf("%.3f", row$Freq_Susceptible), "\n")
  cat("  Difference:", sprintf("%.3f", row$Difference), "\n\n")
}

# Save summary table
output_csv <- paste0(tolower(CHROMOSOME), "_kmer_frequency_summary_corrected.csv")
write_csv(summary_table, output_csv)
cat("✓ Saved:", output_csv, "\n")

cat("\n===============================================================================\n")
cat("✓ COMPLETE!\n")
cat("===============================================================================\n")
cat("\nOutput files:\n")
cat("  1. ", paste0(tolower(CHROMOSOME), "_kmer_frequency_boxplot_corrected.png"), " - Main boxplot comparison\n")
cat("  2. ", paste0(tolower(CHROMOSOME), "_kmer_individual_comparison_corrected.png"), " - Individual k-mer bars\n")
cat("  3. ", paste0(tolower(CHROMOSOME), "_kmer_frequency_summary_corrected.csv"), " - Summary statistics\n")
cat("\nAnalyzed", length(available_kmers), "k-mers across", length(common_samples), "samples\n")
cat("Resistant samples:", length(resistant_common), "\n")
cat("Susceptible samples:", length(susceptible_common), "\n")
cat("Using corrected k-mer sequences and sample groups\n")