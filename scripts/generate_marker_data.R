#!/usr/bin/env Rscript
# ==============================================================================
# GENERATE MARKER_DATA.TXT FROM K-MER AND SNP FILES
# Creates comprehensive marker file for novelty and heritability analyses
# ==============================================================================
setwd("/mnt/data/DanaS/sam_broom_gwas/gwas/heritability_analysis/kmer_novel/try")
library(tidyverse)
library(data.table)
cat("================================================================================\n")
cat("GENERATING MARKER_DATA.TXT FROM SOURCE FILES\n")
cat("================================================================================\n\n")

# ==============================================================================
# DEFINE FILE PATHS
# ==============================================================================

cat("Step 1: Defining file paths...\n\n")

# K-mer files
kmer_files <- c(
  "/mnt/data/DanaS/sam_broom_gwas/gwas/gadot_hel_thresh_05/kmer/get_kmers/ha412/kmers_manhattan_plot_data.txt",
  "/mnt/data/DanaS/sam_broom_gwas/gwas/gadot_hel_thresh_05/kmer/get_kmers/lr1/kmers_manhattan_plot_data.txt",
  "/mnt/data/DanaS/sam_broom_gwas/gwas/gadot_hel_thresh_05/kmer/get_kmers/psc8/kmers_manhattan_plot_data.txt",
  "/mnt/data/DanaS/sam_broom_gwas/gwas/gadot_hel_thresh_05/kmer/get_kmers/xrq/kmers_manhattan_plot_data.txt",
  "/mnt/data/DanaS/sam_broom_gwas/gwas/gadot_nec_thresh_03/kmer/get_kmers/ha412/kmers_manhattan_plot_data.txt",
  "/mnt/data/DanaS/sam_broom_gwas/gwas/gadot_nec_thresh_03/kmer/get_kmers/lr1/kmers_manhattan_plot_data.txt",
  "/mnt/data/DanaS/sam_broom_gwas/gwas/gadot_nec_thresh_03/kmer/get_kmers/psc8/kmers_manhattan_plot_data.txt",
  "/mnt/data/DanaS/sam_broom_gwas/gwas/gadot_nec_thresh_03/kmer/get_kmers/xrq/kmers_manhattan_plot_data.txt",
  "/mnt/data/DanaS/sam_broom_gwas/gwas/yavor_hel_thresh_05/kmer/get_kmers/ha412/kmers_manhattan_plot_data_pvast.txt",
  "/mnt/data/DanaS/sam_broom_gwas/gwas/yavor_hel_thresh_05/kmer/get_kmers/lr1/kmers_manhattan_plot_data_pvast.txt",
  "/mnt/data/DanaS/sam_broom_gwas/gwas/yavor_hel_thresh_05/kmer/get_kmers/psc8/kmers_manhattan_plot_data_pvast.txt",
  "/mnt/data/DanaS/sam_broom_gwas/gwas/yavor_hel_thresh_05/kmer/get_kmers/xrq/kmers_manhattan_plot_data_pvast.txt",
  "/mnt/data/DanaS/sam_broom_gwas/gwas/yavor_nec_thresh_03/kmer/get_kmers/ha412/kmers_manhattan_plot_data_pvast.txt",
  "/mnt/data/DanaS/sam_broom_gwas/gwas/yavor_nec_thresh_03/kmer/get_kmers/lr1/kmers_manhattan_plot_data_pvast.txt",
  "/mnt/data/DanaS/sam_broom_gwas/gwas/yavor_nec_thresh_03/kmer/get_kmers/psc8/kmers_manhattan_plot_data_pvast.txt",
  "/mnt/data/DanaS/sam_broom_gwas/gwas/yavor_nec_thresh_03/kmer/get_kmers/xrq/kmers_manhattan_plot_data_pvast.txt"
)

# SNP files
# snp_files <- c(
#   "/mnt/data/DanaS/sam_broom_gwas/gwas/gadot_hel_thresh_05/work_lift/ann/HA412.snps.bed",
#   "/mnt/data/DanaS/sam_broom_gwas/gwas/gadot_hel_thresh_05/work_lift/ann/LR1.snps.bed",
#   "/mnt/data/DanaS/sam_broom_gwas/gwas/gadot_hel_thresh_05/work_lift/ann/PSC8.snps.bed",
#   "/mnt/data/DanaS/sam_broom_gwas/gwas/gadot_hel_thresh_05/work_lift/ann/XRQ2.snps.bed",
#   "/mnt/data/DanaS/sam_broom_gwas/gwas/gadot_hel_thresh_05/work_lift/ann/XRQ.snps.bed",
#   "/mnt/data/DanaS/sam_broom_gwas/gwas/gadot_nec_thresh_03/work_lift/ann/HA412.snps.bed",
#   "/mnt/data/DanaS/sam_broom_gwas/gwas/gadot_nec_thresh_03/work_lift/ann/LR1.snps.bed",
#   "/mnt/data/DanaS/sam_broom_gwas/gwas/gadot_nec_thresh_03/work_lift/ann/PSC8.snps.bed",
#   "/mnt/data/DanaS/sam_broom_gwas/gwas/gadot_nec_thresh_03/work_lift/ann/XRQ2.snps.bed",
#   "/mnt/data/DanaS/sam_broom_gwas/gwas/yavor_hel_thresh_05/work_lift/ann/HA412.snps.bed",
#   "/mnt/data/DanaS/sam_broom_gwas/gwas/yavor_hel_thresh_05/work_lift/ann/LR1.snps.bed",
#   "/mnt/data/DanaS/sam_broom_gwas/gwas/yavor_hel_thresh_05/work_lift/ann/PSC8.snps.bed",
#   "/mnt/data/DanaS/sam_broom_gwas/gwas/yavor_hel_thresh_05/work_lift/ann/XRQ2.snps.bed",
#   "/mnt/data/DanaS/sam_broom_gwas/gwas/yavor_nec_thresh_03/work_lift/ann/HA412.snps.bed",
#   "/mnt/data/DanaS/sam_broom_gwas/gwas/yavor_nec_thresh_03/work_lift/ann/LR1.snps.bed",
#   "/mnt/data/DanaS/sam_broom_gwas/gwas/yavor_nec_thresh_03/work_lift/ann/PSC8.snps.bed",
#   "/mnt/data/DanaS/sam_broom_gwas/gwas/yavor_nec_thresh_03/work_lift/ann/XRQ2.snps.bed"
# )

snp_files <- c(
"/mnt/data/DanaS/sam_broom_gwas/gwas/gadot_nec_thresh_03/work_lift_bwa100/HanLR1r0.9-20211115.genome.lift.raw.tsv",
"/mnt/data/DanaS/sam_broom_gwas/gwas/gadot_nec_thresh_03/work_lift_bwa100/HanPSC8r1.0-20181105.genome.lift.raw.tsv",
"/mnt/data/DanaS/sam_broom_gwas/gwas/gadot_nec_thresh_03/work_lift_bwa100/HanXRQr2.0-SUNRISE-2.1.genome.lift.raw.tsv",
"/mnt/data/DanaS/sam_broom_gwas/gwas/gadot_hel_thresh_05/work_lift_bwa100/HanLR1r0.9-20211115.genome.lift.raw.tsv",
"/mnt/data/DanaS/sam_broom_gwas/gwas/gadot_hel_thresh_05/work_lift_bwa100/HanPSC8r1.0-20181105.genome.lift.raw.tsv",
"/mnt/data/DanaS/sam_broom_gwas/gwas/gadot_hel_thresh_05/work_lift_bwa100/HanXRQr2.0-SUNRISE-2.1.genome.lift.raw.tsv",
"/mnt/data/DanaS/sam_broom_gwas/gwas/yavor_nec_thresh_03/work_lift_bwa100/HanLR1r0.9-20211115.genome.lift.raw.tsv",
"/mnt/data/DanaS/sam_broom_gwas/gwas/yavor_nec_thresh_03/work_lift_bwa100/HanPSC8r1.0-20181105.genome.lift.raw.tsv",
"/mnt/data/DanaS/sam_broom_gwas/gwas/yavor_nec_thresh_03/work_lift_bwa100/HanXRQr2.0-SUNRISE-2.1.genome.lift.raw.tsv",
"/mnt/data/DanaS/sam_broom_gwas/gwas/yavor_hel_thresh_05/work_lift_bwa100/HanLR1r0.9-20211115.genome.lift.raw.tsv",
"/mnt/data/DanaS/sam_broom_gwas/gwas/yavor_hel_thresh_05/work_lift_bwa100/HanPSC8r1.0-20181105.genome.lift.raw.tsv",
"/mnt/data/DanaS/sam_broom_gwas/gwas/yavor_hel_thresh_05/work_lift_bwa100/HanXRQr2.0-SUNRISE-2.1.genome.lift.raw.tsv")

snp_files_ha412 <- c(
  "/mnt/data/DanaS/sam_broom_gwas/gwas/gadot_hel_thresh_05/work_lift_bwa100/sig.1bp.meta.bed",
  "/mnt/data/DanaS/sam_broom_gwas/gwas/yavor_hel_thresh_05/work_lift_bwa100/sig.1bp.meta.bed",
  "/mnt/data/DanaS/sam_broom_gwas/gwas/yavor_nec_thresh_03/work_lift_bwa100/sig.1bp.meta.bed"
)

# Check which files exist
kmer_exists <- file.exists(kmer_files)
snp_exists <- file.exists(snp_files)

cat("K-mer files found:", sum(kmer_exists), "of", length(kmer_files), "\n")
cat("SNP files found:", sum(snp_exists), "of", length(snp_files), "\n\n")

# Filter to existing files
kmer_files <- kmer_files[kmer_exists]
snp_files <- snp_files[snp_exists]

if(length(kmer_files) == 0 && length(snp_files) == 0) {
  stop("ERROR: No input files found! Check file paths.")
}

# ==============================================================================
# HELPER FUNCTION: PARSE METADATA FROM FILEPATH
# ==============================================================================

parse_meta_from_path <- function(file) {
  # Extract source (gadot / yavor)
  source <- ifelse(grepl("gadot", file, ignore.case=TRUE), "gadot", "yavor")
  
  # Extract trait (hel / nec)
  trait <- ifelse(grepl("_nec_", file, ignore.case=TRUE), "nec", "hel")
  
  # Extract reference genome
  reference <- case_when(
    
    grepl("sig.1bp", file, ignore.case=TRUE) ~ "ha412",
    grepl("HA412|ha412", file, ignore.case=TRUE) ~ "ha412",
    grepl("LR1|lr1", file, ignore.case=TRUE) ~ "lr1",
    grepl("PSC8|psc8", file, ignore.case=TRUE) ~ "psc8",
    grepl("XRQ|xrq", file, ignore.case=TRUE) ~ "xrq",
    TRUE ~ NA_character_
  )
  
  list(source = source, trait = trait, reference = reference)
}

# ==============================================================================
# LOAD K-MER DATA
# ==============================================================================

cat("================================================================================\n")
cat("Step 2: Loading K-mer files...\n")
cat("================================================================================\n\n")

kmer_data_list <- list()

for(file in kmer_files) {
  cat(sprintf("Loading: %s\n", basename(file)))
  
  # Parse metadata
  meta <- parse_meta_from_path(file)
  
  # Read file
  # Expected columns: chr, pos, pval (and possibly others)
  tryCatch({
    data <- read.delim(file, header = TRUE, stringsAsFactors = FALSE)
    
    # Check if required columns exist
    if(!all(c("chr", "pos", "mapped_marker", "pval") %in% colnames(data))) {
      cat("  WARNING: Missing required columns. Trying alternative names...\n")
      # Try common alternatives
      if("chromosome" %in% colnames(data)) colnames(data)[colnames(data) == "chromosome"] <- "chr"
      if("position" %in% colnames(data)) colnames(data)[colnames(data) == "position"] <- "pos"
      if("p_value" %in% colnames(data)) colnames(data)[colnames(data) == "p_value"] <- "pval"
      if("P" %in% colnames(data)) colnames(data)[colnames(data) == "P"] <- "pval"
      if("gwas_sequence" %in% colnames(data)) colnames(data)[colnames(data) == "gwas_sequence"] <- "mapped_marker"
      
    }
    
    # Add metadata
    data$source <- meta$source
    data$type <- meta$trait
    data$reference <- meta$reference
    data$mapped_marker <- data$mapped_marker
    data$marker <- "kmers"
    
    # Select and standardize columns
    data <- data %>%
      select(chr, pos, source, type, reference, marker, mapped_marker, pval) %>%
      filter(!is.na(chr) & !is.na(pos) & !is.na(pval))
    
    cat(sprintf("  Loaded %d k-mer positions\n", nrow(data)))
    
    kmer_data_list[[length(kmer_data_list) + 1]] <- data
    
  }, error = function(e) {
    cat(sprintf("  ERROR loading file: %s\n", e$message))
  })
}

cat("\n")

# Combine all k-mer data
if(length(kmer_data_list) > 0) {
  kmer_data <- bind_rows(kmer_data_list)
  cat(sprintf("Total k-mer positions loaded: %d\n\n", nrow(kmer_data)))
} else {
  kmer_data <- data.frame()
  cat("WARNING: No k-mer data loaded!\n\n")
}

# ==============================================================================
# LOAD SNP DATA
# ==============================================================================

cat("================================================================================\n")
cat("Step 3: Loading SNP files...\n")
cat("================================================================================\n\n")

snp_data_list <- list()

for(file in snp_files) {
  cat(sprintf("Loading: %s\n", basename(file)))
  
  # Parse metadata
  meta <- parse_meta_from_path(file)
  
  # Read BED file
  # BED format: chr, start, end, [optional: name, score, strand, ...]
  tryCatch({
    data <- fread(file, header = FALSE, stringsAsFactors = FALSE)
    
    # BED files typically have no header
    # Assume: col1=chr, col2=start, col3=end
    if(ncol(data) >= 3) {
      colnames(data)[1:6] <- c("chr", "start", "mapped_412marker", "alt", "ref", "strand")
      
      # Use start position as SNP position
      data$pos <- data$start
      data$mapped_marker <- paste0(data$chr, "_", data$start, "_", data$alt, "_", data$ref)

      # If there's a p-value column, try to find it
      # Often in column 5 or 6
      if(ncol(data) >= 5) {
        # Try to identify p-value column (contains small numbers)
        for(i in 4:ncol(data)) {
          if(is.numeric(data[[i]]) && all(data[[i]][!is.na(data[[i]])] <= 1)) {
            data$pval <- data[[i]]
            break
          }
        }
      }
      
      # If no p-value found, set to NA (will need to add separately)
      if(!"pval" %in% colnames(data)) {
        cat("  WARNING: No p-value column detected. Setting pval = NA\n")
        data$pval <- NA
      }
      
    } else {
      cat("  ERROR: BED file has fewer than 3 columns\n")
      next
    }
    
    # Add metadata
    data$source <- meta$source
    data$type <- meta$trait
    data$reference <- meta$reference
    data$marker <- "snp"
    data$mapped_marker <- data$mapped_marker
    
    # Select and standardize columns
    data <- data %>%
      select(chr, pos, source, type, reference, marker, mapped_marker, pval) %>%
      filter(!is.na(chr) & !is.na(pos))
    
    cat(sprintf("  Loaded %d SNP positions\n", nrow(data)))
    
    snp_data_list[[length(snp_data_list) + 1]] <- data
    
  }, error = function(e) {
    cat(sprintf("  ERROR loading file: %s\n", e$message))
  })
}


snp_data_list_412 <- list()

for(file in snp_files_ha412) {
  cat(sprintf("Loading: %s\n", basename(file)))
  
  # Parse metadata
  meta <- parse_meta_from_path(file)
  
  # Read BED file
  # BED format: chr, start, end, [optional: name, score, strand, ...]
  tryCatch({
    data <- read.delim(file, header = FALSE, stringsAsFactors = FALSE)
    
    # BED files typically have no header
    # Assume: col1=chr, col2=start, col3=end
    if(ncol(data) >= 3) {
      
      colnames(data)[1:6] <- c("chr", "start", "end", "mapped_marker", "alt", "ref")
      
      # Use start position as SNP position
      data$pos <- data$start
      data$mapped_marker <- data$mapped_marker
      # If there's a p-value column, try to find it
      # Often in column 5 or 6
      if(ncol(data) >= 5) {
        # Try to identify p-value column (contains small numbers)
        for(i in 4:ncol(data)) {
          if(is.numeric(data[[i]]) && all(data[[i]][!is.na(data[[i]])] <= 1)) {
            data$pval <- data[[i]]
            break
          }
        }
      }
      
      # If no p-value found, set to NA (will need to add separately)
      if(!"pval" %in% colnames(data)) {
        cat("  WARNING: No p-value column detected. Setting pval = NA\n")
        data$pval <- NA
      }
      
    } else {
      cat("  ERROR: BED file has fewer than 3 columns\n")
      next
    }
    
    # Add metadata
    data$source <- meta$source
    data$type <- meta$trait
    data$reference <- meta$reference
    data$marker <- "snp"
    data$mapped_marker <- data$mapped_marker
    
    # Select and standardize columns
    data_412 <- data %>%
      select(chr, pos, source, type, reference, marker, mapped_marker, pval) %>%
      filter(!is.na(chr) & !is.na(pos))
    
    cat(sprintf("  Loaded %d SNP positions\n", nrow(data_412)))
    
    snp_data_list_412[[length(snp_data_list_412) + 1]] <- data_412
    
  }, error = function(e) {
    cat(sprintf("  ERROR loading file: %s\n", e$message))
  })
}




cat("\n")

# Combine all SNP data
# if(length(snp_data_list) > 0) {
  snp_data <- bind_rows(snp_data_list)
#   cat(sprintf("Total SNP positions loaded: %d\n\n", nrow(snp_data)))
# } else {
#   snp_data <- data.frame()
#   cat("WARNING: No SNP data loaded!\n\n")
# }

snp_data_412 <- bind_rows(snp_data_list_412)

snp_data <- rbind(snp_data, snp_data_412)

# ==============================================================================
# COMBINE AND FILTER
# ==============================================================================

cat("================================================================================\n")
cat("Step 4: Combining and filtering markers...\n")
cat("================================================================================\n\n")

# Combine k-mers and SNPs
if(nrow(kmer_data) > 0 && nrow(snp_data) > 0) {
  all_markers <- bind_rows(kmer_data, snp_data)
} else if(nrow(kmer_data) > 0) {
  all_markers <- kmer_data
  cat("WARNING: Only k-mer data available (no SNPs)\n")
} else if(nrow(snp_data) > 0) {
  all_markers <- snp_data
  cat("WARNING: Only SNP data available (no k-mers)\n")
} else {
  stop("ERROR: No data to combine!")
}

cat(sprintf("Combined markers: %d\n", nrow(all_markers)))

# Filter by significance threshold (if p-values available)
if(!all(is.na(all_markers$pval))) {
  threshold <- 5e-8
  significant <- all_markers %>% filter(pval < threshold | is.na(pval))
  
  cat(sprintf("Significant markers (p < %.2e): %d\n", threshold, 
              sum(!is.na(significant$pval))))
  cat(sprintf("Markers without p-value: %d\n", sum(is.na(significant$pval))))
  
  all_markers <- significant
} else {
  cat("WARNING: No p-values available, keeping all markers\n")
}

# ==============================================================================
# SUMMARY STATISTICS
# ==============================================================================

cat("\n================================================================================\n")
cat("Step 5: Summary Statistics\n")
cat("================================================================================\n\n")

cat("Overall Summary:\n")
cat("----------------\n")
cat(sprintf("Total markers: %d\n", nrow(all_markers)))
cat(sprintf("  K-mers: %d (%.1f%%)\n", 
            sum(all_markers$marker == "kmers"),
            sum(all_markers$marker == "kmers") / nrow(all_markers) * 100))
cat(sprintf("  SNPs: %d (%.1f%%)\n", 
            sum(all_markers$marker == "snp"),
            sum(all_markers$marker == "snp") / nrow(all_markers) * 100))
cat("\n")

cat("By Source (Race):\n")
print(table(all_markers$source, all_markers$marker))
cat("\n")

cat("By Trait:\n")
print(table(all_markers$type, all_markers$marker))
cat("\n")

cat("By Reference:\n")
print(table(all_markers$reference, all_markers$marker))
cat("\n")

cat("By Source-Trait-Marker:\n")
summary_table <- all_markers %>%
  group_by(source, type, marker) %>%
  summarize(
    N = n(),
    N_chr = n_distinct(chr),
    Min_pval = min(pval, na.rm = TRUE),
    .groups = "drop"
  )
print(summary_table, n = 100)

# ==============================================================================
# SAVE OUTPUT
# ==============================================================================

cat("\n================================================================================\n")
cat("Step 6: Saving output...\n")
cat("================================================================================\n\n")

# Sort by chromosome and position
all_markers <- all_markers %>%
  arrange(chr, pos)
ha412_markers <- all_markers %>%
  filter(reference == "ha412") %>%
  arrange(chr, pos)
# Save to file
output_file <- "marker_data.txt"
write.table(all_markers, output_file, 
            sep = "\t", row.names = FALSE, quote = FALSE)
ha412_markers_file <- "marker_data_ha412.txt"
write.table(ha412_markers, ha412_markers_file, 
            sep = "\t", row.names = FALSE, quote = FALSE)
cat(sprintf("✓ Saved marker_data.txt\n"))
cat(sprintf("  File: %s\n", output_file))
cat(sprintf("  Size: %d rows, %d columns\n", nrow(all_markers), ncol(all_markers)))
cat(sprintf("  Columns: %s\n", paste(colnames(all_markers), collapse = ", ")))

# Also save to outputs directory if it exists
if(dir.exists("./")) {
  output_file2 <- "./marker_data.txt"
  write.table(all_markers, output_file2, 
              sep = "\t", row.names = FALSE, quote = FALSE)
  cat(sprintf("✓ Also saved to: %s\n", output_file2))
}

# ==============================================================================
# QUICK VALIDATION
# ==============================================================================

cat("\n================================================================================\n")
cat("Step 7: Quick validation...\n")
cat("================================================================================\n\n")

# Check for missing values
cat("Missing values:\n")
cat(sprintf("  chr: %d\n", sum(is.na(all_markers$chr))))
cat(sprintf("  pos: %d\n", sum(is.na(all_markers$pos))))
cat(sprintf("  source: %d\n", sum(is.na(all_markers$source))))
cat(sprintf("  type: %d\n", sum(is.na(all_markers$type))))
cat(sprintf("  reference: %d\n", sum(is.na(all_markers$reference))))
cat(sprintf("  marker: %d\n", sum(is.na(all_markers$marker))))
cat(sprintf("  pval: %d\n", sum(is.na(all_markers$pval))))
cat("\n")

# Check chromosome names
cat("Unique chromosomes:\n")
chr_counts <- sort(table(all_markers$chr), decreasing = TRUE)
print(head(chr_counts, 20))
cat(sprintf("  Total unique: %d\n", length(unique(all_markers$chr))))
cat("\n")

# Check position ranges by chromosome
cat("Position ranges (top 5 chromosomes):\n")
pos_ranges <- all_markers %>%
  group_by(chr) %>%
  summarize(
    N = n(),
    Min_pos = min(pos),
    Max_pos = max(pos),
    Range_Mb = (max(pos) - min(pos)) / 1e6,
    .groups = "drop"
  ) %>%
  arrange(desc(N))

print(head(pos_ranges, 5))
cat("\n")

# ==============================================================================
# COMPLETION
# ==============================================================================

cat("================================================================================\n")
cat("✓ MARKER_DATA.TXT GENERATION COMPLETE!\n")
cat("================================================================================\n\n")

cat("Next steps:\n")
cat("  1. Review marker_data.txt\n")
cat("  2. Run k-mer novelty analysis:\n")
cat("     Rscript kmer_novelty_complete_analysis.R\n")
cat("  3. Run heritability analysis:\n")
cat("     Rscript heritability_snp_kmer_combined.R\n\n")

cat("Files ready for:\n")
cat("  ✓ K-mer novelty analysis\n")
cat("  ✓ Heritability comparison\n")
cat("  ✓ LD analysis\n")
cat("  ✓ Reviewer response\n\n")

cat("================================================================================\n")
