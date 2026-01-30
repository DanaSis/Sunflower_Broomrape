#!/usr/bin/env Rscript
# ==============================================================================
# COMBINE ALL INTERSECTED FILES - ROBUST VERSION
# ==============================================================================

library(tidyverse)

setwd("/mnt/data/DanaS/sam_broom_gwas/gwas/heritability_analysis/kmer_novel/try")

cat("================================================================================\n")
cat("COMBINING INTERSECTED FILES INTO ALL_REF\n")
cat("================================================================================\n\n")

# ==============================================================================
# SETUP
# ==============================================================================

input_dir <- "/mnt/data/DanaS/sam_broom_gwas/significants_with_snp_mapping"

cat(sprintf("Input directory: %s\n\n", input_dir))

# ==============================================================================
# LOAD FILES ONE BY ONE
# ==============================================================================

cat("Loading files...\n\n")

all_data_list <- list()
file_count <- 0
empty_count <- 0
error_count <- 0

# Define all combinations
projects <- c("gadot", "yavor")
traits <- c("hel", "nec")
markers <- c("kmer", "snp")
refs <- c("HA412", "LR1", "PSC8", "XRQ2")

for(project in projects) {
  for(trait in traits) {
    for(marker in markers) {
      for(ref in refs) {
        
        # Construct filename
        filename <- sprintf("%s_%s_%s_%s.intersected250kbed_prod.tsv",
                           project, trait, marker, ref)
        filepath <- file.path(input_dir, filename)
        
        # Check if file exists
        if(!file.exists(filepath)) {
          cat(sprintf("SKIP: %s (file not found)\n", filename))
          next
        }
        
        # Try to load
        cat(sprintf("Loading: %s... ", filename))
        
        tryCatch({
          
          # Read with explicit settings
          data <- read.table(filepath,
                            header = TRUE,
                            sep = "\t",
                            quote = "",
                            comment.char = "",
                            stringsAsFactors = FALSE,
                            fill = TRUE,
                            check.names = FALSE)
          
          # Check if empty
          if(nrow(data) == 0) {
            cat("EMPTY\n")
            empty_count <- empty_count + 1
            next
          }
          
          cat(sprintf("%d rows\n", nrow(data)))
          
          # Add metadata
          data$project <- project
          data$trait <- trait
          data$marker <- marker
          data$ref <- ref
          
          # Store
          all_data_list[[length(all_data_list) + 1]] <- data
          file_count <- file_count + 1
          
        }, error = function(e) {
          cat(sprintf("ERROR: %s\n", e$message))
          error_count <- error_count + 1
        })
      }
    }
  }
}

cat("\n")
cat(sprintf("Summary: %d files loaded, %d empty, %d errors\n\n", 
            file_count, empty_count, error_count))

# ==============================================================================
# COMBINE
# ==============================================================================

if(length(all_data_list) == 0) {
  stop("ERROR: No files were loaded!")
}

cat("Combining datasets...\n")

all_ref <- bind_rows(all_data_list)

cat(sprintf("Total rows: %d\n", nrow(all_ref)))
cat(sprintf("Total columns: %d\n", ncol(all_ref)))
cat("\n")

# ==============================================================================
# CLEAN UP COLUMN NAMES
# ==============================================================================

cat("Checking and standardizing column names...\n\n")

cat("Current columns:\n")
print(colnames(all_ref))
cat("\n")

# Standardize marker column
if("marker" %in% colnames(all_ref)) {
  # This is the metadata column we added
  # Rename it to avoid confusion with data columns
  colnames(all_ref)[colnames(all_ref) == "marker"] <- "marker_type"
}

cat("After standardization:\n")
print(colnames(all_ref))
cat("\n")

# ==============================================================================
# SUMMARY
# ==============================================================================

cat("================================================================================\n")
cat("SUMMARY\n")
cat("================================================================================\n\n")

cat("By marker type:\n")
print(table(all_ref$marker_type))
cat("\n")

cat("By project:\n")
print(table(all_ref$project))
cat("\n")

cat("By trait:\n")
print(table(all_ref$trait))
cat("\n")

cat("By reference:\n")
print(table(all_ref$ref))
cat("\n")

# ==============================================================================
# SAVE
# ==============================================================================

cat("Saving output...\n")

output_file <- "ALL_REFS.all_markers.intersected250kbed_prod.tsv"

write.table(all_ref, output_file,
            sep = "\t", 
            row.names = FALSE, 
            quote = FALSE)

cat(sprintf("✓ Saved: %s\n", output_file))
cat(sprintf("  Rows: %d\n", nrow(all_ref)))
cat(sprintf("  Columns: %d\n", ncol(all_ref)))
cat(sprintf("  Size: %.2f MB\n\n", file.size(output_file) / 1024^2))

# ==============================================================================
# VALIDATION
# ==============================================================================

cat("================================================================================\n")
cat("VALIDATION\n")
cat("================================================================================\n\n")

cat("Checking for required columns...\n")

required <- c("mapped_k_chr", "mapped_k_pos", "mapped_k_gwas_sequence", 
              "project", "trait", "marker_type", "ref")

for(col in required) {
  if(col %in% colnames(all_ref)) {
    cat(sprintf("  ✓ %s\n", col))
  } else {
    cat(sprintf("  ✗ %s MISSING\n", col))
  }
}

cat("\nSample data:\n")
print(head(all_ref %>% 
           select(project, trait, marker_type, ref, 
                  mapped_k_chr, mapped_k_pos, mapped_k_gwas_sequence) %>%
           mutate(seq_preview = substr(mapped_k_gwas_sequence, 1, 20)), 5))

cat("\n================================================================================\n")
cat("✓ COMPLETE!\n")
cat("================================================================================\n")
