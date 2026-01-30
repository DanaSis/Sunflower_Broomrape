#!/usr/bin/env Rscript
# ==============================================================================
# STATISTICAL ANALYSIS OF PATHWAY DIFFERENCES BETWEEN RACES
# Addressing Reviewer Comment 1: Mechanistic differences
# R VERSION
# ==============================================================================
#setwd("/mnt/data/DanaS/sam_broom_gwas/gwas/pathway_analysis")
# Required libraries
required_packages <- c("tidyverse", "ggplot2", "reshape2", "pheatmap", "RColorBrewer")

for(pkg in required_packages) {
  if(!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg, repos = "http://cran.us.r-project.org")
    library(pkg, character.only = TRUE)
  }
}

cat("================================================================================\n")
cat("STATISTICAL ANALYSIS: PATHWAY DIFFERENCES BETWEEN RACES\n")
cat("================================================================================\n\n")

# ==============================================================================
# LOAD DATA
# ==============================================================================
setwd("/mnt/data/DanaS/sam_broom_gwas/gwas/pathway_analysis/files/")
cat("Loading data...\n")

results <- read.delim("./files/ALL_REFS_all_markers_intersected250kbed_prod.tsv",
                      header = TRUE, stringsAsFactors = FALSE)
go_enrichment <- read.delim("./files/enriched_GO.txt",
                            skip = 2, header = TRUE, stringsAsFactors = FALSE)
prioritized_genes <- read.delim("./prioritized_candidate_genes.tsv",
                                header = TRUE, stringsAsFactors = FALSE)

cat(sprintf("✓ Loaded %d marker-gene associations\n", nrow(results)))
cat(sprintf("✓ Loaded %d GO terms\n", nrow(go_enrichment)))
cat(sprintf("✓ Loaded %d prioritized genes\n", nrow(prioritized_genes)))

# ==============================================================================
# ANALYSIS 1: GO TERM ENRICHMENT COMPARISON
# ==============================================================================

cat("\n================================================================================\n")
cat("ANALYSIS 1: GO TERM ENRICHMENT COMPARISON\n")
cat("================================================================================\n\n")

# Filter for Biological Process
go_by_race <- go_enrichment %>%
  filter(ONTOLOGY == "BP") %>%
  select(Broomrape, Trait, ONTOLOGY, Description, Count, pvalue, p.adjust)

# Separate by race
yavor_go <- go_by_race %>% filter(Broomrape == "Yavor")
gadot_go <- go_by_race %>% filter(Broomrape == "Gadot")

cat(sprintf("Yavor enriched GO terms: %d\n", nrow(yavor_go)))
cat(sprintf("Gadot enriched GO terms: %d\n", nrow(gadot_go)))

# Find race-specific vs shared terms
yavor_terms <- unique(yavor_go$Description)
gadot_terms <- unique(gadot_go$Description)

yavor_specific_terms <- setdiff(yavor_terms, gadot_terms)
gadot_specific_terms <- setdiff(gadot_terms, yavor_terms)
shared_terms <- intersect(yavor_terms, gadot_terms)

cat(sprintf("\nYavor-specific GO terms: %d\n", length(yavor_specific_terms)))
cat(sprintf("Gadot-specific GO terms: %d\n", length(gadot_specific_terms)))
cat(sprintf("Shared GO terms: %d\n", length(shared_terms)))

# Chi-square test
contingency <- matrix(c(
  length(yavor_specific_terms), length(shared_terms),
  length(gadot_specific_terms), length(shared_terms)
), nrow = 2, byrow = TRUE)

chi_result <- chisq.test(contingency)

cat(sprintf("\n*** Chi-square test for GO term distribution ***\n"))
cat(sprintf("χ² = %.2f, df = %d, p = %.4f\n",
            chi_result$statistic, chi_result$parameter, chi_result$p.value))

if(chi_result$p.value < 0.05) {
  cat("Result: GO term distributions are SIGNIFICANTLY different between races\n")
} else {
  cat("Result: GO term distributions are not significantly different\n")
}

# ==============================================================================
# ANALYSIS 2: FUNCTIONAL CATEGORY ANALYSIS
# ==============================================================================

cat("\n================================================================================\n")
cat("ANALYSIS 2: FUNCTIONAL CATEGORY ANALYSIS\n")
cat("================================================================================\n\n")

# Define functional categories
functional_categories <- list(
  'Defense/Immunity' = c('defense', 'resistance', 'immune', 'pathogen', 'disease'),
  'Signal Transduction' = c('signal', 'kinase', 'phosphorylation', 'receptor'),
  'Transcription' = c('transcription', 'DNA-binding', 'gene expression'),
  'Cell Wall' = c('cell wall', 'cellulose', 'pectin', 'xyloglucan'),
  'Hormone Response' = c('hormone', 'auxin', 'jasmonate', 'salicylic', 'ethylene', 'abscisic'),
  'Stress Response' = c('stress', 'oxidative', 'reactive oxygen'),
  'Metabolism' = c('metabolic', 'biosynthetic', 'catabolic'),
  'Transport' = c('transport', 'transmembrane')
)

# Function to categorize GO terms
categorize_go_term <- function(description) {
  if(is.na(description)) return("Other")
  desc_lower <- tolower(as.character(description))
  
  categories <- c()
  for(cat_name in names(functional_categories)) {
    keywords <- functional_categories[[cat_name]]
    if(any(sapply(keywords, function(kw) grepl(kw, desc_lower, fixed = TRUE)))) {
      categories <- c(categories, cat_name)
    }
  }
  
  if(length(categories) == 0) return("Other")
  return(categories)
}

# Categorize GO terms
yavor_go$categories <- lapply(yavor_go$Description, categorize_go_term)
gadot_go$categories <- lapply(gadot_go$Description, categorize_go_term)

# Expand to one row per category
yavor_expanded <- yavor_go %>%
  unnest(categories)
gadot_expanded <- gadot_go %>%
  unnest(categories)

# Count by category
yavor_cat_counts <- table(yavor_expanded$categories)
gadot_cat_counts <- table(gadot_expanded$categories)

# Combine into comparison table
all_categories <- sort(unique(c(names(yavor_cat_counts), names(gadot_cat_counts))))

category_comparison <- data.frame(
  Category = all_categories,
  Yavor = sapply(all_categories, function(cat) 
    ifelse(cat %in% names(yavor_cat_counts), yavor_cat_counts[cat], 0)),
  Gadot = sapply(all_categories, function(cat) 
    ifelse(cat %in% names(gadot_cat_counts), gadot_cat_counts[cat], 0)),
  stringsAsFactors = FALSE
)

# Calculate proportions
category_comparison$Yavor_pct <- round(category_comparison$Yavor / sum(category_comparison$Yavor) * 100, 1)
category_comparison$Gadot_pct <- round(category_comparison$Gadot / sum(category_comparison$Gadot) * 100, 1)
category_comparison$Difference <- round(category_comparison$Yavor_pct - category_comparison$Gadot_pct, 1)

cat("Functional category enrichment by race:\n")
print(category_comparison, row.names = FALSE)
cat("\n")

# Fisher's exact tests for each category
cat("\n*** Fisher's exact tests for each functional category ***\n")
cat(paste(rep("-", 80), collapse = ""), "\n")

fisher_results <- data.frame()

for(cat in all_categories) {
  yavor_in_cat <- category_comparison$Yavor[category_comparison$Category == cat]
  yavor_out_cat <- sum(category_comparison$Yavor) - yavor_in_cat
  gadot_in_cat <- category_comparison$Gadot[category_comparison$Category == cat]
  gadot_out_cat <- sum(category_comparison$Gadot) - gadot_in_cat
  
  if(yavor_in_cat > 0 || gadot_in_cat > 0) {
    contingency <- matrix(c(
      yavor_in_cat, yavor_out_cat,
      gadot_in_cat, gadot_out_cat
    ), nrow = 2, byrow = TRUE)
    
    fisher_result <- fisher.test(contingency)
    
    fisher_results <- rbind(fisher_results, data.frame(
      Category = cat,
      Odds_Ratio = fisher_result$estimate,
      P_value = fisher_result$p.value,
      Significant = ifelse(fisher_result$p.value < 0.05, "Yes", "No"),
      stringsAsFactors = FALSE
    ))
    
    sig_marker <- ifelse(fisher_result$p.value < 0.001, "***",
                        ifelse(fisher_result$p.value < 0.01, "**",
                              ifelse(fisher_result$p.value < 0.05, "*", "")))
    
    cat(sprintf("%-25s OR=%6.2f, p=%.4f %s\n",
                cat, fisher_result$estimate, fisher_result$p.value, sig_marker))
  }
}

# Bonferroni correction
fisher_results <- fisher_results %>% arrange(P_value)
n_tests <- nrow(fisher_results)
fisher_results$Bonferroni <- pmin(fisher_results$P_value * n_tests, 1.0)
fisher_results$Bonferroni_Significant <- fisher_results$Bonferroni < 0.05

cat("\nAfter Bonferroni correction:\n")
sig_cats <- fisher_results %>% filter(Bonferroni_Significant == TRUE)
if(nrow(sig_cats) > 0) {
  print(sig_cats[, c("Category", "Odds_Ratio", "P_value", "Bonferroni")], row.names = FALSE)
} else {
  cat("No categories remain significant after Bonferroni correction\n")
}

# Save results
write.table(category_comparison, "./functional_category_comparison_R.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(fisher_results, "./fisher_tests_by_category_R.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)

# ==============================================================================
# ANALYSIS 3: R-GENE SUBCLASS ANALYSIS
# ==============================================================================

cat("\n================================================================================\n")
cat("ANALYSIS 3: R-GENE SUBCLASS DISTRIBUTION\n")
cat("================================================================================\n\n")

# Define R-gene subclasses
r_gene_classes <- list(
  'NBS-LRR' = c('NBS', 'NB-ARC', 'LRR', 'leucine-rich repeat', 'leucine rich repeat'),
  'RLK/RLP' = c('RLK', 'RLP', 'receptor-like kinase', 'receptor like kinase', 'receptor-like protein'),
  'Wall-associated' = c('WAK', 'wall-associated'),
  'PRR' = c('PRR', 'pattern recognition', 'PAMP'),
  'Other kinases' = c('kinase', 'serine/threonine', 'protein kinase')
)

# Function to classify R genes
classify_r_gene <- function(product) {
  if(is.na(product)) return("Other")
  
  product_lower <- tolower(as.character(product))
  
  # Check in order of specificity
  for(r_class in names(r_gene_classes)) {
    keywords <- r_gene_classes[[r_class]]
    if(any(sapply(keywords, function(kw) grepl(kw, product_lower, fixed = TRUE)))) {
      return(r_class)
    }
  }
  
  return("Other")
}

# Classify R genes
r_genes <- prioritized_genes %>%
  filter(is_r_gene == TRUE)

r_genes$r_gene_class <- sapply(r_genes$product, classify_r_gene)

# Count by race and R-gene class
r_class_counts <- r_genes %>%
  group_by(category, r_gene_class) %>%
  summarize(count = n(), .groups = "drop") %>%
  pivot_wider(names_from = r_gene_class, values_from = count, values_fill = 0)

cat("R-gene subclass distribution by race:\n")
print(r_class_counts, width = 100)
cat("\n")

# Calculate proportions
r_class_prop <- r_class_counts
r_class_prop[, -1] <- round(r_class_counts[, -1] / rowSums(r_class_counts[, -1]) * 100, 1)

cat("R-gene subclass proportions (%):\n")
print(r_class_prop, width = 100)
cat("\n")

# Fisher's exact tests for each R-gene class
cat("\n*** Chi-square tests for R-gene subclass distribution ***\n")
cat(paste(rep("-", 80), collapse = ""), "\n")

r_class_names <- setdiff(colnames(r_class_counts), c("category", "Other"))

for(r_class in r_class_names) {
  yavor_count <- r_class_counts[[r_class]][r_class_counts$category == "Yavor-specific"]
  gadot_count <- r_class_counts[[r_class]][r_class_counts$category == "Gadot-specific"]
  
  yavor_total <- sum(r_class_counts[r_class_counts$category == "Yavor-specific", -1])
  gadot_total <- sum(r_class_counts[r_class_counts$category == "Gadot-specific", -1])
  
  if(length(yavor_count) > 0 && length(gadot_count) > 0 && yavor_total > 0 && gadot_total > 0) {
    contingency <- matrix(c(
      yavor_count, yavor_total - yavor_count,
      gadot_count, gadot_total - gadot_count
    ), nrow = 2, byrow = TRUE)
    
    fisher_result <- fisher.test(contingency)
    
    sig_marker <- ifelse(fisher_result$p.value < 0.001, "***",
                        ifelse(fisher_result$p.value < 0.01, "**",
                              ifelse(fisher_result$p.value < 0.05, "*", "")))
    
    cat(sprintf("%-20s OR=%6.2f, p=%.4f %s\n",
                r_class, fisher_result$estimate, fisher_result$p.value, sig_marker))
  }
}

# Save results
write.table(r_class_counts, "./r_gene_subclass_counts_R.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(r_class_prop, "./r_gene_subclass_proportions_R.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)

# ==============================================================================
# ANALYSIS 4: TOP PATHWAYS
# ==============================================================================

cat("\n================================================================================\n")
cat("ANALYSIS 4: TOP 10 PATHWAYS PER RACE\n")
cat("================================================================================\n\n")

# Get top 10 for each race
yavor_top10 <- yavor_go %>%
  arrange(pvalue) %>%
  head(10) %>%
  select(Description, Count, pvalue, p.adjust)

gadot_top10 <- gadot_go %>%
  arrange(pvalue) %>%
  head(10) %>%
  select(Description, Count, pvalue, p.adjust)

cat("Top 10 pathways for YAVOR:\n")
print(yavor_top10, row.names = FALSE)
cat("\n")

cat("Top 10 pathways for GADOT:\n")
print(gadot_top10, row.names = FALSE)
cat("\n")

write.table(yavor_top10, "./yavor_top10_pathways_R.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(gadot_top10, "./gadot_top10_pathways_R.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)

# ==============================================================================
# CREATE FIGURES
# ==============================================================================

cat("\n================================================================================\n")
cat("CREATING FIGURES\n")
cat("================================================================================\n\n")

# Figure 1: Functional category comparison
pdf("./Fig_functional_category_analysis_R.pdf", width = 14, height = 6)
par(mfrow = c(1, 2), mar = c(10, 5, 4, 2))

# Panel A: Absolute counts
x_pos <- barplot(rbind(category_comparison$Yavor, category_comparison$Gadot),
                 beside = TRUE,
                 col = c("#E69F00", "#56B4E9"),
                 border = "black",
                 las = 2,
                 names.arg = category_comparison$Category,
                 ylab = "Number of Enriched GO Terms",
                 main = "A) Pathway Enrichment by Functional Category",
                 cex.names = 0.8)

legend("topright", legend = c("Yavor", "Gadot"),
       fill = c("#E69F00", "#56B4E9"), bty = "n")

# Panel B: Proportions
categories_for_plot <- category_comparison[category_comparison$Category != "Other", ]

barplot(t(as.matrix(categories_for_plot[, c("Yavor_pct", "Gadot_pct")])),
        names.arg = categories_for_plot$Category,
        horiz = TRUE,
        col = c("#E69F00", "#56B4E9"),
        border = "black",
        las = 1,
        xlab = "Proportion of Enriched GO Terms (%)",
        main = "B) Relative Pathway Distribution",
        cex.names = 0.7)

legend("bottomright", legend = c("Yavor", "Gadot"),
       fill = c("#E69F00", "#56B4E9"), bty = "n")

dev.off()

# Figure 2: R-gene subclass distribution
if(nrow(r_class_counts) > 0) {
  pdf("./Fig_R_gene_subclass_distribution_R.pdf", width = 10, height = 6)
  par(mar = c(10, 5, 4, 8), xpd = TRUE)
  
  r_class_matrix <- as.matrix(r_class_counts[, -1])
  rownames(r_class_matrix) <- r_class_counts$category
  
  barplot(r_class_matrix,
          beside = TRUE,
          col = c("#E69F00", "#56B4E9"),
          border = "black",
          las = 2,
          ylab = "Number of R genes",
          main = "R-gene Subclass Distribution by Race",
          cex.names = 0.8)
  
  legend(par("usr")[2], par("usr")[4],
         legend = r_class_counts$category,
         fill = c("#E69F00", "#56B4E9"),
         bty = "n")
  
  dev.off()
}

# Figure 3: Heatmap
pdf("./Fig_pathway_heatmap_R.pdf", width = 10, height = 8)

heatmap_data <- category_comparison[, c("Category", "Yavor", "Gadot")]
rownames(heatmap_data) <- heatmap_data$Category
heatmap_data$Category <- NULL
heatmap_data <- as.matrix(heatmap_data)

# Normalize by row
heatmap_norm <- heatmap_data / rowSums(heatmap_data) * 100

heatmap(heatmap_norm,
        Colv = NA,
        Rowv = NA,
        scale = "none",
        col = colorRampPalette(c("blue", "white", "red"))(50),
        margins = c(8, 12),
        main = "Functional Category Enrichment Heatmap\n(% within each category)")

dev.off()

cat("✓ Created Fig_functional_category_analysis_R.pdf\n")
cat("✓ Created Fig_R_gene_subclass_distribution_R.pdf\n")
cat("✓ Created Fig_pathway_heatmap_R.pdf\n")

# ==============================================================================
# SUMMARY
# ==============================================================================

cat("\n================================================================================\n")
cat("ANALYSIS COMPLETE - SUMMARY\n")
cat("================================================================================\n\n")

yavor_defense_terms <- sum(yavor_expanded$categories == "Defense/Immunity")
gadot_defense_terms <- sum(gadot_expanded$categories == "Defense/Immunity")

cat(sprintf("Yavor: %d/%d defense-related pathways (%.1f%%)\n",
            yavor_defense_terms, nrow(yavor_go),
            yavor_defense_terms/nrow(yavor_go)*100))

cat(sprintf("Gadot: %d/%d defense-related pathways (%.1f%%)\n",
            gadot_defense_terms, nrow(gadot_go),
            gadot_defense_terms/nrow(gadot_go)*100))

cat("\n================================================================================\n")
cat("OUTPUT FILES CREATED:\n")
cat("================================================================================\n")
cat("  1. functional_category_comparison_R.tsv\n")
cat("  2. fisher_tests_by_category_R.tsv\n")
cat("  3. r_gene_subclass_counts_R.tsv\n")
cat("  4. r_gene_subclass_proportions_R.tsv\n")
cat("  5. yavor_top10_pathways_R.tsv\n")
cat("  6. gadot_top10_pathways_R.tsv\n")
cat("  7. Fig_functional_category_analysis_R.pdf\n")
cat("  8. Fig_R_gene_subclass_distribution_R.pdf\n")
cat("  9. Fig_pathway_heatmap_R.pdf\n")
cat("\nAnalysis complete!\n")
