# conda deactivate - deactivate base
# R

library(RColorBrewer)
library(readr)      # For read_tsv if needed
library(DESeq2)
library(tidyverse)
library(ggplot2)
library(pheatmap)
#library(tximeta)
library(data.table)
library(ggrepel)
library(EnhancedVolcano)
library(tximport)
library(grid)

library(gridExtra)
# ------------------------------
# 1. Create Metadata Data Frame
# ------------------------------
setwd("/mnt/data/DanaS/rna/realign_ha412/deseq")
# Define factors for experimental design:
accessions <- c("s", "r", "e")
names <- c("s", "r", "emek3")
# Conditions: a = time 0 (no treatment), b = time 5 untreated, c = time 5 treated.
conditions <- c("a", "b", "c")
replicates <- 1:3  # three replicates per combination

# Create a data frame with all combinations
metadata <- expand.grid(
  accession = accessions,
  condition = conditions,
  replicate = replicates,
  stringsAsFactors = FALSE
)

# Order the metadata for readability (optional)
metadata <- metadata[order(metadata$accession, metadata$condition, metadata$replicate), ]

# Create sample names by concatenating accession, condition, and replicate number:
metadata$sample <- with(metadata, paste0(accession, condition, replicate))

# Add time and treatment columns based on the condition letter:
metadata$time <- ifelse(metadata$condition == "a", 0, 5)
metadata$treatment <- ifelse(metadata$condition == "a", "no_treatment",
                      ifelse(metadata$condition == "b", "untreated", "treated"))
metadata$name <- ifelse(metadata$accession == "e", "emek3",
                      ifelse(metadata$accession == "r", "r", "s"))                      

#metadata$treatment <- ifelse(metadata$condition == "b", "untreated", "treated")
# Reorder columns (optional)
metadata <- metadata[, c("sample", "accession", "condition", "time", "treatment", "replicate", "name")]

# View the metadata data frame:
print(metadata)

# ------------------------------
# 2. Import Quantification Data with tximport
# ------------------------------
# Assume that for each sample, tximport produced a file named "sample.genes.results"
# stored in a directory. Adjust the path as necessary.
# Here we assume that the files are in the "output/results/rsem/" directory.

# Get the vector of sample names from the metadata
samples <- metadata$sample

# Define the directory where your quantification files are stored
tx_dir <- "/mnt/data/DanaS/rna/realign_ha412/output/results/output/output-data/"  # change if necessary

# Create a named vector of file paths
files <- file.path(tx_dir, paste0("rsem_star_", samples, ".genes.results"))
names(files) <- samples

# (Optional) Check that files exist:
if(!all(file.exists(files))){
  stop("Some quantification files are missing. Please check your file paths!")
}

metadata$accession <- as.factor(metadata$accession)
metadata$treatment <- as.factor(metadata$treatment)
metadata$condition <- as.factor(metadata$condition)
metadata$sample <- as.factor(metadata$sample)
metadata$name <- as.factor(metadata$name)
# Import the quantification files with tximport.
# For RSEM, the file usually contains a column named "expected_count" (or similar).
# You may need to adjust the column names via the tximport parameter "txOut" or by renaming
# columns later. Here we assume gene-level counts are directly provided.

txi.rsem <- tximport(files, type = "rsem", txIn = FALSE, txOut = FALSE)
head(txi.rsem$counts)
tail(txi.rsem$counts)
# fix "Error: all(lengths > 0) is not TRUE" error
txi.rsem$length[txi.rsem$length == 0] = 0.01
# ------------------------------
# 3. Create DESeq2 Dataset Using tximport Data
# ------------------------------

# Create DESeq2 metadata: the row names must match the column names in txi.rsem.
# Here we assume that the order of samples in txi.rsem is the same as in metadata.
rownames(metadata) <- metadata$sample

ddstxi <- DESeqDataSetFromTximport(txi.rsem,
                                colData = metadata,
                                design = ~ accession + treatment)
dds <- DESeq(ddstxi)
resultsNames(dds)
dds$group <- factor(paste0(dds$name, dds$treatment))
design(dds) <- ~ group
dds <- DESeq(dds)
resultsNames(dds)

res <- results(dds, contrast=c("group", "emek3treated", "emek3untreated"))
res_p05 <- res %>% as.data.frame() %>% 
            na.omit(.)  %>% 
            filter(padj < 0.05) %>% 
            mutate(gene = rownames(.))
marker_data <- fread("/mnt/data/DanaS/sam_broom_gwas/significants1/ha412_marker_data.txt", 
                     header = TRUE)




res_have <- res %>% as.data.frame() %>% 
            filter(!grepl("gene:Ha412HOChr00c", rownames(res))) %>% 
            na.omit(.)  %>% 
            filter(padj < 0.05) %>% 
            mutate(gene = rownames(.))

res_have <- res_have[c("gene", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")]

res_have$gene <- gsub("gene:", "", res_have$gene)

deg_gwas_com_chr_9 <- left_join(res_have, marker_data, by = c("gene" = "ID")) %>% 
            filter(chromosome == "Ha412HOChr09") %>% 
            select(gene, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj, chromosome, start, end, product) %>% 
            filter(start >= 167075448 | end <= 167691607) 

deg_gwas_com_chr_9_resistance <- deg_gwas_com_chr_9[c(1,4:7),] %>% 
            select(gene, product) %>%
            unique(.)

# Ha412HOChr09:167075448-167691607
write.table(as.data.frame(res_have),
            file = "/mnt/data/DanaS/rna/realign_ha412/deseq/emek_treatment_vs_untreated_new.txt",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)

marker_data <- fread("/mnt/data/DanaS/sam_broom_gwas/significants_with_snp_mapping/combined_by_ref/HA412.all_markers.intersected250kbed_prod.tsv", 
                     header = TRUE)
res <- fread("/mnt/data/DanaS/rna/realign_ha412/deseq/group_euntreated_vs_treatment.txt", 
                     header = TRUE)



rownames(res) <- gsub("gene:", "", rownames(res))
rownames(dds) <- gsub("gene:", "", rownames(dds))
# Make sure marker_data$ID matches rownames of ec_eb
marker_gene_ids <- marker_data$gene_id
#deg_filtered <- rownames(res) %in% marker_gene_ids
deg_filtered <- as.data.frame(res) %>% 
  filter(rownames(res) %in% marker_gene_ids) %>% 
  filter(padj < 0.05)
vsd <- vst(dds, blind = FALSE)  # or rlog() if needed
rownames(vsd) <- gsub("gene:", "", rownames(vsd))
norm_mat <- assay(vsd)[rownames(vsd) %in% rownames(deg_filtered), ]
norm_mat <- norm_mat[,4:9]
ann_col <- as.data.frame(colData(vsd)[, c("treatment", "name")])
an_col <-ann_col[4:9,]
rownames(an_col) <- c("ctrl1","ctrl2","ctrl3","trt1","trt2","trt3")
an_col$trt = ifelse(an_col$treatment=="treated","infected", "non-infected")
ann_col <- c("ctrl1","ctrl2","ctrl3","trt1","trt2","trt3")
colnames(norm_mat) <- c("ctrl1","ctrl2","ctrl3","trt1","trt2","trt3")

jpeg("gwas_deg_col.jpeg", width=7, height=7, units = "in", res = 300)
pheatmap(norm_mat,
         annotation_col = an_col[(c="trt")],
         scale = "row",                # standardize genes
         cluster_cols = TRUE,
         cluster_rows = TRUE,
         show_rownames = TRUE,
         treeheight_col = 0,
         treeheight_row = 5,
         fontsize = 16,
         fontsize_row = 14,          # row label font size
         fontsize_col = 12)



dev.off()



# Your target gene to highlight
introg_ha412 <- fread("/mnt/data/DanaS/sam_broom_gwas/significants_with_snp_mapping/combined_by_ref/HA412_near_by_introgressions.tsv", 
                     header = TRUE)
introg_ha412 <- introg_ha412 %>% 
            filter(introg != "no_match")

introg_gene_id <- introg_ha412$gene_id

write.table(introg_gene_id,
            file = "/mnt/data/DanaS/rna/realign_ha412/deseq/introg_gene_id.txt",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)

write.table(marker_gene_ids,
            file = "/mnt/data/DanaS/rna/realign_ha412/deseq/marker_gene_ids.txt",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)

res_have$gene <- gsub("gene:", "", res_have$gene)
res_p05_gene_ids <- res_have$gene
res_p05_gene_ids <- res_p05_gene_ids %>% 
        as.data.frame()

write.table(res_p05_gene_ids,
            file = "/mnt/data/DanaS/rna/realign_ha412/deseq/res_p05_gene_idss.txt",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)

res_p05_gene_ids <- fread("/mnt/data/DanaS/rna/realign_ha412/deseq/group_euntreated_vs_treatment_sep14.txt", 
                     header = FALSE)
marker_gene_ids <- fread("/mnt/data/DanaS/rna/realign_ha412/deseq/marker_gene_ids.txt", 
                     header = FALSE)
introg_gene_id <- fread("/mnt/data/DanaS/rna/realign_ha412/deseq/introg_gene_id.txt", 
                     header = FALSE)

library(VennDiagram)
    deg = res_p05_gene_ids$V1
    gwas = marker_gene_ids$V1
    Introg = introg_gene_id$V1 
x = list(deg, gwas, Introg)
venn.diagram(
  x,
  filename = "gwas_deg_venn14.jpeg",
  category.names = c("DEG","GWAS","Introgressions"),
  imagetype = "png",
  #height = 2000,
  #width = 2000,
  #resolution = 300,
  fill = c("#ead919d5", "#56B4E9", "#009E73"),
  lwd = 1,
  output=TRUE  
)



jpeg("gwas_deg_venn14.jpeg", width = 500, height = 500, units = "px", quality = 100)
venn.diagram(merged_genes, filename = "gwas_deg_venn.jpeg")
display_venn <- function(x, ...){
  library(VennDiagram)
  grid.newpage()
  venn_object <- venn.diagram(x, filename = NULL, ...)
  grid.draw(venn_object)
}

install.packages("ggvenn")
library(ggvenn)


deg = res_p05_gene_ids$V1
gwas = marker_gene_ids$V1
Introg = introg_gene_id$V1 

merged_genes <- list(
  DEG = deg,
  GWAS = gwas,
  Introgressions = Introg
)

ggvenn(
  merged_genes, 
  fill_color = c("#E69F00", "#56B4E9", "#009E73"),
  stroke_size = 0.5, set_name_size = 4
)




display_venn(merged_genes,
  category.names = c("DEG","GWAS","Introgressions"),
  fill = c("#E69F00", "#56B4E9", "#009E73"), 
  )
dev.off()
target_gene <- deg_filtered %>% 
  filter(rownames(deg_filtered) %in% introg_gene_id)
target_gene <- marker_data[,c(4,10)] %>% as.data.frame() %>%
  filter(introg_gene_id %in% rownames(deg_filtered))
rownames(target_gene) <- target_gene$gene_id
rownames_map <- setNames(target_gene$product, target_gene$gene_id)
norm_mat1 <- norm_mat
rownames(norm_mat1) <- rownames_map[rownames(norm_mat1)]
# Replace rownames of df1
rownames(norm_mat) <- rownames_map[rownames(norm_mat)]
target_gene <- c("Ha412HOChr09g0413391","Ha412HOChr09g0413751")

# Create the pheatmap and save object
# Create a vector of row label colors
gene_labels <- rownames(norm_mat)
label_colors <- ifelse(gene_labels == target_gene, "red", "black")

# Plot
pheatmap_obj <- Heatmap(norm_mat, # Exclude target gene rows for separate handling
        name = "Expression",
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        show_row_names = TRUE,
        row_names_gp = gpar(col = label_colors, fontsize = 8),
        show_column_names = TRUE,
        column_names_gp = gpar(fontsize = 7),
        top_annotation = HeatmapAnnotation(df = ann_col[4:9,]))


jpeg(filename = "heatmap_col1.jpeg", width = 500, height = 500, units = "px", quality = 100)
print(pheatmap_obj)
dev.off()

jpeg("gwas_deg_col.jpeg", width = 500, height = 700, units = "px", quality = 100)
make_bold_names <- function(mat, rc_fun, rc_names) {
  bold_names <- rc_fun(mat)
  ids <- rc_names %>% match(rc_fun(mat))
  ids %>%
    walk(
      function(i)
        bold_names[i] <<-
        bquote(bold(.(rc_fun(mat)[i]))) %>%
        as.expression()
    )
  bold_names
}

wrapped_rownames <- str_wrap(rownames(norm_mat), width = 70, indent = 0)
jpeg(filename = "heatmap_col_names_sep14.jpeg", width=15.5, height=9, units = "in", res = 500)
pheatmap(norm_mat1[-c(1,6,17,25,29:31), ],
         annotation_col = an_col[(c="trt")],
         #annotation_row = target_gene[(c="product")],
         scale = "row",                # standardize genes
         cluster_cols = TRUE,
         cluster_rows = TRUE,
         treeheight_col = 0,
         treeheight_row = 0,
         #width=4.5, height=5,
           cellhight = 60,
           cellwidth = 40,
          fontsize = 15,
          fontsize_row = 16.2,          # row label font size
          fontsize_col = 15,
          #labels_row = wrapped_rownames
         #labels_row = make_bold_names(norm_mat, rownames, target_gene)
)

dev.off()


pheatmap(
  norm_mat,
  labels_row = make_bold_names(norm_mat, rownames, c("Ha412HOChr04g0177021", "Ha412HOChr17g0810781")),
  )

jpeg(filename = "heatmap_col_names.jpeg", width = 500, height = 500, units = "px", quality = 100)
print(pheatmap_obj)
dev.off()


library(pheatmap)
library(grid)
library(gridExtra)

# Target gene to highlight
target_gene <- "Ha412HOChr04g0177021"

# Run pheatmap in silent mode and store object
pheatmap_obj <- pheatmap(norm_mat,
                         annotation_col = ann_col[4:9, ],
                         scale = "row",
                         cluster_cols = TRUE,
                         cluster_rows = TRUE,
                         show_rownames = TRUE,
                         treeheight_col = 0,
                         treeheight_row = 5,
                         fontsize = 12,
                         fontsize_col = 7,
                         fontsize_row = 8,
                         silent = TRUE)

# Create vector of row name colors
row_colors <- ifelse(rownames(norm_mat) == target_gene, "red", "black")

# Locate the row_names grob
rowname_grob_index <- which(sapply(pheatmap_obj$gtable$grobs, function(x) x$name) == "row_names")

# Apply colors to row names
pheatmap_obj$gtable$grobs[[rowname_grob_index]]$gp$col <- row_colors

# Optional: create legend
rowname_legend <- legendGrob(c("Target gene", "Other genes"),
                             nrow = 2,
                             pch = 15,
                             gp = gpar(fontsize = 10, col = c("red", "black")))

# Draw and save
jpeg("heatmap_colored_rownames.jpeg", width = 1000, height = 1000, res = 150)
grid.draw(arrangeGrob(pheatmap_obj$gtable, rowname_legend, ncol = 2, widths = c(5, 1)))
dev.off()


# Note: Modify the design formula as appropriate.
# For example, if you wish to test for treatment effects accounting for time and accession,
# you might use: design = ~ accession + time + treatment
# or if you want to include interaction terms, adjust accordingly.





# Pre-filter low count genes (optional)
smallestGroupSize <- 3
keep <- rowSums(counts(ddstxi) >= 10) >= smallestGroupSize
ddstxi <- ddstxi[keep,]

ddstxi$treatment <- factor(ddstxi$treatment, levels = c("untreated","treated", "no_treatment"))
ddstxi$treatment <- relevel(ddstxi$treatment, ref = "treated")
# ------------------------------
# 4. Run Differential Expression Analysis with DESeq2
# ------------------------------
dds <- DESeq(ddstxi)
resultsNames(dds) # lists the coefficients

#a much simpler approach would give all the results tables that are desired. We will explain this approach first, because it is much simpler to perform. 
#If the comparisons of interest are, for example, the effect of a condition for different sets of samples, a simpler approach than adding interaction terms explicitly to the design formula is to perform the following steps:
#combine the factors of interest into a single factor with all combinations of the original factors
#change the design to include just this factor, e.g. ~ group
dds$group <- factor(paste0(dds$accession, dds$treatment))
#dds$group <- relevel(dds$group, ref = "euntreated")
design(dds) <- ~ group
dds <- DESeq(dds)
resultsNames(dds)
res <- results(dds, contrast=c("group", "euntreated","etreated"))
# Get results; for example, to compare treated vs untreated at time 5 you might use:
# Here, we assume that the "treatment" factor has levels "no_treatment", "untreated", and "treated"
# and that you want to compare "treated" to "untreated" at time 5.
res <- results(dds)
head(res)
res <- results(dds, name="treatment_untreated_vs_treated", test="Wald")
#res <- results(dds, contrast = c("b", "a"))

# Summarize the results
summary(res)

res_have <- res %>% as.data.frame() %>% 
            filter(!grepl("gene:Ha412HOChr00c", rownames(res))) %>% 
            na.omit(.)  %>% 
            filter(padj < 0.05) %>% 
            mutate(gene = rownames(.))

res_have <- res_have[c("gene", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")]

write.table(as.data.frame(res_have),
            file = "/mnt/data/DanaS/rna/realign_ha412/deseq/group_euntreated_vs_treatment_sep14.txt",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)


select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:200]
df <- as.data.frame(colData(dds)[,c("condition","accession")])
d <- pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)
jpeg(filename = "heatmap_col.jpeg", width = 500, height = 500, units = "px", quality = 100)
print(pheatmap_obj)
dev.off()
         
res05 <- results(dds, alpha=0.05)
summary(res05)

p <- plotCounts(dds, gene=which.max(res05$log2FoldChange), intgroup="condition", 
                returnData=TRUE)
ggplot(p, aes(x=condition, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,400))
png("plotCounts", width=7, height=7, units = "in", res = 300)
print(d)
dev.off()

head(res05)
resOrdered <- res[order(res$pvalue),]
# (Optional) Shrink log fold changes for visualization
resLFC <- lfcShrink(dds, coef = "condition_c_vs_b", type = "apeglm")

# multi factor design
dds$group <- factor(paste0(dds$name, dds$treatment))
design(dds) <- ~ group
design(dds) <- ~ accession + condition + accession:condition
dds <- DESeq(dds)
dds$treatment <- relevel(dds$treatment, ref = "treated")
dds <- DESeq(dds)
resultsNames(dds)
# condition treated vs untreated, tells you that the estimates 
#are of the logarithmic fold change log2(treated/untreated). 
eb_vs_eb <- results(dds, name="condition_b_vs_c", test="Wald")

res_p05_lfg1 <- res %>% as.data.frame() %>% 
            na.omit(.)  %>% 
            filter(padj < 0.05) %>% 
            filter(log2FoldChange >= 1 | log2FoldChange <= -0.01 ) %>%
            mutate(gene = rownames(.)) %>% 
            filter(!grepl('Ha412HOChr00c', gene))

res_p05_lfg1 <- res_p05_lfg1[c("gene", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")]

write.table(as.data.frame(res_p05_lfg1),
            file = "/mnt/data/DanaS/rna/realign_ha412/deseq2/emek_c_vs_b_res_p05_lfg1.txt",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)


rc_rb <- results(dds, contrast=c("group", "rc", "rb"))
ec_eb <- results(dds, contrast=c("group", "ec", "eb"))
sc_sb <- results(dds, contrast=c("group", "sc", "sb"))


res_rc_vs_ra <- results(dds,
  contrast = list(c("condition_c_vs_a","accessionr.conditionb"))
)



res_have <- ec_eb %>% as.data.frame() %>% 
            filter(!grepl("gene:Ha412HOChr00c", rownames(ec_eb))) %>% 
            na.omit(.)  %>% 
            filter(padj < 0.05) %>% 
            mutate(gene = rownames(.))

res_have <- res_have[c("gene", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")]

write.table(as.data.frame(res_have),
            file = "/mnt/data/DanaS/rna/realign_ha412/deseq/ec_eb_new.txt",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)


summary(rc_rb)
res_have <- rc_rb %>% as.data.frame() %>% 
            na.omit(.)  %>% 
            filter(padj < 0.05) %>% 
            filter(log2FoldChange >= 1) %>%
            mutate(gene = rownames(.))
res_have <- res_have[c("gene", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")]

write.table(as.data.frame(res_have),
            file = "/mnt/data/DanaS/rna/realign_ha412/deseq2/res_rc_vs_rb.txt",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)


# with interaction
design(dds) <- ~ genotype + condition + genotype:condition

res <- lfcShrink(dds, res=ec_eb, type = "apeglm")

volcano <- EnhancedVolcano(res,
    lab = rownames(res),
    x = 'log2FoldChange',
    y = 'pvalue')

res <- results(dds, contrast=c("treatment","treated","untreated"))
  volcano <- EnhancedVolcano(res,
    lab = rownames(res),
    x = 'log2FoldChange',
    y = 'padj',
    xlim = c(-15, 15),
    ylim = c(0, 16),
    title = 'treated vs. untreated',
    pCutoff = 0.05,
    FCcutoff = 2,
    pointSize = 1.0, labSize = 2.0,)
png("DGE_Volcano_treated_vs_untreated.png", width=7, height=7, units = "in", res = 300)
print(volcano)
dev.off()

# Define clearly the cutoffs
pCutoff <- 0.05
FCcutoff <- 1.0

p <- EnhancedVolcano(res,
    lab = rownames(res),
    x = 'log2FoldChange',
    y = 'padj',
    xlim = c(-15, 15),
    ylim = c(0, 16),
    pCutoff = pCutoff,       # use consistent naming
    FCcutoff = FCcutoff,     # use consistent naming
    pointSize = 1.0,
    labSize = 2.0,
    title = "treated vs. untreated",
    subtitle = "",
    caption = paste0('log2 FC cutoff: ', 
                     FCcutoff, '; p-value cutoff: ', 
                     pCutoff, '\nTotal = ', nrow(res), ' variables'),
    legendLabels = c('NS','Log2 FC','Adjusted p-value', 'Adjusted p-value & Log2 FC'), # Note: changed from 'legend'
    legendPosition = 'bottom',
    legendLabSize = 14,
    legendIconSize = 5.0)

png("DGE_Volcano_treated_vs_untreated.png", width=7, height=7, units = "in", res = 300)
print(p)  # <- you must explicitly print the plot inside png()
dev.off()

# ------------------------------
# 5. Export and Visualize Results
# ------------------------------
# Write the results to a CSV file
write.csv(as.data.frame(resLFC), file = "DESeq2_results.csv")

# Generate an MA plot
plotMA(resLFC, main = "MA-Plot", ylim = c(-5, 5))

# Optionally, inspect and plot other diagnostics such as PCA plots:
rld <- rlog(dds)
plotPCA(rld, intgroup = c("accession", "treatment"))

# End of script.

##############################################################################################################



# ------------------------------
# 3. Create DESeq2 Dataset Using tximport Data
# ------------------------------

# Create DESeq2 metadata: the row names must match the column names in txi.rsem.
# Here we assume that the order of samples in txi.rsem is the same as in metadata.
# identify biological replicates

coldata$accession <- as.factor(coldata$accession)
coldata$time <- factor(coldata$time)
coldata$condition <- factor(coldata$condition)
coldata$sample <- factor(coldata$sample)

rownames(coldata) <- colnames(data)

dds <- DESeqDataSetFromMatrix(countData = data,
                              colData = coldata,
                              design = ~ accession + condition + accession:condition)
# Note: Modify the design formula as appropriate.
# For example, if you wish to test for treatment effects accounting for time and accession,
# you might use: design = ~ accession + time + treatment
# or if you want to include interaction terms, adjust accordingly.

# Pre-filter low count genes (optional)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

smallestGroupSize <- min(table(coldata$condition))
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]

# ------------------------------
# 4. Run Differential Expression Analysis with DESeq2
# ------------------------------
dds <- DESeq(dds)
dds$accession <- factor(dds$accession, levels = c("e","r","s"))
dds$condition <- factor(dds$condition, levels = c("a","b","c"))

levels(dds$accession)
## [1] "e" "r" "s"

levels(dds$condition)
## [1] "a" "b" "c"

design(dds) <- ~ accession + condition + accession:condition
dds <- DESeq(dds)

resultsNames(dds)
res <- results(dds)
normalized_counts <- counts(dds, normalized = T)
normalized_counts <- rownames_to_column(as.data.frame(normalized_counts), var = "ensembl_id")
# filter out rows containing gene:Ha412HOChr00c and write to file
normalized_counts <- normalized_counts %>% 
  filter(!grepl("gene:Ha412HOChr00c", ensembl_id))
write.table(normalized_counts, sep = "\t", file = "normalized_counts.txt", row.names = F)
# in bash use the py script to get the annotations
# normalized_counts=normalized_counts.txt
# gff=/mnt/data/sunflower/HA412/Ha412HOv2.0-20181130.gff3
# python3 /mnt/data/DanaS/rna/get_norm_counts_annot.py $gff $normalized_counts
# this will output annotated_counts.txt
annotation <- read.table("/mnt/data/DanaS/rna/realign_ha412/output/annotated_counts.txt", 
                     sep = "\t", header = T, stringsAsFactors = F)

annotated_data <- right_join(annotation, normalized_counts, by = c("Gene_ID" = "ensembl_id"))
head(annotated_data)


install.packages("ashr")
library(ashr)
resAshT <- lfcShrink(dds, coef="condition_b_vs_a", type="ashr", lfcThreshold=1)
resAshT
sum(resAshT$padj < 0.1, na.rm=TRUE)

jpeg(filename = "MA.jpeg", width = 500, height = 500, units = "px", quality = 100)
plotMA(res, ylim=c(-3,3))
dev.off()
# distance plot - Euclidean distance between the expression values for each individual sample
vsd <- vst(dds, blind = TRUE)
plotDists = function (vsd.obj) {
  sampleDists <- dist(t(assay(vsd.obj)))
  sampleDistMatrix <- as.matrix( sampleDists )
  rownames(sampleDistMatrix) <- paste(vsd.obj$sample )
  colors <- colorRampPalette( rev(RColorBrewer::brewer.pal(9, "Blues")) )(255)
  pheatmap::pheatmap(sampleDistMatrix,
           clustering_distance_rows = sampleDists,
           clustering_distance_cols = sampleDists,
           col = colors)
}
jpeg(filename = "dist.jpeg", width = 500, height = 500, units = "px", quality = 100)
plotDists(vsd)
dev.off()

sampleDists <- dist(t(assay(vsd)))

sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- NULL
colnames(sampleDistMatrix) <- paste(vsd$sample)
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
jpeg(filename = "dist1.jpeg", width = 500, height = 500, units = "px", quality = 100)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dev.off()


rld <- rlog(dds)
sampleDists <- dist(t(assay(vsd)))

sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- NULL
colnames(sampleDistMatrix) <- paste(vsd$sample)
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
jpeg(filename = "dist1.jpeg", width = 500, height = 500, units = "px", quality = 100)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dev.off()

pcaData <- plotPCA(vsd, intgroup=c("accession", "condition"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
pca_groups <- ggplot(pcaData, aes(PC1, PC2, color=condition, shape=accession)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()
ggsave("pca_groups.jpeg", plot = pca_groups, width = 7, height = 5)
# Get results; for example, to compare treated vs untreated at time 5 you might use:
# Here, we assume that the "treatment" factor has levels "no_treatment", "untreated", and "treated"
# and that you want to compare "treated" to "untreated" at time 5.

res <- results(dds, contrast = c("condition", "c", "b"))

# Summarize the results
summary(res)
head(res)
# (Optional) Shrink log fold changes for visualization
resLFC <- lfcShrink(dds, coef = "treatment_treated_vs_untreated", type = "apeglm")
groups <- levels(dds$group)
groups
# [1] "ea" "eb" "ec" "ra" "rb" "rc" "sa" "sb" "sc"

# Create an empty list to hold results
res_list <- list()

# Double loop over all unique pairs
for (i in seq_along(groups)) {
  for (j in seq_along(groups)) {
    if (i < j) {
      grp1 <- groups[i]
      grp2 <- groups[j]
      
      # Create a name like "ea_vs_eb"
      contrast_name <- paste0(grp1, "_vs_", grp2)
      
      # Run the contrast: "group" is the factor in the design,
      # grp1 is numerator, grp2 is denominator
      res_list[[contrast_name]] <- results(dds, 
        contrast = c("group", grp1, grp2)
      )
    }
  }
}

res <- results(dds, contrast = c("sample", "sb", "sc"))
res_samp <- results(dds, name="condition_c_vs_a", test="Wald")
# ------------------------------
# 5. Export and Visualize Results
# ------------------------------
# Write the results to a file
res_rc_vs_ra <- results(dds,
  contrast = list(c("condition_c_vs_a","accessionr.conditionb"))
)



res_have <- res_rc_vs_rb %>% as.data.frame() %>% 
            filter(!grepl("gene:Ha412HOChr00c", rownames(res_rc_vs_rb))) %>% 
            na.omit(.)  %>% 
            filter(padj < 0.05) %>% 
            mutate(gene = rownames(.))

res_have <- res_have[c("gene", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")]

write.table(as.data.frame(res_have),
            file = "/mnt/data/DanaS/rna/realign_ha412/deseq/res_rc_vs_rb.txt",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)
# Generate an MA plot
plotMA(res, main = "MA-Plot", ylim = c(-5, 5))

# Optionally, inspect and plot other diagnostics such as PCA plots:
rld <- rlog(dds)
pca <- plotPCA(rld, intgroup = "group")
jpeg(filename = "pca_rld.jpeg", width = 500, height = 500, units = "px", quality = 100)
plotPCA(rld, intgroup = "group")
dev.off()


count_table_path <- "/mnt/data/DanaS/rna/realign_ha412/output/raw_counts.txt"
data <- read.table(count_table_path, 
                     header = TRUE, sep = "\t", 
                     row.names = "ensembl_id")
data <- data[,sort(colnames(data))]
head(data)
colSums(data)


# identify biological replicates
condition <- c(rep("ea", 3), rep("eb", 3), rep("ec", 3),
                rep("ra", 3), rep("rb", 3), rep("rc", 3),
                rep("sa", 3), rep("sb", 3), rep("sc", 3))

# assign replicates to each sample name to construct colData
my_colData <- as.data.frame(coldata)
rownames(colData) <- colnames(data)
my_colData
design(ddsMF) <- formula(~ accession + time + treatment)
ddsMF <- DESeq(ddsMF)
dds <- DESeqDataSetFromMatrix(countData = data,
                              colData = metadata,
                              design = ~ accession + time + treatment)



select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:100]
df <- as.data.frame(colData(dds)[,c("condition","accession")])
jpeg(filename = "dist2.jpeg", width = 500, height = 500, units = "px", quality = 100)
pheatmap(assay(vsd)[select,], cluster_rows=TRUE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)
dev.off()





d <- plotCounts(dds, gene="gene:Ha412HOChr07g0323951", intgroup="treatment", returnData=TRUE)

pca_gene <- ggplot(d, aes(x=treatment, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0), aes (color = treatment)) + 
  scale_y_log10(breaks=c(25,100,400))
ggsave("pca_groups_Ha412HOChr07g0323951.jpeg", plot = pca_gene, width = 7, height = 5)


res <- results(dds)
min_p <- plotCounts(dds, gene=which.min(res$padj), intgroup="sample", 
                returnData=TRUE)
minp <- ggplot(min_p, aes(x=sample, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,400))
ggsave("pca_groups_min_p.jpeg", plot = minp, width = 7, height = 5)



vsd <- vst(dds, blind = TRUE)
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds)[,c("condition","accession")])
pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)


# multi factor design

ddsMF <- dds
design(ddsMF) <- formula(~ accession + treatment)
ddsMF <- DESeq(ddsMF)
resMF <- results(ddsMF)

min_p <- plotCounts(ddsMF, gene=which.min(res$padj), intgroup="sample", 
                returnData=TRUE)
minp <- ggplot(min_p, aes(x=sample, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,400))
ggsave("pca_groups_min_p.jpeg", plot = minp, width = 7, height = 5)


# time series
# Combine time and treatment into one factor
coldata$group <- factor(paste0("t", coldata$time, "_", coldata$treatment))
# Optionally, set the reference level (here, assuming t0_untreated0 is your baseline)
coldata$group <- relevel(coldata$group, ref = "t0_untreated0")

ddsTS <- DESeqDataSetFromMatrix(countData = data,
                                colData = coldata,
                                design = ~ accession + group + accession:group)

ddsT <- DESeq(ddsTS)
ddsTC <- DESeq(ddsTS, test="LRT", reduced = ~ accession + group)
resTC <- results(ddsTC)

fiss <- plotCounts(ddsTC, which.min(resTC$padj), 
                   intgroup = c("group","accession"), returnData = TRUE)
fiss$group <- as.numeric(as.character(fiss$group))
f <- ggplot(fiss,
  aes(x = group, y = count, color = accession, group = accession)) + 
  geom_point() + stat_summary(fun.y=mean, geom="line") +
  scale_y_log10()
ggsave("pca_groups_min_p.jpeg", plot = minp, width = 7, height = 5)

resultsNames(ddsTC)
res5 <- results(ddsTC, name="accessionr.groupt5_treated5", test="Wald")
ordered <- res5[order(res5$padj),]
head(ordered)

res_samp <- results(dds, name="condition_c_vs_a", test="Wald")
ordered_samp <- res_samp[order(res_samp$padj),]
head(ordered_samp)
summary(res)
normalized_counts <- counts(dds, normalized = T)
normalized_counts <- rownames_to_column(as.data.frame(normalized_counts), var = "ensembl_id")
write.table(normalized_counts, sep = "\t", file = "normalized_counts.txt", row.names = F)
DE_gene_heatmap <- function(res, padj_cutoff = 0.05, ngenes = 20) {
  # generate the color palette
  brewer_palette <- "RdBu"
  ramp <- colorRampPalette(RColorBrewer::brewer.pal(11, brewer_palette))
  mr <- ramp(256)[256:1]
  # obtain the significant genes and order by log2FoldChange
  significant_genes <- res %>% filter(padj < padj_cutoff) %>% arrange (desc(log2FoldChange)) %>% head (ngenes)
  heatmap_values <- as.matrix(significant_genes[,-c(1:8)])
  rownames(heatmap_values) <- significant_genes$Gene.name
  # plot the heatmap using pheatmap
  pheatmap::pheatmap(heatmap_values, color = mr, scale = "row", fontsize_col = 10, fontsize_row = 200/ngenes, fontsize = 5, border_color = NA)
}

jpeg(filename = "deg_try.jpeg", width = 500, height = 500, units = "px", quality = 100)
DE_gene_heatmap(res)
dev.off()


annotation <- read.table("/mnt/data/DanaS/rna/realign_ha412/output/annotated_counts.txt", 
                     sep = "\t", header = T, stringsAsFactors = F)

head(annotation)

annotated_data <- right_join(annotation, normalized_counts, by = c("Gene_ID" = "ensembl_id"))
head(annotated_data)

variable_gene_heatmap <- function (vsd.obj, num_genes = 50, annotation, title = "") {
  brewer_palette <- "RdBu"
  # Ramp the color in order to get the scale.
  ramp <- colorRampPalette( RColorBrewer::brewer.pal(11, brewer_palette))
  mr <- ramp(256)[256:1]
  # get the stabilized counts from the vsd object
  stabilized_counts <- assay(vsd.obj)
  # calculate the variances by row(gene) to find out which genes are the most variable across the samples.
  row_variances <- rowVars(as.matrix(stabilized_counts), useNames = TRUE)
  # get the top most variable genes
  top_variable_genes <- stabilized_counts[order(row_variances, decreasing=T)[1:num_genes],]
  # subtract out the means from each row, leaving the variances for each gene
  top_variable_genes <- top_variable_genes - rowMeans(top_variable_genes, na.rm=T)
  # replace the ensembl ids with the gene names
  gene_names <- annotation$Gene.name[match(rownames(top_variable_genes), annotation$Gene_ID)]
  rownames(top_variable_genes) <- gene_names
  # reconstruct colData without sizeFactors for heatmap labeling
  coldata <- as.data.frame(vsd.obj@colData)
  coldata <- as.data.frame(vsd.obj@colData)
  coldata$sizeFactor <- NULL  # remove sizeFactor if present
  
  # Subset coldata to include only the 'treatment' column.
  if ("condition" %in% colnames(coldata)) {
    coldata <- coldata[, "condition", drop = FALSE]
  } else {
    warning("Column 'treatment' not found in colData; using all available annotations.")
  }
  # draw heatmap using pheatmap
  pheatmap::pheatmap(top_variable_genes, color = mr, annotation_col = coldata, fontsize_col = 8, fontsize_row = 250/num_genes, border_color = NA, main = title)
}
jpeg(filename = "deg.jpeg", width = 500, height = 500, units = "px", quality = 100)

variable_gene_heatmap(vsd, num_genes = 500, annotation = annotation)
dev.off()



generate_DE_results <- function (dds, comparisons, padjcutoff = 0.001, log2cutoff = 0.5, cpmcutoff = 2) {
  # generate average counts per million metric from raw count data 
  raw_counts <- counts(dds, normalized = F)
  cpms <- enframe(rowMeans(edgeR::cpm(raw_counts)))
  colnames(cpms) <- c("ensembl_id", "avg_cpm")
  
  # extract DESeq results between the comparisons indicated
  res <- results(dds, contrast = c("treatment", comparisons[1], comparisons[2]))[,-c(3,4)]
  
  # annotate the data with gene name and average counts per million value
  res <- as_tibble(res, rownames = "ensembl_id")

  res <- left_join(res, annotated_data, by = c("ensembl_id" = "Gene_ID"))
  # append the average cpm value to the results data
  res <- left_join(res, cpms, by = c("ensembl_id" = "ensembl_id"))
  
  # combine normalized counts with entire DE list
  normalized_counts <- round(counts(dds, normalized = TRUE),3)
  pattern <- str_c(comparisons[1], "|", comparisons[2])
  combined_data <- as_tibble(cbind(res, normalized_counts[,grep(pattern, colnames(normalized_counts))] ))
  combined_data <- combined_data[order(combined_data$log2FoldChange, decreasing = T),]
  
  # make ordered rank file for GSEA, selecting only protein coding genes
  res_prot_ranked <- res[order(res$log2FoldChange, decreasing = T),c("Gene.name", "log2FoldChange")]
  res_prot_ranked <- na.omit(res_prot_ranked)
  res_prot_ranked$Gene.name <- str_to_upper(res_prot_ranked$Gene.name)
  
  # generate sorted lists with the indicated cutoff values
  res <- res[order(res$log2FoldChange, decreasing=TRUE ),]
  de_genes_padj <- res[which(res$padj < padjcutoff),]
  de_genes_log2f <- res[which(abs(res$log2FoldChange) > log2cutoff & res$padj < padjcutoff),]
  de_genes_cpm <- res[which(res$avg_cpm > cpmcutoff & res$padj < padjcutoff),]
  
  # write output to files
  write.csv (de_genes_padj, file = paste0(comparisons[1], "_vs_", comparisons[2], "_padj_cutoff.csv"), row.names =F)
  write.csv (de_genes_log2f, file = paste0(comparisons[1], "_vs_", comparisons[2], "_log2f_cutoff.csv"), row.names =F)
  write.csv (de_genes_cpm, file = paste0(comparisons[1], "_vs_", comparisons[2], "_cpm_cutoff.csv"), row.names =F)
  write.csv (combined_data, file = paste0(comparisons[1], "_vs_", comparisons[2], "_allgenes.csv"), row.names =F)
  write.table (res_prot_ranked, file = paste0(comparisons[1], "_vs_", comparisons[2], "_rank.rnk"), sep = "\t", row.names = F, quote = F)
  
  writeLines( paste0("For the comparison: ", comparisons[1], "_vs_", comparisons[2], ", out of ", nrow(combined_data), " genes, there were: \n", 
               nrow(de_genes_padj), " genes below padj ", padjcutoff, "\n",
               nrow(de_genes_log2f), " genes below padj ", padjcutoff, " and above a log2FoldChange of ", log2cutoff, "\n",
               nrow(de_genes_cpm), " genes below padj ", padjcutoff, " and above an avg cpm of ", cpmcutoff, "\n",
               "Gene lists ordered by log2fchange with the cutoffs above have been generated.") )
  gene_count <- tibble (cutoff_parameter = c("padj", "log2fc", "avg_cpm" ), 
                        cutoff_value = c(padjcutoff, log2cutoff, cpmcutoff), 
                        signif_genes = c(nrow(de_genes_padj), nrow(de_genes_log2f), nrow(de_genes_cpm)))
  invisible(gene_count)
}




library(dplyr)
library(VennDiagram)
#install.packages("VennDiagram")
get_sig_genes <- function(file, padj_cutoff = 0.05) {
  df <- read.table(file, header=TRUE, sep="\t", stringsAsFactors=FALSE)
  # Adjust column names if needed (e.g. `df$padj` or `df$FDR` if you used edgeR, etc.)
  df %>%
    filter(padj < padj_cutoff) %>%
    pull(gene) %>%
    unique()
}

# E
genes_ea_eb <- get_sig_genes("/mnt/data/DanaS/rna/realign_ha412/deseq/DESeq2_results_ea_eb.txt")
genes_ea_ec <- get_sig_genes("/mnt/data/DanaS/rna/realign_ha412/deseq/DESeq2_results_ea_ec.txt")
genes_eb_ec <- get_sig_genes("/mnt/data/DanaS/rna/realign_ha412/deseq/DESeq2_results_eb_ec.txt")

# R
genes_ra_rb <- get_sig_genes("/mnt/data/DanaS/rna/realign_ha412/deseq/DESeq2_results_ra_rb.txt")
genes_ra_rc <- get_sig_genes("/mnt/data/DanaS/rna/realign_ha412/deseq/DESeq2_results_ra_rc.txt")
genes_rb_rc <- get_sig_genes("/mnt/data/DanaS/rna/realign_ha412/deseq/DESeq2_results_rb_rc.txt")

# S
genes_sa_sb <- get_sig_genes("/mnt/data/DanaS/rna/realign_ha412/deseq/DESeq2_results_sa_sb.txt")
genes_sa_sc <- get_sig_genes("/mnt/data/DanaS/rna/realign_ha412/deseq/DESeq2_results_sa_sc.txt")
genes_sb_sc <- get_sig_genes("/mnt/data/DanaS/rna/realign_ha412/deseq/DESeq2_results_sb_sc.txt")



set_E_ac <- genes_ea_ec
set_R_ac <- genes_ra_rc
set_S_ac <- genes_sa_sc

venn.plot1 <- venn.diagram(
  x = list(
    E_ac = set_E_ac,
    R_ac = set_R_ac,
    S_ac = set_S_ac
  ),
  filename = "01_Venn_E_R_S_a_vs_c.png",
  main = "Venn 1: E vs R vs S (pre vs post w/ O.cumana)",
  imagetype = "png",
  height = 2000,
  width = 2000,
  resolution = 300,
  fill = c("#FF9999","#99CC99","#9999CC")
)

set_E_ab <- genes_ea_eb
set_R_ab <- genes_ra_rb
set_S_ab <- genes_sa_sb
genes_ea_eb <- rownames(res_ea_eb)
genes_ra_rb <- rownames(res_ra_rb)
genes_sa_sb <- rownames(res_sa_sb)

venn.diagram(
  x = list(
    E_ab = genes_ea_eb,
    R_ab = genes_ra_rb,
    S_ab = genes_sa_sb
  ),
  filename = "02_Venn_E_R_S_a_vs_b.png",
  main = "Venn 2: E, R, S (pre vs post w/o O. cumana)",
  imagetype = "png",
  height = 2000,
  width = 2000,
  resolution = 300,
  fill = c("#FF9999","#99CC99","#9999CC")
)





filter_and_write <- function(res, out_file, pattern="gene:Ha412HOChr00c") {
  df_filt <- res %>%
    as.data.frame() %>%
    # remove genes containing pattern in row name
    filter(!grepl(pattern, rownames(res))) %>%
    # remove rows with NA
    na.omit() %>%
    # keep only padj < 0.05
    filter(padj < 0.05) %>%
    # add gene column
    mutate(gene = rownames(.)) %>%
    # re-arrange columns
    select(gene, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj)
  
  # Write out
  write.table(
    df_filt,
    file = out_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
}

###
# GENOTYPE E (reference)
###

# 1) (e,b) vs (e,a) --> ea_eb
res_ea_eb <- results(dds, name="condition_b_vs_a")
filter_and_write(res_ea_eb, "results_ea_eb.txt")

# 2) (e,c) vs (e,a) --> ea_ec
res_ea_ec <- results(dds, name="condition_c_vs_a")
filter_and_write(res_ea_ec, "results_ea_ec.txt")

# 3) (e,c) vs (e,b) --> eb_ec
#   = (condition_c_vs_a) - (condition_b_vs_a)
res_eb_ec <- results(dds,
  contrast = list("condition_c_vs_a", "condition_b_vs_a"),
  listValues = c(1, -1)
)
filter_and_write(res_eb_ec, "results_eb_ec.txt")


###
# GENOTYPE R
###

# 4) (r,b) vs (r,a) --> ra_rb
#   = condition_b_vs_a + accessionr.conditionb
res_ra_rb <- results(dds,
  contrast = list(c("condition_b_vs_a", "accessionr.conditionb"))
)
filter_and_write(res_ra_rb, "results_ra_rb.txt")

# 5) (r,c) vs (r,a) --> ra_rc
#   = condition_c_vs_a + accessionr.conditionc
res_ra_rc <- results(dds,
  contrast = list(c("condition_c_vs_a", "accessionr.conditionc"))
)
filter_and_write(res_ra_rc, "results_ra_rc.txt")

# 6) (r,c) vs (r,b) --> rb_rc
#   = [c+a(r:c)] - [b+a(r:b)] in main + interaction terms
res_rb_rc <- results(dds,
  contrast = list(
    c("condition_c_vs_a", "accessionr.conditionc"),
    c("condition_b_vs_a", "accessionr.conditionb")
  ),
  listValues = c(1, -1)
)
filter_and_write(res_rb_rc, "results_rb_rc.txt")


###
# GENOTYPE S
###

# 7) (s,b) vs (s,a) --> sa_sb
#   = condition_b_vs_a + accessions.conditionb
res_sa_sb <- results(dds,
  contrast = list(c("condition_b_vs_a", "accessions.conditionb"))
)
filter_and_write(res_sa_sb, "results_sa_sb.txt")

# 8) (s,c) vs (s,a) --> sa_sc
#   = condition_c_vs_a + accessions.conditionc
res_sa_sc <- results(dds,
  contrast = list(c("condition_c_vs_a", "accessions.conditionc"))
)
filter_and_write(res_sa_sc, "results_sa_sc.txt")

# 9) (s,c) vs (s,b) --> sb_sc
res_sb_sc <- results(dds,
  contrast = list(
    c("condition_c_vs_a", "accessions.conditionc"),
    c("condition_b_vs_a", "accessions.conditionb")
  ),
  listValues = c(1, -1)
)
filter_and_write(res_sb_sc, "results_sb_sc.txt")


library(data.table)
res_p05_gene_ids <- fread("/mnt/data/DanaS/rna/realign_ha412/deseq/group_euntreated_vs_treatment_sep14.txt")

str(an_col)
