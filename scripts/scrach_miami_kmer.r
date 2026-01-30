
# Rscript path/to/sript gwas_res1 gwas_res2 simplem_value traits_yavur_gadot

library(miamiplot)
library(dplyr)
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)
gwas_res1=args[1]
simplem_up=args[2] 
bonfer_up=args[3] 
gwas_res2=args[4]
simplem_down=args[5] 
bonfer_down=args[6] 
Tag1=args[7]
Tag2=args[8]  
highlight1=args[9]
highlight2=args[10]
k_thresh_up=args[11]  
k_thresh_down=args[12]        # output perfix => trait name
#kmers=args[3]
print(args)

print("Begin Miami plot....")
# read the data
Raw1 <- read.table(gwas_res1,header=TRUE)
Raw2 <- read.table(gwas_res2,header=TRUE)
highlight_data1 <- read.table(highlight1,header=TRUE)
highlight_data2 <- read.table(highlight2,header=TRUE)
# Add a new column to distinguish the datasets
gwas1 <- Raw1 %>% mutate(study = "Yavur")
gwas2 <- Raw2 %>% mutate(study = "Gadot")
k1 <- highlight_data1 %>% mutate(study = "kYavur")
colnames(k1) <- c("chr","ps","p_wald","study")
k2 <- highlight_data2 %>% mutate(study = "kGadot")
colnames(k2) <- c("chr","ps","p_wald","study")

# Combine the datasets
colnames(gwas1)
colnames(gwas2)
colnames(k1)
colnames(k2)
#cols.num <- c("chr","ps","p_wald")
# gwas1[cols.num] <- sapply(gwas1[cols.num],as.numeric)
# gwas2[cols.num] <- sapply(gwas2[cols.num],as.numeric)
# k1[cols.num] <- sapply(k1[cols.num],as.numeric)
# k2[cols.num] <- sapply(k2[cols.num],as.numeric)
gwas1$chr <- as.numeric(gwas1$chr)
gwas2$chr <- as.numeric(gwas2$chr)
k1$chr <- as.numeric(k1$chr)
k2$chr <- as.numeric(k2$chr)
gwas_results <- bind_rows(gwas1, gwas2, k1, k2)
#gwas_results <- bind_rows(gwas1, gwas2)
head(gwas_results)
colnames(gwas_results) <- c("chr","pos","pval","study")
gwas_results$chr <- as.numeric(gwas_results$chr)
head(gwas_results)
sim_up = as.numeric(simplem_up)
bonf_up = as.numeric(bonfer_up)
bon_up = 0.05/bonf_up
kmer_up = as.numeric(k_thresh_up)
sim_down = as.numeric(simplem_down)
bonf_down = as.numeric(bonfer_down)
bon_down = 0.05/bonf_down
kmer_down = as.numeric(k_thresh_down)
#Say you wanted to plot the results from study "A" in the upper plot and
#study "B" in the lower plot.
    
plot_df <- gwas_results %>% 
  group_by(chr) %>% 
  # Compute chromosome size
  summarise(chrlength = max(pos)) %>%  
  # Calculate cumulative position of each chromosome
  mutate(cumulativechrlength = cumsum(as.numeric(chrlength))-chrlength) %>% 
  select(-chrlength) %>%
  # Temporarily add the cumulative length of each chromosome to the initial 
  # dataset 
  left_join(gwas_results, ., by=c("chr"="chr")) %>%
  # Sort by chr then position 
  arrange(chr, pos) %>%
  # Add the position to the cumulative chromosome length to get the position of 
  # this probe relative to all other probes
  mutate(rel_pos = pos + cumulativechrlength) %>%
  # Calculate the logged p-value too
  mutate(logged_p = -log10(pval)) %>%
  select(-cumulativechrlength)


axis_df <- plot_df %>% 
  group_by(chr) %>% 
  summarize(chr_center=(max(rel_pos) + min(rel_pos)) / 2)


max_y_limit <- max(plot_df$logged_p)


upper_plot <- ggplot() + 
  geom_point(data = plot_df[which(plot_df$study == "Yavur"),], 
             aes(x = rel_pos, y = logged_p, color = as.factor(chr)), 
             size = 2, alpha=0.6) +
  geom_point(data = plot_df[which(plot_df$study == "kYavur"),], 
             aes(x = rel_pos, y = logged_p), 
             size = 2.3,
             color = "#de2b0f") +  # Change size as needed for emphasis
  scale_color_manual(values = rep(c("cadetblue", "bisque4"), nrow(axis_df))) +
  scale_x_continuous(labels = axis_df$chr, 
                     breaks = axis_df$chr_center, 
                     expand = expansion(mult = 0.01)) +
  scale_y_continuous(limits = c(0, 12), 
                     expand = expansion(mult = c(0.02, 0))) + 
  geom_hline(yintercept = -log10(bon_up), color = "#216b07e9", linetype = "dashed", 
             linewidth=1.2) +
  geom_hline(yintercept = -log10(sim_up), color = "blue", linetype = "solid", 
             linewidth=0.7) +
  # geom_hline(yintercept = -log10(kmer_up), color = "#55ff00", linetype = "twodash", 
  #            size = 0.3) +
  labs(x = "", y = 'Yavor') + #quote(atop('Yavor', '-log'[10]*'(p)'))) + 
  theme_classic() +
  theme(legend.position = "none", 
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 35),
        axis.text.y = element_text(size = 35),
        axis.title.y = element_text(size = 41, margin = unit(c(0,1.5,0,0), "cm")),
        plot.margin = unit(c(1,0.1,0,0.1), "cm")) #+ 

lower_plot <- ggplot() + 
  geom_point(data = plot_df[which(plot_df$study == "Gadot"),], 
             aes(x = rel_pos, y = logged_p, color = as.factor(chr)), 
             size = 2, alpha=0.6) +
  geom_point(data = plot_df[which(plot_df$study == "kGadot"),], 
             aes(x = rel_pos, y = logged_p), 
             size = 2.3,
             color = "#de2b0f") +
  scale_color_manual(values = rep(c("cadetblue", "bisque4"), nrow(axis_df))) +
  scale_x_continuous(breaks = axis_df$chr_center, 
                     position = "top",
                     expand = expansion(mult = 0.01)) +
  scale_y_reverse(limits = c(12, 0), 
                     expand = expansion(mult = c(0, 0.02))) + 
  geom_hline(yintercept = -log10(bon_down), color = "#216b07e9", linetype = "dashed", 
             linewidth=1) +
  geom_hline(yintercept = -log10(sim_down), color = "blue", linetype = "solid", 
             linewidth=0.7) +
  # geom_hline(yintercept = -log10(kmer_down), color = "#55ff00", linetype = "twodash", 
  #            size = 0.3) +
  labs(x = "", y = 'Gadot') + #bquote(atop('Gadot', '-log'[10]*'(p)'))) + 
  theme_classic() +
  theme(legend.position = "none", 
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 35),
        axis.title.y = element_text(size = 41, margin = unit(c(0,1.5,0,0), "cm")),
        axis.title.x = element_blank(),
        plot.margin = unit(c(0,0.1,1,0.1), "cm")) #+ 




print("ALMOST DONE")
p <- gridExtra::grid.arrange(upper_plot, lower_plot, nrow = 2)
ggsave(p, filename = paste0(Tag1," ",Tag2,".png"), device = "png", 
width = 26, height = 9, dpi = 300)


print("Miami plot - DONE !")

#scale_color_manual(values = rep(c("darkkhaki","darkmagenta"), nrow(axis_df))) +