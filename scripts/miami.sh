# miami plot
# gemma out look like this:

# chr     rs      ps      n_miss  allele1 allele0 af      beta    se      logl_H1 l_remle l_mle   p_wald  p_lrt   p_score
# 01      Ha412HOChr01_210696_C_T 210696  4       T       C       0.059   2.754329e-02    1.388374e-01    -9.560921e+01 6.537449e-02     7.073131e-02    8.430536e-01    8.476450e-01    8.498864e-01
# 01      Ha412HOChr01_210723_C_T 210723  1       T       C       0.076   -7.521499e-02   1.238212e-01    -9.544745e+01 9.693202e-02     1.018898e-01    5.446129e-01    5.482662e-01    5.609404e-01

# data should look like (first and last 2 rows)

# rsid chr pos        beta         se     tstat      pval study
# 1  rs1   1   1 -0.00445690 0.00976958 -0.456202 0.6482490     A
# 2  rs2   1   2 -0.00274477 0.01417810 -0.193592 0.8464970     A
# 59995 rs29995  22 121  0.00374490 0.00753038  0.4973060 0.6189770     B
# 59996 rs29996  22 122  0.00088521 0.00963444  0.0918797 0.9267940     B


conda activate r_4.1
# get the wanted cols in this order

up=/mnt/data/DanaS/zero_one/filtered_by_NEC/z_25/out_gemma/preplot/n_NEC_mean_zero_one.assoc.txt
down=/mnt/data/DanaS/may_2024/gadot/by_n2r/mean/NEC/thresh_03/out_gemma/n_NEC_n2r_max_zero_one.assoc.txt
simple_up=$(awk 'NR==2 {print $2}' /mnt/data/DanaS/zero_one/filtered_by_NEC/z_25/*_simplem.txt)
bonfer_up=$(awk 'NR==2 {print $1}' /mnt/data/DanaS/zero_one/filtered_by_NEC/z_25/*_simplem.txt)
simple_down=$(awk 'NR==2 {print $2}' /mnt/data/DanaS/may_2024/gadot/by_n2r/mean/NEC/thresh_025/*_simplem.txt)
bonfer_down=$(awk 'NR==2 {print $1}' /mnt/data/DanaS/may_2024/gadot/by_n2r/mean/NEC/thresh_025/*_simplem.txt)


awk '{print $2, $1, $3, $8, $9, $10, $13}' $up > yavur_miami.txt
awk '{print $2, $1, $3, $8, $9, $10, $13}' $down > gadot_miami.txt

# Rscript path/to/sript gwas_res1 gwas_res2 thresh.p_value file_name

Rscript /mnt/data/DanaS/r_scripts/miami_plot.r yavur_miami.txt gadot_miami.txt $simple n_hel_zeroone_thresh05_yavur_gadot

# from scrash
Rscript /mnt/data/DanaS/r_scripts/scrach_miami.r yavur_miami.txt $simple_up $bonfer_up gadot_miami.txt $simple_down $bonfer_down n_NEC_n2r_Yavor_mean_thresh05_zero_one n_NEC_n2r_Gadot_max_thresh025_bw



cd /mnt/data/DanaS/may_2024/miami
up=/mnt/data/DanaS/may_2024/yavor/by_n2r/mean/NEC/thresh_03/out_gemma/n_NEC_n2r_max_zero_one.assoc.txt
simple_up=$(awk 'NR==2 {print $2}' /mnt/data/DanaS/may_2024/yavor/by_n2r/mean/NEC/thresh_03/yavor_NEC_mean_thresh03_simplem.txt)
bonfer_up=$(awk 'NR==2 {print $1}' /mnt/data/DanaS/may_2024/yavor/by_n2r/mean/NEC/thresh_03/yavor_NEC_mean_thresh03_simplem.txt)
kmer_up=/mnt/data/DanaS/sam_broom_gwas/gwas/yavor_nec_thresh_03/kmer/get_kmers/ha412/kmers_manhattan_plot_data_pvast.txt
down=/mnt/data/DanaS/may_2024/gadot/by_n2r/mean/NEC/thresh_03/out_gemma/n_NEC_n2r_max_zero_one.assoc.txt
simple_down=$(awk 'NR==2 {print $2}' /mnt/data/DanaS/may_2024/gadot/by_n2r/mean/NEC/thresh_03/gadot_NEC_mean_thresh03_simplem.txt)
bonfer_down=$(awk 'NR==2 {print $1}' /mnt/data/DanaS/may_2024/gadot/by_n2r/mean/NEC/thresh_03/gadot_NEC_mean_thresh03_simplem.txt)
kmer_down=/mnt/data/DanaS/sam_broom_gwas/gwas/gadot_nec_thresh_03/kmer/get_kmers/ha412/kmers_manhattan_plot_data.txt


awk '{print $1, $3, $13}' $up > snpup.txt
awk '{print $1, $3, $13}' $down > snpdown.txt
awk '{print $1, $2, $6}' $kmer_up > kup.txt 
awk '{print $1, $2, $3}' $kmer_down > kdown.txt 

head -n2 snpup.txt
sed -i 's/^Ha412HOChr//g' snpdown.txt
sed -i 's/^Ha412HOChr//g' snpup.txt 
Rscript /mnt/data/DanaS/r_scripts/scrach_miami.r \
/mnt/data/DanaS/may_2024/miami/snpup.txt $simple_up $bonfer_up \
/mnt/data/DanaS/may_2024/miami/snpdown.txt $simple_down $bonfer_down \
"miami" n_NEC03_n2r_zero_one



up=snp_file_upper.assoc.txt
simple_up=simplem_val_upper
bonfer_up=bonferoni_val_upper
kmer_up=kmers_manhattan_plot_data.txt
down=snp_file_bottom.assoc.txt
simple_down=simplem_val_bottom
bonfer_down=bonferoni_val_bottom
kmer_down=kmers_manhattan_plot_data.txt

awk '{print $1, $3, $13}' $up > snpup.txt
awk '{print $1, $3, $13}' $down > snpdown.txt
awk '{print $1, $2, $6}' $kmer_up > kup.txt 
awk '{print $1, $2, $3}' $kmer_down > kdown.txt 
sed -i 's/^Ha412HOChr//g' snpup.txt snpdown.txt kup.txt kdown.txt

Rscript /mnt/data/DanaS/r_scripts/scrach_miami_kmer.r \
snpup.txt $simple_up $bonfer_up \
snpdown.txt $simple_down $bonfer_down \
"miami" trait_name \
kup.txt kdown.txt




# #####
#
#  NecAD 
#       Yavor (up)
cd /mnt/data/DanaS/may_2024/miami
up=/mnt/data/DanaS/sam_broom_gwas/gwas/yavor_nec_thresh_03/snp/out_gemma/n_NEC_n2r_max_zero_one.assoc.txt
simple_up=$(awk 'NR==2 {print $2}' /mnt/data/DanaS/sam_broom_gwas/gwas/yavor_nec_thresh_03/snp/yavor_NEC_mean_thresh03_simplem.txt)
bonfer_up=$(awk 'NR==2 {print $1}' /mnt/data/DanaS/sam_broom_gwas/gwas/yavor_nec_thresh_03/snp/yavor_NEC_mean_thresh03_simplem.txt)
kmer_up=/mnt/data/DanaS/sam_broom_gwas/gwas/yavor_hel_thresh_05/kmer/get_kmers/ha412/kmers_manhattan_plot_data_pvast.txt
#      Gadot (down)
down=/mnt/data/DanaS/sam_broom_gwas/gwas/gadot_nec_thresh_03/snp/out_gemma/n_NEC_n2r_max_zero_one.assoc.txt
simple_down=$(awk 'NR==2 {print $2}' /mnt/data/DanaS/sam_broom_gwas/gwas/gadot_nec_thresh_03/snp/gadot_NEC_mean_thresh03_simplem.txt)
bonfer_down=$(awk 'NR==2 {print $1}' /mnt/data/DanaS/sam_broom_gwas/gwas/gadot_nec_thresh_03/snp/gadot_NEC_mean_thresh03_simplem.txt)
kmer_down=/mnt/data/DanaS/sam_broom_gwas/gwas/gadot_hel_thresh_05/kmer/get_kmers/ha412/kmers_manhattan_plot_data.txt

awk '{print $1, $3, $13}' $up > snpup.txt
awk '{print $1, $3, $13}' $down > snpdown.txt
awk '{print $1, $2, $6}' $kmer_up > kup.txt 
awk '{print $1, $2, $3}' $kmer_down > kdown.txt 
sed -i 's/^Ha412HOChr//g' snpup.txt snpdown.txt kup.txt kdown.txt
head -n2 snpdown.txt


Rscript /mnt/data/DanaS/r_scripts/scrach_miami_kmer.r \
snpup.txt $simple_up $bonfer_up \
snpdown.txt $simple_down $bonfer_down \
"miami" NecAD.plot \
kup.txt kdown.txt



# #####
#
#  HealAD
#       Yavor (up)
cd /mnt/data/DanaS/may_2024/miami
up=/mnt/data/DanaS/sam_broom_gwas/gwas/yavor_hel_thresh_05/snp/out_gemma/n_HEALTHY_n2r_max_zero_one.assoc.txt
simple_up=$(awk 'NR==2 {print $2}' /mnt/data/DanaS/sam_broom_gwas/gwas/yavor_hel_thresh_05/snp/yavor_mean_thresh05_simplem.txt)
bonfer_up=$(awk 'NR==2 {print $1}' /mnt/data/DanaS/sam_broom_gwas/gwas/yavor_hel_thresh_05/snp/yavor_mean_thresh05_simplem.txt)
kmer_up=/mnt/data/DanaS/sam_broom_gwas/gwas/yavor_hel_thresh_05/kmer/get_kmers/ha412/kmers_manhattan_plot_data_pvast.txt
#      Gadot (down)
down=/mnt/data/DanaS/sam_broom_gwas/gwas/gadot_hel_thresh_05/snp/out_gemma/n_HEALTHY_n2r_max_zero_one.assoc.txt
simple_down=$(awk 'NR==2 {print $2}' /mnt/data/DanaS/sam_broom_gwas/gwas/gadot_hel_thresh_05/snp/gadot_mean_thresh05_simplem.txt)
bonfer_down=$(awk 'NR==2 {print $1}' /mnt/data/DanaS/sam_broom_gwas/gwas/gadot_hel_thresh_05/snp/gadot_mean_thresh05_simplem.txt)
kmer_down=/mnt/data/DanaS/sam_broom_gwas/gwas/gadot_hel_thresh_05/kmer/get_kmers/ha412/kmers_manhattan_plot_data.txt

awk '{print $1, $3, $13}' $up > snpup.txt
awk '{print $1, $3, $13}' $down > snpdown.txt
awk '{print $1, $2, $6}' $kmer_up > kup.txt 
awk '{print $1, $2, $3}' $kmer_down > kdown.txt 
head -n2 snpdown.txt
sed -i 's/^Ha412HOChr//g' snpup.txt snpdown.txt kup.txt kdown.txt

Rscript /mnt/data/DanaS/r_scripts/scrach_miami_kmer.r \
snpup.txt $simple_up $bonfer_up \
snpdown.txt $simple_down $bonfer_down \
"miami" HealAD.plot \
kup.txt kdown.txt





python3 /mnt/data/DanaS/sam_broom_gwas/scripts/plot_chr_region_ld_sorted_heatmap.py \
  --raw  /mnt/data/DanaS/sam_broom_gwas/gwas/gadot_hel_thresh_05/snp/snps_chr9.raw \
  --bim  /mnt/data/DanaS/sam_broom_gwas/gwas/gadot_hel_thresh_05/snp/chr9.bim \
  --pav  /mnt/data/DanaS/sam_broom_gwas/gwas/pav_chr9.tsv \
  --allrefs /mnt/data/DanaS/sam_broom_gwas/significants_with_snp_mapping/ALL_REFS.all_markers.intersected250kbed_prod.tsv \
  --ref HA412 --project gadot --trait hel \
  --chr-bim 9 \
  --chr-pattern Ha412HOChr09 \
  --start 166000000 --end 168000000 \
  --pheno /mnt/data/DanaS/sam_broom_gwas/gwas/gadot_hel_thresh_05/snp/PHENO_FOR_GWAS.txt \
  --sample-col File \
  --pheno-col n_HEALTHY_n2r_mean \
  --ld-matrix /path/to/gwaslab_chr9_ld_matrix.tsv \
  --intro-bed /path/to/introgressions_HA412.bed \
  --intro-chr Ha412HOChr09 \
  --out-png /mnt/data/DanaS/sam_broom_gwas/gwas/gadot_hel_thresh_05/plots/chr9_ld_sorted_heatmap.png \
  --out-meta /mnt/data/DanaS/sam_broom_gwas/gwas/gadot_hel_thresh_05/plots/chr9_ld_sorted_heatmap.meta.tsv
