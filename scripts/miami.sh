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

Rscript /mnt/data/DanaS/r_scripts/scrach_miami_kmer.R \
snpup.txt $simple_up $bonfer_up \
snpdown.txt $simple_down $bonfer_down \
"miami" HealAD.plot \
kup.txt kdown.txt



