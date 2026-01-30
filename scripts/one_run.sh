
############################################# filter  ##############################
PHENO=/mnt/data/DanaS/may_2024/yavor/by_n2r/mean/thresh_05/filtered_by_n2r_HEALTHY_mean_thresh_0.5.txt
OUTDIR=/mnt/data/DanaS/may_2024/yavor/by_n2r/mean/thresh_05/xrq
VCF=/mnt/data/sunflower/SAM_SNPs/Annuus.tranche90.snp.fullsam.90.bi.remappedxrqv2.vcf.gz
PRE=yavor_HEL_mean_thresh05
#####################################################################################


conda deactivate
thrsh=0.5
threshfldr=thresh_05

PHENO=/mnt/data/DanaS/may_2024/yavor/by_n2r/mean/thresh_05/filter_by_total/$threshfldr/total_tubes_n2r_mean_yavor_"$trsh".txt
OUTDIR=/mnt/data/DanaS/may_2024/yavor/by_n2r/mean/thresh_05/filter_by_total/"$threshfldr"/
VCF=/mnt/data/sunflower/SNP_HA412/SAM.HA412.maf5Miss20Tranch99.vcf.gz
PRE="$threshfldr"
#add "SAM" and three leading zeros
awk '{printf "SAM%03d", $1; for (i=2; i<=NF; i++) printf "\t%s", $i; print ""}' "$PHENO" > temp_file && mv temp_file "$PHENO"
cd $OUTDIR
# remove 272 and 85
awk -i inplace '!($1=="SAM272"|| $1=="SAM085")' $PHENO

# create samples2keep file
awk -F "\t" 'NR>1 { print $1}' $PHENO > samples2keep.txt
keep=samples2keep.txt

# filters
#MAF=0.05
#MISS=0.8
filtered_vcf="$PRE"_maf5.miss20
# filter vcf and assign id for the snps
vcftools --gzvcf $VCF --keep $keep --maf 0.05 --max-missing 0.8 --recode --recode-INFO-all --out $filtered_vcf
filtered_vcf="$PRE"_maf5.miss20.recode.vcf
vcf_id="$PRE"_maf5.miss10.with_snp_id
bcftools annotate --set-id +'%CHROM\_%POS\_%REF\_%FIRST_ALT' $filtered_vcf > $vcf_id
rm $filtered_vcf

vcf_id="$PRE"_maf5.miss10.with_snp_id
plink --vcf $vcf_id --make-bed --allow-extra-chr --out "$PRE"
mkdir tmp
mv "$PRE".fam ./tmp

# merge the phenotype file with the .fam file
Rscript /mnt/data/DanaS/r_scripts/MergeFamFile.R ./tmp/"$PRE".fam $PHENO $OUTDIR "$PRE"

# reletedness matrix
gemma -bfile "$PRE" -gk 1 -o "$PRE"

# PCA 
#GO TO pca DIRECTORY and run the Rscript thats located in /mnt/data/DanaS/GWAS2/VCF/PCA_IBS.R, specify the vcf file and  population.file
#Rscript /mnt/data/DanaS/r_scripts/PCA_IBS_ONE_RUN.R $vcf_id pop_code.txt "$PRE"
pop_code=/mnt/data/DanaS/pop_code.txt
awk 'BEGIN {OFS="\t"} {print $1, $2}' $pop_code > POP_FOR_GWAS.txt

awk 'NR==FNR { keep[$1]; next } FNR > 1 && $1 in keep { print $2 }' POP_FOR_GWAS.txt samples2keep.txt > POP_FOR_PCA.txt

#this will create eigval.txt and eigvec.txt files in /pca directory the first col needs to be chnged to 1
#change first col values to 1 | and Delete first (headers) line:

#awk '{$1=1 ; print ;}' ./output/"$PRE".eigVec.txt | sed '1d' > "$PRE".eigVec_ones_NO_HEAD.txt
mkdir pca
plink --bfile "$PRE" --allow-extra-chr --keep-allele-order --pca 10 header --out ./pca/pca_"$PRE"
cd pca
awk 'NR > 1 { $1=""; $2 = ""; print 1,$0 }' ./pca_"$PRE".eigenvec | tr -s " " > pca.txt
cut -d" " -f1-2 pca.txt > PC1.txt
cut -d" " -f1-3 pca.txt > PC1-2.txt
cut -d" " -f1-4 pca.txt > PC1-3.txt
cut -d" " -f1-5 pca.txt > PC1-4.txt
cut -d" " -f1-6 pca.txt > PC1-5.txt
head PC1.txt
pca=./pca/PC1-2.txt
cd ..
plink --vcf $VCF --allow-extra-chr --keep-allele-order --pca 10 header --out pca
# simplem
plink --vcf $vcf_id --allow-extra-chr --vcf-half-call missing --maf 0.05 --max-maf 0.49 --recode 12 --transpose --out OUT_"$PRE" # WAS USED BY HOD
cut -d' ' -f7- OUT_"$PRE".tped > OUT_"$PRE".tab
#this will create a file with bonferoni in first col and simplem in second col
Rscript /mnt/data/DanaS/r_scripts/simplem_4_one_run.R OUT_"$PRE".tab "$PRE"

mkdir out_gemma
conda deactivate # if base is activate
conda activate snakemake
cp /mnt/data/DanaS/may_2024/yavor/by_n2r/mean/thresh_05/filter_by_total/thresh_01/snakefile ./


snakemake --cores 6 --snakefile ../snakefile

mkdir plots
cp /mnt/data/DanaS/may_2024/yavor/by_n2r/mean/thresh_05/filter_by_total/thresh_01/snake_manhattan.smk ./
snakemake --cores 6 --snakefile ../snake_manhattan.smk
mv *.jpeg ./plots

# gemma
mkdir out_gemma
OUTGEMMA=./out_gemma
relat=./output/*.cXX.txt

while read -r col trait 
do
    gemma -bfile "$PRE" \
    -k $relat \
    -lmm 4 \
    -n $col -o $trait \
    -outdir $OUTGEMMA
done < TraitMap.txt

mkdir $OUTGEMMA/preplot
cp ./out_gemma/*.assoc.txt ./out_gemma/preplot/

# if only HanXRQChr than do this 
sed -i 's/^Ha412HOChr//g' ./out_gemma/preplot/*.assoc.txt #for Ha412 

conda activate r_4.1
bonfer=$(awk 'NR==2 {print $1}' *_simplem.txt)
simple=$(awk 'NR==2 {print $2}' *_simplem.txt)
for trait in $(awk '{print $2}' TraitMap.txt); do
Rscript /mnt/data/DanaS/r_scripts/Plot.Manhattan.QQ.R ./out_gemma/$trait.assoc.txt $trait $bonfer $simple
done

mkdir plots
mv *.jpeg ./plots
conda deactivate


# to skip the first 10 lines and then use while read -r :
while read -r col trait 
do
    gemma -bfile "$PRE" \
    -k $relat \
    -lmm 4 \
    -n $col -o $trait \
    -outdir $OUTGEMMA
done

#to remove chrachters (Ha412HOChr/HanXRQChr) from CHR col so it will be numeric - FOR ALL THE FILES IN THE DIRECTORY
# if chr col has boath HanXRQChr and HanXRQChr00c do this:
mkdir ./out_gemma/preplot
for file in ./out_gemma/preplot/*.assoc.txt; do
    sed -i '/^HanXRQChr00c/d' "$file"
    sed -i 's/^HanXRQChr//g' "$file"
done

#________________________________________________________
# make sure all chr are numbers only
cut -f1 ./out_gemma/preplot/*.assoc.txt | sort | uniq


#  plot manhattan and QQ ## bonferoni=numLoci,simplem


screen -dmS plot_man_$PRE
conda activate r_4.1
bonfer=$(awk 'NR==2 {print $1}' *_simplem.txt)
simple=$(awk 'NR==2 {print $2}' *_simplem.txt)
for trait in $(awk '{print $2}' TraitMap.txt); do
Rscript /mnt/data/DanaS/r_scripts/Plot.Manhattan.QQ.R ./out_gemma/preplot/$trait.assoc.txt $trait $bonfer $simple
done

mkdir plots
mv *.jpeg ./plots

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

# get the wanted cols in this order
down=/mnt/data/DanaS/zero_one/gadot/z_25/out_gemma/preplot/n_HEALTHY_mean_zero_one.assoc.txt
up=/mnt/data/DanaS/zero_one/z_25/snps/gemma/out_gemma/preplot/n_HEALTHY_mean_zero_one.assoc.txt
simple=$(awk 'NR==2 {print $2}' /mnt/data/DanaS/zero_one/gadot/z_25/*_simplem.txt)

down=/mnt/data/DanaS/zero_one/gadot/filter_by_NEC/out_gemma/preplot/n_NEC_mean_zero_one.assoc.txt
up=/mnt/data/DanaS/zero_one/filtered_by_NEC/z_25/out_gemma/preplot/n_NEC_mean_zero_one.assoc.txt


awk '{print $2, $1, $3, $8, $9, $10, $13}' $down > gadot_miami.txt
awk '{print $2, $1, $3, $8, $9, $10, $13}' $up > yavur_miami.txt

Rscript /mnt/data/DanaS/r_scripts/miami_plot.r yavur_miami.txt gadot_miami.txt $simple n_hel_zeroone_yavur_gadot

Rscript /mnt/data/DanaS/r_scripts/taitmap.R $PHENO

assocfile=/mnt/data/DanaS/zero_one/gadot/out_gemma/preplot/n_N_max.assoc.txt
sed -i 's/^Ha412HOChr//g' $assocfile
Rscript /mnt/data/DanaS/r_scripts/Plot.Manhattan.QQ.R $assocfile n_N_max $bonfer $simple
# bayesian plot manhattan 
mkdir preplot
cp *.param.txt ./preplot

cd preplot
# while in preplot
# delete last line from col_trait (cuase it doesnt have gemma.asoc result)
sed '$d' /mnt/data/DanaS/gwas_gadot_ha412/run_gwas/stg_test/stg_col_trt.txt > ./preplot/stg_col_trt.txt
#remove chrachters (Ha412HOChr) from CHR col so it will be numeric - FOR ALL THE FILES IN THE DIRECTORY
sed -i 's/^Ha412HOChr//g' *.param.txt 

for trait in $(awk '{print $2}' stg_col_trt.txt); do
Rscript /mnt/data/DanaS/r_scripts/QQMan_bslmm.R \
/mnt/data/DanaS/gwas_gadot_ha412/run_gwas/stg_test/stg_blsmm/preplot/$trait.param.txt \
$trait
done

mkdir plots
mv *.jpeg ./plots
ls ./plots/Man* | wc -l
ls ./out_gemma/preplot/* | wc -l
screen -ls



# Variables
gwas_file=/mnt/data/DanaS/zero_one/gadot/out_gemma/n_HEALTHY_mean_zero_one.assoc.txt
chromosome="13"
num_results=10

# Extract and sort by p-value
awk -v chr="$chromosome" '$1 == chr' "$gwas_file" | sort -k13,13nr | head -n "$num_results"


dos2unix $PHENO

PRE=gadot_z_025
PHENO=/mnt/data/DanaS/may_2024/gadot/zero_one/by_n2r/z_025/HEALTHY_MAX_N2R_Z_025.txt
OUTDIR=/mnt/data/DanaS/may_2024/gadot/zero_one/by_n2r/z_025
VCF=/mnt/data/sunflower/SNP_HA412/SAM.HA412.maf5Miss20Tranch99.vcf.gz


# create samples2keep file
awk -F "\t" 'NR>1 { print $1}' $PHENO > samples2keep.txt
keep=samples2keep.txt

# filters
MAF=0.05
MISS=0.8
filtered_vcf="$PRE"_maf5.miss20
# filter vcf and assign id for the snps
vcftools --gzvcf $VCF --keep $keep --maf $MAF --max-missing $MISS --recode --recode-INFO-all --out $filtered_vcf
filtered_vcf="$PRE"_maf5.miss20.recode.vcf
vcf_id="$PRE"_maf5.miss10.with_snp_id
bcftools annotate --set-id +'%CHROM\_%POS\_%REF\_%FIRST_ALT' $filtered_vcf > $vcf_id

plink --vcf $vcf_id --make-bed --allow-extra-chr --out "$PRE"
mkdir tmp
mv "$PRE".fam ./tmp


# merge the phenotype file with the .fam file
Rscript /mnt/data/DanaS/r_scripts/MergeFamFile.R ./tmp/"$PRE".fam $PHENO $OUTDIR "$PRE"

# reletedness matrix
gemma -bfile "$PRE" -gk 1 -o "$PRE"

relat=./output/z_25.cXX.txt

# PCA 
#GO TO pca DIRECTORY and run the Rscript thats located in /mnt/data/DanaS/GWAS2/VCF/PCA_IBS.R, specify the vcf file and  population.file
#Rscript /mnt/data/DanaS/r_scripts/PCA_IBS_ONE_RUN.R $vcf_id pop_code.txt "$PRE"


#this will create eigval.txt and eigvec.txt files in /pca directory the first col needs to be chnged to 1
#change first col values to 1 | and Delete first (headers) line:

#awk '{$1=1 ; print ;}' ./output/"$PRE".eigVec.txt | sed '1d' > "$PRE".eigVec_ones_NO_HEAD.txt

# Simple.M
#CREATE PLINK .tped FILE FROM VCF FILE AND CONVERT IT TO .tab FILE:
plink --vcf $vcf_id --allow-extra-chr --vcf-half-call missing --maf 0.05 --max-maf 0.49 --recode 12 --transpose --out OUT_"$PRE" # WAS USED BY HOD
cut -d' ' -f7- OUT_"$PRE".tped > OUT_"$PRE".tab
Rscript /mnt/data/DanaS/r_scripts/simplem_4_one_run.R OUT_"$PRE".tab "$PRE" #this will create a file with bonferoni in first col and simplem in second col


# gemma
mkdir out_gemma
OUTGEMMA=./out_gemma

screen -S "$PRE"
OUTGEMMA=./out_gemma
relat=./output/*.cXX.txt
PRE=gadot_z_025
while read -r col trait 
do
    gemma -bfile "$PRE" \
    -k $relat \
    -lmm 4 \
    -n $col -o $trait \
    -outdir $OUTGEMMA
done < TraitMap.txt


mkdir $OUTGEMMA/preplot
cp ./out_gemma/*.assoc.txt ./out_gemma/preplot/


# if only HanXRQChr than do this 
sed -i 's/^Ha412HOChr//g' ./out_gemma/preplot/*.assoc.txt #for Ha412 

conda activate r_4.1
bonfer=$(awk 'NR==2 {print $1}' *_simplem.txt)
simple=$(awk 'NR==2 {print $2}' *_simplem.txt)
for trait in $(awk '{print $2}' TraitMap.txt); do
Rscript /mnt/data/DanaS/r_scripts/Plot.Manhattan.QQ.R ./out_gemma/preplot/$trait.assoc.txt $trait $bonfer $simple
done
k=/mnt/data/DanaS/sam_broom_gwas/gwas/aeg/kmer/get_kmers/tot/ha412/kmers_manhattan_plot_data.txt

Rscript /mnt/data/DanaS/r_scripts/Plot.Manhattan.QQ.R /mnt/data/DanaS/sam_broom_gwas/gwas/aeg/out_gemma/total_n2r_max.assoc.txt "n_tot_n2r_max" $bonfer $simple $k
mkdir plots
mv *.jpeg ./plots

# calculate LD Decay with plink
sed -i 's/^Ha412HOChr//g' /mnt/data/DanaS/may_2024/gadot/by_n2r/mean/thresh_05/gadot_mean_thresh05.bed
sed -i 's/^Ha412HOChr//g' /mnt/data/DanaS/may_2024/gadot/by_n2r/mean/thresh_05/gadot_mean_thresh05.bim
sed -i 's/^Ha412HOChr//g' /mnt/data/DanaS/may_2024/gadot/by_n2r/mean/thresh_05/gadot_mean_thresh05.fam

bfiles=gadot_mean_thresh05
vcf=/mnt/data/DanaS/may_2024/yavor/by_n2r/mean/thresh_05/yavor_mean_thresh05_maf5.miss10.with_snp_id
out=ld_decay_thresh05
plink \
    --bfile $bfiles --r2 gz \
    --ld-window-kb 1000 --ld-window-r2 0 --ld-window 100 --out ld_decay
# Calculating average LD across set distances
# run python script
conda deactivate
conda activate python_env

python /mnt/data/DanaS/r_scripts/ld_decay_calc.py -i ld_decay.ld.gz -o ld_decay
#plot
conda activate r_4.1
Rscript /mnt/data/DanaS/r_scripts/ld_decay_plot.r /mnt/data/DanaS/may_2024/gadot/by_n2r/mean/thresh_05/ld_decay.ld_decay_bins

# calculate LD Decay with popladecay
/mnt/data/Install/PopLDdecay/PopLDdecay -InVCF $vcf -MaxDist 1000 -OutStat $out

perl /mnt/data/Install/PopLDdecay/bin/Plot_OnePop.pl -inFile ld_decay_thresh05.stat.gz -output ld_decay_thresh05

grep "Ha412HOChr07" /mnt/data/DanaS/may_2024/yavor/by_n2r/mean/thresh_05/out_gemma/n_HEALTHY_n2r_mean_zero_one.assoc.txt | awk -F'\t' '$3 == 157328855 {print $13}'
grep "^Ha412HOChr07" ./out_gemma/n_HEALTHY_n2r_mean_zero_one.assoc.txt | awk -F'\t' '$3 >= 143267372 && $3 <= 143298531 {print $13}'
grep "^Ha412HOChr07" ./out_gemma/n_HEALTHY_n2r_mean_zero_one.assoc.txt | awk -F'\t' '$13 >= 8.06047289182362111e-08 && $13 <= 8.06047289182362111e-07 {print $2}'

grep "^Ha412HOChr04" /mnt/data/DanaS/may_2024/gadot/by_n2r/mean/thresh_05/out_gemma/n_HEALTHY_norm2root_max.assoc.txt | awk -F'\t' '$13 >= 8.06047289182362111e-08 && $13 <= 8.06047289182362111e-07 {print $2}'
# get significant snps
awk -v simple="$simple" 'BEGIN{FS=OFS="\t"}; $13 < simple' ../out_gemma/n_HEALTHY_n2r_max_zero_one.assoc.txt > ./sig_snps/sig_snps.txt
awk 'BEGIN{FS=OFS="\t"}; $13 < 8.02458733559627e-07' ./out_gemma/n_NEC_n2r_max_zero_one.assoc.txt > ./sig_snp/sig_snps.txt 

 # annotate
gff_file=/mnt/data/sunflower/HA412/Ha412HOv2.0-20181130.gff3.gz
xrq_gff=/mnt/data/DanaS/XRQ2_functionalAnnotationForDana.gtf
ha412_gff=/mnt/data/DanaS/HAN412_Eugene_curated_v1_1.gff3
zgrep -i "HanXRQChr02g0050061" $gff_file
grep -i "HanXRQr2_Chr07g0314271" $xrq_gff
grep -i "HanXRQChr07g0314271" $ha412_gff
HanXRQChr07g0314341 # looks like this in ha412_gff
HanXRQr2_Chr07g0314271 # looks like this in xrq_gff

#  xrq_gff looks like this:
# HanXRQCP        EuGene  exon    463     1524    .       -       .       transcript_id "mRNA:HanXRQr2_CPg0835831"; gene_id "gene:HanXRQr2_CPg0835831"; gene_name "HanXRQr2CPg0835831"; Ontology_term "SO:0005845"; est_cons "30.7"; est_incons "0.0"; Name "HanXRQr2CPg0835831"; locus_tag "HanXRQr2_CPg0835831"; product "Putative%20photosystem%20II%20protein%20D1%2FD2%20superfamily"; Ontology_term "GO:0009523,GO:0009535,GO:0010287,GO:0016021,GO:0005506,GO:0005515,GO:0016168,GO:0016682,GO:0045156,GO:0009635,GO:0009772,GO:0018298"; ec_number "1.10.3.9";

# h412_gff looks like this:
# Ha412HOChr02	EuGene	gene	18236122	18237527	.	+	.	Name=Ha412HOChr02g0052541;ID=gene:Ha412HOChr02g0052541;Expression_Levels=786;PSC8_Gene=HanPSC8Chr02g0049571-HanPSC8Chr02g0049551-HanPSC8Chr02g0049431-HanPSC8Chr02g0049491-HanPSC8Chr02g0049481-HanPSC8Chr02g0049531;XRQv2_Gene=HanXRQChr02g0050201-HanXRQChr02g0050181-HanXRQChr02g0050001-HanXRQChr02g0050071-HanXRQChr02g0050061-HanXRQChr02g0050111	Ha412HOChr02	18231906	18236906	785

# need to get gene id from ha412 and cross it with xrq
# create 2 tables with all relevant info: chr, start pos, stop pos,  gene id, gene name, GO term from both files

# so first 
################ working #####################333
mkdir -p ./annotated_snps
mkdir -p ./annotated_genes
mkdir -p ./proces_snp

for file in ./sig_snps/*.txt; do
    # Define the output filename
    output_file="./annotated_snps/$(basename "$file" .txt).txt"
    
    # Use awk to process each file
    awk 'BEGIN {FS="\t"; OFS="\t"} NR>1 {gsub("_", "\t", $2); print $2}' "$file" > "$output_file"
done

snps=./sig_snps/sig_snps.txt
awk 'BEGIN {FS="\t"; OFS="\t"} {print $1,  $3, int($3)-25000, int($3)+25000}' sig_snps.txt | awk '{print $1 "\t" $3 "\t" $4}' > ./snp_ld_intervl.txt

awk '{print $1 "\t" $2-1 "\t" $3}' snp_ld_intervl.txt > snp_coordinates.bed
bedtools intersect -a $ha412_gff -b snp_coordinates.bed -wa > output_412.gff
bedtools intersect -a $ha412_gff -b snp_coordinates.bed -wa | awk '$3 == "gene"' > snp2gene_412.gff
bedtools intersect -a $ha412_gff -b snp_coordinates.bed -wa | awk '$3 == "mRNA"' > snp2mRNA_412.gff

HanXRQChr07g0314341 # looks like this in ha412_gff
HanXRQr2_Chr07g0314271 # looks like this in xrq_gff

#so add "r2_" and gerp in xrq_gff
sed 's/old_string/new_string/g' input_file > output_file
sed 's/HanXRQChr/HanXRQr2_Chr/g' output_412.gff > snp2gene_412_4xrq.gff
# do the grep
#grep the patterns
grep "HanXRQr2_Chr" snp2gene_412_4xrq.gff > patterns_to_extrct_from_xrq.txt

# Extract patterns from patterns_to_extrct_from_xrq.txt
awk -F'XRQv2_Gene=' '{print $2}' patterns_to_extrct_from_xrq.txt | awk -F';' '{print $1}' > patterns.txt

# Output file for matching entries
output_file="matching_entries.txt"
rm -f $output_file  # Clear the output file if it exists

# Search for each pattern in xrq_gff and append matches to the output file
while IFS= read -r pattern; do
    grep "$pattern" $xrq_gff >> $output_file
done < patterns.txt

echo "Matching entries have been written to $output_file"
###########################

for FN in ./annotated_snps/*.txt; do

  output_fn="./proces_snp/$(basename "$FN" .txt).txt"
  awk 'BEGIN{FS=OFS="\t"}{print $1, int($2)-25000, int($2)+25000}' "$FN" > "$output_fn"
done

for SN in ./proces_snp/*.txt; do
  intersect_output="./annotated_genes/$(basename "$SN" .txt).txt"
  bedtools intersect -wo -a "$gff_file" -b "$SN" | awk '$3 == "gene"' > "$intersect_output"
done 

rm -r ./proces_snp
 ### END OF WORKING  CODE ######  




#________ EMMAX _________
conda activate emmax


wget https://csg.sph.umich.edu//kang/emmax/download/index.html
plink --vcf $VCF --keep $SAMPLES --set-missing-var-ids @:# --maf 0.05 --allow-extra-chr --recode12 --output-missing-genotype 0 --transpose --out $NAME
mkdir emmax
cd emmax

NAME=yavor_thresh025_emmax
VCF=/mnt/data/DanaS/may_2024/yavor/by_n2r/mean/thresh_025/yavor_mean_thresh025_maf5.miss10.with_snp_id
PHENO=/mnt/data/DanaS/may_2024/yavor/by_n2r/mean/thresh_025/filtered_by_n2r_HEALTHY_mean_thresh_0.25.txt
mkdir phenotypes 

# Case-control Phenotypes
# If you encode case/control to 2/1, you will be able to run case-control analysis. 
# Because EMMAX is based on linear mixed model rather than generalized mixed model, 
# the effect size (beta) would not be meaningful, but the p-values should be reliable 
# (unless case/control counts are highly imbalanced).

# Get the header line (column names)
header=$(head -1 "$PHENO")
# Convert the header line into an array of column names
IFS=$'\t' read -r -a col_names <<< "$header"
# Loop through each column index and name
for i in "${!col_names[@]}"; do
    COL_NAME="${col_names[$i]}"
    # Create a file for each column
    awk -F "\t" -v col=$((i+1)) 'NR>1 { print $1, $1, $col }' "$PHENO" > ./phenotypes/"${COL_NAME}.txt"
done

for i in "${!col_names[@]}"; do
    COL_NAME="${col_names[$i]}"
    # Create a file for each column
    awk -F "\t" -v col=$((i+1)) 'NR>1 { print $1, $1, $col }' "$PHENO" > ./phenotypes/"${COL_NAME}.txt"
done

# plink
plink --vcf $VCF --set-missing-var-ids @:# --maf 0.05 --allow-extra-chr --recode12 --output-missing-genotype 0 --transpose --out "$NAME"
# Kinship Matrix
sed -i 's/^Ha412HOChr//g' ${NAME}.tped
sed -i 's/^Ha412HOChr//g' ${NAME}.map
emmax-kin-intel64 -v -s -d 10 $NAME

# EMMAX Association
mkdir results
for PHENOFILE in ./phenotypes/*.txt; do
    emmax-intel64 -v -d 10 -t $NAME -p $PHENOFILE -k yavor_thresh025_emmax.aIBS.kinf -o "$PHENOFILE"_res
done


for FILE in ./phenotypes/*.ps; do
    Rscript /mnt/data/DanaS/r_scripts/manhattan_emmax.R "$FILE" ./plots $bonfer $simple
done

bonfer=$(awk 'NR==2 {print $1}' /mnt/data/DanaS/may_2024/yavor/by_n2r/mean/thresh_025/*_simplem.txt)
simple=$(awk 'NR==2 {print $2}' /mnt/data/DanaS/may_2024/yavor/by_n2r/mean/thresh_025/*_simplem.txt)
for trait in $(awk '{print $2}' TraitMap.txt); do
Rscript /mnt/data/DanaS/r_scripts/Plot.Manhattan.QQ.R ./out_gemma/preplot/$trait.assoc.txt $trait $bonfer $simple
done

# get significant snps
#!/bin/bash

# Ensure the output directory exists
mkdir -p ./significants

# Loop through all .assoc.txt files
for ASSOC in ./gwas/*/snp/out_gemma/*.assoc.txt; do
    # Extract the base directory for the current assoc file
    DIR=$(dirname "$ASSOC")
    BASE_DIR=$(dirname "$DIR") # Go one level up to find the directory containing *_simplem.txt
    POP_NAME=$(basename "$BASE_DIR") # Extract population name (e.g., gadot_nec_thresh_03)

    # Locate the corresponding *_simplem.txt file
    SIMPLE_FILE=$(find "$BASE_DIR" -maxdepth 1 -name "*_simplem.txt")

    # Check if the corresponding simplem file exists
    if [[ -f "$SIMPLE_FILE" ]]; then
        # Extract the 'simple' threshold value
        simple=$(awk 'NR==2 {print $2}' "$SIMPLE_FILE")
    else
        echo "No *_simplem.txt file found for $ASSOC. Skipping..."
        continue
    fi

    # Prepare the output file
    OUTPUT="./significants/${POP_NAME}_$(basename "$ASSOC" .assoc.txt)_sig_snps.txt"
    
    # Initialize the output file with a header
    echo -e "chr\tpos\tpval" > "$OUTPUT"
    
    # Extract significant SNPs and append them to the output file
    awk -v simple="$simple" 'BEGIN{FS=OFS="\t"} $13 < simple {print $1, $2, $13}' "$ASSOC" >> "$OUTPUT"
done

echo "Significant SNPs extracted to ./significants/ for all populations and assoc files."


mkdir -p ./max/{thresh_001,thresh_002,thresh_003,thresh_004,thresh_005}

trsh="05"
conda deactivate
############################################ aeg try ##############################
PHENO=/mnt/data/DanaS/sam_broom_gwas/gwas/aeg/thresh/max/thresh_0"$trsh"/filtered_by_n_HEALTHY_max_thresh_0."$trsh".txt
OUTDIR=/mnt/data/DanaS/sam_broom_gwas/gwas/aeg/thresh/max/thresh_0"$trsh"
VCF=/mnt/data/sunflower/SNP_HA412/SAM.HA412.maf5Miss20Tranch99.vcf.gz
PRE=aeg_HEALTHY_max_thresh_0"$trsh"
#########################################################################################
cd $OUTDIR
awk '{printf "SAM%03d", $1; for (i=2; i<=NF; i++) printf "\t%s", $i; print ""}' "$PHENO" > temp_file && mv temp_file "$PHENO"
# remove 272 and 85
awk -i inplace '!($1=="SAM272"|| $1=="SAM085")' $PHENO

# create samples2keep file
awk -F "\t" 'NR>1 { print $1}' $PHENO > samples2keep.txt
keep=samples2keep.txt
head $keep
# filters
#MAF=0.05
#MISS=0.8
filtered_vcf="$PRE"_maf5.miss20
# filter vcf and assign id for the snps
vcftools --gzvcf $VCF --keep $keep --maf 0.05 --max-missing 0.8 --recode --recode-INFO-all --out $filtered_vcf
filtered_vcf="$PRE"_maf5.miss20.recode.vcf
vcf_id="$PRE"_maf5.miss10.with_snp_id
bcftools annotate --set-id +'%CHROM\_%POS\_%REF\_%FIRST_ALT' $filtered_vcf > $vcf_id
rm $filtered_vcf

vcf_id="$PRE"_maf5.miss10.with_snp_id
plink --vcf $vcf_id --make-bed --allow-extra-chr --out "$PRE"
mkdir tmp
mv "$PRE".fam ./tmp

# merge the phenotype file with the .fam file
Rscript /mnt/data/DanaS/r_scripts/MergeFamFile.R ./tmp/"$PRE".fam $PHENO $OUTDIR "$PRE"

# reletedness matrix
gemma -bfile "$PRE" -gk 1 -o "$PRE"

# simplem
plink --vcf $vcf_id --allow-extra-chr --vcf-half-call missing --maf 0.05 --max-maf 0.49 --recode 12 --transpose --out OUT_"$PRE" # WAS USED BY HOD
cut -d' ' -f7- OUT_"$PRE".tped > OUT_"$PRE".tab
#this will create a file with bonferoni in first col and simplem in second col
Rscript /mnt/data/DanaS/r_scripts/simplem_4_one_run.R OUT_"$PRE".tab "$PRE"

# gemma
mkdir out_gemma
OUTGEMMA=./out_gemma
relat=./output/*.cXX.txt
while read -r col trait 
do
    gemma -bfile "$PRE" \
    -k $relat \
    -lmm 4 \
    -n $col -o $trait \
    -outdir $OUTGEMMA
done < TraitMap.txt


conda activate r_4.1
bonfer=$(awk 'NR==2 {print $1}' *_simplem.txt)
simple=$(awk 'NR==2 {print $2}' *_simplem.txt)
for trait in $(awk '{print $2}' TraitMap.txt); do
Rscript /mnt/data/DanaS/r_scripts/Plot.Manhattan.QQ.R ./out_gemma/$trait.assoc.txt $trait $bonfer $simple
done

mkdir plots
mv *.jpeg ./plots

#for file in /mnt/data/DanaS/sam_broom_gwas/gwas/aeg/thresh/max/*.txt; do
#wc -l $file;
#done
trsh="04"
