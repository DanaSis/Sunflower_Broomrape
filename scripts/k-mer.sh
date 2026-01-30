###################### k-mer anlysis ########################

# k-mer table path/prefix

kmer_table=/mnt/data/sunflower/sunflowerKmer/kmers_table

## kmers_table.table is the binary file with all the k-mers presence/absence information.


# k-mer based kinship path

/mnt/data/sunflower/sunflowerKmer/kmers_table.kinship


/mnt/data/sunflower/sunflowerKmer/src/kmers_table_to_bed.cpp # this if we want to run the kmers in gemman like we did with SNPs

######################################################### pheno file #####################################################

# Format of phenotype file: the phenotype file should be with two columns separated by tabs (\t), with “accession_id” and “phenotype_value”
# so copy the original pheno file and create a new onw with the wanted colunms and col names:

cp /mnt/data/DanaS/gwas_remappedHa412/pheno/pheno_030822.txt /mnt/data/DanaS/kmer/pheno # copy the original

awk 'BEGIN{FS=OFS="\t"; print "accession_id", "phenotype_value"} NR>1{print $1, $2}' pheno.txt > pheno_for_kmer.txt # create a new onw with the wanted colunms and col names:

# or, can do the same with awk + sed: (but the previos one is one line and more elegant)
awk '{ print $1,  $2 }' pheno.txt > pheno_for_kmer.txt # create a new one with the first and second col
sed '1{ s/accesion/accession_id/; s/n_HEL/phenotype_value/; }' pheno_for_kmer.txt # chnge the headers

######################################################## Run k-mers-based GWAS ##########################################################

python2.7 <KMERS-GWAS-PATH>/kmers_gwas.py --pheno phenotype.pheno --kmers_table kmers_table -l 31 -p 8 --outdir output_dir

wc -l /mnt/data/DanaS/kmer/broomrape/kmers/pass_threshold_5per /mnt/data/DanaS/kmer/broomrape/kmers/pass_threshold_10per

# python2.7 /mnt/data/sunflower/sunflowerKmer/kmers_gwas.py --help
#Parameters used:
# --pheno - a file with numerical phenotypic information (format
# described in “k-mers table conversion to PLINK binary format“
# section)
# --kmers_table - path to k-mers table - only 
# -l - k-mers length
# -p - maximum number of threads to use
# --outdir - path to a directory for output
# -k, --kmers_number - number of k-mers to filter from first step (defualt is 10001)

#Regarding positioning the k-mers, there is no one way to do this. You can look at the examples in our paper, where we tried 
#different things. To summarise I would say, first option is to map k-mers, second would be to find the reads from which the k-mer 
#of interest originated and map the reads, and finally if the previous two didn't work, you can assemble all reads which are linked 
#to the same k-mer and then try to find the assembled fragment position. As this method is not reference dependent in some case the 
#k-mers could just not map to the reference, which are the most interesting cases.


# use filter_table for textual output of k-mers presence/absence pattern. this will output a .txt file with all accessions and 
# presence/absence 1/0
#filter_kmers -t kmers_table -k kmers_list.txt -o output.txt

#Parameters used:
#-t - k-mers table prefix path
#-k - file with k-mers, each k-mer in a seperate line
#-o - output file

## process the kmer file (pass_threshold_10per) for alignment:

# remove all col but the second
awk 'BEGIN{FS=OFS="\t"; print "kmerseq"} NR>1{print $2}' /mnt/data/DanaS/kmer/broomrape/kmers/pass_threshold_10per > /mnt/data/DanaS/kmer/broomrape/kmers/pass_threshold_10per_kmers_list.txt
# remove enithing after the  underscort _
sed 's/\_.*//' /mnt/data/DanaS/kmer/broomrape/kmers/pass_threshold_10per_kmers_list.txt > /mnt/data/DanaS/kmer/broomrape/kmers/pass_threshold_10per_kmers_list_no_under.txt
# remove the header
sed '1d' /mnt/data/DanaS/kmer/broomrape/kmers/pass_threshold_10per_kmers_list_no_under.txt > /mnt/data/DanaS/kmer/broomrape/kmers/pass_threshold_10per_kmers_list_no_under_no_header.txt


        
#convert the kmers_list.txt to .fasts file
awk '{ print ">"NR"\n"$0 }' /mnt/data/DanaS/kmer/broomrape/kmers/pass_threshold_10per_kmers_list_no_under_no_header.txt > /mnt/data/DanaS/kmer/broomrape/kmers/align_kmers/10per_kmers.fas

#-------------------------------------------------------------------------------------------------------------------------------------
# alinghment of the k-mers with bowyie 1 (For relatively short reads (e.g. less than 50 bp) Bowtie 1 is sometimes faster and/or more sensitive)

bowtie [options]* <ebwt> {-1 <m1> -2 <m2> | --12 <r> | --interleaved <i> | <s>} [<hit>]

./bowtie -a --best --strata
<s> A comma-separated list of files containing unpaired reads to be aligned


-a/--all           report all alignments per read
--best             hits guaranteed best stratum; ties broken by quality
--strata           hits in sub-optimal strata arent reported (requires --best)

# Specifying --strata in addition to -a and --best causes bowtie to report only those alignments in the best alignment “stratum”. 
# The alignments in the best stratum are those having the least number of mismatches (or mismatches just in the “seed” portion of 
# the alignment in the case of -n mode). Note that if --strata is specified, --best must also be specified.
#---------------------------------------------------------------------------------------------------------------------------------------

################# align with bwa (cause allredy has index of HA412) #################

bwa aln [-n maxDiff] [-o maxGapO] [-e maxGapE] [-d nDelTail] [-i nIndelEnd] [-k maxSeedDiff] [-l seedLen] 
[-t nThrds] [-cRN] [-M misMsc] [-O gapOsc] [-E gapEsc] [-q trimQual] <in.db.fasta> <in.query.fq> > <out.sai>


bwa aln [optiont] <index_prefix> <input_reads.fasta> > bwa_aln_alignments.sai

## HA412 index file psth/prefix

# run:
index_prefix=/mnt/data/sunflower/HA412/Ha412HOv2.0-20181130.genome.fasta
input_reads=/mnt/data/DanaS/kmer/broomrape/kmers/align_kmers/10per_kmers.fas

bwa aln -o 0 -n 0 $index_prefix $input_reads > /mnt/data/DanaS/kmer/broomrape/kmers/align_kmers/bwa_aln_10per_kmers_zero_dif.sai

# alignment output is .sai file, we need to make a .sam file from it  (using bwa samse)
in_sai=/mnt/data/DanaS/kmer/broomrape/kmers/align_kmers/bwa_aln_10per_kmers_zero_dif.sai

bwa samse $index_prefix $in_sai $input_reads > /mnt/data/DanaS/kmer/broomrape/kmers/align_kmers/bwa_aln_10per_kmers_zero_dif.sam


#To obtain more information from the alignment files we need to convert the sam format to bam using samtools

samtools view -bS /mnt/data/DanaS/kmer/broomrape/kmers/align_kmers/bwa_aln_10per_kmers_zero_dif.sam > /mnt/data/DanaS/kmer/broomrape/kmers/align_kmers/bwa_aln_10per_kmers_zero_dif.bam

#Then sort alignments by coordinates in the reference genome:

samtools sort /mnt/data/DanaS/kmer/broomrape/kmers/align_kmers/bwa_aln_10per_kmers_zero_dif.bam -o /mnt/data/DanaS/kmer/broomrape/kmers/align_kmers/bwa_aln_10per_kmers_zero_dif_sorted.bam

#And last, we need Samtools to index the BAM file:
samtools index /mnt/data/DanaS/kmer/broomrape/kmers/align_kmers/bwa_aln_10per_kmers_zero_dif_sorted.bam > /mnt/data/DanaS/kmer/broomrape/kmers/align_kmers/bwa_aln_10per_kmers_zero_dif_sorted_indx.bam

#Now that we’ve generated the files, we can view the output with tview
samtools tview /mnt/data/DanaS/kmer/broomrape/kmers/align_kmers/bwa_aln_10per_kmers_sorted.bam

samtools view /mnt/data/DanaS/kmer/broomrape/kmers/align_kmers/bwa_aln_10per_kmers_zero_dif_sorted_indx.bam Ha412HOChr14

# -----------------------------------------------------------------------------------
## see stuff in the sam /bam file:
cat /mnt/data/DanaS/kmer/broomrape/kmers/align_kmers/bwa_aln_10per_kmers.sam | less 
# see when X0 (Number of best hits) is 1
grep -c X0:i:1 /mnt/data/DanaS/kmer/broomrape/kmers/align_kmers/bwa_aln_10per_kmers.sam
#count uniqecaly mapped reads
grep -c XT:A:U /mnt/data/DanaS/kmer/broomrape/kmers/align_kmers/bwa_aln_10per_kmers_zero_dif.sam
#see uniqecaly mapped reads
grep -c XT:A:U /mnt/data/DanaS/kmer/broomrape/kmers/align_kmers/bwa_aln_10per_kmers_zero_dif.sam | less
# Count UNmapped reads
samtools view -f4 -c /mnt/data/DanaS/kmer/broomrape/kmers/align_kmers/bwa_aln_10per_kmers_zero_dif_sorted.bam
# Count the number of alignments (reads mapping to multiple locations counted multiple times):
samtools view -F 0x04 -c /mnt/data/DanaS/kmer/broomrape/kmers/align_kmers/bwa_aln_10per_kmers_zero_dif.sam
#Count number of records (unmapped reads + each aligned location per mapped read) in a bam file:
samtools view -c /mnt/data/DanaS/kmer/broomrape/kmers/align_kmers/bwa_aln_10per_kmers_zero_dif_sorted.bam
##to count alignments with mapping quality score >30
samtools view -q 20 -c /mnt/data/DanaS/kmer/broomrape/kmers/align_kmers/bwa_aln_10per_kmers_zero_dif_sorted.bam
#------------------------------------------------------

# get the reads that are mapped and have map quality > 20 in a new bam file 
samtools view -q 20 -b /mnt/data/DanaS/kmer/broomrape/kmers/align_kmers/bwa_aln_10per_kmers_zero_dif_sorted.bam > /mnt/data/DanaS/kmer/broomrape/kmers/align_kmers/aligned_kmers_over_q20.bam

# get chr names from the sorted bam with idxstats:
samtools idxstats /mnt/data/DanaS/kmer/broomrape/kmers/align_kmers/bwa_aln_10per_kmers_sorted.bam | head -n 17

#the chrom name are "Ha412HOChr01" "Ha412HOChr02" "Ha412HOChr03" exc..
samtools view -b $bamfile Ha412HOChr01 > Ha412HOChr01.bam
# now get the chr startpos end pos and the kmer in a tab file
bamfile=aligned_kmers_over_q20.bam
chrom_file=aligned_kmers_over_q20.tab

samtools view $bamfile | awk '{print $3 "\t" $4 "\t" $4+length($10)-1 "\t" $10}' > $chrom_file
#or by sequence
samtools view $bamfile | grep 'CGTAAGAGTCTTAGGCAAACATGGCCTTTGA' | awk '{print $3 "\t" $4 "\t" $4+length($10)-1}' > CGTAAGAGTCTTAGGCAAACATGGCCTTTGA.tab

samtools view -h /mnt/data/DanaS/kmer/broomrape/kmers/align_kmers/bwa_aln_10per_kmers_sorted.bam Ha412HOChr01 | head


conda deactivate
conda activate gemma098
--kmers_for_no_perm_phenotype

/usr/local/bin/gemma -bfile /mnt/data/DanaS/may_2024/yavor/by_n2r/mean/NEC/thresh_03/kmer/out/kmers/pheno.0.phenotype_value -lmm 2 -k /mnt/data/DanaS/may_2024/yavor/by_n2r/mean/NEC/thresh_03/kmer/out/pheno.kinship -outdir /mnt/data/DanaS/may_2024/yavor/by_n2r/mean/NEC/thresh_03/kmer/kmer_no_perm -o no_perm -maf 0.050000 -miss 0.5 

out_kmer=/mnt/data/DanaS/may_2024/yavor/by_n2r/mean/NEC/thresh_03/kmer_no_perm
pheno=/mnt/data/DanaS/may_2024/yavor/by_n2r/mean/NEC/thresh_03/n_NEC_n2r_mean_thresh03_yavor.txt
kmer_table=/mnt/data/sunflower/sunflowerKmer/kmers_table
#kmer_kinship=/mnt/data/sunflower/sunflowerKmer/kmers_table.kinship
# k-mer table path:
table=/mnt/data/DanaS/kmer/kmers_table
# path to a directory for output (dont use an existing dir i.e it will crearte a dir for the output)

screen -S kmers 
python2.7 /mnt/data/sunflower/sunflowerKmer/kmers_gwas.py --pheno $pheno --kmers_table $kmer_table -l 31 -p 20 --kmers_number 1000001 --outdir $out_kmer
python2.7 /mnt/data/DanaS/kmer/kmers_gwas.py --pheno $pheno --kmers_table $table --gemma_path /usr/local/bin/gemma -l 31 -p 20 --kmers_number 10000001 --permutations 20 --outdir $out_kmer

python2.7 $gwas_path --pheno $pheno --kmers_table $table --gemma_path /usr/local/bin/gemma -l 31 -p 20 --kmers_number 1000001 --outdir $out_kmer


#to plot a QQ plot you will need the accurate p-values for all variants, you can do this using the 
#kmers_table_to_bed functionality and then running it with an outside GWAS program (i.e GEMMA)

<KMERS-GWAS-PATH>/bin/kmers_table_to_bed -t kmers_table -k 31 -p phenotypes.pheno --maf 0.05 --mac 5 -b 10000000 -o output_file

# Parameters used:
# -t - k-mers table prefix path
# -k - length of k-mers
# -p - phenotype file, condense output only to individuals with a phenotype
# --maf - minor allele frequency
# --mac - minor allele count
# -b - maximal number of variants in each PLINK bed file (separate to many files)
# -o - prefix for output files
# -u - output only unique presence/absence patterns (optional)


# SAM003 SAM006 SAM085 and SAM272 are not in the kmer_table so need to delet them from pheno file
# in my file only sam003 and sam005 (in lines 4 and 6) apearing so remove them
sed '4d;6d' ../Total_No_of_attach_NEC_and_HEALTHY_norm2root.txt > pheno4qq.txt
pheno=/mnt/data/sunflower/sunflowerKmer/Branching_Gwas/pheno.phenotypes

/mnt/data/sunflower/sunflowerKmer/bin/kmers_table_to_bed -t /mnt/data/sunflower/sunflowerKmer/kmers_table \
-u -k 31 -p $pheno \
--maf 0.05 --mac 5 -b 500000000 -o table_to_bed_kmers

plink --bfile /mnt/data/DanaS/general_kmer_analysis/table_to_bed_kmers.0 --recode vcf --out vcf_kmers0
plink --bfile table_to_bed_kmers --freq counts
mv /mnt/data/DanaS/kmer/broomrape/kmers/QQ/pheno4qq.txt /mnt/data/DanaS/kmer/broomrape/kmers/qq/pheno4qq.txt

# kmer PAV
kmer_list_pval=/mnt/data/DanaS/may_2024/yavor/by_n2r/mean/NEC/thresh_025/kmers/get_kmers/aa_pass_bonfer_with_pval.txt
Rscript /mnt/data/DanaS/r_scripts/subset_kmer.r $kmer_list_pval
topkmers=/mnt/data/DanaS/may_2024/yavor/by_n2r/mean/NEC/thresh_03/kmer/get_kmers/top_100_kmers
table=/mnt/data/DanaS/kmer/kmers_table
kmers=/mnt/data/DanaS/may_2024/yavor/by_n2r/mean/NEC/thresh_025/kmers/get_kmers/top_100_kmers_with_pval
awk 'BEGIN{FS=OFS="\t"} NR>1{print $1}' $kmers > list100
kmer_list=/mnt/data/DanaS/may_2024/yavor/by_n2r/mean/NEC/thresh_025/kmers/get_kmers/top_100_kmers_list
kmer_list=/mnt/data/DanaS/may_2024/yavor/by_n2r/mean/NEC/thresh_03/kmer/get_kmers/kmers_list.txt
/mnt/data/DanaS/kmer/bin/filter_kmers -t $table --kmers_file $kmer_list -o output_pav.txt
kmer_list=/mnt/data/DanaS/sam_broom_gwas/gwas/gadot_nec_thresh_03/kmer/get_kmers/get_kmer/kmers_list.txt
pass_bonfer=/mnt/data/DanaS/may_2024/gadot/by_n2r/mean/NEC/thresh_03/kmers/n_NEC_n2r_max_zero_one/get_kmers/aa_pass_bonfer_with_pval.txt
pass_bonfer=/mnt/data/DanaS/may_2024/gadot/by_n2r/mean/NEC/thresh_01/kmers/processing_kmer/sig_kemers.txt
pass_10per=/mnt/data/DanaS/may_2024/gadot/by_n2r/mean/thresh_05/kmer/out/kmers/pass_threshold_10per
sed 's/_.*//' $pass_bonfer > kmers_list.txt
awk '{print $2, $10}' $pass_10per > kmers_list_pval.txt
awk 'NR>1{print $2, $10}' $pass_10per | sed 's/_.*//' > kmers_list.txt
kmer_list=/mnt/data/DanaS/sam_broom_gwas/gwas/yavor_nec_thresh_03/kmer/get_kmers/kmers_list.txt
kmer_list=/mnt/data/DanaS/sam_broom_gwas/gwas/yavor_hel_thresh_05/kmer/get_kmers/kmers_list.txt
/mnt/data/DanaS/kmer/bin/filter_kmers -t $table --kmers_file $kmer_list -o output_pav.txt

awk 'BEGIN{FS=OFS="\t"} NR>1{print $1}' /mnt/data/DanaS/may_2024/yavor/by_n2r/mean/kmer/thresh_05/n_HEALTHY_n2r_max_zero_one/processing_kmer/top_100_kmers > list100
pav=/mnt/data/DanaS/may_2024/gadot/by_n2r/mean/thresh_05/pav/output_pav.txt
pav=/mnt/data/DanaS/may_2024/gadot/by_n2r/mean/thresh_05/processing_kmer/output_pav.txt 
order=/mnt/data/DanaS/may_2024/gadot/pheno_means_gadot.txt
traits="n_HEALTHY_n2r_max"
traits="total_tubes_n2r_max"
broom=Gadot_HEL

pav=/mnt/data/DanaS/may_2024/gadot/by_n2r/mean/NEC/thresh_01/kmers/processing_kmer/output_pav.txt
pav=/mnt/data/DanaS/may_2024/gadot/by_n2r/mean/NEC/thresh_03/pav/output_pav.txt
pav=/mnt/data/DanaS/may_2024/gadot/by_n2r/mean/NEC/thresh_03/kmers/n_NEC_n2r_max_zero_one/get_kmers/output_pav.txt
order=/mnt/data/DanaS/may_2024/gadot/pheno_means_gadot.txt
order=/mnt/data/DanaS/may_2024/gadot/gadot_summary_all.txt
traits="n_NEC_n2r_max"
traits="prcnt_N_max"
broom=Gadot_NEC


pav=/mnt/data/DanaS/may_2024/yavor/by_n2r/mean/kmer/thresh_05/n_HEALTHY_n2r_max_zero_one/processing_kmer/output_pav.txt
order=/mnt/data/DanaS/may_2024/yavor/means_max_healthy_yavor.txt
traits="n_HEALTHY_n2r_max"
broom=CLASS_Yavor_HEL_

pav=/mnt/data/DanaS/may_2024/yavor/by_n2r/mean/NEC/thresh_025/kmers/get_kmers/output_pav100.txt
pav=/mnt/data/DanaS/may_2024/yavor/by_n2r/mean/NEC/thresh_025/kmers/processing_kmer/output_pav.txt
pav=/mnt/data/DanaS/may_2024/yavor/by_n2r/mean/NEC/thresh_03/kmer/get_kmers/output_pav.txt
order=/mnt/data/DanaS/may_2024/yavor/by_n2r/mean/NEC/means_max_nec_yavor_n2r.txt
order=/mnt/data/DanaS/may_2024/yavor/yavor_summary_all.txt
traits="pecnt_N_max"
traits="n_NEC_n2r_max"
broom=CLASS_Yavor_NEC
class="oil_NonOIL"

Class   class   class_others    HA_RHA_others   oil_NonOIL_others

Rscript /mnt/data/DanaS/r_scripts/pav_LD_plot.r $pav $order $traits $broom 70
Rscript /mnt/data/DanaS/r_scripts/pav_clust_by_infection_val.r $pav $order $traits $broom 50
Rscript /mnt/data/DanaS/r_scripts/ggalign.r $pav $order $traits $broom
Rscript /mnt/data/DanaS/r_scripts/pav_cluster_by_ld.r $pav $order $traits $broom 50
Rscript /mnt/data/DanaS/r_scripts/pav_cluster_by_ld_class.r $pav $order $traits $broom $class 50
Rscript /mnt/data/DanaS/r_scripts/pav_clust_by_hclust_kmers_then_kmeans_infection.r $pav $order $traits $broom $class 50
Rscript /mnt/data/DanaS/r_scripts/pav_clust_by_ld_kmer_infect_pop.r $pav $order $traits 

introg_list=/mnt/data/DanaS/sam_broom_gwas/gwas/haplotyps/introg_kmers.txt
/mnt/data/sunflower/sunflowerKmer/bin/filter_kmers \
-t /mnt/data/sunflower/sunflowerKmer/kmers_table \
-k $introg_list -o ./pav_table_introg.txt

introg_list=/mnt/data/DanaS/sam_broom_gwas/gwas/haplotyps/introg_kmer09s.txt
/mnt/data/sunflower/sunflowerKmer/bin/filter_kmers \
-t /mnt/data/sunflower/sunflowerKmer/kmers_table \
-k $introg_list -o ./pav_table_introg09.txt


introg_list="/mnt/data/DanaS/sam_broom_gwas/gwas/haplotyps/introg_kmers07_deb.txt"
/mnt/data/sunflower/sunflowerKmer/bin/filter_kmers \
-t /mnt/data/sunflower/sunflowerKmer/kmers_table \
-k $introg_list -o ./pav_table_introg07_deb.txt

pav=./pav_table_introg.txt
pav=./pav_table_introg07_deb.txt
order_full=/mnt/data/DanaS/may_2024/yavor/yavor_summary_all.txt
order_sub=/mnt/data/DanaS/sam_broom_gwas/gwas/yavor_hel_thresh_05/snp/filtered_by_n2r_HEALTHY_mean_thresh_0.5.txt
traits="n_HEALTHY_n2r_mean"
broom=Yavor_HealAD_ch07_deb
titele=HealAD
Rscript /mnt/data/DanaS/sam_broom_gwas/scripts/pav_clust_ld_new.r $pav $order_full $traits $broom $order_sub $titele

pav=./pav_table_introg09.txt
order_full=/mnt/data/DanaS/may_2024/gadot/gadot_summary_all.txt
order_sub=/mnt/data/DanaS/sam_broom_gwas/gwas/gadot_hel_thresh_05/snp/PHENO_FOR_GWAS.txt
traits="n_HEALTHY_n2r_mean"
broom=Gadot_HealAD_chr09
titele=HealAD




kmer_list=/mnt/data/DanaS/sam_broom_gwas/gwas/aeg/kmer/get_kmers/nec/kmers_list.txt
/mnt/data/sunflower/sunflowerKmer/bin/filter_kmers \
-t /mnt/data/sunflower/sunflowerKmer/kmers_table \
-k $kmer_list -o /mnt/data/DanaS/sam_broom_gwas/gwas/aeg/kmer/n_NEC_n2r_max/pav_table.txt

kmer_list=/mnt/data/DanaS/sam_broom_gwas/gwas/aeg/kmer/get_kmers/nec/kmers_list.txt
/mnt/data/sunflower/sunflowerKmer/bin/filter_kmers \
-t /mnt/data/sunflower/sunflowerKmer/kmers_table \
-k $kmer_list -o /mnt/data/DanaS/sam_broom_gwas/gwas/aeg/kmer/n_NEC_n2r_max/pav_table.txt

pav=./output_pav.txt
order_full=/mnt/data/DanaS/may_2024/gadot/gadot_summary_all.txt
order_sub=/mnt/data/DanaS/sam_broom_gwas/gwas/gadot_nec_thresh_03/kmer/n_NEC_n2r_mean_thresh03_gadot.txt
traits="n_NEC_n2r_max"
broom=Gadot_NecAD
titele=NecAD

pav=./output_pav.txt
order_full=/mnt/data/DanaS/may_2024/gadot/gadot_summary_all.txt
order_sub=/mnt/data/DanaS/sam_broom_gwas/gwas/gadot_hel_thresh_05/snp/PHENO_FOR_GWAS.txt
traits="n_HEALTHY_n2r_mean"
broom=Gadot_HealAD
titele=HealAD

pav=/mnt/data/DanaS/sam_broom_gwas/gwas/aeg/kmer/n_NEC_n2r_max/pav_table.txt
order_full=/mnt/data/DanaS/sam_broom_gwas/gwas/aeg/aeg_summary_all.txt
order_sub=/mnt/data/DanaS/sam_broom_gwas/gwas/aeg/aeg_summary_all.txt
traits="n_NEC_n2r_max"
broom=aeg_NecAD
titele=NecAD

Rscript /mnt/data/DanaS/sam_broom_gwas/scripts/pav_clust_ld_new.r $pav $order_full $traits $broom $order_sub $titele


pav=./output_pav.txt
order_full=/mnt/data/DanaS/may_2024/yavor/yavor_summary_all.txt
order_sub=/mnt/data/DanaS/sam_broom_gwas/gwas/yavor_hel_thresh_05/snp/filtered_by_n2r_HEALTHY_mean_thresh_0.5.txt
traits="n_HEALTHY_n2r_max"
broom=Yavor_HealAD_ch07
titele=HealAD

pav=./output_pav.txt
order_full=/mnt/data/DanaS/may_2024/yavor/yavor_summary_all.txt
order_sub=/mnt/data/DanaS/sam_broom_gwas/gwas/yavor_nec_thresh_03/snp/n_NEC_n2r_mean_thresh03_yavor.txt
traits="n_NEC_n2r_max"
broom=Yavor_NecAD
titele=NecAD

## pvast
cd /mnt/data/DanaS/sam_broom_gwas/gwas/yavor_nec_thresh_03/kmer/get_kmers
best_p=9.22374
cd /mnt/data/DanaS/sam_broom_gwas/gwas/yavor_hel_thresh_05/kmer/get_kmers
best_p=9.09233
for data in ./*/kmers_manhattan_plot_data.txt
do
    dir_path=$(dirname "$data")
    file="${data}"
    out_prefix="kmers_manhattan_plot_data_pvast.txt"
    out="${dir_path}/${out_prefix}"
    best_p=9.09233
    Rscript /mnt/data/DanaS/sam_broom_gwas/scripts/pvast.r $file $best_p $out
done
cd /mnt/data/DanaS/sam_broom_gwas/gwas/gadot_nec_thresh_03/kmer/get_kmers
Rscript /mnt/data/DanaS/sam_broom_gwas/scripts/pvast.r $file $best_p $out

pav=./output_pav.txt
order_full=/mnt/data/DanaS/may_2024/gadot/gadot_summary_all.txt
order_sub=/mnt/data/DanaS/sam_broom_gwas/gwas/gadot_hel_thresh_05/snp/PHENO_FOR_GWAS.txt
traits="n_HEALTHY_n2r_max"
broom=Gadot_HealAD
titele=HealADbroom=Yavor_HEL_thresh_05

order=/mnt/data/DanaS/may_2024/yavor/yavor_summary_all.txt
broom=Yavor_HEL_fullSAM
#prcnt_H_n2r_max
pav=pav_table.txt
traits="n_NEC_n2r_max"
titele=NEC-AD

order=/mnt/data/DanaS/sam_broom_gwas/gwas/yavor_nec_thresh_03/snp/n_NEC_n2r_mean_thresh03_yavor.txt
broom=Yavor_NEC_thresh_03

order=/mnt/data/DanaS/may_2024/yavor/yavor_summary_all.txt
broom=Yavor_NEC_fullSAM
#prcnt_N_n2r_max


cd /mnt/data/DanaS/sam_broom_gwas/gwas/gadot_hel_thresh_05/kmer/pav
pav="pav_table.txt"
traits="n_HEALTHY_n2r_mean"
titele="HEAL-AD"

order="/mnt/data/DanaS/may_2024/gadot/gadot_summary_all.txt"
broom="Gadot_HEL_fullSAM"

order="/mnt/data/DanaS/sam_broom_gwas/gwas/gadot_hel_thresh_05/snp/PHENO_FOR_GWAS.txt"
broom="Gadot_HEL_thresh_05"
#prcnt_H_n2r_max

pav=pav_table.txt
traits="n_NEC_n2r_mean"
titele=NEC-AD
order="/mnt/data/DanaS/may_2024/gadot/gadot_summary_all.txt"
broom=Gadot_NEC_fullSAM

order=/mnt/data/DanaS/sam_broom_gwas/gwas/gadot_nec_thresh_03/snp/n_NEC_n2r_mean_thresh03_gadot.txt
broom=Gadot_NEC_thresh_03
#prcnt_N_n2r_max
Rscript /mnt/data/DanaS/sam_broom_gwas/scripts/pav_cluster_by_ld.r $pav $order $traits $broom 50 $titele



order=/mnt/data/DanaS/sam_broom_gwas/gwas/yavor_nec_thresh_03/snp/n_NEC_n2r_mean_thresh03_yavor.txt
Rscript /mnt/data/DanaS/sam_broom_gwas/scripts/pav_anlysis.r $pav $order $traits $broom 50 $titele
Rscript /mnt/data/DanaS/sam_broom_gwas/scripts/pav_clust_ld_new.r $pav $order_full $traits $broom $order_sub $titele
pop="/mnt/data/DanaS/Pop_Code.txt"
grouping="Class"

for grouping in Class ha_rha oil_nonoil
do
    Rscript /mnt/data/DanaS/sam_broom_gwas/scripts/pav_candidate_rkmers.r $pav $order $traits $broom 50 $titele $pop $grouping
done


Rscript /mnt/data/DanaS/sam_broom_gwas/scripts/pav_candidate_rkmers.r $pav $order $traits $broom 50 $titele $pop $grouping


# Initialize analyzer
analyzer = PAVTraitAnalysis(threshold=0.05)

# Load and analyze data
pav_matrix, trait_values = analyzer.load_and_prepare_data(
    "pav_table.txt",
    "yavor_summary_all.txt",
    "n_HEALTHY_n2r_max"
)

# Run analysis
results = analyzer.analyze_kmers(pav_matrix, trait_values)
significant = analyzer.identify_significant_kmers(results)

# Create visualization
analyzer.plot_results(results, significant, "n_HEALTHY_n2r_max", "kmer_associations.png")


for class in  class_others HA_RHA_others oil_NonOIL_others 
do
    Rscript /mnt/data/DanaS/r_scripts/pav_cluster_by_ld_class.r $pav $order $traits $broom 50
done

for file in  do

Rscript /mnt/data/DanaS/r_scripts/pav_clust_by_infection_val.r $pav $order $traits $broom 50

for file in /mnt/data/DanaS/r_scripts/pav_clust_by*; do
    Rscript "$file" $pav $order $traits $broom 25
done

Rscript /mnt/data/DanaS/r_scripts/pav_clust_by_kmeans_combine_both_infection_and_kmer.r $pav $order $traits $broom 25



cp $chrom_file /mnt/data/DanaS/kmer/broomrape/kmers/man_p/aligned_kmers_over_q20.tab
#get rid of whats after the _ in the kmer pass_threshold_10per file that came from the assoc analysis
awk 'BEGIN{FS=OFS="\t"} {sub(/_.*/,"",$2)} 1' /mnt/data/DanaS/kmer/broomrape/kmers/pass_threshold_10per > /mnt/data/DanaS/kmer/broomrape/kmers/man_p/pass_threshold_10per_kmers_list.txt

#vlookup:
# take only sequence and p value col
awk 'BEGIN{FS=OFS="\t"; print "seq", "pvalue"} NR>1{print $2, $9}' pass_threshold_5per_kmers_list.txt > pval.txt
sed '1d' pval.txt > p_val.txt #remove header
file2=p_val.txt
file1=aligned_kmers_over_q20.tab

awk 'NR == FNR{a[$1] = $2;next}; {print $0, $4 in a?a[$4]: "NA"}' $file2 $file1 > file3.txt

#LD - to calculate LD between k-mers you just need the population presence/absence patterns of your associated k-mers. You can get 
#those by using the /bin/filter_kmers functionality. Then you can use any LD measure.
# You basically take your SNP of interest and extract their population pattern from the 1001G SNP matrix 
# and your k-mers presence/absence from the k-mers table and calculate LD between them.
