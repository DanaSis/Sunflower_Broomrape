#!/bin/bash

# --- 1. Configuration ---

# Define the genomes and their index paths in an associative array
declare -A GENOMES=(
    ["ha412"]="/mnt/data/sunflower/HA412/Ha412HOv2.0-20181130.genome.fasta"
    ["xrq"]="/mnt/data/sunflower/XRQ2/HanXRQr2.0-SUNRISE-2.1.genome.fasta"
    ["psc8"]="/mnt/data/sunflower/PSC8/HanPSC8r1.0-20181105.genome.fasta"
    ["lr1"]="/mnt/data/DanaS/sam_broom_gwas/20211123_HanLR1r0.9-20211115_FunctionalAnnotation/fasta/HanLR1r0.9-20211115.genome.fasta"
)

# Define the traits and base directory
TRAITS=("n_HEL_n2r_max" "n_NEC_n2r_max" "total_n2r_max")
BASE_DIR="/mnt/data/DanaS/sam_broom_gwas/gwas/aeg/kmer"

# Define the thresholds to process
THRESHOLDS=("10per" "5per")

# R script location
RSCRIPT="/mnt/data/DanaS/r_scripts/match_blast.r"

# --- 2. Function to Process Each Genome ---
# This function encapsulates the entire mapping and analysis pipeline
process_genome() {
    local genome_id="$1"
    local index_prefix="$2"
    local input_kmers_fas="$3"
    local input_gwas_data="$4"
    local output_dir="$5" # The directory: .../map_kmers/$threshold
    
    # Create and move into the genome-specific directory within the output path
    mkdir -p "${output_dir}/${genome_id}"
    cd "${output_dir}/${genome_id}"
    
    echo -e "\n### Starting analysis for ${genome_id} in ${output_dir}/${genome_id} ###"
    
    # Define files with genome prefix for clarity and uniqueness
    local genome_prefix="${genome_id}_"
    local chrom_file="${genome_prefix}aligned_kmers_over_q20.tab"
    local for_match_blast="for_match_blast.txt"

    # Alignment and processing steps (using BWA, samtools)
    # The commands below are the standard steps for k-mer mapping
    bwa aln -o 0 -n 0 "$index_prefix" "$input_kmers_fas" > "${genome_prefix}bwa_aln_kmers_0_mismatc.sai"
    bwa samse "$index_prefix" "${genome_prefix}bwa_aln_kmers_0_mismatc.sai" "$input_kmers_fas" > "${genome_prefix}bwa_aln_kmers_0_mismatc.sam"
    samtools view -bS "${genome_prefix}bwa_aln_kmers_0_mismatc.sam" > "${genome_prefix}bwa_aln_kmers_0_mismatc.bam"
    samtools sort "${genome_prefix}bwa_aln_kmers_0_mismatc.bam" -o "${genome_prefix}bwa_aln_kmers_0_mismatc_sorted.bam"
    samtools index "${genome_prefix}bwa_aln_kmers_0_mismatc_sorted.bam"
    samtools view -q 20 -b "${genome_prefix}bwa_aln_kmers_0_mismatc_sorted.bam" > "${genome_prefix}aligned_kmers_over_q20.bam"
    
    # Extract coordinates (Chr, Start, End, Kmer sequence)
    echo "1. Extracting aligned coordinates (Q>20)..."
    samtools view "${genome_prefix}aligned_kmers_over_q20.bam" | \
        awk '{print $3 "\t" $4 "\t" $4+length($10)-1 "\t" $10}' > "$chrom_file"
    
    # Prepare for BLAST
    echo "2. Running BLAST comparison..."
    local gwas_kmers_fasta="gwas_kmers_with_ids.fasta"
    local fasta_align_over_q20="over_q20.fas"

    # Prepare GWAS data (ID suffix for BLAST header)
    awk 'NR>1 {split($1, id, "_"); print ">"id[2]"\n"substr($1, 1, length($1)-length(id[2])-1)}' "$input_gwas_data" > "$gwas_kmers_fasta"
    # Prepare Aligned data (for BLAST subject)
    awk '{ print ">"NR"\n"$4 }' "$chrom_file" > "$fasta_align_over_q20"
    
    # Run BLAST
    blastn -query "$gwas_kmers_fasta" \
        -subject "$fasta_align_over_q20" \
        -outfmt "6 qseqid qseq sseqid sseq pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
        -out blast_results.txt -strand both
    
    # Get only the first hit for each query kmer and add headers
    awk 'BEGIN{FS=OFS="\t"; print "qseqid", "qseq", "sseqid", "sseq", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"} {print}' blast_results.txt | \
        awk '{FS=OFS="\t"} !seen[$1]++' > blast_results_with_seq_h_firsts.txt
    
    # Run R script for final table
    echo "3. Running R script for final table..."
    
    # Prepare GWAS data (kmer and p-value) with correct headers for R script
    awk 'BEGIN{FS=OFS="\t"; print "rs", "p_lrt"} {print $1, $2}' "$input_gwas_data" > "$for_match_blast"

    Rscript "$RSCRIPT" "$for_match_blast" "blast_results_with_seq_h_firsts.txt" "$chrom_file"

    # Move back to the trait directory (2 levels up)
    cd ../.. 
}

# --- 3. Main Pipeline Execution ---

echo "Starting Kmer Mapping Pipeline for AEG Traits..."

# Loop 1: Traits
for trait in "${TRAITS[@]}"
do
    TRAIT_PATH="${BASE_DIR}/${trait}"
    
    if [ ! -d "$TRAIT_PATH" ]; then
        echo "Error: Trait directory not found: $TRAIT_PATH"
        continue
    fi

    echo -e "\n========================================================"
    echo "STARTING TRAIT: $trait"
    echo "========================================================"
    
    # Loop 2: Thresholds (5per and 10per)
    for threshold in "${THRESHOLDS[@]}"
    do
        PASS_FILE="${TRAIT_PATH}/out/kmers/pass_threshold_${threshold}"
        
        if [ ! -f "$PASS_FILE" ]; then
            echo "Skipping threshold ${threshold} for $trait: Input file not found: $PASS_FILE"
            continue
        fi

        # Define the output directory path for the current threshold
        MAP_BASE_DIR="${TRAIT_PATH}/map_kmers"
        OUTPUT_DIR="${MAP_BASE_DIR}/${threshold}"
        
        echo -e "\n--- Processing $trait/$threshold ---"
        
        # --- Pre-processing steps (Kmer extraction and FASTA creation) ---
        
        # 1. Extract Kmer and P-value into a temporary file
        TEMP_GWAS_DATA="${TRAIT_PATH}/gwas_data_${threshold}_temp.txt"
        awk 'NR>1{ print $2, $10 }' "$PASS_FILE" > "$TEMP_GWAS_DATA"
        
        # 2. Extract Kmers list and convert to FASTA format for BWA
        KMERS_LIST="${TRAIT_PATH}/kmers_list_${threshold}.txt"
        KMERS_FASTA="${TRAIT_PATH}/${threshold}_kmers.fas"
        
        awk '{ print $1 }' "$TEMP_GWAS_DATA" | sed 's/_.*//' > "$KMERS_LIST"
        awk '{ print ">"NR"\n"$0 }' "$KMERS_LIST" > "$KMERS_FASTA"
        
        # Loop 3: Genomes
        for genome_id in "${!GENOMES[@]}"
        do
            index_path="${GENOMES[${genome_id}]}"
            
            # Execute the main analysis function
            process_genome "$genome_id" "$index_path" "$KMERS_FASTA" "$TEMP_GWAS_DATA" "$OUTPUT_DIR"
        done

        # Clean up temporary files created for the current threshold
        rm "$TEMP_GWAS_DATA" "$KMERS_LIST" "$KMERS_FASTA"

    done # End Threshold Loop

done # End Trait Loop

echo -e "\nPipeline finished for AEG traits. Results are in $BASE_DIR/*/map_kmers/*"