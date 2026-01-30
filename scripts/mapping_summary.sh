#!/bin/bash

# --- Configuration ---
YEARS="2021 2023"
THRESHOLDS="5per 10per"
GENOMES="ha412 xrq psc8 lr1"
BASE_PATH="/mnt/data/DanaS/seeds/HA_only/"
INT_BASE_PATH="/mnt/data/DanaS/ANNOTATION_FILES/introgressions/"
EXTEND_BP=250000 
LOG_FILE="summary_introgression_ancestry_log.txt"

# Define the output file for the summary table
SUMMARY_FILE="kmer_mapping_summary_with_introgression_$(date +%Y%m%d_%H%M%S).txt"

# --- Function to Check Introgression Status (CONFIRMED) ---
check_introgression_status() {
    local genome_id="$1"
    local manhattan_file="$2"
    local run_identifier="$3" 
    local int_file="${INT_BASE_PATH}${genome_id}_introgression.txt"
    
    if [ ! -f "$manhattan_file" ] || [ ! -f "$int_file" ]; then
        echo "0\tNA"  
        return
    fi
    
    # 1. Prepare Kmer Mapped Regions (BED format - Zero-based)
    awk 'BEGIN{FS=OFS="\t"} 
        NR > 1 { 
            # pos ($2) is 1-based. BED Start is 0-based.
            start = ($2 - '"$EXTEND_BP"' - 1) > 0 ? ($2 - '"$EXTEND_BP"' - 1) : 0;
            end = $2 + '"$EXTEND_BP"'; 
            print $1, start, end, $2, $3  # Field 4 ($2) is the original position
        }' "$manhattan_file" > kmer_regions.bed
        
    # --- Function to Check Introgression Status (FINAL FINAL REVISION) ---
check_introgression_status() {
    local genome_id="$1"
    local manhattan_file="$2"
    local run_identifier="$3" 
    local int_file="${INT_BASE_PATH}${genome_id}_introgression.txt"
    
    if [ ! -f "$manhattan_file" ] || [ ! -f "$int_file" ]; then
        echo "0\tNA"  
        return
    fi
    
    # ... (Kmer Mapped Regions preparation remains unchanged and is correct) ...
    awk 'BEGIN{FS=OFS="\t"} 
        NR > 1 { 
            start = ($2 - '"$EXTEND_BP"' - 1) > 0 ? ($2 - '"$EXTEND_BP"' - 1) : 0;
            end = $2 + '"$EXTEND_BP"'; 
            print $1, start, end, $2, $3
        }' "$manhattan_file" > kmer_regions.bed
        
    # 2. Prepare Introgression Data (BED format - Zero-based, FIXED ANCESTRY)
    # The tr -d ',' removes commas.
    # The subsequent awk command prints columns 1, 2-1, 3, AND everything from column 4 onwards ($4 to NF)
    tail -n +2 "$int_file" | tr -d ',' | \
        awk 'BEGIN{FS=OFS="\t"} {
            # $1, $2, $3 are Chromosome, Start, End. 
            # The rest of the line, starting from the 4th field, is stored in a variable:
            full_ancestry = $4; 
            for (i=5; i<=NF; i++) full_ancestry = full_ancestry " " $i;
            
            # Print Chromosome, Start-1 (0-based), End, and the reconstructed Ancestry string
            print $1, $2 - 1, $3, full_ancestry 
        }' > introgression.bed
        
    # ... (Rest of the function logic remains unchanged and is now correct) ...
    
    # 3. Use BEDTools intersect to find overlaps
    OVERLAP_DATA=$(bedtools intersect -a kmer_regions.bed -b introgression.bed -wao | \
                   awk '$NF > 0')

    # Count unique Kmer positions ($4) that overlap
    OVERLAP_COUNT=$(echo "$OVERLAP_DATA" | awk 'NF > 0 { print $4 }' | sort -u | wc -l)
    
    # 4. Get unique Ancestries found (Ancestor is now field $8 in the intersect output)
    if [ "$OVERLAP_COUNT" -gt 0 ]; then
        UNIQUE_ANCESTRY=$(echo "$OVERLAP_DATA" | awk 'NF > 0 { print $8 }' | sort -u | tr '\n' ',' | sed 's/,$//')

        # --- Detailed Log ---
        echo -e "\n### ${run_identifier} / ${genome_id} ###" >> "$LOG_FILE"
        echo "Ancestries found: ${UNIQUE_ANCESTRY}" >> "$LOG_FILE"
        
        # Report back the count and the unique ancestry list
        echo "${OVERLAP_COUNT}\t${UNIQUE_ANCESTRY}"
    else
        echo "0\tNone"
    fi
    
    # Clean up temporary files
    rm kmer_regions.bed introgression.bed
}
}

# --- Initialization ---
echo -e "--- Kmer Introgression Ancestry Log ---" > "$LOG_FILE"

# Print header for the summary table
echo -e "Year\tTrait\tThreshold\tTotal_Passing_Kmers\tGenome\tMapped_Kmers_Q20\tManhattan_Plot_Rows\tKmers_in_Introgression_Region\tIntrogression_Ancestry" > "$SUMMARY_FILE"

# --- Main Analysis Loop ---
for year in $YEARS; do
    YEAR_DIR="${BASE_PATH}${year}"
    for trait_path in "$YEAR_DIR"/*/ ; do
        if [ ! -d "$trait_path" ]; then continue; fi
        trait=$(basename "$trait_path")
        KMER_INPUT_DIR="${trait_path}out/kmers"
        MAP_BASE_DIR="${trait_path}mapp_kmers"
        if [ ! -d "$MAP_BASE_DIR" ]; then continue; fi
        
        for threshold in $THRESHOLDS; do
            PASS_FILE="${KMER_INPUT_DIR}/pass_threshold_${threshold}"
            OUTPUT_DIR="${MAP_BASE_DIR}/${threshold}"
            
            # Get Initial Kmer Count
            if [ -f "$PASS_FILE" ]; then TOTAL_KMERS=$(($(wc -l < "$PASS_FILE") - 1)); else TOTAL_KMERS=0; fi
            if [ ! -d "$OUTPUT_DIR" ]; then continue; fi

            for genome_id in $GENOMES; do
                GENOME_DIR="${OUTPUT_DIR}/${genome_id}"
                MANHATTAN_FILE="${GENOME_DIR}/kmers_manhattan_plot_data.txt"
                
                # Check for file existence/line counts
                Q20_TAB_FILE="${GENOME_DIR}/${genome_id}_aligned_kmers_over_q20.tab"
                if [ -f "$Q20_TAB_FILE" ]; then MAPPED_Q20=$(wc -l < "$Q20_TAB_FILE"); else MAPPED_Q20=0; fi
                if [ -f "$MANHATTAN_FILE" ]; then MANHATTAN_ROWS=$(($(wc -l < "$MANHATTAN_FILE") - 1)); else MANHATTAN_ROWS=0; fi

                # Check Introgression & Ancestry
                RUN_ID="${year}/${trait}/${threshold}"
                if [ "$MANHATTAN_ROWS" -gt 0 ]; then
                    # Call function and capture two outputs: Count and Ancestry list
                    INT_RESULT=$(check_introgression_status "$genome_id" "$MANHATTAN_FILE" "$RUN_ID")
                    INT_COUNT=$(echo "$INT_RESULT" | awk '{print $1}')
                    ANCESTRY_LIST=$(echo "$INT_RESULT" | awk '{print $2}')
                else
                    INT_COUNT=0
                    ANCESTRY_LIST="None"
                fi
                
                # Print the results row
                echo -e "${year}\t${trait}\t${threshold}\t${TOTAL_KMERS}\t${genome_id}\t${MAPPED_Q20}\t${MANHATTAN_ROWS}\t${INT_COUNT}\t${ANCESTRY_LIST}" >> "$SUMMARY_FILE"
            done
        done
    done
done

echo -e "\nSummary table generated in ${SUMMARY_FILE}"
echo -e "Detailed ancestry information is in ${LOG_FILE}"