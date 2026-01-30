#!/usr/bin/env bash
# RUN FROM the GWAS dir (e.g. /mnt/data/DanaS/sam_broom_gwas/gwas/gadot_hel_thresh_05)
# chmod +x /mnt/data/DanaS/sam_broom_gwas/scripts/lift_snps_between_genomes.sh
# /mnt/data/DanaS/sam_broom_gwas/scripts/lift_snps_between_genomes.sh
set -euo pipefail

# ------------------------ USER INPUTS ------------------------
# Source (Ha412HOv2) FASTA .fai and targets’ .fai — ".fai" will be stripped automatically
SOURCE_FA_FAI="/mnt/data/DanaS/ANNOTATION_FILES/genomes/HA412/Ha412HOv2.0-20181130.genome.fasta.fai"
TARGETS_FAI=(
  "/mnt/data/sunflower/XRQ2/HanXRQr2.0-SUNRISE-2.1.genome.fasta.fai"
  "/mnt/data/DanaS/sam_broom_gwas/20211123_HanLR1r0.9-20211115_FunctionalAnnotation/fasta/HanLR1r0.9-20211115.genome.fasta.fai"
  "/mnt/data/sunflower/PSC8/HanPSC8r1.0-20181105.genome.fasta.fai"
)

# GEMMA assoc table
#ASSOC="/mnt/data/DanaS/may_2024/yavor/by_n2r/mean/NEC/thresh_03/out_gemma/n_NEC_n2r_mean_zero_one.assoc.txt"
ASSOC="/mnt/data/DanaS/sam_broom_gwas/gwas/aeg/thresh/nec/mean/thresh_001/out_gemma/n_NEC_n2r_mean_zero_one.assoc.txt"
# Choose p-value column & threshold (edit as needed)
# Common columns: p_wald, p_lrt, p_score
PVAL_COL="p_wald"
#PVAL_THRESH="8.02458733559627e-07"   # SimpleM threshold (edit if needed)
PVAL_THRESH=7.788465282916e-08
# Flank size (±bp around SNP on source for mapping)
FLANK=50

# Mapping/filter
THREADS=8
MAPQ=20

# Work/output dirs (created under current directory)
WORK="work_lift"
OUTDIR="lifted_by_genome"
# -------------------------------------------------------------

# ------------------------ PRECHECKS --------------------------
need() { command -v "$1" >/dev/null 2>&1 || { echo "ERROR: '$1' not found in PATH" >&2; exit 1; }; }
need awk; need samtools; need bedtools; need minimap2; need bcftools

mkdir -p "$WORK" "$OUTDIR"

SOURCE_FA="${SOURCE_FA_FAI%.fai}"
mapfile -t TARGETS < <(for t in "${TARGETS_FAI[@]}"; do echo "${t%.fai}"; done)

[[ -f "$SOURCE_FA" ]] || { echo "ERROR: SOURCE fasta not found: $SOURCE_FA" >&2; exit 1; }
[[ -f "${SOURCE_FA}.fai" ]] || samtools faidx "$SOURCE_FA"
for TGT in "${TARGETS[@]}"; do
  [[ -f "$TGT" ]] || { echo "ERROR: TARGET fasta not found: $TGT" >&2; exit 1; }
  [[ -f "${TGT}.fai" ]] || samtools faidx "$TGT"
done
[[ -f "$ASSOC" ]] || { echo "ERROR: assoc file not found: $ASSOC" >&2; exit 1; }

# ----------------- 1) pick significant SNPs -----------------
# Keep: chr, ps(1-based), rs, allele1, allele0, pval
awk -v col="$PVAL_COL" -v thr="$PVAL_THRESH" '
BEGIN{FS=OFS="\t"}
NR==1{
  for(i=1;i<=NF;i++) h[$i]=i;
  if(!(col in h)){print "ERROR: column "col" not found in assoc file." > "/dev/stderr"; exit 1}
  next
}
{
  p=$h[col]+0;
  if (p <= thr) print $h["chr"], $h["ps"], $h["rs"], $h["allele1"], $h["allele0"], p
}' "$ASSOC" > "$WORK/sig.tsv"

if [[ ! -s "$WORK/sig.tsv" ]]; then
  echo "No SNPs pass ${PVAL_COL} <= ${PVAL_THRESH}. Nothing to lift."
  exit 0
fi

# ----------------- 2) flanks on source ----------------------
# Build slopped BED with meta columns: chr start end rs a1 a0
awk -v OFS="\t" '{print $1,$2-1,$2,$3,$4,$5}' "$WORK/sig.tsv" > "$WORK/sig.1bp.meta.bed"
bedtools slop -i "$WORK/sig.1bp.meta.bed" -g <(cut -f1,2 "${SOURCE_FA}.fai") -b "$FLANK" > "$WORK/sig.slop.bed"

# Make a named BED where col4 is a rich header that minimap2 carries as QNAME
# qname format: CHR:START-END|A1=..|A0=..|RS=..
awk 'BEGIN{OFS="\t"}{
  chr=$1; st=$2; en=$3; rs=$4; a1=$5; a0=$6;
  name=chr ":" st "-" en "|A1=" a1 "|A0=" a0 "|RS=" rs;
  print chr, st, en, name
}' "$WORK/sig.slop.bed" > "$WORK/sig.slop.named.bed"

# Extract flanks with the rich QNAME
bedtools getfasta -fi "$SOURCE_FA" -bed "$WORK/sig.slop.named.bed" -name > "$WORK/sig.flanks.fa"

# ----------------- 3) map to each target & emit VCF ---------
for TGT in "${TARGETS[@]}"; do
  TNAME=$(basename "$TGT" | sed 's/\.fa.*$//')
  echo ">> Lifting to $TNAME"

  # Minimizer index (once)
  MMI="$WORK/${TNAME}.mmi"
  [[ -f "$MMI" ]] || minimap2 -d "$MMI" "$TGT"

  # Align flanks; keep unique/confident mappings
  BAM="$WORK/${TNAME}.bam"
  minimap2 -x sr -t "$THREADS" -a "$MMI" "$WORK/sig.flanks.fa" \
    | samtools view -b -q "$MAPQ" \
    | samtools sort -o "$BAM"
  samtools index "$BAM"

  # Convert to BED6: tChr tStart tEnd qName score strand
  bedtools bamtobed -i "$BAM" > "$WORK/${TNAME}.bed"

  # Optional: warn if not all flanks mapped uniquely
  n_reads=$(grep -c '^>' "$WORK/sig.flanks.fa" || true)
  n_aln=$(wc -l < "$WORK/${TNAME}.bed")
  if [ "$n_aln" -lt "$n_reads" ]; then
    echo "Note: $((n_reads - n_aln)) of $n_reads flanks did not align uniquely to $TNAME (MAPQ>=$MAPQ)." >&2
  fi

  # Parse A1/A0/RS from qName and compute target SNP position (offset = FLANK)
  awk -v OFS="\t" -v FL="$FLANK" '
  {
    # qName in $4 is "CHR:START-END|A1=..|A0=..|RS=.."
    split($4, P, /\|/);         # P[1]=CHR:START-END, P[2]=A1=.., P[3]=A0=.., P[4]=RS=..
    a1=P[2]; sub("^A1=","",a1)
    a0=P[3]; sub("^A0=","",a0)
    rs=P[4]; sub("^RS=","",rs)

    if ($6 == "+") pos0 = $2 + FL; else pos0 = $3 - FL - 1;
    pos1 = pos0 + 1
    print $1, pos1, rs, a1, a0, $6
  }' "$WORK/${TNAME}.bed" > "$WORK/${TNAME}.lift.raw.tsv"
  # Columns: TCHR TPOS RS A1 A0 STRAND

  # Fetch target REF base safely; determine REF/ALT, mark mismatches
  awk -v OFS="\t" '{print $1,$2,$3,$4,$5,$6}' "$WORK/${TNAME}.lift.raw.tsv" \
  | while IFS=$'\t' read -r tchr tpos rs a1 a0 strand; do
      tref=$(samtools faidx "$TGT" "${tchr}:${tpos}-${tpos}" | awk 'NR>1{printf toupper($0)}')
      A1=$(echo "$a1" | tr a-z A-Z); A0=$(echo "$a0" | tr a-z A-Z)
      status="OK"
      if   [[ "$A1" != "$tref" && "$A0" == "$tref" ]]; then ref="$tref"; alt="$A1"
      elif [[ "$A0" != "$tref" && "$A1" == "$tref" ]]; then ref="$tref"; alt="$A0"
      elif [[ "$A0" != "$tref" && "$A1" != "$tref" ]]; then status="MISMATCH"; ref="$tref"; alt="$A1"
      else status="TIE"; ref="$tref"; alt="$A1"; fi
      printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$tchr" "$tpos" "$rs" "$ref" "$alt" "$strand" "$status"
    done > "$WORK/${TNAME}.lift.checked.tsv"

  # Emit VCF, dropping MISMATCH/TIE rows (they remain in lift.checked.tsv for QC)
  VCF="$OUTDIR/${TNAME}.lifted.vcf.gz"
  {
    echo "##fileformat=VCFv4.2"
    echo "##source=flanklift_minimap2"
    echo "##reference=$TGT"
    echo '##INFO=<ID=LIFT_STRAND,Number=1,Type=String,Description="Strand of flank alignment on target (+ or -)">'
    echo '##INFO=<ID=STATUS,Number=1,Type=String,Description="Allele check: OK|MISMATCH|TIE">'
    echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
    awk -v OFS="\t" '$7!="MISMATCH" && $7!="TIE" {print $1,$2,$3,$4,$5,".","PASS","LIFT_STRAND="$6";STATUS="$7}' "$WORK/${TNAME}.lift.checked.tsv"
  } | bgzip > "$VCF"
  bcftools index -f "$VCF"

  echo "Summary ($TNAME):"
  awk '{c[$7]++} END{for(k in c) printf "%-12s\t%d\n", k, c[k]}' "$WORK/${TNAME}.lift.checked.tsv" | sort
  echo "VCF => $VCF"
  echo
done

echo "Done. Per-target VCFs are in: $OUTDIR/"





python3 /mnt/data/DanaS/sam_broom_gwas/scripts/merge_plot_chr9_all_markers.py \
  --snps /mnt/data/DanaS/sam_broom_gwas/gwas/gadot_hel_thresh_05/snp/snps_chr9.raw \
  --pav  /mnt/data/DanaS/sam_broom_gwas/gwas/pav_chr9.kmerXsample.tsv \
  --pheno /mnt/data/DanaS/sam_broom_gwas/gwas/gadot_hel_thresh_05/snp/PHENO_FOR_GWAS.txt \
  --sample-col sample \
  --pheno-col NecAD_mean \
  --out /mnt/data/DanaS/sam_broom_gwas/gwas/gadot_hel_thresh_05/plots/chr9_all_markers.NecAD_heatmap.png \
  --max-snps 400 \
  --max-kmers 400


