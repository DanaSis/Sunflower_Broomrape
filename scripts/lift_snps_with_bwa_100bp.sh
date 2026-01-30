#!/usr/bin/env bash
# Map ±50bp flanks (100bp total) with BWA-MEM and emit per-target lift.raw.tsv
# Run from the GWAS dir (where work_lift/ exists): 
#   /mnt/data/DanaS/sam_broom_gwas/scripts/lift_snps_with_bwa_100bp.sh
set -euo pipefail

# ---------- inputs ----------
SOURCE_FA_FAI="/mnt/data/DanaS/ANNOTATION_FILES/genomes/HA412/Ha412HOv2.0-20181130.genome.fasta.fai"
TARGETS_FAI=(
  "/mnt/data/sunflower/XRQ2/HanXRQr2.0-SUNRISE-2.1.genome.fasta.fai"
  "/mnt/data/DanaS/sam_broom_gwas/20211123_HanLR1r0.9-20211115_FunctionalAnnotation/fasta/HanLR1r0.9-20211115.genome.fasta.fai"
  "/mnt/data/sunflower/PSC8/HanPSC8r1.0-20181105.genome.fasta.fai"
)
# ----------------------------
# GEMMA assoc table
#ASSOC="/mnt/data/DanaS/sam_broom_gwas/gwas/yavor_hel_thresh_05/snp/out_gemma/n_HEALTHY_n2r_max_zero_one.assoc.txt"
#PVAL_THRESH="8.06047289182362e-08" 
#ASSOC="/mnt/data/DanaS/sam_broom_gwas/gwas/yavor_nec_thresh_03/snp/out_gemma/n_NEC_n2r_max_zero_one.assoc.txt"
#PVAL_THRESH="8.02458733559627e-07"  # not real thresh.. just to see top ones in chr10 and chr03
#ASSOC="/mnt/data/DanaS/sam_broom_gwas/gwas/gadot_hel_thresh_05/snp/out_gemma/n_HEALTHY_n2r_max_zero_one.assoc.txt"
#PVAL_THRESH="8.37993976499297e-08" 


ASSOC="/mnt/data/DanaS/sam_broom_gwas/gwas/aeg/thresh/nec/mean/thresh_001/out_gemma/n_NEC_n2r_mean_zero_one.assoc.txt"
PVAL_THRESH=7.788465282916e-08
# Choose p-value column & threshold (edit as needed)
# Common columns: p_wald, p_lrt, p_score
PVAL_COL="p_wald"
  # SimpleM threshold (edit if needed)

# Flank size (±bp around SNP on source for mapping)
FLANK=100

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
                  # 50 bp each side => 100 bp total

WORK="work_lift_bwa100"
OUTDIR="lifted_by_genome_bwa100"
ASSOC="work_lift/sig.tsv" 
mkdir -p "$WORK" "$OUTDIR"

need(){ command -v "$1" >/dev/null || { echo "ERROR: '$1' not in PATH" >&2; exit 1; }; }
need awk; need bedtools; need samtools; need bwa

SOURCE_FA="${SOURCE_FA_FAI%.fai}"
mapfile -t TARGETS < <(for t in "${TARGETS_FAI[@]}"; do echo "${t%.fai}"; done)

# 0) make sure we have sig.tsv (chr pos rs a1 a0 p)
[[ -s "$ASSOC" ]] || { echo "Missing $ASSOC (expected from your main run)"; exit 1; }

# 1) Build ±50 bp flanks FASTA (100 bp total), with A1/A0/RS in the header
#    Re-uses your established pattern and names.map trick
awk -v OFS="\t" '{print $1,$2-1,$2,$3,$4,$5}' "$ASSOC" > "$WORK/sig.1bp.meta.bed"
bedtools slop -i "$WORK/sig.1bp.meta.bed" -g <(cut -f1,2 "${SOURCE_FA}.fai") -b "$FLANK" > "$WORK/sig.slop50.bed"

# names like: chrom:start-end|A1=…|A0=…|RS=…
paste <(awk -v OFS="" '{print $1,":",$2,"-",$3}' "$WORK/sig.slop50.bed") \
      <(cut -f4 "$WORK/sig.slop50.bed") \
      <(cut -f5 "$WORK/sig.slop50.bed") \
      <(cut -f6 "$WORK/sig.slop50.bed") \
| awk -v OFS="|" '{print $1,"A1="$3,"A0="$4,"RS="$2}' > "$WORK/names.map"

bedtools getfasta -fi "$SOURCE_FA" -bed "$WORK/sig.slop50.bed" -name > "$WORK/sig.flanks50.fa"

# 2) For each target: BWA index (once), map flanks, and emit lift.raw.tsv
for TGT in "${TARGETS[@]}"; do
  TNAME=$(basename "$TGT" | sed 's/\.fa.*$//')
  echo ">> BWA lifting to $TNAME (±${FLANK}bp)"

  # index if missing
  if [[ ! -f "${TGT}.bwt" ]]; then
    echo "Indexing $TGT with BWA..."
    bwa index "$TGT"
  fi

  BAM="$WORK/${TNAME}.bwa.bam"
  # map: for 100bp queries, bwa mem is appropriate
  bwa mem -t "$THREADS" "$TGT" "$WORK/sig.flanks50.fa" \
    | samtools view -b -q "$MAPQ" \
    | samtools sort -o "$BAM"
  samtools index "$BAM"

  # BAM -> BED (qname preserved = header with A1/A0/RS)
  bedtools bamtobed -i "$BAM" > "$WORK/${TNAME}.bwa.bed"

  # Compute SNP 1-based POS from aligned flank and our FLANK size; parse RS/A1/A0 from name
  # BED: chrom start end name score strand
  awk -v OFS="\t" -v FL="$FLANK" '
    {
      rs="."; a1="."; a0=".";
      if (match($4,/RS=([^|]+)/,m)) rs=m[1];
      if (match($4,/A1=([^|]+)/,m)) a1=m[1];
      if (match($4,/A0=([^|]+)/,m)) a0=m[1];
      strand=$6;
      pos = (strand=="+") ? $2+FL+1 : $3-FL;  # center = FL + 1 from left end on +
      print $1, pos, rs, a1, a0, strand
    }' "$WORK/${TNAME}.bwa.bed" > "$WORK/${TNAME}.lift.raw.tsv"

  # Summarize quick
  echo "Wrote raw lifts: $WORK/${TNAME}.lift.raw.tsv"
  echo "Example:"
  head -n3 "$WORK/${TNAME}.lift.raw.tsv" || true

  # Optionally copy to OUTDIR with a clear name
  cp "$WORK/${TNAME}.lift.raw.tsv" "$OUTDIR/${TNAME}.bwa100.lift.raw.tsv"
done

echo "Done. Per-target raw lifts in: $OUTDIR/"

