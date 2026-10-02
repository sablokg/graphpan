#!/usr/bin/env bash
###############################################################################
# GRAPH PANGENOME PIPELINE — build -> normalize -> variant call -> annotate
#
# Stack: pggb (graph construction) + odgi (graph manipulation/stats)
#        + vg (variant extraction) + bcftools + snpEff/VEP (annotation)
#
# Alternative graph builders (swap Stage 2 if preferred):
#   - minigraph-cactus  (better for large, divergent genomes / SVs)
#   - minigraph alone   (reference-based, faster, less complete graph)
#
# USAGE:
#   ./pangenome_pipeline.sh -i genomes/ -r reference.fa -o results -g annotation.gff3 -t 16
###############################################################################

set -euo pipefail

# ---------------------------- 0. ARGUMENTS ----------------------------------
usage() {
  echo "Usage: $0 -i <input_genomes_dir> -r <reference.fa> -o <outdir> -g <ref_annotation.gff3> [-t threads] [-p ploidy]"
  exit 1
}

THREADS=8
PLOIDY=2

while getopts "i:r:o:g:t:p:h" opt; do
  case $opt in
    i) INDIR=$OPTARG ;;
    r) REF=$OPTARG ;;
    o) OUTDIR=$OPTARG ;;
    g) REF_GFF=$OPTARG ;;
    t) THREADS=$OPTARG ;;
    p) PLOIDY=$OPTARG ;;
    h) usage ;;
    *) usage ;;
  esac
done

[[ -z "${INDIR:-}" || -z "${REF:-}" || -z "${OUTDIR:-}" || -z "${REF_GFF:-}" ]] && usage

mkdir -p "$OUTDIR"/{00_input,01_graph,02_odgi,03_variants,04_annotation,05_stats,logs}

echo "[INFO] Genome input dir : $INDIR"
echo "[INFO] Reference genome : $REF"
echo "[INFO] Reference GFF3   : $REF_GFF"
echo "[INFO] Output dir       : $OUTDIR"
echo "[INFO] Threads          : $THREADS"

# ---------------------------- 1. ENVIRONMENT ---------------------------------
# Recommended: create a dedicated conda/mamba env with pinned tool versions.
# Uncomment to auto-create it.
#
# mamba create -n pangenome -c bioconda -c conda-forge \
#     pggb odgi vg bcftools samtools vcfbub snpeff seqkit panacus -y
# conda activate pangenome

command -v pggb   >/dev/null || { echo "[ERROR] pggb not found in PATH"; exit 1; }
command -v odgi   >/dev/null || { echo "[ERROR] odgi not found in PATH"; exit 1; }
command -v vg     >/dev/null || { echo "[ERROR] vg not found in PATH"; exit 1; }
command -v bcftools >/dev/null || { echo "[ERROR] bcftools not found in PATH"; exit 1; }

# ---------------------------- 2. PREPARE INPUT --------------------------------
# All input genomes + reference combined into a single FASTA, indexed.
# PGGB requires all sequences merged, bgzip-compressed and samtools-faidx'd.

echo "[STAGE 1] Merging & indexing input genomes"

MERGED="$OUTDIR/00_input/all_genomes.fa"
cat "$REF" "$INDIR"/*.fa "$INDIR"/*.fasta 2>/dev/null > "$MERGED" || true

# Deduplicate sequence names if needed (PanSN-spec naming is strongly recommended:
#   sample#haplotype#contig  e.g.  HG002#1#chr1
# Rename here if your FASTAs don't already follow this convention.
seqkit rename "$MERGED" > "$OUTDIR/00_input/all_genomes.renamed.fa" 2>/dev/null || \
  cp "$MERGED" "$OUTDIR/00_input/all_genomes.renamed.fa"

bgzip -f "$OUTDIR/00_input/all_genomes.renamed.fa"
samtools faidx "$OUTDIR/00_input/all_genomes.renamed.fa.gz"

INPUT_FA="$OUTDIR/00_input/all_genomes.renamed.fa.gz"
N_SEQS=$(grep -c '^>' <(zcat "$INPUT_FA") || true)
echo "[INFO] $N_SEQS sequences merged into pangenome input"

# ---------------------------- 3. BUILD GRAPH (PGGB) ---------------------------
echo "[STAGE 2] Building pangenome graph with pggb"

# -n = expected number of haplotypes/genomes (approx)
# -p = minimum sequence identity for mapping (%), tune per divergence
# -s = segment length for mapping
pggb \
  -i "$INPUT_FA" \
  -o "$OUTDIR/01_graph" \
  -n "$N_SEQS" \
  -p 90 \
  -s 5000 \
  -t "$THREADS" \
  --temp-dir "$OUTDIR/logs" \
  2>&1 | tee "$OUTDIR/logs/01_pggb.log"

GRAPH_GFA=$(find "$OUTDIR/01_graph" -maxdepth 1 -name "*.smooth.final.gfa" | head -n1)
[[ -z "$GRAPH_GFA" ]] && { echo "[ERROR] pggb output GFA not found"; exit 1; }
echo "[INFO] Graph GFA: $GRAPH_GFA"

# ---------------------------- 4. GRAPH -> ODGI (normalize + stats) ------------
echo "[STAGE 3] Converting to odgi format & computing graph stats"

ODGI_GRAPH="$OUTDIR/02_odgi/pangenome.og"
odgi build -g "$GRAPH_GFA" -o "$ODGI_GRAPH" -t "$THREADS"

odgi stats -i "$ODGI_GRAPH" -S > "$OUTDIR/05_stats/graph_stats.tsv"
odgi viz   -i "$ODGI_GRAPH" -o "$OUTDIR/05_stats/graph_viz.png" -x 4000 -y 800
odgi paths -i "$ODGI_GRAPH" -L > "$OUTDIR/05_stats/path_list.txt"

# Optional: core/accessory/pangenome openness plot (needs panacus)
if command -v panacus >/dev/null; then
  panacus hist "$GRAPH_GFA" > "$OUTDIR/05_stats/pangenome_growth.tsv"
fi

# ---------------------------- 5. VARIANT CALLING FROM GRAPH -------------------
echo "[STAGE 4] Converting to vg format and calling variants against reference path"

VG_GRAPH="$OUTDIR/03_variants/pangenome.vg"
vg convert -g "$GRAPH_GFA" -p > "$VG_GRAPH"

vg index -x "$OUTDIR/03_variants/pangenome.xg" "$VG_GRAPH" -t "$THREADS"

REF_PATH_NAME=$(basename "$REF" | sed 's/\.[^.]*$//')

# vg deconstruct: extract variants of all paths relative to reference path(s)
vg deconstruct \
  -p "$REF_PATH_NAME" \
  -t "$THREADS" \
  "$OUTDIR/03_variants/pangenome.xg" \
  > "$OUTDIR/03_variants/pangenome.raw.vcf"

# Normalize/clean multiallelic sites and nested bubbles
bcftools sort "$OUTDIR/03_variants/pangenome.raw.vcf" \
  -Oz -o "$OUTDIR/03_variants/pangenome.sorted.vcf.gz"
tabix -p vcf "$OUTDIR/03_variants/pangenome.sorted.vcf.gz"

# Optional: collapse nested/overlapping bubbles with vcfbub
if command -v vcfbub >/dev/null; then
  vcfbub -l 0 -a 100000 \
    -i "$OUTDIR/03_variants/pangenome.sorted.vcf.gz" \
    > "$OUTDIR/03_variants/pangenome.clean.vcf"
  bgzip -f "$OUTDIR/03_variants/pangenome.clean.vcf"
  tabix -p vcf "$OUTDIR/03_variants/pangenome.clean.vcf.gz"
  FINAL_VCF="$OUTDIR/03_variants/pangenome.clean.vcf.gz"
else
  FINAL_VCF="$OUTDIR/03_variants/pangenome.sorted.vcf.gz"
fi

bcftools stats "$FINAL_VCF" > "$OUTDIR/05_stats/vcf_stats.txt"

# ---------------------------- 6. ANNOTATION ------------------------------------
echo "[STAGE 5] Annotating variants"

# 6a. Build/use a snpEff database for the reference genome
SNPEFF_DB_NAME="pangenome_ref"
SNPEFF_DATA_DIR="$OUTDIR/04_annotation/snpeff_data"
mkdir -p "$SNPEFF_DATA_DIR/$SNPEFF_DB_NAME"

cp "$REF" "$SNPEFF_DATA_DIR/$SNPEFF_DB_NAME/sequences.fa"
cp "$REF_GFF" "$SNPEFF_DATA_DIR/$SNPEFF_DB_NAME/genes.gff"

cat > "$OUTDIR/04_annotation/snpeff.config" <<EOF
data.dir = $SNPEFF_DATA_DIR
$SNPEFF_DB_NAME.genome : $SNPEFF_DB_NAME
EOF

snpEff build -gff3 -c "$OUTDIR/04_annotation/snpeff.config" \
  -v "$SNPEFF_DB_NAME" -noCheckCds -noCheckProtein \
  2>&1 | tee "$OUTDIR/logs/05_snpeff_build.log"

# 6b. Annotate the final variant callset (functional consequence: missense,
#     synonymous, UTR, intergenic, frameshift, etc.)
snpEff -c "$OUTDIR/04_annotation/snpeff.config" \
  "$SNPEFF_DB_NAME" "$FINAL_VCF" \
  > "$OUTDIR/04_annotation/pangenome.annotated.vcf" \
  2> "$OUTDIR/logs/06_snpeff_annotate.log"

bgzip -f "$OUTDIR/04_annotation/pangenome.annotated.vcf"
tabix -p vcf "$OUTDIR/04_annotation/pangenome.annotated.vcf.gz"

# 6c. Project reference gene annotation (GFF3) onto the pangenome graph itself,
#     so gene coordinates are queryable directly on graph nodes/paths.
vg annotate \
  -x "$OUTDIR/03_variants/pangenome.xg" \
  -f "$REF_GFF" \
  -p > "$OUTDIR/04_annotation/pangenome.graph_annotated.gam" \
  2>> "$OUTDIR/logs/06_snpeff_annotate.log" || \
  echo "[WARN] vg annotate step skipped/failed — check vg version supports -f gff mode"

# ---------------------------- 7. SUMMARY ---------------------------------------
echo "[DONE] Pipeline complete."
echo "  Graph (GFA)         : $GRAPH_GFA"
echo "  Graph (odgi)         : $ODGI_GRAPH"
echo "  Graph stats           : $OUTDIR/05_stats/graph_stats.tsv"
echo "  Graph visualization    : $OUTDIR/05_stats/graph_viz.png"
echo "  Raw variants (VCF)     : $OUTDIR/03_variants/pangenome.sorted.vcf.gz"
echo "  Annotated variants     : $OUTDIR/04_annotation/pangenome.annotated.vcf.gz"
echo "  snpEff HTML summary    : snpEff_summary.html (current dir)"
echo "  Graph-level annotation : $OUTDIR/04_annotation/pangenome.graph_annotated.gam"
