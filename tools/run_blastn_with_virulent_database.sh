#!/bin/bash
# ===================================================
# Usage:
#   bash run_blastn_with_virulent_database.sh \
#     -i <input_fasta> \
#     --prefix <output_prefix> \
#     [--outdir <output_dir>] \
#     [-db <database_path>]
#
# Example:
#   bash run_blastn_with_virulent_database.sh \
#     -i barcode07_trimmed.fasta \
#     --prefix barcode07_vfdb \
#     --outdir 06_blast_vfdb
#
# Output files (in outdir):
#   <prefix>.blast.tsv
#   <prefix>.best.tsv
#   <prefix>.sseqid_counts.tsv   (columns: sseqid, n_qseqid_detected, RPM)
#   <prefix>.species_counts.tsv  (columns: species, n_qseqid_detected, RPM)
#   RPM = n_qseqid_detected / total_reads_in_input_fasta * 1,000,000
# ===================================================

set -euo pipefail

# ==== 1️⃣ Default parameters ====
DB_DEFAULT="/mnt/d/YiLun/ncbi_db/virulent_database_from_paper_vfdb/virulent_database_from_paper_vfdb.fasta"
THREADS=10
EVALUE=1e-50
IDENTITY=95
MAX_TARGET=5

# ==== 2️⃣ Parse arguments ====
while [[ $# -gt 0 ]]; do
  case "$1" in
    -i|--input)
      QUERY="$2"; shift 2 ;;
    --prefix)
      PREFIX="$2"; shift 2 ;;
    --outdir)
      OUTDIR="$2"; shift 2 ;;
    -db|--database)
      DB="$2"; shift 2 ;;
    -h|--help)
      echo "Usage: $0 -i <input_fasta> --prefix <output_prefix> [--outdir <dir>] [-db <database_path>]"
      exit 0 ;;
    *)
      echo "[ERROR] Unknown parameter: $1"
      echo "Usage: $0 -i <input_fasta> --prefix <output_prefix> [--outdir <dir>] [-db <database_path>]"
      exit 1 ;;
  esac
done

# ==== 3️⃣ Check required arguments ====
if [[ -z "${QUERY:-}" || -z "${PREFIX:-}" ]]; then
  echo "[ERROR] -i and --prefix are required"
  echo "Usage: $0 -i <input_fasta> --prefix <output_prefix> [--outdir <dir>] [-db <database_path>]"
  exit 1
fi

DB="${DB:-$DB_DEFAULT}"
DB_EXPANDED=$(eval echo "$DB")

OUTDIR="${OUTDIR:-.}"
mkdir -p "$OUTDIR"

# ==== 4️⃣ Count total reads in input FASTA (for RPM) ====
echo "[INFO] Counting total reads in input FASTA..."
TOTAL_READS=$(grep -c '^>' "$QUERY")
echo "[INFO] Total reads: ${TOTAL_READS}"

# ==== 5️⃣ Output filenames from prefix ====
OUT="${OUTDIR}/${PREFIX}.blast.tsv"
OUT_BEST="${OUTDIR}/${PREFIX}.best.tsv"
OUT_SSEQID_SUM="${OUTDIR}/${PREFIX}.sseqid_counts.tsv"
OUT_SPECIES_SUM="${OUTDIR}/${PREFIX}.species_counts.tsv"

# ==== 6️⃣ Temp files (avoid collision in parallel runs) ====
TMP_BLAST="$(mktemp "${OUTDIR}/tmp_blast.XXXXXX.tsv")"
TMP_WITH_HEADER="$(mktemp "${OUTDIR}/tmp_with_header.XXXXXX.tsv")"

cleanup() {
  rm -f "$TMP_BLAST" "$TMP_WITH_HEADER"
}
trap cleanup EXIT

# ==== 6️⃣ Run blastn ====

echo "[INFO] Running blastn..."
blastn -task megablast \
  -query "$QUERY" \
  -db "$DB_EXPANDED" \
  -evalue $EVALUE \
  -outfmt "6 qseqid sseqid pident length qlen slen qcovs evalue bitscore" \
  -perc_identity $IDENTITY \
  -max_target_seqs $MAX_TARGET \
  -num_threads $THREADS \
  > "$TMP_BLAST"

# ==== 7️⃣ Calculate scovs + filter + sort ====
SCOVS_MIN=0  # ← 這裡改你要的 scovs 門檻 (%)

echo "[INFO] Calculating scovs (>= ${SCOVS_MIN}%) and sorting..."

awk -v min_scovs="$SCOVS_MIN" 'BEGIN{
  OFS="\t";
  print "qseqid","sseqid","pident","length","qlen","slen","qcovs","evalue","bitscore","scovs"
}
{
  scovs=($4/$6)*100;
  if (scovs >= min_scovs) {
    printf "%s\t%s\t%.2f\t%d\t%d\t%d\t%.2f\t%s\t%s\t%.2f\n",
           $1,$2,$3,$4,$5,$6,$7,$8,$9,scovs
  }
}' "$TMP_BLAST" > "$TMP_WITH_HEADER"

{
  head -n 1 "$TMP_WITH_HEADER"
  tail -n +2 "$TMP_WITH_HEADER" \
    | LC_ALL=C sort -t $'\t' -k10,10nr -k3,3nr
} > "$OUT"


# ==== 8️⃣ Best hit per qseqid (cap scovs > 100 => 100, then pident) ====
echo "[INFO] Selecting best hit per qseqid..."
awk -F'\t' 'BEGIN{OFS="\t"}
NR==1{print $0,"\tscovs_capped"; next}
{
  sc=$10+0; if(sc>100) sc=100;
  print $0, sprintf("%.2f", sc)
}' "$OUT" \
| {
    read -r header
    echo -e "$header"
    LC_ALL=C sort -t $'\t' -k1,1 -k11,11nr -k3,3nr \
    | awk -F'\t' 'BEGIN{OFS="\t"} !seen[$1]++'
  } > "$OUT_BEST"

# ==== 9️⃣ sseqid summary (sort first, then add header) ====
echo "[INFO] Summarizing by sseqid..."
awk -F'\t' -v total="$TOTAL_READS" 'NR>1{cnt[$2]++}
END{
  for(k in cnt) printf "%s\t%d\t%.4f\n", k, cnt[k], (cnt[k]/total)*1e6
}' "$OUT_BEST" \
| LC_ALL=C sort -t $'\t' -k2,2nr -k1,1 \
| {
    echo -e "sseqid\tn_qseqid_detected\tRPM"
    cat
  } > "$OUT_SSEQID_SUM"

# ==== 🔟 Species (Genus_species) summary (sort first, then add header) ====
echo "[INFO] Summarizing by species..."
awk -F'\t' -v total="$TOTAL_READS" 'NR>1{
  split($2,a,"_");
  sp=a[1]"_"a[2];
  cnt[sp]++
}
END{
  for(k in cnt) printf "%s\t%d\t%.4f\n", k, cnt[k], (cnt[k]/total)*1e6
}' "$OUT_BEST" \
| LC_ALL=C sort -t $'\t' -k2,2nr -k1,1 \
| {
    echo -e "species\tn_qseqid_detected\tRPM"
    cat
  } > "$OUT_SPECIES_SUM"

echo "[DONE]"
echo "  - $OUT"
echo "  - $OUT_BEST"
echo "  - $OUT_SSEQID_SUM"
echo "  - $OUT_SPECIES_SUM"