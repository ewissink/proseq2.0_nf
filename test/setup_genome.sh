#!/usr/bin/env bash
#
# setup_genome.sh — download a genome + build the BWA index and chromInfo for testing.
#
# Usage: test/setup_genome.sh {mm10|hg38} [OUTDIR]
#
# NOTE: this is the heavy, one-time step (mammalian `bwa index` takes ~1 h and
# several GB). If you already have an index, skip this and point run_concordance.sh
# at it instead.
#
# rDNA caveat: proseq2.0 recommends an index that INCLUDES a copy of the rDNA
# (GenBank U13369.1) so rRNA reads don't misalign; the plain UCSC genome below does
# not. That only affects rRNA *filtering* / biological interpretation — it does NOT
# affect the old-vs-new concordance test, since both pipelines use the same index.
# For real analysis, append the rDNA sequence to the FASTA before indexing.
#
set -euo pipefail

GENOME="${1:-}"
OUTDIR="${2:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/ref}"
[ -n "$GENOME" ] || { echo "Usage: $0 {mm10|hg38} [OUTDIR]" >&2; exit 1; }

case "$GENOME" in
  mm10) BASE="https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips" ;;
  hg38) BASE="https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips" ;;
  *) echo "ERROR: genome must be mm10 or hg38" >&2; exit 1 ;;
esac

for tool in wget bwa; do
  command -v "$tool" >/dev/null 2>&1 || { echo "ERROR: '$tool' not on PATH (activate test env)" >&2; exit 1; }
done

mkdir -p "$OUTDIR"; cd "$OUTDIR"

echo ">> Downloading $GENOME FASTA + chrom.sizes ..."
[ -s "${GENOME}.fa.gz" ]      || wget -q --show-progress "$BASE/${GENOME}.fa.gz"
[ -s "${GENOME}.chrom.sizes" ] || wget -q --show-progress "$BASE/${GENOME}.chrom.sizes"
cp -f "${GENOME}.chrom.sizes" "${GENOME}.chromInfo"

echo ">> Building BWA index (this is the slow part) ..."
if [ ! -s "${GENOME}.fa.gz.bwt" ]; then
  bwa index "${GENOME}.fa.gz"
fi

echo
echo "Done. Pass these to run_concordance.sh:"
echo "  --bwa-index  $OUTDIR/${GENOME}.fa.gz"
echo "  --chrom-info $OUTDIR/${GENOME}.chromInfo"
