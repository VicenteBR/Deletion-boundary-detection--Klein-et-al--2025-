#!/usr/bin/env bash
set -euo pipefail

# Make sure the output directory exists
OUTDIR="filtered_bams"
mkdir -p "$OUTDIR"

# Loop over all .bam files in this folder
for f in *.bam; do
  # Derive a base name (e.g. sample.bam -> sample.spliced.bam)
  base=$(basename "$f" .bam)
  out_bam="$OUTDIR/${base}.spliced.bam"
  sorted_bam="$OUTDIR/${base}.spliced.sorted.bam"

  echo "Filtering spliced reads with MAPQ >= 60 from $f → $out_bam"

  # Extract header + only CIGARs with 'N' and MAPQ >= 60, write to unsorted BAM
  {
    samtools view -H "$f"
    samtools view "$f" \
      | awk '$6 ~ /N/ && $5 >= 60'
  } | samtools view -b -o "$out_bam" -

  # Sort & index
  samtools sort -@12 -o "$sorted_bam" "$out_bam"
  samtools index "$sorted_bam"

  # Optionally remove the unsorted BAM to save space:
  rm "$out_bam"

  echo "  → sorted and indexed as $sorted_bam"
done

echo "All done. Filtered BAMs are in $OUTDIR/" 
