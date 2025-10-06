Scripts utilized to identify, map, and quantify deletion boundaries for the publication "Engineered Type IV-A and Type I-Fv CRISPR-Cas effectors function in reduction of the _Escherichia coli_ genome (Klein et al, 2025)."

### SPLICED ALIGNMENTS FILTERING
### splice_extractor.sh
Given one or more BAM files in the current directory, the script:
- Keeps only spliced alignments (CIGAR = N)
- Applies mapping-quality filtering (MAPQ >= 60)
- Writes a BAM file with filtered reads and index it
- These files are  the inputs for the detect_junctions.py script

### Requirements
- bash
- samtools >= 1.1
- awk

### JUNCTION DETECTION
### detect_junctions.py
Given one or more coordinate-sorted, indexed BAMs, a reference FASTA, and a GFF3 annotation, the script:
- Detects splice/skipped-region junctions from each BAM (CIGAR N operations).
- Clusters nearly identical junctions using a user-defined tolerance (in bp).
- Filters to a user-specified sequence region.
- Annotates each junction with overlapping Left Feature (at the start site) and Right Feature (at the end site) based on GFF3 attributes (locus_tag or gene).

Outputs per-BAM:
- an interactive HTML report (Plotly + DataTables),
- a TSV table of junctions,
- a BED track (for IGV),
- four FASTA files with sequence context around start/end sites.

### Requirements
- Python 3.8+
- Python packages: pysam, pandas, numpy, plotly, biopython

### Input files
- BAM files: coordinate-sorted and indexed (output from splice_extractor.sh).
- Reference FASTA (indexed)
- GFF3 annotation file

### USAGE
python detect_junctions.py \
  --bamdir path/to/bams \
  --outdir path/to/output \
  --gff annotations.gff3 \
  --fasta genome.fa \
  --region chr:start-end \
  [--tol 3] \
  [--flank 50]

## Arguments
--bamdir

Directory containing one or more *.bam files (each must have a .bai index in the same directory).

--outdir

Output directory (created if it does not exist).

--gff

GFF3 annotation. The script reads seqid, start, end, attributes and looks inside attributes for locus_tag= or gene= to annotate Left Feature/Right Feature.

--fasta

Reference FASTA (must have .fai index).

--region

Genomic window of interest in the form chrom:start-end.
Junctions are included if either the start or the end falls inside this window.

--tol (default: 3)

Clustering tolerance in nucleotides for both junction start and end. Junctions within ±tol bp of a cluster’s mean start and end are merged.

--flank (default: 50)

Number of nucleotides to extract upstream (left) and downstream (right) around the start and end sites for FASTA outputs.

### Generated outputs per BAM files
- sample.tsv - Tab-separated table with columns:
junction_id (chrom:start-end), chrom, start, end, length, count, Left Feature, Right Feature, left_start, right_start, left_end, right_end (sequence strings; --flank nt each)

- sample.bed - BED track with columns that can be loaded in IGV to visualize junction spans.
chrom start end junction_id count strand(.)

Fasta files for each flanking region:
- sample.left_start.fa
- sample.right_start.fa
- sample.left_end.fa
- sample.right_end.fa

- sample.html - Interactive HTML report with:
    - A Plotly horizontal bar visualization (each bar = one junction; start, end, length, count, and features in hover),
    - A DataTables table of all columns,
    - Paths to the generated TSV/BED/FASTA files.


### Column information
- count: number of reads (in that BAM) supporting the clustered junction.

- length: junction span = end - start.

- Left/Right Feature: locus_tag or gene covering the start or end site (inclusive overlap). Useful to see whether a junction connects two exons/genes or occurs within a feature.

- left_/right_ (start/end): flanking sequence context.

- left_*: upstream of the site

- right_*: downstream of the site
