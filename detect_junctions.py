#!/usr/bin/env python3
"""
Script to detect and summarize splice junctions from BAM files,
cluster them, and generate:
  1) An interactive HTML report per region
  2) A TSV and a BED track for IGV visualization
  3) Four FASTA files with sequence context:
     - left_start.fa, right_start.fa, left_end.fa, right_end.fa

Features now include Left Feature and Right Feature (locus_tag or gene name).

Requires: pysam, pandas, numpy, plotly, Biopython
Usage:
    python detect_junctions.py \
      --bamdir path/to/bams --outdir path/to/output \
      --gff annotations.gff3 --fasta genome.fa \
      --region chr:start-end [--tol 3] [--flank 50]
"""
import os
import argparse
import pysam
import numpy as np
import pandas as pd
import plotly.graph_objects as go
import plotly.io as pio


def extract_junctions(bam_file):
    juncs = []
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        for read in bam.fetch():
            if read.cigartuples:
                pos = read.reference_start
                for op, length in read.cigartuples:
                    if op == 3:  # N (skipped region) = junction
                        juncs.append((read.reference_name, pos, pos + length))
                        pos += length
                    elif op in (0, 2, 7, 8):  # M, D, =, X move along reference
                        pos += length
    return juncs


def cluster_junctions(juncs, tol):
    clusters = []
    for chrom, s, e in juncs:
        for cl in clusters:
            if cl['chrom'] == chrom and abs(s - cl['start_mean']) <= tol and abs(e - cl['end_mean']) <= tol:
                cl['count'] += 1
                cl['start_mean'] = (cl['start_mean'] * (cl['count'] - 1) + s) / cl['count']
                cl['end_mean']   = (cl['end_mean']   * (cl['count'] - 1) + e) / cl['count']
                break
        else:
            clusters.append({'chrom': chrom, 'start_mean': s, 'end_mean': e, 'count': 1})
    out = []
    for cl in clusters:
        s = int(round(cl['start_mean']))
        e = int(round(cl['end_mean']))
        out.append({
            'junction_id': f"{cl['chrom']}:{s}-{e}",
            'chrom': cl['chrom'],
            'start': s,
            'end': e,
            'count': cl['count'],
            'length': e - s
        })
    return out


def read_gff(path):
    cols = ['seqid','source','type','start','end','score','strand','phase','attributes']
    df = pd.read_csv(path, sep='\t', comment='#', names=cols)
    df['start'] = pd.to_numeric(df['start'], errors='coerce')
    df['end']   = pd.to_numeric(df['end'], errors='coerce')
    df = df.dropna(subset=['start','end']).copy()
    df['start'] = df['start'].astype(int)
    df['end']   = df['end'].astype(int)
    return df


def extract_flanks_handle(fasta_handle, chrom, pos, flank):
    """Extract left (upstream) and right (downstream) flanks around a genomic pos."""
    left  = fasta_handle.fetch(chrom, max(0, pos - flank), pos)
    right = fasta_handle.fetch(chrom, pos, pos + flank)
    return left, right


def parse_attribute(attr_str, keys):
    for part in attr_str.split(';'):
        kv = part.strip().split('=')
        if len(kv) == 2 and kv[0] in keys:
            return kv[1]
    return None


def assign_features(row, gff_df):
    chrom, row_s, row_e = row['chrom'], row['start'], row['end']
    sub = gff_df[gff_df['seqid'] == chrom]
    # left feature at start, right feature at end
    left  = sub[(sub.start <= row_s) & (sub.end >= row_s)]
    right = sub[(sub.start <= row_e) & (sub.end >= row_e)]
    lf = left.attributes.map(lambda x: parse_attribute(x, ['locus_tag', 'gene'])).dropna().unique()
    rf = right.attributes.map(lambda x: parse_attribute(x, ['locus_tag', 'gene'])).dropna().unique()
    return (','.join(lf) if len(lf) > 0 else 'None',
            ','.join(rf) if len(rf) > 0 else 'None')


def write_fasta(df, seq_col, out_path, site_label):
    """
    Write sequences from df[seq_col] to a FASTA file.
    Header format:
      >{junction_id}|{site_label}|count={count}|chrom={chrom}|pos={pos}
    where pos is start or end depending on site_label.
    """
    with open(out_path, 'w') as fh:
        for _, r in df.iterrows():
            pos = r['start'] if 'start' in site_label else r['end']
            header = f">{r['junction_id']}|{site_label}|count={r['count']}|chrom={r['chrom']}|pos={pos}"
            seq = str(r[seq_col]).strip()
            fh.write(header + "\n")
            fh.write(seq + "\n")


def generate_html_and_tracks(clusters, region, gff_df, fasta_path, flank, outdir, base):
    chrom, rs, re = region
    df = pd.DataFrame(clusters)

    # include junctions that start inside or end inside the region
    df = df[(df['chrom'] == chrom) & ((df['start'].between(rs, re)) | (df['end'].between(rs, re)))].reset_index(drop=True)

    # open FASTA once
    fa = pysam.FastaFile(fasta_path)

    # sequence context for start site (donor)
    df[['left_start','right_start']] = df.apply(
        lambda r: pd.Series(extract_flanks_handle(fa, chrom, int(r.start), flank)),
        axis=1
    )
    # sequence context for end site (acceptor)
    df[['left_end','right_end']] = df.apply(
        lambda r: pd.Series(extract_flanks_handle(fa, chrom, int(r.end), flank)),
        axis=1
    )

    # assign left/right features (based on overlap at start and end)
    feats = df.apply(lambda r: assign_features(r, gff_df), axis=1)
    df['Left Feature']  = feats.map(lambda x: x[0])
    df['Right Feature'] = feats.map(lambda x: x[1])

    # save TSV
    tsv = os.path.join(outdir, f"{base}.tsv")
    df.to_csv(tsv, sep='\t', index=False)

    # save BED
    bed = os.path.join(outdir, f"{base}.bed")
    bed_df = df[['chrom','start','end','junction_id','count']].copy()
    bed_df['strand'] = '.'
    bed_df.to_csv(bed, sep='\t', header=False, index=False)

    # write four FASTAs
    left_start_fa  = os.path.join(outdir, f"{base}.left_start.fa")
    right_start_fa = os.path.join(outdir, f"{base}.right_start.fa")
    left_end_fa    = os.path.join(outdir, f"{base}.left_end.fa")
    right_end_fa   = os.path.join(outdir, f"{base}.right_end.fa")

    write_fasta(df, 'left_start',  left_start_fa,  'start_left')
    write_fasta(df, 'right_start', right_start_fa, 'start_right')
    write_fasta(df, 'left_end',    left_end_fa,    'end_left')
    write_fasta(df, 'right_end',   right_end_fa,   'end_right')

    fa.close()

    # Plotly visualization: bars per junction
    fig = go.Figure()
    y_positions = [i * 0.2 for i in range(len(df))]
    for idx, row in df.iterrows():
        y = y_positions[idx]
        length = row['length']
        start = row['start']
        gid = row['junction_id']
        cnt = row['count']
        lf = row['Left Feature']
        rf = row['Right Feature']
        fig.add_trace(go.Bar(
            x=[length],
            y=[y],
            base=[start],
            orientation='h',
            marker=dict(color='steelblue'),
            showlegend=False,
            hovertemplate=(
                f"ID: {gid}<br>"
                f"Start: {start} | End: {row['end']}<br>"
                f"Length: {length}<br>"
                f"Count: {cnt}<br>"
                f"Left Feature: {lf}<br>"
                f"Right Feature: {rf}"
            )
        ))
    maxy = y_positions[-1] if y_positions else 0
    fig.update_layout(
        title=f"Junctions {chrom}:{rs}-{re}",
        xaxis_title='Genomic pos',
        yaxis=dict(showticklabels=False, range=[-0.5, maxy + 0.2]),
        bargap=0.2,
        template='plotly_white'
    )

    # HTML
    html = os.path.join(outdir, f"{base}.html")
    cols = [
        'junction_id','start','end','length','count',
        'Left Feature','Right Feature',
        'left_start','right_start','left_end','right_end'
    ]
    header = ''.join(f"<th>{c}</th>" for c in cols)
    rows = ''.join(
        '<tr>' + ''.join(f'<td>{r[c]}</td>' for c in cols) + '</tr>'
        for _, r in df.iterrows()
    )
    fasta_links = (
        f"<ul>"
        f"<li>Left (start): {left_start_fa}</li>"
        f"<li>Right (start): {right_start_fa}</li>"
        f"<li>Left (end): {left_end_fa}</li>"
        f"<li>Right (end): {right_end_fa}</li>"
        f"</ul>"
    )
    content = f"""
<html><head><meta charset='utf-8'>
<script src='https://cdn.plot.ly/plotly-latest.min.js'></script>
<link rel='stylesheet' href='https://cdn.datatables.net/1.10.24/css/jquery.dataTables.min.css'>
<script src='https://code.jquery.com/jquery-3.5.1.js'></script>
<script src='https://cdn.datatables.net/1.10.24/js/jquery.dataTables.min.js'></script>
</head><body>
<h1>{base} {chrom}:{rs}-{re}</h1>
<div>{pio.to_html(fig, include_plotlyjs=False, full_html=False)}</div>
<h2>Outputs</h2>
<p>TSV: {tsv}<br>BED: {bed}</p>
<h3>FASTA files</h3>
{fasta_links}
<h2>Details</h2>
<table id='tbl'><thead><tr>{header}</tr></thead><tbody>
{rows}
</tbody></table>
<script>$(function(){{$('#tbl').DataTable();}});</script>
</body></html>
"""
    with open(html, 'w') as f:
        f.write(content)
    print(f"Saved TSV: {tsv}, BED: {bed}, HTML: {html}")
    print(f"Saved FASTAs: {left_start_fa}, {right_start_fa}, {left_end_fa}, {right_end_fa}")


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--bamdir', required=True)
    p.add_argument('--outdir', required=True)
    p.add_argument('--gff', required=True)
    p.add_argument('--fasta', required=True)
    p.add_argument('--region', required=True)
    p.add_argument('--tol', type=int, default=3)
    p.add_argument('--flank', type=int, default=50)
    args = p.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    gff_df = read_gff(args.gff)

    chrom, coords = args.region.split(':')
    rs, re = map(int, coords.split('-'))

    for bam in sorted(os.listdir(args.bamdir)):
        if not bam.endswith('.bam'):
            continue
        base = os.path.splitext(bam)[0]
        j = extract_junctions(os.path.join(args.bamdir, bam))
        cl = cluster_junctions(j, args.tol)
        generate_html_and_tracks(cl, (chrom, rs, re), gff_df, args.fasta, args.flank, args.outdir, base)


if __name__ == '__main__':
    main()
