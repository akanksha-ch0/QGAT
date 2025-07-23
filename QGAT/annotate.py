import pandas as pd
import argparse
import re
from pathlib import Path
import os  #  Added for summary file paths

def sort_and_load_gtf(gtf_path, is_ncbi=False):
    try:
        gtf_df = pd.read_csv(gtf_path, sep="\t", comment="#", header=None, names=[
            "chromosome", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"
        ], dtype=str)
    except Exception as e:
        raise FileNotFoundError(f" GTF file not found or unreadable: {e}")

    if gtf_df.empty:
        raise ValueError(" GTF file is empty after parsing.")

    if is_ncbi:
        gtf_df["gene_id"] = gtf_df["attribute"].str.extract(r'db_xref "GeneID:([^"]+)"')[0]
        gtf_df["gene_name"] = gtf_df["attribute"].str.extract(r'gene "([^"]+)"')[0]
    else:
        gtf_df["gene_id"] = gtf_df["attribute"].str.extract(r'gene_id "([^"]+)"')[0]
        gtf_df["gene_name"] = gtf_df["attribute"].str.extract(r'gene_name "([^"]+)"')[0]

    gtf_df["gene_biotype"] = gtf_df["attribute"].str.extract(r'gene_biotype "([^"]+)"')[0]

    # Keep only gene-like rows
    gtf_df = gtf_df[gtf_df["feature"].isin(["gene", "transcript", "exon"])]

    gtf_df = gtf_df.dropna(subset=["chromosome", "start", "end", "gene_id"])
    gtf_df["start"] = pd.to_numeric(gtf_df["start"], errors="coerce")
    gtf_df["end"] = pd.to_numeric(gtf_df["end"], errors="coerce")
    gtf_df = gtf_df.dropna(subset=["start", "end"])

    gtf_df["width"] = gtf_df["end"] - gtf_df["start"] + 1

    return gtf_df[[
        "chromosome", "start", "end", "width", "strand",
        "gene_id", "gene_name", "gene_biotype"
    ]]


def annotate_regions(input_file, gtf_df, is_ncbi=False):
    try:
        input_df = pd.read_csv(input_file, sep="\t", dtype=str)
    except Exception as e:
        raise FileNotFoundError(f" Input file could not be read: {e}")

    if input_df.empty:
        raise ValueError(" Input file is empty.")

    required_cols = {"chromosome", "start", "end"}
    if not required_cols.issubset(set(input_df.columns)):
        raise ValueError(f" Input file must contain columns: {required_cols}")

    input_df["start"] = pd.to_numeric(input_df["start"], errors="coerce")
    input_df["end"] = pd.to_numeric(input_df["end"], errors="coerce")
    input_df = input_df.dropna(subset=["start", "end"])

    overlaps = []
    for _, region in input_df.iterrows():
        chrom, start, end = region["chromosome"], region["start"], region["end"]
        match = gtf_df[
            (gtf_df["chromosome"] == chrom) &
            (gtf_df["end"] >= start) &
            (gtf_df["start"] <= end)
        ]
        for _, gene in match.iterrows():
            overlaps.append({
                "chromosome": chrom,
                "start": start,
                "end": end,
                "gene_chr": gene["chromosome"],
                "gene_start_pos": gene["start"],
                "gene_end_pos": gene["end"],
                "width": gene["width"],
                "strand": gene["strand"],
                "gene_id": gene["gene_id"],
                "gene_name": gene["gene_name"],
                "gene_biotype": gene["gene_biotype"]
            })

    # NEW: robust mismatch check
    if not overlaps:
        input_chrs = set(input_df["chromosome"].unique())
        gtf_chrs = set(gtf_df["chromosome"].unique())
        common_chrs = input_chrs.intersection(gtf_chrs)

        if not common_chrs:
            raise ValueError(
                f"\n No overlapping genes found AND no matching chromosome names.\n"
                f"  Input file chromosomes: {sorted(input_chrs)}\n"
                f"  GTF file chromosomes: {sorted(gtf_chrs)}\n"
                f"  → Possible cause: input uses plain numbers (1,2,3) but GTF uses NCBI IDs (e.g., NC_037328.1).\n"
                f"  → Fix: Adjust input chromosome names to match GTF, or use --ncbi if not already.\n"
            )
        else:
            print(" No overlapping genes found but chromosome names DO match.")
        return pd.DataFrame()

    result_df = pd.DataFrame(overlaps)
    result_df["gene_name"] = result_df["gene_name"].fillna("NA").replace("", "NA")
    return result_df

# Add this run_annotate implementation at the end:
def run_annotate(args):
    print("\n🔄 Sorting and parsing GTF file...")
    gtf_df = sort_and_load_gtf(args.gtf, is_ncbi=args.ncbi)

    print(" Annotating regions...")
    annotated_df = annotate_regions(args.input, gtf_df, is_ncbi=args.ncbi)

    if annotated_df.empty:
        print("  No overlapping genes found.")
    else:
        annotated_df.to_csv(args.output, sep='\t', index=False)

    total_regions = pd.read_csv(args.input, sep='\t').shape[0]
    total_overlaps = annotated_df.shape[0]

    print("\n Annotated output saved to {}".format(args.output))
    print("\n Total input regions: {}".format(total_regions))
    print(" Total overlapping genes: {}".format(total_overlaps))

    if not annotated_df.empty and 'gene_name' in annotated_df.columns:
        unique_genes_df = annotated_df.drop_duplicates(subset=['gene_name'])
        total_unique_genes = unique_genes_df.shape[0]

        print(" Total unique genes: {}".format(total_unique_genes))

        # Write summary text file
        summary_file = os.path.splitext(args.output)[0] + "_summary.txt"
        with open(summary_file, 'w') as f:
            f.write("QGAT Annotation Summary\n")
            f.write("=======================\n")
            f.write(f"Input file: {args.input}\n")
            f.write(f"GTF file: {args.gtf}\n")
            f.write(f"Output file: {args.output}\n")
            f.write(f"Total input regions: {total_regions}\n")
            f.write(f"Total overlapping genes: {total_overlaps}\n")
            f.write(f"Total unique genes: {total_unique_genes}\n")

        print(f"\n Summary saved to {summary_file}")

        # Save unique genes table (all columns, no duplicates)
        unique_file = os.path.splitext(args.output)[0] + "_unique.tsv"
        unique_genes_df.to_csv(unique_file, sep='\t', index=False)
        print(f" Unique genes table saved to {unique_file}")

    elif not annotated_df.empty:
        print(" Warning: 'gene_name' column not found. Cannot compute unique genes.")
