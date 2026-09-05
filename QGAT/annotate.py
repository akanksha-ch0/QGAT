import pandas as pd
import argparse
import re
from pathlib import Path
import os
import sys

from QGAT.utils import parse_input_regions


def sort_and_load_gtf(gtf_path, is_ncbi=False, features=None):
    """
    Parse a GTF/GFF file and return a DataFrame of gene-level features.

    Parameters
    ----------
    gtf_path : str
        Path to the GTF file.
    is_ncbi : bool
        If True, parse NCBI-style attributes (db_xref, gene).
    features : list of str or None
        GTF feature types to keep (e.g., ['gene', 'transcript']).
        Defaults to ['gene', 'transcript', 'exon'].

    Returns
    -------
    pd.DataFrame with columns:
        chromosome, start, end, width, strand, gene_id, gene_name, gene_biotype
    """
    if features is None:
        features = ["gene", "transcript", "exon"]

    try:
        gtf_df = pd.read_csv(gtf_path, sep="\t", comment="#", header=None, names=[
            "chromosome", "source", "feature", "start", "end",
            "score", "strand", "frame", "attribute"
        ], dtype=str)
    except Exception as e:
        raise FileNotFoundError(f"GTF file not found or unreadable: {e}")

    if gtf_df.empty:
        raise ValueError("GTF file is empty after parsing.")

    # --- Extract gene attributes -----------------------------------------
    if is_ncbi:
        gtf_df["gene_id"] = gtf_df["attribute"].str.extract(r'db_xref "GeneID:([^"]+)"')[0]
        # Try 'gene' attribute first, then fall back to 'gene_id' for NCBI
        gtf_df["gene_name"] = gtf_df["attribute"].str.extract(r'gene "([^"]+)"')[0]
        # If gene_name is still NA, try product
        mask = gtf_df["gene_name"].isna()
        if mask.any():
            gtf_df.loc[mask, "gene_name"] = gtf_df.loc[mask, "attribute"].str.extract(
                r'product "([^"]+)"'
            )[0]
    else:
        gtf_df["gene_id"] = gtf_df["attribute"].str.extract(r'gene_id "([^"]+)"')[0]
        gtf_df["gene_name"] = gtf_df["attribute"].str.extract(r'gene_name "([^"]+)"')[0]

    gtf_df["gene_biotype"] = gtf_df["attribute"].str.extract(r'gene_biotype "([^"]+)"')[0]

    # Also try 'gene_type' (used by GENCODE)
    mask_biotype = gtf_df["gene_biotype"].isna()
    if mask_biotype.any():
        gtf_df.loc[mask_biotype, "gene_biotype"] = gtf_df.loc[mask_biotype, "attribute"].str.extract(
            r'gene_type "([^"]+)"'
        )[0]

    # --- Filter to requested features ------------------------------------
    gtf_df = gtf_df[gtf_df["feature"].isin(features)]
    print(f"  [GTF] Keeping features: {features} ({len(gtf_df):,} rows)")

    gtf_df = gtf_df.dropna(subset=["chromosome", "start", "end", "gene_id"])
    gtf_df["start"] = pd.to_numeric(gtf_df["start"], errors="coerce")
    gtf_df["end"] = pd.to_numeric(gtf_df["end"], errors="coerce")
    gtf_df = gtf_df.dropna(subset=["start", "end"])
    gtf_df["start"] = gtf_df["start"].astype(int)
    gtf_df["end"] = gtf_df["end"].astype(int)

    gtf_df["width"] = gtf_df["end"] - gtf_df["start"] + 1

    # Report stats
    n_genes = gtf_df["gene_id"].nunique()
    n_chrs = gtf_df["chromosome"].nunique()
    print(f"  [GTF] Loaded {n_genes:,} unique gene IDs across {n_chrs} chromosome(s).")

    return gtf_df[[
        "chromosome", "start", "end", "width", "strand",
        "gene_id", "gene_name", "gene_biotype"
    ]]


def annotate_regions(input_file, gtf_df, is_ncbi=False, window=0,
                     upstream=0, downstream=0):
    """
    Annotate genomic regions with overlapping genes from a GTF DataFrame.

    Parameters
    ----------
    input_file : str
        Path to input file with genomic regions.
    gtf_df : pd.DataFrame
        Parsed GTF data from sort_and_load_gtf().
    is_ncbi : bool
        Whether NCBI format is used.
    window : int
        Window size for position-only input (bp).
    upstream : int
        Extend search upstream of each region (bp).
    downstream : int
        Extend search downstream of each region (bp).

    Returns
    -------
    pd.DataFrame with annotated overlaps.
    """
    # --- Use the smart input parser --------------------------------------
    input_df = parse_input_regions(input_file, window=window)

    if input_df.empty:
        raise ValueError("Input file is empty after parsing.")

    # --- Apply upstream/downstream extension -----------------------------
    if upstream > 0 or downstream > 0:
        print(
            f"  [Annotate] Extending search: "
            f"upstream={upstream:,} bp, downstream={downstream:,} bp"
        )
        input_df["start"] = (input_df["start"] - upstream).clip(lower=0)
        input_df["end"] = input_df["end"] + downstream

    # --- Find overlaps ---------------------------------------------------
    overlaps = []
    total = len(input_df)
    report_interval = max(1, total // 10)  # report ~10 times

    for idx, (_, region) in enumerate(input_df.iterrows(), 1):
        if total > 50 and idx % report_interval == 0:
            pct = int((idx / total) * 100)
            print(f"  [Progress] Annotating region {idx:,}/{total:,} ({pct}%)...")

        chrom = region["chromosome"]
        start = region["start"]
        end = region["end"]

        match = gtf_df[
            (gtf_df["chromosome"] == str(chrom)) &
            (gtf_df["end"] >= start) &
            (gtf_df["start"] <= end)
        ]

        for _, gene in match.iterrows():
            # Calculate overlap statistics
            overlap_start = max(start, gene["start"])
            overlap_end = min(end, gene["end"])
            overlap_bp = max(0, overlap_end - overlap_start + 1)
            region_len = end - start + 1
            gene_len = gene["end"] - gene["start"] + 1
            overlap_pct_region = round((overlap_bp / region_len) * 100, 2) if region_len > 0 else 0.0
            overlap_pct_gene = round((overlap_bp / gene_len) * 100, 2) if gene_len > 0 else 0.0

            # Distance to gene TSS (transcription start site)
            tss = gene["start"] if gene.get("strand", "+") == "+" else gene["end"]
            region_mid = (start + end) // 2
            distance_to_tss = region_mid - tss

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
                "gene_biotype": gene["gene_biotype"],
                "overlap_bp": overlap_bp,
                "overlap_pct_region": overlap_pct_region,
                "overlap_pct_gene": overlap_pct_gene,
                "distance_to_tss": distance_to_tss,
            })

    # --- Handle no-overlap case ------------------------------------------
    if not overlaps:
        input_chrs = set(input_df["chromosome"].astype(str).unique())
        gtf_chrs = set(gtf_df["chromosome"].astype(str).unique())
        common_chrs = input_chrs.intersection(gtf_chrs)

        if not common_chrs:
            raise ValueError(
                f"\n  No overlapping genes found AND no matching chromosome names.\n"
                f"    Input file chromosomes (sample): {sorted(input_chrs)[:10]}\n"
                f"    GTF file chromosomes (sample):   {sorted(gtf_chrs)[:10]}\n"
                f"    Possible cause: input uses plain numbers (1,2,3) but GTF uses "
                f"NCBI IDs (e.g., NC_037328.1).\n"
                f"    Fix: Adjust input chromosome names to match GTF, or use -ncbi if "
                f"not already."
            )
        else:
            print(
                f"  [Info] No overlapping genes found, but chromosome names DO match.\n"
                f"    Common chromosomes: {sorted(common_chrs)[:10]}\n"
                f"    Your regions may fall outside annotated gene boundaries.\n"
                f"    Try increasing --upstream / --downstream or --window."
            )
        return pd.DataFrame()

    result_df = pd.DataFrame(overlaps)
    result_df["gene_name"] = result_df["gene_name"].fillna("NA").replace("", "NA")

    print(f"  [Annotate] Found {len(result_df):,} gene overlaps.")

    return result_df


# ---------------------------------------------------------------------------
# Standalone run_annotate (also called from main.py)
# ---------------------------------------------------------------------------
def run_annotate(args):
    """Execute the annotate subcommand."""
    print("\n  Sorting and parsing GTF file...")

    # Parse feature filter
    features = None
    if hasattr(args, "feature") and args.feature:
        features = [f.strip() for f in args.feature.split(",")]
        print(f"  [Config] Feature filter: {features}")

    gtf_df = sort_and_load_gtf(args.gtf, is_ncbi=args.ncbi, features=features)

    print("  Annotating regions...")

    # Get upstream/downstream
    upstream = getattr(args, "upstream", 0) or 0
    downstream = getattr(args, "downstream", 0) or 0
    window = getattr(args, "window", 0) or 0

    annotated_df = annotate_regions(
        args.input, gtf_df,
        is_ncbi=args.ncbi,
        window=window,
        upstream=upstream,
        downstream=downstream,
    )

    if annotated_df.empty:
        print("  No overlapping genes found.")
        return

    annotated_df.to_csv(args.output, sep='\t', index=False)

    total_regions = parse_input_regions(args.input, window=window).shape[0]
    total_overlaps = annotated_df.shape[0]

    print(f"\n  Annotated output saved to {args.output}")
    print(f"\n  Total input regions: {total_regions}")
    print(f"  Total overlapping genes: {total_overlaps}")

    if 'gene_name' in annotated_df.columns:
        unique_genes_df = annotated_df.drop_duplicates(subset=['gene_name'])
        total_unique_genes = unique_genes_df.shape[0]

        print(f"  Total unique genes: {total_unique_genes}")

        # Biotype breakdown
        if 'gene_biotype' in annotated_df.columns:
            biotype_counts = annotated_df['gene_biotype'].value_counts()
            print(f"\n  Gene biotype distribution:")
            for bt, count in biotype_counts.head(10).items():
                print(f"    {bt}: {count}")

        # Write summary text file
        summary_file = os.path.splitext(args.output)[0] + "_summary.txt"
        with open(summary_file, 'w') as f:
            f.write("QGAT Annotation Summary\n")
            f.write("=======================\n")
            f.write(f"Input file: {args.input}\n")
            f.write(f"GTF file: {args.gtf}\n")
            f.write(f"Output file: {args.output}\n")
            if window > 0:
                f.write(f"Window size: {window:,} bp\n")
            if upstream > 0:
                f.write(f"Upstream extension: {upstream:,} bp\n")
            if downstream > 0:
                f.write(f"Downstream extension: {downstream:,} bp\n")
            if features:
                f.write(f"Feature filter: {', '.join(features)}\n")
            f.write(f"Total input regions: {total_regions}\n")
            f.write(f"Total overlapping genes: {total_overlaps}\n")
            f.write(f"Total unique genes: {total_unique_genes}\n")
            if 'gene_biotype' in annotated_df.columns:
                f.write(f"\nGene biotype distribution:\n")
                for bt, count in annotated_df['gene_biotype'].value_counts().items():
                    f.write(f"  {bt}: {count}\n")

        print(f"\n  Summary saved to {summary_file}")

        # Save unique genes table
        unique_file = os.path.splitext(args.output)[0] + "_unique.tsv"
        unique_genes_df.to_csv(unique_file, sep='\t', index=False)
        print(f"  Unique genes table saved to {unique_file}")

    elif not annotated_df.empty:
        print("  Warning: 'gene_name' column not found. Cannot compute unique genes.")
