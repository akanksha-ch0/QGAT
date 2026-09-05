import argparse
import sys
import pandas as pd
import re
import matplotlib.pyplot as plt
from pathlib import Path
import os
import gzip

from QGAT import utils
from QGAT.annotate import sort_and_load_gtf, annotate_regions


def open_auto(path):
    """
    Open a file, using .gz if needed. Always latin-1 for text.
    """
    if path.exists():
        return open(path, 'r', encoding='latin-1')
    gz_path = str(path) + ".gz"
    if os.path.exists(gz_path):
        print(f"  [File] Opening compressed file: {gz_path}")
        return gzip.open(gz_path, 'rt', encoding='latin-1')
    raise FileNotFoundError(f"Neither {path} nor {gz_path} found.")


def clean_qtl_name(name):
    """Remove QTL ID in brackets from QTL name, e.g., 'Milk fat QTL (1234)' -> 'Milk fat QTL'"""
    return re.sub(r"\s*\(.*?\)", "", name).strip()


# ---------------------------------------------------------------------------
# QTL subcommand
# ---------------------------------------------------------------------------
def run_qtl(args):
    base_dir = Path(__file__).resolve().parent
    qtl_dir = base_dir / "QTLdb"

    bed_files = {
        "cattle": "QTLdb_cattleARS_UCD1.bed",
        "chicken": "QTLdb_chickenGG4.bed",
        "goat": "QTLdb_goatCHIR_1.bed",
        "sheep": "QTLdb_sheepOAR3.bed",
        "horse": "QTLdb_horseEC2.bed",
        "pig": "QTLdb_pigMARC1.bed"
    }

    gff_files = {
        "cattle": "QTLdb_cattleARS_UCD1.gff",
        "chicken": "QTLdb_chickenGG4.gff",
        "goat": "QTLdb_goatCHIR_1.gff",
        "sheep": "QTLdb_sheepOAR3.gff",
        "horse": "QTLdb_horseEC2.gff",
        "pig": "QTLdb_pigMARC1.gff"
    }

    print(f"\n  [QTL] Species: {args.species}")
    print(f"  [QTL] Mode: {'GFF (trait-level)' if args.trait else 'BED (standard)'}")

    # Use the new smart input parser with QTL-specific Chr. prefixing
    window = getattr(args, "window", 0) or 0
    if window > 0:
        print(f"  [QTL] Window: +/- {window:,} bp")

    input_df = utils.parse_input_regions_qtl(args.input, window=window)

    if args.trait:
        gff_path = qtl_dir / gff_files[args.species]
        with open_auto(gff_path) as gff_file:
            qtl_df = utils.parse_qtl_gff(gff_file)
    else:
        bed_path = qtl_dir / bed_files[args.species]
        with open_auto(bed_path) as bed_file:
            qtl_df = utils.parse_qtl_bed(bed_file)

    print(f"  [QTL] Database loaded: {len(qtl_df):,} QTL entries")

    result_df = utils.find_overlapping_qtls(input_df, qtl_df)
    result_df.to_csv(args.output, sep='\t', index=False)
    print(f"\n  Results saved to {args.output}\n")

    # --- Summary ---------------------------------------------------------
    print("  Summary Report:")
    print(f"    Species: {args.species}")
    print(f"    Total input regions: {len(input_df)}")
    print(f"    Total overlapping QTLs: {len(result_df)}")

    if not result_df.empty:
        regions_with_hits = result_df[['chromosome', 'region_start', 'region_end']].drop_duplicates()
        print(f"    Regions with at least 1 QTL: {len(regions_with_hits)}/{len(input_df)}")

        if not args.trait:
            result_df["qtl_trait"] = result_df["qtl_name"].apply(clean_qtl_name)
            top_qtls = result_df["qtl_trait"].value_counts().head(10)
            print("    Top 10 QTLs (by trait):")
            for name, count in top_qtls.items():
                print(f"      {name} ({count})")


# ---------------------------------------------------------------------------
# Annotate subcommand
# ---------------------------------------------------------------------------
def run_annotate(args):
    print("\n  Sorting and parsing GTF file...")

    # Parse feature filter
    features = None
    if hasattr(args, "feature") and args.feature:
        features = [f.strip() for f in args.feature.split(",")]
        print(f"  [Config] Feature filter: {features}")

    gtf_df = sort_and_load_gtf(args.gtf, is_ncbi=args.ncbi, features=features)

    print("  Annotating regions...")

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

    total_regions = utils.parse_input_regions(args.input, window=window).shape[0]
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


# ---------------------------------------------------------------------------
# Plot subcommand
# ---------------------------------------------------------------------------
def run_plot(args):
    if not Path(args.input).exists():
        raise FileNotFoundError("Input file not found: {}".format(args.input))

    df = pd.read_csv(args.input, sep='\t')
    if df.empty:
        raise ValueError("Input file is empty.")
    if 'trait' not in df.columns:
        raise ValueError("Input file must contain a 'trait' column from GFF (3rd column).")

    trait_counts = df['trait'].value_counts()

    num_traits = len(trait_counts)
    fig_height = max(6, num_traits * 0.4)

    plt.figure(figsize=(12, fig_height))
    bars = plt.barh(trait_counts.index[::-1], trait_counts.values[::-1], color='steelblue', edgecolor='black')

    for bar in bars:
        width = bar.get_width()
        plt.text(width + 1, bar.get_y() + bar.get_height() / 2, str(width), va='center', fontsize=8)

    plt.xlabel("QTL Count", fontsize=12)
    plt.ylabel("Trait", fontsize=12)
    plt.title("QTL Trait Distribution", fontsize=14, weight='bold')
    plt.tight_layout()
    plt.savefig(args.output, dpi=300)
    print("\n  Trait plot saved to: {}".format(args.output))


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------
def main():
    from QGAT import __version__

    parser = argparse.ArgumentParser(
        description="QGAT: QTL and Gene Annotation Tool",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # QTL search with range input (chr, start, end):
  QGAT qtl -i regions.tsv -species cattle -o results.tsv

  # QTL search with position-only input (chr, position) + window:
  QGAT qtl -i snps.tsv -species cattle -o results.tsv --window 50000

  # Gene annotation with GTF:
  QGAT annotate -i regions.tsv -gtf genome.gtf -o annotated.tsv

  # Gene annotation with position-only input + window:
  QGAT annotate -i snps.tsv -gtf genome.gtf -o annotated.tsv --window 100000

  # Gene annotation with upstream/downstream extension:
  QGAT annotate -i regions.tsv -gtf genome.gtf -o annotated.tsv --upstream 5000 --downstream 5000

  # Filter to only gene features:
  QGAT annotate -i regions.tsv -gtf genome.gtf -o annotated.tsv --feature gene

  # Trait-level QTL search + plot:
  QGAT qtl -i regions.tsv -species cattle -trait -o traits.tsv
  QGAT plot -i traits.tsv -o trait_plot.png

Input Format:
  Your input file can use EITHER format:
    Format 1 (ranges):    chromosome  start    end
    Format 2 (positions): chromosome  position
  Column names are case-insensitive. Aliases like chr/CHR/chrom/pos/bp/BP all work.
  When using Format 2, you MUST provide --window SIZE.
        """
    )
    parser.add_argument(
        "--version", action="version",
        version=f"QGAT v{__version__}"
    )

    subparsers = parser.add_subparsers(dest='command', required=True)

    # --- QTL subcommand --------------------------------------------------
    qtl_parser = subparsers.add_parser(
        "qtl", help="Find overlapping QTLs from Animal QTLdb",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    qtl_parser.add_argument(
        "-i", "--input", required=True,
        help="Input file with genomic regions (TSV/CSV). "
             "Columns: chromosome+start+end OR chromosome+position (with --window)."
    )
    qtl_parser.add_argument(
        "-species", required=True,
        choices=["cattle", "chicken", "goat", "sheep", "horse", "pig"],
        help="Animal species"
    )
    qtl_parser.add_argument("-o", "--output", required=True, help="Output TSV file")
    qtl_parser.add_argument(
        "-trait", action="store_true",
        help="Use GFF file to extract trait-level QTLs (includes trait type and P-value)"
    )
    qtl_parser.add_argument(
        "-w", "--window", type=int, default=0,
        help="Window size in bp around each position. "
             "Use when input has 'position' column instead of 'start'/'end'. "
             "Creates regions: [position - window, position + window]. "
             "Can also extend existing start/end regions."
    )

    # --- Annotate subcommand ---------------------------------------------
    annot_parser = subparsers.add_parser(
        "annotate", help="Annotate genomic regions with genes using a GTF file",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    annot_parser.add_argument(
        "-i", "--input", required=True,
        help="Input file with genomic regions (TSV/CSV). "
             "Columns: chromosome+start+end OR chromosome+position (with --window)."
    )
    annot_parser.add_argument("-gtf", required=True, help="GTF file for gene annotation")
    annot_parser.add_argument("-o", "--output", required=True, help="Output TSV file with gene annotations")
    annot_parser.add_argument(
        "-ncbi", action="store_true",
        help="Indicate that the GTF file is in NCBI format "
             "(uses db_xref for gene_id, gene for gene_name)"
    )
    annot_parser.add_argument(
        "-w", "--window", type=int, default=0,
        help="Window size in bp around each position. "
             "Use when input has 'position' column instead of 'start'/'end'. "
             "Creates regions: [position - window, position + window]. "
             "Can also extend existing start/end regions."
    )
    annot_parser.add_argument(
        "--feature", type=str, default=None,
        help="Comma-separated list of GTF feature types to keep. "
             "Default: gene,transcript,exon. Example: --feature gene "
             "or --feature gene,transcript"
    )
    annot_parser.add_argument(
        "--upstream", type=int, default=0,
        help="Extend search upstream of each region by this many bp. "
             "Useful for finding nearby genes. Example: --upstream 5000"
    )
    annot_parser.add_argument(
        "--downstream", type=int, default=0,
        help="Extend search downstream of each region by this many bp. "
             "Example: --downstream 5000"
    )

    # --- Plot subcommand -------------------------------------------------
    plot_parser = subparsers.add_parser(
        "plot", help="Plot trait distribution from GFF-based QTL extraction"
    )
    plot_parser.add_argument(
        "-i", "--input", required=True,
        help="TSV file output from -trait-based QTL extraction"
    )
    plot_parser.add_argument(
        "-o", "--output", required=True,
        help="Output plot image file (PNG)"
    )

    args = parser.parse_args()

    try:
        if args.command == "qtl":
            run_qtl(args)
        elif args.command == "annotate":
            run_annotate(args)
        elif args.command == "plot":
            run_plot(args)
    except FileNotFoundError as e:
        print(f"\n  [Error] File not found: {e}", file=sys.stderr)
        sys.exit(1)
    except ValueError as e:
        print(f"\n  [Error] {e}", file=sys.stderr)
        sys.exit(1)
    except KeyboardInterrupt:
        print("\n  [Interrupted] Operation cancelled by user.")
        sys.exit(130)
    except Exception as e:
        print(f"\n  [Error] Unexpected error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
