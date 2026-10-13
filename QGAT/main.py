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
    Open a file, properly decompressing .gz if the path ends in .gz
    or if only the .gz version exists on disk. Always latin-1 for text.
    """
    p = Path(path)
    if str(p).endswith('.gz'):
        if p.exists():
            return gzip.open(p, 'rt', encoding='latin-1')
        raise FileNotFoundError(f"Compressed file not found: {p}")

    if p.exists():
        return open(p, 'r', encoding='latin-1')

    gz_path = Path(str(p) + ".gz")
    if gz_path.exists():
        print(f"  [File] Opening compressed file: {gz_path}")
        return gzip.open(gz_path, 'rt', encoding='latin-1')

    raise FileNotFoundError(f"Neither {p} nor {gz_path} found.")


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

    custom_db = getattr(args, "custom_db", None) or getattr(args, "db", None)

    if not args.species and not custom_db:
        raise ValueError(
            "Missing database source. You must specify EITHER:\n"
            "  (a) -species <cattle|pig|chicken|goat|sheep|horse> (for bundled Animal QTLdb)\n"
            "  (b) -db <path_to_file.bed|.gff> (to upload or supply your own custom QTL database)"
        )

    # 1. Load Database
    if custom_db:
        custom_path = Path(custom_db)
        if not custom_path.exists():
            gz_path = Path(str(custom_path) + ".gz")
            if gz_path.exists():
                custom_path = gz_path
            else:
                raise FileNotFoundError(f"Custom QTL database file not found: {custom_db}")

        species_label = f"Custom ({custom_path.name})"
        is_gff = args.trait or ".gff" in custom_path.name.lower()
        print(f"\n  [QTL] Mode: Custom Database -> {custom_path}")
        print(f"  [QTL] Format: {'GFF (trait-level)' if is_gff else 'BED (standard)'}")

        if is_gff:
            with open_auto(custom_path) as gff_file:
                qtl_df = utils.parse_qtl_gff(gff_file)
        else:
            with open_auto(custom_path) as bed_file:
                qtl_df = utils.parse_qtl_bed(bed_file)
    else:
        species_label = args.species
        print(f"\n  [QTL] Species: {args.species}")
        print(f"  [QTL] Mode: {'GFF (trait-level)' if args.trait else 'BED (standard)'}")

        if args.trait:
            gff_path = qtl_dir / gff_files[args.species]
            with open_auto(gff_path) as gff_file:
                qtl_df = utils.parse_qtl_gff(gff_file)
        else:
            bed_path = qtl_dir / bed_files[args.species]
            with open_auto(bed_path) as bed_file:
                qtl_df = utils.parse_qtl_bed(bed_file)

    print(f"  [QTL] Database loaded: {len(qtl_df):,} QTL entries")

    # 2. Parse input regions with smart chromosome alignment
    window = getattr(args, "window", 0) or 0
    if window > 0:
        print(f"  [QTL] Window: +/- {window:,} bp")

    input_df = utils.parse_input_regions_qtl(args.input, window=window, qtl_df=qtl_df)

    result_df = utils.find_overlapping_qtls(input_df, qtl_df)
    result_df.to_csv(args.output, sep='\t', index=False)
    print(f"\n  Results saved to {args.output}\n")

    # --- Summary ---------------------------------------------------------
    print("  Summary Report:")
    print(f"    Database: {species_label}")
    print(f"    Total input regions: {len(input_df)}")
    print(f"    Total overlapping QTLs: {len(result_df)}")

    if not result_df.empty:
        regions_with_hits = result_df[['chromosome', 'region_start', 'region_end']].drop_duplicates()
        print(f"    Regions with at least 1 QTL: {len(regions_with_hits)}/{len(input_df)}")

        if not args.trait and 'qtl_name' in result_df.columns:
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
# Plot subcommand — Publication-Ready Modern Aesthetics
# ---------------------------------------------------------------------------
def run_plot(args):
    if not Path(args.input).exists():
        raise FileNotFoundError(f"Input file not found: {args.input}")

    df = pd.read_csv(args.input, sep='\t')
    if df.empty:
        raise ValueError("Input file is empty.")
    if 'trait' not in df.columns:
        raise ValueError("Input file must contain a 'trait' column from GFF (3rd column).")

    # Aggregate counts
    trait_counts = df['trait'].value_counts()
    total_annotations = len(df)
    unique_traits = len(trait_counts)

    # Optional top N filter if user specified or default to top 20 for crisp visualization
    top_n = getattr(args, "top", None) or 20
    if len(trait_counts) > top_n:
        top_counts = trait_counts.head(top_n)
        other_count = trait_counts.iloc[top_n:].sum()
        if other_count > 0:
            top_counts["Other Traits"] = other_count
        plot_counts = top_counts
    else:
        plot_counts = trait_counts

    num_bars = len(plot_counts)
    fig_height = max(6.0, num_bars * 0.48 + 2.0)

    # Modern color palette: polished gradient from deep navy/indigo to teal
    import numpy as np
    colors = plt.cm.viridis(np.linspace(0.25, 0.85, num_bars))[::-1]

    fig, ax = plt.subplots(figsize=(11, fig_height), dpi=300)
    fig.patch.set_facecolor("#ffffff")
    ax.set_facecolor("#fafbfc")

    y_pos = np.arange(num_bars)
    # Clean formatted labels (replace underscores, title case)
    clean_labels = [str(l).replace("_", " ").title() for l in plot_counts.index[::-1]]
    bar_values = plot_counts.values[::-1]

    # Draw bars with slight rounded appearance, crisp borders
    bars = ax.barh(
        y_pos, bar_values,
        color=colors,
        edgecolor="#ffffff",
        linewidth=1.2,
        height=0.68,
        zorder=3
    )

    # Despine: clean modern scientific aesthetic
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_color("#94a3b8")
    ax.spines['bottom'].set_linewidth(1.0)

    # Light subtle vertical grid behind bars
    ax.grid(axis='x', color="#e2e8f0", linestyle='--', linewidth=0.7, alpha=0.9, zorder=0)

    # Numerical & percentage value annotations right next to each bar
    max_val = max(bar_values) if len(bar_values) > 0 else 1
    for i, bar in enumerate(bars):
        width = bar.get_width()
        pct = (width / total_annotations) * 100 if total_annotations > 0 else 0
        ax.text(
            width + (max_val * 0.015),
            bar.get_y() + bar.get_height() / 2,
            f"{width:,}  ({pct:.1f}%)",
            va='center',
            ha='left',
            fontsize=9.5,
            fontweight='600',
            color='#1e293b'
        )

    # Set axis limits & labels
    ax.set_xlim(0, max_val * 1.25)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(clean_labels, fontsize=10.5, fontweight='500', color='#0f172a')
    ax.tick_params(axis='y', length=0, pad=8)
    ax.tick_params(axis='x', colors='#475569', labelsize=10)

    # Clear, descriptive typography
    ax.set_xlabel("Number of Overlapping QTLs", fontsize=11.5, fontweight='600', color='#334155', labelpad=10)

    # Title & Subtitle
    plt.suptitle(
        "Distribution of Quantitative Trait Loci (QTLs) by Trait Class",
        fontsize=14.5,
        fontweight='700',
        color='#0f172a',
        x=0.08,
        y=0.98,
        ha='left'
    )
    ax.set_title(
        f"Total Associations: N = {total_annotations:,} | Distinct Trait Categories: {unique_traits:,} | Dataset: {Path(args.input).name}",
        fontsize=10.0,
        fontweight='400',
        color='#64748b',
        pad=16,
        loc='left'
    )

    plt.tight_layout()
    plt.subplots_adjust(top=0.90)
    plt.savefig(args.output, dpi=300, facecolor=fig.get_facecolor(), bbox_inches='tight')
    plt.close()
    print(f"\n  [Plot] Publication-ready trait distribution saved to: {args.output}")


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
        "qtl", help="Find overlapping QTLs from Animal QTLdb or custom database",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    qtl_parser.add_argument(
        "-i", "--input", required=True,
        help="Input file with genomic regions (TSV/CSV). "
             "Columns: chromosome+start+end OR chromosome+position (with --window)."
    )
    qtl_parser.add_argument(
        "-species", required=False, default=None,
        choices=["cattle", "chicken", "goat", "sheep", "horse", "pig"],
        help="Animal species (for pre-bundled Animal QTLdb files)"
    )
    qtl_parser.add_argument(
        "-db", "--db", "--qtl-db", dest="custom_db", default=None,
        help="Path to custom QTL database file (BED or GFF format, supports .gz). "
             "Allows uploading or supplying your own QTL database."
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
    plot_parser.add_argument(
        "--top", type=int, default=20,
        help="Display top N most frequent traits (default: 20, grouping remainder into 'Other Traits')"
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
