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
        print(f"Opening compressed file: {gz_path}")
        return gzip.open(gz_path, 'rt', encoding='latin-1')
    raise FileNotFoundError(f"Neither {path} nor {gz_path} found.")

def clean_qtl_name(name):
    """Remove QTL ID in brackets from QTL name, e.g., 'Milk fat QTL (1234)' → 'Milk fat QTL'"""
    return re.sub(r"\s*\(.*?\)", "", name).strip()

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

    input_df = utils.parse_input_regions(args.input)

    if args.trait:
        gff_path = qtl_dir / gff_files[args.species]
        with open_auto(gff_path) as gff_file:
            qtl_df = utils.parse_qtl_gff(gff_file)
    else:
        bed_path = qtl_dir / bed_files[args.species]
        with open_auto(bed_path) as bed_file:
            qtl_df = utils.parse_qtl_bed(bed_file)

    result_df = utils.find_overlapping_qtls(input_df, qtl_df)
    result_df.to_csv(args.output, sep='\t', index=False)
    print(f"\n Results saved to {args.output}\n")

    print(" Summary Report:")
    print(f" Total input regions: {len(input_df)}")
    print(f" Total overlapping QTLs: {len(result_df)}")

    if not result_df.empty and not args.trait:
        result_df["qtl_trait"] = result_df["qtl_name"].apply(clean_qtl_name)
        top_qtls = result_df["qtl_trait"].value_counts().head(10)
        print(" Top 10 QTLs (by trait):")
        for name, count in top_qtls.items():
            print(f"  {name} ({count})")



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
    print("\n📊 Trait plot saved to: {}".format(args.output))


def main():
    parser = argparse.ArgumentParser(description="QGAT: QTL and Gene Annotation Tool")
    subparsers = parser.add_subparsers(dest='command', required=True)

    qtl_parser = subparsers.add_parser("qtl", help="Find overlapping QTLs")
    qtl_parser.add_argument("-i", "--input", required=True, help="Input TSV file with regions")
    qtl_parser.add_argument("--species", required=True, choices=["cattle", "chicken", "goat", "sheep", "horse", "pig"], help="Animal species")
    qtl_parser.add_argument("-o", "--output", required=True, help="Output TSV file")
    qtl_parser.add_argument("--trait", action="store_true", help="Use GFF file to extract trait-level QTLs")

    annot_parser = subparsers.add_parser("annotate", help="Annotate genomic regions using GTF")
    annot_parser.add_argument("-i", "--input", required=True, help="Input TSV file with regions")
    annot_parser.add_argument("--gtf", required=True, help="GTF file for gene annotation")
    annot_parser.add_argument("-o", "--output", required=True, help="Output TSV file with gene annotations")
    annot_parser.add_argument("--ncbi", action="store_true", help="Indicate if the GTF file is from NCBI format")

    plot_parser = subparsers.add_parser("plot", help="Plot trait distribution from GFF-based QTL extraction")
    plot_parser.add_argument("-i", "--input", required=True, help="TSV file output from --trait-based QTL extraction")
    plot_parser.add_argument("-o", "--output", required=True, help="Output plot image file (PNG)")

    args = parser.parse_args()

    try:
        if args.command == "qtl":
            run_qtl(args)
        elif args.command == "annotate":
            run_annotate(args)
        elif args.command == "plot":
            run_plot(args)
    except Exception as e:
        print("\n {}".format(e))
        sys.exit(1)


if __name__ == "__main__":
    main()
