import os
import gzip
import pandas as pd
import csv
import re
from collections import Counter


def parse_input_regions(input_path):
    df = pd.read_csv(input_path, sep='\t', dtype=str)
    expected_cols = ['chromosome', 'start', 'end']
    if not all(col in df.columns for col in expected_cols):
        raise ValueError(f"Input file must contain columns: {expected_cols}")
    df['chromosome'] = df['chromosome'].apply(
        lambda x: f"Chr.{x}" if not str(x).startswith("Chr.") else x
    )
    df['start'] = df['start'].astype(int)
    df['end'] = df['end'].astype(int)
    return df


def parse_qtl_bed(bed_source):
    """
    Parse a BED file.
    Accepts either:
      - a file path (str, Path, bytes)
      - or an open file-like object
    """
    qtl_data = []

    if isinstance(bed_source, (str, bytes, os.PathLike)):
        if str(bed_source).endswith('.gz'):
            f = gzip.open(bed_source, 'rt', encoding='latin-1')
        else:
            f = open(bed_source, 'r', encoding='latin-1')
        should_close = True
    else:
        f = bed_source
        should_close = False

    try:
        reader = csv.reader(f, delimiter='\t')
        for row in reader:
            if len(row) < 4:
                continue
            chrom = row[0]
            start = int(row[1])
            end = int(row[2])

            qtl_name_parts = []
            qtl_number = None
            for val in row[3:]:
                qtl_name_parts.append(val)
                match = re.match(r"\((\d+)\)", val)
                if match:
                    qtl_number = match.group(1)
                    break

            qtl_name = " ".join(qtl_name_parts).strip()
            qtl_data.append({
                'chromosome': chrom,
                'qtl_start': start,
                'qtl_end': end,
                'qtl_name': qtl_name,
                'qtl_number': qtl_number
            })
    finally:
        if should_close:
            f.close()

    return pd.DataFrame(qtl_data)


def parse_qtl_gff(gff_source):
    """
    Parse a GFF file.
    Accepts either:
      - a file path (str, Path, bytes)
      - or an open file-like object
    Supports both plain .gff and .gff.gz if given as a path.
    """
    gff_data = []

    # If user passed a path, auto-handle .gz if needed:
    if isinstance(gff_source, (str, bytes, os.PathLike)):
        if str(gff_source).endswith('.gz'):
            f = gzip.open(gff_source, 'rt', encoding='latin-1')
        else:
            f = open(gff_source, 'r', encoding='latin-1')
        should_close = True
    else:
        f = gff_source
        should_close = False

    try:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.strip().split('\t')
            if len(parts) != 9:
                continue

            chrom = parts[0]
            trait_type = parts[2]
            start = int(parts[3])
            end = int(parts[4])
            info = parts[8]

            info_dict = dict(
                item.split('=') for item in info.strip(';').split(';') if '=' in item
            )

            gff_entry = {
                'chromosome': chrom,
                'qtl_start': start,
                'qtl_end': end,
                'qtl_name': info_dict.get("Name", ""),
                'qtl_number': info_dict.get("QTL_ID", ""),
                'trait': trait_type,
                'p_value': info_dict.get("P-value", "")
            }

            gff_data.append(gff_entry)
    finally:
        if should_close:
            f.close()

    return pd.DataFrame(gff_data)


def find_overlapping_qtls(input_df, qtl_df):
    results = []

    for _, region in input_df.iterrows():
        overlaps = qtl_df[
            (qtl_df['chromosome'] == region['chromosome']) &
            (qtl_df['qtl_end'] >= region['start']) &
            (qtl_df['qtl_start'] <= region['end'])
        ]
        for _, qtl in overlaps.iterrows():
            result = {
                'chromosome': region['chromosome'],
                'region_start': region['start'],
                'region_end': region['end'],
                'qtl_start': qtl['qtl_start'],
                'qtl_end': qtl['qtl_end'],
                'qtl_name': qtl['qtl_name'],
                'qtl_number': qtl['qtl_number']
            }
            for col in qtl.index:
                if col not in result:
                    result[col] = qtl[col]
            results.append(result)

    return pd.DataFrame(results)


def generate_summary(input_df, result_df):
    total_regions = len(input_df)
    total_qtls = len(result_df)
    qtl_counter = Counter()
    for qtl in result_df['qtl_name']:
        if isinstance(qtl, str):
            qtl_counter[qtl] += 1

    top_qtls = qtl_counter.most_common(10)
    print("\n📊 Summary Report:")
    print(f"Total input regions: {total_regions}")
    print(f"Total overlapping QTLs: {total_qtls}")
    print("Top 10 QTLs:")
    for name, count in top_qtls:
        print(f"  {name} ({count})")
