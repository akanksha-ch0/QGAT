import os
import gzip
import pandas as pd
import csv
import re
import sys
from collections import Counter
from io import StringIO

# ---------------------------------------------------------------------------
# Column-name aliases (case-insensitive matching)
# ---------------------------------------------------------------------------
COLUMN_ALIASES = {
    "chromosome": [
        "chromosome", "chr", "chrom", "#chrom", "#chr",
        "contig", "seqname", "seqnames",
    ],
    "start": [
        "start", "begin", "chromstart", "region_start", "startpos",
        "start_pos", "start_position",
    ],
    "end": [
        "end", "stop", "chromend", "region_end", "endpos",
        "end_pos", "end_position",
    ],
    "position": [
        "position", "pos", "bp", "loc", "location",
        "snp_pos", "snp_position", "base_pair", "basepair",
    ],
}


def _build_alias_map():
    """Return a dict mapping every alias (lower-cased) → canonical name."""
    alias_map = {}
    for canonical, aliases in COLUMN_ALIASES.items():
        for alias in aliases:
            alias_map[alias.lower()] = canonical
    return alias_map


_ALIAS_MAP = _build_alias_map()


# ---------------------------------------------------------------------------
# Separator detection
# ---------------------------------------------------------------------------
def detect_separator(filepath, n_lines=10):
    """
    Sniff the first *n_lines* of a file and return the most likely
    column separator: tab, comma, or whitespace.
    """
    with open(filepath, "r", encoding="utf-8-sig") as fh:
        head = ""
        for i, line in enumerate(fh):
            if i >= n_lines:
                break
            head += line

    # Count occurrences
    tabs = head.count("\t")
    commas = head.count(",")
    # Rough heuristic – tabs almost always win in genomic files
    if tabs >= commas and tabs > 0:
        return "\t"
    if commas > tabs and commas > 0:
        return ","
    # Fallback: whitespace
    return r"\s+"


# ---------------------------------------------------------------------------
# Column normalisation
# ---------------------------------------------------------------------------
def normalize_columns(columns):
    """
    Map a list of raw column names to their canonical equivalents.
    Returns (mapped_columns, mapping_dict) where mapping_dict records
    what was renamed.

    Raises ValueError with a helpful message if no chromosome column
    is found.
    """
    mapped = []
    mapping = {}
    seen_canonical = set()

    for col in columns:
        clean = col.strip().lstrip("\ufeff")  # strip BOM if present
        key = clean.lower()
        canonical = _ALIAS_MAP.get(key)
        if canonical and canonical not in seen_canonical:
            mapped.append(canonical)
            mapping[col] = canonical
            seen_canonical.add(canonical)
        else:
            mapped.append(clean)

    if "chromosome" not in seen_canonical:
        raise ValueError(
            f"Could not identify a chromosome column in your input file.\n"
            f"  Detected columns: {columns}\n"
            f"  Accepted names (case-insensitive): "
            f"{COLUMN_ALIASES['chromosome']}\n"
            f"  Please rename your chromosome column to one of the above."
        )

    return mapped, mapping


# ---------------------------------------------------------------------------
# Main input parser
# ---------------------------------------------------------------------------
def parse_input_regions(input_path, window=0):
    """
    Read a genomic-region file and return a DataFrame with columns
    [chromosome, start, end] plus any extra columns from the original file.

    Supported input styles
    ----------------------
    1. **Range input** – columns: chromosome, start, end
    2. **Position input** – columns: chromosome, position
       Requires *window > 0*; creates start = pos - window, end = pos + window.

    Parameters
    ----------
    input_path : str or Path
        Path to the input file (TSV / CSV / space-delimited).
    window : int
        Half-window size in bp.  Only used when the input has a *position*
        column but no *start* / *end*.

    Returns
    -------
    pd.DataFrame  with at least [chromosome, start, end].
    """
    input_path = str(input_path)

    # --- detect separator ------------------------------------------------
    sep = detect_separator(input_path)

    # --- read file -------------------------------------------------------
    try:
        df = pd.read_csv(
            input_path,
            sep=sep,
            dtype=str,
            encoding="utf-8-sig",   # handles BOM
            engine="python",        # needed for regex sep
            skipinitialspace=True,
            comment="#",            # skip comment lines
        )
    except Exception as e:
        raise FileNotFoundError(
            f"Could not read input file '{input_path}': {e}"
        )

    if df.empty:
        raise ValueError(f"Input file '{input_path}' is empty.")

    # Strip whitespace from column names
    df.columns = [c.strip() for c in df.columns]

    # --- normalise column names ------------------------------------------
    original_cols = list(df.columns)
    new_cols, mapping = normalize_columns(original_cols)
    df.columns = new_cols

    if mapping:
        renamed = ", ".join(f"'{k}' -> '{v}'" for k, v in mapping.items())
        print(f"  [Input] Auto-mapped columns: {renamed}")

    # --- determine input mode --------------------------------------------
    has_start = "start" in df.columns
    has_end = "end" in df.columns
    has_position = "position" in df.columns

    if has_start and has_end:
        # ---- Range mode -------------------------------------------------
        mode = "range"
        df["start"] = pd.to_numeric(df["start"], errors="coerce")
        df["end"] = pd.to_numeric(df["end"], errors="coerce")

        bad_rows = df[df["start"].isna() | df["end"].isna()]
        if not bad_rows.empty:
            print(
                f"  [Warning] Dropped {len(bad_rows)} row(s) with "
                f"non-numeric start/end values."
            )
            df = df.dropna(subset=["start", "end"])

        df["start"] = df["start"].astype(int)
        df["end"] = df["end"].astype(int)

        # If window is provided in range mode, extend the regions
        if window > 0:
            print(
                f"  [Input] Extending regions by +/- {window:,} bp "
                f"(--window applied to range input)."
            )
            df["start"] = (df["start"] - window).clip(lower=0)
            df["end"] = df["end"] + window

    elif has_position:
        # ---- Position mode ----------------------------------------------
        mode = "position"
        df["position"] = pd.to_numeric(df["position"], errors="coerce")

        bad_rows = df[df["position"].isna()]
        if not bad_rows.empty:
            print(
                f"  [Warning] Dropped {len(bad_rows)} row(s) with "
                f"non-numeric position values."
            )
            df = df.dropna(subset=["position"])

        df["position"] = df["position"].astype(int)

        if window <= 0:
            raise ValueError(
                "Your input has a 'position' column but no 'start'/'end'.\n"
                "  You must provide --window SIZE (in bp) to create regions\n"
                "  around each position.  Example:  --window 50000\n"
                "  This creates regions [position - 50000, position + 50000]."
            )

        print(
            f"  [Input] Position mode: creating regions of +/- {window:,} bp "
            f"around each position."
        )
        df["start"] = (df["position"] - window).clip(lower=0)
        df["end"] = df["position"] + window

    else:
        raise ValueError(
            f"Cannot determine genomic coordinates from your input.\n"
            f"  Detected columns: {original_cols}\n\n"
            f"  Your file must have EITHER:\n"
            f"    (a) 'chromosome', 'start', 'end' columns   (range input)\n"
            f"    (b) 'chromosome', 'position' columns        (position input, requires --window)\n\n"
            f"  Accepted column names (case-insensitive):\n"
            f"    Chromosome: {COLUMN_ALIASES['chromosome']}\n"
            f"    Start:      {COLUMN_ALIASES['start']}\n"
            f"    End:        {COLUMN_ALIASES['end']}\n"
            f"    Position:   {COLUMN_ALIASES['position']}"
        )

    # --- validate coordinates --------------------------------------------
    invalid = df[df["start"] > df["end"]]
    if not invalid.empty:
        print(
            f"  [Warning] {len(invalid)} region(s) have start > end. "
            f"Swapping start and end for those rows."
        )
        mask = df["start"] > df["end"]
        df.loc[mask, ["start", "end"]] = df.loc[mask, ["end", "start"]].values

    negative = df[df["start"] < 0]
    if not negative.empty:
        print(
            f"  [Warning] {len(negative)} region(s) have negative start. "
            f"Clipping to 0."
        )
        df["start"] = df["start"].clip(lower=0)

    # --- deduplicate with warning ----------------------------------------
    n_before = len(df)
    df = df.drop_duplicates(subset=["chromosome", "start", "end"])
    n_dropped = n_before - len(df)
    if n_dropped > 0:
        print(
            f"  [Warning] Removed {n_dropped} duplicate region(s)."
        )

    # --- chromosome name cleanup -----------------------------------------
    df["chromosome"] = df["chromosome"].astype(str).str.strip()

    print(
        f"  [Input] Loaded {len(df)} region(s) in {mode} mode "
        f"from '{os.path.basename(input_path)}'."
    )

    return df


# ---------------------------------------------------------------------------
# QTL input parser — uses the new smart parser but adds Chr. prefix for QTLdb
# ---------------------------------------------------------------------------
def parse_input_regions_qtl(input_path, window=0):
    """
    Parse input regions for the QTL subcommand.
    Adds the 'Chr.' prefix required by QTLdb BED/GFF files.
    """
    df = parse_input_regions(input_path, window=window)
    df["chromosome"] = df["chromosome"].apply(
        lambda x: f"Chr.{x}" if not str(x).startswith("Chr.") else x
    )
    return df


# ---------------------------------------------------------------------------
# QTL BED parser
# ---------------------------------------------------------------------------
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
            try:
                start = int(row[1])
                end = int(row[2])
            except ValueError:
                continue  # skip malformed rows

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


# ---------------------------------------------------------------------------
# QTL GFF parser
# ---------------------------------------------------------------------------
def parse_qtl_gff(gff_source):
    """
    Parse a GFF file.
    Accepts either:
      - a file path (str, Path, bytes)
      - or an open file-like object
    Supports both plain .gff and .gff.gz if given as a path.
    """
    gff_data = []

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
            try:
                start = int(parts[3])
                end = int(parts[4])
            except ValueError:
                continue  # skip malformed rows

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


# ---------------------------------------------------------------------------
# QTL overlap finder
# ---------------------------------------------------------------------------
def find_overlapping_qtls(input_df, qtl_df):
    """Find overlapping QTLs between input regions and QTL database."""
    results = []

    total = len(input_df)
    for idx, (_, region) in enumerate(input_df.iterrows(), 1):
        if total > 100 and idx % 500 == 0:
            print(f"  [Progress] Processed {idx:,}/{total:,} regions...")

        overlaps = qtl_df[
            (qtl_df['chromosome'] == region['chromosome']) &
            (qtl_df['qtl_end'] >= region['start']) &
            (qtl_df['qtl_start'] <= region['end'])
        ]
        for _, qtl in overlaps.iterrows():
            # Calculate overlap statistics
            overlap_start = max(region['start'], qtl['qtl_start'])
            overlap_end = min(region['end'], qtl['qtl_end'])
            overlap_bp = max(0, overlap_end - overlap_start + 1)
            region_len = region['end'] - region['start'] + 1
            overlap_pct = round((overlap_bp / region_len) * 100, 2) if region_len > 0 else 0.0

            result = {
                'chromosome': region['chromosome'],
                'region_start': region['start'],
                'region_end': region['end'],
                'qtl_start': qtl['qtl_start'],
                'qtl_end': qtl['qtl_end'],
                'qtl_name': qtl['qtl_name'],
                'qtl_number': qtl['qtl_number'],
                'overlap_bp': overlap_bp,
                'overlap_pct': overlap_pct,
            }
            for col in qtl.index:
                if col not in result:
                    result[col] = qtl[col]
            results.append(result)

    if not results:
        # Diagnostic: check chromosome name mismatch
        input_chrs = set(input_df["chromosome"].unique())
        qtl_chrs = set(qtl_df["chromosome"].unique())
        common = input_chrs & qtl_chrs
        if not common:
            print(
                f"\n  [Warning] No overlapping QTLs found AND no matching chromosome names.\n"
                f"    Input chromosomes: {sorted(input_chrs)[:10]}\n"
                f"    QTL DB chromosomes: {sorted(qtl_chrs)[:10]}\n"
                f"    Possible cause: chromosome naming mismatch.\n"
                f"    The QTL database uses 'Chr.X' format. Your input names are auto-prefixed,\n"
                f"    but if you used full NCBI accessions, they won't match."
            )

    return pd.DataFrame(results)


# ---------------------------------------------------------------------------
# Summary generator
# ---------------------------------------------------------------------------
def generate_summary(input_df, result_df):
    """Print a summary of QTL overlap results."""
    total_regions = len(input_df)
    total_qtls = len(result_df)
    qtl_counter = Counter()
    for qtl in result_df['qtl_name']:
        if isinstance(qtl, str):
            qtl_counter[qtl] += 1

    top_qtls = qtl_counter.most_common(10)
    print("\n Summary Report:")
    print(f"  Total input regions: {total_regions}")
    print(f"  Total overlapping QTLs: {total_qtls}")
    if top_qtls:
        print("  Top 10 QTLs:")
        for name, count in top_qtls:
            print(f"    {name} ({count})")
