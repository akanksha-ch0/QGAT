QGAT: QTL and Gene Annotation Tool
==================================

**Web Application (Live)**: [https://akanksha-ch0.github.io/QGAT/](https://akanksha-ch0.github.io/QGAT/)

QGAT is a Python-based (version 3.7 or higher) command-line tool and interactive web application designed to identify Quantitative Trait Loci (QTLs) and annotate genomic regions with gene information. It supports major livestock species and allows users to input genomic coordinates and obtain QTL or gene overlaps using internal or user-provided data.

Features
--------
- QTL overlap detection from Animal QTLdb (for six species).
- Support for both BED-based and GFF-based QTL lookup (`-trait` flag).
- Gene annotation using user-supplied GTF files.
- Support for both Ensembl and NCBI GTF formats.
- **NEW: `--window` flag** — use position-only input (chr + position) by specifying a window size.
- **NEW: Flexible input parsing** — auto-detects column names (case-insensitive), file separators (TSV/CSV), and handles Excel BOM.
- **NEW: `--feature` filter** — restrict annotation to specific GTF feature types (e.g., `gene` only).
- **NEW: `--upstream` / `--downstream`** — extend search beyond exact region boundaries.
- **NEW: Overlap statistics** — overlap_bp, overlap_pct, distance_to_tss in output.
- Automatic GTF sorting and formatting.
- Trait-wise bar plot generation.
- Clean tabular output and summary report.
- Subcommand interface for clear usage.

Supported Species for QTL Search:
---------------------------------
- Cattle
- Goat
- Sheep
- Pig
- Chicken
- Horse

Installation
------------
1. Download: https://github.com/akanksha-ch0/QGAT.git
    
2. Install:  
Open your terminal, navigate to the downloaded project folder, and install the package using `pip`: pip install ./QGAT-main

3. OR directly install in command-line using : pip install git+https://github.com/akanksha-ch0/QGAT.git

4. Check version:

       QGAT --version


Note:
-----
The `QTLdb/` folder containing species-specific `.bed` files is included in the package and required for QTL search (except when using `-trait` which uses GFF).

The `input/` and `output/` folders are also included for organizing example data and results. The input file can be in TXT, TSV and CSV format.

Input Format
------------
QGAT accepts **two input formats**. Column names are **case-insensitive** and many aliases are accepted.

### Format 1: Range Input (chromosome + start + end)

    chromosome	start	end
    1	22834818	23381127
    2	38590	3831532

### Format 2: Position-Only Input (chromosome + position)

    chromosome	position
    1	22834818
    2	38590

When using Format 2, you **must** provide `--window SIZE` to create regions around each position:

    QGAT annotate -i snps.tsv -gtf genome.gtf -o out.tsv --window 50000

This creates regions: `[position - 50000, position + 50000]`

### Accepted Column Names

| Canonical | Accepted Aliases (case-insensitive) |
|-----------|-------------------------------------|
| chromosome | chr, chrom, CHR, #CHROM, contig, seqname |
| start | START, begin, chromStart, region_start |
| end | END, stop, chromEnd, region_end |
| position | pos, POS, bp, BP, loc, location, snp_pos |

Usage
-----
Use `QGAT --help` to see available commands and flags.

### Subcommand: QTL 

Find overlapping QTLs for genomic regions using BED files.

    QGAT qtl -i path/to/input.tsv -species goat -o path/to/output_qtls.tsv

**With position-only input and window:**

    QGAT qtl -i path/to/snps.tsv -species goat -o path/to/output_qtls.tsv --window 50000

Arguments:
- `-i`, `--input`     : Input file with columns: `chromosome`, `start`, `end` (or `chromosome`, `position` with `--window`).
- `-species`         : One of `cattle`, `goat`, `sheep`, `pig`, `chicken`, `horse`.
- `-o`, `--output`    : Output TSV file.
- `-w`, `--window`    : Window size (bp) for position-only input. Creates regions `[pos - window, pos + window]`.

### Trait-Specific QTL Search (from GFF)

use: --trait flag

    QGAT qtl -i path/to/input.tsv -species cattle -trait -o path/to/trait_qtl_output.tsv

- The output includes full trait names and significance (P-value) from the GFF.
- All columns from the GFF entry will be retained in the output.


### Subcommand: Plot

Plot trait frequency (only for GFF-based QTL search using `-trait`):

    QGAT plot -i path/to/trait_qtl_output.tsv -o path/to/trait_plot.png

Arguments:
- `-i`, `--input`     : Output file from `-trait` QTL run.
- `-o`, `--output`    : Output image file (e.g. PNG).


### Subcommand: Annotate

Annotate input regions with overlapping genes from a GTF file.

    QGAT annotate -i path/to/input.tsv -gtf path/to/genome.gtf -o path/to/annotated.tsv

**With position-only input:**

    QGAT annotate -i path/to/snps.tsv -gtf path/to/genome.gtf -o path/to/annotated.tsv --window 100000

**With upstream/downstream extension:**

    QGAT annotate -i path/to/input.tsv -gtf path/to/genome.gtf -o path/to/annotated.tsv --upstream 5000 --downstream 5000

**Filter to specific GTF features:**

    QGAT annotate -i path/to/input.tsv -gtf path/to/genome.gtf -o path/to/annotated.tsv --feature gene

For NCBI GTF format, add:

    QGAT annotate -i path/to/input.tsv -gtf path/to/ncbi_genome.gtf -ncbi -o path/to/annotated.tsv

Arguments:
- `-i`, `--input`     : Input file with `chromosome`, `start`, `end` columns (or `chromosome`, `position` with `--window`).
- `-gtf`              : GTF file path (Ensembl or NCBI format).
- `-ncbi`             : Optional flag to indicate NCBI format.
- `-o`, `--output`    : Output file with gene annotations.
- `-w`, `--window`    : Window size (bp) for position-only input.
- `--feature`         : Comma-separated GTF feature types to keep (default: gene,transcript,exon). Example: `--feature gene`
- `--upstream`        : Extend search upstream of each region (bp). Example: `--upstream 5000`
- `--downstream`      : Extend search downstream of each region (bp). Example: `--downstream 5000`


### Gene Ontology & Functional Enrichment Pipeline

Extract unique candidate gene symbols from `annotated.tsv` for downstream functional enrichment analysis:

    # Extract unique gene symbols
    cut -f 10 path/to/annotated.tsv | sort -u > candidate_genes.txt

These candidate genes can be directly piped into functional enrichment tools (**DAVID**, **g:Profiler**, **PANTHER**, **ShinyGO**) or submitted with 1-click via the **QGAT Web Application**.


Output Files
------------
- `output_qtls.tsv`     : Overlapping QTLs from BED-based QTLdb (includes overlap_bp, overlap_pct).
- `trait_qtl_output.tsv`: GFF-based trait-overlap results.
- `trait_plot.png`      : Trait frequency bar plot.
- `annotated.tsv`       : Gene annotations from GTF files (includes overlap_bp, overlap_pct_region, overlap_pct_gene, distance_to_tss).
- `annotated_summary.txt`: Summary statistics.
- `annotated_unique.tsv` : Deduplicated gene list.


Output Columns (Annotate)
--------------------------
| Column | Description |
|--------|-------------|
| chromosome | Input region chromosome |
| start | Input region start |
| end | Input region end |
| gene_chr | Gene chromosome |
| gene_start_pos | Gene start position |
| gene_end_pos | Gene end position |
| width | Gene width (bp) |
| strand | Gene strand (+/-) |
| gene_id | Gene identifier |
| gene_name | Gene name/symbol |
| gene_biotype | Gene biotype (e.g., protein_coding) |
| overlap_bp | Overlap length in bp |
| overlap_pct_region | % of input region overlapping the gene |
| overlap_pct_gene | % of gene overlapping the input region |
| distance_to_tss | Distance from region midpoint to gene TSS |


Important Notes
---------------
### Chromosome Naming Consistency

Ensure that chromosome names in your input file match those used in the GTF/GFF/QTLdb files.

Examples:
- If GTF uses `NC_037328.1`, use the same in your input.
- If GTF uses `1`, input should use `1` (not `Chr.1` or other variations).

If mismatched, either:
- No overlaps will be found, or
- The program will raise a warning with diagnostic information showing which chromosome names were found in each file.

Help Command
------------
Use the help flag to see available options:

    QGAT --help

    QGAT qtl --help
    QGAT annotate --help
    QGAT plot --help

License
-------
NBAGR
