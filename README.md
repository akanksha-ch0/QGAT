QGAT: QTL and Gene Annotation Tool
==================================

QGAT is a Python-based command-line tool designed to identify Quantitative Trait Loci (QTLs) and annotate genomic regions with gene information. It supports major livestock species and allows users to input genomic coordinates and obtain QTL or gene overlaps using internal or user-provided data.

Features
--------
- QTL overlap detection from Animal QTLdb (for six species).
- Support for both BED-based and GFF-based QTL lookup (`--trait` flag).
- Gene annotation using user-supplied GTF files.
- Support for both Ensembl and NCBI GTF formats.
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


Note:
-----
The `QTLdb/` folder containing species-specific `.bed` files is included in the package and required for QTL search (except when using `--trait` which uses GFF).

The `input/` and `output/` folders are also included for organizing example data and results.

Usage
-----
Use `QGAT --help` to see available commands and flags.

### Subcommand: QTL 

Find overlapping QTLs for genomic regions using BED files.

    QGAT qtl -i path/to/input.tsv --species goat -o path/to/output_qtls.tsv

Arguments:
- `-i`, `--input`     : Input file (TXT,TSV and CSV) with columns: `chromosome`, `start`, `end`.
- `--species`         : One of `cattle`, `goat`, `sheep`, `pig`, `chicken`, `horse`.
- `-o`, `--output`    : Output TSV file.

### Trait-Specific QTL Search (from GFF)

If you have a GFF file (with trait information), use:

    QGAT qtl -i path/to/input.tsv --species cattle --trait -o path/to/trait_qtl_output.tsv

- The output includes full trait names and significance (P-value) from the GFF.
- All columns from the GFF entry will be retained in the output.

### Subcommand: Annotate

Annotate input regions with overlapping genes from a GTF file.

    QGAT annotate -i path/to/input.tsv --gtf path/to/genome.gtf -o path/to/annotated.tsv

For NCBI GTF format, add:

    QGAT annotate -i path/to/input.tsv --gtf path/to/ncbi_genome.gtf --ncbi -o path/to/annotated.tsv

Arguments:
- `-i`, `--input`     : Input TSV with `chromosome`, `start`, `end` columns.
- `--gtf`             : GTF file path (Ensembl or NCBI format).
- `--ncbi`            : Optional flag to indicate NCBI format.
- `-o`, `--output`    : Output file with gene annotations.

### Subcommand: Plot

Plot trait frequency (only for GFF-based QTL search using `--trait`):

    QGAT plot -i path/to/trait_qtl_output.tsv -o path/to/trait_plot.png

Arguments:
- `-i`, `--input`     : Output file from `--trait` QTL run.
- `-o`, `--output`    : Output image file (e.g. PNG).

Output Files
------------
- `output_qtls.tsv`     : Overlapping QTLs from BED-based QTLdb.
- `trait_qtl_output.tsv`: GFF-based trait-overlap results.
- `annotated.tsv`       : Gene annotations from GTF file.
- `trait_plot.png`      : Trait frequency bar plot.

Important Notes
---------------
### ⚠ Chromosome Naming Consistency

Ensure that chromosome names in your input file match those used in the GTF/GFF/QTLdb files.

Examples:
- If GTF uses `NC_037328.1`, use the same in your input.
- If GTF uses `1`, input should use `1` (not `Chr.1` or other variations).

If mismatched, either:
- No overlaps will be found, or
- The program will raise a warning or error.

Help Command
------------
Use the help flag to see available options:

    QGAT --help

    QGAT qtl --help
    QGAT annotate --help
    QGAT plot --help

License
-------
NBAGR_KU
