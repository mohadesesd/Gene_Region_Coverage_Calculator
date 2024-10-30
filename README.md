# Gene Region Coverage Calculator

This repository provides a Python script to calculate the coverage percentage of exon and intron regions of specified genes based on a BED file. It uses a GFF file for gene annotations and a list of gene symbols to analyze.

---

## Setup

### Requirements

- **Python**: Version 3.6 or higher
- **Packages**: Install using the following command:
  ```bash
  pip install pandas pybedtools argparse
## Requirements

### BEDTools
BEDTools must be installed and accessible in your system PATH.

### File Requirements

- **GFF file**: A GFF file for gene annotations (e.g., Gencode).
- **BED file**: A BED file with regions of interest (e.g., `output.bed`).
- **Gene List**: A list of gene symbols (e.g., `TP53`, `GAPDH`, `CDH1`).

---

## Usage

### Clone the repository and navigate to it:
```bash
git clone https://github.com/yourusername/repository-name.git
cd repository-name
```
### Run the script:

```bash
python calculate_coverage.py --bed_file path/to/output.bed --gff_file path/to/gencode.v47lift37.basic.annotation.gff3 --gene_list TP53,GAPDH,CDH1
```
## Output

The script will generate a CSV file named `coverage_results.csv` with coverage data for exons and introns in each specified gene.

## Arguments

| Argument      | Description                    | Example Value                                    |
|---------------|--------------------------------|--------------------------------------------------|
| `--bed_file`  | Path to the BED file           | `path/to/output.bed`                             |
| `--gff_file`  | Path to the GFF annotation file| `path/to/gencode.v47lift37.basic.annotation.gff3`|
| `--gene_list` | Comma-separated list of gene symbols | `TP53,GAPDH,CDH1`                               |

## Example Output

The output CSV file `coverage_results.csv` contains:

| gene  | region_type | region_start | region_end | coverage_percent |
|-------|-------------|--------------|------------|-------------------|
| TP53  | exon        | 7590694      | 7590808    | 95.3              |
| GAPDH | intron      | 1005632      | 1005890    | 88.0              |
## Notes

- Ensure that the BED and GFF files align to the same reference genome (e.g., hg19).
- The script only considers canonical transcripts from the GFF file.

