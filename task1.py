import pandas as pd
import argparse
from pybedtools import BedTool
import re

def parse_gene_list(gene_list_str):
    # Replace commas and tabs with spaces, then split by whitespace
    return [gene.strip() for gene in gene_list_str.replace(',', ' ').replace('\t', ' ').split()]

def parse_gff(gff_file, gene_list):
    gene_data = []
    with open(gff_file, 'r') as file:
        for line in file:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            if fields[2] == "exon":
                # Extract attributes using the provided logic
                attributes = {item.split('=')[0]: item.split('=')[1] for item in fields[8].split(';') if '=' in item}
                gene_name = attributes.get("gene_name") or attributes.get("Name")
                if (gene_name in gene_list) and ('canonical' in attributes.get("tag")):
                    gene_data.append({
                        "chrom": fields[0],
                        "start": int(fields[3]) - 1,  # GFF is 1-based, BED is 0-based
                        "end": int(fields[4]),
                        "type": fields[2],
                        "gene": gene_name
                    })
    return pd.DataFrame(gene_data)

def calculate_introns(exons_df):
    intron_data = []

    # Group by gene
    grouped_exons = exons_df.groupby('gene')

    for gene_name, group in grouped_exons:
        group = group.sort_values('start')
        # Calculate introns by subtracting the end of one exon from the start of the next
        for i in range(len(group) - 1):
            intron_start = group.iloc[i]['end']
            intron_end = group.iloc[i + 1]['start']
            if intron_start < intron_end:  # Ensure valid intron
                intron_data.append({
                    'chrom': group.iloc[i]['chrom'],
                    'start': intron_start,
                    'end': intron_end,
                    'type': 'intron',
                    'gene': gene_name  # Keep the gene name for identification
                })

    return pd.DataFrame(intron_data)

def calculate_coverage(bed_file, regions_df):
    bed = BedTool(bed_file)

    coverage_results = []

    for _, region in regions_df.iterrows():
        # Create a BedTool object for the current region
        region_bed = BedTool([[region['chrom'], region['start'], region['end']]]).saveas()

        # Calculate coverage
        covered = region_bed.coverage(bed)
        coverage_fraction = 0
        #for i in covered:
        coverage_fraction += float(covered[0][-1]) * 100  # Convert to percentage
        coverage_results.append({
            'gene': region['gene'],
            'region_type': region['type'],
            'region_start': region['start'],
            'region_end': region['end'],
            'coverage_percent': min(coverage_fraction, 100)
        })

    return pd.DataFrame(coverage_results)

def main():
    parser = argparse.ArgumentParser(description='Calculate gene region coverage from BED and GFF files.')
    parser.add_argument('--bed_file', required=True, help='Path to the BED file.')
    parser.add_argument('--gff_file', required=True, help='Path to the GFF file.')
    parser.add_argument('--gene_list', required=True, help="List of gene symbols (e.g., TP53,CDH1,GAPDH).")
    args = parser.parse_args()
    gene_list = parse_gene_list(args.gene_list)

    # Parse GFF file to get exons
    exons_df = parse_gff(args.gff_file, gene_list)

    # Calculate introns from exons
    introns_df = calculate_introns(exons_df)

    # Calculate coverage for exons and introns
    exon_coverage = calculate_coverage(args.bed_file, exons_df)
    intron_coverage = calculate_coverage(args.bed_file, introns_df)

    # Combine results
    coverage_df = pd.concat([exon_coverage, intron_coverage], ignore_index=True)

    # Save DataFrame to a CSV file
    coverage_df.to_csv("coverage_results.csv", index=False, sep='\t')
    print("Coverage results saved to 'coverage_results.csv'")

if __name__ == "__main__":
    main()
