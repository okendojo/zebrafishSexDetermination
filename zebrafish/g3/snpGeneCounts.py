import pandas as pd

# Load haplotype information into a DataFrame
hap_df = pd.read_csv('haplotype_info.txt', sep='\t', header=None, names=['CHROM', 'POS', 'HAP', 'REF', 'ALT'])

# Assuming you have a separate file with gene information (gene_bed_file.bed)
gene_bed = pd.read_csv('gene_bed_file.bed', sep='\t', header=None, names=['CHROM', 'START', 'END', 'GENE'])

# Merge haplotype and gene information based on chromosome
merged_df = pd.merge(hap_df, gene_bed, on='CHROM')

# Filter SNPs within gene boundaries
gene_snps = merged_df[(merged_df['POS'] >= merged_df['START']) & (merged_df['POS'] <= merged_df['END'])]

# Count SNPs per haplotype and gene
snp_counts = gene_snps.groupby(['GENE', 'HAP']).size().reset_index(name='SNP_COUNT')

print(snp_counts)

