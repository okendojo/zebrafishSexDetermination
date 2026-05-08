import pandas as pd
from scipy.stats import pearsonr

# Load gene expression data
gene_expression = pd.read_csv('gene_expression.csv', index_col=0)

# Load annotated SNP data
annotated_snps = pd.read_csv('annotated_snps.vcf', delimiter='\t', comment='#', header=None)

# Merge the datasets based on gene IDs or positions
merged_data = pd.merge(gene_expression, annotated_snps, how='inner', left_on='gene_id', right_on='gene_id')

# Perform correlation analysis
correlation_coefficient, p_value = pearsonr(merged_data['gene_expression'], merged_data['snp_genotype'])

print(f'Correlation Coefficient: {correlation_coefficient}')
print(f'P-value: {p_value}')

