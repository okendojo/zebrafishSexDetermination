


gatk SelectVariants  -R hg38.fa.  --variant dna_variants.vcf.bgz --restrict-alleles-to BIALLELIC -select 'vc.getHetCount()==1' --select-type-to-include SNP -O dna_variants.selected.vcf.bgz

bcftools norm --rm-dup all dna_variants.selected.vcf.bgz | bgzip > out.vcf.gz
