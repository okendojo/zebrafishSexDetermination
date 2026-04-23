#script used to extract gene intervals on sar4 block in the q-arm of chr4

GTF=GRCz12tu_genomic.gtf
awk -F'\t' 'BEGIN{OFS="\t"}
$1=="chr4" && $3=="gene" && $5 >= 72316977 && $4 <= 76164240 {
    gene="."
    if (match($9, /gene_name "([^"]+)"/, m)) gene=m[1];
    else if (match($9, /gene_id "([^"]+)"/, m)) gene=m[1];
    print $1, $4-1, $5, gene
}' "$GTF" | sort -k2,2n > chr4_sar4block_genes.bed
