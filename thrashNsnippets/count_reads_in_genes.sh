#!/bin/bash

source activate emma_env

fragments_file=$1

prefix=$(echo $fragments_file | sed s/fragments.tsv//)

sort -k1,1 -k2,2 -k 3,3 $fragments_file > ${prefix}fragments.sort.bed

# bedmap --delim '\t' --echo --echo-map-id ${prefix}.fragments.sort.bed annotations/hg38_genes.bed > ${prefix}.genes_bc.bed
# bedmap --delim '\t' --echo --echo-map-id ${prefix}.fragments.sort.bed annotations/hg38_genes.bed > ${prefix}.promoters_bc.bed

bedtools intersect -a annotations/hg38_genes_2.bed -b ${prefix}fragments.sort.bed -wb | bedtools groupby -g 4 -c 8 -o collapse > ${prefix}genes_bc.bed
bedtools intersect -a annotations/hg38_promoters_500flank_2.bed -b ${prefix}fragments.sort.bed -wb | bedtools groupby -g 4 -c 8 -o collapse > ${prefix}promoters_bc.bed

