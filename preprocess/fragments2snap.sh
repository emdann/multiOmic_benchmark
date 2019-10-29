#!/bin/bash

fragments_tsv=$1
sample_prefix=$(echo $fragments_tsv | sed 's/.fragments.tsv//')

source activate emma_env

## Sort tsv file by barcode
sort -k4,4 $fragments_tsv > ${sample_prefix}.fragments.bed
## Compress bed file
gzip ${sample_prefix}.fragments.bed
## Run snaptools
snaptools snap-pre \
	--input-file=${sample_prefix}.fragments.bed.gz \
	--output-snap=${sample_prefix}.snap \
	--genome-name=hg38 \
	--genome-size=hg38.chrom.sizes \
	--min-mapq=30  \
	--min-flen=50  \
	--max-flen=1000  \
	--keep-chrm=TRUE  \
	--keep-single=FALSE  \
	--keep-secondary=False  \
	--overwrite=True  \
	--max-num=20000  \
	--verbose=True

snaptools snap-add-bmat \
	--snap-file ${sample_prefix}.snap \
	--bin-size-list 5000 10000 \
	--verbose=True

