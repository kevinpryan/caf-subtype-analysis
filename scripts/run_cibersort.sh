#!/bin/bash

while getopts "m:r:p:" flag 
do
    case "${flag}" in
        m) mixture=${OPTARG};;
	r) ref=${OPTARG};;
	p) pheno=${OPTARG};;
    esac
done

echo "mixture file: $mixture";
echo "reference file: $ref";
echo "phenotype classes file: $pheno"

docker run -v /home/kevin/Documents/PhD/subtypes/caf-subtype-analysis/intermediate_files/cibersort:/src/data \
	   -v /home/kevin/Documents/PhD/subtypes/caf-subtype-analysis/intermediate_files/cibersort/outputs_2022-10-07:/src/outdir \
	   cibersortx/fractions \
	   --username k.ryan45@nuigalway.ie \
	   --token b7f03b943ade9b4146dc2126b4ac9d19 \
	   --refsample $ref \
	   --mixture $mixture \
	   --rmbatchBmode FALSE \
           --phenoclasses $pheno \
	   --perm 500 \
	   --verbose TRUE \
	   --G.min 300 \
           --G.max 500
