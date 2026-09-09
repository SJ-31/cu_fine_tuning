#!/usr/bin/env bash

dbsnp="/data/project/stemcell/shannc/reference/variants/dbSNP_renamed.vcf.gz"
exprs="(INFO/CLNSIG ~ '2' || INFO/CLNSIG ~ '3') && (INFO/CLNREVSTAT ~ 'mult' || INFO/CLNREVSTAT ~ 'exp' || INFO/CLNREVSTAT ~ 'practice') && INFO/VC == 'SNV'"
output="/data/project/stemcell/shannc/repos/evo2_fine_tune/raw/dbSNP_benign.tsv"
query="%INFO/VC\t%INFO/GENEINFO\t%INFO/FREQ\t%INFO/CLNHGVS\t%INFO/CLNSIG"

echo -e "vc\tgeneinfo\tfreq\thgvs\tclinsig" >"${output}"
bcftools filter -i "${exprs}" "${dbsnp}" |
	bcftools query -f "${query}" >>"${output}"
