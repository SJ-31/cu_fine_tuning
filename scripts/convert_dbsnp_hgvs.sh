#!/usr/bin/env bash

from_api=data/processed/dbSNP_benign_formatted.csv
if [[ ! -e "${from_api}" ]]; then
	tmp=data/raw/dbSNP_hgvsg.txt
	xan select hgvs data/raw/dbSNP_benign.tsv |
		xan enum | xan explode hgvs --sep , |
		xan filter "'>' in hgvs" |
		xan groupby index 'first(hgvs)' |
		xan select 'first(hgvs)' |
		xan behead >"${tmp}"

	variant_validator_api.py -i "${tmp}" \
		-o data/processed/dbSNP_benign_formatted.csv
	rm "${tmp}"

fi

af=data/processed/dbSNP_benign_af.csv
if [[ ! -e "${af}" ]]; then
	xan select id,freq data/raw/dbSNP_benign.tsv |
		xan explode freq |
		xan separate freq ":" --into study,freq |
		xan map "max(slice(filter(split(freq, ','), n => len(n) > 1), 1, 10)) as alt_af" |
		xan drop freq |
		xan fill -v 0 -s alt_af |
		xan filter 'alt_af >0.05' |
		xan groupby id 'first(alt_af) as af' >"${af}"
fi
