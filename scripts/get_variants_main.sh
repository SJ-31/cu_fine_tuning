#!/usr/bin/env bash

# Call from project root

python main.py \
	--seqrepo data/seqrepo/2024-12-20 \
	--database data/all_seqs.db \
	--hgvs_column hgvs \
	--id_col transcript_id \
	--input data/processed/passing_variants.csv \
	--alias_spec scripts/aliases.yml \
	--load_sequences scripts/sequences.yml \
	--output_passed data/processed/generated_passed.csv \
	--output_failed data/processed/generated_failed.csv
