#!/usr/bin/env ipython


import os
from pathlib import Path

import polars as pl
from pyhere import here

raw: Path = here("data", "raw")
remote: Path = here("data")


# TODO: [2026-08-03 Mon] this script should generate the
# HGVSc codes for all the ProteinGym samples,
# Also generate a metadata file with
# extract ClinVar
# metadata including clinical significance,
# the source of the variant


sources = {
    "ProteinGym-snps-clinvar": here(
        raw, "ProteinGym", "substitutions_preprocessed.csv"
    ),
    "ProteinGym-indels-gnomAD": here(
        raw, "ProteinGym", "indels_preprocessed_gnomad.csv"
    ),  # TODO: replace ENST ids with MANE
    "ProteinGym-indels-clinvar": here(
        raw, "ProteinGym", "indels_preprocessed_clinvar.csv"
    ),  # TODO: extract HGVSc from name column
    "ClinGen": here(raw, "clingen_erepo-tabbed.tsv"),
    "CiVIC": here(raw, "CIViC", "nightly-civic_accepted_civic_2026-08-11.txt"),
}

# [2026-08-11 Tue] format each of these separately

civic_columns = [
    "Allele",
    "Consequence",
    "SYMBOL",
    "Entrez Gene ID",
    "Feature_type",
    "Feature",
    "HGVSc",
    "HGVSp",
    "CIViC Variant Name",
    "CIViC Variant ID",
    "CIViC Variant Aliases",
    "CIViC Variant URL",
    "CIViC Molecular Profile Name",
    "CIViC Molecular Profile ID",
    "CIViC Molecular Profile Aliases",
    "CIViC Molecular Profile URL",
    "CIViC HGVS",
    "Allele Registry ID",
    "ClinVar IDs",
    "CIViC Molecular Profile Score",
    "CIViC Entity Type",
    "CIViC Entity ID",
    "CIViC Entity URL",
    "CIViC Entity Source",
    "CIViC Entity Variant Origin",
    "CIViC Entity Status",
    "CIViC Entity Significance",
    "CIViC Entity Direction",
    "CIViC Entity Disease",
    "CIViC Entity Therapies",
    "CIViC Entity Therapy Interaction Type",
    "CIViC Evidence Phenotypes",
    "CIViC Evidence Level",
    "CIViC Evidence Rating",
    "CIViC Assertion ACMG Codes",
    "CIViC Assertion AMP Category",
    "CIViC Assertion ClinGen Codes",
    "CIViC Assertion NCCN Guideline",
    "CIViC Assertion Regulatory Approval",
    "CIViC Assertion FDA Companion Test",
]

# Get NM_ hgvsc from the CIViC HGVS column
# civic = pl.read_csv(, separator = "|", new_columns = civic_columns, truncate_ragged_lines = True)
