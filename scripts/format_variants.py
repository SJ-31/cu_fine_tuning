#!/usr/bin/env ipython

import os
from collections.abc import Callable
from pathlib import Path
from typing import Literal

import polars as pl
from pyhere import here

RAW: Path = here("data", "raw")


# TODO: [2026-08-03 Mon] this script should generate the
# HGVSc codes for all the ProteinGym samples,
# Also generate a metadata file with
# extract ClinVar
# metadata including clinical significance,
# the source of the variant

# Columns are source, transcript_id, variant_class, symbol,
# disease, consequence, hgvs, clinsig

MAPPING: pl.DataFrame = (
    pl.read_csv(here("data", "mart_2026-08-03_filtered.csv"))
    .select(
        [
            "Transcript stable ID version",
            "RefSeq match transcript (MANE Select)",
            "Gene name",
        ]
    )
    .rename(
        {
            "Gene name": "symbol",
            "RefSeq match transcript (MANE Select)": "transcript_id",
            "Transcript stable ID version": "ensembl_transcript_id",
        }
    )
)


def map_ids(
    df: pl.DataFrame,
    on: Literal["transcript_id", "ensembl_transcript_id"],
    target: Literal["transcript_id", "symbol"],
) -> pl.DataFrame:
    assert on in df.columns, f"Column used for `on` {on} must be in df"
    mapping = MAPPING.select([on, target]).unique(on)
    return df.join(mapping, how="left", on=on)


def convert_ensembl_hgvsc(df: pl.DataFrame, col: str = "hgvs") -> pl.DataFrame:
    struct_names = ["ensembl_transcript_id", "hgvs_tmp"]
    df = (
        df.with_columns(
            pl.col(col)
            .str.split_exact(":", 1)
            .struct.rename_fields(struct_names)
            .alias("fields")
        )
        .unnest("fields")
        .drop(col)
    )
    if "transcript_id" not in df.columns:
        df = map_ids(df, on="ensembl_transcript_id", target="transcript_id")
    df = (
        df.filter(pl.col("transcript_id").is_not_null())
        .with_columns(hgvs=pl.col("transcript_id") + ":" + pl.col("hgvs_tmp"))
        .drop("hgvs_tmp")
    )
    return df


SOURCES = {
    "ProteinGym-snps-clinvar": here(
        RAW, "ProteinGym", "substitutions_preprocessed.csv"
    ),
    "ProteinGym-indels-gnomAD": here(
        RAW, "ProteinGym", "indels_preprocessed_gnomad.csv"
    ),  # TODO: replace ENST ids with MANE
    "ProteinGym-indels-clinvar": here(
        RAW, "ProteinGym", "indels_preprocessed_clinvar.csv"
    ),  # TODO: extract HGVSc from name column
    "ClinGen": here(RAW, "clingen_erepo-tabbed.tsv"),
    "CiVIC": here(RAW, "CIViC", "nightly-civic_accepted_civic_2026-08-11.txt"),
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


def format_proteingym_snps_clinvar(file) -> pl.DataFrame:
    df: pl.DataFrame = pl.read_csv(file, infer_schema_length=None)
    df = (
        df.filter(pl.col("CLNREVSTAT") != "criteria_provided,_single_submitter")
        .select(["HGVSc", "clinsig", "Consequence", "Feature", "CLNDN"])
        .rename(
            {
                "Consequence": "consequence",
                "Feature": "transcript_id",
                "HGVSc": "hgvs",
                "CLNDN": "disease",
            }
        )
        .with_columns(
            pl.col("consequence").str.replace_all(",", ";"),
            pl.col("disease")
            .str.split("|")
            .list.filter(pl.element() != "not_provided")
            .list.filter(pl.element() != "not_specified"),
        )
    )
    df = map_ids(df, "transcript_id", "symbol")
    return df


def format_proteingym_indels_gnomad(file) -> pl.DataFrame:
    df: pl.DataFrame = pl.read_csv(file, infer_schema_length=None)
    df = (
        df.select(["SYMBOL", "HGVSc"])
        .with_columns(
            pl.lit(None).alias("consequence"),
            pl.lit("benign").alias("clinsig"),
        )
        .rename({"SYMBOL": "symbol", "HGVSc": "hgvs"})
    )
    df = convert_ensembl_hgvsc(df)
    return df

def format_proteingym_indels_clinvar(file) -> pl.DataFrame:
    df: pl.DataFrame = pl.read_csv(file, infer_schema_length=None).filter(
        ~pl.col("Review status").is_in(
            ["criteria provided, single submitter", "no assertion criteria provided"]
        )
    )
    df = (
        df.select(
            ["Name", "refseq_id", "gene", "Clinical significance", "Condition(s)"]
        )
        .rename(
            {
                "Condition(s)": "disease",
                "Clinical significance": "clinsig",
                "gene": "symbol",
                "refseq_id": "transcript_id",
            }
        )
        .with_columns(
            pl.col("Name")
            .str.split_exact(":", 1)
            .struct.rename_fields(["_tmp", "hgvs"])
            .alias("fields")
        )
        .unnest("fields")
        .drop("_tmp")
        .with_columns(
            pl.col("hgvs").str.replace(" \\(p.*\\)", ""),
            pl.col("disease").str.split("|"),
            pl.lit(None).alias("consequence"),
        )
        .with_columns(hgvs=pl.col("transcript_id") + ":" + pl.col("hgvs"))
        .drop("Name")
    )
    return df

