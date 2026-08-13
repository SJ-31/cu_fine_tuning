#!/usr/bin/env ipython

from collections.abc import Callable
from pathlib import Path
from typing import Literal

import polars as pl
from hgvs.exceptions import HGVSParseError
from hgvs.parser import Parser
from pyhere import here

RAW: Path = here("data", "raw")


# Columns are source, transcript_id, variant_class, symbol,
# disease, consequence, hgvs, clinsig

MAPPING: pl.DataFrame = (
    pl.read_csv(here("data", "mart_2026-08-03_filtered.csv"))
    .select(
        [
            "Transcript stable ID version",
            "RefSeq match transcript (MANE Select)",
            "Gene name",
            "transcript_len",
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
MAPPING = pl.concat(
    [
        MAPPING,
        pl.read_csv(
            here("data", "mane_transcript_mapping_116_2026-08-12.tsv"), separator="\t"
        )
        .select(
            [
                "ensembl_transcript_id_version",
                "transcript_mane_select",
                "hgnc_symbol",
                "transcript_len",
            ]
        )
        .rename(
            {
                "hgnc_symbol": "symbol",
                "ensembl_transcript_id_version": "ensembl_transcript_id",
                "transcript_mane_select": "transcript_id",
            }
        )
        .filter(
            (pl.col("transcript_id").is_not_null())
            & (pl.col("ensembl_transcript_id").is_not_null())
        ),
    ],
    how="vertical_relaxed",
).unique("transcript_id")


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
        .drop("ensembl_transcript_id")
    )
    return df


def extract_transcript_id(col: str, df: pl.DataFrame) -> pl.DataFrame:
    struct_names = ["transcript_id", "hgvs_tmp"]
    df = (
        df.with_columns(
            pl.col(col)
            .str.split_exact(":", 1)
            .struct.rename_fields(struct_names)
            .alias("fields")
        )
        .unnest("fields")
        .drop("hgvs_tmp")
    )
    return df


SOURCES = {
    "ProteinGym-snps-clinvar": here(
        RAW, "ProteinGym", "substitutions_preprocessed.csv"
    ),
    "ProteinGym-indels-gnomAD": here(
        RAW, "ProteinGym", "indels_preprocessed_gnomad.csv"
    ),
    "ProteinGym-indels-clinvar": here(
        RAW, "ProteinGym", "indels_preprocessed_clinvar.csv"
    ),
    "ClinGen": here(RAW, "clingen_erepo-tabbed.tsv"),
    "CIViC": here(RAW, "CIViC", "nightly-civic_accepted_civic_2026-08-11.txt"),
    "COSMIC_resistance": here(
        RAW, "COSMIC", "Cosmic_ResistanceMutations_v104_GRCh38.tsv.gz"
    ),
    "COSMIC_census": here(RAW, "COSMIC", "Cosmic_MutantCensus_v104_GRCh38.tsv.gz"),
}


CIVIC_COLS = [
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


def consequence_from_hgvsp(
    df: pl.DataFrame, col: str = "prot", tmp_prefix: bool = True
) -> pl.DataFrame:
    parser = Parser()

    def get_consequence(hgvsp: str) -> str | None:
        try:
            var = parser.parse(hgvsp)
            if "alt" in dir(var.posedit.edit):
                type = var.posedit.edit.type
                ref = var.posedit.pos.start.aa
                alt = var.posedit.edit.alt
                if type == "sub" and alt != ref and alt != "*":
                    return "missense_variant"
                elif type == "sub" and alt == "*":
                    return "stop_gained"
                elif type == "identity":
                    return "synonymous_variant"
            return
        except HGVSParseError:
            return

    if tmp_prefix:
        df = df.with_columns(("TMP" + ":" + pl.col(col)).alias(col))
    df = df.with_columns(
        pl.col(col)
        .map_elements(get_consequence, return_dtype=pl.String)
        .alias("consequence")
    ).drop(col)
    return df


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
            .list.filter(pl.element() != "not_specified")
            .list.join(";"),
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
    df = extract_transcript_id("hgvs", df)
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
            pl.col("hgvs").str.extract("\\((p.*)\\)$").alias("prot"),
            pl.col("hgvs").str.replace(" \\(p.*\\)", ""),
            pl.col("disease").str.split("|").list.join(";"),
        )
        .with_columns(hgvs=pl.col("transcript_id") + ":" + pl.col("hgvs"))
        .drop("Name")
    )
    df = consequence_from_hgvsp(df)
    return df


def format_clingen_erepo(file) -> pl.DataFrame:
    df: pl.DataFrame = pl.read_csv(file, separator="\t", ignore_errors=True)
    df = (
        df.select(["HGNC Gene Symbol", "Disease", "Assertion", "Variation"])
        .rename(
            {"Disease": "disease", "Assertion": "clinsig", "HGNC Gene Symbol": "symbol"}
        )
        .with_columns(
            pl.col("Variation")
            .str.split_exact(":", 1)
            .struct.rename_fields(["transcript_id", "hgvs"])
        )
        .unnest("Variation")
        .with_columns(
            pl.col("hgvs").str.extract("\\((p.*)\\)$").alias("prot"),
            pl.col("transcript_id").str.replace("\\(.*\\)$", ""),
            pl.col("hgvs").str.replace(" \\(p.*\\)$", ""),
        )
        .with_columns(hgvs=pl.col("transcript_id") + ":" + pl.col("hgvs"))
    )
    df = consequence_from_hgvsp(df)
    return df


def format_civic(file) -> pl.DataFrame:
    df: pl.DataFrame = pl.read_csv(
        file, separator="|", new_columns=CIVIC_COLS, truncate_ragged_lines=True
    )
    df = (
        df.select(
            ["SYMBOL", "Consequence", "HGVSc", "CIViC Entity Disease", "CIViC HGVS"]
        )
        .rename(
            {
                "SYMBOL": "symbol",
                "Consequence": "consequence",
                "CIViC Entity Disease": "disease",
            }
        )
        .with_columns(
            pl.col("CIViC HGVS")
            .str.split("&")
            .list.filter(pl.element().str.starts_with("NM_"))
            .list.first()
            .alias("hgvs")
        )
    )
    has_hgvs = (
        df.filter(pl.col("hgvs").is_not_null())
        .with_columns(pl.col("hgvs").str.extract("(NM_.*):").alias("transcript_id"))
        .drop("HGVSc")
    )
    no_hgvs = df.filter(pl.col("hgvs").is_null()).drop("hgvs")
    no_hgvs = extract_transcript_id("HGVSc", no_hgvs).rename({"HGVSc": "hgvs"})
    return (
        pl.concat([has_hgvs, no_hgvs], how="diagonal_relaxed")
        .drop("CIViC HGVS")
        .with_columns(pl.col("consequence").str.replace_all("&", ";"))
    )


def format_cosmic(file) -> pl.DataFrame:
    cosmic_samples = (
        pl.read_csv(
            here(RAW, "COSMIC", "Cosmic_Sample_v104_GRCh38.tsv.gz"),
            separator="\t",
            infer_schema_length=None,
        )
        .select(["COSMIC_SAMPLE_ID", "TUMOUR_REMARK"])
        .rename({"TUMOUR_REMARK": "disease"})
    )
    df: pl.DataFrame = pl.read_csv(
        file, separator="\t", infer_schema_length=None
    ).filter(pl.col("MUTATION_SOMATIC_STATUS") != "Variant of unknown origin")
    selection = ["COSMIC_SAMPLE_ID", "HGVSC", "GENE_SYMBOL", "TRANSCRIPT_ACCESSION"]
    rename = {
        "GENE_SYMBOL": "symbol",
        "HGVSC": "hgvs",
        "TRANSCRIPT_ACCESSION": "transcript_id",
    }
    if "MUTATION_DESCRIPTION" in df.columns:
        selection.append("MUTATION_DESCRIPTION")
        rename["MUTATION_DESCRIPTION"] = "consequence"
    df = df.select(selection).rename(rename)
    df = df.join(cosmic_samples, on="COSMIC_SAMPLE_ID").drop("COSMIC_SAMPLE_ID")
    if "consequence" in df.columns:
        df = df.with_columns(pl.col("consequence").str.replace_all(",", ";"))
    return df


def get_variant_class(parser: Parser, val: str) -> str | None:
    try:
        var = parser.parse(val)
        return var.posedit.edit.type
    except HGVSParseError:
        return None


def main():
    parser = Parser()
    formatters: dict[str, Callable[[str], pl.DataFrame]] = {
        "ProteinGym-snps-clinvar": format_proteingym_snps_clinvar,
        "ProteinGym-indels-gnomAD": format_proteingym_indels_gnomad,
        "ProteinGym-indels-clinvar": format_proteingym_indels_clinvar,
        "ClinGen": format_clingen_erepo,
        "CIViC": format_civic,
        "COSMIC_census": format_cosmic,
        "COSMIC_resistance": format_cosmic,
    }
    dfs = [
        read_fn(SOURCES[source]).with_columns(pl.lit(source).alias("source"))
        for source, read_fn in formatters.items()
    ]
    combined = (
        pl.concat(dfs, how="diagonal_relaxed")
        .filter(pl.col("hgvs").is_not_null())
        .unique("hgvs")
    ).with_columns(
        pl.col("hgvs")
        .map_elements(lambda x: get_variant_class(parser, x), return_dtype=pl.String)
        .alias("variant_class"),
        pl.col("clinsig").str.to_lowercase(),
        pl.col("consequence").str.to_lowercase().str.replace_all(" ", "_"),
    )
    combined = combined.join(
        MAPPING.unique("transcript_id").select(["transcript_id", "transcript_len"]),
        on="transcript_id",
        how="left",
    )
    failed = combined.filter(pl.col("variant_class").is_null())
    passed = combined.filter(pl.col("variant_class").is_not_null())
    return passed, failed


if __name__ == "__main__":
    passed, failed = main()
    passed.write_csv(here("data", "processed", "passing_variants.csv"))
    failed.write_csv(here("data", "processed", "failed_variants.csv"))
