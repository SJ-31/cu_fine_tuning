#!/usr/bin/env ipython

import sys
from pathlib import Path

import polars as pl
from pyhere import here

sys.path.append(str(here()))
from search_ontology import SearchOLS, SearchOntology


def lookup_safe(x: str, ont: SearchOntology) -> str | None:
    try:
        return ont.lookup(x)
    except KeyError:
        pass


tmp_file = here("data", "processed", "mondo_disease_mapping_tmp.csv")
failure_file: Path = here("data", "processed", "ols_failed_lookups.txt")

fill_custom = here("data", "processed", "mondo_disease_mapping_review.csv")

if failure_file.exists():
    failures: list = failure_file.read_text().splitlines()
else:
    failures = []
if tmp_file.exists():
    mappings = pl.read_csv(tmp_file).filter(
        (pl.col("disease").is_not_null()) & (pl.col("mondo").is_not_null())
    )
    if fill_custom.exists():
        mappings = pl.concat([mappings, pl.read_csv(fill_custom)])
else:
    mappings = pl.DataFrame({"disease": [], "mondo": []})
df = (
    pl.scan_csv(here("data", "processed", "passing_variants.csv"))
    .select("disease")
    .collect()
    .filter(pl.col("disease").is_not_null())
)
excluded = [
    "not specified",
    "not provided",
    "See cases",
    "non-transition zone (non-transition zone)",
    "no recurrence (no recurrence)",
    ".",
]
diseases = (
    df["disease"]
    .str.split(";")
    .explode()
    .str.replace_all("_", " ")
    .value_counts()
    .sort("count")
    .filter(
        (pl.col("disease").str.len_chars() > 0)
        & (~pl.col("disease").is_in(excluded))
        & (~pl.col("disease").is_in(failures))
        & (~pl.col("disease").is_in(mappings["disease"]))
    )
)

ontology = SearchOntology(
    str(here("data", "mondo.owl")), "MONDO", "http://purl.obolibrary.org/obo/MONDO_"
)
from_ontology = diseases.with_columns(
    pl.Series([lookup_safe(x, ontology) for x in diseases["disease"]]).alias("mondo")
).select(["disease", "mondo"])
to_join = from_ontology.filter(~pl.col("mondo").is_null())


lookups_left = from_ontology.filter(pl.col("mondo").is_null())
ols = SearchOLS(["mondo"], cache=here("data", "ols_cache.db"))
results = []
for d in lookups_left["disease"]:
    try:
        df = ols.search(d)
        results.append(df)
    except Exception as e:
        print(f"Exception when looking up {d}: {e}")
        failures.append(d)

ols_df = pl.concat(results).with_columns(pl.col("description").list.first())
case_same = (
    ols_df.filter(
        (pl.col("label").str.to_lowercase() == pl.col("query").str.to_lowercase())
        | (
            pl.col("query")
            .str.to_lowercase()
            .is_in(pl.col("exact_synonyms").list.eval(pl.element().str.to_lowercase()))
        )
    )
    .select(["query", "obo_id"])
    .rename({"query": "disease", "obo_id": "mondo"})
)
ols_df = ols_df.filter(~pl.col("query").is_in(case_same["disease"]))

# Write all files
#
ols_df.with_columns(pl.col("exact_synonyms").list.join(";")).write_csv(
    here("data", "processed", "mondo_ols_lookups.csv")
)
pl.concat([mappings, to_join, case_same], how="diagonal_relaxed").write_csv(tmp_file)
lookups_left.filter(
    (~pl.col("disease").is_in(failures))
    & (~pl.col("disease").is_in(case_same["disease"]))
).write_csv(fill_custom)
failure_file.write_text("\n".join(failures))
