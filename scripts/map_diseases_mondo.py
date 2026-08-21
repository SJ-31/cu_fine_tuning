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


mapped_file = here("data", "processed", "mondo_disease_mapping.csv")
failure_file: Path = here("data", "processed", "ols_failed_lookups.txt")
if failure_file.exists():
    failures: list = failure_file.read_text().splitlines()
else:
    failures = []
if mapped_file.exists():
    mappings = pl.read_csv(mapped_file).filter(
        (pl.col("disease").is_not_null()) & (pl.col("mondo").is_not_null())
    )
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
        & (~pl.col("disease").is_in(mappings["disease"]))
    )
)

ontology = SearchOntology(
    str(here("data", "mondo.owl")), "MONDO", "http://purl.obolibrary.org/obo/MONDO_"
)
from_ontology = diseases.with_columns(
    pl.Series([lookup_safe(x, ontology) for x in diseases["disease"]]).alias("mondo")
).select(["disease", "mondo"])
from_ontology.write_csv(mapped_file)

lookups_left = from_ontology.filter(pl.col("mondo").is_null())
ols = SearchOLS(["mondo"], cache=here("data", "ols_cache.db"))
results = []
for d in lookups_left["disease"]:
    try:
        df = ols.search(d)
        results.append(df)
    except:
        failures.append(d)
pl.concat(results).write_csv(here("data", "processed", "mondo_ols_lookups.csv"))
failure_file.write_text("\n".join(failures))
