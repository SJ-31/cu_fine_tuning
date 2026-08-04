#!/usr/bin/env ipython


import os
from pathlib import Path

import polars as pl
from pyhere import here

raw: Path = here("data", "raw", "ProteinGym")
remote: Path = here("data")


# TODO: [2026-08-03 Mon] this script should generate the
# HGVSc codes for all the ProteinGym samples,
# Also generate a metadata file with
# extract ClinVar
# metadata including clinical significance,
# the source of the variant

clinvar: pl.DataFrame = pl.scan_csv(here("data", "raw", "clinvar_hgvsg.csv"))

indels = pl.scan_csv(raw / "indels_preprocessed_clinvar.csv").select("VariationID")
indels_control = pl.scan_csv(raw / "indels_preprocessed_gnomad.csv")
substitutions = pl.scan_csv(raw / "substitutions_preprocessed.csv").select(
    ["clinsig", "HGVSc"]
)
