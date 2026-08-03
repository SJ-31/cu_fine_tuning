#!/usr/bin/env ipython

from pathlib import Path

import polars as pl
from pyhere import here

raw: Path = here("data", "raw", "ProteinGym")

# TODO: [2026-08-03 Mon] this script should generate the
# HGVSc codes for all the ProteinGym samples,
# Also generate a metadata file with
# extract ClinVar
# metadata including clinical significance,
# the source of the variant

indels = pl.scan_csv(raw / "indels_preprocessed_clinvar.csv")
indels_control = pl.scan_csv(raw / "indels_preprocessed_gnomad.csv")["VariationID"]
substitutions = pl.read_csv()
