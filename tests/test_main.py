#!/usr/bin/env python3

import sys

import gffutils
import polars as pl
import pyfaidx
import pytest
from biocommons.seqrepo import SeqRepo
from hgvs.parser import Parser
from pyhere import here

sys.path.append(str(here()))
import main as m

DATA = here("data")

SR: SeqRepo = SeqRepo(DATA / "seqrepo/2024-12-20/")
DB = DATA / "GCF_000001405.40.db"

# TODO: record creation as
# db = gffutils.create_db(
#     gff, data / "GCF_000001405.40.db", merge_strategy="create_unique", force=True
# )

HP = Parser()


"""
RefSeq         gene_name    strand    seq
NR_023317.1     RNU7-1     1          CAGTGTTACAGCTCTTTTAGAATTTGTCTAGTAGGCTTTCTGGCTTTTTACCGGAAAGCCCCT
NR_104088.1     RNU6-8     -1         GTGCTCGCTTCGGCAGCACATATACTAAAATTGGAACGATACAGAGAAGATTAGCATGGCCCCTGCGCAAGGATGACACGCAAATTCGTGAAGCGTTCCATATTTTG

"""

snps = pl.read_csv(here("tests", "snps.csv"))

downloads = pl.scan_csv(here("tests", "download_test.csv"))


@pytest.mark.parametrize("id", list(downloads.collect()["id"]))
def test_download(id):
    print(id)
    lookup = downloads.collect().rows_by_key("id", unique=True, named=True)[id]
    fp, tp, cds = lookup["fp"], lookup["tp"], lookup["cds"]
    result = m.SeqDB.download(id)
    assert fp == result["5p_utr"]
    assert tp == result["3p_utr"]
    assert cds == result["cds"]


def test_seq_delete():
    seq = m.Sequence.new("ATGAGACTAGACAGTGA", "fiveprime", "threeprime", "start")
    old_stop = seq.stop_codon
    old_start = seq.start_codon
    assert seq[1] == "A"
    del seq[2]
    assert old_start == seq.start_codon
    assert seq[2] == "G"
    assert seq.stop_codon == old_stop - 1
    del seq[2:5]
    assert seq.stop_codon == old_stop - 1 - 3
    assert str(seq[2:5]) == "ACT"


def test_seq():
    seq = m.Sequence.new("ATGAGACTAGACAGTGA", "fiveprime", "threeprime", "start")
    assert seq[-1] == "e"
    assert seq[1] == "A"
    assert seq[3] == "G"
    seq[1] = "G"
    assert seq[1] == "G"
    assert seq[1:3] == "GT"
    seq.relative_to = "stop"
    assert seq[1] == "t"
    assert seq[2] == "h"
    seq.insert(3, "inserted")
    assert seq[3:11] == "inserted"
    del seq[2:5]


# TODO: [2026-08-06 Thu] will need to download the sequences with
# NCBI's datasets, cause you can't risk finding the ORF manually.


# @pytest.mark.parametrize("gene,hgvs,alt,supported", list(snps.iter_rows()))
# def test_sub(gene, hgvs: str, alt: str, supported: bool):
#     G = m.VariantGenerator(sr=SR, parser=HP, seqtype="dna")
#     if supported:
#         generated = G.gen(gene, hgvs)
#         assert alt == generated
#     else:
#         with pytest.raises(m.VariantUnsupportedError):
#             G.gen(gene, hgvs)
