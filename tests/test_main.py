#!/usr/bin/env python3

import sys

import polars as pl
import pytest
from biocommons.seqrepo import SeqRepo
from hgvs.parser import Parser
from pyhere import here

sys.path.append(str(here()))
import main as m

DATA = here("data")

SR: SeqRepo = SeqRepo(DATA / "seqrepo/2024-12-20/")

sr_dir = "/home/shannc/Bio_SDD/stem_synology/chula_mount/shannc/repos/evo2_fine_tune/seqrepo/2024-12-20/"

HP = Parser()


"""
RefSeq         gene_name    strand    seq
NR_023317.1     RNU7-1     1          CAGTGTTACAGCTCTTTTAGAATTTGTCTAGTAGGCTTTCTGGCTTTTTACCGGAAAGCCCCT
NR_104088.1     RNU6-8     -1         GTGCTCGCTTCGGCAGCACATATACTAAAATTGGAACGATACAGAGAAGATTAGCATGGCCCCTGCGCAAGGATGACACGCAAATTCGTGAAGCGTTCCATATTTTG

"""

snps = pl.read_csv(here("tests", "snps.csv"))


@pytest.mark.parametrize("gene,hgvs,alt", list(snps.iter_rows()))
def test_sub(gene, hgvs: str, alt: str):
    G = m.VariantGenerator(sr=SR, parser=HP, seqtype="dna")
    generated = G.gen(gene, hgvs)
    assert alt == generated
