#!/usr/bin/env python3

import sys

import numpy as np
import polars as pl
import pytest
from Bio import Align
from biocommons.seqrepo import SeqRepo
from hgvs.parser import Parser
from pyhere import here

sys.path.append(str(here()))
import main as m

DATA = here("data")

SR: SeqRepo = SeqRepo(DATA / "seqrepo/2024-12-20/")
DB = DATA / "GCF_000001405.40.db"
HP = Parser()


"""
RefSeq         gene_name    strand
NR_023317.1     RNU7-1     1
NR_104088.1     RNU6-8     -1

"""
EXTRA = {
    "NR_023317.1": "CAGTGTTACAGCTCTTTTAGAATTTGTCTAGTAGGCTTTCTGGCTTTTTACCGGAAAGCCCCT",
    "NR_104088.1": "GTGCTCGCTTCGGCAGCACATATACTAAAATTGGAACGATACAGAGAAGATTAGCATGGCCCCTGCGCAAGGATGACACGCAAATTCGTGAAGCGTTCCATATTTTG",
}
DOWNLOADS = pl.scan_csv(here("tests", "data", "download_test.csv"))
IDS = [
    "NR_023317.1",
    "NR_104088.1",
    "NM_052951.3",
    "NM_001172655.1",
    "fail",
    "NR_fail",
    "NM_fail",
]


def edit_distance(a: str, b: str) -> int:
    """
    Calculate the edit distance between a, b using the Wagner-Fischer algorithm
    https://en.wikipedia.org/wiki/Wagner%E2%80%93Fischer_algorithm

    Returns
    -------
    Levenshtein edit distance from string a to b
    """
    m, n = len(a), len(b)
    d = np.zeros((m + 1, n + 1))
    for i in range(1, m + 1):
        d[i, 0] = i + 1
    for j in range(1, n + 1):
        d[0, j] = j + 1
    for j in range(1, n + 1):
        for i in range(1, m + 1):
            if a[i - 1] == b[j - 1]:
                sub_cost = 0
            else:
                sub_cost = 1
            d[i, j] = min(
                d[i - 1, j] + 1,  # deletion
                d[i, j - 1] + 1,  # insertion
                d[i - 1, j - 1] + sub_cost,  # substitution
            )
    return d[m - 1, n - 1]


@pytest.fixture
def default_db(tmp_path):
    db_file = tmp_path / "seqs.db"
    db = m.SeqDB(file=db_file)
    db.set_aliases(
        pl.read_csv(here("data", "mart_2026-08-03_filtered.csv")),
        id_col="RefSeq match transcript (MANE Select)",
        alias_col="Transcript stable ID version",
        namespace="ensembl",
    )
    return db


def test_add_tabular(default_db):
    db: m.SeqDB = default_db
    db.add_refseq_transcripts(IDS, sr=SR)
    db.add_transcripts_tabular(
        id_col="ensembl_transcript_id_version",
        file=here("tests", "data", "test_tabular_seqs.tsv"),
    )
    assert "ENST00000002829.8" in db.seen
    assert "NR_104088.1" in db.seen
    assert db.fetch("ENST00000002829.8")["full"] is not None
    assert db.fetch("NR_104088.1")["full"] is not None
    assert db.fetch("ENST00000216019.11")["5p_utr"] is None
    with pytest.raises(KeyError):
        db.fetch("FOOO")


def test_seqdb(default_db):
    db: m.SeqDB = default_db
    lookup = DOWNLOADS.collect().rows_by_key("id", unique=True, named=True)
    failed = db.add_refseq_transcripts(IDS, sr=SR)
    assert "fail" in failed
    assert "NR_fail" in failed
    assert "NM_fail" in failed
    assert db.fetch(IDS[0])["full"] == EXTRA[IDS[0]]
    assert db.fetch(IDS[1])["full"] == EXTRA[IDS[1]]
    for key in ("5p_utr", "3p_utr", "cds"):
        assert db.fetch(IDS[2])[key] == lookup[IDS[2]][key]
        assert db.fetch(IDS[3])[key] == lookup[IDS[3]][key]
    add, _ = db.add_refseq_transcript(IDS[0], sr=SR, mapping_key="ensembl")
    assert not add
    assert db.seen == set(IDS)


@pytest.mark.parametrize("id", list(DOWNLOADS.collect()["id"]))
def test_download(id):
    lookup = DOWNLOADS.collect().rows_by_key("id", unique=True, named=True)[id]
    fp, tp, cds = lookup["5p_utr"], lookup["3p_utr"], lookup["cds"]
    result = m.SeqDB.download(id)
    assert fp == result["5p_utr"]
    assert tp == result["3p_utr"]
    assert cds == result["cds"]


def test_seq_delete():
    seq = m.Transcript.new("ATGAGACTAGACAGTGA", "fiveprime", "threeprime", "start")
    old_stop = seq.end
    old_start = seq.start
    assert seq[1] == "A"
    del seq[2]
    assert old_start == seq.start
    assert seq[2] == "G"
    assert seq.end == old_stop - 1
    del seq[2:5]
    assert seq.end == old_stop - 1 - 3
    assert str(seq[2:5]) == "ACT"


def test_seq():
    seq = m.Transcript.new("ATGAGACTAGACAGTGA", "fiveprime", "threeprime", "start")
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


def test_seq_insert():
    seq = m.Transcript.new("ATG", "fp", "tp", "start")
    seq.insert(2, "FOO")
    assert str(seq) == "fpAFOOTGtp"
    s2 = m.Transcript.new("ATG", is_cds=False)
    s2.insert(2, "FOO")
    assert str(s2) == "AFOOTG"


# [2026-08-10 Mon] TODO: should add c. type tests for delins, dup, too
# repeats, inv


@pytest.mark.parametrize(
    "vtype", ["snps", "ins", "del", "delins", "dup", "repeat", "inv"]
)
def test_from_csv(vtype, default_db, subtests):
    file = here("tests", "data", f"{vtype}.csv")
    db = default_db
    G = m.VariantGenerator(db=db, parser=HP, seqtype="dna", sr=SR)
    aligner = Align.PairwiseAligner()
    df = pl.read_csv(
        file,
        schema={
            "gene": pl.String,
            "hgvs": pl.String,
            "alt": pl.String,
            "supported": pl.Boolean,
            "pos_check": pl.String,
            "char": pl.String,
        },
    )
    for gene, hgvs, alt, supported, pos_check, char in df.iter_rows():
        with subtests.test(msg=f"Testing {hgvs}"):
            if supported:
                generated = G.gen(gene, hgvs)
                if alt:
                    print(aligner.align(alt, generated)[0])
                    assert alt == generated
                else:
                    seq = SR.fetch(gene)
                    score = edit_distance(seq, generated)
                    try:
                        alignments = aligner.align(seq, generated)
                        a1 = alignments[0]
                        print(f"Alignment, score {a1.score} (max {len(seq)})\n{a1}")
                    except OverflowError:
                        pass
                    if vtype == "snps":
                        expected_score = 1
                    elif vtype in "ins":
                        expected_score = len(generated) - len(seq)
                    elif vtype == "del":
                        expected_score = len(seq) - len(generated)
                    assert score == expected_score
                if pos_check and char:
                    if "-" not in pos_check:
                        assert generated[int(pos_check)] == char
                    else:
                        start, end = pos_check.split("-")
                        assert generated[int(start) : int(end)] == char
            else:
                with pytest.raises(m.VariantUnsupportedError):
                    G.gen(gene, hgvs)
