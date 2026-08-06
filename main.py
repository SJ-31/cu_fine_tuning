#!/usr/bin/env python3

from __future__ import annotations

import os
from pathlib import Path
from typing import Literal

import duckdb
import polars as pl
from attrs import Factory, define, field
from Bio.Seq import MutableSeq, Seq
from biocommons.seqrepo import SeqRepo
from hgvs.exceptions import HGVSParseError
from hgvs.location import BaseOffsetPosition
from hgvs.parser import Parser
from hgvs.sequencevariant import SequenceVariant


@define
class Sequence:
    """Class representing a mutable transcript, which indexes
    relative to the start codon, following the HGVS numbering system
    """

    s: MutableSeq
    five_p: int | None
    three_p: int | None
    start_codon: int  # 0-based index of the first base of the start codon
    stop_codon: int  # 0-based index of the last base of the stop codon
    relative_to: Literal["start", "stop"] = "start"
    # Index relative to the start, stop codons. Use for
    # Negative indices relative to start index into the 5' UTR
    # Positive indices relative to stop index into 3' UTR

    def _shift_index(self, i: int | slice) -> int | slice:
        if i == 0:
            raise ValueError("Base index of 0 is not defined")
        elif (
            (isinstance(i, int) and i <= 0) or (isinstance(i, slice) and i.start <= 0)
        ) and self.relative_to == "stop":
            raise ValueError(
                "Negative indexing is not defined when relative to the stop codon"
            )
        offset = (
            self.start_codon if self.relative_to == "start" else self.stop_codon + 1
        )
        if isinstance(i, int):
            if i < 0:
                i += 1
            return offset + i - 1
        start, stop = i.start, i.stop
        if start < 0:
            start += 1
        return slice(offset + start - 1, offset + stop - 1, i.step)

    def __setitem__(self, i: int, base: str):
        self.s[self._shift_index(i)] = base

    def __getitem__(self, i: int | slice) -> Seq | str:
        return self.s[self._shift_index(i)]

    def _adjust_indices(self, pos_change: int, change: int):
        if pos_change < self.start_codon:
            if self.start_codon != 0:
                self.start_codon += change
            self.stop_codon += change
        elif pos_change < self.stop_codon:
            self.stop_codon += change

    def insert(self, i, s):
        def insert_char(i, c):
            shifted = self._shift_index(i)
            self._adjust_indices(shifted, 1)
            self.s.insert(shifted, c)

        if len(s) == 1:
            insert_char(i, s)
        for char in s[::-1]:
            insert_char(i, char)

    def __delitem__(self, i: int | slice):
        def d(idx: int):
            shifted: int = self._shift_index(idx)
            del self.s[shifted]
            self._adjust_indices(shifted, -1)

        if isinstance(i, slice):
            if i.stop <= i.start:
                raise ValueError("Cannot delete with reversed indices")
            for _ in range(i.start, i.stop):
                d(i.start)
        else:
            d(i)

    @classmethod
    def new(
        cls,
        cds: str,
        five_p: str | None = None,
        three_p: str | None = None,
        relative_to: Literal["start", "stop"] = "start",
    ):
        if cds[:3].upper() != "ATG":
            print("WARNING: CDS has no start codon")
        if cds[-3:].upper() not in {"TAG", "TGA", "TAA"}:
            print("WARNING: CDS has no stop codon")
        fp, tp = None, None
        start_codon, stop_codon = 0, len(cds) - 1
        if not five_p and not three_p:
            seq = MutableSeq(cds)
        elif not five_p and three_p:
            seq = MutableSeq(cds + three_p)
            tp = len(cds)
        elif not three_p and five_p:
            seq = MutableSeq(five_p + cds)
            fp = 0
            start_codon = len(five_p)
            stop_codon += len(five_p)
        elif three_p and five_p:
            seq = MutableSeq(five_p + cds + three_p)
            fp, start_codon, tp = 0, len(five_p), len(cds)
            stop_codon += len(five_p)
        return cls(
            s=seq,
            five_p=fp,
            three_p=tp,
            start_codon=start_codon,
            stop_codon=stop_codon,
            relative_to=relative_to,
        )

    def __str__(self) -> str:
        return str(self.s)


@define
class SeqDB:
    file: Path
    aliases: dict[str, dict[str, str]] = field(factory=dict)
    db: duckdb.DuckDBPyConnection = field(
        init=False, default=Factory(lambda x: duckdb.connect(x.file), takes_self=True)
    )

    def set_aliases(
        self, mapping: pl.DataFrame, id_col: str, alias_col: str, namespace: str
    ):
        mapping = mapping.filter(
            (pl.col(id_col).is_not_null()) & (pl.col(alias_col).is_not_null())
        )
        self.aliases[namespace] = dict(zip(mapping[id_col], mapping[alias_col]))

    def fetch(self, id: str, namespace: str | None = None) -> Sequence:
        if namespace is not None:
            id = self.aliases[namespace][id]
        five_p, three_p, cds = self.db.execute(
            """
        SELECT  5p_utr, 3p_utr, cds FROM t WHERE id == ?
        """,
            [id],
        )

    def __attrs_post_init__(self):
        if not self.file.exists():
            self.db.sql("""
            CREATE TABLE t (id PRIMARY_KEY VARCHAR,
            5p_utr VARCHAR,
            3p_utr VARCHAR,
            cds VARCHAR
            )
            """)


class VariantUnsupportedError(Exception):
    pass


def get_pos(p: BaseOffsetPosition) -> int:
    offset = p.offset
    if offset == 0:
        # HGVS syntax is 1-indexed
        return p.base - 1
    raise VariantUnsupportedError(
        "Can only generate intron variants with g. definitions"
    )


def ends(v: SequenceVariant) -> tuple[int, int]:
    return get_pos(v.posedit.pos.start), get_pos(v.posedit.pos.end)


# TODO: check out
# https://github.com/biocommons/hgvs/blob/82500b8f5c9f08a44094096dac9457606735205b/src/hgvs/utils/altseqbuilder.py
# the class is mainly used for HGVSc to HGVSp conversion, so you would
# need to modify it. But good to reference it to check for edge cases

# TODO: cause the RefSeq sequences contain the UTRs, and HGVS numbers
# from the start codon, need to identify them


@define
class VariantGenerator:
    sr: SeqRepo
    parser: Parser = field(factory=Parser)
    seqtype: Literal["aa", "dna"] = "dna"

    def lookup(self, name: str) -> MutableSeq:
        namespace = None
        if name.startswith("ENS"):
            namespace = "ensembl"
        return MutableSeq(self.sr.fetch(name, namespace=namespace))

    def _validate_var(self, v: SequenceVariant) -> None:
        if v.type == "g":
            raise VariantUnsupportedError("Cannot generate from HGVSg")
        if v.type not in {"c", "g", "n"} and self.seqtype == "dna":
            raise VariantUnsupportedError(
                "Can only generate DNA variants from HGVSg or HGVSc strings"
            )

    def extract_gene(self, edited: str, gene: str) -> str:
        """For HGVSg variants, extract the gene subsequence"""
        raise NotImplementedError()

    def gen_repeat(self, id: str, hgvs: str):
        raise NotImplementedError()

    def gen_del(self, id: str, v: SequenceVariant) -> str:
        pass

    def gen_dup(self, id: str, v: SequenceVariant) -> str:
        pass

    def gen_ins(self, id: str, v: SequenceVariant) -> str:
        seq: MutableSeq = self.lookup(id)
        pos = ends(v)
        pass

    def _check_ref(
        self, seq: MutableSeq, pos: tuple[int, int], v: SequenceVariant, type: str
    ):
        v_ref = v.posedit.edit.ref
        if type == "sub":
            ref = seq[pos[0]]
        else:
            ref = seq[pos[0] : pos[1]]
        if ref != v_ref:
            print(
                f"WARNING: ref {v_ref} in variant `{v}` doesn't match ref {ref} in sequence"
            )

    def gen_sub(self, id: str, v: SequenceVariant) -> str:
        seq: MutableSeq = self.lookup(id)
        pos = ends(v)
        self._check_ref(seq, pos, v, "sub")
        seq[pos[0]] = v.posedit.edit.alt
        return str(seq)

    def gen(self, id: str, hgvs: str) -> str:
        try:
            v: SequenceVariant = self.parser.parse(hgvs)
            self._validate_var(v)
            vtype = v.posedit.edit.type
            if vtype == "sub":
                edited = self.gen_sub(id, v)
            elif vtype == "del":
                edited = self.gen_del(id, v)
            elif vtype == "del":
                edited = self.gen_del(id, v)
            elif vtype == "dup":
                edited = self.gen_dup(id, v)
            else:
                raise NotImplementedError(
                    f"Can't yet generate variants of type {vtype}"
                )
            if v.type == "g":
                return self.extract_gene(edited, id)
            return edited
        except HGVSParseError:
            if "[" in hgvs and "]" in hgvs:
                return self.gen_repeat(id, hgvs)
            return ""


# def main():


def parse_args():
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-s", "--seqrepo", default=None, help="seqrepo directory", required=True
    )
    # parser.add_argument("-c", "--", default = , help = "", action = "store")
    parser.add_argument("-o", "--output")
    parser.add_argument(
        "-h",
        "--hgvs_column",
        default="hgvs",
        help="Column in input containing HGVS strings",
        action="store",
    )
    parser.add_argument(
        "-i", "--input_file", default=False, help="Test", action="store_true"
    )
    args = vars(parser.parse_args())
    return args


if __name__ == "__main__":
    args = parse_args()
    sr: SeqRepo = SeqRepo(args["seqrepo"])
