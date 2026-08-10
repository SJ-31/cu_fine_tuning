#!/usr/bin/env python3

from __future__ import annotations

import os
import subprocess as sp
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Literal
from zipfile import ZipFile

import duckdb
import polars as pl
from attrs import Factory, define, field
from Bio import SeqIO
from Bio.Seq import MutableSeq, Seq
from biocommons.seqrepo import SeqRepo
from hgvs.exceptions import HGVSParseError
from hgvs.location import BaseOffsetPosition
from hgvs.parser import Parser
from hgvs.sequencevariant import SequenceVariant


@define
class Transcript:
    """Class representing a mutable transcript, which indexes
    relative to the start codon, following the HGVS numbering system

    Parameters
    ----------
    start : int
        If `is_cds`, 0-based index of the first base of the start codon
    Otherwise, 0 for the first base
    end : int
        If `is_cds`, the 0-based index of the last base of the stop codon
    Otherwise, the index at the last base of the sequence
    relative_to : Literal["start", "stop", None]
    How to interpret indices, following HGVS numbering conventions
    - Negative indices relative to start index into the 5' UTR
    - Positive indices relative to stop index into 3' UTR
    - Setting to None indexes directly into `self.s`
    is_cds : bool
        Whether `s` is a coding sequence
    """

    s: MutableSeq
    five_p: int | None
    three_p: int | None
    start: int
    end: int
    is_cds: bool
    relative_to: Literal["start", "stop", None] = "start"

    def _shift_index(self, i: int | slice) -> int | slice:
        if self.relative_to == "none":
            return i
        if i == 0:
            raise ValueError("Base index of 0 is not defined")
        elif (
            (isinstance(i, int) and i <= 0) or (isinstance(i, slice) and i.start <= 0)
        ) and self.relative_to == "stop":
            raise ValueError(
                "Negative indexing is not defined when relative to the stop codon"
            )
        offset = self.start if self.relative_to == "start" else self.end + 1
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
        if pos_change < self.start:
            if self.start != 0:
                self.start += change
            self.end += change
        elif pos_change < self.end:
            self.end += change

    def insert(self, i, s):
        def insert_char(i, c):
            shifted = self._shift_index(i)
            if self.is_cds:
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
            if self.is_cds:
                self._adjust_indices(shifted, -1)
            else:
                self.end = len(self.s) - 1

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
        s: str,
        five_p: str | None = None,
        three_p: str | None = None,
        relative_to: Literal["start", "stop", None] = "start",
        is_cds: bool = True,
    ):
        if not is_cds:
            return cls(
                s=MutableSeq(s),
                start=0,
                end=len(s) - 1,
                five_p=None,
                three_p=None,
                is_cds=False,
                relative_to=relative_to,
            )
        if s[:3].upper() != "ATG":
            print("WARNING: CDS has no start codon")
        if s[-3:].upper() not in {"TAG", "TGA", "TAA"}:
            print("WARNING: CDS has no stop codon")
        fp, tp = None, None
        start_codon, stop_codon = 0, len(s) - 1
        if not five_p and not three_p:
            seq = MutableSeq(s)
        elif not five_p and three_p:
            seq = MutableSeq(s + three_p)
            tp = len(s)
        elif not three_p and five_p:
            seq = MutableSeq(five_p + s)
            fp = 0
            start_codon = len(five_p)
            stop_codon += len(five_p)
        elif three_p and five_p:
            seq = MutableSeq(five_p + s + three_p)
            fp, start_codon, tp = 0, len(five_p), len(s)
            stop_codon += len(five_p)
        return cls(
            s=seq,
            five_p=fp,
            three_p=tp,
            start=start_codon,
            end=stop_codon,
            relative_to=relative_to,
            is_cds=True,
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
    seen: set = field(init=False, factory=set)

    def set_aliases(
        self, mapping: pl.DataFrame, id_col: str, alias_col: str, namespace: str
    ):
        mapping = mapping.filter(
            (pl.col(id_col).is_not_null()) & (pl.col(alias_col).is_not_null())
        )
        self.aliases[namespace] = dict(zip(mapping[id_col], mapping[alias_col]))

    @staticmethod
    def download(id: str, take_first: bool = True) -> dict[str, str]:
        """Download RefSeq sequence given by `id` with NCBI datasets.

        If `id` can't be found in the fasta files provided by the data package and `take_first`,
        retrieve the first sequence in the fasta file
        """
        tmp = {}
        with TemporaryDirectory() as d:
            dir = Path(d)
            zip_file = dir / "dataset.zip"
            command = [
                "datasets",
                "download",
                "gene",
                "accession",
                id,
                "--include",
                "gene,cds,3p-utr,5p-utr",
                "--no-progressbar",
                "--filename",
                str(zip_file),
            ]
            proc = sp.run(" ".join(command), shell=True, capture_output=True)
            proc.check_returncode()
            with ZipFile(zip_file, "r") as z:
                contents = z.namelist()
                for seqtype in ("gene", "cds", "3p_utr", "5p_utr"):
                    extract_to = dir / seqtype
                    path = f"ncbi_dataset/data/{seqtype}.fna"
                    if path in contents:
                        z.extract(member=path, path=extract_to)
                        file = extract_to / path
                        seqs = [seq for seq in SeqIO.parse(file, "fasta")]
                        tmp[seqtype] = seqs or []
        result = {}
        for k, seqlist in tmp.items():
            for seq in seqlist:
                if seq.id.startswith(id) and ":" in seq.id:
                    result[k] = str(seq.seq)
                    break
            if not result.get(k) and take_first and seqlist:
                result[k] = str(seqlist[0].seq)
        return result

    def add_refseq_transcript(
        self,
        id: str,
        sr: SeqRepo,
        mapping_key: str | None = "ensembl",
        lookup: dict | None = None,
    ) -> tuple[bool, str | None]:
        """
        Helper function to add RefSeq transcript into db

        Returns
        -------
        tuple[bool, str]
        The first element indicates whether `id` was added into `db`
        If this is False, the second element is the reason why
        """
        if id not in self.seen:
            self.seen.add(id)
        else:
            return False, "already in db"
        if not lookup and mapping_key:
            lookup = self.aliases.get(mapping_key)
        if id.startswith("ENST") and lookup:
            if id not in lookup:
                return False, "could not map from Ensembl to RefSeq"
            else:
                id = lookup[id]
        elif id.startswith("NR_"):
            fp, tp, cds = "", "", ""
            try:
                full = sr.fetch(id)
            except KeyError:
                return False, "not found in SeqRepo"
        elif id.startswith("NM_"):
            try:
                downloaded = self.download(id)
                fp = downloaded.get("5p_utr", "")
                tp = downloaded.get("3p_utr", "")
                cds = downloaded.get("cds")
                full = ""
                if not cds:
                    return False, "no CDS could be downloaded with datasets"
            except sp.CalledProcessError:
                return False, "datasets raised CalledProcessError"
        else:
            return False, "unsupported prefix"
        to_insert = [id, fp, tp, cds, full]
        self.db.execute("INSERT INTO t VALUES (?, ?, ?, ?, ?)", to_insert)
        return True, None

    def __contains__(self, val: str) -> bool:
        return val in self.seen

    def add_refseq_transcripts(
        self, ids: list[str], sr: SeqRepo, mapping_key: str | None = "ensembl"
    ) -> dict[str, str]:
        """
        Add RefSeq transcripts in `ids` to db.

        Parameters
        ----------
        ids : list[str]
            List of RefSeq or Ensembl transcript ids (MANE).
            The latter will be mapped
            using `self.aliases` by `mapping_key` if provided
            ids must be prefixed with `NM_`, `NR_`, or `ENST`

        BUG: because datasets doesn't provide the NR_ accessions in its
            data packages for non-coding transcripts, must use SeqRepo instead

        Returns
        -------
        Dict of failed ids > failure reason
        """
        if mapping_key:
            lookup: dict | None = self.aliases.get(mapping_key)
            if not lookup:
                print(
                    f"WARNING: {mapping_key} not in `aliases`, cannot map Ensembl transcripts"
                )
        else:
            print("WARNING: no mapping key provided, cannot map Ensembl transcripts")
            lookup = {}
        failed = {}
        ids = list(set(ids))
        for id in ids:  # Save each id individually in case of network errors
            if id.startswith("ENST") and lookup:
                if id not in lookup:
                    failed[id] = "could not map from Ensembl to RefSeq"
                    self.seen.add(id)
                    continue
                else:
                    id = lookup[id]
            if id in self.seen:
                continue
            added, comment = self.add_refseq_transcript(id, sr=sr, lookup=lookup)
            if not added:
                failed[id] = comment
        return failed

    def fetch_transcript(self, id: str, namespace: str | None = None) -> Transcript:
        res = self.fetch(id, namespace)
        if id.startswith("NR_") and res.get("full"):
            return Transcript.new(s=res.get("full", ""), is_cds=False)
        return Transcript.new(
            s=res.get("cds", ""),
            five_p=res.get("5p_utr"),
            three_p=res.get("3p_utr"),
            is_cds=True,
        )

    def fetch(self, id: str, namespace: str | None = None) -> dict:
        if namespace is not None:
            id = self.aliases[namespace][id]
        fp, tp, cds, full = self.db.execute(
            """
        SELECT  fp_utr, tp_utr, cds, full_seq FROM t WHERE id == ?
        """,
            [id],
        ).fetchone()
        return {"5p_utr": fp, "3p_utr": tp, "cds": cds, "full": full}

    def __attrs_post_init__(self):
        self.db.sql("""
        CREATE TABLE IF NOT EXISTS t (
            id VARCHAR PRIMARY KEY,
            fp_utr VARCHAR,
            tp_utr VARCHAR,
            cds VARCHAR,
            full_seq VARCHAR
        );
            """)

        self.seen |= {p[0] for p in self.db.execute("SELECT id FROM t").fetchall()}


class VariantUnsupportedError(Exception):
    pass


def get_pos(p: BaseOffsetPosition) -> int:
    offset = p.offset
    if offset == 0:
        # WARNING: HGVS syntax is 1-indexed, but this is accounted for
        # by Transcript class
        return p.base
    raise VariantUnsupportedError(
        "Can only generate intron variants with g. definitions"
    )


def ends(v: SequenceVariant) -> tuple[int, int]:
    return get_pos(v.posedit.pos.start), get_pos(v.posedit.pos.end)


@define
class VariantGenerator:
    db: SeqDB
    sr: SeqRepo
    parser: Parser = field(factory=Parser)
    seqtype: Literal["aa", "dna"] = "dna"

    def lookup(self, name: str, v: SequenceVariant) -> Transcript:
        namespace = None
        if name.startswith("ENS"):
            namespace = "ensembl"
        if name not in self.db:
            added, reason = self.db.add_refseq_transcript(name, sr=self.sr)
            if not added:
                raise ValueError(f"Sequence for `{name}` unavailable. Reason: {reason}")
        transcript = self.db.fetch_transcript(name, namespace=namespace)
        if (
            "datum" in dir(v.posedit.pos.start)
            and transcript.is_cds
            and v.posedit.pos.start.datum.name == "CDS_END"
        ):
            transcript.relative_to = "stop"
        return transcript

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
        seq: Transcript = self.lookup(id)
        pos = ends(v)
        pass

    def _check_ref(
        self, seq: Transcript, pos: tuple[int, int], v: SequenceVariant, type: str
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
        seq: Transcript = self.lookup(id, v)
        pos = ends(v)
        self._check_ref(seq, pos, v, "sub")
        seq[pos[0]] = v.posedit.edit.alt
        return str(seq)

    def gen_safe(
        self, id: str, hgvs: str, allow_uncertain: bool = False
    ) -> tuple[bool, str]:
        """
        Wrapper around variant generation to ignore unsupported or failed
        variants

        Returns
        -------
        tuple of
            - bool : whether the variant was successfully generated
            - str : variant sequence if first element was true,
            otherwise reason indicating why generation failed
        """
        try:
            v: SequenceVariant = self.parser.parse(hgvs)
            if v.posedit.pos.uncertain and not allow_uncertain:
                return False, "Variant position uncertain"
            self._validate_var(v)
            result = self.gen(id, hgvs=v)
            return True, result
        except (VariantUnsupportedError, ValueError) as u:
            return False, str(u)
        except HGVSParseError as e:
            return False, f"HGVS library failed to parse: {str(e)}"

    def gen(self, id: str, hgvs: str | SequenceVariant) -> str:
        """
        Generate variant from HGVS string

        Parameters
        ----------
        id : str
            Sequence identifier for variant
        """
        try:
            if isinstance(hgvs, str):
                v: SequenceVariant = self.parser.parse(hgvs)
            else:
                v = hgvs
            vtype = v.posedit.edit.type
            if vtype == "sub":
                edited = self.gen_sub(id, v)
            elif vtype == "del":
                edited = self.gen_del(id, v)
            elif vtype == "ins":
                edited = self.gen_ins(id, v)
            elif vtype == "dup":
                edited = self.gen_dup(id, v)
            else:
                raise NotImplementedError(f"Can't generate variants of type {vtype}")
            if v.type == "g":
                return self.extract_gene(edited, id)
            return edited
        except HGVSParseError as e:
            if "[" in hgvs and "]" in hgvs:
                return self.gen_repeat(id, hgvs)
            raise e


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
