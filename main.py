#!/usr/bin/env python3
from __future__ import annotations

try:
    from icecream import ic

    ic.configureOutput("dbg: ", includeContext=True)
except ImportError:  # Graceful fallback if IceCream isn't installed.
    ic = lambda *a: None if not a else (a[0] if len(a) == 1 else a)  # noqa

import argparse
import re
import subprocess as sp
from collections.abc import Callable
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Literal
from zipfile import ZipFile

import duckdb
import polars as pl
import yaml
from attrs import Factory, define, field
from Bio import SeqIO
from Bio.Seq import MutableSeq, Seq
from biocommons.seqrepo import SeqRepo
from hgvs.exceptions import HGVSParseError
from hgvs.location import BaseOffsetPosition, SimplePosition
from hgvs.parser import Parser
from hgvs.sequencevariant import SequenceVariant


@define
class ReferenceSeq:
    """Class representing a mutable reference sequence, which indexes
    following the HGVS numbering system

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
        shifted = self._shift_index(i)
        try:
            return self.s[shifted]
        except IndexError:
            raise VariantUnsupportedError(
                f"Index {shifted} is not defined for sequence of length {len(self.s)}"
            )

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

        if len(s) == 1 and isinstance(s, str):
            insert_char(i, s)
        elif len(s) == 1:
            insert_char(i, str(s))
        else:
            for char in s[::-1]:
                insert_char(i, char)

    def window(
        self,
        pos: int,
        bound: int = 1000,
        upstream: int | None = None,
        downstream: int | None = None,
    ) -> str:
        """Return a string centered on `pos`,
        optionally with upstream and downstream sequence content
        `pos` is interpreted with hgvs numbering

        Parameters
        ----------
        pos : int
            Position to center window on, which will be translated to
            hgvs numbering
        bound : int
            Bound for the length of upstream and downstream sequence to
            include. Both can be set specifically by their parameters

        Returns
        -------
        String of self[pos - upstream:pos] + self[pos] + self[pos + 1: downstream]
        with length upstream + downstream + 1

        """
        pos = self._shift_index(pos)
        up = pos - (bound if upstream is None else upstream)
        up = 0 if up < 0 else up
        down = pos + (bound if downstream is None else downstream)
        s = str(self)
        return s[up:pos] + s[pos] + s[pos + 1 : down + 1]

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
                relative_to="start",
            )
        if not s:
            raise ValueError("CDS is missing")
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
        self,
        mapping: pl.DataFrame | str | Path,
        id_col: str,
        alias_col: str,
        namespace: str,
    ):
        if not isinstance(mapping, pl.DataFrame):
            mapping = pl.read_csv(mapping)
        mapping = mapping.filter(
            (pl.col(id_col).is_not_null()) & (pl.col(alias_col).is_not_null())
        )
        if namespace not in self.aliases:
            self.aliases[namespace] = dict(zip(mapping[id_col], mapping[alias_col]))
        else:
            self.aliases[namespace].update(
                dict(zip(mapping[id_col], mapping[alias_col]))
            )

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
        elif id.startswith("NR_") or id.startswith("NC_"):
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

    def add_transcripts_tabular(
        self,
        id_col: str,
        file: str | Path,
        cds: str = "coding",
        full: str = "cdna",
        fp_utr: str = "5utr",
        tp_utr: str = "3utr",
        missing_val: str | None = "Sequence unavailable",
    ) -> None:
        """
        Add transcripts from a tabular file
        """
        file = Path(file) if isinstance(file, str) else file
        if file.name.endswith(".csv"):
            df = pl.read_csv(file)
        else:
            df = pl.read_csv(file, separator="\t")
        df = df.select([id_col, fp_utr, tp_utr, cds, full])
        if missing_val is not None:
            df = df.with_columns(
                *[
                    pl.col(c).replace(missing_val, None)
                    for c in (cds, full, fp_utr, tp_utr)
                ]
            )
        df = df.filter(~pl.col(id_col).is_in(self.seen))
        self.db.execute("INSERT INTO t SELECT * FROM df")
        self.seen |= set(df[id_col])

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

    def fetch_transcript(self, id: str, namespace: str | None = None) -> ReferenceSeq:
        res = self.fetch(id, namespace)
        if (id.startswith("NR_") or id.startswith("NC_")) and res.get("full"):
            return ReferenceSeq.new(s=res.get("full", ""), is_cds=False)
        return ReferenceSeq.new(
            s=res.get("cds", ""),
            five_p=res.get("5p_utr"),
            three_p=res.get("3p_utr"),
            is_cds=True,
        )

    def fetch(self, id: str, namespace: str | None = None) -> dict:
        if namespace is not None:
            if id in self.aliases[namespace]:
                id = self.aliases[namespace][id]
        res = self.db.execute(
            """
        SELECT fp_utr, tp_utr, cds, full_seq FROM t WHERE id = ?
        """,
            [id],
        ).fetchone()
        if res is None:
            raise KeyError(f"id {id} is not present in database")
        fp, tp, cds, full = res
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


def get_pos(p: BaseOffsetPosition | SimplePosition) -> int:
    if isinstance(p, SimplePosition):
        return p.base
    offset = p.offset
    if offset == 0:
        # WARNING: HGVS syntax is 1-indexed, but this is accounted for
        # by Transcript class
        return p.base
    raise VariantUnsupportedError(
        "Can only generate intron variants with g. definitions"
    )


def ends(v: SequenceVariant) -> tuple[int, int]:
    start, end = get_pos(v.posedit.pos.start), get_pos(v.posedit.pos.end)
    if (end < 0) and (end < start) and (v.type == "c"):
        start, end = end, start
    return start, end


@define
class VariantGenerator:
    db: SeqDB
    sr: SeqRepo
    parser: Parser = field(factory=Parser)
    seqtype: Literal["aa", "dna"] = "dna"
    bounds: tuple[int, int] = (1000, 1000)  # Bounds for genomic variants (type g)

    def lookup(self, name: str, v: SequenceVariant | None = None) -> ReferenceSeq:
        namespace = None
        if name.startswith("ENS"):
            namespace = "ensembl"
        if name not in self.db and not name.startswith("ENS"):
            added, reason = self.db.add_refseq_transcript(name, sr=self.sr)
            if not added:
                raise ValueError(f"Sequence for `{name}` unavailable. Reason: {reason}")
        transcript = self.db.fetch_transcript(name, namespace=namespace)
        if (
            v is not None
            and "datum" in dir(v.posedit.pos.start)
            and transcript.is_cds
            and v.posedit.pos.start.datum.name == "CDS_END"
        ):
            transcript.relative_to = "stop"
        return transcript

    def _validate_var(self, v: SequenceVariant | str) -> None:
        vtype = v.type if isinstance(v, SequenceVariant) else v
        if vtype not in {"c", "g", "n"} and self.seqtype == "dna":
            raise VariantUnsupportedError(
                "Can only generate DNA variants from HGVSg or HGVSc strings"
            )
        if vtype == "n" and isinstance(v, SequenceVariant):
            start, end = ends(v)
            if start <= 0 or end <= 0:
                raise VariantUnsupportedError(
                    "UTR indexing is not defined for non-coding sequences"
                )

    def _convert_string(
        self,
        seq: ReferenceSeq,
        v: SequenceVariant | tuple[str, int],
        extra_upstream: int = 0,
        extra_downstream: int = 0,
    ) -> str:
        if isinstance(v, tuple):
            pos = v[1]
            type = v[0]
        else:
            pos = v.posedit.pos.start.base
            type = v.type
        if type == "g":
            return seq.window(
                pos=pos,
                upstream=self.bounds[0] + extra_upstream,
                downstream=self.bounds[1] + extra_downstream,
            )
        return str(seq)

    def gen_repeat(self, id: str, hgvs: str) -> str:
        seq: ReferenceSeq = self.lookup(id)
        if "*" in hgvs:
            seq.relative_to = "stop"
        if "(" in hgvs:
            raise VariantUnsupportedError("Cannot generate from uncertain sequences")
        if ";" in hgvs:
            raise VariantUnsupportedError("Multi-allelic repeats not supported")
        _, hgvs = hgvs.split(":")
        type = hgvs[0]
        if hgvs[1] != ".":
            raise HGVSParseError("Expected hgvs type")
        hgvs = hgvs[2:]
        check_mixed = len(re.findall("\\[[0-9]+\\]", hgvs)) > 1

        def parse_helper(with_underscore: str, without: str) -> list:
            if "_" in hgvs:
                m = re.findall(with_underscore, hgvs)
            else:
                m = re.findall(without, hgvs)
            if not m:
                raise HGVSParseError("Failed to parse repeat")
            return m

        if check_mixed:
            match = parse_helper("([*-]?[0-9]+)_[*-]?[0-9]+.*", "([*-]?[0-9]+).*")
            repeat_start = match[0]
            tmp = [
                g[0] * int(g[1]) for g in re.findall("([a-zA-Z]+)\\[([0-9]+)\\]", hgvs)
            ]
            to_insert = "".join(tmp)
        else:
            match = parse_helper(
                "([*-]?[0-9]+)_[*-]?[0-9]+([a-zA-Z]+)\\[([0-9]+)\\]",
                "([*-]?[0-9]+)([a-zA-Z]+)\\[([0-9]+)\\]",
            )
            repeat_start, unit, times = match[0]
            to_insert = unit * int(times)
        if "_" in hgvs:
            i1, i2 = re.findall("([*-]?[0-9]+)_([*-]?[0-9]+)", hgvs)[0]
            i1 = int(i1[1:]) if i1.startswith("*") else int(i1)
            i2 = int(i2[1:]) if i2.startswith("*") else int(i2)
            to_insert = to_insert[: i2 - i1]
        if repeat_start.startswith("*"):
            repeat_start = int(repeat_start[1:])
        else:
            repeat_start = int(repeat_start)
        seq.insert(repeat_start + 1, to_insert)
        return self._convert_string(
            seq, (type, repeat_start), extra_downstream=len(to_insert)
        )

    def gen_delins(self, id: str, v: SequenceVariant) -> str:
        seq: ReferenceSeq = self.lookup(id, v)
        pos = ends(v)
        if pos[1] == pos[0]:
            del seq[pos[0]]
        else:
            for _ in range(pos[1] - pos[0] + 1):
                del seq[pos[0]]
        seq.insert(pos[0], v.posedit.edit.alt)
        return self._convert_string(seq, v)

    def gen_del(self, id: str, v: SequenceVariant) -> str:
        """Generate del variant

        Ranges are inclusive, following HGVS syntax
        """
        seq: ReferenceSeq = self.lookup(id, v)
        pos = ends(v)
        try:
            if pos[1] == pos[0]:
                del seq[pos[0]]
            else:
                for _ in range(pos[1] - pos[0] + 1):
                    del seq[pos[0]]
        except IndexError:
            raise ValueError(f"Indexing error in deletion for variant {v}")
        return self._convert_string(seq, v)

    def gen_inv(self, id: str, v: SequenceVariant) -> str:
        seq: ReferenceSeq = self.lookup(id, v)
        pos = ends(v)
        to_insert = str(seq[pos[0] : pos[1] + 1])[::-1]
        for _ in range(pos[1] - pos[0] + 1):
            del seq[pos[0]]
        if not to_insert:
            raise ValueError(f"Range to invert for variant {v} is empty")
        seq.insert(pos[0], to_insert)
        return self._convert_string(seq, v)

    def gen_dup(self, id: str, v: SequenceVariant) -> str:
        """Generate dup variant"""
        seq: ReferenceSeq = self.lookup(id, v)
        pos = ends(v)
        if pos[0] == pos[1]:
            to_dup = seq[pos[0]]
            seq.insert(pos[1], to_dup)
        else:
            to_dup = seq[pos[0] : pos[1] + 1]  # inclusive range following HGVS
            if not to_dup:
                raise ValueError("Duplicated region is empty")
            seq.insert(pos[1] + 1, to_dup)
        return self._convert_string(seq, v)

    def gen_ins(self, id: str, v: SequenceVariant) -> str:
        """
        Generate simple, non self-referential insertion.
        Complex insertions (see https://hgvs-nomenclature.org/stable/recommendations/DNA/insertion/)
        currently unsupported
        """
        seq: ReferenceSeq = self.lookup(id, v)
        pos = ends(v)
        if not v.posedit.pos.uncertain and pos[1] - pos[0] > 1:
            raise ValueError("Insertion range must be adjacent")
        seq.insert(pos[1], v.posedit.edit.alt)
        return self._convert_string(seq, v, extra_downstream=len(v.posedit.edit.alt))

    def _check_ref(
        self, seq: ReferenceSeq, pos: tuple[int, int], v: SequenceVariant, type: str
    ):
        v_ref = v.posedit.edit.ref
        if type == "sub":
            ref = seq[pos[0]]
        else:
            ref = seq[pos[0] : pos[1]]
        if ref != v_ref:
            raise ValueError(
                f"ref {v_ref} in variant `{v}` doesn't match ref {ref} in sequence"
            )

    def gen_sub(self, id: str, v: SequenceVariant) -> str:
        seq: ReferenceSeq = self.lookup(id, v)
        pos = ends(v)
        self._check_ref(seq, pos, v, "sub")
        seq[pos[0]] = v.posedit.edit.alt
        return self._convert_string(seq, v)

    def safe_gen(
        self, id: str, hgvs: str, allow_uncertain: bool = False, as_dict: bool = True
    ) -> tuple[bool, str] | dict:
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
            success, alt = True, result
        except (
            VariantUnsupportedError,
            ValueError,
            KeyError,
            NotImplementedError,
        ) as u:
            success, alt = False, str(u)
        except HGVSParseError as e:
            if ("[" in hgvs) and ("]" in hgvs):
                return False, self.gen_repeat(id, hgvs)
            success, alt = False, f"HGVS library failed to parse: {str(e)}"
        if not as_dict:
            return success, alt
        return {"gen_success": success, "alt_seq": alt}

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
            elif vtype == "inv":
                edited = self.gen_inv(id, v)
            elif vtype == "delins":
                edited = self.gen_delins(id, v)
            else:
                raise NotImplementedError(f"Can't generate variants of type {vtype}")
            return edited
        except HGVSParseError as e:
            if isinstance(hgvs, str) and ("[" in hgvs) and ("]" in hgvs):
                return self.gen_repeat(id, hgvs)
            raise e


def spec_helper(file: str, fn: Callable, spec_name: str, file_key: str = "file"):
    with open(file, "r") as f:
        lst = yaml.safe_load(f)
        if not isinstance(lst, list):
            print(f"WARNING: {spec_name} specification is not a list. Ignoring...")
            return
        for i, group in enumerate(lst):
            if not isinstance(group, dict):
                print(f"WARNING: item {i} of spec is not a mapping")
                continue
            file = group[file_key]
            if not Path(file).exists():
                raise ValueError(f"file in item {i} of mapping does not exist")
            fn(**group)


# * CLI entry


def gen_batch(
    generator: VariantGenerator, args: dict, batch_df: pl.DataFrame, write_prefix: int
) -> tuple[pl.DataFrame, pl.DataFrame]:
    wd: Path = Path(args["workdir"])
    result = batch_df.with_columns(
        pl.Series(
            [
                generator.safe_gen(id, hgvs, as_dict=True)
                for id, hgvs in zip(
                    batch_df[args["id_column"]], batch_df[args["hgvs_column"]]
                )
            ]
        ).alias("fields")
    ).unnest("fields")
    failed = result.filter(~pl.col("gen_success")).rename({"alt_seq": "fail_reason"})
    passed = result.filter(pl.col("gen_success"))
    f_write, p_write = (
        wd / f"{write_prefix}_failed.csv",
        wd / f"{write_prefix}_passed.csv",
    )
    if f_write.exists():
        raise FileExistsError(f"File for failed entries {f_write} should not exist")
    if p_write.exists():
        raise FileExistsError(f"File for passing entries {p_write} should not exist")
    passed.write_csv(p_write)
    failed.write_csv(f_write)
    return passed, failed


def main(args: dict):
    wd: Path = Path(args["workdir"])
    if not wd.exists():
        wd.mkdir()
    seqdb = SeqDB(file=args["database"])
    if args["alias_spec"] is not None:
        spec_helper(args["alias_spec"], seqdb.set_aliases, "alias", "mapping")
    if args["load_sequences"] is not None:
        spec_helper(
            args["load_sequences"], seqdb.add_transcripts_tabular, "tabular sequence"
        )
    parser = Parser()
    sr = SeqRepo(args["seqrepo"])
    df: pl.DataFrame = pl.read_csv(args["input"])
    generator = VariantGenerator(db=seqdb, parser=parser, seqtype="dna", sr=sr)
    previous = list(wd.glob("*csv"))
    if not previous:
        start_index = 0
    else:
        start_index = max([int(file.stem.split("_")[0]) for file in previous]) + 1
    hgvs_col: str = args["hgvs_column"]
    if previous:
        attempted_hgvs = pl.concat([pl.read_csv(f).select(hgvs_col) for f in previous])[
            hgvs_col
        ].to_list()
        df = df.filter(~pl.col(hgvs_col).is_in(attempted_hgvs))
        failed_tmp = [pl.scan_csv(f) for f in wd.glob("*_failed.csv")]
        passed_tmp = [pl.scan_csv(f) for f in wd.glob("*_passed.csv")]
    else:
        failed_tmp, passed_tmp = [], []
    for batch in df.iter_slices(args["save_interval"]):
        passed, failed = gen_batch(
            batch_df=batch, args=args, generator=generator, write_prefix=start_index
        )
        failed_tmp.append(failed.lazy())
        passed_tmp.append(passed.lazy())
        start_index += 1
    return pl.concat(passed_tmp, how="diagonal_relaxed").collect(), pl.concat(
        failed_tmp, how="diagonal_relaxed"
    ).collect()


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-s", "--seqrepo", default=None, help="seqrepo directory", required=True
    )
    parser.add_argument(
        "-d",
        "--database",
        default=None,
        help="Sequence database file",
        action="store",
        required=True,
    )
    parser.add_argument("-p", "--output_passed")
    parser.add_argument("-f", "--output_failed")
    parser.add_argument(
        "-g",
        "--hgvs_column",
        default="hgvs",
        help="Column in input containing HGVS strings",
        action="store",
    )
    parser.add_argument(
        "-a",
        "--alias_spec",
        default=None,
        help="""
        YAML file containing sequence aliases.
        The format is a list of mappings with four keys:
        - file: the path to the mapping file
        - id_col: column in the mapping file with identifiers 
        - alias_col: column in the mapping file with identifier aliases
        - namespace: namespace key in sequence database to use e.g. `ensembl`
        """,
        action="store",
    )
    parser.add_argument(
        "-l",
        "--load_sequences",
        default=None,
        help="""YAML file specifying tabular files to initially load into SeqDB""",
        action="store",
    )
    parser.add_argument(
        "-w", "--workdir", help="Working directory to cache results", action="store"
    )
    parser.add_argument(
        "-v",
        "--save_interval",
        default=3000,
        help="Number of sequences to generate before saving a batch to the working directory",
        action="store",
        type=int,
    )
    parser.add_argument(
        "-t",
        "--id_column",
        default="transcript_id",
        help="Column containing transcript identifiers",
        action="store",
    )
    parser.add_argument("-i", "--input", default=False, help="Test", action="store")
    args = vars(parser.parse_args())
    return args


if __name__ == "__main__":
    args = parse_args()
    passed, failed = main(args)
    passed.write_csv(args["output_passed"])
    failed.write_csv(args["output_failed"])
