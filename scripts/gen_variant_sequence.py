#!/usr/bin/env ipython

import os
from typing import Literal

from attrs import define, field
from Bio.Seq import MutableSeq
from biocommons.seqrepo import SeqRepo
from hgvs.exceptions import HGVSParseError
from hgvs.parser import Parser
from hgvs.sequencevariant import SequenceVariant

sr_dir = "/home/shannc/Bio_SDD/stem_synology/chula_mount/shannc/repos/evo2_fine_tune/seqrepo/2024-12-20/"


def ends(v: SequenceVariant) -> tuple[int, int]:
    return v.posedit.pos.start, v.posedit.pos.end


# TODO: check out
# https://github.com/biocommons/hgvs/blob/82500b8f5c9f08a44094096dac9457606735205b/src/hgvs/utils/altseqbuilder.py
# the class is mainly used for HGVSc to HGVSp conversion, so you would
# need to modify it. But good to reference it to check for edge cases


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
        if v.type not in {"c", "g"} and self.seqtype == "dna":
            raise ValueError(
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
        # BUG: you need to check the indexing for these
        if type == "sub":
            # HGVS are 1-indexed
            ref = seq[pos[0] + 1]
        else:
            ref = seq[pos[0] - 1 : pos[1]]
        if ref != v_ref:
            print(
                f"WARNING: ref {v_ref} in variant `{v}` doesn't match ref {ref} in sequence"
            )

    def gen_sub(self, id: str, v: SequenceVariant) -> str:
        seq: MutableSeq = self.lookup(id)
        pos = ends(v)
        self._check_ref(seq, pos, v, "sub")
        seq[pos[0]] = v.posedit.edit.alt
        return str(id)

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
