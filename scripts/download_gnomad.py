#!/usr/bin/env python3

import sys
from pathlib import Path

import polars as pl
from gql import Client, GraphQLRequest, gql
from gql.transport.aiohttp import AIOHTTPTransport
from hgvs.parser import Parser
from pyhere import here

sys.path.append(str(here()))

import main as m

# Retrieve high-quality benign variants from gnomAD

OUT: Path = here("data", "raw", "gnomAD_bg")
OUT.mkdir(exist_ok=True)

MIN_POPMAX: float = 0.05

DB = m.SnpSpace(here("data", "all_snps.db"))
SDB: m.SeqDB = m.SeqDB(here("data", "all_seqs.db"))
SDB.set_aliases(
    pl.read_csv(
        here("data", "mane_transcript_mapping_116_2026-08-12.tsv"), separator="\t"
    ),
    id_col="transcript_mane_select",
    alias_col="ensembl_transcript_id_version",
    namespace="ensembl",
)
MAPPING: dict = SDB.aliases["ensembl"]

TRANSPORT = AIOHTTPTransport(url="https://gnomad.broadinstitute.org/api")
CLIENT = Client(transport=TRANSPORT, fetch_schema_from_transport=True)


# AF is given by AC / AN
#
def create_query(symbol: str) -> GraphQLRequest:
    query = """
    {
    gene(gene_symbol: "_TMP_", reference_genome: GRCh38) {
        variants(dataset: gnomad_r4) {
          pos
          transcript_consequence {
            refseq_id
            refseq_version
            major_consequence
            hgvsc
            hgvs
            gene_symbol
            transcript_id
            transcript_version
            is_mane_select
          }
          exome {
            ac
            an
            populations {
              id
              ac
              an
            }
            filters
            flags
          }
          genome {
            ac
            an
            populations {
              id
              ac
              an
            }
            filters
            flags
          }
        }
      }
    }
    """
    return gql(query.replace("_TMP_", symbol))


# TODO: make the calls more efficient


def get_filter_variants(symbol: str) -> tuple[bool, pl.DataFrame]:
    tmp = {"transcript_id": [], "hgvs": [], "af": []}
    query = create_query(symbol)
    result: dict = CLIENT.execute(query)
    try:
        variants = result["gene"]["variants"]
        for v in variants:
            if not v["exome"]:
                continue
            popmax = max([v["ac"] / v["an"] for v in v["exome"]["populations"]])
            af = v["exome"]["ac"] / v["exome"]["an"]
            if popmax <= MIN_POPMAX:
                continue
            match v:
                case {"transcript_consequence": cons}:
                    refseq = cons.get("refseq_id")
                    version = cons.get("refseq_version")
                    hgvs = cons.get("hgvsc") or cons.get("hgvs")
                    if not hgvs:
                        continue
                    if not refseq:
                        ens = cons.get("transcript_id")
                        version = cons.get("transcript_version")
                        if ens not in SDB:
                            continue
                        transcript = f"{ens}.{version}"
                    else:
                        transcript = f"{refseq}.{version}"
                    tmp["transcript_id"].append(transcript)
                    hgvs = f"{transcript}:{hgvs}"
                    tmp["hgvs"].append(hgvs)
                    tmp["af"].append(af)
                case _:
                    continue
        return True, pl.DataFrame(tmp)
    except KeyError as e:
        print(f"WARNING: KeyError {e} while looking up variants for gene {symbol}")
        return False, pl.DataFrame()


def main():
    generated = pl.scan_csv(here("data", "processed", "generated_passed.csv"))
    from_generated = (
        generated.filter(
            (
                pl.col("clinsig")
                .str.replace_all(" ", "_")
                .is_in(["likely_benign", "benign"])
            )
            & (pl.col("variant_class") == "sub")
        )
        .collect()
        .select(["transcript_id", "hgvs"])
    )
    symbols_to_get: list = (
        generated.filter(pl.col("symbol").is_not_null())
        .select("symbol")
        .collect()["symbol"]
        .unique()
        .to_list()
    )
    failed_file = OUT / "failed.txt"
    if failed_file.exists():
        failed_symbols = (OUT / "failed.txt").read_text().splitlines()
    else:
        failed_symbols = []
    dfs = []
    for symbol in symbols_to_get:
        file = OUT / f"{symbol}.csv"
        if not file.exists():
            passed, df = get_filter_variants(symbol)
            if passed:
                df.write_csv(file)
            else:
                failed_symbols.append(symbol)
                failed_file.write_text("\n".join(failed_symbols))
    all_hgvs = pl.concat([from_generated] + dfs, how="diagonal_relaxed")
    parser = Parser()
    for group, hgvs_df in all_hgvs.group_by("transcript_id"):
        id = group[0]
        hgvs_df["hgvs"]
        DB.add(
            id,
            hgvs=hgvs_df["hgvs"].to_list(),
            parser=parser,
            sdb=SDB,
            af=hgvs_df["af"].to_list(),
        )
