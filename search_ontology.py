#!/usr/bin/env ipython


from pathlib import Path

import duckdb
import polars as pl
import pronto
import requests
from attrs import Factory, define, field

# obo = pronto.Ontology(file)


@define
class SearchOLS:
    """
    Simple class interfacing with Ontology Lookup Service.
    Caches results in a duckdb database for future lookups.

    The central database table has columns for "queryFields" + "obo_id",
    a separate table stores previous lookups

    Parameters
    ----------
    fields : dict[str, str]
        Additional fields to store from the look up. `obo_id` and `label`
        are always included. A dictionary mapping the name of the key
        (see [1] for options) to the SQL storage type (default VARCHAR)
    query_fields : tuple
        Fields to search for in the query
        WARNING: searching in the `description` field takes a really
        long time

    References
    ----------
    [1] https://www.ebi.ac.uk/ols4/api-docs
    """

    ontologies: list
    cache: Path = field(converter=Path)
    fields: dict[str, str] = field(
        factory=lambda: {"description": "VARCHAR[]", "exact_synonyms": "VARCHAR[]"}
    )
    query_fields: tuple[str, ...] = ("label", "synonym")
    endpoint: str = "https://www.ebi.ac.uk/ols4"
    timeout: int | None = None
    db: duckdb.DuckDBPyConnection = field(
        init=False, default=Factory(lambda x: duckdb.connect(x.cache), takes_self=True)
    )
    seen: set[str] = field(factory=set)
    curies: set[str] = field(factory=set)
    columns: tuple[str, ...] = field(
        init=False,
        default=Factory(
            lambda x: tuple(["obo_id", "label"] + list(x.fields.keys())),
            takes_self=True,
        ),
    )

    def __getitem__(self, curie: str) -> pl.DataFrame:
        return self.db.execute("SELECT * FROM t WHERE obo_id = ?", [curie]).pl()

    def pl(self) -> pl.DataFrame:
        return self.db.execute("SELECT * FROM t").pl()

    def __attrs_post_init__(self):
        column_str = ",".join([f"{k} {v}" for k, v in self.fields.items()])
        self.db.sql(f"""
        CREATE TABLE IF NOT EXISTS t (
        obo_id VARCHAR PRIMARY KEY,
        label VARCHAR,
        {column_str}
        )
        """)
        self.db.sql(f"""
        CREATE TABLE IF NOT EXISTS lookups (
        query VARCHAR PRIMARY KEY,
        obo_id VARCHAR[],
        )
        """)
        self.seen |= set(
            self.db.execute("SELECT query FROM lookups").fetchnumpy()["query"]
        )
        self.curies |= set(
            self.db.execute("SELECT obo_id FROM t").fetchnumpy()["obo_id"]
        )

    def search(self, query: str) -> pl.DataFrame:
        """Search for query

        See the docs for `/api/search`
        """
        if query in self.seen:
            return self.db.execute(
                """
            SELECT * FROM t
            JOIN (SELECT query, unnest(obo_id) AS obo_id FROM lookups) AS r
            ON t.obo_id = r.obo_id
            WHERE r.query = ?
            """,
                [query],
            ).pl()
        req = requests.get(
            f"{self.endpoint}/api/search",
            {
                "q": query,
                "ontology": self.ontologies,
                "fieldList": list(self.fields.keys()) + ["obo_id", "label"],
                "queryFields": self.query_fields,
            },
        )
        req.raise_for_status()
        result = req.json()["response"]["docs"]
        df: pl.DataFrame = pl.DataFrame(result).select(self.columns)
        lookups_df = (
            df.with_columns(pl.lit(query).alias("query"))
            .group_by("query")
            .agg(pl.col("obo_id"))
        ).select(["query", "obo_id"])
        df = df.filter(~pl.col("obo_id").is_in(self.curies))
        self.seen.add(query)
        if not df.is_empty():
            self.db.execute("INSERT INTO t SELECT * FROM df")
        self.db.execute("INSERT INTO lookups SELECT * FROM lookups_df")
        return df
