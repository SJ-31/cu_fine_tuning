#!/usr/bin/env ipython
from pathlib import Path

import duckdb
import polars as pl
import requests
from attrs import Factory, define, field
from pyhornedowl import PyIndexedOntology, open_ontology
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry


@define
class SearchOntology:
    """
    Class to search ontologies by exact label matching, from Owl files

    Modified from
    https://incenp.org/notes/2025/comparing-python-ontology-libraries.html
    """

    ont: PyIndexedOntology = field(
        converter=lambda x: x if isinstance(x, PyIndexedOntology) else open_ontology(x)
    )
    prefix_name: str  # e.g. MONDO
    prefix: str  # e.g. http://purl.obolibrary.org/obo/MONDO_
    synonym2iri: dict[str, str] = field(factory=dict)

    def __attrs_post_init__(self):
        self.ont.add_prefix_mapping(self.prefix_name, self.prefix)
        self.ont.add_prefix_mapping("rdfs", "http://www.w3.org/2000/01/rdf-schema#")
        self.ont.add_prefix_mapping(
            "oio", "http://www.geneontology.org/formats/oboInOwl#"
        )
        self.ont.build_indexes()
        klasses = [c for c in self.ont.get_classes() if c.startswith(self.prefix)]
        if not klasses:
            raise ValueError(
                f"No class in the ontology starts with prefix `{self.prefix}`"
            )
        for klass in klasses:
            for synonym in self.ont.get_annotations(klass, "oio:hasExactSynonym"):
                self.synonym2iri[synonym] = klass
        if not self.synonym2iri:
            print("WARNING: no synonyms available")

    @property
    def annotation_properties(self):
        """
        Return all available annotation properties in the ontology.

        These can be used with `get_annotations` after the prefix
        is registered
        """
        return self.ont.get_annotation_properties()

    def lookup(self, s: str, as_iri: bool = False) -> str:
        iri = self.ont.get_iri_for_label(s) or self.synonym2iri.get(s)
        curie = self.ont.get_id_for_iri(iri)
        if not iri or not curie:
            raise KeyError(f"`{s}` doesn't exist as a synonym or label in the ontology")
        if as_iri:
            return iri
        return curie


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
    timeout = 10
    session: requests.Session = field(factory=requests.Session)
    columns: tuple[str, ...] = field(
        init=False,
        default=Factory(
            lambda x: tuple(
                ["obo_id", "label", "ontology_name"] + list(x.fields.keys())
            ),
            takes_self=True,
        ),
    )

    def __getitem__(self, curie: str) -> dict:
        return (
            self.db.execute("SELECT * FROM t WHERE obo_id = ?", [curie])
            .pl()
            .rows_by_key("obo_id", named=True, unique=True)[curie]
        )

    def pl(self) -> pl.DataFrame:
        return self.db.execute("SELECT * FROM t").pl()

    def __attrs_post_init__(self):
        column_str = ",".join([f"{k} {v}" for k, v in self.fields.items()])
        self.db.sql(f"""
        CREATE TABLE IF NOT EXISTS t (
        obo_id VARCHAR PRIMARY KEY,
        label VARCHAR,
        ontology_name VARCHAR,
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
        retry = Retry(
            total=5,
            backoff_factor=0.5,
            backoff_jitter=0.1,
            status_forcelist=(429, 500, 502, 503, 504),
            allowed_methods=frozenset(["GET"]),
            raise_on_status=False,
        )
        adapter = HTTPAdapter(max_retries=retry)
        self.session.mount("https://", adapter)
        self.session.mount("http://", adapter)

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
        req = self.session.get(
            f"{self.endpoint}/api/search",
            {
                "q": query,
                "ontology": self.ontologies,
                "fieldList": list(self.fields.keys())
                + ["obo_id", "label", "ontology_name"],
                "queryFields": self.query_fields,
            },
            timeout=self.timeout,
        )
        req.raise_for_status()
        result = req.json()["response"]["docs"]
        df: pl.DataFrame = (
            pl.DataFrame(result)
            .select(self.columns)
            .filter(pl.col("ontology_name").is_in(self.ontologies))
        )
        lookups_df = (
            df.with_columns(pl.lit(query).alias("query"))
            .group_by("query")
            .agg(pl.col("obo_id"))
        ).select(["query", "obo_id"])
        df = df.filter(~pl.col("obo_id").is_in(self.curies))
        self.seen.add(query)
        if not df.is_empty():
            self.curies |= set(df["obo_id"])
            self.db.execute("INSERT INTO t SELECT * FROM df")
        self.db.execute("INSERT INTO lookups SELECT * FROM lookups_df")
        return df
