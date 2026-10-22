"""Build the Allen lateralisation DuckDB database from downloaded CSV files."""

from __future__ import annotations

import argparse
from pathlib import Path

from allen_lateralisation import write_raw_csv_tables, write_to_duckdb
from allen_transform import (
    build_fact_expression_from_csv,
    build_gene_dimension,
    build_region_pairs_from_ontology,
    load_raw_csvs,
    run_diagnostics,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--raw-dir",
        default="raw",
        help="Directory containing MicroarrayExpression.csv, Probes.csv, SampleAnnot.csv, and Ontology.csv.",
    )
    parser.add_argument(
        "--database",
        default="allen_lateralisation.duckdb",
        help="Output DuckDB database path.",
    )
    parser.add_argument(
        "--max-probes",
        type=int,
        default=None,
        help="Optional probe-row limit for a quick smoke test.",
    )
    parser.add_argument(
        "--donor-id",
        default=None,
        help="Donor ID to use when SampleAnnot.csv has no donor ID column.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    raw_dir = Path(args.raw_dir)

    raw_tables = load_raw_csvs(raw_dir)
    probes = raw_tables["probes"]
    samples = raw_tables["samples"]
    ontology = raw_tables["ontology"]

    region_pairs = build_region_pairs_from_ontology(ontology)
    genes = build_gene_dimension(probes)
    fact_expression = build_fact_expression_from_csv(
        raw_dir,
        probes,
        samples,
        max_probes=args.max_probes,
        donor_id=args.donor_id,
    )

    diagnostics = run_diagnostics(fact_expression)
    write_to_duckdb(
        args.database,
        fact_expression=fact_expression,
        region_pairs=region_pairs,
        genes=genes,
    )
    write_raw_csv_tables(args.database, raw_tables=raw_tables)

    print(f"Wrote {args.database}")
    print(f"Region pairs: {len(region_pairs):,}")
    print(f"Fact rows: {len(fact_expression):,}")
    print(f"Diagnostics: {diagnostics}")


if __name__ == "__main__":
    main()
