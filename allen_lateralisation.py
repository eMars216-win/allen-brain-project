"""Warehouse helpers for the Allen lateralisation pipeline."""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from typing import Iterable

import duckdb
import numpy as np
import pandas as pd


LEFT_HEMISPHERE_ID = 1
RIGHT_HEMISPHERE_ID = 2


@dataclass(frozen=True)
class HemispherePair:
    parent_structure_id: int
    pair_index: int
    left_structure_id: int
    right_structure_id: int
    left_name: str
    right_name: str


def build_hemisphere_pairs(structures: Iterable[dict]) -> list[HemispherePair]:
    """Pair left/right structures by parent and sorted graph order."""
    by_parent: dict[int, dict[int, list[dict[str, Any]]]] = defaultdict(
        lambda: defaultdict(list)
    )

    for structure in structures:
        parent_id = structure.get("parent_structure_id")
        hemisphere_id = structure.get("hemisphere_id")
        if parent_id is None or hemisphere_id not in {
            LEFT_HEMISPHERE_ID,
            RIGHT_HEMISPHERE_ID,
        }:
            continue
        by_parent[int(parent_id)][int(hemisphere_id)].append(structure)

    pairs: list[HemispherePair] = []
    for parent_id, hemispheres in sorted(by_parent.items()):
        left = sorted(
            hemispheres.get(LEFT_HEMISPHERE_ID, []),
            key=lambda row: row.get("graph_order", 0),
        )
        right = sorted(
            hemispheres.get(RIGHT_HEMISPHERE_ID, []),
            key=lambda row: row.get("graph_order", 0),
        )
        if len(left) != len(right):
            continue

        for pair_index, (left_row, right_row) in enumerate(zip(left, right), start=1):
            pairs.append(
                HemispherePair(
                    parent_structure_id=parent_id,
                    pair_index=pair_index,
                    left_structure_id=int(left_row["id"]),
                    right_structure_id=int(right_row["id"]),
                    left_name=str(left_row.get("name", "")),
                    right_name=str(right_row.get("name", "")),
                )
            )

    return pairs


def compute_lateralisation_index(expression_by_structure: pd.DataFrame) -> pd.DataFrame:
    """Compute LI from left/right expression rows already aligned by pair."""
    required = {
        "gene_id",
        "donor_id",
        "region_pair_id",
        "left_expression",
        "right_expression",
    }
    missing = required - set(expression_by_structure.columns)
    if missing:
        raise ValueError(f"Missing required columns: {sorted(missing)}")

    df = expression_by_structure.copy()
    denominator = df["left_expression"] + df["right_expression"]
    df["lateralisation_index"] = np.where(
        denominator != 0,
        (df["left_expression"] - df["right_expression"]) / denominator,
        np.nan,
    )
    return df


def write_to_duckdb(
    database_path: str,
    *,
    fact_expression: pd.DataFrame,
    region_pairs: pd.DataFrame,
    genes: pd.DataFrame,
) -> None:
    """Write core analytical tables and views to DuckDB."""
    with duckdb.connect(database_path) as con:
        con.register("fact_expression_df", fact_expression)
        con.register("region_pairs_df", region_pairs)
        con.register("genes_df", genes)

        con.execute("create or replace table fact_expression as select * from fact_expression_df")
        con.execute("create or replace table dim_region_pairs as select * from region_pairs_df")
        con.execute("create or replace table dim_gene as select * from genes_df")

        con.execute(
            """
            create or replace view v_expression_by_structure as
            select
                gene_id,
                donor_id,
                structure_id,
                avg(expression_value) as mean_expression,
                count(*) as measurement_count
            from fact_expression
            group by gene_id, donor_id, structure_id
            """
        )


def write_raw_csv_tables(
    database_path: str,
    *,
    raw_tables: dict[str, pd.DataFrame],
) -> None:
    """Write raw downloaded CSVs into DuckDB staging tables."""
    with duckdb.connect(database_path) as con:
        for table_name, frame in raw_tables.items():
            registered_name = f"{table_name}_df"
            con.register(registered_name, frame)
            con.execute(
                f"create or replace table raw_{table_name} as select * from {registered_name}"
            )

        con.execute(
            """
            create or replace view v_lateralisation as
            select
                rp.region_pair_id,
                rp.left_structure_id,
                rp.right_structure_id,
                l.gene_id,
                l.donor_id,
                l.mean_expression as left_expression,
                r.mean_expression as right_expression,
                case
                    when l.mean_expression + r.mean_expression = 0 then null
                    else (l.mean_expression - r.mean_expression)
                         / (l.mean_expression + r.mean_expression)
                end as lateralisation_index
            from dim_region_pairs rp
            join v_expression_by_structure l
              on l.structure_id = rp.left_structure_id
            join v_expression_by_structure r
              on r.structure_id = rp.right_structure_id
             and r.gene_id = l.gene_id
             and r.donor_id = l.donor_id
            """
        )

        con.execute(
            """
            create or replace view v_lateralisation_summary as
            select
                gene_id,
                region_pair_id,
                avg(lateralisation_index) as mean_lateralisation_index,
                stddev_samp(lateralisation_index) as donor_li_stddev,
                count(*) as donor_count
            from v_lateralisation
            group by gene_id, region_pair_id
            """
        )
