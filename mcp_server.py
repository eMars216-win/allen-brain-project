"""Read-only MCP tools for the Allen Brain Atlas lateralisation warehouse."""

from __future__ import annotations

import math
import os
from datetime import date, datetime
from decimal import Decimal
from pathlib import Path
from typing import Any

import duckdb
from mcp.server import MCPServer


SERVER_DIR = Path(__file__).resolve().parent
DATABASE_PATH = Path(
    os.environ.get("ALLEN_DUCKDB_PATH", SERVER_DIR / "allen_lateralisation.duckdb")
).expanduser()

server = MCPServer(
    "Allen Brain Atlas Lateralisation",
    instructions=(
        "Read-only access to an Allen Human Brain Atlas lateralisation DuckDB "
        "warehouse. Build the database with build_duckdb_from_csv.py before querying."
    ),
)


def _connect() -> duckdb.DuckDBPyConnection:
    path = DATABASE_PATH.resolve()
    if not path.is_file():
        raise FileNotFoundError(
            f"DuckDB database not found at {path}. Build it with: "
            "python build_duckdb_from_csv.py --raw-dir raw --database "
            "allen_lateralisation.duckdb"
        )
    return duckdb.connect(str(path), read_only=True)


def _json_value(value: Any) -> Any:
    if value is None:
        return None
    if isinstance(value, (datetime, date)):
        return value.isoformat()
    if isinstance(value, Decimal):
        return float(value)
    if isinstance(value, bytes):
        return value.decode("utf-8", errors="replace")
    if isinstance(value, float) and not math.isfinite(value):
        return None
    if hasattr(value, "item"):
        return _json_value(value.item())
    return value


def _fetch(sql: str, parameters: list[Any] | None = None) -> list[dict[str, Any]]:
    with _connect() as connection:
        cursor = connection.execute(sql, parameters or [])
        columns = [column[0] for column in cursor.description]
        return [
            {column: _json_value(value) for column, value in zip(columns, row)}
            for row in cursor.fetchall()
        ]


def _check_limit(limit: int) -> int:
    if not 1 <= limit <= 100:
        raise ValueError("limit must be between 1 and 100")
    return limit


@server.tool()
def dataset_status() -> dict[str, Any]:
    """Show whether the warehouse exists and which tables and views it contains."""
    if not DATABASE_PATH.expanduser().is_file():
        return {
            "database_path": str(DATABASE_PATH.expanduser().resolve()),
            "ready": False,
            "build_command": (
                "python build_duckdb_from_csv.py --raw-dir raw "
                "--database allen_lateralisation.duckdb"
            ),
        }

    objects = _fetch(
        """
        select table_name, table_type
        from information_schema.tables
        where table_schema = 'main'
        order by table_type, table_name
        """
    )
    names = {row["table_name"] for row in objects}
    required = {
        "fact_expression",
        "dim_region_pairs",
        "dim_gene",
        "v_expression_by_structure",
        "v_lateralisation",
        "v_lateralisation_summary",
    }
    return {
        "database_path": str(DATABASE_PATH.expanduser().resolve()),
        "ready": required.issubset(names),
        "missing_objects": sorted(required - names),
        "objects": objects,
    }


@server.tool()
def find_genes(query: str, limit: int = 20) -> list[dict[str, Any]]:
    """Find genes by Allen ID, Entrez ID, symbol, name, or probe name."""
    limit = _check_limit(limit)
    query = query.strip()
    if not query or len(query) > 120:
        raise ValueError("query must contain between 1 and 120 characters")
    pattern = f"%{query}%"
    return _fetch(
        """
        with gene_annotations as (
            select
                gene_id,
                min(try_cast(entrez_id as bigint)) as entrez_id,
                any_value(nullif(gene_symbol, '')) filter (
                    where try_cast(entrez_id as bigint) is not null
                ) as gene_symbol,
                any_value(nullif(gene_name, '')) filter (
                    where try_cast(entrez_id as bigint) is not null
                ) as gene_name,
                string_agg(distinct probe_name, ', ' order by probe_name) as probe_names
            from raw_probes
            group by gene_id
        )
        select gene_id, entrez_id, gene_symbol, gene_name, probe_names
        from gene_annotations
        where cast(gene_id as varchar) = ?
           or cast(entrez_id as varchar) = ?
           or gene_symbol ilike ?
           or gene_name ilike ?
           or probe_names ilike ?
        order by gene_symbol nulls last, gene_id
        limit ?
        """,
        [query, query, pattern, pattern, pattern, limit],
    )


@server.tool()
def top_lateralised_genes(
    direction: str = "any",
    min_donors: int = 2,
    limit: int = 20,
    mapped_only: bool = True,
) -> list[dict[str, Any]]:
    """Rank mapped genes by average absolute LI across paired regions.

    direction may be "any", "left", or "right". Positive LI means
    left-lateralised; negative LI means right-lateralised.
    Set mapped_only=False to include Allen probe features without an Entrez ID.
    """
    limit = _check_limit(limit)
    direction = direction.lower().strip()
    if direction not in {"any", "left", "right"}:
        raise ValueError('direction must be "any", "left", or "right"')
    if min_donors < 1:
        raise ValueError("min_donors must be at least 1")

    sign_filter = {
        "any": "true",
        "left": "mean(s.mean_lateralisation_index) > 0",
        "right": "mean(s.mean_lateralisation_index) < 0",
    }[direction]
    return _fetch(
        f"""
        with per_gene as (
            select
                s.gene_id,
                avg(s.mean_lateralisation_index) as mean_lateralisation_index,
                avg(abs(s.mean_lateralisation_index)) as mean_absolute_li,
                avg(s.donor_li_stddev) as mean_donor_li_stddev,
                count(*) as region_pair_count,
                min(s.donor_count) as minimum_donor_count
            from v_lateralisation_summary s
            group by s.gene_id
            having min(s.donor_count) >= ?
               and {sign_filter}
        )
        , gene_annotations as (
            select
                gene_id,
                min(try_cast(entrez_id as bigint)) as entrez_id,
                any_value(nullif(gene_symbol, '')) filter (
                    where try_cast(entrez_id as bigint) is not null
                ) as gene_symbol,
                any_value(nullif(gene_name, '')) filter (
                    where try_cast(entrez_id as bigint) is not null
                ) as gene_name,
                string_agg(distinct probe_name, ', ' order by probe_name) as probe_names
            from raw_probes
            group by gene_id
        )
        select
            p.*,
            a.entrez_id,
            a.gene_symbol,
            a.gene_name,
            a.probe_names
        from per_gene p
        left join gene_annotations a on a.gene_id = p.gene_id
        where (? = false or a.entrez_id is not null)
        order by p.mean_absolute_li desc, p.gene_id
        limit ?
        """,
        [min_donors, mapped_only, limit],
    )


@server.tool()
def region_lateralisation(
    region: str, gene_query: str = "", limit: int = 30
) -> list[dict[str, Any]]:
    """Return gene lateralisation summaries for region-pair names matching text."""
    limit = _check_limit(limit)
    region = region.strip()
    if not region or len(region) > 120:
        raise ValueError("region must contain between 1 and 120 characters")
    region_pattern = f"%{region}%"
    gene_pattern = f"%{gene_query.strip()}%"
    return _fetch(
        """
        select
            rp.region_pair_id,
            rp.parent_structure_id,
            rp.left_name,
            rp.right_name,
            s.gene_id,
            g.acronym,
            g.name,
            s.mean_lateralisation_index,
            s.donor_li_stddev,
            s.donor_count
        from v_lateralisation_summary s
        join dim_region_pairs rp on rp.region_pair_id = s.region_pair_id
        left join (
            select gene_id, any_value(acronym) as acronym,
                   any_value(name) as name
            from dim_gene group by gene_id
        ) g on g.gene_id = s.gene_id
        where (rp.left_name ilike ? or rp.right_name ilike ?)
          and (? = '' or cast(s.gene_id as varchar) = ?
               or g.acronym ilike ? or g.name ilike ?)
        order by abs(s.mean_lateralisation_index) desc, s.gene_id
        limit ?
        """,
        [
            region_pattern,
            region_pattern,
            gene_query.strip(),
            gene_query.strip(),
            gene_pattern,
            gene_pattern,
            limit,
        ],
    )


@server.tool()
def donor_variability(region: str = "", limit: int = 20) -> list[dict[str, Any]]:
    """Rank region pairs by average across-gene donor variability in LI."""
    limit = _check_limit(limit)
    region_pattern = f"%{region.strip()}%"
    return _fetch(
        """
        select
            rp.region_pair_id,
            rp.parent_structure_id,
            rp.left_name,
            rp.right_name,
            avg(s.donor_li_stddev) as mean_gene_li_stddev,
            max(s.donor_li_stddev) as maximum_gene_li_stddev,
            count(*) as gene_count
        from v_lateralisation_summary s
        join dim_region_pairs rp on rp.region_pair_id = s.region_pair_id
        where (? = '' or rp.left_name ilike ? or rp.right_name ilike ?)
        group by rp.region_pair_id, rp.parent_structure_id,
                 rp.left_name, rp.right_name
        order by mean_gene_li_stddev desc nulls last
        limit ?
        """,
        [region.strip(), region_pattern, region_pattern, limit],
    )


if __name__ == "__main__":
    server.run()
