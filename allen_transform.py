"""Transformation helpers for downloaded Allen CSV files."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd


CSV_ALIASES = {
    "expression": ("MicroarrayExpression.csv", "microarray_expression.csv"),
    "probes": ("Probes.csv", "probes.csv"),
    "samples": ("SampleAnnot.csv", "sample_annot.csv", "samples.csv"),
    "ontology": ("Ontology.csv", "ontology.csv", "structures.csv"),
}


def normalize_column_name(column: Any) -> str:
    """Normalize Allen CSV headers for predictable downstream references."""
    return str(column).strip().lower().replace(" ", "_").replace("-", "_")


def find_csv(raw_dir: str | Path, logical_name: str) -> Path:
    """Find a known Allen CSV by logical name."""
    raw_path = Path(raw_dir)
    for filename in CSV_ALIASES[logical_name]:
        candidate = raw_path / filename
        if candidate.exists():
            return candidate
    expected = ", ".join(CSV_ALIASES[logical_name])
    raise FileNotFoundError(f"Could not find {logical_name} CSV. Expected one of: {expected}")


def read_metadata_csv(path: str | Path) -> pd.DataFrame:
    """Read an Allen metadata CSV and normalize its column names."""
    frame = pd.read_csv(path)
    frame.columns = [normalize_column_name(column) for column in frame.columns]
    return frame


def first_existing(columns: list[str], candidates: tuple[str, ...]) -> str | None:
    """Return the first candidate column present in a list of columns."""
    for candidate in candidates:
        if candidate in columns:
            return candidate
    return None


def require_column(frame: pd.DataFrame, candidates: tuple[str, ...], label: str) -> str:
    """Resolve a required column from several known Allen header variants."""
    column = first_existing(list(frame.columns), candidates)
    if column is None:
        raise ValueError(
            f"Could not identify {label}. Tried {candidates}; found {list(frame.columns)}"
        )
    return column


def load_raw_csvs(raw_dir: str | Path) -> dict[str, pd.DataFrame]:
    """Load the downloaded Allen metadata CSVs into normalized DataFrames."""
    return {
        "probes": read_metadata_csv(find_csv(raw_dir, "probes")),
        "samples": read_metadata_csv(find_csv(raw_dir, "samples")),
        "ontology": read_metadata_csv(find_csv(raw_dir, "ontology")),
    }


def build_region_pairs_from_ontology(ontology: pd.DataFrame) -> pd.DataFrame:
    """Build dim_region_pairs from a downloaded ontology CSV."""
    structure_id = require_column(ontology, ("id", "structure_id"), "structure id")
    parent_id = require_column(
        ontology,
        ("parent_structure_id", "parent_id"),
        "parent structure id",
    )
    hemisphere_id = require_column(
        ontology, ("hemisphere_id", "hemisphere"), "hemisphere id"
    )
    graph_order = require_column(ontology, ("graph_order",), "graph order")
    name = first_existing(list(ontology.columns), ("name", "structure_name", "label"))

    structures = ontology.copy()
    structures = structures.rename(
        columns={
            structure_id: "id",
            parent_id: "parent_structure_id",
            graph_order: "graph_order",
        }
    )
    if hemisphere_id != "hemisphere_id":
        hemisphere_values = structures[hemisphere_id].astype("string").str.strip().str.upper()
        hemisphere_mapping = {
            "L": 1,
            "LEFT": 1,
            "R": 2,
            "RIGHT": 2,
            "B": 3,
            "BILATERAL": 3,
            "M": 3,
            "MIDLINE": 3,
        }
        numeric_hemispheres = pd.to_numeric(hemisphere_values, errors="coerce")
        structures["hemisphere_id"] = numeric_hemispheres.fillna(
            hemisphere_values.map(hemisphere_mapping)
        ).astype("Int64")
    else:
        structures["hemisphere_id"] = pd.to_numeric(
            structures[hemisphere_id], errors="coerce"
        ).astype("Int64")
    if name and name != "name":
        structures = structures.rename(columns={name: "name"})
    if "name" not in structures.columns:
        structures["name"] = structures["id"].astype(str)

    left = structures[structures["hemisphere_id"].astype("Int64") == 1].copy()
    right = structures[structures["hemisphere_id"].astype("Int64") == 2].copy()
    left = left.sort_values(["parent_structure_id", "graph_order"])
    right = right.sort_values(["parent_structure_id", "graph_order"])
    left["pair_index"] = left.groupby("parent_structure_id").cumcount() + 1
    right["pair_index"] = right.groupby("parent_structure_id").cumcount() + 1

    pairs = left.merge(
        right,
        on=["parent_structure_id", "pair_index"],
        suffixes=("_left", "_right"),
    )
    return pd.DataFrame(
        {
            "region_pair_id": range(1, len(pairs) + 1),
            "parent_structure_id": pairs["parent_structure_id"],
            "pair_index": pairs["pair_index"],
            "left_structure_id": pairs["id_left"],
            "right_structure_id": pairs["id_right"],
            "left_name": pairs["name_left"],
            "right_name": pairs["name_right"],
        }
    )


def build_gene_dimension(probes: pd.DataFrame) -> pd.DataFrame:
    """Build dim_gene from the downloaded Probes CSV."""
    probe_id = require_column(probes, ("probe_id", "id"), "probe id")
    gene_id = first_existing(list(probes.columns), ("gene_id", "entrez_id"))
    acronym = first_existing(list(probes.columns), ("gene_symbol", "acronym", "symbol"))
    name = first_existing(list(probes.columns), ("gene_name", "name"))

    dim = pd.DataFrame()
    dim["probe_id"] = probes[probe_id]
    dim["gene_id"] = probes[gene_id] if gene_id else probes[probe_id]
    dim["acronym"] = probes[acronym] if acronym else pd.NA
    dim["name"] = probes[name] if name else pd.NA
    return dim.drop_duplicates()


def build_fact_expression_from_csv(
    raw_dir: str | Path,
    probes: pd.DataFrame,
    samples: pd.DataFrame,
    *,
    max_probes: int | None = None,
    donor_id: str | int | None = None,
) -> pd.DataFrame:
    """Convert Allen's expression matrix CSV into long fact-expression rows.

    Allen's matrix has probe IDs in its first column and one expression column
    per sample after that. SampleAnnot.csv lists samples in the same order as
    those expression columns, so the source column position is retained as a
    temporary key and joined to the annotation row position. Donor ID is taken
    from SampleAnnot when present, or from the explicit `donor_id` argument.
    """
    expression_path = find_csv(raw_dir, "expression")
    expression = pd.read_csv(expression_path, header=None, nrows=max_probes)

    probe_id = require_column(probes, ("probe_id", "id"), "probe id")
    gene_id = first_existing(list(probes.columns), ("gene_id", "entrez_id"))
    sample_id = require_column(samples, ("sample_id", "well_id", "id"), "sample id")
    structure_id = require_column(
        samples,
        ("structure_id", "ontology_structure_id", "structure_id_path"),
        "sample structure id",
    )
    sample_donor_column = first_existing(list(samples.columns), ("donor_id", "donor_name"))
    if sample_donor_column is None and donor_id is None:
        raise ValueError(
            "SampleAnnot.csv has no donor ID column. Supply donor_id explicitly "
            "(for example, donor_id='9861')."
        )

    expression_sample_count = expression.shape[1] - 1
    if expression_sample_count != len(samples):
        raise ValueError(
            "Expression/sample annotation mismatch: MicroarrayExpression.csv "
            f"contains {expression_sample_count:,} sample columns after its probe ID, "
            f"but SampleAnnot.csv contains {len(samples):,} rows. "
            "Use the matching complete SampleAnnot.csv for this donor."
        )

    # Allen's first expression-matrix field is the probe ID, not a measurement.
    # Use it as the join key so probe ordering differences cannot silently mislabel rows.
    expression = expression.rename(columns={0: "probe_id"})
    measurement_columns = list(expression.columns[1:])
    expression = expression.melt(
        id_vars=["probe_id"],
        value_vars=measurement_columns,
        var_name="sample_column",
        value_name="expression_value",
    )
    expression["sample_column"] = expression["sample_column"].astype(int) - 1

    probe_columns = [probe_id] + ([gene_id] if gene_id else [])
    probe_lookup = probes[probe_columns].drop_duplicates(subset=[probe_id]).copy()
    probe_lookup = probe_lookup.rename(columns={probe_id: "probe_id"})
    if not gene_id:
        probe_lookup["gene_id"] = probe_lookup["probe_id"]
    elif gene_id != "gene_id":
        probe_lookup = probe_lookup.rename(columns={gene_id: "gene_id"})
    expression = expression.merge(probe_lookup, on="probe_id", how="left", validate="many_to_one")

    sample_lookup = samples.reset_index(names="sample_column").copy()
    sample_lookup = sample_lookup.rename(
        columns={
            sample_id: "sample_id",
            structure_id: "structure_id",
        }
    )
    if sample_donor_column is not None:
        sample_lookup = sample_lookup.rename(columns={sample_donor_column: "donor_id"})
    else:
        sample_lookup["donor_id"] = str(donor_id)

    long = expression.merge(
        sample_lookup[["sample_column", "sample_id", "structure_id", "donor_id"]],
        on="sample_column",
        how="left",
        validate="many_to_one",
    )
    long["expression_value"] = pd.to_numeric(long["expression_value"], errors="coerce")
    return long[
        ["donor_id", "probe_id", "gene_id", "sample_id", "structure_id", "expression_value"]
    ]


def run_diagnostics(fact_expression: pd.DataFrame) -> dict[str, Any]:
    """Return basic grain and completeness diagnostics."""
    required = {"donor_id", "probe_id", "gene_id", "sample_id", "structure_id", "expression_value"}
    missing_columns = sorted(required - set(fact_expression.columns))
    if missing_columns:
        return {"missing_columns": missing_columns, "row_count": len(fact_expression)}

    duplicate_grain = fact_expression.duplicated(
        ["donor_id", "probe_id", "sample_id", "structure_id"]
    ).sum()

    return {
        "row_count": len(fact_expression),
        "donor_count": fact_expression["donor_id"].nunique(dropna=True),
        "gene_count": fact_expression["gene_id"].nunique(dropna=True),
        "structure_count": fact_expression["structure_id"].nunique(dropna=True),
        "null_expression_count": fact_expression["expression_value"].isna().sum(),
        "duplicate_grain_count": int(duplicate_grain),
    }
