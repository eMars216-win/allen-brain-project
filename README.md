# Allen Brain Atlas Lateralisation Pipeline

A data engineering project exploring gene expression asymmetry across
the left and right hemispheres of the human brain, using the Allen
Human Brain Atlas API as the primary data source.

Built at the intersection of data engineering, neuroscience, and
philosophy of mind.

---

## Motivation

The human brain is roughly bilaterally symmetric in structure but
demonstrably asymmetric in function. Language production is
left-lateralised. Holistic, contextual attention tends to be
right-lateralised. The anatomical and functional bases of this
asymmetry are well-documented at the level of brain imaging and
lesion studies, but the molecular underpinnings — which genes are
differentially expressed between hemispheres, in which regions, and
with what functional implications — are less systematically explored.

This project approaches that question as a data engineering problem:
build a pipeline that ingests gene expression data from the Allen Human
Brain Atlas, computes lateralisation indices across all hemisphere-paired
brain regions, and produces an analytical layer queryable in SQL.

The intellectual motivation draws on Iain McGilchrist's hemispheric
framework (*The Master and His Emissary*, 2009), which argues that the
left and right hemispheres differ not merely in task specialisation but
in fundamental mode of attention — analytic and instrumental vs holistic
and contextual. If this framework has molecular correlates, they should
be visible in the gene expression asymmetry data. This pipeline is a
partial empirical test of that hypothesis.

---

## Research Questions

**Primary:**
- Which genes show the strongest and most consistent lateralisation
  across donors in the Allen Human Brain Atlas?
- Does the pattern of gene expression asymmetry respect anatomical
  boundaries, or do expression-similar regions cross anatomical lines?
- Do right-lateralised genes cluster into functional categories
  associated with McGilchrist's characterisation of right hemisphere
  function (social cognition, embodied processing, contextual attention)?

**Secondary:**
- Does individual variation in lateralisation indices correlate across
  genes — suggesting co-regulated lateralisation mechanisms?
- Which brain regions show the highest inter-donor variability in
  lateralisation — pointing toward individually variable asymmetry?

**Planned extension:**
- Phase 2: Replace scalar lateralisation index with kernel-based
  distributional comparison (Maximum Mean Discrepancy) to detect
  distributional differences beyond mean expression shift.

---

## Data Source

**Allen Human Brain Atlas**
- 6 post-mortem human brain donors
- ~900 brain regions per donor
- ~20,000 genes measured via microarray
- REST API: api.brain-map.org
- Free, no authentication required
- Documentation: help.brain-map.org

**Key API finding (discovered during development):**
The Human Brain Atlas encodes hemisphere via the `hemisphere_id` field
(1=left, 2=right, 3=bilateral/midline), not via name string suffixes.
Hemisphere pairs are matched by shared `parent_structure_id` and sorted
by `graph_order` within each parent group. 313 valid hemisphere pairs
identified across 1839 total structures. All 313 shared parents have
exactly equal left and right child counts — zero mismatches, confirming
positional pairing by graph_order is reliable.

---

## Pipeline Architecture

```
Allen Brain Atlas API
        │
        ▼
  [Ingestion Layer]
  Python / requests
  Paginated API calls
  Retry + backoff handling
        │
        ▼
  [Raw Storage]
  MongoDB (document store)
  Raw JSON responses preserved
  Idempotent — safe to re-run
        │
        ▼
  [Transformation Layer]
  Python / Pandas
  Flatten nested JSON
  Resolve positional alignment
  (probe × sample × donor)
        │
        ▼
  [Analytical Warehouse]
  DuckDB
  Star schema:
    fact_expression
    dim_region_pairs
    dim_gene
  Analytical views:
    v_expression_by_structure
    v_lateralisation
    v_lateralisation_summary
        │
        ▼
  [Analysis]
  SQL queries
  Lateralisation index (LI)
  Gene ontology enrichment
  McGilchrist framework mapping
```

**Why MongoDB + DuckDB rather than one or the other:**
MongoDB stores raw API responses as documents before any transformation.
This preserves provenance — if the transformation logic changes, raw
data can be reprocessed without re-hitting the API. DuckDB provides
the analytical layer: columnar storage, fast aggregation, standard SQL
including window functions and CTEs. The two tools do different jobs
and are not in competition.

---

## Lateralisation Index

For each gene × brain region pair × donor:

```
LI = (left_expression - right_expression) / (left_expression + right_expression)
```

- LI = +1: exclusively left-hemisphere expressed
- LI = -1: exclusively right-hemisphere expressed
- LI =  0: perfectly symmetric

LI is computed per donor then averaged across donors.
Standard deviation across donors flags regions of high individual
variability — analytically interesting in their own right.

---

## Key Schema Design Decisions

**Hemisphere pairing via parent_structure_id + graph_order:**
The Allen ontology is a tree. Each brain region has a parent. Left and
right homologous structures share the same parent. Within each parent,
children are sorted by graph_order and matched positionally. This
emerged from direct API exploration — the initial name-string matching
approach was incorrect and would have silently dropped data.

**Bridge table for ontology hierarchy:**
The Allen Brain Atlas brain region ontology is a tree 8+ levels deep.
Querying "all descendants of the hippocampus" requires traversing this
tree. A bridge table (ancestor_id, descendant_id, depth) pre-computes
all ancestor-descendant relationships at build time, converting
recursive tree traversal into a simple SQL join at query time.
Trade-off: build-time cost for query-time simplicity. Correct choice
for an ontology that changes rarely.

**Grain of fact_expression:**
One row = one probe × one sampling site × one donor measurement.
Aggregation (mean across probes, mean across sites, mean across donors)
is deferred to SQL views rather than baked into the transformation.
This keeps the raw data honest and makes aggregation decisions
transparent and changeable.

---

## Repository Structure

```
allen-lateralisation/
│
├── README.md                    ← this file
├── requirements.txt             ← Python dependencies
│
├── allen_lateralisation.py      ← ingestion pipeline
│                                   fetch_structures()
│                                   build_hemisphere_pairs()
│                                   fetch_donors()
│                                   fetch_probes()
│                                   fetch_and_store_expression()
│                                   compute_lateralisation_index()
│                                   write_to_duckdb()
│
├── allen_transform.py           ← transformation layer
│                                   flatten_expression_doc()
│                                   flatten_all_documents()
│                                   write_fact_expression()
│                                   run_diagnostics()
│
├── exploration/                 ← Jupyter notebooks
│   ├── 01_api_exploration.ipynb ← API structure investigation
│   └── 02_hemisphere_pairs.ipynb← pairing logic development
│
└── queries/                     ← saved SQL analyses
    ├── top_lateralised_genes.sql
    ├── region_asymmetry.sql
    └── donor_variability.sql
```

---

## Current Status

| Component | Status |
|-----------|--------|
| API connectivity | Confirmed working |
| Structure ontology fetch | Working — 1839 structures |
| Hemisphere pair identification | Working — 313 pairs confirmed |
| Pairing logic (hemisphere_id + graph_order) | Implemented |
| Human microarray expression endpoint | In progress |
| MongoDB ingestion | Implemented, awaiting expression data |
| DuckDB transformation | Implemented, awaiting expression data |
| Lateralisation index computation | Implemented, awaiting expression data |
| SQL analytical views | Implemented, awaiting expression data |
| Initial gene analysis (BDNF, FOXP2) | Pending |
| McGilchrist GO enrichment analysis | Planned — Phase 2 |
| Kernel MMD distributional analysis | Planned — Phase 2 |

---

## Intellectual Context

This project sits at the intersection of several threads:

**Data engineering:** The pipeline architecture follows standard DE
patterns — raw ingestion to document store, transformation to
relational warehouse, analytical layer in SQL. The Allen Atlas's
hierarchical ontology and multi-modal data structure provide genuine
schema design challenges beyond tutorial-level exercises.

**Neuroscience:** The lateralisation question engages active research
on the molecular basis of hemispheric asymmetry. The planned kernel
extension (Phase 2) connects to topological data analysis methods
increasingly used in transcriptomics.

**Philosophy of mind:** McGilchrist's framework is the interpretive
lens through which the analytical results will be read. The pipeline
is a partial empirical test of whether the hemispheric functional
asymmetry he describes has correlates at the gene expression level.
The results will be interesting regardless of direction — confirmation,
disconfirmation, and anomaly are all informative.

---

## Dependencies

```
requests      # Allen API HTTP calls
pymongo       # MongoDB document storage
duckdb        # Analytical warehouse
pandas        # Tabular transformation
numpy         # Numerical operations
```                                                                                                                                                                          

Install:
```cmd
pip install requests pymongo duckdb pandas numpy
```

---

## References

- Allen Human Brain Atlas: https://human.brain-map.org
- API Documentation: http://help.brain-map.org/display/api/Allen+Brain+Atlas+API
- Hawrylycz et al., *An anatomically comprehensive atlas of the adult
  human brain transcriptome*, Nature 2012
- McGilchrist, *The Master and His Emissary*, Yale University Press 2009
- Gretton et al., *A Kernel Two-Sample Test*, JMLR 2012

---

*This project is a work in progress. The pipeline architecture and
analytical questions will evolve as the data reveals its structure.*

---

## Known Issues

### API Service
human_microarray_expression API service intermittently down (confirmed via Allen community forum, recurring as of June 2026)
