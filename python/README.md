# CIRCUST in Python

Python implementation of the CIRCUST algorithm as described in Larriba et al. (2023) and its supplementary material. It replicates the four stages of the R version in this repository (root folder).

The package is self-contained: it lives in this `python/` folder and does not depend on the R files.

## What it computes

CIRCUST does not estimate a single sample order but **K synchronized circular orders**, and summarizes the fits of the TOP genes using the median of the K parameters.

| Stage | What it does |
|---|---|
| 1. Preprocessing | Removes genes with more than 30% zeros and normalizes each gene to [-1, 1]. |
| 2. Preliminary order | CPCA with the 12 core genes. Detects and removes outliers. Synchronizes the order: ARNTL peak at π and direction fixed with DBP. |
| 3. TOP genes | ORI filter (R²_ORI ≥ 0.5), split of the circumference into intervals, FMM fit in each one (ω > 0.1, R²_FMM > 0.5) and selection of the best ones. Core genes are always included. |
| 4. Robust estimation | K random subsets of 2/3 of the TOP genes, one CPCA and one synchronization per subset, FMM fit of all TOP genes in each of the K orders, and median of the K parameters. |

## Installation

Python 3.10 or later.

```bash
cd python
python3 -m venv .venv && source .venv/bin/activate     # Windows: .venv\Scripts\activate
pip install -e .                                        # installs the circust package
pip install -r requirements.txt                        # includes Numba (optional but recommended)
```

`requirements.txt` includes Numba, which speeds up the FMM fit. Without Numba the result is the same but slower.

## Usage

```bash
python scripts/run_circust.py --data data/matrixIn_mini.parquet -o results
```

The bundled example (412 genes × 120 samples) takes about 10 seconds. The full `matrixIn` (56,200 genes × 479 samples) takes about 1 minute with 24 cores.

The input matrix has genes in rows and samples in columns, with gene symbols as identifiers. CSV, TSV, Parquet (`gene_id` column), Excel and GCT are supported.

### Arguments

| Argument | Default | Description |
|---|---|---|
| `--data` | required | Expression matrix. |
| `--gene-column` | `gene_id` for Parquet, first column otherwise | Column containing the gene symbols. |
| `--core-genes` | `circust` | Preset (`circust`, `zhang`, `ruben`) or comma-separated list. Requires ARNTL and DBP. CRY1 is recommended. |
| `-o`, `--output` | `output_circust` | Output folder. |
| `--n-reps` | 5 | Number K of orders. |
| `--seed` | 42 | Seed for the subset sampling. With the same seed and the same data the result is reproducible. |
| `--values` | `normalized` | Values of the ordered matrices: `normalized` (to [-1, 1]) or `original` (as provided). |
| `--matrix-format` | `csv` | Format of the full matrices: `csv` or `parquet`. With `matrixIn` each matrix takes about 90 MB in Parquet and considerably more in CSV. |
| `--n-intervals` | 8 | Number of intervals the circumference is split into in Stage 3. |
| `--genes-per-interval` | 50 | Number of genes admitted per interval in Stage 3. |
| `--jobs` | -1 | Parallel processes (-1 uses all). |
| `--quiet` | | Do not print progress. |

## Outputs

For each k = 1, …, K:

| File | Content |
|---|---|
| `sample_order_k<k>.csv` | Order k of the samples: `sample_name`, `position`, `circular_phase` (radians in [0, 2π)). |
| `ordered_matrix_k<k>` | **Full** matrix (all genes that pass preprocessing, without the removed samples) with columns in order k. |
| `fitted_top_k<k>.csv` | FMM fitted values of each TOP gen, with columns in order k. |

In addition:

| File | Content |
|---|---|
| `circust_atlas.csv` | Tissue atlas: one row per TOP gene with the **median of the K** parameters `M, A, alpha, beta, omega, t_U, t_L` and `r2_fmm`, and `classification` (`day` if the peak is in [0, π), `night` otherwise). |
| `parameters_by_k.csv` | The parameters of each gene in each of the K fits (K rows per gene). |
| `outliers.csv` | Removed samples and the criterion that removed them. |
| `top_genes.csv` | TOP genes with their interval and their values in the preliminary order. |
| `manifest.json` | Configuration, seed, orientation of each order and timings. |

Phases are given in radians and are not converted to hours, because that would require assuming a 24 h period that the method does not verify. The angles in the atlas (`alpha`, `beta`, `t_U`, `t_L`) are summarized with the **circular median** and the rest with the ordinary median, because the ordinary median is not valid across the cut between 2π and 0.


## Using it from Python

```python
from circust.core_genes import CoreGeneSelector
from circust.cpca import CPCA
from circust.preprocessing import Preprocessor, load_expression_matrix
from circust.robust_order import RobustOrderEstimator
from circust.synchronizer import CircularSynchronizer
from circust.top_genes import TopGeneSelector

raw  = load_expression_matrix("data/matrixIn_mini.parquet", gene_column="gene_id")
prep = Preprocessor(verbose=False).run(raw)
core = CoreGeneSelector.from_string(None, verbose=False).select(prep.expr_norm).genes

cpca = CPCA(core_genes=core, verbose=False).run(prep.expr_norm)           # Stage 2.1
sync = CircularSynchronizer(verbose=False).run(                           # Stage 2.2
    cpca.expr_norm_final, cpca.circular_scale,
    cpca.core_genes_found, cpca.fmm_fits_final)
top = TopGeneSelector(verbose=False).run(                                 # Stage 3
    sync.expr_ordered, sync.circular_scale, cpca.core_genes_found)
robust = RobustOrderEstimator(n_reps=5, seed=1, verbose=False).run(       # Stage 4
    top.candidate_matrix, top.circular_scale,
    cpca.core_genes_found, forced_genes=top.forced_genes)

robust.atlas                # median of the K parameters per TOP gene
robust.parameters           # K rows per gene
robust.ordered_names(0)     # samples in order 1
robust.fitted[0]            # fitted values of the TOP genes in order 1
```


## References

- Larriba, Y., Mason, I. C., Saxena, R., Scheer, F. A. J. L. and Rueda, C. (2023). CIRCUST: A novel methodology for temporal order reconstruction of molecular rhythms. *PLOS Computational Biology*, 19(9), e1011510. https://doi.org/10.1371/journal.pcbi.1011510
- Rueda, C., Larriba, Y. and Peddada, S. D. (2019). Frequency Modulated Möbius model accurately predicts rhythmic signals in biological and physical sciences. *Scientific Reports*, 9, 18701.
- Larriba, Y., Rueda, C., Fernández, M. A. and Peddada, S. D. (2019). Order restricted inference in chronobiology. *Statistics in Medicine*.
