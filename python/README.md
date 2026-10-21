# CIRCUST — Python Implementation

Python implementation of CIRCUST for reconstructing circadian sample orderings from gene expression matrices.

This implementation was developed in collaboration with Samuel Ortega Bustos.

## 1. Repository Structure

The Python implementation is located in the `python/` directory of the [CIRCUST GitHub repository](https://github.com/yolandalago/CIRCUST).

```text
python/
├── README.md
├── requirements.txt
├── pyproject.toml
├── LICENSE
├── .gitignore
├── data/
│   └── .gitkeep
|   └── matrixIn.parquet
├── scripts/
│   └── run_circust.py
└── src/
    └── circust/
        ├── preprocessing.py
        ├── core_genes.py
        ├── cpca.py
        ├── synchronizer.py
        ├── top_genes.py
        ├── robust_order.py
        ├── outputs.py
        └── fitting/
```

**Important:** After downloading and extracting the GitHub repository, navigate to the `python/` directory containing `pyproject.toml`, `requirements.txt`, `scripts/`, and `src/`. All commands below must be executed **from this directory**, not from the root directory of the GitHub repository.

The original R implementation remains available in the repository's root directory.

## 2. Installation on Windows (PowerShell)

Install Python 3.10 or later. The version used for testing was Python 3.14.4 on Windows; not all intermediate versions have been tested.

Open PowerShell in the `python/` directory of the downloaded repository. For example:

```powershell
cd "C:\Users\YourUsername\Downloads\CIRCUST-main\python"
```

**Replace this path with the actual location of your downloaded repository.** Verify that you are in the correct directory:

```powershell
Test-Path .\pyproject.toml
Test-Path .\scripts\run_circust.py
```

**Both commands should return `True`.** If either returns `False`, check your working directory before proceeding.

Create a virtual environment and install the package and its dependencies:

```powershell
py -m venv .venv
.\.venv\Scripts\python.exe -m pip install --upgrade pip
.\.venv\Scripts\python.exe -m pip install -e .
.\.venv\Scripts\python.exe -m pip install -r requirements.txt
```

Running `pip install -e .` is necessary to ensure that `from circust...` imports work when executing `scripts/run_circust.py`. The `requirements.txt` file also installs Numba to accelerate computations. Activating the virtual environment is not required.

Verify the installation:

```powershell
.\.venv\Scripts\python.exe -c "import circust; print('CIRCUST imported successfully')"
.\.venv\Scripts\python.exe .\scripts\run_circust.py --help
```

If PowerShell reports that `py` is not recognized, try `python -m venv .venv` and continue using the commands beginning with `.\.venv\Scripts\python.exe`.

If you encounter `ModuleNotFoundError: No module named 'circust'`, verify that you are in the directory containing `pyproject.toml` and rerun:

```powershell
.\.venv\Scripts\python.exe -m pip install -e .
```

## 3. Input Data Preparation

### 3.1. Reference dataset: `matrixIn.parquet`

The reference dataset, `data/matrixIn.parquet`, is included in this repository. It contains the same gene expression matrix used in the original R implementation of CIRCUST, allowing users to reproduce the analysis in Python and compare the results obtained with both implementations.

CIRCUST takes a gene expression matrix as input, with **genes in rows and samples in columns**.

The reference dataset is `matrixIn.parquet`, which contains the same expression matrix used in the original R implementation of CIRCUST. It is provided to facilitate reproducibility and comparison between both implementations.

The expected data structure is:

| gene_id | Sample_1 | Sample_2 | ... | Sample_n |
|---|---|---|---|---|
| Gene_1 | 5.21 | 3.47 | ... | 4.18 |
| Gene_2 | 2.36 | 6.12 | ... | 3.95 |
| ... | ... | ... | ... | ... |
| Gene_G | 4.72 | 5.08 | ... | 2.61 |

- **Rows (genes):** each row represents a gene. Gene identifiers should be provided in a dedicated column named `gene_id` (e.g., `PER1`, `PER2`, `ARNTL`). If a different column name is used, specify it with `--gene-column COLUMN_NAME`.
- **Columns (samples):** each remaining column represents one biological sample. Column headers are used directly as sample identifiers; no specific naming convention is required. Sample identifiers should be unique and non-empty.
- **Expression values:** all sample columns must contain numeric gene expression values. Do not include sample metadata (e.g., age, sex, tissue, or batch) as additional columns.
- **Circadian core genes:** the input matrix must contain the core clock genes required by the CIRCUST pipeline, using gene identifiers compatible with those expected by the implementation.

No separate sample identifier column or sample metadata file is required. Sample identifiers are read directly from the column names.

The default gene identifier column is `gene_id`. If another name is used, specify it using `--gene-column COLUMN_NAME`.

CIRCUST performs its own preprocessing and normalization. **Do not manually normalize expression values to [-1, 1].**

### 3.2. Running the reference dataset

Place `matrixIn.parquet` in the `data/` directory and execute:

```powershell
.\.venv\Scripts\python.exe .\scripts\run_circust.py --data data/matrixIn.parquet --r-compatible --n-reps 5 --jobs -1 -o resultados_circust
```

This runs the complete Python workflow using the same input expression matrix as the original R implementation.

### 3.3. Other supported input formats (optional)

CIRCUST also supports **CSV, TSV, and Excel** files, provided they follow the same structure: gene identifiers in the first column and numeric expression values in the remaining sample columns.

For example, a standard CSV file should have the following structure:

```csv
gene_id,Sample_1,Sample_2,Sample_3
PER1,5.21,4.18,3.76
PER2,3.47,5.12,4.09
ARNTL,6.14,5.83,4.92
```

To run CIRCUST on a CSV file, replace the input path:

```powershell
.\.venv\Scripts\python.exe .\scripts\run_circust.py --data data/my_expression_matrix.csv --r-compatible --n-reps 5 --jobs -1 -o resultados_circust
```

CSV files should use commas as field separators and periods as decimal separators. Files using semicolons and decimal commas should be converted to standard CSV format before execution.

For Excel files, the additional dependency `openpyxl` may be required:

```powershell
.\.venv\Scripts\python.exe -m pip install openpyxl
```

**R data files:** `.rds` and `.RData` files are not read directly by the Python script. Expression matrices stored in these formats must first be exported from R to CSV or another supported format, preserving gene identifiers as the first column and sample identifiers as column names.

## 4. Running CIRCUST

Place your expression matrix in the `data/` directory.

For example:

```powershell
.\.venv\Scripts\python.exe .\scripts\run_circust.py --data data/matrixIn.csv --r-compatible --n-reps 5 --jobs -1 -o resultados_circust
```

The `--r-compatible` option enables the FMM fitting procedure designed to reproduce the methodological choices of the original R implementation.

The `--n-reps 5` option specifies five repetitions, while `--jobs -1` enables parallel execution using the available computational resources.

All calculations are performed from the input expression matrix, without relying on previously computed results.

## 5. Output Files

With K = 5, **a successfully completed execution generates 20 output files**:

- `sample_order_k1.csv` to `sample_order_k5.csv`: sample identifiers, positions, and estimated circular phases.
- `ordered_matrix_k1.csv` to `ordered_matrix_k5.csv`: observed expression values of the preprocessed genes, reordered according to the estimated sample ordering (CSV or Parquet, depending on the selected output format).
- `fitted_top_k1.csv` to `fitted_top_k5.csv`: FMM-fitted expression values for the TOP genes.
- `parameters_by_k.csv`: fitted parameters for each TOP gene and each repetition.
- `circust_atlas.csv`: summary based on the medians of the K estimates.
- `top_genes.csv`: selected TOP genes and preliminary statistics.
- `outliers.csv`: excluded samples.
- `manifest.json`: execution settings, selection information, and computational timings.

To inspect the execution manifest:

```powershell
Get-Content .\resultados_circust\manifest.json
```

## 6. Computational Performance

In a previous execution using the MatrixIn dataset (Windows, Python 3.14.4), CIRCUST processed **56,200 genes and 479 samples**.

After preprocessing, 21,149 genes were retained, 22 samples were excluded, and 65 TOP genes were selected. Five sample orderings were obtained.

The reported execution time was **48.4 seconds**, including output generation.

This result corresponds to a single execution and should be considered an illustrative performance reference rather than a formal benchmark against the R implementation.

## 7. GitHub Repository and Data Management

The Python implementation is maintained in the `python/` directory of the main CIRCUST repository.

Source code, configuration files, and documentation should be maintained within this directory, while the original R implementation remains in the repository root.

Do not upload virtual environments (`.venv/`), generated results, temporary files, or third-party materials without appropriate authorization.

Expression datasets should only be distributed when their licenses and applicable data-sharing conditions permit redistribution.
