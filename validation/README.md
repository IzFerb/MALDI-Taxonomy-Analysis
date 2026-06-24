# Validation workflow: R vs MATLAB

This folder holds side-by-side exports from both pipelines for comparison.

## Step 1 — Prepare sample data

Extract `Gums files.zip` (if not already done):

```r
source("r/setup_sample_data.R")
```

Expected layout in `sample_data/`:

| Group ID | Folder | Role |
|----------|--------|------|
| 1 | Blanks | Blank spectra (first 4 rows for PQN) |
| 2 | Historical samples | Historical Egyptian samples |
| 3 | Pigmented references | Pigmented reference gums |
| 4 | Reference samples | Acacia reference library |

## Step 2 — Run the R pipeline

```bash
Rscript r/validate_pipeline.R
```

This writes reference outputs to `validation/r/`:

- `similarity_matrix.csv` — Acacia Pei1 similarity matrix (primary comparison)
- `acacia_feature_matrix.csv` — peak feature matrix after blank removal
- `sample_manifest.csv` — file order and group mapping
- `summary.txt` — dimensions and basic stats

## Step 3 — Run MATLAB export

In MATLAB (requires the **icoshift** toolbox):

```matlab
cd('/path/to/art_institute_code')
Analysis_export
```

This writes matching files to `validation/matlab/`.

## Step 4 — Compare

If MATLAB exports exist:

```bash
Rscript r/validate_pipeline.R --compare-only
```

The report is written to `validation/r/comparison_report.txt`.

### What to look for

| Metric | Good match |
|--------|------------|
| Matrix dimensions | Identical (55×55 for Acacia after exclusions) |
| Max absolute difference | < 0.0001 (exact match unlikely due to icoshift/peak differences) |
| Pearson correlation | > 0.99 |
| PCA plots | Visually similar group separation |

Small differences are expected where R and MATLAB use different peak-finding or icoshift implementations. A correlation above 0.99 with similar PCA structure indicates a successful port.

## Troubleshooting

**Dimension mismatch** — Compare `validation/r/sample_manifest.csv` and `validation/matlab/sample_manifest.csv`. File order must match (alphabetical folders and files).

**Large differences** — Check that both pipelines use `sample_data/` as the data source and the same group ID mapping (see `r/config/config_sample_data.yaml`).
