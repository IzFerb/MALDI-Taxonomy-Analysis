# MALDI Gum Taxonomy Analysis — R Guide

This folder contains an R port of the original MATLAB `Analysis.m` pipeline from the Art Institute of Chicago gum taxonomy study. It is designed so researchers **without coding experience** can run the full workflow.

## Quick start (recommended)

### Option A — Point-and-click app (easiest)

1. Install [R](https://cran.r-project.org/) and optionally [RStudio](https://posit.co/download/rstudio-desktop/).
2. Open this project folder in RStudio.
3. Run once to install dependencies:
   ```r
   source("r/install_packages.R")
   ```
4. Launch the app:
   ```r
   source("r/launch_app.R")
   ```
5. In the app, confirm the **Project folder** path, then click **Run analysis**.
6. View PCA plots in the **PCA plots** tab. Results are saved to `r/output/`.

### Option B — One-click script

From Terminal (in the project root):

```bash
Rscript r/run_analysis.R
```

Or in R:

```r
source("r/run_analysis.R")
```

### Option C — Edit settings without code

Open `r/config/config.yaml` in any text editor. Change values like m/z range, alignment settings, or group IDs, then run Option A or B.

---

## Preparing your data

The pipeline expects the same layout as the MATLAB script:

```
project_folder/
├── Blank 1.txt          # blank spectra (project root)
├── Blank 2.txt
├── Blank 3.txt
├── Blank 4.txt
├── Acacia/              # one subfolder per sample group
│   ├── sample1.txt
│   └── sample2.txt
├── Historical/
│   └── ...
└── Pigmented/
    └── ...
```

Each `.txt` file has two tab-separated columns: **m/z** and **intensity**.

Folder order defines numeric group IDs (1st folder = group 1, 2nd = group 2, …). Match the original study by placing folders in the same order, or update `reference_groups`, `acacia_group_id`, etc. in `config.yaml`.

---

## What the pipeline does

| Step | MATLAB equivalent | R function |
|------|-------------------|------------|
| Load spectra | `dir()` loop | `load_spectra()` |
| Interpolate m/z grid | `interp1(..., 'spline')` | `interpolate_spectra()` |
| Align spectra | `icoshift('median', ..., 'whole', 34)` | `icoshift()` |
| Standardize | row-wise SNV | `snv_rows()` |
| Peak picking | `findpeaks` | `find_and_align_peaks()` |
| Blank removal | PQN step | `remove_blank_peaks()` |
| Similarity | `simbin(..., 'Pei1')` | `simbin()` |
| PCA plots | `plot_pca_metabolomics` | `plot_pca_metabolomics()` |

---

## Output files

After a successful run, check `r/output/`:

- `01_aligned_spectra.png` — aligned spectra overlay
- `02_snv_spectra.png` — standardized spectra
- `acacia_pca_*.png` — PCA score plots (PC pairs)
- `historical_pca.png` — historical samples projected onto Acacia model
- `pigmented_pca.png` — pigmented references projected onto Acacia model
- `acacia_similarity_matrix.csv` — pairwise Pei1 similarity
- `acacia_sample_classes.csv` — sample class labels

---

## Troubleshooting

**"No sample .txt files found"**  
Create subfolders with `.txt` spectra. Files in the project root alone (except blanks) are not treated as samples.

**Package errors**  
Run `source("r/install_packages.R")` again.

**Results differ slightly from MATLAB**  
Small numerical differences can come from interpolation, peak-finding, or icoshift implementation. The R icoshift uses the same median-reference, whole-spectrum alignment with a configurable max shift.

---

## Validating against MATLAB

Sample data is in `Gums files.zip`. Full instructions are in **[validation/README.md](../validation/README.md)**.

Quick version:

```bash
# 1. Extract sample data and run R pipeline (~2 min)
Rscript r/validate_pipeline.R

# 2. In MATLAB (requires icoshift toolbox):
#    Analysis_export

# 3. Compare results
Rscript r/validate_pipeline.R --compare-only
```

The primary check is the **Acacia similarity matrix** (`validation/r/similarity_matrix.csv` vs `validation/matlab/similarity_matrix.csv`). A correlation above 0.99 indicates a successful match.

---

## Support

Original MATLAB authors: G.H. Tinnevelt et al. — gerjen.tinnevelt@ru.nl

For the R port, open an issue or contact your local analysis support team.

## License

CC-BY 4.0 (same as original project)
