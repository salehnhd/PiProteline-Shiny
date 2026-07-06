# PiProteline Shiny App

A graphical interface for the **PiProteline** proteomics pipeline:
preprocessing, quantitative analysis, protein–protein interaction (PPI) network
analysis, and functional enrichment — without writing any R code.

> **On the CoPPIs dependency.** This app builds the human interactome using
> `CoPPIs::interactome.hs` and `CoPPIs::filter_interactome()`. CoPPIs is a
> published R package (see *Citation* below); it provides a fixed, citable
> interactome snapshot, which keeps the network step reproducible. The score
> thresholds used here are `experimental = 150` and `database = 300`. CoPPIs is
> required and is listed in the install steps below.

---

## 1. Requirements

* **R ≥ 4.3** (RStudio recommended).
* The PiProteline backend and its dependencies.
* The CoPPIs package (interactome).

### Install the packages

Some dependencies live on **Bioconductor**, so install BiocManager first.

```r
# Bioconductor dependencies used by PiProteline
install.packages("BiocManager")
BiocManager::install(c("ggtree", "ggtreeExtra"))

# PiProteline backend and CoPPIs interactome
install.packages("remotes")
remotes::install_github("lomi95/PiProteline")
remotes::install_github("lomi95/CoPPIs")

# Shiny app dependencies
install.packages(c(
  "shiny", "shinyWidgets", "DT", "openxlsx",
  "igraph", "dplyr", "ggplot2"
))
```

> If `install_github` fails on `ggtree`/`ggtreeExtra`, install those from
> Bioconductor first (command above), then re-run the PiProteline install.

---

## 2. Launch the app

Put all the app files in one folder, then either:

* open `app.R` in RStudio and click **Run App**, or
* run from the R console:

```r
shiny::runApp("path/to/PiProteline-Shiny")
```

Files in this repository:

```
app.R
helpers.R
preprocessing_tab.R
quantitative_analysis_tab.R
network_analysis_tab.R
functional_analysis_tab.R
```

---

## 3. Input data format

Upload an `.xlsx`, `.xls`, or `.csv` file shaped like this:

| GeneName | Control_1 | Control_2 | Treated_1 | Treated_2 |
|----------|-----------|-----------|-----------|-----------|
| TP53     | 1203      | 1185      | 540       | 533       |
| ...      | ...       | ...       | ...       | ...       |

Rules:

* **Excel files with several sheets:** the app shows a **Worksheet** selector and
  defaults to the sheet with the most columns (your data), so a "Legend" sheet
  will not be loaded by mistake. Pick the correct data sheet if the guess is wrong.
* **One identifier column** (gene or protein). You pick which column this is in
  the app ("Select gene column").
* **Intensity / abundance columns**, one per sample.
* **Groups are read from the sample column names.** You type the group tags
  (e.g. `Control,Treated`) and the app finds every column whose name contains
  that tag. So `Control_1`, `Control_2`, ... all belong to group `Control`.
  ⚠️ Because matching is by substring, avoid group names where one is
  contained inside another (e.g. `A` and `AB`, or `WT` and `WT2`) — a column
  like `AB_1` would then match *both* `A` and `AB`. Use tags that are clearly
  distinct, e.g. `Ctrl` / `Treat` instead of `A` / `AB`.

### Reserved words (important)

Because groups are matched by text inside column names, a few words are
**not allowed** anywhere in your group names or sample column names — not even
as part of a word. The backend will stop with an error if it sees them:

```
fc, specific, centrality, weighted, number_of_genes, p_value, fdr, input_genes
```

For example, a group called `FCa` fails because it contains `fc`. Rename such
columns before uploading.

---

## 4. Workflow (run the tabs in order)

1. **Preprocessing**
   * Upload your file.
   * Select the gene column.
   * Type the group names (comma-separated).
   * Choose a normalization. *Tip:* `TotSigNorm` keeps values positive;
     `Znorm`, `Robust`, and `MinMax` can produce negative values that break the
     log2 fold change used later.
   * Click **Preprocess data**. (Optional: **Descriptive statistics**.)
2. **Quantitative Analysis**
   * Click **Run Quantitative Analysis**.
   * View the MANOVA table, and pick any pairwise **Volcano** or **MDS**
     contrast from the dropdowns.
3. **Network Analysis**
   * Set the critical-node quantile and a save prefix.
   * Click **Run network + functional analysis** (this reuses the preprocessed
     data; it does not redo preprocessing).
   * Browse centralities, critical nodes, and PPI graphs per group.
4. **Functional Analysis**
   * Explore enrichment tables and plots produced in step 3.

---

## 5. Quick test

To confirm your install works, upload any dataset in the format shown in
section 3 (one identifier column plus sample columns whose names contain the
group tags) and click through the four tabs in order. If every tab produces
output, the app is working.

---

## Notes

* The interactome here is **human** (`tax_ID = 9606`), matching
  `CoPPIs::interactome.hs`. Other species are not supported by this snapshot.
* Network analysis can take several minutes on large datasets.
* Results can be saved to disk via the **Save results as (prefix)** field in the
  Network tab.

---


