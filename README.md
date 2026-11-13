
# Sparse Vector Autoregression Simulation

## Project Overview

This project investigates how Gaussian innovation covariance structures affect estimation and support recovery in a sparse vector autoregressive model of order one (VAR(1)):

$$
Y_t = A Y_{t-1} + \varepsilon_t, \qquad
\varepsilon_t \overset{\mathrm{i.i.d.}}{\sim} N_d(0, \Sigma_\varepsilon).
$$

The simulation compares a fixed-penalty LASSO estimator (`glmnet`, with $\lambda = 0.05$) and a time-slice-tuned LASSO estimator (`sparsevar`). The innovation covariance structures are **Block-I**, **Block-II**, and **Toeplitz**. Accuracy is assessed using relative Frobenius error, mean squared error, precision, recall, and F1-score.

## Repository Structure

| File | Description |
| --- | --- |
| `00_requirements.R` | Loads R packages used in the project |
| `010_download_genomics.R` | Downloads the optional GSE130811 gene-expression dataset |
| `01_generate_A.R` | Generates a stable, sparse transition matrix $A$ |
| `02_generate_errors.R` | Generates Gaussian innovations with Block-I, Block-II, or Toeplitz covariance |
| `03_generate_VAR.R` | Simulates VAR(1) observations using $A$ and the innovations |
| `04_estimate_A.R` | Estimates $\hat A$ using fixed and time-slice-tuned LASSO |
| `05_evaluate_A.R` | Calculates coefficient-estimation and support-recovery metrics |
| `06_run_simulation.R` | Runs the simulation and saves results |
| `07_SimAnalysis.Rmd` | Analyzes and visualizes the simulation results |
| `08_DataCleaning.Rmd` | Explores and preprocesses the optional gene-expression dataset |
| `FinalReport.Rmd` | Source for the written report |
| `FinalReport.pdf` | Rendered report |
| `results/` | Saved simulation results and coefficient matrix |
| `data/` | Raw and processed gene-expression data |

## Reproducing the Reported Simulation

Run the commands below **from the root of the cloned repository** in R or RStudio.

### 1. Install required packages

On a fresh R installation, kindly install the CRAN dependencies:

```r
install.packages(c("glmnet", "sparsevar", "dplyr", "ggplot2", "tidyr", "rmarkdown"))
```

> **Note:** The current `00_requirements.R` also checks for and loads `GEOquery`, and attempts to install it through Bioconductor if absent. This is needed for the optional genomics work  (not for the VAR simulation itself).

### 2. Run the experiment used in the report

```r
source("00_requirements.R")
source("06_run_simulation.R")

sim <- runSim(
  T = 300,
  d = 500,
  spar = 1 - (5 / 499),
  rho_A = 0.8,
  lambda = 0.05,
  burn_in = 200,
  n_rep = 1,
  rhoVal = 0.7,
  seedA = 293,
  seedSim = 123,
  saveResults = TRUE,
  keepEstimates = TRUE
)

print(sim$summary)
```

### 3. Inspect saved outputs

When `saveResults = TRUE`, the simulation writes the following files to `results/`:

- `simulation_results.csv` — metrics for each estimator and covariance structure
- `simulation_summary.csv` — condition-level summaries
- `paired_comparison.csv` — paired differences between the estimators
- `true_A.csv` — the true transition matrix

To run the analysis notebook, open `07_SimAnalysis.Rmd` in RStudio and knit it. Note that this notebook runs the 500-variable simulation again; its saved results files will be overwritten when the simulation completes.

## Optional Genomic Data Exploration

`08_DataCleaning.Rmd` processes the GSE130811 *C. elegans* gene-expression time series. Its download step is implemented in `010_download_genomics.R`. This analysis is separate from the simulated VAR(1) coefficient-recovery experiment.

## Reference

Basu, S., & Michailidis, G. (2015). Regularized estimation in sparse high-dimensional time series models. *The Annals of Statistics, 43*(4), 1535–1567. https://doi.org/10.1214/15-AOS1315
