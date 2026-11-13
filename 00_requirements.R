## Packages to install for simulation!
# CRAN packages
cran_packages <- c(
  "glmnet",
  "sparsevar",
  "dplyr",
  "ggplot2",
  "tidyr",
  "rmarkdown"
)

missing_packages <- cran_packages[
  !vapply(cran_packages, requireNamespace,
          logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  install.packages(
    missing_packages,
    repos = "https://cloud.r-project.org"
  )
}

# Load packages
library(glmnet)
library(sparsevar)
library(dplyr)
library(ggplot2)
library(tidyr)