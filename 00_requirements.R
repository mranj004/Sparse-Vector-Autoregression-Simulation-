## Packages to install for simulation!
library(glmnet)
library(sparsevar)

## Packages for data analysis
library(dplyr)
library(ggplot2)

## Package for genomic data retrieval
if (!requireNamespace("GEOquery", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  
  BiocManager::install("GEOquery", ask = FALSE, update = FALSE)
}

library(GEOquery)