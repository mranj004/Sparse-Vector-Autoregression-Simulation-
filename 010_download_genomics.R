# 010_download_genomics

source("00_requirements.R")
data_dir <- "data"
dir.create(data_dir, showWarnings = FALSE)

# Download GSE130811 processed expression data
GEOquery::getGEOSuppFiles(
  "GSE130811",
  baseDir = data_dir,
  filter_regex = "expr_mRNA_CE10_coding"
)