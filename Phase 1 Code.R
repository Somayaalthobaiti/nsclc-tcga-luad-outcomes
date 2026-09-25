# ==============================================================================
# Phase 1: TCGA-LUAD Data Discovery & Metadata Retrieval
# ==============================================================================

# Install required packages if not already present
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("TCGAbiolinks", quietly = TRUE)) BiocManager::install("TCGAbiolinks")
if (!requireNamespace("tidyverse", quietly = TRUE)) install.packages("tidyverse")

library(TCGAbiolinks)
library(tidyverse)

# ------------------------------------------------------------------------------
# Step 1: Fetch Indexed Clinical Metadata via GDC API
# ------------------------------------------------------------------------------
cat("--- Fetching GDC Indexed Clinical Data ---\n")
luad_indexed <- GDCquery_clinic(project = "TCGA-LUAD", type = "clinical")

# Display dimensions (rows = patients, columns = indexed variables)
cat("\nIndexed Clinical Data Dimensions (Rows x Columns):\n")
print(dim(luad_indexed))

# Display column names of the indexed clinical file
cat("\nVariables present in Indexed Clinical Data:\n")
print(colnames(luad_indexed))


# ------------------------------------------------------------------------------
# Step 2: Query and Download BCR Biotab Clinical Supplements
# ------------------------------------------------------------------------------
cat("\n--- Querying GDC BCR Biotab Files ---\n")
query_biotab <- GDCquery(
  project = "TCGA-LUAD",
  data.category = "Clinical",
  data.type = "Clinical Supplement",
  data.format = "BCR Biotab"
)

# Inspect query manifest
biotab_manifest <- getResults(query_biotab)
cat("\nBCR Biotab Manifest Summary:\n")
print(biotab_manifest[, c("file_name", "file_size", "data_type")])

# Download Biotab files to local directory
GDCdownload(query_biotab)

# Prepare/extract Biotab tables into a named list of data frames
biotab_tables <- GDCprepare(query_biotab)

cat("\nExtracted Biotab Table Names:\n")
print(names(biotab_tables))