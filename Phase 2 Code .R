# ==============================================================================
# Phase 2: Systematic Variable Audit Script
# ==============================================================================

library(tidyverse)

# ------------------------------------------------------------------------------
# Step 1: Isolate Data Tables
# ------------------------------------------------------------------------------
patient_df  <- biotab_tables$clinical_patient_luad
drug_df     <- biotab_tables$clinical_drug_luad
rad_df      <- biotab_tables$clinical_radiation_luad
follow_df   <- biotab_tables$clinical_follow_up_v1.0_luad

# ------------------------------------------------------------------------------
# Step 2: Search Candidate Column Names Across Indexed & Biotab Datasets
# ------------------------------------------------------------------------------
search_terms <- c("age", "gender", "sex", "race", "ethnicity", "smoke", "tobacco", "pack",
                  "vital", "status", "death", "follow", "stage", "pathologic", "grade", 
                  "residual", "margin", "therapy", "treatment", "radiation", "drug")

cat("=== Candidate Columns in Indexed Data (`luad_indexed`) ===\n")
indexed_matches <- grep(paste(search_terms, collapse = "|"), colnames(luad_indexed), ignore.case = TRUE, value = TRUE)
print(indexed_matches)

cat("\n=== Candidate Columns in Patient Biotab (`clinical_patient_luad`) ===\n")
biotab_matches <- grep(paste(search_terms, collapse = "|"), colnames(patient_df), ignore.case = TRUE, value = TRUE)
print(biotab_matches)


# ------------------------------------------------------------------------------
# Step 3: Inspect Key Survival & Exposure Fields in `luad_indexed`
# ------------------------------------------------------------------------------
cat("\n=== Summary of Key Candidate Fields in Indexed Data ===\n")
indexed_key_vars <- c("patient_id", "age_at_index", "gender", "race", "ethnicity", 
                      "vital_status", "days_to_death", "days_to_last_follow_up", 
                      "paper_tobacco_smoking_history", "ajcc_pathologic_stage",
                      "ajcc_pathologic_t", "ajcc_pathologic_n", "ajcc_pathologic_m",
                      "prior_malignancy", "synchronous_malignancy")

indexed_present <- intersect(indexed_key_vars, colnames(luad_indexed))

for (var in indexed_present) {
  cat("\n--- Indexed Field:", var, "---\n")
  print(table(luad_indexed[[var]], useNA = "ifany") %>% head(10))
}


# ------------------------------------------------------------------------------
# Step 4: Inspect Key Survival, Clinical & Staging Fields in Biotab Patient Data
# ------------------------------------------------------------------------------
# Note: Biotab tables often contain CDE descriptions in rows 1 and 2. We inspect row 3 onwards for data values.
cat("\n=== Value Distributions of Candidate Fields in Biotab Patient Table ===\n")
biotab_key_vars <- c("bcr_patient_barcode", "age_at_initial_pathologic_diagnosis", 
                     "gender", "race", "ethnicity", "tobacco_smoking_history", 
                     "number_pack_years_smoked", "vital_status", "days_to_death", 
                     "days_to_last_known_alive", "days_to_last_followup", 
                     "ajcc_pathologic_tumor_stage", "ajcc_pathologic_t", 
                     "ajcc_pathologic_n", "ajcc_pathologic_m", 
                     "histologic_grade", "residual_tumor")

biotab_present <- intersect(biotab_key_vars, colnames(patient_df))

for (var in biotab_present) {
  cat("\n--- Biotab Field:", var, "---\n")
  # Filter out header description lines if present
  vals <- patient_df[[var]][-c(1:2)]
  print(table(vals, useNA = "ifany") %>% head(12))
}


# ------------------------------------------------------------------------------
# Step 5: Audit Treatment Records (Drug & Radiation Tables)
# ------------------------------------------------------------------------------
cat("\n=== Treatment Data Availability Summary ===\n")
cat("Drug Table Dimensions:", dim(drug_df), "\n")
cat("Radiation Table Dimensions:", dim(rad_df), "\n")

if ("pharmaceutical_therapy_drug_name" %in% colnames(drug_df)) {
  cat("\nTop Reported Drugs in Drug Table:\n")
  print(table(drug_df$pharmaceutical_therapy_drug_name[-c(1:2)], useNA = "ifany") %>% head(10))
}

if ("radiation_therapy_type" %in% colnames(rad_df)) {
  cat("\nRadiation Therapy Types:\n")
  print(table(rad_df$radiation_therapy_type[-c(1:2)], useNA = "ifany"))
}