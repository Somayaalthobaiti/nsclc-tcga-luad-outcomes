# ==============================================================================
# Phase 4: Data Cleaning, Variable Harmonization & Cohort Derivation (Fixed Lengths)
# ==============================================================================

library(tidyverse)

# ------------------------------------------------------------------------------
# Step 1: Load and Clean Raw Patient Biotab & Indexed Data
# ------------------------------------------------------------------------------
patient_raw <- biotab_tables$clinical_patient_luad
indexed_raw <- luad_indexed

# Strip Biotab header description rows (rows 1 and 2)
patient_clean <- patient_raw[-c(1:2), ]

# Standardize TCGA missing value strings to NA
recode_missing <- function(x) {
  x_clean <- trimws(as.character(x))
  missing_strings <- c("[Not Available]", "[Not Evaluated]", "[Unknown]", 
                       "[Discrepancy]", "[Not Applicable]", "null", "")
  ifelse(x_clean %in% missing_strings, NA_character_, x_clean)
}

patient_clean <- patient_clean %>% mutate(across(everything(), recode_missing))

# ------------------------------------------------------------------------------
# Step 2: Merge Data Sources FIRST to Standardize Row Length (N = 585)
# ------------------------------------------------------------------------------
# Extract patient barcode identifier column dynamically from indexed dataset
id_col_idx <- grep("submitter_id|barcode|patient_id", colnames(indexed_raw), ignore.case = TRUE, value = TRUE)[1]

indexed_clean <- indexed_raw %>%
  mutate(bcr_patient_barcode = substr(.data[[id_col_idx]], 1, 12)) %>%
  distinct(bcr_patient_barcode, .keep_all = TRUE)

# Left-join Biotab data onto the complete indexed dataset
df_merged <- indexed_clean %>%
  left_join(patient_clean, by = "bcr_patient_barcode", suffix = c("_idx", "_bio"))

# ------------------------------------------------------------------------------
# Step 3: Extract and Coalesce Fields Directly from Merged Dataframe
# ------------------------------------------------------------------------------
# Helper function to find column in df_merged and return character vector of length nrow(df)
get_field <- function(df, pattern1, pattern2 = NULL) {
  cols1 <- grep(pattern1, colnames(df), ignore.case = TRUE, value = TRUE)
  v1 <- if (length(cols1) > 0) as.character(df[[cols1[1]]]) else rep(NA_character_, nrow(df))
  
  if (!is.null(pattern2)) {
    cols2 <- grep(pattern2, colnames(df), ignore.case = TRUE, value = TRUE)
    v2 <- if (length(cols2) > 0) as.character(df[[cols2[1]]]) else rep(NA_character_, nrow(df))
    return(coalesce(v1, v2))
  }
  return(v1)
}

# Extract coalesced character vectors (all guaranteed to have length = nrow(df_merged) = 585)
days_death_vec  <- get_field(df_merged, "days_to_death.*_bio", "days_to_death")
days_follow_vec <- get_field(df_merged, "days_to_last_follow.*_bio", "days_to_last_follow")
vital_stat_vec  <- get_field(df_merged, "vital_status.*_bio", "vital_status")
age_vec         <- get_field(df_merged, "age_at_initial_pathologic_diagnosis", "age")
gender_vec      <- get_field(df_merged, "gender.*_bio", "gender|sex")
stage_vec       <- get_field(df_merged, "ajcc_pathologic_tumor_stage", "ajcc_pathologic_stage|stage")
res_tumor_vec   <- get_field(df_merged, "residual_tumor")
tobacco_vec     <- get_field(df_merged, "tobacco_smoking_history")
pack_yrs_vec    <- get_field(df_merged, "pack")
prior_mal_vec   <- get_field(df_merged, "prior_malignancy")

# Assemble clean harmonized dataset
df_harmonized <- tibble(
  bcr_patient_barcode = df_merged$bcr_patient_barcode,
  days_to_death       = as.numeric(days_death_vec),
  days_to_followup    = as.numeric(days_follow_vec),
  vital_status_raw    = vital_stat_vec,
  age                 = as.numeric(age_vec),
  pack_years          = as.numeric(pack_yrs_vec),
  gender_raw          = gender_vec,
  stage_str           = stage_vec,
  residual_str        = res_tumor_vec,
  tobacco_str         = tobacco_vec,
  prior_malignancy    = prior_mal_vec
) %>%
  mutate(
    # Vital Status (1 = Dead, 0 = Alive)
    vital_status_clean = case_when(
      toupper(vital_status_raw) == "DEAD" ~ 1,
      toupper(vital_status_raw) == "ALIVE" ~ 0,
      TRUE ~ NA_real_
    ),
    
    # Gender
    gender = case_when(
      toupper(gender_raw) %in% c("FEMALE", "F") ~ "Female",
      toupper(gender_raw) %in% c("MALE", "M") ~ "Male",
      TRUE ~ NA_character_
    ),
    gender = factor(gender, levels = c("Female", "Male")),
    
    # Smoking History
    smoking_status = case_when(
      tobacco_str == "1" ~ "Lifelong Non-Smoker",
      tobacco_str %in% c("2", "3", "4", "5") ~ "Ever Smoker",
      TRUE ~ NA_character_
    ),
    smoking_status = factor(smoking_status, levels = c("Lifelong Non-Smoker", "Ever Smoker")),
    
    # Stage Grouping
    stage_grouped = case_when(
      grepl("Stage I[AB]?$", stage_str, ignore.case = TRUE) ~ "Stage I",
      grepl("Stage II[AB]?$", stage_str, ignore.case = TRUE) ~ "Stage II",
      grepl("Stage III[AB]?$", stage_str, ignore.case = TRUE) ~ "Stage III",
      grepl("Stage IV$", stage_str, ignore.case = TRUE) ~ "Stage IV",
      TRUE ~ NA_character_
    ),
    stage_grouped = factor(stage_grouped, levels = c("Stage I", "Stage II", "Stage III", "Stage IV")),
    
    # Surgical Residual Margin
    residual_margin = case_when(
      residual_str == "R0" ~ "R0 (No Residual)",
      residual_str %in% c("R1", "R2", "RX") ~ "R1/R2/RX (Residual/Unknown)",
      TRUE ~ "Unknown/Unrecorded"
    ),
    residual_margin = factor(residual_margin, levels = c("R0 (No Residual)", "R1/R2/RX (Residual/Unknown)", "Unknown/Unrecorded")),
    
    # Primary Overall Survival Time (Days)
    os_days = case_when(
      vital_status_clean == 1 ~ days_to_death,
      vital_status_clean == 0 ~ days_to_followup,
      TRUE ~ NA_real_
    ),
    
    # OS Time Transformations
    os_months = os_days / 30.4375,
    os_years  = os_days / 365.25
  )

# ------------------------------------------------------------------------------
# Step 4: Apply Inclusion / Exclusion Criteria & Build Attrition Table
# ------------------------------------------------------------------------------
n0 <- nrow(df_harmonized)

# Inclusion Step 1: Valid vital status and non-negative follow-up duration
step1_df <- df_harmonized %>%
  filter(!is.na(vital_status_clean), !is.na(os_days), os_days >= 0)
n1 <- nrow(step1_df)

# Inclusion Step 2: Exclude documented prior malignancy
step2_df <- step1_df %>%
  filter(is.na(prior_malignancy) | prior_malignancy != "yes")
n2 <- nrow(step2_df)

# Inclusion Step 3: Confirmed AJCC Pathologic Stage
step3_df <- step2_df %>%
  filter(!is.na(stage_grouped))
n3 <- nrow(step3_df)

# Final Analytic Cohort Assignment
analytic_cohort <- step3_df

# Construct Attrition Flow Table
attrition_table <- tibble(
  Selection_Step = c(
    "1. Initial TCGA-LUAD Cohort (Indexed GDC Manifest)",
    "2. Exclude missing vital status or invalid follow-up duration (< 0 days)",
    "3. Exclude documented prior malignancy",
    "4. Exclude missing/unclassified baseline AJCC stage"
  ),
  Patients_Remaining = c(n0, n1, n2, n3),
  Excluded = c(0, n0 - n1, n1 - n2, n2 - n3)
)

cat("=== Cohort Attrition Flow Table ===\n")
print(attrition_table)

cat("\n=== Final Analytic Cohort Metrics ===\n")
cat("Total Analytic Cohort Size (N):", nrow(analytic_cohort), "\n")
cat("Total Observed Deaths (Events):", sum(analytic_cohort$vital_status_clean == 1), "\n")
cat("Median Follow-up Duration (Months):", round(median(analytic_cohort$os_months, na.rm = TRUE), 2), "\n")