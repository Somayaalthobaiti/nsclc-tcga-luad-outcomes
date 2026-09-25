# ==============================================================================
# Phase 5 Table 1 Fix & Phase 6: Survival Analysis & Diagnostic Pipeline
# ==============================================================================

if (!requireNamespace("survival", quietly = TRUE)) install.packages("survival")
if (!requireNamespace("survminer", quietly = TRUE)) install.packages("survminer")

library(tidyverse)
library(survival)
library(survminer)
library(gtsummary)

# Ensure analytic_cohort is present
if (!exists("analytic_cohort")) {
  stop("Error: 'analytic_cohort' object not found. Please run Phase 4 script.")
}

# ------------------------------------------------------------------------------
# Step 1: Fixed Publication Table 1 (Chi-squared Test Fix)
# ------------------------------------------------------------------------------
table1_clean <- analytic_cohort %>%
  select(
    age,
    gender,
    smoking_status,
    pack_years,
    stage_grouped,
    residual_margin,
    vital_status_clean,
    os_months
  ) %>%
  tbl_summary(
    by = stage_grouped,
    label = list(
      age ~ "Age at Diagnosis (Years)",
      gender ~ "Sex",
      smoking_status ~ "Smoking Status",
      pack_years ~ "Smoking Pack-Years",
      residual_margin ~ "Surgical Residual Margin",
      vital_status_clean ~ "Vital Status (Deceased)",
      os_months ~ "Follow-up Duration (Months)"
    ),
    statistic = list(
      all_continuous() ~ "{median} ({p25}, {p75})",
      age ~ "{mean} ({sd})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    digits = all_continuous() ~ 1,
    missing = "ifany",
    missing_text = "Missing / Unrecorded"
  ) %>%
  add_overall(last = FALSE) %>%
  add_p(test = list(
    all_continuous() ~ "kruskal.test",
    all_categorical() ~ "chisq.test"
  )) %>%
  bold_labels()

cat("=== Fixed Table 1: Baseline Characteristics Stratified by Stage ===\n")
print(table1_clean)


# ------------------------------------------------------------------------------
# Step 2: Kaplan-Meier Survival Analysis (Overall & Stratified)
# ------------------------------------------------------------------------------
# Construct primary survival object (Time in Months, Event = Deceased)
surv_obj <- Surv(time = analytic_cohort$os_months, event = analytic_cohort$vital_status_clean)

# Overall Cohort KM Fit
km_overall <- survfit(surv_obj ~ 1, data = analytic_cohort)

cat("\n=== Overall Cohort Survival Summary ===\n")
print(km_overall)

# KM Fit Stratified by AJCC Stage
km_stage <- survfit(surv_obj ~ stage_grouped, data = analytic_cohort)

cat("\n=== Kaplan-Meier Median Survival by AJCC Stage ===\n")
print(km_stage)

# Log-Rank Test across Stage Strata
logrank_stage <- survdiff(surv_obj ~ stage_grouped, data = analytic_cohort)
cat("\n=== Log-Rank Test (AJCC Stage) ===\n")
print(logrank_stage)

# KM Fit Stratified by Residual Margin Status
km_margin <- survfit(surv_obj ~ residual_margin, data = analytic_cohort)
logrank_margin <- survdiff(surv_obj ~ residual_margin, data = analytic_cohort)

cat("\n=== Log-Rank Test (Surgical Margins) ===\n")
print(logrank_margin)


# ------------------------------------------------------------------------------
# Step 3: Univariable & Multivariable Cox Regression Modeling
# ------------------------------------------------------------------------------
# Univariable Cox Models
uv_cox_table <- analytic_cohort %>%
  select(
    os_months,
    vital_status_clean,
    age,
    gender,
    smoking_status,
    stage_grouped,
    residual_margin
  ) %>%
  tbl_uvregression(
    method = coxph,
    y = Surv(os_months, vital_status_clean),
    exponentiate = TRUE,
    label = list(
      age ~ "Age (per year)",
      gender ~ "Sex",
      smoking_status ~ "Smoking Status",
      stage_grouped ~ "AJCC Pathologic Stage",
      residual_margin ~ "Surgical Margin Status"
    )
  ) %>%
  bold_labels() %>%
  bold_p(t = 0.05)

cat("\n=== Univariable Cox Proportional Hazards Models ===\n")
print(uv_cox_table)

# Multivariable Cox Model
mv_cox_fit <- coxph(
  Surv(os_months, vital_status_clean) ~ age + gender + smoking_status + stage_grouped + residual_margin,
  data = analytic_cohort
)

mv_cox_table <- tbl_regression(
  mv_cox_fit,
  exponentiate = TRUE,
  label = list(
    age ~ "Age (per year)",
    gender ~ "Sex",
    smoking_status ~ "Smoking Status",
    stage_grouped ~ "AJCC Pathologic Stage",
    residual_margin ~ "Surgical Margin Status"
  )
) %>%
  bold_labels() %>%
  bold_p(t = 0.05)

cat("\n=== Multivariable Adjusted Cox Proportional Hazards Model ===\n")
print(mv_cox_table)


# ------------------------------------------------------------------------------
# Step 4: Proportional Hazards Assumption Diagnostics (Schoenfeld Residuals)
# ------------------------------------------------------------------------------
ph_diagnostics <- cox.zph(mv_cox_fit)

cat("\n=== Proportional Hazards Test (Scaled Schoenfeld Residuals) ===\n")
print(ph_diagnostics)