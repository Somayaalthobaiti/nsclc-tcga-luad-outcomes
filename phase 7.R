# ==============================================================================
# Phase 7: Stratified Cox Modeling, Publication KM Curves & Sensitivity Analysis
# ==============================================================================

library(tidyverse)
library(survival)
library(survminer)
library(gtsummary)

if (!exists("analytic_cohort")) {
  stop("Error: 'analytic_cohort' object not found. Please run Phase 4 script first.")
}

# ------------------------------------------------------------------------------
# Step 1: Remediation — Stratified Cox Proportional Hazards Model
# ------------------------------------------------------------------------------
strat_cox_fit <- coxph(
  Surv(os_months, vital_status_clean) ~ age + gender + smoking_status + residual_margin + strata(stage_grouped),
  data = analytic_cohort
)

strat_cox_table <- tbl_regression(
  strat_cox_fit,
  exponentiate = TRUE,
  label = list(
    age ~ "Age (per year)",
    gender ~ "Sex",
    smoking_status ~ "Smoking Status",
    residual_margin ~ "Surgical Margin Status"
  )
) %>%
  bold_labels() %>%
  bold_p(t = 0.05)

cat("=== Stratified Multivariable Cox Model [strata(stage_grouped)] ===\n")
print(strat_cox_table)

# Post-Remediation Schoenfeld Diagnostics
ph_diag_strat <- cox.zph(strat_cox_fit)
cat("\n=== Post-Stratification Proportional Hazards Diagnostic ===\n")
print(ph_diag_strat)


# ------------------------------------------------------------------------------
# Step 2: Publication-Grade Kaplan-Meier Survival Curves
# ------------------------------------------------------------------------------
km_stage <- survfit(Surv(os_months, vital_status_clean) ~ stage_grouped, data = analytic_cohort)

# Figure 1: KM Curve by AJCC Pathologic Stage
km_plot_stage <- ggsurvplot(
  km_stage,
  data = analytic_cohort,
  pval = TRUE,
  pval.coord = c(120, 0.85),
  pval.size = 4.5,
  conf.int = FALSE,
  risk.table = TRUE,
  risk.table.col = "strata",
  risk.table.height = 0.25,
  risk.table.y.text = FALSE,
  ggtheme = theme_classic(base_size = 12),
  palette = c("#2b5c8f", "#41b6c4", "#e6ab02", "#d95f02"),
  legend.title = "AJCC Stage",
  legend.labs = c("Stage I", "Stage II", "Stage III", "Stage IV"),
  title = "Overall Survival by AJCC Pathologic Stage (TCGA-LUAD)",
  xlab = "Follow-up Duration (Months)",
  ylab = "Overall Survival Probability"
)

# Figure 2: KM Curve by Surgical Residual Margin
km_margin <- survfit(Surv(os_months, vital_status_clean) ~ residual_margin, data = analytic_cohort)

km_plot_margin <- ggsurvplot(
  km_margin,
  data = analytic_cohort,
  pval = TRUE,
  pval.coord = c(120, 0.85),
  pval.size = 4.5,
  conf.int = FALSE,
  risk.table = TRUE,
  risk.table.col = "strata",
  risk.table.height = 0.25,
  risk.table.y.text = FALSE,
  ggtheme = theme_classic(base_size = 12),
  palette = c("#1b9e77", "#d95f02", "#7570b3"),
  legend.title = "Surgical Margin",
  title = "Overall Survival by Surgical Residual Margin Status",
  xlab = "Follow-up Duration (Months)",
  ylab = "Overall Survival Probability"
)

print(km_plot_stage)
print(km_plot_margin)


# ------------------------------------------------------------------------------
# Step 3: Sensitivity Analysis — Complete Case vs. Imputed/Subgroup Stability
# ------------------------------------------------------------------------------
# Sensitivity Model 1: Exclude Stage IV Patients (Early-Stage Resectable Subcohort)
early_stage_df <- analytic_cohort %>% filter(stage_grouped %in% c("Stage I", "Stage II", "Stage III"))

sens_fit_early <- coxph(
  Surv(os_months, vital_status_clean) ~ age + gender + smoking_status + residual_margin + strata(stage_grouped),
  data = early_stage_df
)

sens_table_early <- tbl_regression(
  sens_fit_early,
  exponentiate = TRUE,
  label = list(
    age ~ "Age (per year)",
    gender ~ "Sex",
    smoking_status ~ "Smoking Status",
    residual_margin ~ "Surgical Margin Status"
  )
) %>%
  bold_labels()

cat("\n=== Sensitivity Model: Resectable Cohort (Stages I–III, N = ", nrow(early_stage_df), ") ===\n", sep = "")
print(sens_table_early)