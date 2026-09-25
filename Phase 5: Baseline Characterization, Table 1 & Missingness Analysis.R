# ==============================================================================
# Phase 5: Baseline Characterization, Table 1 & Missingness Analysis
# ==============================================================================

# Install required packages if not present
if (!requireNamespace("gtsummary", quietly = TRUE)) install.packages("gtsummary")
if (!requireNamespace("patchwork", quietly = TRUE)) install.packages("patchwork")

library(tidyverse)
library(gtsummary)
library(patchwork)

# Ensure analytic_cohort from Phase 4 is loaded in memory
if (!exists("analytic_cohort")) {
  stop("Error: 'analytic_cohort' not found. Please run the Phase 4 script first.")
}

# ------------------------------------------------------------------------------
# Step 1: Missingness Audit Table
# ------------------------------------------------------------------------------
missing_summary <- analytic_cohort %>%
  select(
    `Age at Diagnosis` = age,
    `Sex` = gender,
    `Smoking Status` = smoking_status,
    `Pack-Years` = pack_years,
    `AJCC Stage` = stage_grouped,
    `Residual Margin` = residual_margin,
    `Vital Status` = vital_status_clean,
    `Overall Survival (Months)` = os_months
  ) %>%
  summarise(across(everything(), list(
    Missing_N = ~sum(is.na(.)),
    Missing_Pct = ~round(mean(is.na(.)) * 100, 2)
  ))) %>%
  pivot_longer(
    cols = everything(),
    names_to = c("Variable", ".value"),
    names_pattern = "(.*)_(Missing_N|Missing_Pct)"
  ) %>%
  arrange(desc(Missing_Pct))

cat("=== Comprehensive Variable Missingness Audit ===\n")
print(missing_summary)

# ------------------------------------------------------------------------------
# Step 2: Generate Publication-Quality Table 1 (Stratified by AJCC Stage)
# ------------------------------------------------------------------------------
table1 <- analytic_cohort %>%
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
    all_categorical() ~ "fisher.test"
  )) %>%
  bold_labels() %>%
  modify_header(label ~ "**Clinical Characteristic**")

cat("\n=== Table 1 Summary Text Representation ===\n")
print(table1)

# ------------------------------------------------------------------------------
# Step 3: Diagnostic Distribution Plots (Continuous Variables)
# ------------------------------------------------------------------------------
p1 <- ggplot(analytic_cohort, aes(x = age)) +
  geom_histogram(fill = "#2b5c8f", color = "white", bins = 25, alpha = 0.8) +
  geom_density(aes(y = after_stat(count) * 2.5), color = "#d95f02", linewidth = 1) +
  theme_minimal() +
  labs(
    title = "Age Distribution at Diagnosis",
    x = "Age (Years)",
    y = "Patient Count"
  )

p2 <- ggplot(analytic_cohort, aes(x = pack_years)) +
  geom_histogram(fill = "#41b6c4", color = "white", bins = 25, alpha = 0.8) +
  geom_density(aes(y = after_stat(count) * 5), color = "#d95f02", linewidth = 1) +
  theme_minimal() +
  labs(
    title = "Smoking Pack-Years Distribution",
    x = "Pack-Years",
    y = "Patient Count"
  )

p3 <- ggplot(analytic_cohort, aes(x = os_months, fill = factor(vital_status_clean))) +
  geom_histogram(position = "identity", alpha = 0.6, bins = 30, color = "white") +
  scale_fill_manual(
    values = c("0" = "#2ca25f", "1" = "#de2d26"),
    labels = c("Censored (Alive)", "Event (Deceased)"),
    name = "Status"
  ) +
  theme_minimal() +
  labs(
    title = "Overall Survival Follow-up Time Distribution",
    x = "Follow-up (Months)",
    y = "Patient Count"
  )

combined_plots <- (p1 | p2) / p3
print(combined_plots)