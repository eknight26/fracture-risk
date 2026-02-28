
# FRACTURE RISK PREDICTION - NHANES 2017-2020
# Script 02: T-Score Calculation
# Author: Ernest Caballero
# Date: 2024
# Description: Calculates femur neck and Ward's Triangle T-scores for the
# 2017-2020 cohort (50+) using a young adult reference group
# (NHANES 2010-2011, age 20-30), stratified by gender and race.

# WHO defines T-score as: (patient BMD - young adult mean BMD) / young adult SD
# Interpretation:
#    T-score >= -1.0            = Normal bone density
#    -1.0 > T-score > -2.5      = Osteopenia (low bone mass)
#    T-score   <=    -2.5       = Osteoporosis
# Reference: https://pmc.ncbi.nlm.nih.gov/articles/PMC2827823/

rm(list = ls())
options(scipen = 999)

library(tidyverse)
library(haven)


# -----------------------------------------------------------
# SECTION 1: YOUNG ADULT REFERENCE GROUP (NHANES 2010-2011)
# -----------------------------------------------------------
# Reference group: age 20-30, stratified by gender (RIAGENDR) and race (RIDRETH1)
# This mirrors standard clinical DXA T-score methodology

demo_ref <- read_xpt("data/raw/DEMO_F.XPT") %>%
  filter(RIDAGEYR >= 20 & RIDAGEYR <= 30) %>%
  select(SEQN, RIAGENDR, RIDAGEYR, RIDRETH1)

dxa_femur_ref <- read_xpt("data/raw/DXXFEM_F.XPT") %>%
  select(SEQN, DXXNKBMD, DXXWDBMD)

# Compute reference means and SDs per gender x race stratum
# These are used as the "young normal" baseline for T-score calculation
reference_stats <- demo_ref %>%
  inner_join(dxa_femur_ref, by = "SEQN") %>%
  group_by(RIAGENDR, RIDRETH1) %>%
  summarise(
    mean_ref_nk = mean(DXXNKBMD, na.rm = TRUE),
    sd_ref_nk   = sd(DXXNKBMD,   na.rm = TRUE),
    mean_ref_wd = mean(DXXWDBMD,  na.rm = TRUE),
    sd_ref_wd   = sd(DXXWDBMD,   na.rm = TRUE),
    n_ref       = n(),           # useful to inspect reference group sizes
    .groups = "drop"
  )

cat("Reference group strata:\n")
print(reference_stats)


# =============================================================================
# SECTION 2: PREPARE STUDY COHORT (NHANES 2017-2020, age 50+)
# =============================================================================

demo_study <- read_xpt("data/raw/P_DEMO.XPT") %>%
  filter(RIDAGEYR >= 50) %>%
  select(SEQN, RIAGENDR, RIDAGEYR, RIDRETH1)

dxa_femur_study <- read_xpt("data/raw/P_DXXFEM.XPT") %>%
  filter(DXAFMRST == 1) %>%          # valid scans only
  select(SEQN, DXXNKBMD, DXXWDBMD)


# =============================================================================
# SECTION 3: CALCULATE T-SCORES
# =============================================================================
# Join study cohort with reference stats by matching gender x race stratum,
# then apply WHO T-score formula

t_scores <- demo_study %>%
  inner_join(dxa_femur_study, by = "SEQN") %>%
  left_join(reference_stats, by = c("RIAGENDR", "RIDRETH1")) %>%
  mutate(
    # WHO T-score formula: (individual BMD - young adult mean) / young adult SD
    T_score_nk = (DXXNKBMD - mean_ref_nk) / sd_ref_nk,   # femur neck
    T_score_wd = (DXXWDBMD - mean_ref_wd) / sd_ref_wd,   # Ward's Triangle

    # Clinical classification derived from T-scores (WHO 1994)
    osteo_class_nk = case_when(
      T_score_nk >= -1.0              ~ "Normal",
      T_score_nk > -2.5               ~ "Osteopenia",
      T_score_nk <= -2.5              ~ "Osteoporosis",
      TRUE                            ~ NA_character_
    )
  ) %>%
  select(SEQN, T_score_nk, T_score_wd, osteo_class_nk)


# =============================================================================
# SECTION 4: QUALITY CHECK & SAVE
# =============================================================================

cat("\nT-score summary:\n")
summary(t_scores[, c("T_score_nk", "T_score_wd")])

cat("\nClinical classification (femur neck):\n")
print(table(t_scores$osteo_class_nk, useNA = "ifany"))

# Save for use in 01_data_prep.R
write_csv(t_scores, "data/processed/dexa_t_scores.csv")
cat("Saved: data/processed/dexa_t_scores.csv\n")
