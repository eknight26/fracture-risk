# =============================================================================
# FRACTURE RISK PREDICTION - NHANES 2017-2020
# Script 03: Missing Data Imputation
# Author: Ernest Caballero
# Description: Imputes missing values using MICE (Multiple Imputation by
#              Chained Equations). A single completed dataset is produced
#              (m=1 imputation) suitable for ML modelling. Survey weights
#              and ID variables are excluded from imputation but retained.
#
# Method rationale:
#   MICE is preferred over simple median/mean imputation because it preserves
#   relationships between variables. For a clinical dataset with structured
#   missingness (DEXA scans only on eligible participants), this is important.
#   m=1 is used here for ML pipelines; for inferential statistics use m=5+.
# =============================================================================

rm(list = ls())
options(scipen = 999)

library(tidyverse)
library(mice)      # multiple imputation by chained equations


# =============================================================================
# SECTION 1: LOAD DATA
# =============================================================================

main_data <- read_csv("data/processed/main_data_with_tscores.csv",
                      show_col_types = FALSE)

cat("Dimensions on load:", dim(main_data), "\n")
cat("Total missing values:", sum(is.na(main_data)), "\n")


# =============================================================================
# SECTION 2: DEFINE FEATURE GROUPS
# =============================================================================
# Three clinically meaningful groups reflecting how fracture risk is assessed:
#   Group 1 — Demography & Lifestyle: who the patient is and how they live
#   Group 2 — Laboratory: what their blood reveals about bone metabolism
#   Group 3 — Bone Scans: what DEXA imaging and T-scores directly show

group_demo_lifestyle <- c(
  # Demographics
  "age", "gender", "ethnicity", "income_ratio",
  # Body measurements & vitals
  "bmi", "mean_systolic", "mean_diastolic",
  # Comorbidities
  "high_bp", "high_choles", "diabetes_status", "a1c_level",
  "arthritis", "chf", "chd", "heart_attack", "stroke",
  "liver", "thyroid", "copd", "overweight", "asthma", "cancer",
  # Osteoporosis-related history
  "had_osteoporosis", "other_fractures", "prednisone",
  "mother_fx_hip", "father_fx_hip",
  # Lifestyle
  "alc_drinks_per_day", "cigarettes_per_day", "sedentary_mins_day"
)

group_laboratory <- c(
  # Bone metabolism markers
  "alp", "phosphate", "calcium",
  # Nutritional marker
  "serum_folate",
  # Environmental exposures linked to bone loss
  "blood_lead", "blood_cadmium", "blood_mercury",
  "blood_selenium", "blood_manganese",
  # Glycaemic markers
  "ldl_level"
)

group_bone_scans <- c(
  # Femur BMD (site-specific)
  "bmd_femur_total", "bmd_femur_neck", "bmd_femur_troch",
  "bmd_femur_inter", "bmd_femur_ward",
  # Spine BMD (vertebral levels)
  "bmd_spine_total", "bmd_l1", "bmd_l2", "bmd_l3", "bmd_l4",
  # WHO T-scores (derived in script 02)
  "T_score_nk", "T_score_wd"
)

# Variables never imputed — identifiers, weights, target
exclude_from_imputation <- c(
  "SEQN", "SDMVPSU", "SDMVSTRA",
  "samp_wt_intrvw", "samp_wt_mec", "folate_wt",
  "hx_fracture"
)

# All features to impute
all_features <- c(group_demo_lifestyle, group_laboratory, group_bone_scans)

cat("\nFeature group sizes:\n")
cat("  Demo & Lifestyle:", length(group_demo_lifestyle), "features\n")
cat("  Laboratory:      ", length(group_laboratory), "features\n")
cat("  Bone Scans:      ", length(group_bone_scans), "features\n")
cat("  Total:           ", length(all_features), "features\n")


# =============================================================================
# SECTION 3: PREPARE DATA FOR IMPUTATION
# =============================================================================

# Separate columns to impute from those to exclude
data_to_impute <- main_data %>% select(all_of(all_features))
data_excluded  <- main_data %>% select(all_of(exclude_from_imputation))

# Convert binary/categorical variables to factors for MICE
# MICE uses logistic regression for factors, PMM for numeric — important distinction
binary_vars <- c(
  "gender", "high_bp", "high_choles", "diabetes_status",
  "arthritis", "chf", "chd", "heart_attack", "stroke",
  "liver", "thyroid", "copd", "overweight", "asthma", "cancer",
  "had_osteoporosis", "other_fractures", "prednisone",
  "mother_fx_hip", "father_fx_hip"
)

data_to_impute <- data_to_impute %>%
  mutate(across(all_of(binary_vars), as.factor))

cat("\nMissing data per group before imputation:\n")
cat("  Demo & Lifestyle:", sum(is.na(select(data_to_impute,
                                             all_of(group_demo_lifestyle)))), "\n")
cat("  Laboratory:      ", sum(is.na(select(data_to_impute,
                                             all_of(group_laboratory)))), "\n")
cat("  Bone Scans:      ", sum(is.na(select(data_to_impute,
                                             all_of(group_bone_scans)))), "\n")


# =============================================================================
# SECTION 4: RUN MICE IMPUTATION
# =============================================================================
# m = 1       : single imputed dataset (sufficient for ML, use 5+ for inference)
# method      : "pmm" for continuous (predictive mean matching)
#               "logreg" for binary factors (auto-detected)
# maxit = 10  : 10 iterations of the chained equations
# seed        : set for reproducibility

cat("\nRunning MICE imputation — this may take a few minutes...\n")

set.seed(42)
imputed <- mice(
  data_to_impute,
  m      = 1,
  maxit  = 10,
  method = "pmm",   # PMM works for both numeric and handles factors gracefully
  seed   = 42,
  printFlag = FALSE  # suppress iteration output; set TRUE to monitor progress
)

# Extract the single completed dataset
data_imputed <- complete(imputed, action = 1)

cat("Imputation complete.\n")
cat("Missing values after imputation:", sum(is.na(data_imputed)), "\n")


# =============================================================================
# SECTION 5: REASSEMBLE FULL DATASET
# =============================================================================

# Reattach excluded columns (IDs, weights, target) to imputed features
main_data_imputed <- bind_cols(data_excluded, data_imputed)

# Sanity check
stopifnot("Row count changed after imputation!" =
            nrow(main_data_imputed) == nrow(main_data))

cat("\nFinal imputed dataset dimensions:", dim(main_data_imputed), "\n")


# =============================================================================
# SECTION 6: SAVE FEATURE GROUP DATASETS
# =============================================================================
# Save both the full imputed set and pre-split group datasets
# Group datasets are used directly in 05_modelling.R

id_target_weights <- c("SEQN", "hx_fracture", "SDMVPSU", "SDMVSTRA",
                        "samp_wt_mec", "folate_wt", "samp_wt_intrvw")

# Full imputed dataset
write_csv(main_data_imputed, "data/processed/main_data_imputed.csv")

# Group 1: Demography & Lifestyle
main_data_imputed %>%
  select(all_of(id_target_weights), all_of(group_demo_lifestyle)) %>%
  write_csv("data/processed/group1_demo_lifestyle.csv")

# Group 2: Laboratory
main_data_imputed %>%
  select(all_of(id_target_weights), all_of(group_laboratory)) %>%
  write_csv("data/processed/group2_laboratory.csv")

# Group 3: Bone Scans
main_data_imputed %>%
  select(all_of(id_target_weights), all_of(group_bone_scans)) %>%
  write_csv("data/processed/group3_bone_scans.csv")

cat("\nSaved:\n")
cat("  data/processed/main_data_imputed.csv\n")
cat("  data/processed/group1_demo_lifestyle.csv\n")
cat("  data/processed/group2_laboratory.csv\n")
cat("  data/processed/group3_bone_scans.csv\n")
cat("Script 03 complete. Run 04_eda_visualisation.R next.\n")
