# -----------------------------------------------------------------------------
# FRACTURE RISK PREDICTION - NHANES 2017-2020
# SCRIPT 02: T-Score Calculation
# AUTHOR: Ernest Caballero
# DESCRIPTION: Calculates femur neck and Ward's Triangle T-scores for the
# 2017-2020 cohort (50+) using a young adult reference group
# (NHANES 2009-2010, age 20-30), stratified by gender and race.
# WHO defines T-score as: (patient BMD - young adult mean BMD) / young adult SD
# Classification criteria of osteoporosis T-score: 
# https://pmc.ncbi.nlm.nih.gov/articles/PMC10721754/
#    T-score >= -1.0            = Normal bone density
#    -2.5 < T-score < -1.0      = Osteopenia (low bone mass)
#    T-score   <=    -2.5       = Osteoporosis
# Reference: https://pmc.ncbi.nlm.nih.gov/articles/PMC2827823/
#------------------------------------------------------------------------------


# -----------------------------------------------------------
# SECTION 1: YOUNG ADULT REFERENCE GROUP (NHANES 2009-2010)
# -----------------------------------------------------------
# Reference group: age 20-30, stratified by gender (`RIAGENDR`) and harmonised
# race/ethnicity. Survey weights (`WTMEC2YR`) applied when computing reference
# means and SDs because NHANES oversamples minorities — unweighted stats would
# bias T-scores for all groups.
# WTMEC2YR = MEC examination weight, correct for DEXA physical exam data.
#
# Race/ethnicity harmonisation:
# RIDRETH1 (2009-2010) and RIDRETH3 (2017-2020) use different coding schemes.
# RIDRETH3 added Non-Hispanic Asian as a separate category (code 6) in 2011.
# Solution: collapse both variables to 5 comparable categories.
#
# Harmonised categories:
#   1 = Mexican American
#   2 = Other Hispanic
#   3 = Non-Hispanic White
#   4 = Non-Hispanic Black
#   5 = Other / Non-Hispanic Asian / Multiracial
# Source: CDC NHANES analytic guidelines for cross-cycle comparisons


# LIBRARIES
library(tidyverse)
library(survey)
library(foreign)

# Handle lonely PSUs (strata with single sampling unit after cohort subsetting)
# 'adjust' centres contribution around grand mean — standard practice for NHANES
options(survey.lonely.psu = 'adjust')


# --- Load demographic reference data (2009-2010) ---
download.file("https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2009/DataFiles/DEMO_F.XPT",
              tf <- tempfile(), mode = "wb")
demo_ref <- foreign::read.xport(tf)[, c("SEQN", "RIDAGEYR", "RIAGENDR",
                                        "RIDRETH1", "WTMEC2YR",
                                        "SDMVPSU", "SDMVSTRA")] %>%
  filter(RIDAGEYR >= 20 & RIDAGEYR <= 30) %>%
  # Retain weights and design variables — required for weighted reference stats
  select(SEQN, RIAGENDR, RIDAGEYR, RIDRETH1, WTMEC2YR, SDMVPSU, SDMVSTRA)

# Recode RIDRETH1 to harmonised race variable
demo_ref <- demo_ref %>%
  mutate(race_harmonised = case_when(
    RIDRETH1 == 1 ~ 1L,   # Mexican American
    RIDRETH1 == 2 ~ 2L,   # Other Hispanic
    RIDRETH1 == 3 ~ 3L,   # Non-Hispanic White
    RIDRETH1 == 4 ~ 4L,   # Non-Hispanic Black
    RIDRETH1 == 5 ~ 5L,   # Other / Asian / Multiracial (grouped in this cycle)
    TRUE          ~ NA_integer_
  ))

cat("Reference group size (age 20-30):", nrow(demo_ref), "\n")
cat("Missing race/ethnicity in reference:", sum(is.na(demo_ref$race_harmonised)), "\n")



# --- Load DEXA femur reference data (2009-2010) ---
download.file("https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2009/DataFiles/DXXFEM_F.XPT",
              tf <- tempfile(), mode = "wb")

dxa_femur_ref <- foreign::read.xport(tf)[, c("SEQN", "DXAFMRST", "DXXNKBMD")] %>%
  filter(DXAFMRST == 1) %>%    # valid scans only
  select(SEQN, DXXNKBMD)


# --- Merge demographic and DEXA reference data ---
ref_merged <- demo_ref %>%
  inner_join(dxa_femur_ref, by = "SEQN") %>%
  # Remove rows with missing BMD or weight — cannot compute weighted stats otherwise
  filter(!is.na(DXXNKBMD), !is.na(WTMEC2YR), WTMEC2YR > 0)

cat("Reference group after merge and filter:", nrow(ref_merged), "\n")


# --- Survey design for reference group ---
# NHANES complex survey design applied to the young adult reference sample.
# Accounts for clustering, stratification, and unequal sampling probabilities.
# https://wwwn.cdc.gov/nchs/nhanes/tutorials/SampleDesign.aspx
# https://wwwn.cdc.gov/nchs/nhanes/tutorials/Weighting.aspx
ref_design <- svydesign(
  id      = ~SDMVPSU,    # primary sampling unit
  strata  = ~SDMVSTRA,   # stratification variable
  weights = ~WTMEC2YR,   # MEC examination weight, correct for DEXA data
  data    = ref_merged,
  nest    = TRUE         # PSU IDs nested within strata, not globally unique
)


# --- Compute weighted reference means and SDs per gender & race stratum ---
# Split the young adult reference group by gender and race, 
# compute the survey-weighted mean femur neck BMD for each subgroup, 
# convert the SE into a SD using the stratum sample size

reference_stats <- svyby(
  formula = ~DXXNKBMD,                     # femur neck
  by      = ~RIAGENDR + race_harmonised,   # stratify by gender and harmonised race
  design  = ref_design,
  FUN     = svymean,                       # svymean returns weighted mean and SE per stratum
  na.rm   = TRUE
) %>%
  as_tibble() %>%
  rename(
    mean_ref_nk = DXXNKBMD,
    se_ref_nk   = se) %>%
  
  # Join unweighted stratum counts to calculate SD from SE
  left_join(
    ref_merged %>%
      group_by(RIAGENDR, race_harmonised) %>%
      summarise(n_ref = n(), .groups = "drop"),
    by = c("RIAGENDR", "race_harmonised")
  ) %>%
  mutate(
    # SD back-calculated from SE
    sd_ref_nk = se_ref_nk * sqrt(n_ref)   # SD = SE × sqrt(n), converts survey standard error to population SD
  ) %>%
  select(RIAGENDR, race_harmonised, mean_ref_nk, sd_ref_nk, n_ref)

cat("\nWeighted reference statistics by gender & race stratum:\n")
print(reference_stats)



# Sanity check — SD must be positive and non-NA for all strata
# Zero or NA SD means too few respondents in that stratum for reliable reference
problem_strata <- reference_stats %>%
  filter(is.na(sd_ref_nk) | sd_ref_nk == 0)

if (nrow(problem_strata) > 0) {
  cat("\nWARNING: Following strata have unusable SD — too few respondents:\n")
  print(problem_strata %>% select(RIAGENDR, race_harmonised, n_ref))
  cat("Consider collapsing race categories for affected strata.\n")
} else {
  cat("\nAll strata have valid weighted means and SDs.\n")
}



# -------------------------------------------------------------
# SECTION 2: PREPARE STUDY COHORT (NHANES 2017-2020, age 50+)
# -------------------------------------------------------------

download.file("https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2017/DataFiles/P_DEMO.XPT",
              tf <- tempfile(), mode = "wb")

demo_study <- foreign::read.xport(tf)[, c("SEQN", "RIDAGEYR", "RIAGENDR",
                                          "RIDRETH3", "INDFMPIR", "WTINTPRP",
                                          "WTMECPRP", "SDMVPSU", "SDMVSTRA")] %>%
  filter(RIDAGEYR >= 50) %>%
  # FIX: Retain weights and design variables for downstream survey analysis
  # in 01_data_prep.R — previous version dropped these in a second select()
  select(SEQN, RIAGENDR, RIDAGEYR, RIDRETH3, INDFMPIR,
         WTINTPRP, WTMECPRP, SDMVPSU, SDMVSTRA)

cat("\nStudy cohort size (age 50+):", nrow(demo_study), "\n")

# Recode RIDRETH3 → harmonised race variable
# RIDRETH3 codes: 1=Mexican American, 2=Other Hispanic, 3=NH White,
#                 4=NH Black, 6=NH Asian, 7=Other/Multiracial
# NOTE: RIDRETH3 does not use code 5 — it was skipped when Asian (6) was added
# FIX: Previous version had RIDRETH3 %in% c(6,7) but the comment said
#      code 5 was missed — confirmed: code 5 is simply unused in RIDRETH3
demo_study <- demo_study %>%
  mutate(race_harmonised = case_when(
    RIDRETH3 == 1          ~ 1L,   # Mexican American
    RIDRETH3 == 2          ~ 2L,   # Other Hispanic
    RIDRETH3 == 3          ~ 3L,   # Non-Hispanic White
    RIDRETH3 == 4          ~ 4L,   # Non-Hispanic Black
    RIDRETH3 %in% c(6, 7)  ~ 5L,   # Non-Hispanic Asian + Other/Multiracial
    TRUE                   ~ NA_integer_
  ))

cat("Missing race/ethnicity in study cohort:",
    sum(is.na(demo_study$race_harmonised)), "\n")


# --- Load study cohort DEXA femur data ---
download.file("https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2017/DataFiles/P_DXXFEM.XPT",
              tf <- tempfile(), mode = "wb")

dexa_femur_study <- foreign::read.xport(tf)[, c("SEQN", "DXAFMRST", "DXXNKBMD")] %>%
  filter(DXAFMRST == 1) %>%    # valid scans only
  select(SEQN, DXXNKBMD)

cat("Study participants with valid DEXA:", nrow(dexa_femur_study), "\n")



# -------------------------------
# SECTION 3: CALCULATE T-SCORES
# -------------------------------
# Join study cohort with reference stats by matching gender & race stratum, then apply T-score formula.
# Classification based on femur neck T-score https://pubmed.ncbi.nlm.nih.gov/17223616/

n_before_dexa <- nrow(demo_study)

t_scores <- demo_study %>%
  inner_join(dexa_femur_study, by = "SEQN") %>%   # retain only those with DEXA
  left_join(reference_stats,
            by = c("RIAGENDR", "race_harmonised")) %>%   # match to reference stratum
  mutate(
    # T-score formula: (individual BMD - young adult reference mean) / reference SD
    T_score_nk = (DXXNKBMD - mean_ref_nk) / sd_ref_nk,
    
    # Classification
    osteo_class = case_when(
      T_score_nk >= -1.0  ~ "Normal",
      T_score_nk > -2.5   ~ "Osteopenia",
      T_score_nk <= -2.5  ~ "Osteoporosis",
      TRUE                ~ NA_character_   # uncomputable if reference SD missing
    )
  ) %>%
  select(SEQN, T_score_nk, osteo_class)

# Participant flow report
cat("\nParticipant flow:\n")
cat("Study cohort (age 50+):", n_before_dexa, "\n")
cat("After DEXA join (valid scan only):", nrow(t_scores), "\n")
cat("Missing (no scan performed):", n_before_dexa - nrow(t_scores), "\n")
cat("T-score NA (unmatched reference):", sum(is.na(t_scores$T_score_nk)), "\n")



# ---------------------------------
# SECTION 4: QUALITY CHECK & SAVE
# ---------------------------------

cat("\nT-score summary:\n")
print(summary(t_scores[, "T_score_nk"]))

cat("\nClassification — femur neck T-score:\n")
print(table(t_scores$osteo_class, useNA = "ifany"))

# Save
write_csv(t_scores, "data/processed/dexa_t_scores.csv")
cat("\nSaved: data/processed/dexa_t_scores.csv\n")



# ----------------------------------------
# SECTION 5: MERGE WITH MAIN DATA & SAVE
# ----------------------------------------

main_data_clean <- read_csv("data/processed/main_data_clean.csv",
                            show_col_types = FALSE)

cat("Rows before T-score join:", nrow(main_data_clean), "\n")

# Sanity check: verify all T-score SEQNs exist in main_data_clean before joining
# If any SEQN in t_scores is not in main_data_clean, something went wrong upstream
unmatched_seqn <- sum(!t_scores$SEQN %in% main_data_clean$SEQN)
cat("T-score SEQNs not found in main data:", unmatched_seqn, "\n")

# Merge T-scores into main dataset
# left_join: keeps ALL participants from main dataset, including those without a valid DEXA scan. 
# T-score columns (T_score_nk, osteo_class) will be NA for those participants — MICE will impute these alongside other missing values.
# Stratum matching (gender & race) was already applied, so joining on SEQN alone is correct here.
main_data_final <- main_data_clean %>%
  left_join(t_scores, by = "SEQN")

cat("Rows after T-score join:", nrow(main_data_final), "\n")
cat("Participants missing T-score (NA):", sum(is.na(main_data_final$T_score_nk)), "\n")
cat("Total missing values:", sum(is.na(main_data_final)), "\n")
cat("Dimensions of the merged dataset:", dim(main_data_final), "\n")

# Save
write_csv(main_data_final, "data/processed/main_data_final.csv")
cat("\nSaved: data/processed/main_data_final.csv\n")

cat("\nScript 02 complete. Run script 03 next. \n")
