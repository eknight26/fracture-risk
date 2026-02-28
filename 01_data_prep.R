
# FRACTURE RISK PREDICTION - NHANES 2017-2020
# Script 01: Data Preparation & Merging
# Author: Ernest Caballero
# Date: 2024
# Description: Loads, cleans, and merges all NHANES XPT modules into a single
# analysis-ready dataset. Target variable is history of fracture
# at hip, wrist, or spine (hx_fracture).


# Clear environment and set options
rm(list = ls())
options(scipen = 999)

# LIBRARIES
# tidyverse covers: ggplot2, dplyr, tidyr, readr, purrr, stringr, forcats
library(dplyr)
library(survey)
library(tidyverse)
# library(haven)     # read XPT (SAS transport) files
# library(janitor)   # clean_names(), tabyl()
# library(naniar)    # missing data visualisation
# library(pander)    # formatted table output
# library(knitr)     # knitting reports

# Handle lonely PSUs (strata with single sampling unit) after cohort subsetting
# 'adjust' centres contribution around grand mean — standard practice for NHANES
options(survey.lonely.psu = 'adjust')

# Display Version Information
cat("R package versions:\n")
for (p in c("base", "survey","dplyr")) { 
  cat(p, ": ", as.character(packageVersion(p)), "\n")
}

# Create output folders if they don't exist
dir.create("outputs",       showWarnings = FALSE)
dir.create("outputs/plots", showWarnings = FALSE, recursive = TRUE)
dir.create("data/processed", showWarnings = FALSE, recursive = TRUE)



# ---------------------------------------------- 
# SECTION 1: LOAD & PREPARE INDIVIDUAL MODULES
# ----------------------------------------------

# --- Demographics ---
# SEQN = participant sequence number, RIDAGEYR = age, RIAGENDR = gender, 
# RIDRETH3 = race/ethnicity, INDFMPIR = income-to-poverty ratio (capped at 5.0)
# WTINTPRP = full sample interview weight, WTMECPRP = full sample MEC exam weight
# SDMVPSU  = Masked variance unit pseudo-PSU variable for variance estimation
# SDMVSTRA = Masked variance unit pseudo-stratum variable for variance estimation

download.file("https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2017/DataFiles/P_DEMO.XPT", tf <- tempfile(), mode="wb")
demo <- foreign::read.xport(tf)[,c("SEQN", "RIDAGEYR", "RIAGENDR", 
                                   "RIDRETH3", "INDFMPIR", "WTINTPRP", 
                                   "WTMECPRP", "SDMVPSU", "SDMVSTRA")] %>%
  select(SEQN, RIDAGEYR, RIAGENDR, RIDRETH3, INDFMPIR, WTINTPRP, 
         WTMECPRP, SDMVPSU, SDMVSTRA) %>%
  rename(
    age            = RIDAGEYR,
    gender         = RIAGENDR,
    ethnicity      = RIDRETH3,
    income_ratio   = INDFMPIR,
    samp_wt_intrvw = WTINTPRP,
    samp_wt_mec    = WTMECPRP
  )


# --- Osteoporosis (Target Variable) ---
# Self-reported information regarding fractured hip, wrist, or spine
# OSQ010A = fractured a hip, OSQ010B = fractured a wrist, OSQ010C = fractured spine
# hx_fracture = 1 if any fracture at hip, wrist, or spine; else 0

download.file("https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2017/DataFiles/P_OSQ.XPT", tf <- tempfile(), mode="wb")
osteo <- foreign::read.xport(tf)[,c("SEQN", "OSQ060", "OSQ080", "OSQ130", 
                                    "OSQ170", "OSQ200", "OSQ010A", "OSQ010B", 
                                    "OSQ010C")] %>%
  select(SEQN, OSQ060, OSQ080, OSQ130, OSQ170, OSQ200,
         OSQ010A, OSQ010B, OSQ010C) %>%
  rename(
    had_osteoporosis = OSQ060,    # ever told you had osteoporosis?
    other_fractures  = OSQ080,    # doctor ever told any other fractures?
    prednisone       = OSQ130,    # ever taken prednisone daily?
    mother_fx_hip    = OSQ170,    # did mother ever fracture a hip
    father_fx_hip    = OSQ200,    # did father ever fracture a hip
    fx_hip           = OSQ010A,   # broken/fractured a hip
    fx_wrist         = OSQ010B,   # broken/fractured a wrist
    fx_spine         = OSQ010C    # broken/fractured spine
    
  ) %>%
  mutate(
    hx_fracture = if_else(
      rowSums(across(c(fx_hip, fx_wrist, fx_spine)) == 1, na.rm = TRUE) > 0,
      1L, 0L
    )
  ) %>%
  select(-fx_hip, -fx_wrist, -fx_spine)  # drop source columns after deriving target


# --- DEXA Femur (Bone Mineral Density) ---
# DXAFMRST == 1 = valid scan; selected site-specific BMD columns

download.file("https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2017/DataFiles/P_DXXFEM.XPT", tf <- tempfile(), mode="wb")
dexa_femur <- foreign::read.xport(tf)[,c("SEQN", "DXAFMRST", "DXXOFBMD", 
                                         "DXXNKBMD", "DXXTRBMD", "DXXINBMD", 
                                         "DXXWDBMD")] %>%
  filter(DXAFMRST == 1) %>%
  select(SEQN, DXXOFBMD, DXXNKBMD, DXXTRBMD, DXXINBMD, DXXWDBMD) %>%
  rename(
    bmd_femur_total    = DXXOFBMD,
    bmd_femur_neck     = DXXNKBMD,
    bmd_femur_troch    = DXXTRBMD,
    bmd_femur_inter    = DXXINBMD,
    bmd_femur_ward     = DXXWDBMD
  )


# --- DEXA Spine (Bone Mineral Density) ---
# DXASPNST == 1 = complete valid scan
# DXASPNST == 2 = complete but invalid (included; see note below)
# NOTE: Invalid scans (code 2) are retained because the main causes —
# degenerative disease, severe scoliosis, spinal fusion — are prevalent
# in the 50+ population and exclusion would bias the sample.

download.file("https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2017/DataFiles/P_DXXSPN.XPT", tf <- tempfile(), mode="wb")
dexa_spine <- foreign::read.xport(tf)[,c("SEQN", "DXASPNST", "DXXOSBMD", "DXXL1BMD", 
                                         "DXXL2BMD", "DXXL3BMD", "DXXL4BMD")] %>%
  filter(DXASPNST %in% c(1, 2)) %>%
  select(SEQN, DXXOSBMD, DXXL1BMD, DXXL2BMD, DXXL3BMD, DXXL4BMD) %>%
  rename(
    bmd_spine_total = DXXOSBMD,
    bmd_l1          = DXXL1BMD,
    bmd_l2          = DXXL2BMD,
    bmd_l3          = DXXL3BMD,
    bmd_l4          = DXXL4BMD
  )


# --- Body Measurements ---
# BMDSTATS == 1 = valid body measurement exam
download.file("https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2017/DataFiles/P_BMX.XPT", tf <- tempfile(), mode="wb")
bmx <- foreign::read.xport(tf)[,c("SEQN", "BMDSTATS", "BMXBMI")] %>%
  filter(BMDSTATS == 1) %>%
  select(SEQN, BMXBMI) %>%
  rename(bmi = BMXBMI)


# --- Blood Pressure ---
# Mean systolic and diastolic BP calculated from 3 readings (BPXOSY1-3)
download.file("https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2017/DataFiles/P_BPXO.XPT", tf <- tempfile(), mode="wb")
bp <- foreign::read.xport(tf)[,c("SEQN", "BPXOSY1", "BPXOSY2", "BPXOSY3", 
                                 "BPXODI1", "BPXODI2", "BPXODI3")] %>%
  mutate(mean_systolic = rowMeans(select(., BPXOSY1, BPXOSY2, BPXOSY3),
                                     na.rm = TRUE)) %>%
  mutate(mean_diastolic = rowMeans(select(., BPXODI1, BPXODI2, BPXODI3),
                                   na.rm = TRUE)) %>%
  select(SEQN, mean_systolic, mean_diastolic)


# --- BP and Cholesterol (Hypertension Diagnosis) ---
download.file("https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2017/DataFiles/P_BPQ.XPT", tf <- tempfile(), mode="wb")
bp_chol <- foreign::read.xport(tf)[,c("SEQN", "BPQ020", "BPQ080")] %>%
  select(SEQN, BPQ020, BPQ080) %>%
  rename(
    high_bp     = BPQ020,
    high_choles = BPQ080
  )



# --- Alcohol Intake ---
# ALQ111 = ever drank alcohol (1 = yes); ALQ130 = avg drinks/day last 12 months
# alcohol_consumed = drinks/day if drinker, else 0

download.file("https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2017/DataFiles/P_ALQ.XPT", tf <- tempfile(), mode="wb")
alc_intake <- foreign::read.xport(tf)[,c("SEQN", "ALQ111", "ALQ130")] %>%
  select(SEQN, ALQ111, ALQ130) %>%
  mutate(alc_drinks_per_day = if_else(ALQ111 == 1, ALQ130, 0)) %>%
  select(SEQN, alc_drinks_per_day)



# --- Diabetes ---
# DIQ010: 1 = yes diabetic, 2 = no, 3 = borderline
# DIQ280: aic level, 777 = refused, 9999 = don't know

download.file("https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2017/DataFiles/P_DIQ.XPT", tf <- tempfile(), mode="wb")
diabetes <- foreign::read.xport(tf)[,c("SEQN", "DIQ010", "DIQ280", "DID320")] %>%
  select(SEQN, DIQ010, DIQ280, DID320) %>%
  rename(
    diabetes_status = DIQ010,
    a1c_level       = DIQ280,
    ldl_level       = DID320
  )



# --- Medical Conditions ---
# Comorbidities: arthritis, CHF, coronary heart disease, stroke, emphysema,
# thyroid, COPD, osteoporosis, overweight, asthma, cancer
# 1 = yes, 2 = no, 7 = refused, 9 = don't know

download.file("https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2017/DataFiles/P_MCQ.XPT", tf <- tempfile(), mode="wb")
med_cond <- foreign::read.xport(tf)[,c("SEQN", "MCQ160A", "MCQ160B", "MCQ160C", 
                                       "MCQ160E", "MCQ160F", "MCQ160L", 
                                       "MCQ160M", "MCQ160P", "MCQ080", 
                                       "MCQ010", "MCQ220")] %>%
  select(SEQN, MCQ160A, MCQ160B, MCQ160C, MCQ160E, MCQ160F,
         MCQ160L, MCQ160M, MCQ160P, MCQ080, MCQ010, MCQ220) %>%
  rename(
    arthritis    = MCQ160A,
    chf          = MCQ160B,
    chd          = MCQ160C,
    heart_attack = MCQ160E,
    stroke       = MCQ160F,
    liver        = MCQ160L,
    thyroid      = MCQ160M,
    copd         = MCQ160P,
    overweight   = MCQ080,
    asthma       = MCQ010,
    cancer       = MCQ220
  )



# --- Physical Activity (Sedentary Behaviour) ---
# PAD680 = minutes of sedentary activity per day

download.file("https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2017/DataFiles/P_PAQ.XPT", tf <- tempfile(), mode="wb")
activity <- foreign::read.xport(tf)[,c("SEQN", "PAD680")] %>%
  select(SEQN, PAD680) %>%
  rename(sedentary_mins_day = PAD680)



# --- Smoking ---
# SMQ681 = smoked in last 5 days (1 = yes); SMQ720 = cigarettes/day last 5 days
# 95 = if more than 95 cigarettes/day (ceiling)

download.file("https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2017/DataFiles/P_SMQRTU.XPT", tf <- tempfile(), mode="wb")
smoking <- foreign::read.xport(tf)[,c("SEQN", "SMQ681", "SMQ720")] %>%
  select(SEQN, SMQ681, SMQ720) %>%
  mutate(cigarettes_per_day = if_else(SMQ681 == 1, SMQ720, 0)) %>%
  select(SEQN, cigarettes_per_day)



# --- Blood Tests: Folate ---
# WTFOLPRP: Folate Weight pre-pandemic - weight adjustment due to oversampling, 
# LBDFOTSI: serum total folate (SI unit nmol/L)

download.file("https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2017/DataFiles/P_FOLFMS.XPT", tf <- tempfile(), mode="wb")
folate <- foreign::read.xport(tf)[,c("SEQN", "WTFOLPRP", "LBDFOTSI")] %>%
  select(SEQN, LBDFOTSI, WTFOLPRP) %>%
  rename(
    serum_folate     = LBDFOTSI,
    folate_wt        = WTFOLPRP     
  )



# --- Blood Tests: Biochemistry (ALP, Phosphate, Calcium) ---
# Selected for bone metabolism relevance:
# ALP = bone/liver activity marker; phosphate & calcium = bone mineralisation

download.file("https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2017/DataFiles/P_BIOPRO.XPT", tf <- tempfile(), mode="wb")
biochem <- foreign::read.xport(tf)[,c("SEQN", "LBXSAPSI", "LBDSPHSI", "LBDSCASI")] %>%
  select(SEQN, LBXSAPSI, LBDSPHSI, LBDSCASI) %>%
  rename(
    alp        = LBXSAPSI,
    phosphate  = LBDSPHSI,
    calcium    = LBDSCASI
  )



# --- Blood Tests: Heavy Metals (Lead, Cadmium, Mercury, Selenium, Manganese) ---
# Environmental exposures linked to bone density loss

download.file("https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2017/DataFiles/P_PBCD.XPT", tf <- tempfile(), mode="wb")
metals <- foreign::read.xport(tf)[,c("SEQN", "LBDBPBSI", "LBDBCDSI", "LBDTHGSI", 
                                     "LBDBSESI", "LBDBMNSI")] %>%
  select(SEQN, LBDBPBSI, LBDBCDSI, LBDTHGSI, LBDBSESI, LBDBMNSI) %>%
  rename(
    blood_lead      = LBDBPBSI,
    blood_cadmium   = LBDBCDSI,
    blood_mercury   = LBDTHGSI,
    blood_selenium  = LBDBSESI,
    blood_manganese = LBDBMNSI
  )



# ------------------------------
# SECTION 2: MERGE ALL MODULES
# ------------------------------

# Base: demographics inner joined with osteoporosis (to ensure only respondents
# with both demographic and fracture data are retained)
main_data <- demo %>%
  inner_join(osteo, by = "SEQN")

# All supplementary datasets
# left_join preserves all rows from main_data; unmatched rows get NA
supplementary <- list(
  dexa_femur, dexa_spine,
  bmx, bp, bp_chol,
  folate, biochem, metals,
  alc_intake, diabetes,
  med_cond, activity, smoking
)

# purrr::reduce applies left_join sequentially — cleaner than 13 manual joins
main_data <- reduce(supplementary, ~left_join(.x, .y, by = "SEQN"),
                    .init = main_data)

cat("Dataset dimensions after merge:", dim(main_data), "\n")


# ----------------------------------------
# SECTION 3: RECODE SPECIAL VALUES TO NA
# ----------------------------------------
# NHANES encodes refused (7/77/777) and don't know (9/99/999/9999) as numeric.
# Must be set to NA before imputation to prevent corruption of data

# Binary/categorical variables (1=yes, 2=no, 7=refused, 9=don't know)
binary_vars <- c(
  "had_osteoporosis", "other_fractures", "prednisone",
  "mother_fx_hip", "father_fx_hip",
  "high_bp", "high_choles", "diabetes_status",
  "arthritis", "chf", "chd", "heart_attack", "stroke",
  "liver", "thyroid", "copd", "overweight", "asthma", "cancer"
)

main_data <- main_data %>%
  mutate(across(all_of(binary_vars), ~na_if(., 7))) %>%   # refused
  mutate(across(all_of(binary_vars), ~na_if(., 9))) %>%   # don't know
  # Recode binary: 1 = yes, 2 = no → convert 2 to 0 for modelling
  mutate(across(all_of(binary_vars), ~if_else(. == 2, 0L, as.integer(.))))

# Continuous variables with special codes (777/7777 = refused, 999/9999 = don't know)
main_data <- main_data %>%
  mutate(
    sedentary_mins_day = na_if(sedentary_mins_day, 9999),  # don't know
    sedentary_mins_day = na_if(sedentary_mins_day, 7777),  # refused
    a1c_level          = na_if(a1c_level, 999),
    a1c_level          = na_if(a1c_level, 777),
    ldl_level          = na_if(ldl_level, 5555),           # never heard of LDL
    ldl_level          = na_if(ldl_level, 6666),           # never had cholesterol test
    ldl_level          = na_if(ldl_level, 7777),
    ldl_level          = na_if(ldl_level, 9999),
    alc_drinks_per_day = na_if(alc_drinks_per_day, 999),
    alc_drinks_per_day = na_if(alc_drinks_per_day, 777),
    cigarettes_per_day = na_if(cigarettes_per_day, 999),   
    cigarettes_per_day = na_if(cigarettes_per_day, 777)
  )



# ----------------------------
# SECTION 4: QUALITY CHECKS
# ----------------------------

# Check for duplicates
stopifnot("Duplicate SEQNs detected!" =
            length(unique(main_data$SEQN)) == nrow(main_data))
cat("Check passed: No duplicate respondents.\n")

# Remove single row where target variable is NA (unresolvable)
main_data <- main_data %>% drop_na(hx_fracture)
cat("Rows after removing NA target:", nrow(main_data), "\n")

# Check class label distribution — expect imbalance (~85% no fracture, ~15% fracture)
cat("\nClass label distribution (hx_fracture):\n")
print(table(main_data$hx_fracture))
cat("Class imbalance ratio:", round(sum(main_data$hx_fracture == 0) /
                                      sum(main_data$hx_fracture == 1), 1), ":1\n")

# Convert target to factor
main_data <- main_data %>%
  mutate(hx_fracture = factor(hx_fracture, levels = c(0, 1),
                              labels = c("No_Fracture", "Fracture")))



# ----------------------------------------
# SECTION 5: MISSING DATA SUMMARY & PLOT
# ----------------------------------------

na_summary <- main_data %>%
  summarise(across(everything(), ~sum(is.na(.)))) %>%
  # rotates output so each column becomes its own row
  pivot_longer(everything(), names_to = "feature", values_to = "n_missing") %>%
  # Adds new column `pct_missing` to get a percentage of missing data for each feature
  mutate(pct_missing = round(n_missing / nrow(main_data) * 100, 1)) %>%
  arrange(desc(pct_missing))

write_csv(na_summary, "outputs/na_summary.csv")
cat("\nTop 20 features with most missing data:\n")
print(head(na_summary, 20))

# Plot missing data — only features with any missing values
na_summary %>%
  filter(pct_missing > 0) %>%
  ggplot(aes(x = reorder(feature, pct_missing), y = pct_missing)) +
  geom_col(fill = "steelblue4", alpha = 0.8) +
  geom_text(aes(label = paste0(pct_missing, "%")),
            hjust = -0.1, size = 3, colour = "black") +
  coord_flip() +
  labs(
    x        = NULL,
    y        = "Missing Values (%)",
    title    = "Missing Data by Feature",
    subtitle = "NHANES 2017-2020 Pre-Pandemic | Fracture Risk Cohort",
    caption  = "Source: CDC/NCHS NHANES"
  ) +
  theme_minimal(base_size = 12) +
  theme(axis.text = element_text(colour = "black")) +
  scale_y_continuous(limits = c(0, 100), expand = expansion(mult = c(0, 0.1)))

ggsave("outputs/plots/01_missing_data.png", width = 10, height = 8, dpi = 300)



# -----------------
# SECTION 6: SAVE
# -----------------

write_csv(main_data, "data/processed/main_data_clean.csv")
cat("\nSaved: data/processed/main_data_clean.csv\n")
cat("Script 01 complete. Run 02_t_score_calc.R next.\n")

