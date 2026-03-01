# -----------------------------------------------------------------------------
# FRACTURE RISK PREDICTION - NHANES 2017-2020
# SCRIPT 03: Missing Data Imputation
# AUTHOR: Ernest Caballero
# DESCRIPTION: Imputes missing values using MICE (Multiple Imputation by
# Chained Equations). A single completed dataset is produced (m=1 imputation) 
# suitable for ML modelling. Survey weights and ID variables are excluded 
# from imputation but retained.
# METHOD RATIONALE: MICE is preferred over simple median/mean imputation 
# because it preserves relationships between variables. 
# For a clinical dataset with structured missingness (DEXA scans only on eligible participants), 
# this is important. m=1 is used here for ML pipelines; for inferential statistics use m=5+.
# -----------------------------------------------------------------------------

# LIBRARIES
library(tidyverse)
library(mice)      # multiple imputation by chained equations (MICE)


# ---------------------------------
# SECTION 1: LOAD PROCESSED DATA
# ---------------------------------

main_data_imp <- read_csv("data/processed/main_data_final.csv",
                          show_col_types = FALSE)

cat("Dimensions:", dim(main_data_imp), "\n")
cat("Total missing values:", sum(is.na(main_data_imp)), "\n")
print(head(main_data_imp))



# ----------------------------------
# SECTION 2: DEFINE FEATURE GROUPS
# ----------------------------------
# We define three clinically meaningful groups to reflect how fracture risk is assessed:
# Group 1 — Demography & Lifestyle: who the patient is and how they live
# Group 2 — Laboratory: what their blood reveals about bone metabolism
# Group 3 — Bone Scans: what DEXA imaging and T-scores directly show

group_demo_lifestyle <- c(
  # Demographics
  "age", "gender", "ethnicity", "income_ratio",
  
  # Body measurements & vitals
  "bmi", "mean_systolic", "mean_diastolic",
  
  # Comorbidities
  "high_bp", "high_choles", "diabetes_status", 
  "arthritis", "chf", "chd", "heart_attack", 
  "stroke", "liver", "thyroid", "copd", 
  "overweight", "asthma", "cancer",
  
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
  "ldl_level", "a1c_level"
)

group_bone_scans <- c(
  # Femur BMD
  "bmd_femur_total", "bmd_femur_neck", "bmd_femur_troch",
  "bmd_femur_inter", "bmd_femur_ward",
  
  # Spine BMD
  "bmd_spine_total", "bmd_l1", "bmd_l2", "bmd_l3", "bmd_l4",
  
  # T-scores (derived in script 02)
  "T_score_nk", "osteo_class"
)

# Variables never imputed — identifiers, weights, target
exclude_from_imputation <- c("SEQN", "SDMVPSU", "SDMVSTRA", 
                             "samp_wt_intrvw", "samp_wt_mec", "folate_wt",
                             "hx_fracture")

# All features to impute
all_features <- c(group_demo_lifestyle, group_laboratory, group_bone_scans)

cat("\nFeature group sizes:\n")
cat("Demo & Lifestyle:", length(group_demo_lifestyle), "features\n")
cat("Laboratory:", length(group_laboratory), "features\n")
cat("Bone Scans:", length(group_bone_scans), "features\n")
cat("Total:", length(all_features), "features\n")



# ----------------------------------------
# SECTION 3: PREPARE DATA FOR IMPUTATION
# ----------------------------------------

# Select columns
data_to_impute <- main_data_imp %>% select(all_of(all_features))
data_excluded  <- main_data_imp %>% select(all_of(exclude_from_imputation))

# Define imputation method per feature type
# Binary variables (0/1) — logistic regression
binary_vars <- c("had_osteoporosis", "other_fractures", "prednisone",
                 "mother_fx_hip", "father_fx_hip", "high_bp", "high_choles", 
                 "arthritis", "chf", "chd", "heart_attack", "stroke", "liver", 
                 "thyroid", "copd", "overweight", "asthma", "cancer")

# Unordered categorical (3+ categories) — polytomous logistic regression
poly_vars <- c("diabetes_status", "osteo_class")

# Continuous numeric — predictive mean matching
pmm_vars <- c( "income_ratio", "bmd_femur_total", "bmd_femur_neck", 
               "bmd_femur_troch", "bmd_femur_inter", "bmd_femur_ward", 
               "bmd_spine_total", "bmd_l1", "bmd_l2", "bmd_l3", "bmd_l4", 
               "bmi", "mean_systolic", "mean_diastolic", "serum_folate", 
               "alp", "phosphate", "calcium", "blood_lead", "blood_cadmium", 
               "blood_mercury", "blood_selenium", "blood_manganese",
               "ldl_level", "a1c_level", "alc_drinks_per_day", 
               "cigarettes_per_day", "sedentary_mins_day", "T_score_nk")

data_to_impute <- data_to_impute %>%
  mutate(gender = factor(gender,
                         levels = c(1, 2),
                         labels = c("Male", "Female"))) %>%
  mutate(ethnicity = factor(ethnicity,
                            levels = c(1, 2, 3, 4, 6, 7),
                            labels = c("Mexican American", "Other Hispanic",
                                       "Non-Hispanic White", "Non-Hispanic Black",
                                       "Non-Hispanic Asian", "Other/Multiracial"))) %>%
  mutate(across(all_of(binary_vars),
                ~ factor(.x, levels = c(0, 1), labels = c("No", "Yes")))) %>%
  mutate(diabetes_status = factor(diabetes_status,
                                  levels = c(1, 2, 3),        
                                  labels = c("Yes", "No", "Borderline"))) %>% 
  mutate(osteo_class = factor(osteo_class,
                              levels = c("Normal", "Osteopenia", "Osteoporosis")))

cat("\nSpot check factor levels:\n")
cat("Gender:\n"); print(levels(data_to_impute$gender))
cat("Diabetes status:\n"); print(levels(data_to_impute$diabetes_status))
cat("Osteo class:\n"); print(levels(data_to_impute$osteo_class))
cat("Arthritis:\n"); print(levels(data_to_impute$arthritis))


# --- Missing data counts for each column --- 
na_check <- data_to_impute %>%
  summarise(across(everything(), ~ sum(is.na(.)))) %>%
  pivot_longer(everything(),
               names_to  = "feature",
               values_to = "n_missing") %>%
  mutate(pct_missing = round(n_missing / nrow(data_to_impute) * 100, 2)) %>%
  arrange(desc(n_missing))

cat("Columns with missing values:\n")
print(na_check %>% filter(n_missing > 0))


# --- Columns to DROP from imputation (>80% missing) ---
# Too little observed data — imputed values would be unreliable
drop_from_imputation <- c("ldl_level", "a1c_level")

# Remove columns from data_to_impute and pmm_vars
data_to_impute <- data_to_impute %>%
  select(-all_of(drop_from_imputation))
pmm_vars <- setdiff(pmm_vars, drop_from_imputation)
group_laboratory <- setdiff(group_laboratory, drop_from_imputation)

cat("Dropped from imputation (>80% missing):", drop_from_imputation, "\n")
cat("Remaining columns:", ncol(data_to_impute), "\n")


# --- Columns to FLAG (50-80% missing) ---
# MICE will impute these (to document high uncertainty)
high_missing_flag <- c("bmd_spine_total", "serum_folate")

cat("\nHigh missingness columns (imputed with caution):\n")
cat("bmd_spine_total: 57.5% missing\n")
cat("serum_folate: 55.9% missing\n")


# Check missing data numbers per group
cat("\nMissing data per group before imputation:\n")
cat("Demo & Lifestyle:", sum(is.na(select(data_to_impute, all_of(group_demo_lifestyle)))), "\n")
cat("Laboratory:", sum(is.na(select(data_to_impute, all_of(group_laboratory)))), "\n")
cat("Bone Scans:", sum(is.na(select(data_to_impute, all_of(group_bone_scans)))), "\n")

# Build method vector
method_vector <- rep("", ncol(data_to_impute))
names(method_vector) <- names(data_to_impute)

method_vector[names(method_vector) %in% pmm_vars] <- "pmm"
method_vector[names(method_vector) %in% binary_vars] <- "logreg"
method_vector[names(method_vector) %in% poly_vars] <- "polyreg"

cat("Imputation methods assigned:\n")
print(table(method_vector))



# --------------------------------
# SECTION 4: RUN MICE IMPUTATION
# --------------------------------
# m = 5: 5 imputed datasets — needed for proper inference (Rubin's rules)
# also good practice for ML to average predictions across datasets
# maxit = 10: 10 iterations of chained equations per imputation
# seed = 42 : reproducibility

cat("\nRunning MICE imputation ......\n")

imputed <- mice(
  data_to_impute,
  m         = 5,          # 5 datasets for both inference and ML
  maxit     = 30,
  method    = method_vector,
  seed      = 42,
  printFlag = TRUE        # set TRUE to monitor iteration progress
)

cat("Imputation complete.\n")


# --- Inspect convergence ---
# Trace plots show the mean and SD of imputed values across iterations mix well (no clear trend upward/downward),
# which indicates the algorithm has converged and imputations are stable.
plot(imputed)


# --- Check imputation quality ---
# Imputed values follow a similar shape to observed data
densityplot(imputed)


# --- Extract completed datasets ---
# For ML: use a single completed dataset (action = 1) to train your model
# For inference: pool results across all 5 datasets using Rubin's rules:
# fit model on each, then pool with mice::pool()

# Dataset for ML
data_imputed_ml <- complete(imputed, action = 1)
# Dataset for inference, for pooling
data_imputed_all <- complete(imputed, action = "long", include = FALSE) 

cat("Missing values after imputation (dataset 1):", sum(is.na(data_imputed_ml)), "\n")
cat("Rows in stacked dataset (5 x n):", nrow(data_imputed_all), "\n")




# ------------------------------------
# SECTION 5: REASSEMBLE FULL DATASET
# ------------------------------------

# Reattach excluded columns to imputed features
main_data_imputed <- bind_cols(data_excluded, data_imputed)

# Sanity check
stopifnot("Row count changed after imputation!" =
            nrow(main_data_imputed) == nrow(main_data))

cat("\nFinal imputed dataset dimensions:", dim(main_data_imputed), "\n")


# ----------------------------------------
# SECTION 6: SAVE FEATURE GROUP DATASETS
# ----------------------------------------
# Save both the full imputed set and pre-split group datasets
# Group datasets will be used in script 05_modelling.R

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
cat("data/processed/main_data_imputed.csv\n")
cat("data/processed/group1_demo_lifestyle.csv\n")
cat("data/processed/group2_laboratory.csv\n")
cat("data/processed/group3_bone_scans.csv\n")
cat("Script 03 complete. Run 04 next.\n")
