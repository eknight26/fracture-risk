# -----------------------------------------------------------------------------
# SCRIPT 03: Missing Data Imputation
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

main_data_imp <- read_csv("data/processed/main_data_for_imputation.csv",
                          show_col_types = FALSE)

cat("Dimensions:", dim(main_data_imp), "\n")
cat("Total missing values:", sum(is.na(main_data_imp)), "\n")
print(head(main_data_imp))

# Correct datatypes
main_data_imp <- main_data_imp %>%
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
  "alc_drinks_per_day", "ever_smoked", "sedentary_mins_day"
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

# All features to impute
all_features <- c(group_demo_lifestyle, group_laboratory, group_bone_scans)

cat("\nFeature group sizes:\n")
cat("Demo & Lifestyle:", length(group_demo_lifestyle), "features\n")
cat("Laboratory:", length(group_laboratory), "features\n")
cat("Bone Scans:", length(group_bone_scans), "features\n")
cat("Total:", length(all_features), "features\n")

# Variables to exclude in imputation
# Identifiers, weights, target, and complete cols
exclude_from_imputation <- c("SEQN", "SDMVPSU", "SDMVSTRA", 
                             "samp_wt_intrvw", "samp_wt_mec", "folate_wt",
                             "hx_fracture")



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
                 "thyroid", "copd", "overweight", "asthma", "cancer",
                 "ever_smoked")

# Unordered categorical (3+ categories) — polytomous logistic regression
poly_vars <- c("diabetes_status", "osteo_class")

# Continuous numeric — predictive mean matching
pmm_vars <- c("income_ratio", "bmd_femur_neck",
              "bmd_femur_troch", "bmd_femur_inter", "bmd_femur_ward",
              "bmd_l1", "bmd_l2", "bmd_l3", "bmd_l4",
              "bmi", "mean_systolic", "mean_diastolic", "serum_folate",
              "alp", "phosphate", "calcium", "blood_lead", "blood_cadmium",
              "blood_mercury", "blood_selenium", "blood_manganese",
              "alc_drinks_per_day", "sedentary_mins_day", "T_score_nk")
# Note: ldl_level, a1c_level, bmd_femur_total, bmd_spine_total removed —
# either >80% missing or composite columns causing multicollinearity


# Missing data checks for each column 
# Age, gender, and ethnicity does not have missing values
na_check <- data_to_impute %>%
  summarise(across(everything(), ~ sum(is.na(.)))) %>%
  pivot_longer(everything(),
               names_to  = "feature",
               values_to = "n_missing") %>%
  mutate(pct_missing = round(n_missing / nrow(data_to_impute) * 100, 2)) %>%
  arrange(desc(n_missing))

cat("Columns with missing values:\n")
print(na_check %>% filter(n_missing > 0))

cat("\nSpot check factor levels:\n")
cat("Gender:\n"); print(levels(data_to_impute$gender))
cat("Had osteoporosis:\n"); print(levels(data_to_impute$had_osteoporosis))
cat("Diabetes status:\n"); print(levels(data_to_impute$diabetes_status))
cat("Osteo class:\n"); print(levels(data_to_impute$osteo_class))


# --- Columns to DROP from imputation (>80% missing and BMD composites) ---
# Too little observed data — imputed values would be unreliable
drop_missing_flag <- na_check %>%
  filter(pct_missing >= 80)

drop_from_imputation <- drop_missing_flag %>% pull(feature)
drop_from_bmd <- c("bmd_femur_total", "bmd_spine_total")   # highly correlated BMD composites

cat("\nVery high missingness columns dropped (>80% missing):\n")
print(drop_missing_flag %>% select(feature, n_missing, pct_missing))

# Remove from data_to_impute and variable lists
data_to_impute <- data_to_impute %>%
  select(-all_of(drop_from_imputation)) %>%
  select(-all_of(drop_from_bmd))

cat("Dropped (>80% missing):", drop_from_imputation, "\n")
cat("Dropped (BMD composites):", drop_from_bmd, "\n")
cat("Remaining columns:", ncol(data_to_impute), "\n")


# --- Update group to remove dropped columns ---
pmm_vars <- setdiff(pmm_vars, c(drop_from_imputation, drop_from_bmd))
group_laboratory <- setdiff(group_laboratory, drop_from_imputation)
group_bone_scans <- setdiff(group_bone_scans, drop_from_bmd)

# Rebuild all_features after dropping
all_features <- c(group_demo_lifestyle, group_laboratory, group_bone_scans)
cat("All features after dropping:", length(all_features), "\n")


# --- Columns to FLAG (50-80% missing) ---
# MICE will impute these but high uncertainty should be noted in methods section
high_missing_flag <- na_check %>%
  filter(pct_missing >= 50 & pct_missing < 80)

cat("\nHigh missingness columns (imputed with caution):\n")
print(high_missing_flag %>% select(feature, n_missing, pct_missing))


# --- Missing data per group after all column changes ---
cat("\nMissing data per group before imputation:\n")
cat("Demo & Lifestyle:", sum(is.na(select(data_to_impute, any_of(group_demo_lifestyle)))), "\n")
cat("Laboratory:", sum(is.na(select(data_to_impute, any_of(group_laboratory)))), "\n")
cat("Bone Scans:", sum(is.na(select(data_to_impute, any_of(group_bone_scans)))), "\n")


# --- Build method vector ---
method_vector <- rep("", ncol(data_to_impute))
names(method_vector) <- names(data_to_impute)

method_vector[names(method_vector) %in% pmm_vars] <- "pmm"
method_vector[names(method_vector) %in% binary_vars] <- "logreg"
method_vector[names(method_vector) %in% poly_vars] <- "polyreg"

cat("\nImputation methods assigned:\n")
print(table(method_vector))

# Data leakage check - verify target is not in data_to_impute
stopifnot("hx_fracture found in imputation data — data leakage risk!" =
            !"hx_fracture" %in% names(data_to_impute))
cat("Leakage check passed: hx_fracture not in data_to_impute.\n")


# --------------------------------
# SECTION 4: RUN MICE IMPUTATION
# --------------------------------
# m = 5     : 5 imputed datasets for inference (Rubin's rules) and ML
# maxit = 30: increased from 10 — needed for convergence with correlated BMD features
# seed = 42 : reproducibility
# quickpred(): custom predictor matrix to control which variables predict which
# mincor = 0.1  : keeps weak but clinically meaningful predictors
# (important for heavy metals, folate → bone density relationships)
# minpuc = 0.25 : predictor must have observed data for >= 25% of cases
# where the target is missing. Filters near-empty predictors, which is providing 
# very little information to the imputation model, adds computational cost and noise without meaningful benefit

pred_matrix <- quickpred(data_to_impute, mincor = 0.1, minpuc = 0.25)

# Inspect predictor counts — sanity check before running
# What you want to see:
# Most columns have 5-15 predictors (filtered down from ~48)
# Columns showing 0 predictors (MICE has no information to impute it from)
# Columns showing ~48 predictors (close to total columns), almost every column passed the filters for that variable, quickpred had no effect
# If ALL columns show ~48 predictors, quickpred filtered almost nothing overall, mincor too low
cat("\nNumber of predictors assigned per column:\n")
print(sort(rowSums(pred_matrix)))

# --- Manually force predictors for calcium and blood_selenium ---
# Both columns had 0 predictors assigned causing erratic density plots (maxit = 30)

# Calcium — clinically related to phosphate, ALP, and bone density
pred_matrix["calcium", "phosphate"]      <- 1
pred_matrix["calcium", "alp"]            <- 1
pred_matrix["calcium", "bmd_femur_neck"] <- 1
pred_matrix["calcium", "age"]            <- 1
pred_matrix["calcium", "bmi"]            <- 1

# Blood selenium level — related to other blood metals
pred_matrix["blood_selenium", "blood_lead"]      <- 1
pred_matrix["blood_selenium", "blood_cadmium"]   <- 1
pred_matrix["blood_selenium", "blood_mercury"]   <- 1
pred_matrix["blood_selenium", "blood_manganese"] <- 1
pred_matrix["blood_selenium", "age"]             <- 1

# Verify
cat("Predictors for calcium (", sum(pred_matrix["calcium", ]), "):\n")
print(names(pred_matrix["calcium", pred_matrix["calcium", ] == 1]))
cat("\nPredictors for blood_selenium (", sum(pred_matrix["blood_selenium", ]), "):\n")
print(names(pred_matrix["blood_selenium", pred_matrix["blood_selenium", ] == 1]))


# --- Run the MICE ---
cat("\nRunning MICE imputation — this may take a few minutes...\n")

imputed <- mice(
  data_to_impute,
  m               = 5,
  maxit           = 30,
  method          = method_vector,
  predictorMatrix = pred_matrix,
  seed            = 42,
  printFlag       = TRUE
)

cat("Imputation complete.\n")


# --------------------------------
# SECTION 5: CONVERGENCE CHECK
# --------------------------------
# Trace plots: mean and SD of imputed values across iterations
# Good convergence: lines mix well with no clear upward/downward trend
# Initial runs without quickpred() showed upward trends in BMD features - removed BMD composites
# quickpred() with mincor = 0.1 resolved convergence for all features
plot(imputed)


# --- Separate numeric and factor columns for plotting ---
numeric_cols <- names(data_to_impute)[sapply(data_to_impute, is.numeric)]
factor_cols  <- names(data_to_impute)[sapply(data_to_impute, is.factor)]

cat("Numeric columns:", length(numeric_cols), "\n")
cat("Factor columns: ", length(factor_cols),  "\n")


# --- Exclude columns that cannot be density plotted ---
# age, gender, ethnicity  — 0% missing, no imputed values to plot
# these columns stay in data_to_impute as predictors for other columns
# but densityplot() has nothing to show for them
skip_plot <- c("age", "gender", "ethnicity")
numeric_cols_plot <- setdiff(numeric_cols, skip_plot)
factor_cols_plot <- setdiff(factor_cols, skip_plot)

cat("Columns excluded from density plot (0% missing):", skip_plot, "\n")
cat("Columns included in density plot:", length(numeric_cols_plot), "\n")


# --- Density plots for numeric columns ---
# Blue line = observed data distribution
# Red lines = imputed values across m = 5 datasets
# Good: red lines follow similar shape to blue
# Bad:  red lines are completely different shape — misspecified imputation model

cat("\nPlotting density plots for numeric columns...\n")

for (col in numeric_cols_plot) {
  cat("Plotting:", col, "\n")
  tryCatch({
    print(
      densityplot(imputed,
                  as.formula(paste("~", col)),
                  main = paste("Density Plot —", col))
    )
  }, error = function(e) {
    cat("FAILED:", col, "—", conditionMessage(e), "\n")
  })
}

# --- Bar plots for factor columns ---
# Shows category proportions in observed vs imputed data
# Good: similar proportions between observed and imputed
# Bad:  large shifts in category proportions after imputation

cat("\nPlotting bar plots for factor columns...\n")

for (col in factor_cols_plot) {
  cat("Plotting:", col, "\n")
  tryCatch({
    observed <- data_to_impute[[col]]
    
    # Get imputed values across all 5 datasets
    imputed_vals <- lapply(1:imputed$m, function(i) {
      complete(imputed, action = i)[[col]]
    })
    
    # Build a dataframe comparing observed vs imputed proportions
    obs_table <- prop.table(table(observed)) %>%
      as.data.frame() %>%
      mutate(source = "Observed")
    
    imp_table <- lapply(seq_len(imputed$m), function(i) {
      prop.table(table(imputed_vals[[i]])) %>%
        as.data.frame() %>%
        mutate(source = paste("Imputed", i))
    }) %>% bind_rows()
    
    # Combine observed and imputed
    colnames(obs_table)[1] <- "category"
    colnames(imp_table)[1] <- "category"
    plot_data <- bind_rows(obs_table, imp_table)
    
    # Plot
    p <- ggplot(plot_data, aes(x = category, y = Freq, fill = source)) +
      geom_bar(stat = "identity", position = "dodge", alpha = 0.8) +
      scale_fill_manual(
        values = c("Observed" = "steelblue4",
                   setNames(rep("firebrick3", imputed$m),
                            paste("Imputed", seq_len(imputed$m))))
      ) +
      labs(
        title = paste("Imputed vs Observed —", col),
        subtitle = "Blue = Observed | Red = Imputed datasets 1-5",
        x = col,
        y = "Proportion",
        fill = "Source"
      ) +
      theme_minimal(base_size = 12) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    print(p)
    cat("OK\n")
    
  }, error = function(e) {
    cat("FAILED:", col, "—", conditionMessage(e), "\n")
  })
}

cat("\nImputation quality plots complete.\n")


# --- Extract completed datasets ---
# For ML: use a single completed dataset (action = 1) to train your model
# For inference: pool results across all 5 datasets using Rubin's rules —
# fit model on each dataset separately, then pool with mice::pool()

# Single dataset for ML
data_imputed_ml <- complete(imputed, action = 1)

# All 5 datasets stacked for inference pooling
data_imputed_all <- complete(imputed, action = "long", include = FALSE)

cat("Missing values after imputation (dataset 1):", sum(is.na(data_imputed_ml)), "\n")
cat("Rows in stacked dataset (5 x n):", nrow(data_imputed_all), "\n")


# ------------------------------------
# SECTION 5: REASSEMBLE FULL DATASET
# ------------------------------------

# Reattach excluded columns (SEQN, weights, target, etc) to imputed features
main_data_imputed <- bind_cols(data_excluded, data_imputed_ml)

# Sanity checks
# BUG FIX: was main_data — correct variable is main_data_imp
stopifnot("Row count changed after imputation!" =
            nrow(main_data_imputed) == nrow(main_data_imp))

stopifnot("Duplicate SEQNs in imputed dataset!" =
            length(unique(main_data_imputed$SEQN)) == nrow(main_data_imputed))

cat("Check passed: row count unchanged after imputation.\n")
cat("Check passed: no duplicate SEQNs.\n")
cat("\nFinal imputed dataset dimensions:", dim(main_data_imputed),"\n")
cat("Missing values in final dataset: ", sum(is.na(main_data_imputed)), "\n")


# ----------------------------------------------
# SECTION 7: SAVE MAIN & FEATURE GROUP DATASETS
# ----------------------------------------------
# Save full imputed dataset and pre-split group datasets
# Group datasets used in script 05_modelling.R

id_target_weights <- c("SEQN", "hx_fracture", "SDMVPSU", "SDMVSTRA",
                       "samp_wt_mec", "folate_wt", "samp_wt_intrvw")

# Full imputed dataset ML
write_csv(main_data_imputed, "data/processed/main_data_imputed.csv")
cat("\nSaved: data/processed/main_data_imputed.csv\n")

# Save a copy of un-imputed dataset
# Reattach excluded columns first
main_data_not_imputed <- bind_cols(data_excluded, data_imputed_ml)
write_csv(main_data_not_imputed, "data/processed/main_data_not_imputed.csv")
cat("\nSaved: data/processed/main_data_not_imputed.csv\n")


# Group 1: Demography & Lifestyle
main_data_imputed %>%
  select(any_of(id_target_weights), any_of(group_demo_lifestyle)) %>%
  write_csv("data/processed/group1_demo_lifestyle.csv")

# Group 2: Laboratory
main_data_imputed %>%
  select(any_of(id_target_weights), any_of(group_laboratory)) %>%
  write_csv("data/processed/group2_laboratory.csv")

# Group 3: Bone Scans
main_data_imputed %>%
  select(any_of(id_target_weights), any_of(group_bone_scans)) %>%
  write_csv("data/processed/group3_bone_scans.csv")

cat("\nSaved:\n")
cat("data/processed/main_data_imputed.csv\n")
cat("data/processed/group1_demo_lifestyle.csv\n")
cat("data/processed/group2_laboratory.csv\n")
cat("data/processed/group3_bone_scans.csv\n")
cat("\nScript 03 complete. Run script 04 next.\n")

