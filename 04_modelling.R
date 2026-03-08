# -----------------------------------------------------------------------------
# FRACTURE RISK PREDICTION - NHANES 2017-2020
# SCRIPT 04: Machine Learning Modelling
# AUTHOR: Ernest Caballero
# DESCRIPTION: Trains Logistic Regression, Random Forest, and XGBoost models
# on three feature group datasets and one combined dataset.
# Class imbalance handled with SMOTE. Hyperparameters tuned with 5-fold CV.
# Models evaluated on a held-out test set using ROC AUC, sensitivity,
# specificity, F1, and precision.
# -----------------------------------------------------------------------------

# LIBRARIES
library(tidyverse)
library(tidymodels)
library(themis)    # step_smote() for class imbalance
library(vip)       # variable importance plots
library(patchwork) # combine multiple ggplots

set.seed(42)

# "Fracture" is the second factor level — tell yardstick it is the positive class
# This ensures sens/spec/F1/precision are calculated correctly for fracture prediction
options(yardstick.event_first = FALSE)

# Output folders
dir.create("outputs/models", showWarnings = FALSE, recursive = TRUE)
dir.create("outputs/plots/models", showWarnings = FALSE, recursive = TRUE)


# ----------------------
# SECTION 1: LOAD DATA
# ----------------------

# Survey weights (samp_wt_mec, samp_wt_intrvw, folate_wt) and design variables
# (SDMVPSU, SDMVSTRA) are NOT used as model features.
# These columns are dropped before training below.

group1 <- read_csv("data/processed/group1_demo_lifestyle.csv", show_col_types = FALSE)
group2 <- read_csv("data/processed/group2_laboratory.csv", show_col_types = FALSE)
group3 <- read_csv("data/processed/group3_bone_scans.csv", show_col_types = FALSE)
combined <- read_csv("data/processed/main_data_imputed.csv", show_col_types = FALSE)

cat("Datasets loaded:\n")
cat("Group 1 (Demo & Lifestyle):", nrow(group1), "rows,", ncol(group1), "cols\n")
cat("Group 2 (Laboratory):", nrow(group2), "rows,", ncol(group2), "cols\n")
cat("Group 3 (Bone Scans):", nrow(group3), "rows,", ncol(group3), "cols\n")
cat("Combined:", nrow(combined), "rows,", ncol(combined), "cols\n")


# ------------------------------------
# SECTION 2: DATA PREPARATION HELPER
# ------------------------------------
# prepare_data() restores the correct types and drops non-feature columns

# Columns that are never model features
drop_cols <- c("SEQN", "SDMVPSU", "SDMVSTRA",
               "samp_wt_intrvw", "samp_wt_mec", "folate_wt")

# Binary features encoded as 0/1 integers in script 03
binary_vars <- c("had_osteoporosis", "other_fractures", "prednisone",
                 "mother_fx_hip", "father_fx_hip", "high_bp", "high_choles",
                 "arthritis", "chf", "chd", "heart_attack", "stroke", "liver",
                 "thyroid", "copd", "overweight", "asthma", "cancer", 
                 "ever_smoked")

prepare_data <- function(data) {
  data %>%
    select(-any_of(drop_cols)) %>%
    mutate(
      hx_fracture = factor(hx_fracture,
                           levels = c("No_Fracture", "Fracture")),
      
      gender    = factor(gender,
                         levels = c("Male", "Female")),
      
      ethnicity = factor(ethnicity,
                         levels = c("Mexican American", "Other Hispanic",
                                    "Non-Hispanic White", "Non-Hispanic Black",
                                    "Non-Hispanic Asian", "Other/Multiracial")),
      
      diabetes_status = factor(diabetes_status,
                               levels = c("Yes", "No", "Borderline"))
    ) %>%
    mutate(across(
      any_of(binary_vars), ~ factor(.x, levels = c("No", "Yes"))
    )) %>%
    { if ("osteo_class" %in% names(.))
      mutate(., osteo_class = factor(osteo_class,
                                     levels = c("Normal", "Osteopenia", "Osteoporosis")))
      else . }
}

# Apply preparation to all datasets
group1 <- prepare_data(group1)
group2 <- prepare_data(group2)
group3 <- prepare_data(group3)
combined <- prepare_data(combined)

cat("\nClass distribution (combined dataset):\n")
print(table(combined$hx_fracture))

cat("Imbalance ratio:", round(sum(combined$hx_fracture == "No_Fracture") / sum(combined$hx_fracture == "Fracture"), 1), ":1\n")



# -----------------------------------------------
# SECTION 3: MODELLING FUNCTION
# -----------------------------------------------
# train_and_evaluate() runs the full ML pipeline for one dataset:
#   1. Stratified 80/20 train/test split
#   2. 5-fold cross-validation for hyperparameter tuning
#   3. Preprocessing (dummy encoding, normalisation, SMOTE)
#   4. Tune and fit: Logistic Regression, Random Forest, XGBoost
#   5. Evaluate on held-out test set
#   6. Save ROC curve and variable importance plots
# Returns a list of metrics, plots, and fitted model objects.

train_and_evaluate <- function(data, group_name) {
  
  cat("\n", strrep("=", 60), "\n", sep = "")
  cat("GROUP:", group_name, "\n")
  cat(strrep("=", 60), "\n", sep = "")
  cat("Features:", ncol(data) - 1, "| Rows:", nrow(data), "\n")
  
  # --- 1. Train / Test Split ---
  # Stratified on hx_fracture: preserves the ~85/15 class ratio in both sets
  split <- initial_split(data, prop = 0.8, strata = hx_fracture)
  train <- training(split)
  test <- testing(split)
  
  cat("Train:", nrow(train), "| Test:", nrow(test), "\n")
  
  
  # --- 2. Cross-Validation Folds ---
  cv_folds <- vfold_cv(train, v = 5, strata = hx_fracture)
  
  
  # --- 3. Preprocessing Recipe ---
  # All steps applied inside CV folds to prevent data leakage
  # SMOTE is applied to training folds only, never on validation fold
  recipe <- recipe(hx_fracture ~ ., data = train) %>%
    
    # Convert nominal predictors to dummy variables (0/1 columns); required for LR and XGBoost
    step_dummy(all_nominal_predictors()) %>%
    
    # Standardise numeric predictors
    step_normalize(all_numeric_predictors()) %>%
    
    # SMOTE: creates synthetic fracture cases to balance the ~85/15 class split
    # over_ratio = 0.5 (minority class grows to 50% the size of majority class)
    # Applied after normalisation so synthetic points are on the same scale
    step_smote(hx_fracture, over_ratio = 0.5, seed = 42)
  
  
  # --- 4. Define Models ---
  # Logistic Regression — LASSO penalised (mixture = 1)
  lr_model <- logistic_reg(penalty = tune(), mixture = 1) %>%
    set_engine("glmnet") %>%
    set_mode("classification")
  
  # Random Forest
  # trees: 500 fixed — enough for stable predictions without excessive runtime
  rf_model <- rand_forest(mtry = tune(), min_n = tune(), trees = 500) %>%
    set_engine("ranger", importance = "impurity") %>%
    set_mode("classification")
  
  # XGBoost (Gradient Boosted Trees)
  # trees: number of boosting rounds
  # learn_rate:shrinkage per round — lower = slow but robust
  # tree_depth: maximum depth per tree — controls model complexity
  # min_n: minimum samples to split a node
  xgb_model <- boost_tree(trees = tune(),
                          learn_rate = tune(),
                          tree_depth = tune(),
                          min_n = tune()
                          ) %>%
    set_engine("xgboost") %>%
    set_mode("classification")
  
  
  # --- 5. Build Workflows ---
  # A workflow bundles the preprocessing steps + model spec together
  lr_wf  <- workflow() %>% add_recipe(recipe) %>% add_model(lr_model)
  rf_wf  <- workflow() %>% add_recipe(recipe) %>% add_model(rf_model)
  xgb_wf <- workflow() %>% add_recipe(recipe) %>% add_model(xgb_model)
  
  
  # --- 6. Tune Hyperparameters ---
  # tune_grid() evaluates random hyperparameter combinations on CV folds
  # grid = 10: 10 random combinations
  # Primary metric: ROC AUC
  
  num_grid <- 10
  tune_metrics <- metric_set(roc_auc, f_meas, sens, spec, precision)
  
  cat("\nTuning Logistic Regression...\n")
  lr_tune <- tune_grid(lr_wf,
                       resamples = cv_folds,
                       grid = num_grid,
                       metrics = tune_metrics
                      )
  
  cat("Tuning Random Forest...\n")
  rf_tune <- tune_grid(rf_wf,
                       resamples = cv_folds,
                       grid = num_grid,
                       metrics = tune_metrics
                       )
  
  cat("Tuning XGBoost...\n")
  xgb_tune <- tune_grid(xgb_wf, 
                        resamples = cv_folds,
                        grid = num_grid,
                        metrics = tune_metrics
                        )
  
  
  # --- 7. Select Best Hyperparameters ---
  best_lr  <- select_best(lr_tune, metric = "roc_auc")
  best_rf  <- select_best(rf_tune, metric = "roc_auc")
  best_xgb <- select_best(xgb_tune, metric = "roc_auc")
  
  cat("\nBest hyperparameters (by ROC AUC):\n")
  cat("  LR  — penalty:", round(best_lr$penalty, 5), "\n")
  cat("  RF  — mtry:", best_rf$mtry, "| min_n:", best_rf$min_n, "\n")
  cat("  XGB — trees:", best_xgb$trees,
      "| learn_rate:", round(best_xgb$learn_rate, 4),
      "| depth:", best_xgb$tree_depth, "\n")
  
  
  # --- 8. Finalise and Fit on Full Training Set ---
  # final_workflow() replaces tune() placeholders with the best values
  # last_fit() trains on the full training set and evaluates on the test set
  final_lr_wf <- final_workflow(lr_wf,  best_lr)
  final_rf_wf <- final_workflow(rf_wf,  best_rf)
  final_xgb_wf <- final_workflow(xgb_wf, best_xgb)
  
  lr_fit <- last_fit(final_lr_wf, split, metrics = tune_metrics)
  rf_fit <- last_fit(final_rf_wf, split, metrics = tune_metrics)
  xgb_fit <- last_fit(final_xgb_wf, split, metrics = tune_metrics)
  
  
  # --- 9. Collect Test Set Metrics ---
  add_labels <- function(fit, model_name) {
    collect_metrics(fit) %>%
      mutate(model = model_name, group = group_name)
  }
  
  all_metrics <- bind_rows(
    add_labels(lr_fit, "Logistic Regression"),
    add_labels(rf_fit, "Random Forest"),
    add_labels(xgb_fit, "XGBoost")
  )
  
  cat("\nTest set results:\n")
  all_metrics %>%
    select(model, .metric, .estimate) %>%
    mutate(.estimate = round(.estimate, 3)) %>%
    arrange(model, .metric) %>%
    print()
  
  
  # --- 10. ROC Curve Plot ---
  roc_data <- bind_rows(
    collect_predictions(lr_fit)  %>% mutate(model = "Logistic Regression"),
    collect_predictions(rf_fit)  %>% mutate(model = "Random Forest"),
    collect_predictions(xgb_fit) %>% mutate(model = "XGBoost")
  )
  
  # Compute AUC labels for the legend
  auc_labels <- all_metrics %>%
    filter(.metric == "roc_auc") %>%
    mutate(label = paste0(model, " (AUC = ", round(.estimate, 3), ")")) %>%
    select(model, label)
  
  roc_plot_data <- roc_data %>%
    group_by(model) %>%
    roc_curve(truth = hx_fracture, .pred_Fracture) %>%
    left_join(auc_labels, by = "model")
  
  roc_plot <- ggplot(roc_plot_data,
                     aes(x = 1 - specificity, y = sensitivity, colour = label)) +
    geom_line(linewidth = 1) +
    geom_abline(linetype = "dashed", colour = "grey60") +
    labs(
      title    = paste("ROC Curves —", group_name),
      subtitle = "Evaluated on held-out test set (20%)",
      x        = "False Positive Rate (1 - Specificity)",
      y        = "True Positive Rate (Sensitivity)",
      colour   = NULL,
      caption  = "Dashed line = random classifier (AUC = 0.5)"
    ) +
    scale_colour_brewer(palette = "Set1") +
    theme_minimal(base_size = 12) +
    theme(legend.position = "bottom")
  
  safe_name <- gsub("[^A-Za-z0-9]", "_", group_name)
  ggsave(
    paste0("outputs/plots/models/roc_", safe_name, ".png"),
    roc_plot, width = 8, height = 6, dpi = 300, bg = "white"
  )
  cat("Saved ROC plot:", safe_name, "\n")
  
  
  # --- 11. Variable Importance Plots ---
  # RF: Gini impurity — how much each feature reduces node impurity on average
  # XGB: Gain — how much each feature improves loss when used in a split
  # LR: Absolute coefficient magnitude after LASSO penalisation
  
  rf_vip <- extract_fit_parsnip(rf_fit) %>%
    vip(num_features = 15, aesthetics = list(fill = "steelblue4", alpha = 0.85)) +
    labs(title = "Random Forest", subtitle = "Gini impurity importance") +
    theme_minimal(base_size = 11)
  
  xgb_vip <- extract_fit_parsnip(xgb_fit) %>%
    vip(num_features = 15, aesthetics = list(fill = "firebrick3", alpha = 0.85)) +
    labs(title = "XGBoost", subtitle = "Gain importance") +
    theme_minimal(base_size = 11)
  
  lr_vip <- extract_fit_parsnip(lr_fit) %>%
    vip(num_features = 15, aesthetics = list(fill = "darkgreen", alpha = 0.85)) +
    labs(title = "Logistic Regression", subtitle = "Absolute LASSO coefficient") +
    theme_minimal(base_size = 11)
  
  vip_panel <- lr_vip + rf_vip + xgb_vip +
    plot_annotation(
      title = paste("Variable Importance —", group_name),
      caption = "Top 15 features per model"
    )
  
  ggsave(
    paste0("outputs/plots/models/vip_", safe_name, ".png"),
    vip_panel, width = 14, height = 6, dpi = 300, bg = "white"
  )
  cat("Saved VIP plot:", safe_name, "\n")
  
  
  # Return everything — used for cross-group comparison
  list(
    group = group_name,
    metrics = all_metrics,
    roc_plot = roc_plot,
    lr_fit = lr_fit,
    rf_fit = rf_fit,
    xgb_fit = xgb_fit,
    best_lr = best_lr,
    best_rf = best_rf,
    best_xgb = best_xgb
  )
}


# -----------------------------------------------
# SECTION 4: RUN MODELS FOR ALL GROUPS
# -----------------------------------------------
# imap() passes both the dataset and its name from the named list

datasets <- list(
  "Group 1: Demographics & Lifestyle" = group1,
  "Group 2: Laboratory"               = group2,
  "Group 3: Bone Scans"               = group3,
  "Combined"                          = combined
)

cat("\n\nStarting model training...\n")
results <- imap(datasets, ~ train_and_evaluate(.x, .y))
cat("\nAll models trained.\n")


# -----------------------------------------------
# SECTION 5: COMPARE RESULTS ACROSS GROUPS
# -----------------------------------------------

# Combine all metrics from all groups into one dataframe
all_metrics <- map_dfr(results, ~ .x$metrics)

# ROC AUC summary table — primary comparison metric
auc_summary <- all_metrics %>%
  filter(.metric == "roc_auc") %>%
  select(Group = group, Model = model, ROC_AUC = .estimate) %>%
  mutate(ROC_AUC = round(ROC_AUC, 3)) %>%
  arrange(desc(ROC_AUC))

cat("\n=== ROC AUC Summary — All Groups & Models ===\n")
print(auc_summary, n = Inf)

# Full metrics summary table
metrics_summary <- all_metrics %>%
  mutate(.estimate = round(.estimate, 3)) %>%
  select(Group = group, Model = model, Metric = .metric, Score = .estimate) %>%
  pivot_wider(names_from = Metric, values_from = Score) %>%
  arrange(desc(roc_auc))

cat("\n=== Full Metrics Summary ===\n")
print(metrics_summary, n = Inf)


# --- Comparison Plot: All metrics × All groups × All models ---
comparison_plot <- all_metrics %>%
  filter(.metric %in% c("roc_auc", "f_meas", "sens", "spec")) %>%
  mutate(
    .metric = recode(.metric,
                     "roc_auc" = "ROC AUC",
                     "f_meas"  = "F1 Score",
                     "sens"    = "Sensitivity",
                     "spec"    = "Specificity"),
    # Shorten group names for plot labels
    group = recode(group,
                   "Group 1: Demographics & Lifestyle" = "Group 1\nDemo & Lifestyle",
                   "Group 2: Laboratory"               = "Group 2\nLaboratory",
                   "Group 3: Bone Scans"               = "Group 3\nBone Scans",
                   "Combined"                          = "Combined")
  ) %>%
  ggplot(aes(x = model, y = .estimate, fill = model)) +
  geom_col(alpha = 0.85, width = 0.7) +
  geom_text(aes(label = round(.estimate, 2)),
            vjust = -0.4, size = 2.8, colour = "black") +
  facet_grid(.metric ~ group) +
  scale_fill_brewer(palette = "Set1") +
  scale_y_continuous(limits = c(0, 1.05), expand = expansion(mult = c(0, 0.05))) +
  labs(
    title    = "Model Performance by Feature Group",
    subtitle = "Evaluated on held-out test set (20%) | Positive class = Fracture",
    x        = NULL,
    y        = "Score",
    fill     = "Model",
    caption  = "NHANES 2017-2020 | Class imbalance addressed with SMOTE (over_ratio = 0.5)"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    strip.text   = element_text(size = 9, face = "bold"),
    legend.position = "bottom"
  )

ggsave("outputs/plots/models/comparison_all_groups.png",
       comparison_plot, width = 14, height = 10, dpi = 300, bg = "white")
cat("\nSaved: outputs/plots/models/comparison_all_groups.png\n")


# --- Combined ROC curve: best model from each group side by side ---
# Identify best model per group by ROC AUC
best_per_group <- auc_summary %>%
  group_by(Group) %>%
  slice_max(ROC_AUC, n = 1) %>%
  ungroup()

cat("\nBest model per group:\n")
print(best_per_group)


# -----------------------------------------------
# SECTION 6: SAVE RESULTS
# -----------------------------------------------

write_csv(all_metrics, "outputs/model_metrics.csv")
write_csv(metrics_summary, "outputs/model_metrics_summary.csv")

cat("\nSaved: outputs/model_metrics.csv\n")
cat("Saved: outputs/model_metrics_summary.csv\n")

# Save all model result objects (for downstream use in script 05)
saveRDS(results, "outputs/models/all_model_results.rds")
cat("Saved: outputs/models/all_model_results.rds\n")

cat("\nScript 04 complete. Run script 05 next.\n")
