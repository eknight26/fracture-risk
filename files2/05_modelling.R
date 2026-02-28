# =============================================================================
# FRACTURE RISK PREDICTION - NHANES 2017-2020
# Script 05: Modelling — Three Feature Groups × Three Algorithms
# Author: Ernest Caballero
# Description: Trains and evaluates logistic regression, random forest, and
#              XGBoost on each feature group separately. Class imbalance is
#              handled via class weights. Results saved for comparison in 06.
#
# Feature Groups:
#   Group 1 — Demography & Lifestyle
#   Group 2 — Laboratory (biochem, folate, heavy metals)
#   Group 3 — Bone Scans (DEXA BMD + T-scores)
#
# Models per group:
#   - Logistic Regression (interpretable baseline)
#   - Random Forest       (handles non-linearity, missing patterns)
#   - XGBoost             (gradient boosting, typically best performance)
# =============================================================================

rm(list = ls())
options(scipen = 999)

library(tidyverse)
library(caret)         # unified model training interface
library(randomForest)  # random forest
library(xgboost)       # gradient boosting
library(pROC)          # AUC, ROC curves
library(patchwork)     # combine plots


# =============================================================================
# SECTION 1: LOAD DATA & DEFINE FEATURE GROUPS
# =============================================================================

g1 <- read_csv("data/processed/group1_demo_lifestyle.csv", show_col_types = FALSE)
g2 <- read_csv("data/processed/group2_laboratory.csv",     show_col_types = FALSE)
g3 <- read_csv("data/processed/group3_bone_scans.csv",     show_col_types = FALSE)

# Feature columns per group (exclude IDs, weights, target)
exclude_cols <- c("SEQN", "SDMVPSU", "SDMVSTRA",
                  "samp_wt_intrvw", "samp_wt_mec", "folate_wt", "hx_fracture")

get_features <- function(df) setdiff(names(df), exclude_cols)

groups <- list(
  "Group1_Demo_Lifestyle" = g1,
  "Group2_Laboratory"     = g2,
  "Group3_Bone_Scans"     = g3
)

cat("Feature counts per group:\n")
for (name in names(groups)) {
  cat(" ", name, ":", length(get_features(groups[[name]])), "features\n")
}


# =============================================================================
# SECTION 2: TRAIN/TEST SPLIT
# =============================================================================
# Same seed and split applied to all groups to ensure fair comparison.
# Stratified split preserves class ratio in both sets.

set.seed(42)

# Use group 1 to define the split (same SEQN rows used for all groups)
split_index <- createDataPartition(g1$hx_fracture, p = 0.75, list = FALSE)

split_data <- function(df) {
  features <- get_features(df)
  train <- df[split_index,  c("hx_fracture", features)]
  test  <- df[-split_index, c("hx_fracture", features)]

  # Convert all character/factor predictors to numeric for XGBoost compatibility
  train <- train %>% mutate(across(where(is.character), as.numeric))
  test  <- test  %>% mutate(across(where(is.character), as.numeric))

  # Ensure target is factor with valid R level names
  train$hx_fracture <- factor(train$hx_fracture,
                               levels = c("No_Fracture", "Fracture"))
  test$hx_fracture  <- factor(test$hx_fracture,
                               levels = c("No_Fracture", "Fracture"))
  list(train = train, test = test)
}

split_groups <- lapply(groups, split_data)


# =============================================================================
# SECTION 3: TRAINING CONFIGURATION
# =============================================================================
# 5-fold cross-validation, repeated once — balances reliability and runtime.
# Class weights passed directly to caret for all models.

train_control <- trainControl(
  method          = "cv",
  number          = 5,
  classProbs      = TRUE,          # needed for AUC/ROC
  summaryFunction = twoClassSummary,
  savePredictions = "final"
)

# Calculate class weights to address ~85:15 imbalance
# Weight = total / (n_classes × class_count) — standard inverse frequency
compute_weights <- function(train_df) {
  n_total    <- nrow(train_df)
  n_no_fx    <- sum(train_df$hx_fracture == "No_Fracture")
  n_fx       <- sum(train_df$hx_fracture == "Fracture")
  c(
    "No_Fracture" = n_total / (2 * n_no_fx),
    "Fracture"    = n_total / (2 * n_fx)
  )
}


# =============================================================================
# SECTION 4: MODEL TRAINING FUNCTION
# =============================================================================

train_all_models <- function(group_name, split) {
  cat("\n", rep("=", 60), "\n", sep = "")
  cat("Training:", group_name, "\n")

  train_df <- split$train
  test_df  <- split$test
  features <- setdiff(names(train_df), "hx_fracture")
  weights  <- compute_weights(train_df)

  # Sample weights vector for caret
  sample_wts <- ifelse(train_df$hx_fracture == "Fracture",
                       weights["Fracture"], weights["No_Fracture"])

  results <- list()

  # ---- Logistic Regression ----
  cat("  [1/3] Logistic Regression...\n")
  set.seed(42)
  results$logistic <- train(
    hx_fracture ~ .,
    data        = train_df,
    method      = "glm",
    family      = "binomial",
    trControl   = train_control,
    metric      = "ROC",
    weights     = sample_wts
  )

  # ---- Random Forest ----
  cat("  [2/3] Random Forest...\n")
  rf_grid <- expand.grid(mtry = c(2, floor(sqrt(length(features))),
                                   floor(length(features) / 3)))
  set.seed(42)
  results$rf <- train(
    hx_fracture ~ .,
    data        = train_df,
    method      = "rf",
    trControl   = train_control,
    metric      = "ROC",
    tuneGrid    = rf_grid,
    weights     = sample_wts,
    ntree       = 500,
    importance  = TRUE
  )

  # ---- XGBoost ----
  cat("  [3/3] XGBoost...\n")
  xgb_grid <- expand.grid(
    nrounds          = c(100, 200),
    max_depth        = c(3, 6),
    eta              = 0.1,
    gamma            = 0,
    colsample_bytree = 0.8,
    min_child_weight = 1,
    subsample        = 0.8
  )
  set.seed(42)
  results$xgb <- train(
    hx_fracture ~ .,
    data        = train_df,
    method      = "xgbTree",
    trControl   = train_control,
    metric      = "ROC",
    tuneGrid    = xgb_grid,
    weights     = sample_wts,
    verbose     = 0
  )

  # Return trained models + test set for evaluation in section 5
  list(models = results, test = test_df, group = group_name)
}

# Train all groups — this is the main run (may take 5-15 mins)
cat("\nStarting model training across all groups...\n")
all_trained <- mapply(train_all_models,
                      names(split_groups), split_groups,
                      SIMPLIFY = FALSE)


# =============================================================================
# SECTION 5: EVALUATION FUNCTION
# =============================================================================

evaluate_model <- function(model, test_df, model_name, group_name) {
  pred_class <- predict(model, newdata = test_df)
  pred_prob  <- predict(model, newdata = test_df, type = "prob")[, "Fracture"]

  cm  <- confusionMatrix(pred_class, test_df$hx_fracture, positive = "Fracture")
  roc <- roc(test_df$hx_fracture, pred_prob,
             levels = c("No_Fracture", "Fracture"), quiet = TRUE)

  tibble(
    Group       = group_name,
    Model       = model_name,
    AUC         = round(auc(roc), 4),
    Accuracy    = round(cm$overall["Accuracy"], 4),
    Sensitivity = round(cm$byClass["Sensitivity"], 4),  # recall — critical for clinical
    Specificity = round(cm$byClass["Specificity"], 4),
    PPV         = round(cm$byClass["Pos Pred Value"], 4),
    F1          = round(cm$byClass["F1"], 4)
  )
}

# Evaluate all models across all groups
cat("\nEvaluating models on held-out test set...\n")

results_list <- list()
roc_list     <- list()

for (group_name in names(all_trained)) {
  trained  <- all_trained[[group_name]]
  test_df  <- trained$test
  models   <- trained$models

  model_labels <- c(logistic = "Logistic Regression",
                    rf       = "Random Forest",
                    xgb      = "XGBoost")

  for (m in names(models)) {
    row <- evaluate_model(models[[m]], test_df, model_labels[[m]], group_name)
    results_list[[paste(group_name, m)]] <- row

    # Store ROC for plotting
    pred_prob <- predict(models[[m]], newdata = test_df, type = "prob")[, "Fracture"]
    roc_obj   <- roc(test_df$hx_fracture, pred_prob,
                     levels = c("No_Fracture", "Fracture"), quiet = TRUE)
    roc_list[[paste(group_name, m)]] <- list(
      roc   = roc_obj,
      group = group_name,
      model = model_labels[[m]]
    )
  }
}

results_df <- bind_rows(results_list)
cat("\nAll model results:\n")
print(results_df, n = Inf)


# =============================================================================
# SECTION 6: FEATURE IMPORTANCE PLOTS
# =============================================================================

plot_importance <- function(trained_group, group_label, filename) {
  # Use random forest variable importance (most interpretable)
  rf_model  <- trained_group$models$rf$finalModel
  imp_df    <- importance(rf_model) %>%
    as.data.frame() %>%
    rownames_to_column("Feature") %>%
    arrange(desc(MeanDecreaseGini)) %>%
    head(15)   # top 15 features

  ggplot(imp_df, aes(x = reorder(Feature, MeanDecreaseGini),
                     y = MeanDecreaseGini)) +
    geom_col(fill = "steelblue4", alpha = 0.85) +
    coord_flip() +
    labs(
      title    = paste("Feature Importance —", group_label),
      subtitle = "Random Forest: Mean Decrease in Gini Impurity (top 15)",
      x        = NULL,
      y        = "Mean Decrease Gini"
    ) +
    theme_minimal(base_size = 12) +
    theme(axis.text = element_text(colour = "black"))

  ggsave(filename, width = 10, height = 7, dpi = 300)
  cat("Saved:", filename, "\n")
}

plot_importance(all_trained[["Group1_Demo_Lifestyle"]], "Group 1: Demo & Lifestyle",
                "outputs/plots/08_importance_group1.png")
plot_importance(all_trained[["Group2_Laboratory"]], "Group 2: Laboratory",
                "outputs/plots/09_importance_group2.png")
plot_importance(all_trained[["Group3_Bone_Scans"]], "Group 3: Bone Scans",
                "outputs/plots/10_importance_group3.png")


# =============================================================================
# SECTION 7: ROC CURVES PER GROUP
# =============================================================================

plot_roc_group <- function(group_name, roc_list, filename) {
  # Filter to current group
  group_rocs <- roc_list[grepl(group_name, names(roc_list))]

  colours <- c("Logistic Regression" = "steelblue4",
               "Random Forest"       = "darkorange",
               "XGBoost"             = "firebrick3")

  p <- ggplot()
  for (item in group_rocs) {
    roc_df <- data.frame(
      fpr   = 1 - item$roc$specificities,
      tpr   = item$roc$sensitivities,
      model = paste0(item$model, " (AUC=", round(auc(item$roc), 3), ")")
    )
    p <- p + geom_line(data = roc_df,
                       aes(x = fpr, y = tpr, colour = model), linewidth = 1)
  }

  p +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey60") +
    scale_colour_manual(values = setNames(
      c("steelblue4", "darkorange", "firebrick3"),
      sapply(group_rocs, function(x)
        paste0(x$model, " (AUC=", round(auc(x$roc), 3), ")"))
    )) +
    labs(
      title    = paste("ROC Curves —", gsub("_", " ", group_name)),
      subtitle = "Comparison of logistic regression, random forest, and XGBoost",
      x        = "False Positive Rate (1 - Specificity)",
      y        = "True Positive Rate (Sensitivity)",
      colour   = "Model"
    ) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "bottom")

  ggsave(filename, width = 8, height = 7, dpi = 300)
  cat("Saved:", filename, "\n")
}

plot_roc_group("Group1", roc_list, "outputs/plots/11_roc_group1.png")
plot_roc_group("Group2", roc_list, "outputs/plots/12_roc_group2.png")
plot_roc_group("Group3", roc_list, "outputs/plots/13_roc_group3.png")


# =============================================================================
# SECTION 8: SAVE RESULTS
# =============================================================================

write_csv(results_df, "outputs/model_results.csv")
saveRDS(all_trained, "outputs/all_trained_models.rds")  # save model objects

cat("\nSaved: outputs/model_results.csv\n")
cat("Saved: outputs/all_trained_models.rds\n")
cat("Script 05 complete. Run 06_model_comparison.R next.\n")
