# =============================================================================
# FRACTURE RISK PREDICTION - NHANES 2017-2020
# Script 06: Model Comparison — Feature Groups × Algorithms
# Author: Ernest Caballero
# Description: Loads results from script 05 and produces a formatted comparison
#              table and summary visualisation across all 9 models
#              (3 groups × 3 algorithms). Identifies the best-performing
#              group and model for fracture prediction.
# =============================================================================

rm(list = ls())
options(scipen = 999)

library(tidyverse)
library(patchwork)
library(scales)


# =============================================================================
# SECTION 1: LOAD RESULTS
# =============================================================================

results_df <- read_csv("outputs/model_results.csv", show_col_types = FALSE)

# Clean group labels for display
results_df <- results_df %>%
  mutate(
    Group_Label = case_when(
      grepl("Group1", Group) ~ "Group 1\nDemo & Lifestyle",
      grepl("Group2", Group) ~ "Group 2\nLaboratory",
      grepl("Group3", Group) ~ "Group 3\nBone Scans",
      TRUE ~ Group
    )
  )

cat("Results loaded:", nrow(results_df), "models evaluated\n")
print(results_df %>% select(-Group_Label), n = Inf)


# =============================================================================
# SECTION 2: FORMATTED COMPARISON TABLE
# =============================================================================
# A clean, readable summary of all 9 models across all metrics.
# AUC and Sensitivity are most clinically relevant for fracture screening.

comparison_table <- results_df %>%
  select(Group_Label, Model, AUC, Sensitivity, Specificity,
         Accuracy, PPV, F1) %>%
  rename(
    `Feature Group`  = Group_Label,
    `AUC`            = AUC,
    `Sensitivity`    = Sensitivity,
    `Specificity`    = Specificity,
    `Accuracy`       = Accuracy,
    `PPV`            = PPV,
    `F1 Score`       = F1
  ) %>%
  arrange(`Feature Group`, desc(AUC))

cat("\n=== MODEL COMPARISON TABLE ===\n")
print(comparison_table, n = Inf)

write_csv(comparison_table, "outputs/model_comparison_table.csv")
cat("Saved: outputs/model_comparison_table.csv\n")


# =============================================================================
# SECTION 3: BEST MODEL IDENTIFICATION
# =============================================================================

best_overall <- results_df %>%
  slice_max(AUC, n = 1)

best_per_group <- results_df %>%
  group_by(Group_Label) %>%
  slice_max(AUC, n = 1) %>%
  ungroup()

cat("\n=== BEST MODEL OVERALL ===\n")
print(best_overall %>% select(Group, Model, AUC, Sensitivity, Specificity))

cat("\n=== BEST MODEL PER FEATURE GROUP ===\n")
print(best_per_group %>% select(Group_Label, Model, AUC, Sensitivity, F1))


# =============================================================================
# SECTION 4: VISUALISATION — AUC COMPARISON ACROSS GROUPS & MODELS
# =============================================================================

model_colours <- c(
  "Logistic Regression" = "steelblue4",
  "Random Forest"       = "darkorange",
  "XGBoost"             = "firebrick3"
)

# Grouped bar chart: AUC by feature group, coloured by model
p_auc <- ggplot(results_df,
                aes(x = Group_Label, y = AUC, fill = Model)) +
  geom_col(position = position_dodge(width = 0.75),
           width = 0.65, alpha = 0.88) +
  geom_text(aes(label = round(AUC, 3)),
            position = position_dodge(width = 0.75),
            vjust = -0.4, size = 3.2, fontface = "bold") +
  scale_fill_manual(values = model_colours) +
  scale_y_continuous(limits = c(0, 1),
                     breaks = seq(0, 1, 0.1),
                     labels = number_format(accuracy = 0.1)) +
  geom_hline(yintercept = 0.7, linetype = "dashed",
             colour = "grey50", linewidth = 0.6) +
  annotate("text", x = 0.55, y = 0.72, label = "AUC = 0.70 (acceptable threshold)",
           size = 3, colour = "grey40") +
  labs(
    title    = "Model Comparison: AUC by Feature Group",
    subtitle = "Fracture risk prediction — NHANES 2017-2020",
    x        = NULL,
    y        = "Area Under ROC Curve (AUC)",
    fill     = "Model",
    caption  = "Higher AUC = better discrimination | Dashed line = acceptable performance threshold"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position   = "bottom",
    axis.text         = element_text(colour = "black"),
    plot.title        = element_text(face = "bold")
  )


# =============================================================================
# SECTION 5: SENSITIVITY VS SPECIFICITY TRADEOFF
# =============================================================================
# Sensitivity (recall) = critical in clinical screening — we want to CATCH fractures
# This plot shows the tradeoff across all models

p_sens_spec <- results_df %>%
  select(Group_Label, Model, Sensitivity, Specificity) %>%
  pivot_longer(c(Sensitivity, Specificity), names_to = "Metric", values_to = "Value") %>%
  ggplot(aes(x = Group_Label, y = Value, fill = Model)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.65, alpha = 0.88) +
  scale_fill_manual(values = model_colours) +
  scale_y_continuous(limits = c(0, 1), labels = number_format(accuracy = 0.01)) +
  facet_wrap(~Metric) +
  labs(
    title    = "Sensitivity vs Specificity by Feature Group",
    subtitle = "Sensitivity = clinical priority: correctly identifying fracture cases",
    x        = NULL,
    y        = "Value",
    fill     = "Model"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    axis.text.x     = element_text(colour = "black", size = 9),
    strip.text      = element_text(face = "bold", size = 11)
  )


# =============================================================================
# SECTION 6: F1 SCORE COMPARISON (handles class imbalance better than accuracy)
# =============================================================================

p_f1 <- ggplot(results_df, aes(x = reorder(paste(Group_Label, Model), F1),
                                y = F1, fill = Model)) +
  geom_col(alpha = 0.88) +
  geom_text(aes(label = round(F1, 3)), hjust = -0.15, size = 3.2) +
  coord_flip() +
  scale_fill_manual(values = model_colours) +
  scale_y_continuous(limits = c(0, 1),
                     expand = expansion(mult = c(0, 0.12))) +
  labs(
    title    = "F1 Score — All Models",
    subtitle = "F1 balances precision and recall; robust to class imbalance",
    x        = NULL,
    y        = "F1 Score",
    fill     = "Model"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "none",
    axis.text       = element_text(colour = "black")
  )


# =============================================================================
# SECTION 7: SAVE COMBINED COMPARISON FIGURE
# =============================================================================

combined_comparison <- p_auc / (p_sens_spec | p_f1) +
  plot_annotation(
    title   = "Fracture Risk Prediction: Feature Group & Model Comparison",
    caption = "NHANES 2017-2020 Pre-Pandemic | Ernest Caballero",
    theme   = theme(
      plot.title   = element_text(face = "bold", size = 16),
      plot.caption = element_text(colour = "grey50")
    )
  )

ggsave("outputs/plots/14_model_comparison.png", combined_comparison,
       width = 16, height = 14, dpi = 300)
cat("Saved: outputs/plots/14_model_comparison.png\n")

ggsave("outputs/plots/15_auc_comparison.png", p_auc,
       width = 10, height = 7, dpi = 300)
cat("Saved: outputs/plots/15_auc_comparison.png\n")


# =============================================================================
# SECTION 8: PLAIN-LANGUAGE SUMMARY FOR README / REPORT
# =============================================================================

cat("\n", rep("=", 60), "\n", sep = "")
cat("SUMMARY\n")
cat(rep("=", 60), "\n", sep = "")
cat("Best overall model:\n")
cat("  Group:", best_overall$Group, "\n")
cat("  Algorithm:", best_overall$Model, "\n")
cat("  AUC:", best_overall$AUC, "\n")
cat("  Sensitivity:", best_overall$Sensitivity, "\n")
cat("  F1:", best_overall$F1, "\n")
cat("\nClinical interpretation:\n")
cat("  Sensitivity is the priority metric — a fracture screening tool must\n")
cat("  minimise false negatives (missed fractures) even at the cost of\n")
cat("  some false positives (unnecessary follow-up).\n")
cat(rep("=", 60), "\n", sep = "")

cat("\nScript 06 complete. Pipeline finished.\n")
cat("Check outputs/ for plots and outputs/model_comparison_table.csv\n")
