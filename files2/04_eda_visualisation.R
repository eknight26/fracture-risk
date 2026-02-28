# =============================================================================
# FRACTURE RISK PREDICTION - NHANES 2017-2020
# Script 04: Exploratory Data Analysis & Visualisation
# Author: Ernest Caballero
# Description: Produces visualisations for each feature group to explore
#              distributions, relationships with fracture outcome, and
#              class imbalance. All plots saved as PNG to outputs/plots/.
# =============================================================================

rm(list = ls())
options(scipen = 999)

library(tidyverse)
library(patchwork)   # combine multiple ggplots into one figure
library(scales)      # axis formatting helpers


# =============================================================================
# SECTION 1: LOAD DATA
# =============================================================================

main_data <- read_csv("data/processed/main_data_imputed.csv",
                      show_col_types = FALSE) %>%
  mutate(hx_fracture = factor(hx_fracture,
                               levels = c("No_Fracture", "Fracture")))

dir.create("outputs/plots", showWarnings = FALSE, recursive = TRUE)

# Consistent colour palette across all plots
fracture_colours <- c("No_Fracture" = "steelblue4", "Fracture" = "firebrick3")

cat("Data loaded:", nrow(main_data), "respondents\n")


# =============================================================================
# SECTION 2: CLASS DISTRIBUTION (TARGET VARIABLE)
# =============================================================================

class_counts <- main_data %>%
  count(hx_fracture) %>%
  mutate(pct = round(n / sum(n) * 100, 1),
         label = paste0(n, "\n(", pct, "%)"))

p_class <- ggplot(class_counts, aes(x = hx_fracture, y = n, fill = hx_fracture)) +
  geom_col(width = 0.5, alpha = 0.9) +
  geom_text(aes(label = label), vjust = -0.3, size = 4, fontface = "bold") +
  scale_fill_manual(values = fracture_colours, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(
    title    = "Class Distribution: Fracture History",
    subtitle = "NHANES 2017-2020 Pre-Pandemic Cohort (Age 50+)",
    x        = NULL,
    y        = "Count",
    caption  = "Note: Class imbalance present — addressed in modelling via class weights"
  ) +
  theme_minimal(base_size = 13) +
  theme(axis.text = element_text(colour = "black"))

ggsave("outputs/plots/02_class_distribution.png", p_class,
       width = 6, height = 5, dpi = 300)
cat("Saved: 02_class_distribution.png\n")


# =============================================================================
# SECTION 3: GROUP 1 — DEMOGRAPHY & LIFESTYLE
# =============================================================================

# --- Age distribution by fracture status ---
p_age <- ggplot(main_data, aes(x = age, fill = hx_fracture)) +
  geom_histogram(binwidth = 5, position = "dodge", alpha = 0.85, colour = "white") +
  scale_fill_manual(values = fracture_colours, name = "Fracture History") +
  labs(title = "Age Distribution by Fracture Status",
       x = "Age (years)", y = "Count") +
  theme_minimal(base_size = 12)

# --- BMI boxplot by fracture status ---
p_bmi <- ggplot(main_data, aes(x = hx_fracture, y = bmi, fill = hx_fracture)) +
  geom_boxplot(alpha = 0.8, outlier.size = 0.8, outlier.alpha = 0.4) +
  scale_fill_manual(values = fracture_colours, guide = "none") +
  labs(title = "BMI by Fracture Status", x = NULL, y = "BMI (kg/m²)") +
  theme_minimal(base_size = 12)

# --- Gender breakdown ---
p_gender <- main_data %>%
  mutate(gender_label = if_else(gender == 1, "Male", "Female")) %>%
  count(gender_label, hx_fracture) %>%
  group_by(gender_label) %>%
  mutate(pct = n / sum(n) * 100) %>%
  ggplot(aes(x = gender_label, y = pct, fill = hx_fracture)) +
  geom_col(position = "dodge", alpha = 0.88) +
  scale_fill_manual(values = fracture_colours, name = "Fracture History") +
  labs(title = "Fracture Rate by Gender",
       x = NULL, y = "Proportion (%)") +
  theme_minimal(base_size = 12)

# --- Sedentary minutes ---
p_sedentary <- ggplot(main_data,
                       aes(x = hx_fracture, y = sedentary_mins_day,
                           fill = hx_fracture)) +
  geom_boxplot(alpha = 0.8, outlier.size = 0.8, outlier.alpha = 0.4) +
  scale_fill_manual(values = fracture_colours, guide = "none") +
  labs(title = "Sedentary Time by Fracture Status",
       x = NULL, y = "Minutes/Day") +
  theme_minimal(base_size = 12)

# Combine group 1 plots into one figure
group1_panel <- (p_age | p_bmi) / (p_gender | p_sedentary) +
  plot_annotation(
    title    = "Group 1: Demography & Lifestyle",
    subtitle = "NHANES 2017-2020 | Fracture Risk Cohort",
    caption  = "Source: CDC/NCHS NHANES",
    theme    = theme(plot.title = element_text(face = "bold", size = 15))
  )

ggsave("outputs/plots/03_group1_demo_lifestyle.png", group1_panel,
       width = 13, height = 10, dpi = 300)
cat("Saved: 03_group1_demo_lifestyle.png\n")


# =============================================================================
# SECTION 4: GROUP 2 — LABORATORY VALUES
# =============================================================================

# Helper function: boxplot for lab values — avoids repeating code
lab_boxplot <- function(var, label, units) {
  ggplot(main_data, aes(x = hx_fracture, y = .data[[var]], fill = hx_fracture)) +
    geom_boxplot(alpha = 0.8, outlier.size = 0.7, outlier.alpha = 0.3) +
    scale_fill_manual(values = fracture_colours, guide = "none") +
    labs(title = label, x = NULL, y = units) +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 15, hjust = 1))
}

p_alp       <- lab_boxplot("alp",            "Alkaline Phosphatase",  "U/L")
p_calcium   <- lab_boxplot("calcium",        "Serum Calcium",         "mmol/L")
p_phosphate <- lab_boxplot("phosphate",      "Serum Phosphate",       "mmol/L")
p_folate    <- lab_boxplot("serum_folate",   "Serum Folate",          "nmol/L")
p_lead      <- lab_boxplot("blood_lead",     "Blood Lead",            "µmol/L")
p_cadmium   <- lab_boxplot("blood_cadmium",  "Blood Cadmium",         "nmol/L")
p_mercury   <- lab_boxplot("blood_mercury",  "Blood Mercury",         "nmol/L")
p_selenium  <- lab_boxplot("blood_selenium", "Blood Selenium",        "µmol/L")

group2_panel <- (p_alp | p_calcium | p_phosphate | p_folate) /
                (p_lead | p_cadmium | p_mercury | p_selenium) +
  plot_annotation(
    title    = "Group 2: Laboratory Values by Fracture Status",
    subtitle = "Bone metabolism markers and environmental metal exposures",
    caption  = "Source: CDC/NCHS NHANES 2017-2020",
    theme    = theme(plot.title = element_text(face = "bold", size = 15))
  )

ggsave("outputs/plots/04_group2_laboratory.png", group2_panel,
       width = 16, height = 10, dpi = 300)
cat("Saved: 04_group2_laboratory.png\n")


# =============================================================================
# SECTION 5: GROUP 3 — BONE SCANS & T-SCORES
# =============================================================================

# BMD density plots — overlapping distributions by fracture status
bmd_density <- function(var, label) {
  ggplot(main_data, aes(x = .data[[var]], fill = hx_fracture)) +
    geom_density(alpha = 0.55) +
    scale_fill_manual(values = fracture_colours, name = "Fracture History") +
    labs(title = label, x = "BMD (g/cm²)", y = "Density") +
    theme_minimal(base_size = 11)
}

p_neck  <- bmd_density("bmd_femur_neck",  "Femur Neck BMD")
p_ward  <- bmd_density("bmd_femur_ward",  "Ward's Triangle BMD")
p_spine <- bmd_density("bmd_spine_total", "Spine Total BMD")
p_l1    <- bmd_density("bmd_l1",          "L1 Vertebra BMD")

# T-score distribution with WHO threshold lines
p_tscore_nk <- ggplot(main_data, aes(x = T_score_nk, fill = hx_fracture)) +
  geom_density(alpha = 0.55) +
  geom_vline(xintercept = -1.0, linetype = "dashed",
             colour = "darkorange", linewidth = 0.8) +
  geom_vline(xintercept = -2.5, linetype = "dashed",
             colour = "red3", linewidth = 0.8) +
  annotate("text", x = -1.0, y = Inf, label = "Osteopenia\n(-1.0)",
           vjust = 1.5, hjust = -0.1, size = 3, colour = "darkorange") +
  annotate("text", x = -2.5, y = Inf, label = "Osteoporosis\n(-2.5)",
           vjust = 1.5, hjust = -0.1, size = 3, colour = "red3") +
  scale_fill_manual(values = fracture_colours, name = "Fracture History") +
  labs(title = "Femur Neck T-Score Distribution",
       subtitle = "Dashed lines = WHO diagnostic thresholds",
       x = "T-Score", y = "Density") +
  theme_minimal(base_size = 11)

p_tscore_wd <- ggplot(main_data, aes(x = T_score_wd, fill = hx_fracture)) +
  geom_density(alpha = 0.55) +
  geom_vline(xintercept = -1.0, linetype = "dashed",
             colour = "darkorange", linewidth = 0.8) +
  geom_vline(xintercept = -2.5, linetype = "dashed",
             colour = "red3", linewidth = 0.8) +
  scale_fill_manual(values = fracture_colours, name = "Fracture History") +
  labs(title = "Ward's Triangle T-Score Distribution",
       x = "T-Score", y = "Density") +
  theme_minimal(base_size = 11)

group3_panel <- (p_neck | p_ward | p_spine | p_l1) /
                (p_tscore_nk | p_tscore_wd) +
  plot_annotation(
    title    = "Group 3: Bone Scan BMD & T-Scores by Fracture Status",
    subtitle = "DEXA femur and spine measurements | WHO T-score thresholds shown",
    caption  = "Source: CDC/NCHS NHANES 2017-2020",
    theme    = theme(plot.title = element_text(face = "bold", size = 15))
  )

ggsave("outputs/plots/05_group3_bone_scans.png", group3_panel,
       width = 16, height = 10, dpi = 300)
cat("Saved: 05_group3_bone_scans.png\n")


# =============================================================================
# SECTION 6: CORRELATION HEATMAP PER GROUP
# =============================================================================

plot_correlation <- function(data, vars, title, filename) {
  corr_matrix <- data %>%
    select(all_of(vars)) %>%
    cor(use = "pairwise.complete.obs") %>%
    as.data.frame() %>%
    rownames_to_column("var1") %>%
    pivot_longer(-var1, names_to = "var2", values_to = "correlation")

  ggplot(corr_matrix, aes(x = var1, y = var2, fill = correlation)) +
    geom_tile(colour = "white") +
    geom_text(aes(label = round(correlation, 2)), size = 2.5) +
    scale_fill_gradient2(low = "steelblue4", mid = "white", high = "firebrick3",
                         midpoint = 0, limits = c(-1, 1), name = "r") +
    labs(title = title, x = NULL, y = NULL) +
    theme_minimal(base_size = 10) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text   = element_text(colour = "black"))

  ggsave(filename, width = 12, height = 10, dpi = 300)
  cat("Saved:", filename, "\n")
}

# Continuous lab variables only for correlation (exclude binary)
lab_continuous <- c("alp", "phosphate", "calcium", "serum_folate",
                    "blood_lead", "blood_cadmium", "blood_mercury",
                    "blood_selenium", "blood_manganese", "ldl_level")

bone_continuous <- c("bmd_femur_total", "bmd_femur_neck", "bmd_femur_troch",
                     "bmd_femur_inter", "bmd_femur_ward",
                     "bmd_spine_total", "bmd_l1", "bmd_l2", "bmd_l3", "bmd_l4",
                     "T_score_nk", "T_score_wd")

plot_correlation(main_data, lab_continuous,
                 "Group 2: Laboratory Correlation Matrix",
                 "outputs/plots/06_corr_laboratory.png")

plot_correlation(main_data, bone_continuous,
                 "Group 3: Bone Scan Correlation Matrix",
                 "outputs/plots/07_corr_bone_scans.png")

cat("\nScript 04 complete. Run 05_modelling.R next.\n")
