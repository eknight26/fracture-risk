# =============================================================================
# FRACTURE RISK PREDICTION - NHANES 2017-2020
# Script 01: Data Preparation & Merging (CONTINUATION)
# Author: Ernest Caballero
# Description: Merges all loaded modules, quality checks, missing data summary,
#              and saves analysis-ready dataset. Run after all module loads above.
# =============================================================================

# =============================================================================
# SECTION 2: MERGE ALL MODULES
# =============================================================================

# Base: demographics inner joined with osteoporosis (ensures only respondents
# with both demographic and fracture data are retained)
main_data <- demo %>%
  inner_join(osteo, by = "SEQN")

# All supplementary datasets — listed explicitly for transparency
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


# =============================================================================
# SECTION 3: RECODE SPECIAL VALUES TO NA
# =============================================================================
# NHANES encodes refused (7/77/777) and don't know (9/99/999/9999) as numeric.
# These must be set to NA before imputation — otherwise they corrupt the data.

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

# Continuous variables with special codes
main_data <- main_data %>%
  mutate(
    sedentary_mins_day = na_if(sedentary_mins_day, 9999),  # don't know
    sedentary_mins_day = na_if(sedentary_mins_day, 7777),  # refused
    a1c_level          = na_if(a1c_level, 9999),
    a1c_level          = na_if(a1c_level, 777),
    alc_drinks_per_day = na_if(alc_drinks_per_day, 999),
    alc_drinks_per_day = na_if(alc_drinks_per_day, 777),
    cigarettes_per_day = na_if(cigarettes_per_day, 999)
  )


# =============================================================================
# SECTION 4: OUTLIER DETECTION — NUMERICAL FEATURES
# =============================================================================
# Outlier detection runs AFTER special value recoding (Section 3) so that
# NHANES codes like 9999 are already NA and won't be flagged as extreme values.
# Runs BEFORE imputation so outliers don't bias the imputed values.
#
# Method: Interquartile Range (IQR) rule
#   Lower fence = Q1 - 3 * IQR   (using 3× for clinical data — more conservative
#   Upper fence = Q3 + 3 * IQR    than the standard 1.5×, reduces false flags)
#
# Decision rule: Winsorize (cap) rather than delete outliers.
#   Deletion loses already-scarce fracture cases (minority class).
#   Winsorizing replaces extreme values with the fence value, preserving the row.
#
# Biological implausibles (e.g. BMI of 1 or 200) are genuinely erroneous.
# Extreme-but-plausible values (e.g. very high lead in an exposed individual)
# are real signal — winsorizing caps without deleting them.
# =============================================================================

# Define all continuous numerical features to check
# Excludes: binary variables, IDs, weights, target, and ordinal categoricals
continuous_vars <- c(
  # Demographics / body
  "age", "income_ratio", "bmi", "mean_systolic", "mean_diastolic",
  # Lifestyle
  "alc_drinks_per_day", "cigarettes_per_day", "sedentary_mins_day",
  # Glycaemic
  "a1c_level", "ldl_level",
  # Bone metabolism bloods
  "alp", "phosphate", "calcium", "serum_folate",
  # Heavy metals
  "blood_lead", "blood_cadmium", "blood_mercury",
  "blood_selenium", "blood_manganese",
  # DEXA femur BMD
  "bmd_femur_total", "bmd_femur_neck", "bmd_femur_troch",
  "bmd_femur_inter", "bmd_femur_ward",
  # DEXA spine BMD
  "bmd_spine_total", "bmd_l1", "bmd_l2", "bmd_l3", "bmd_l4"
)

# Filter to only columns that actually exist in main_data at this point
# (T-score columns added later in script 01 Section 6 after t-score merge)
continuous_vars <- intersect(continuous_vars, names(main_data))

# --- Step 4a: Compute outlier summary BEFORE winsorizing ---
outlier_summary <- map_dfr(continuous_vars, function(var) {
  x    <- main_data[[var]]
  x    <- x[!is.na(x)]                        # exclude NA — not outliers
  q1   <- quantile(x, 0.25)
  q3   <- quantile(x, 0.75)
  iqr  <- q3 - q1
  lo   <- q1 - 3 * iqr                        # lower fence (3× IQR)
  hi   <- q3 + 3 * iqr                        # upper fence (3× IQR)
  n_lo <- sum(x < lo)
  n_hi <- sum(x > hi)

  tibble(
    Feature       = var,
    Q1            = round(q1, 3),
    Q3            = round(q3, 3),
    IQR           = round(iqr, 3),
    Lower_Fence   = round(lo, 3),
    Upper_Fence   = round(hi, 3),
    N_Below_Fence = n_lo,
    N_Above_Fence = n_hi,
    Total_Outliers = n_lo + n_hi,
    Pct_Outliers  = round((n_lo + n_hi) / length(x) * 100, 2)
  )
}) %>%
  arrange(desc(Total_Outliers))

cat("\nOutlier summary (3× IQR fences):\n")
print(outlier_summary, n = Inf)

write_csv(outlier_summary, "outputs/outlier_summary.csv")
cat("Saved: outputs/outlier_summary.csv\n")


# --- Step 4b: Winsorize — cap values at 3× IQR fences ---
winsorize_var <- function(x) {
  q1  <- quantile(x, 0.25, na.rm = TRUE)
  q3  <- quantile(x, 0.75, na.rm = TRUE)
  iqr <- q3 - q1
  lo  <- q1 - 3 * iqr
  hi  <- q3 + 3 * iqr
  pmin(pmax(x, lo), hi)   # pmax caps below lo; pmin caps above hi
}

main_data <- main_data %>%
  mutate(across(all_of(continuous_vars), winsorize_var))

cat("\nWinsorizing complete. Extreme values capped at 3× IQR fences.\n")
cat("Rows retained (no deletion):", nrow(main_data), "\n")


# --- Step 4c: Visualise top outlier-affected features ---
# Plot the 6 features with the most outliers before/after is implicitly shown
# via boxplots — winsorized data will have compressed whiskers at fences

top_outlier_vars <- outlier_summary %>%
  filter(Total_Outliers > 0) %>%
  slice_head(n = 6) %>%
  pull(Feature)

if (length(top_outlier_vars) > 0) {

  outlier_plots <- map(top_outlier_vars, function(var) {
    ggplot(main_data, aes(y = .data[[var]])) +
      geom_boxplot(fill = "steelblue4", alpha = 0.7,
                   outlier.colour = "firebrick3",
                   outlier.size   = 1.2,
                   outlier.alpha  = 0.5) +
      labs(title = var, x = NULL, y = NULL) +
      theme_minimal(base_size = 10) +
      theme(axis.text.x = element_blank(),
            axis.ticks.x = element_blank())
  })

  # Combine into one panel using patchwork
  library(patchwork)
  outlier_panel <- wrap_plots(outlier_plots, ncol = 3) +
    plot_annotation(
      title    = "Boxplots: Features with Most Outliers (after winsorizing)",
      subtitle = "Remaining points beyond whiskers are within 3× IQR fences | Red = residual extremes",
      caption  = "NHANES 2017-2020 | Outliers winsorized, not deleted",
      theme    = theme(plot.title = element_text(face = "bold", size = 14))
    )

  ggsave("outputs/plots/00_outlier_boxplots.png", outlier_panel,
         width = 12, height = 8, dpi = 300)
  cat("Saved: outputs/plots/00_outlier_boxplots.png\n")

} else {
  cat("No outliers detected beyond 3× IQR fences.\n")
}


# =============================================================================
# SECTION 5: QUALITY CHECKS
# =============================================================================

# No duplicate respondents — stops script if violated
stopifnot("Duplicate SEQNs detected!" =
            length(unique(main_data$SEQN)) == nrow(main_data))
cat("Check passed: No duplicate respondents.\n")

# Remove single row where target variable is NA (unresolvable)
main_data <- main_data %>% drop_na(hx_fracture)
cat("Rows after removing NA target:", nrow(main_data), "\n")

# Class label distribution — expect imbalance (~85% no fracture, ~15% fracture)
cat("\nClass label distribution (hx_fracture):\n")
print(table(main_data$hx_fracture))
cat("Class imbalance ratio:", round(sum(main_data$hx_fracture == 0) /
                                      sum(main_data$hx_fracture == 1), 1), ":1\n")

# Convert target to factor
main_data <- main_data %>%
  mutate(hx_fracture = factor(hx_fracture, levels = c(0, 1),
                               labels = c("No_Fracture", "Fracture")))


# =============================================================================
# SECTION 6: MISSING DATA SUMMARY & PLOT
# =============================================================================

na_summary <- main_data %>%
  summarise(across(everything(), ~sum(is.na(.)))) %>%
  pivot_longer(everything(), names_to = "feature", values_to = "n_missing") %>%
  mutate(pct_missing = round(n_missing / nrow(main_data) * 100, 1)) %>%
  arrange(desc(pct_missing))

write_csv(na_summary, "outputs/na_summary.csv")
cat("\nTop 10 features with most missing data:\n")
print(head(na_summary, 10))

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


# =============================================================================
# SECTION 7: SAVE
# =============================================================================

dir.create("data/processed", showWarnings = FALSE, recursive = TRUE)
dir.create("outputs/plots", showWarnings = FALSE, recursive = TRUE)

write_csv(main_data, "data/processed/main_data_clean.csv")
cat("\nSaved: data/processed/main_data_clean.csv\n")
cat("Script 01 complete. Run 02_t_score_calc.R next.\n")
