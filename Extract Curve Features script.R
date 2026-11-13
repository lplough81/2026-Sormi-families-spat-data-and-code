# Input: `df` with columns
#   date, plate, well, timepoint, Family, temperature, type, size,
#   resazurin_counts, fluorescence.norm, value.mm

# Assumes helper functions are already loaded:
#   safe_divide, trapezoid_auc, segment_auc, classify_trajectory,
#   extract_curve_features_one
 
library(dplyr)
library(tidyr)
library(tibble)
 
select <- dplyr::select
filter <- dplyr::filter
 
# Settings 
value_col    <- "value.mm"   # or "fluorescence.norm" for unscaled fold change
metric_label <- value_col
min_points   <- 3            # curves with fewer timepoints are dropped
 
# Map this dataset onto the names the feature functions expect 
df_curves <- df %>%
  ungroup() %>%                                  # input is a grouped tibble
  mutate(
    time_hr         = suppressWarnings(as.numeric(as.character(timepoint))),
    family_id_group = Family,
    treatment_group = temperature,
    trace_id        = paste(date, plate, well, sep = "_"),
    value           = suppressWarnings(as.numeric(.data[[value_col]]))
  )
 
 
# Keep sample wells only
print(count(df_curves, type))                    # check what `type` contains
 
df_clean <- df_curves %>%
  filter(
    !grepl("blank", type, ignore.case = TRUE),   # drop blank wells
    !is.na(family_id_group),
    is.finite(time_hr),
    is.finite(value)
  )
 
id_vars <- c("trace_id", "date", "plate", "well",
             "family_id_group", "treatment_group")
 
 
df_clean <- df_clean %>%
  group_by(across(all_of(id_vars)), time_hr) %>%
  summarise(value = mean(value), .groups = "drop")
 
# How many timepoints each curve has
points_per_curve <- df_clean %>% count(across(all_of(id_vars)), name = "n_points")
print(count(points_per_curve, n_points, name = "n_curves"))
 
df_clean <- df_clean %>%
  semi_join(filter(points_per_curve, n_points >= min_points), by = id_vars)
 
# Features for each curve 
features <- df_clean %>%
  group_by(across(all_of(id_vars))) %>%
  group_modify(~ extract_curve_features_one(.x, "value", metric_label)) %>%
  ungroup()
 
if (nrow(features) == 0) {
  stop("No curve produced features; check `points_per_curve` and the filters.")
}
 
feature_cols <- features %>%
  select(where(is.numeric)) %>%
  names() %>%
  setdiff(c(id_vars, "n_timepoints"))   # n_timepoints is bookkeeping, not a feature
 
# Summary by Family x temperature 
summarise_features <- function(dat, group_vars) {
  dat %>%
    pivot_longer(all_of(feature_cols),
                 names_to = "feature", values_to = "feature_value") %>%
    filter(is.finite(feature_value)) %>%
    group_by(across(all_of(group_vars)), feature) %>%
    summarise(
      n_curves = n(),
      mean     = mean(feature_value),
      sd       = sd(feature_value),
      se       = sd / sqrt(n_curves),
      median   = median(feature_value),
      .groups  = "drop"
    )
}
 
family_temp_summary <- summarise_features(
  features, c("family_id_group", "treatment_group"))
family_summary <- summarise_features(features, "family_id_group")  # pooled over temperature
 
family_temp_means_wide <- family_temp_summary %>%
  select(family_id_group, treatment_group, feature, mean) %>%
  pivot_wider(names_from = feature, values_from = mean)

write.csv(family_temp_means_wide,"42.Fam.means.csv")
