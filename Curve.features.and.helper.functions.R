###Helper functions and curve feature extraction functions from Roberts lab github - used to extract features of resazurin curve data. 
### https://github.com/RobertsLab/sormi-assay-development/blob/main/Resazurin/code/03.5-resazurin-curve-features-all-experiments.Rmd

safe_divide <- function(num, den) {
  if_else(is.finite(num) & is.finite(den) & abs(den) > .Machine$double.eps,
          num / den, NA_real_)
}

trapezoid_auc <- function(time_hr, value) {
  ok <- is.finite(time_hr) & is.finite(value)
  t <- time_hr[ok]
  v <- value[ok]
  if (length(t) < 2) return(NA_real_)
  ord <- order(t)
  t <- t[ord]
  v <- v[ord]
  sum(diff(t) * (head(v, -1) + tail(v, -1)) / 2)
}

segment_auc <- function(time_hr, value) {
  ok <- is.finite(time_hr) & is.finite(value)
  t <- time_hr[ok]
  v <- value[ok]
  if (length(t) < 2) {
    return(tibble(t_start = numeric(), t_end = numeric(),
                  t_mid = numeric(), auc = numeric()))
  }
  ord <- order(t)
  t <- t[ord]
  v <- v[ord]
  tibble(
    t_start = head(t, -1),
    t_end   = tail(t, -1),
    t_mid   = (head(t, -1) + tail(t, -1)) / 2,
    auc     = diff(t) * (head(v, -1) + tail(v, -1)) / 2
  )
}

classify_trajectory <- function(slopes, slope_tol = 1e-12) {
  if (length(slopes) == 0 || all(!is.finite(slopes))) return(NA_character_)
  s <- slopes[is.finite(slopes)]
  signs <- case_when(
    s >  slope_tol ~  1L,
    s < -slope_tol ~ -1L,
    TRUE           ~  0L
  )
  nonzero <- signs[signs != 0L]
  if (length(nonzero) == 0) return("stable")
  if (all(nonzero > 0)) return("monotonic_increase")
  if (all(nonzero < 0)) return("monotonic_decrease")
  if (any(nonzero > 0) && any(nonzero < 0)) {
    first_neg <- which(nonzero < 0)[1]
    first_pos <- which(nonzero > 0)[1]
    if (first_pos < first_neg) return("rise_then_depress")
    if (first_neg < first_pos) return("depress_then_recover")
  }
  "complex"
}

extract_curve_features_one <- function(df, value_col, metric_label) {
  df <- df %>%
    filter(is.finite(time_hr), is.finite(.data[[value_col]])) %>%
    arrange(time_hr)

  if (nrow(df) < 2) return(tibble())

  t <- df$time_hr
  v <- df[[value_col]]
  dt <- diff(t)
  dv <- diff(v)
  valid_slope <- is.finite(dt) & dt > 0 & is.finite(dv)
  slopes <- dv[valid_slope] / dt[valid_slope]
  slope_midpoints <- (head(t, -1) + tail(t, -1)) / 2
  slope_midpoints <- slope_midpoints[valid_slope]
  seg_auc <- segment_auc(t, v)
  split_time <- median(seg_auc$t_mid, na.rm = TRUE)
  early_auc <- seg_auc %>%
    filter(t_mid <= split_time) %>%
    summarise(x = sum(auc, na.rm = TRUE), .groups = "drop") %>%
    pull(x)
  late_auc <- seg_auc %>%
    filter(t_mid > split_time) %>%
    summarise(x = sum(auc, na.rm = TRUE), .groups = "drop") %>%
    pull(x)

  if (length(early_auc) == 0) early_auc <- NA_real_
  if (length(late_auc) == 0) late_auc <- NA_real_

  baseline_value <- v[1]
  final_value <- v[length(v)]
  peak_idx <- which.max(v)
  trough_idx <- which.min(v)
  peak_value <- v[peak_idx]
  trough_value <- v[trough_idx]
  peak_after_baseline <- peak_value - baseline_value
  depression_abs <- if (peak_idx < length(v)) {
    max(peak_value - final_value, 0, na.rm = TRUE)
  } else {
    0
  }

  vmax_idx <- if (length(slopes) > 0) which.max(slopes) else NA_integer_
  min_slope_idx <- if (length(slopes) > 0) which.min(slopes) else NA_integer_
  curvature <- diff(slopes)
  inflection_idx <- if (length(curvature) > 0) {
    which.max(abs(curvature)) + 1L
  } else {
    NA_integer_
  }

  tibble(
    value_metric = metric_label,
    baseline_value = baseline_value,
    final_value = final_value,
    peak_value = peak_value,
    trough_value = trough_value,
    time_to_peak = t[peak_idx],
    time_to_trough = t[trough_idx],
    auc_total = trapezoid_auc(t, v),
    auc_early = early_auc,
    auc_late = late_auc,
    delta_auc_late_minus_early = late_auc - early_auc,
    late_early_auc_ratio = safe_divide(late_auc, early_auc),
    initial_slope = slopes[1],
    vmax = if (length(slopes) > 0) max(slopes, na.rm = TRUE) else NA_real_,
    time_to_vmax = if (!is.na(vmax_idx)) slope_midpoints[vmax_idx] else NA_real_,
    min_slope = if (length(slopes) > 0) min(slopes, na.rm = TRUE) else NA_real_,
    time_to_min_slope = if (!is.na(min_slope_idx)) {
      slope_midpoints[min_slope_idx]
    } else {
      NA_real_
    },
    inflection_time = if (!is.na(inflection_idx)) {
      slope_midpoints[inflection_idx]
    } else {
      NA_real_
    },
    metabolic_scope = peak_after_baseline,
    final_delta = final_value - baseline_value,
    depression_abs = depression_abs,
    metabolic_depression_index = safe_divide(depression_abs,
                                             abs(peak_after_baseline)),
    resilience_ratio = safe_divide(final_value, peak_value),
    stability_cv = safe_divide(sd(v, na.rm = TRUE), abs(mean(v, na.rm = TRUE))),
    trajectory_class = classify_trajectory(slopes),
    n_timepoints = nrow(df)
  )
}

compute_curve_features <- function(df, value_col, metric_label) {
  id_vars <- intersect(
    c("experiment", "trace_id", "family_id_group", "treatment_group",
      "round_group", "cup_id_group", "sample_id_group"),
    names(df)
  )

  df %>%
    filter(!is.na(family_id_group), is.finite(.data[[value_col]])) %>%
    group_by(across(all_of(id_vars))) %>%
    group_modify(~ extract_curve_features_one(.x, value_col, metric_label)) %>%
    ungroup()
}

tidy_tukey <- function(model, term) {
  tk <- tryCatch(TukeyHSD(aov(model), which = term)[[term]],
                 error = function(e) NULL)
  if (is.null(tk)) return(tibble())
  as_tibble(tk, rownames = "contrast") %>%
    rename(
      estimate = diff,
      p.value = `p adj`
    ) %>%
    mutate(SE = NA_real_, df = NA_real_, t.ratio = NA_real_) %>%
    select(contrast, estimate, SE, df, t.ratio, p.value)
}

run_feature_family_stats <- function(df) {
  df <- df %>%
    filter(is.finite(feature_value), !is.na(family_id_group)) %>%
    mutate(family = factor(family_id_group))

  if (nrow(df) < 3 || n_distinct(df$family) < 2) return(tibble())

  model <- lm(feature_value ~ family, data = df)
  aov_tab <- anova(model)

  ss_family <- aov_tab["family", "Sum Sq"]
  df_family <- aov_tab["family", "Df"]
  ms_error <- aov_tab["Residuals", "Mean Sq"]
  ss_total <- sum(aov_tab$`Sum Sq`, na.rm = TRUE)
  eta_sq <- ss_family / ss_total
  omega_sq <- (ss_family - df_family * ms_error) / (ss_total + ms_error)

  fam_means <- df %>%
    group_by(family) %>%
    summarise(mean = mean(feature_value), .groups = "drop")

  pair_grid <- tidyr::crossing(fam1 = fam_means$family,
                               fam2 = fam_means$family) %>%
    filter(as.character(fam1) < as.character(fam2)) %>%
    left_join(fam_means, by = c("fam1" = "family")) %>%
    rename(mean1 = mean) %>%
    left_join(fam_means, by = c("fam2" = "family")) %>%
    rename(mean2 = mean) %>%
    mutate(abs_d = abs(mean1 - mean2) / sqrt(ms_error))

  kw <- kruskal.test(feature_value ~ family, data = df)

  tibble(
    n = nrow(df),
    n_family = n_distinct(df$family),
    f_statistic = aov_tab["family", "F value"],
    p_value = aov_tab["family", "Pr(>F)"],
    kruskal_p_value = kw$p.value,
    eta_sq = eta_sq,
    omega_sq = max(omega_sq, 0, na.rm = TRUE),
    max_abs_pairwise_d = max(pair_grid$abs_d, na.rm = TRUE)
  )
}

run_feature_pairwise <- function(df) {
  df <- df %>%
    filter(is.finite(feature_value), !is.na(family_id_group)) %>%
    mutate(family = factor(family_id_group))

  if (nrow(df) < 3 || n_distinct(df$family) < 2) return(tibble())

  model <- lm(feature_value ~ family, data = df)
  tidy_tukey(model, "family")
}

make_palette <- function(n) {
  okabe_ito_7 <- c(
    "#E69F00", "#56B4E9", "#009E73", "#F0E442",
    "#0072B2", "#D55E00", "#CC79A7"
  )
  if (n == 0L) return(character(0))
  if (n <= length(okabe_ito_7)) return(okabe_ito_7[seq_len(n)])
  colorspace::qualitative_hcl(n, palette = "Dynamic")
}
