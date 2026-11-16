#From https://github.com/RobertsLab/sormi-assay-development/blob/main/Resazurin/code/03.5-resazurin-curve-features-all-experiments.Rmd
# Helper Functions

```{r helper-fns}
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
