library(glmnet)
library(dplyr)
library(ggplot2)
library(patchwork)  # for side-by-side plots; install.packages("patchwork") if needed

# ---- Prepare data ----
auc_comp_minslope_vmax <- auc_comp_minslope_vmax %>%
  mutate(
    z_glyc = as.numeric(scale(mean_glyc)),
    z_vmax = as.numeric(scale(mean_Vmax)),
    glyc_vmax_composite = (z_glyc + z_vmax) / 2
  )

predictors <- c("glyc_vmax_composite", "mean_AUC_42", "mean_AUC_21", "mean_time_min_slope")
target <- "composite_score"

x <- as.matrix(auc_comp_minslope_vmax[, predictors])
y <- auc_comp_minslope_vmax[[target]]
n <- nrow(auc_comp_minslope_vmax)

# ---- LASSO with LOOCV ----
set.seed(123)
lasso_cv <- cv.glmnet(x, y, alpha = 1, nfolds = n, standardize = TRUE)
best_lambda <- lasso_cv$lambda.min

lasso_preds <- numeric(n)
for (i in 1:n) {
  fit_i <- glmnet(x[-i, ], y[-i], alpha = 1, lambda = best_lambda, standardize = TRUE)
  lasso_preds[i] <- predict(fit_i, newx = x[i, , drop = FALSE])
}

r2_lasso <- cor(y, lasso_preds)^2
rmse_lasso <- sqrt(mean((y - lasso_preds)^2))
cat("LOOCV R¬≤:", round(r2_lasso, 3), "  RMSE:", round(rmse_lasso, 3), "\n")

pred_df <- data.frame(family = auc_comp_minslope_vmax$family, actual = y, predicted = lasso_preds)

# ---- Relative influence ----
final_coefs <- as.matrix(coef(lasso_cv, s = "lambda.min"))[-1, ]  # drop intercept

influence_df <- data.frame(variable = names(final_coefs), coefficient = final_coefs) %>%
  mutate(abs_coef = abs(coefficient),
         rel_influence_pct = 100 * abs_coef / sum(abs_coef),
         direction = ifelse(coefficient < 0, "Negative", "Positive")) %>%
  arrange(rel_influence_pct)

# ---- Plot 1: Actual vs Predicted ----
p1 <- ggplot(pred_df, aes(x = actual, y = predicted, label = family)) +
  geom_point(size = 3, color = "steelblue") +
  geom_text(vjust = -1, size = 3.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  annotate("text", x = min(pred_df$actual), y = max(pred_df$predicted),
           label = paste0("LOOCV R¬≤ = ", round(r2_lasso, 3)),
           hjust = 0, vjust = 1, size = 4) +
  labs(x = "Actual composite_score", y = "Predicted composite_score",
       title = "LOOCV: Actual vs Predicted") +
  theme_minimal()

# ---- Plot 2: Relative Influence ----
p2 <- ggplot(influence_df, aes(x = reorder(variable, rel_influence_pct), y = rel_influence_pct, fill = direction)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_manual(values = c("Negative" = "indianred", "Positive" = "steelblue")) +
  geom_text(aes(label = paste0(round(rel_influence_pct, 1), "%")), hjust = -0.1, size = 3.5) +
  ylim(0, max(influence_df$rel_influence_pct) * 1.2) +
  labs(x = "Variable", y = "Relative influence (%)",
       title = "Relative Influence on composite_score") +
  theme_minimal()

# ---- Combine side by side ----
p1 + p2
