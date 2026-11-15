predictors <- c("mean_glyc", "auc_total_21C", "auc_total_42C")

# Keep the outcome and predictors; drop any family with missing values
df <- as.data.frame(merged_all_data)[, c("family", "composite_score", predictors)]
df <- na.omit(df)
y  <- df$composite_score
n  <- nrow(df)
model <- reformulate(predictors, response = "composite_score")

#  Fit on all families 
fit <- lm(model, data = df)
summary(fit)          # coefficients, R^2, adjusted R^2, overall F-test p-value
confint(fit)          # 95% confidence intervals

# Collinearity: VIF = 1 / (1 - R^2 of each predictor regressed on the others)
vif <- sapply(predictors, function(v)
  1 / (1 - summary(lm(reformulate(setdiff(predictors, v), v), data = df))$r.squared))
round(vif, 2)

#  Leave-one-family-out cross-validation 

pred     <- numeric(n)
std_coef <- matrix(NA, n, length(predictors), dimnames = list(NULL, predictors))

for (i in 1:n) {
  fit_i   <- lm(model, data = df[-i, ])
  pred[i] <- predict(fit_i, newdata = df[i, ])
  # standardized coefficients (coefficient x predictor SD), for relative influence
  std_coef[i, ] <- coef(fit_i)[predictors] * sapply(df[-i, predictors], sd)
}

q2   <- 1 - sum((y - pred)^2) / sum((y - mean(y))^2)   # cross-validated R^2
rmse <- sqrt(mean((y - pred)^2))                      # cross-validated RMSE
cat(sprintf("LOOCV: Q2 = %.2f, RMSE = %.1f\n", q2, rmse))

results <- data.frame(family = df$family, observed = y,
                      predicted = round(pred, 1), error = round(pred - y, 1))
results

#  Plot: observed vs predicted 
plot(y, pred, pch = 19, xlab = "Observed composite score",
     ylab = "Predicted (left-out family)",
     main = sprintf("OLS leave-one-family-out CV\nQ2 = %.2f, RMSE = %.1f", q2, rmse))
abline(0, 1, lty = 2)                       
text(y, pred, df$family, pos = 3)

#  Plot, relative influence 
# Each predictor's share (%) of the summed absolute standardized coefficients,
# averaged over the cross-validation fits
influence <- colMeans(abs(std_coef))
influence <- sort(100 * influence / sum(influence))
direction <- sign(coef(fit)[names(influence)])

par(mar = c(4.5, 9, 3, 1))
barplot(influence, horiz = TRUE, las = 1, xlim = c(0, 100),
        col = ifelse(direction > 0, "steelblue", "firebrick"),
        xlab = "Relative influence (%)", main = "Relative influence")
legend("bottomright", c("positive", "negative"),
       fill = c("steelblue", "firebrick"), bty = "n")
par(mar = c(5, 4, 4, 2) + 0.1)
round(influence, 1)
