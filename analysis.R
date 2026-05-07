################################################################################
# Replication of Stock & Watson (1999) "Forecasting Inflation"
# Journal of Monetary Economics, 44(2), 293-335
#
# MATH 1220: Mathematical Statistics — Final Project
# 
# This script replicates the core analysis: comparing Phillips curve inflation
# forecasts against AR benchmarks using pseudo out-of-sample evaluation.
# 
# Data: CPI-U, Core CPI, Unemployment Rate from FRED (1959:01 - present)
# Forecast horizon: h = 12 months
################################################################################

# ==============================================================================
# 0. SETUP & PACKAGES
# ==============================================================================

# Install packages if needed (uncomment as necessary)
# install.packages(c("fredr", "forecast", "lmtest", "sandwich", "tseries",
#                     "nortest", "strucchange", "ggplot2"))

library(fredr)      # FRED API access
library(forecast)   # ARIMA, auto.arima
library(lmtest)     # coeftest, bptest, bgtest
library(sandwich)   # HAC standard errors
library(tseries)    # adf.test, jarque.bera.test
library(ggplot2)    # plotting

# -------------------------------------------------------
# You need a free FRED API key from https://fred.stlouisfed.org/docs/api/api_key.html
# Set it here:
# fredr_set_key("YOUR_API_KEY_HERE")
#
# ALTERNATIVE: If you don't want to use the API, you can download CSVs directly
# from FRED and read them in with read.csv(). See Section 0B below.
# -------------------------------------------------------


# ==============================================================================
# 0B. ALTERNATIVE DATA LOADING (no API key needed)
# ==============================================================================
# If you prefer not to use the fredr package, download CSVs from:
#   - https://fred.stlouisfed.org/series/CPIAUCSL  (CPI-U, all items)
#   - https://fred.stlouisfed.org/series/CPILFESL  (Core CPI, less food & energy)
#   - https://fred.stlouisfed.org/series/UNRATE    (Civilian unemployment rate)
# Then uncomment and use:
#
# cpi_raw   <- read.csv("CPIAUCSL.csv")
# core_raw  <- read.csv("CPILFESL.csv")
# unemp_raw <- read.csv("UNRATE.csv")
#
# cpi_raw$date   <- as.Date(cpi_raw$DATE)
# cpi_raw$value  <- as.numeric(cpi_raw$CPIAUCSL)
# core_raw$date  <- as.Date(core_raw$DATE)
# core_raw$value <- as.numeric(core_raw$CPILFESL)
# unemp_raw$date <- as.Date(unemp_raw$DATE)
# unemp_raw$value <- as.numeric(unemp_raw$UNRATE)
#
# Then skip Section 1 and go directly to Section 2, replacing the merge logic.


# ==============================================================================
# 1. DATA: Pull from FRED
# ==============================================================================

cat("Pulling data from FRED...\n")

cpi_all  <- fredr(series_id = "CPIAUCSL",  # CPI-U All Items
                  observation_start = as.Date("1957-01-01"),
                  observation_end   = as.Date("2025-12-01"))

cpi_core <- fredr(series_id = "CPILFESL",  # Core CPI (less food & energy)
                  observation_start = as.Date("1957-01-01"),
                  observation_end   = as.Date("2025-12-01"))

unemp    <- fredr(series_id = "UNRATE",    # Civilian unemployment rate
                  observation_start = as.Date("1957-01-01"),
                  observation_end   = as.Date("2025-12-01"))


# ==============================================================================
# 2. CONSTRUCT VARIABLES
# ==============================================================================
# Stock & Watson define:
#   pi_t     = 1200 * ln(P_t / P_{t-1})          (monthly inflation, annualized)
#   pi^h_t   = (1200/h) * ln(P_t / P_{t-h})      (h-month inflation, annualized)
#   Dependent variable: pi^h_{t+h} - pi_t
#
# We set h = 12 (12-month forecast horizon)

h <- 12  # forecast horizon in months

# Merge datasets
df <- merge(
  data.frame(date = cpi_all$date, cpi = cpi_all$value),
  data.frame(date = unemp$date,  urate = unemp$value),
  by = "date"
)
df <- merge(
  df,
  data.frame(date = cpi_core$date, cpi_core = cpi_core$value),
  by = "date",
  all.x = TRUE
)
df <- df[order(df$date), ]

# Compute inflation measures
df$pi      <- c(NA, 1200 * diff(log(df$cpi)))            # monthly inflation (annualized)
df$pi_core <- c(NA, 1200 * diff(log(df$cpi_core)))       # core monthly inflation

# h-period inflation (annualized)
n <- nrow(df)
df$pi_h <- NA
for (i in (h + 1):n) {
  df$pi_h[i] <- (1200 / h) * log(df$cpi[i] / df$cpi[i - h])
}

# Dependent variable: acceleration of inflation = pi^h_{t+h} - pi_t
# We need to shift pi_h forward by h periods (i.e., align future inflation with current data)
df$y <- NA
for (i in 1:(n - h)) {
  df$y[i] <- df$pi_h[i + h] - df$pi[i]
}

# Change in monthly inflation (for lagged regressors)
df$dpi <- c(NA, diff(df$pi))

# Trim to sample starting 1959:01 (as in the paper)
df <- df[df$date >= as.Date("1959-01-01"), ]
rownames(df) <- NULL

cat("Dataset constructed:", nrow(df), "observations from",
    as.character(min(df$date)), "to", as.character(max(df$date)), "\n\n")


# ==============================================================================
# PART 1: DESCRIPTIVE STATISTICS & NONPARAMETRIC TESTS
# Topics: Weeks 1-4 (probability review, EDF, CLT, histograms, Q-Q plots)
#         Week 8 (K-S test, chi-square goodness of fit)
# ==============================================================================

cat("=== PART 1: Descriptive Statistics ===\n\n")

# Summary statistics for key variables
vars_of_interest <- c("pi", "urate", "dpi")
for (v in vars_of_interest) {
  x <- na.omit(df[[v]])
  cat(sprintf("%-8s  mean=%.3f  sd=%.3f  skew=%.3f  n=%d\n",
              v, mean(x), sd(x),
              mean((x - mean(x))^3) / sd(x)^3,
              length(x)))
}
cat("\n")

# --- Histogram of monthly inflation ---
p1 <- ggplot(df[!is.na(df$pi), ], aes(x = pi)) +
  geom_histogram(aes(y = after_stat(density)), bins = 50, fill = "steelblue", alpha = 0.7) +
  stat_function(fun = dnorm, args = list(mean = mean(df$pi, na.rm = TRUE),
                                          sd = sd(df$pi, na.rm = TRUE)),
                color = "red", linewidth = 0.8) +
  labs(title = "Distribution of Monthly Inflation (annualized)",
       x = expression(pi[t]), y = "Density") +
  theme_minimal()
print(p1)

# --- Q-Q plot for inflation changes ---
p2 <- ggplot(df[!is.na(df$dpi), ], aes(sample = dpi)) +
  stat_qq() + stat_qq_line(color = "red") +
  labs(title = "Q-Q Plot: Changes in Monthly Inflation",
       x = "Theoretical Quantiles", y = "Sample Quantiles") +
  theme_minimal()
print(p2)

# --- Kolmogorov-Smirnov test for normality of dpi ---
dpi_clean <- na.omit(df$dpi)
ks_result <- ks.test(dpi_clean, "pnorm", mean(dpi_clean), sd(dpi_clean))
cat("Kolmogorov-Smirnov test for normality of dpi:\n")
cat(sprintf("  D = %.4f, p-value = %.4f\n\n", ks_result$statistic, ks_result$p.value))

# --- Jarque-Bera test (uses skewness & kurtosis — connects to moment estimators) ---
jb_result <- jarque.bera.test(ts(dpi_clean))
cat("Jarque-Bera test for normality of dpi:\n")
cat(sprintf("  JB = %.2f, p-value = %.6f\n\n", jb_result$statistic, jb_result$p.value))

# --- Chi-square goodness-of-fit test (bin residuals, test against normal) ---
# This manually implements the Pearson chi-square statistic from Week 8
k <- 10  # number of bins
breaks <- qnorm(seq(0, 1, length.out = k + 1), mean(dpi_clean), sd(dpi_clean))
breaks[1] <- -Inf; breaks[k + 1] <- Inf
observed <- table(cut(dpi_clean, breaks = breaks))
expected <- rep(length(dpi_clean) / k, k)
chi2_stat <- sum((as.numeric(observed) - expected)^2 / expected)
chi2_pval <- 1 - pchisq(chi2_stat, df = k - 1 - 2)  # subtract 2 for estimated mu, sigma
cat("Chi-square goodness-of-fit test (k=10 bins, normal null):\n")
cat(sprintf("  chi2 = %.2f, df = %d, p-value = %.4f\n\n", chi2_stat, k - 1 - 2, chi2_pval))


# ==============================================================================
# PART 2: PHILLIPS CURVE ESTIMATION VIA OLS
# Topics: Weeks 10-12 (linear regression, confidence intervals, t-tests, F-tests)
#         Week 3 (delta method for NAIRU CI)
# ==============================================================================

cat("=== PART 2: Phillips Curve OLS Estimation ===\n\n")

# Prepare lagged variables for regression
# Following S&W: use 12 lags of unemployment, 12 lags of dpi
max_lag <- 12

# Create lagged columns
for (j in 0:max_lag) {
  df[[paste0("urate_L", j)]] <- c(rep(NA, j), df$urate[1:(n - j)])[1:nrow(df)]
  if (j > 0) {
    df[[paste0("dpi_L", j)]]  <- c(rep(NA, j), df$dpi[1:(n - j)])[1:nrow(df)]
  }
}

# --- Baseline Phillips curve: 4 lags of u, 4 lags of dpi ---
# (S&W use various lag lengths; 4 is a standard starting point)
pc_formula <- as.formula(
  paste("y ~", paste0("urate_L", 0:3, collapse = " + "), "+",
        paste0("dpi_L", 1:4, collapse = " + "))
)

# Estimate on full available sample (complete cases only)
pc_data <- df[complete.cases(df[, c("y", paste0("urate_L", 0:3), paste0("dpi_L", 1:4))]), ]
pc_ols  <- lm(pc_formula, data = pc_data)

cat("--- Baseline Phillips Curve (4 lags u, 4 lags dpi) ---\n")
summary(pc_ols)

# HAC standard errors (Newey-West) — appropriate for h-step-ahead overlapping errors
cat("\n--- HAC (Newey-West) Standard Errors ---\n")
print(coeftest(pc_ols, vcov = NeweyWest(pc_ols, lag = h - 1)))

# --- F-test: joint significance of unemployment lags ---
# H0: all beta_urate = 0 (unemployment has no predictive power for inflation)
pc_restricted <- lm(
  as.formula(paste("y ~", paste0("dpi_L", 1:4, collapse = " + "))),
  data = pc_data
)
f_test <- anova(pc_restricted, pc_ols)
cat("\n--- F-test: Joint significance of unemployment lags ---\n")
print(f_test)
cat(sprintf("F = %.3f, p-value = %.4f\n\n", f_test$F[2], f_test$`Pr(>F)`[2]))

# --- Confidence intervals for coefficients (95%) ---
cat("--- 95% Confidence Intervals for OLS Coefficients ---\n")
print(confint(pc_ols, level = 0.95))
cat("\n")

# --- Delta method for implied NAIRU ---
# If the Phillips curve is: y = alpha + beta(1)*u + gamma(L)*dpi + e
# Setting y = 0 and dpi = 0: u_nairu = -alpha / sum(beta_j)
# This is a nonlinear function of parameters -> use delta method (Week 3)
coefs <- coef(pc_ols)
alpha_hat <- coefs["(Intercept)"]
beta_sum  <- sum(coefs[grep("urate", names(coefs))])
nairu_hat <- -alpha_hat / beta_sum
cat(sprintf("Implied NAIRU = %.2f%%\n", nairu_hat))

# Delta method SE: g(theta) = -alpha/sum(beta), gradient = (-1/sum(beta), alpha/sum(beta)^2, ...)
# For a rough CI, use bootstrap or numerical delta method:
V <- vcov(pc_ols)
urate_idx <- grep("urate|Intercept", names(coefs))
grad <- rep(0, length(coefs))
grad[1] <- -1 / beta_sum                                       # d/d(alpha)
for (j in grep("urate", names(coefs))) {
  grad[j] <- alpha_hat / beta_sum^2                             # d/d(beta_j)
}
nairu_se <- sqrt(t(grad) %*% V %*% grad)
cat(sprintf("Delta method SE = %.3f\n", nairu_se))
cat(sprintf("95%% CI for NAIRU: [%.2f, %.2f]\n\n",
            nairu_hat - 1.96 * nairu_se, nairu_hat + 1.96 * nairu_se))


# ==============================================================================
# PART 3: MLE UNDER NORMALITY & LIKELIHOOD RATIO TESTS
# Topics: Weeks 5-6 (MLE, sufficiency, Cramer-Rao)
#         Week 7 (LRT, hypothesis testing, decision theory / AIC / BIC)
# ==============================================================================

cat("=== PART 3: MLE & Likelihood Ratio Tests ===\n\n")

# Under Gaussian errors, OLS = MLE. We verify this and extract the log-likelihood.
loglik_full <- logLik(pc_ols)
cat(sprintf("Log-likelihood (full model, 4 lags): %.2f\n", loglik_full))

# --- Restricted model: 2 lags only ---
pc_2lag <- lm(
  as.formula(paste("y ~", paste0("urate_L", 0:1, collapse = " + "), "+",
                   paste0("dpi_L", 1:2, collapse = " + "))),
  data = pc_data
)
loglik_restricted <- logLik(pc_2lag)
cat(sprintf("Log-likelihood (restricted, 2 lags): %.2f\n", loglik_restricted))

# LRT statistic: Lambda = -2 * (l_restricted - l_full) ~ chi-sq(df)
lrt_stat <- -2 * (as.numeric(loglik_restricted) - as.numeric(loglik_full))
lrt_df   <- length(coef(pc_ols)) - length(coef(pc_2lag))
lrt_pval <- 1 - pchisq(lrt_stat, df = lrt_df)
cat(sprintf("LRT statistic = %.3f, df = %d, p-value = %.4f\n", lrt_stat, lrt_df, lrt_pval))
cat("Interpretation: ", ifelse(lrt_pval < 0.05,
    "Reject H0 — additional lags improve the model.",
    "Fail to reject H0 — parsimonious model sufficient."), "\n\n")

# --- Model selection via AIC and BIC (connects to decision theory / loss functions) ---
cat("--- Model Selection (AIC / BIC) ---\n")
cat("AIC and BIC penalize complexity differently. AIC minimizes expected\n")
cat("Kullback-Leibler divergence (a specific loss function from decision theory).\n\n")

lag_specs <- list(
  "2 lags" = list(u = 0:1, dpi = 1:2),
  "4 lags" = list(u = 0:3, dpi = 1:4),
  "6 lags" = list(u = 0:5, dpi = 1:6),
  "8 lags" = list(u = 0:7, dpi = 1:8)
)

aic_bic_table <- data.frame(Spec = character(), AIC = numeric(), BIC = numeric(),
                            stringsAsFactors = FALSE)
for (spec_name in names(lag_specs)) {
  spec <- lag_specs[[spec_name]]
  f <- as.formula(paste("y ~",
    paste0("urate_L", spec$u, collapse = " + "), "+",
    paste0("dpi_L", spec$dpi, collapse = " + ")))
  fit <- lm(f, data = pc_data)
  aic_bic_table <- rbind(aic_bic_table,
    data.frame(Spec = spec_name, AIC = round(AIC(fit), 2), BIC = round(BIC(fit), 2)))
}
print(aic_bic_table, row.names = FALSE)
cat("\n")


# ==============================================================================
# PART 4: AR BENCHMARK & TIME SERIES MODELS
# Topics: Weeks 12-13 (AR, MA, ARIMA, GARCH, ACF/PACF)
# ==============================================================================

cat("=== PART 4: AR Benchmark & Time Series Models ===\n\n")

# --- AR(p) benchmark for inflation changes ---
# Stock & Watson compare Phillips curve against a univariate AR model
# AR order selected by AIC (as in the paper)

y_ts <- ts(na.omit(df$y), frequency = 12)

# Fit AR by AIC
ar_fit <- ar(y_ts, order.max = 12, method = "ols", aic = TRUE)
cat(sprintf("AR benchmark: order selected by AIC = %d\n", ar_fit$order))
cat(sprintf("AR coefficients:\n"))
print(round(ar_fit$ar, 4))
cat("\n")

# --- ARIMA via auto.arima ---
pi_ts <- ts(na.omit(df$pi), frequency = 12, start = c(1959, 2))
arima_fit <- auto.arima(pi_ts, max.p = 6, max.q = 6, max.d = 2,
                         seasonal = FALSE, stepwise = FALSE, approximation = FALSE)
cat("auto.arima selected model:\n")
print(arima_fit)
cat(sprintf("\nAIC = %.2f, BIC = %.2f\n\n", AIC(arima_fit), BIC(arima_fit)))

# --- ACF / PACF of inflation (useful for identifying AR/MA orders) ---
cat("Plotting ACF and PACF of monthly inflation...\n")
par(mfrow = c(1, 2))
acf(na.omit(df$pi), lag.max = 36, main = "ACF of Monthly Inflation")
pacf(na.omit(df$pi), lag.max = 36, main = "PACF of Monthly Inflation")
par(mfrow = c(1, 1))


# ==============================================================================
# PART 5: PSEUDO OUT-OF-SAMPLE FORECAST COMPARISON
# This is the core of Stock & Watson (1999)
# Topics: OLS, time series, MSE (loss function from decision theory)
# ==============================================================================

cat("\n=== PART 5: Pseudo Out-of-Sample Forecast Evaluation ===\n\n")

# Setup: rolling estimation window
# Training starts at 1959:01, first forecast origin at 1970:01
# (gives ~11 years of initial training data)

forecast_start <- as.Date("1970-01-01")
forecast_rows  <- which(df$date >= forecast_start & !is.na(df$y))

# Storage for forecast errors
errors_pc <- numeric(length(forecast_rows))
errors_ar <- numeric(length(forecast_rows))

cat("Running pseudo out-of-sample evaluation...\n")
cat(sprintf("Forecast origins: %d\n", length(forecast_rows)))

for (i in seq_along(forecast_rows)) {
  t_idx <- forecast_rows[i]

  # Training data: all observations before forecast origin
  train <- df[1:(t_idx - 1), ]
  train <- train[complete.cases(train[, c("y", paste0("urate_L", 0:3),
                                           paste0("dpi_L", 1:4))]), ]

  # Actual outcome
  actual <- df$y[t_idx]

  # --- Phillips curve forecast ---
  pc_fit <- tryCatch(
    lm(pc_formula, data = train),
    error = function(e) NULL
  )
  if (!is.null(pc_fit)) {
    pred_pc <- predict(pc_fit, newdata = df[t_idx, , drop = FALSE])
    errors_pc[i] <- actual - pred_pc
  } else {
    errors_pc[i] <- NA
  }

  # --- AR benchmark forecast ---
  # Simple: regress y on 4 lags of y
  train$y_L1 <- c(NA, train$y[-nrow(train)])
  train$y_L2 <- c(NA, NA, train$y[1:(nrow(train) - 2)])
  train$y_L3 <- c(NA, NA, NA, train$y[1:(nrow(train) - 3)])
  train$y_L4 <- c(NA, NA, NA, NA, train$y[1:(nrow(train) - 4)])

  ar_train <- train[complete.cases(train[, c("y", "y_L1", "y_L2", "y_L3", "y_L4")]), ]
  ar_fit_oos <- tryCatch(
    lm(y ~ y_L1 + y_L2 + y_L3 + y_L4, data = ar_train),
    error = function(e) NULL
  )

  if (!is.null(ar_fit_oos)) {
    # Build the AR prediction for the forecast origin
    y_vals <- df$y[(t_idx - 1):(t_idx - 4)]
    if (all(!is.na(y_vals))) {
      newdata_ar <- data.frame(y_L1 = y_vals[1], y_L2 = y_vals[2],
                               y_L3 = y_vals[3], y_L4 = y_vals[4])
      pred_ar <- predict(ar_fit_oos, newdata = newdata_ar)
      errors_ar[i] <- actual - pred_ar
    } else {
      errors_ar[i] <- NA
    }
  } else {
    errors_ar[i] <- NA
  }
}

# Compute MSFEs (mean squared forecast errors)
msfe_pc <- mean(errors_pc^2, na.rm = TRUE)
msfe_ar <- mean(errors_ar^2, na.rm = TRUE)
relative_msfe <- msfe_pc / msfe_ar

cat(sprintf("\n--- Forecast Comparison Results ---\n"))
cat(sprintf("Phillips Curve MSFE:  %.4f\n", msfe_pc))
cat(sprintf("AR(4) Benchmark MSFE: %.4f\n", msfe_ar))
cat(sprintf("Relative MSFE (PC/AR): %.4f\n", relative_msfe))
cat(ifelse(relative_msfe < 1,
    "Phillips curve outperforms the AR benchmark.\n",
    "AR benchmark outperforms the Phillips curve.\n"))
cat("\n")

# --- Diebold-Mariano test for equal predictive accuracy ---
# H0: E[e_pc^2] = E[e_ar^2]
# This is a formal hypothesis test on forecast loss differentials
d <- errors_pc^2 - errors_ar^2
d <- d[!is.na(d)]
dm_stat <- mean(d) / (sd(d) / sqrt(length(d)))
dm_pval <- 2 * (1 - pnorm(abs(dm_stat)))
cat(sprintf("Diebold-Mariano test: DM = %.3f, p-value = %.4f\n", dm_stat, dm_pval))
cat(ifelse(dm_pval < 0.05,
    "Reject H0: forecast accuracy differs significantly.\n\n",
    "Fail to reject H0: no significant difference in accuracy.\n\n"))


# ==============================================================================
# PART 6: RESIDUAL DIAGNOSTICS
# Topics: Week 8 (chi-square, K-S, Q-Q), Weeks 10-12 (regression diagnostics)
# ==============================================================================

cat("=== PART 6: Residual Diagnostics ===\n\n")

resid_pc <- residuals(pc_ols)

# --- Normality tests on residuals ---
cat("Jarque-Bera test on OLS residuals:\n")
jb_resid <- jarque.bera.test(ts(resid_pc))
print(jb_resid)

cat("\nKolmogorov-Smirnov test on OLS residuals:\n")
ks_resid <- ks.test(resid_pc, "pnorm", mean(resid_pc), sd(resid_pc))
print(ks_resid)

# --- Q-Q plot of residuals ---
p3 <- ggplot(data.frame(r = resid_pc), aes(sample = r)) +
  stat_qq() + stat_qq_line(color = "red") +
  labs(title = "Q-Q Plot of Phillips Curve Residuals",
       x = "Theoretical Quantiles", y = "Sample Quantiles") +
  theme_minimal()
print(p3)

# --- Breusch-Godfrey test for serial correlation in residuals ---
cat("\nBreusch-Godfrey test for serial correlation (12 lags):\n")
bg_test <- bgtest(pc_ols, order = 12)
print(bg_test)

# --- Breusch-Pagan test for heteroskedasticity ---
cat("\nBreusch-Pagan test for heteroskedasticity:\n")
bp_test <- bptest(pc_ols)
print(bp_test)


# ==============================================================================
# PART 7: STRUCTURAL STABILITY (Chow test / QLR)
# Topics: Week 7 (F-test), regression diagnostics
# S&W find evidence of parameter instability in the Phillips curve
# ==============================================================================

cat("\n=== PART 7: Structural Stability ===\n\n")

# Simple Chow test: split sample at a candidate break date
# Common choices: 1984:01 (Great Moderation), 1990:01, etc.
break_date <- as.Date("1984-01-01")

pre  <- pc_data[pc_data$date < break_date, ]
post <- pc_data[pc_data$date >= break_date, ]

fit_pre  <- lm(pc_formula, data = pre)
fit_post <- lm(pc_formula, data = post)

# Chow F-statistic
SSR_full <- sum(residuals(pc_ols)^2)
SSR_pre  <- sum(residuals(fit_pre)^2)
SSR_post <- sum(residuals(fit_post)^2)
k <- length(coef(pc_ols))
n_full <- nrow(pc_data)

chow_F <- ((SSR_full - SSR_pre - SSR_post) / k) / ((SSR_pre + SSR_post) / (n_full - 2 * k))
chow_p <- 1 - pf(chow_F, k, n_full - 2 * k)

cat(sprintf("Chow test (break at %s):\n", break_date))
cat(sprintf("  F = %.3f, df1 = %d, df2 = %d, p-value = %.4f\n",
            chow_F, k, n_full - 2 * k, chow_p))
cat(ifelse(chow_p < 0.05,
    "  Reject H0: evidence of structural break.\n\n",
    "  Fail to reject H0: no evidence of structural break.\n\n"))


# ==============================================================================
# PART 8: EXTENSION — UPDATED SAMPLE
# Re-run the forecast comparison on post-2000 data to see how the Phillips
# curve has performed more recently (it has weakened — the "flattening" result)
# ==============================================================================

cat("=== PART 8: Extension — Post-2000 Forecast Performance ===\n\n")

post2000_rows <- which(df$date >= as.Date("2000-01-01") & !is.na(df$y))

if (length(post2000_rows) > 0) {
  post_errors_pc <- errors_pc[forecast_rows %in% post2000_rows]
  post_errors_ar <- errors_ar[forecast_rows %in% post2000_rows]

  # If the indices don't align perfectly, recompute for the subset
  idx_post <- which(forecast_rows %in% post2000_rows)
  if (length(idx_post) > 10) {
    msfe_pc_post <- mean(errors_pc[idx_post]^2, na.rm = TRUE)
    msfe_ar_post <- mean(errors_ar[idx_post]^2, na.rm = TRUE)
    cat(sprintf("Post-2000 Phillips Curve MSFE:  %.4f\n", msfe_pc_post))
    cat(sprintf("Post-2000 AR(4) Benchmark MSFE: %.4f\n", msfe_ar_post))
    cat(sprintf("Post-2000 Relative MSFE (PC/AR): %.4f\n", msfe_pc_post / msfe_ar_post))
    cat("\nThis lets you discuss why the Phillips curve's forecasting power may\n")
    cat("have weakened: anchored inflation expectations, flattening of the\n")
    cat("wage-price relationship, and structural changes in the labor market.\n\n")
  }
}


# ==============================================================================
# SUMMARY: COURSE TOPIC MAPPING
# ==============================================================================

cat("
================================================================================
COURSE TOPIC MAPPING — Where each part connects to MATH 1220
================================================================================

Part 1 (Descriptive Stats)
  -> Weeks 1-4: sample mean, variance, SE, histograms, EDF
  -> Week 8:    K-S test, chi-square goodness-of-fit, Q-Q plots

Part 2 (Phillips Curve OLS)
  -> Weeks 10-12: linear regression, t-tests, F-tests, confidence intervals
  -> Week 3:      delta method (for NAIRU confidence interval)
  -> Week 4:      t and F distributions

Part 3 (MLE & LRT)
  -> Weeks 5-6: MLE, sufficiency (OLS = MLE under normality)
  -> Week 7:    likelihood ratio test, AIC/BIC as decision-theoretic criteria
  -> Week 6:    Cramer-Rao (MLE achieves the bound under regularity conditions)

Part 4 (AR/ARIMA)
  -> Weeks 12-13: AR, MA, ARIMA, ACF/PACF, model identification

Part 5 (Out-of-Sample Evaluation)
  -> Week 6:  MSE as a loss function (decision theory)
  -> Week 7:  hypothesis test on forecast accuracy (Diebold-Mariano)

Part 6 (Diagnostics)
  -> Week 8:  chi-square test, K-S test, Q-Q plots on residuals
  -> Week 10: regression diagnostics, heteroskedasticity, serial correlation

Part 7 (Structural Stability)
  -> Week 7:  Chow test is an F-test; connects to UMP testing framework

Part 8 (Extension)
  -> Original contribution — updated data, economic interpretation
================================================================================
\n")

cat("Done! See plots in your R graphics window.\n")