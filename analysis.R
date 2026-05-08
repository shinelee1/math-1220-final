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

# Run once to install dependencies, then leave commented out:
# install.packages(c("fredr", "forecast", "lmtest", "sandwich", "tseries",
#                    "ggplot2", "modelsummary", "gt"))

library(fredr)        # FRED API access
library(forecast)     # ARIMA, auto.arima, ggAcf, ggPacf
library(lmtest)       # coeftest, bptest, bgtest
library(sandwich)     # HAC standard errors (NeweyWest)
library(tseries)      # adf.test, jarque.bera.test
library(ggplot2)      # plotting
library(gt)           # HTML table rendering (required by modelsummary HTML output)
library(modelsummary) # clean regression / summary tables

# FIX: explicitly load project .Renviron so FRED_API_KEY is available in any
# R session, not just ones launched via RStudio's project loader.
if (file.exists(".Renviron")) readRenviron(".Renviron")

# Set FRED API key from environment variable — never hard-code the key
fredr_set_key(Sys.getenv("FRED_API_KEY"))

# FIX: create output directories up front so ggsave / write.csv never fail
dir.create("outputs/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/tables",  recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# 1. DATA: Pulling data from FRED
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
# 2. VARIABLE CONSTRUCTION
# ==============================================================================
# Stock & Watson define:
#   pi_t     = 1200 * ln(P_t / P_{t-1})          (monthly inflation, annualized)
#   pi^h_t   = (1200/h) * ln(P_t / P_{t-h})      (h-month inflation, annualized)
#   Dependent variable: pi^h_{t+h} - pi_t
#
# We set h = 12 (12-month forecast horizon)

h <- 12  # forecast horizon in months

# Merge datasets on date
df <- merge(
  data.frame(date = cpi_all$date,  cpi      = cpi_all$value),
  data.frame(date = unemp$date,    urate    = unemp$value),
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
df$pi      <- c(NA, 1200 * diff(log(df$cpi)))       # monthly annualized inflation
df$pi_core <- c(NA, 1200 * diff(log(df$cpi_core)))  # core monthly annualized inflation

# FIX: vectorized pi_h construction (avoids loop and is easier to verify).
# diff(log(x), lag = h)[j] = log(x[j+h]) - log(x[j]), so prepending h NAs
# places pi^h_t = (1200/h)*log(P_t/P_{t-h}) at row t, for t = h+1, ..., n.
n      <- nrow(df)   # row count BEFORE trimming — needed for y construction below
df$pi_h <- c(rep(NA, h), (1200 / h) * diff(log(df$cpi), lag = h))

# FIX: vectorized y construction.
# y[t] = pi^h_{t+h} - pi_t requires leading pi_h by h periods.
# The last h rows cannot have a valid y (no future P data), so they get NA.
df$y <- c(df$pi_h[(h + 1):n] - df$pi[1:(n - h)], rep(NA, h))

# Change in monthly inflation (lagged regressor in Phillips curve)
df$dpi <- c(NA, diff(df$pi))

# Trim to paper sample starting 1959:01 and reset row names
df <- df[df$date >= as.Date("1959-01-01"), ]
rownames(df) <- NULL

# FIX: update n AFTER trimming so lag construction below uses the correct length
n <- nrow(df)

cat("Dataset constructed:", n, "observations from",
    as.character(min(df$date)), "to", as.character(max(df$date)), "\n\n")


# ==============================================================================
# PART 1: DESCRIPTIVE STATISTICS & NONPARAMETRIC TESTS
# Topics from class: Weeks 1-4 (probability review, EDF, CLT, histograms, Q-Q plots)
#                    Week 8 (K-S test, chi-square goodness of fit)
# ==============================================================================

cat("=== PART 1: Descriptive Statistics ===\n\n")

vars_of_interest <- c("pi", "urate", "dpi")

# Build summary table with a single lapply — avoids the incremental rbind loop
summary_table <- do.call(rbind, lapply(vars_of_interest, function(v) {
  x <- na.omit(df[[v]])
  data.frame(
    Variable = v,
    Mean     = round(mean(x), 4),
    SD       = round(sd(x), 4),
    Skewness = round(mean((x - mean(x))^3) / sd(x)^3, 4),
    N        = length(x)
  )
}))

print(summary_table, row.names = FALSE)

# FIX: modelsummary() expects model objects (lm, glm, …), not plain data frames.
# datasummary_df() is the correct modelsummary function for tabulating a data frame.
datasummary_df(
  summary_table,
  output = "outputs/tables/descriptive_statistics.html"
)

# Histogram of monthly inflation
p1 <- ggplot(df[!is.na(df$pi), ], aes(x = pi)) +
  geom_histogram(
    aes(y = after_stat(density)),
    bins = 50, fill = "steelblue", alpha = 0.7
  ) +
  stat_function(
    fun  = dnorm,
    args = list(mean = mean(df$pi, na.rm = TRUE), sd = sd(df$pi, na.rm = TRUE)),
    color = "red", linewidth = 0.8
  ) +
  labs(
    title = "Distribution of Monthly Inflation (annualized)",
    x = expression(pi[t]),
    y = "Density"
  ) +
  theme_minimal()

ggsave("outputs/figures/inflation_histogram.png", p1, width = 8, height = 5)

# Q-Q plot for inflation changes
p2 <- ggplot(df[!is.na(df$dpi), ], aes(sample = dpi)) +
  stat_qq() +
  stat_qq_line(color = "red") +
  labs(
    title = "Q-Q Plot: Changes in Monthly Inflation",
    x = "Theoretical Quantiles", y = "Sample Quantiles"
  ) +
  theme_minimal()

ggsave("outputs/figures/inflation_qqplot.png", p2, width = 8, height = 5)

# Kolmogorov-Smirnov test for normality of dpi
dpi_clean <- na.omit(df$dpi)
ks_result <- ks.test(dpi_clean, "pnorm", mean(dpi_clean), sd(dpi_clean))
cat("Kolmogorov-Smirnov test for normality of dpi:\n")
cat(sprintf("  D = %.4f, p-value = %.4f\n\n", ks_result$statistic, ks_result$p.value))

# Jarque-Bera test (uses skewness & kurtosis — connects to moment estimators!)
jb_result <- jarque.bera.test(ts(dpi_clean))
cat("Jarque-Bera test for normality of dpi:\n")
cat(sprintf("  JB = %.2f, p-value = %.6f\n\n", jb_result$statistic, jb_result$p.value))

# Chi-square goodness-of-fit test (Pearson chi-square statistic from Week 8).
# FIX: renamed k -> k_bins to avoid collision with k (# parameters) in Part 7.
k_bins  <- 10
breaks  <- qnorm(seq(0, 1, length.out = k_bins + 1), mean(dpi_clean), sd(dpi_clean))
breaks[1] <- -Inf; breaks[k_bins + 1] <- Inf
observed  <- table(cut(dpi_clean, breaks = breaks))
expected  <- rep(length(dpi_clean) / k_bins, k_bins)
chi2_stat <- sum((as.numeric(observed) - expected)^2 / expected)
chi2_df   <- k_bins - 1 - 2   # subtract 2 for estimated mu, sigma
chi2_pval <- 1 - pchisq(chi2_stat, df = chi2_df)
cat("Chi-square goodness-of-fit test (k=10 bins, normal null):\n")
cat(sprintf("  chi2 = %.2f, df = %d, p-value = %.4f\n\n", chi2_stat, chi2_df, chi2_pval))

# Save normality test results to disk
norm_tests <- data.frame(
  Test      = c("KS (dpi)", "Jarque-Bera (dpi)", "Chi-sq GoF (dpi, k=10)"),
  Statistic = c(round(ks_result$statistic, 4),
                round(jb_result$statistic, 2),
                round(chi2_stat, 2)),
  p_value   = c(round(ks_result$p.value, 4),
                round(jb_result$p.value, 6),
                round(chi2_pval, 4))
)
write.csv(norm_tests, "outputs/tables/normality_tests.csv", row.names = FALSE)


# ==============================================================================
# PART 2: PHILLIPS CURVE ESTIMATION VIA OLS
# Topics from class: Weeks 10-12 (linear regression, confidence intervals, t-tests, F-tests)
#                    Week 3 (delta method for NAIRU CI)
# ==============================================================================

cat("=== PART 2: Phillips Curve OLS Estimation ===\n\n")

# FIX: construct lags using the post-trim n so indices are correct.
# seq_len(n - j) is safer than 1:(n-j) when j could equal n.
max_lag <- 12

for (j in 0:max_lag) {
  df[[paste0("urate_L", j)]] <- c(rep(NA, j), df$urate[seq_len(n - j)])
  if (j > 0) {
    df[[paste0("dpi_L", j)]] <- c(rep(NA, j), df$dpi[seq_len(n - j)])
  }
}

# Baseline Phillips curve: 4 lags of u, 4 lags of dpi
# (S&W use various lag lengths; 4 is a standard starting point)
pc_formula <- as.formula(
  paste("y ~", paste0("urate_L", 0:3, collapse = " + "), "+",
        paste0("dpi_L", 1:4, collapse = " + "))
)

# Estimate on full available sample (complete cases only)
pc_data <- df[complete.cases(df[, c("y", paste0("urate_L", 0:3),
                                        paste0("dpi_L", 1:4))]), ]
pc_ols  <- lm(pc_formula, data = pc_data)

cat("--- Baseline Phillips Curve (4 lags u, 4 lags dpi) ---\n")
# FIX: summary() is invisible in scripts — must wrap with print()
print(summary(pc_ols))

# HAC (Newey-West) standard errors — h-step overlapping errors induce MA(h-1) structure
cat("\n--- HAC (Newey-West) Standard Errors ---\n")
hac_se <- coeftest(pc_ols, vcov = NeweyWest(pc_ols, lag = h - 1))
print(hac_se)

# Save regression table with HAC SEs to HTML
modelsummary(
  list("Phillips Curve (OLS)" = pc_ols),
  vcov   = NeweyWest(pc_ols, lag = h - 1),
  stars  = TRUE,
  output = "outputs/tables/phillips_curve_ols.html"
)

# F-test: joint significance of unemployment lags
# H0: all beta_urate = 0 (unemployment has no predictive power for inflation)
pc_restricted <- lm(
  as.formula(paste("y ~", paste0("dpi_L", 1:4, collapse = " + "))),
  data = pc_data
)
f_test <- anova(pc_restricted, pc_ols)
cat("\n--- F-test: Joint significance of unemployment lags ---\n")
print(f_test)
cat(sprintf("F = %.3f, p-value = %.4f\n\n", f_test$F[2], f_test$`Pr(>F)`[2]))

# Confidence intervals for coefficients (95%)
cat("--- 95% Confidence Intervals for OLS Coefficients ---\n")
print(confint(pc_ols, level = 0.95))
cat("\n")

# Delta method for implied NAIRU
# If the Phillips curve is: y = alpha + beta(u)*u + gamma(L)*dpi + e
# Setting y = 0 and dpi changes = 0: u_nairu = -alpha / sum(beta_j)
# This is a nonlinear function of parameters -> use delta method (Week 3)
coefs     <- coef(pc_ols)
alpha_hat <- coefs["(Intercept)"]
beta_sum  <- sum(coefs[grep("urate", names(coefs))])
nairu_hat <- -alpha_hat / beta_sum
cat(sprintf("Implied NAIRU = %.2f%%\n", nairu_hat))

# Delta method SE: g(theta) = -alpha/sum(beta)
# Gradient: d/d(alpha) = -1/sum(beta); d/d(beta_j) = alpha/sum(beta)^2
V    <- vcov(pc_ols)
grad <- rep(0, length(coefs))
grad[1] <- -1 / beta_sum                            # d/d(alpha)
for (j in grep("urate", names(coefs))) {
  grad[j] <- alpha_hat / beta_sum^2                 # d/d(beta_j)
}
# FIX: t(grad) %*% V %*% grad returns a 1x1 matrix; as.numeric() extracts the scalar
nairu_se <- sqrt(as.numeric(t(grad) %*% V %*% grad))
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
  data = pc_data  # same sample as full model, valid LRT comparison
)
loglik_restricted <- logLik(pc_2lag)
cat(sprintf("Log-likelihood (restricted, 2 lags): %.2f\n", loglik_restricted))

# LRT statistic: Lambda = -2 * (l_restricted - l_full) ~ chi-sq(df)
lrt_stat <- -2 * (as.numeric(loglik_restricted) - as.numeric(loglik_full))
lrt_df   <- length(coef(pc_ols)) - length(coef(pc_2lag))
lrt_pval <- 1 - pchisq(lrt_stat, df = lrt_df)
cat(sprintf("LRT statistic = %.3f, df = %d, p-value = %.4f\n", lrt_stat, lrt_df, lrt_pval))
cat("Interpretation:", ifelse(lrt_pval < 0.05,
    "Reject H0 — additional lags improve the model.",
    "Fail to reject H0 — parsimonious model sufficient."), "\n\n")

# --- Model selection via AIC and BIC ---
cat("--- Model Selection (AIC / BIC) ---\n")
cat("AIC and BIC penalize complexity differently. AIC minimizes expected\n")
cat("Kullback-Leibler divergence (a specific loss function from decision theory).\n\n")

lag_specs <- list(
  "2 lags" = list(u = 0:1, dpi = 1:2),
  "4 lags" = list(u = 0:3, dpi = 1:4),
  "6 lags" = list(u = 0:5, dpi = 1:6),
  "8 lags" = list(u = 0:7, dpi = 1:8)
)

aic_bic_table <- do.call(rbind, lapply(names(lag_specs), function(spec_name) {
  spec <- lag_specs[[spec_name]]
  f <- as.formula(paste("y ~",
    paste0("urate_L", spec$u, collapse = " + "), "+",
    paste0("dpi_L", spec$dpi, collapse = " + ")))
  fit <- lm(f, data = pc_data)
  data.frame(Spec = spec_name, AIC = round(AIC(fit), 2), BIC = round(BIC(fit), 2))
}))

print(aic_bic_table, row.names = FALSE)
write.csv(aic_bic_table, "outputs/tables/aic_bic_comparison.csv", row.names = FALSE)
cat("\n")


# ==============================================================================
# PART 4: AR BENCHMARK & TIME SERIES MODELS
# Topics: Weeks 12-13 (AR, MA, ARIMA, ACF/PACF)
# ==============================================================================

cat("=== PART 4: AR Benchmark & Time Series Models ===\n\n")

# --- AR(p) benchmark for inflation changes ---
# Stock & Watson compare Phillips curve against a univariate AR model
# AR order selected by AIC (as in the paper)

y_ts <- ts(na.omit(df$y), frequency = 12)

ar_fit <- ar(y_ts, order.max = 12, method = "ols", aic = TRUE)
cat(sprintf("AR benchmark: order selected by AIC = %d\n", ar_fit$order))
cat("AR coefficients:\n")
print(round(ar_fit$ar, 4))
cat("\n")

# --- ARIMA via auto.arima ---
pi_ts     <- ts(na.omit(df$pi), frequency = 12, start = c(1959, 2))
arima_fit <- auto.arima(pi_ts, max.p = 6, max.q = 6, max.d = 2,
                         seasonal = FALSE, stepwise = FALSE, approximation = FALSE)
cat("auto.arima selected model:\n")
print(arima_fit)
cat(sprintf("\nAIC = %.2f, BIC = %.2f\n\n", AIC(arima_fit), BIC(arima_fit)))

# FIX: use ggAcf / ggPacf (from the forecast package) so plots can be saved
# with ggsave. The original par(mfrow=...) approach could not be saved to file.
p_acf  <- ggAcf(na.omit(df$pi), lag.max = 36) +
  ggtitle("ACF of Monthly Inflation") + theme_minimal()
p_pacf <- ggPacf(na.omit(df$pi), lag.max = 36) +
  ggtitle("PACF of Monthly Inflation") + theme_minimal()

ggsave("outputs/figures/inflation_acf.png",  p_acf,  width = 8, height = 4)
ggsave("outputs/figures/inflation_pacf.png", p_pacf, width = 8, height = 4)
cat("ACF and PACF plots saved.\n\n")


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

errors_pc <- numeric(length(forecast_rows))
errors_ar <- numeric(length(forecast_rows))

cat(sprintf("Running pseudo out-of-sample evaluation (%d forecast origins)...\n",
            length(forecast_rows)))

for (i in seq_along(forecast_rows)) {
  t_idx  <- forecast_rows[i]
  actual <- df$y[t_idx]

  # Training data: all complete-case observations strictly before the forecast origin
  train <- df[seq_len(t_idx - 1), ]
  train <- train[complete.cases(train[, c("y", paste0("urate_L", 0:3),
                                               paste0("dpi_L", 1:4))]), ]

  # --- Phillips curve forecast ---
  pc_fit <- tryCatch(lm(pc_formula, data = train), error = function(e) NULL)
  if (!is.null(pc_fit)) {
    errors_pc[i] <- actual - predict(pc_fit, newdata = df[t_idx, , drop = FALSE])
  } else {
    errors_pc[i] <- NA
  }

  # --- AR(4) benchmark forecast ---
  # FIX: build lags directly from the y column of the training data frame instead
  # of appending columns with manual c(NA, ...) indexing, which was error-prone.
  y_tr <- train$y   # already NA-free because train is filtered on complete.cases(y)
  m    <- length(y_tr)

  ar_fit_oos <- NULL
  if (m >= 5) {
    ar_df <- data.frame(
      y    = y_tr[5:m],
      y_L1 = y_tr[4:(m - 1)],
      y_L2 = y_tr[3:(m - 2)],
      y_L3 = y_tr[2:(m - 3)],
      y_L4 = y_tr[1:(m - 4)]
    )
    ar_df <- ar_df[complete.cases(ar_df), ]
    ar_fit_oos <- tryCatch(
      lm(y ~ y_L1 + y_L2 + y_L3 + y_L4, data = ar_df),
      error = function(e) NULL
    )
  }

  if (!is.null(ar_fit_oos)) {
    # Most recent 4 y-values at the forecast origin (descending index = correct lag order)
    y_recent <- df$y[(t_idx - 1):(t_idx - 4)]
    if (length(y_recent) == 4 && !anyNA(y_recent)) {
      newdata_ar   <- data.frame(y_L1 = y_recent[1], y_L2 = y_recent[2],
                                  y_L3 = y_recent[3], y_L4 = y_recent[4])
      errors_ar[i] <- actual - predict(ar_fit_oos, newdata = newdata_ar)
    } else {
      errors_ar[i] <- NA
    }
  } else {
    errors_ar[i] <- NA
  }
}

# Compute MSFEs (mean squared forecast errors)
msfe_pc       <- mean(errors_pc^2, na.rm = TRUE)
msfe_ar       <- mean(errors_ar^2, na.rm = TRUE)
relative_msfe <- msfe_pc / msfe_ar

cat(sprintf("\n--- Forecast Comparison Results ---\n"))
cat(sprintf("Phillips Curve MSFE:   %.4f\n", msfe_pc))
cat(sprintf("AR(4) Benchmark MSFE:  %.4f\n", msfe_ar))
cat(sprintf("Relative MSFE (PC/AR): %.4f\n", relative_msfe))
cat(ifelse(relative_msfe < 1,
    "Phillips curve outperforms the AR benchmark.\n",
    "AR benchmark outperforms the Phillips curve.\n"))
cat("\n")

# --- Diebold-Mariano test for equal predictive accuracy ---
# H0: E[e_pc^2] = E[e_ar^2]
# This is a formal hypothesis test on forecast loss differentials
d       <- errors_pc^2 - errors_ar^2
d       <- d[!is.na(d)]
dm_stat <- mean(d) / (sd(d) / sqrt(length(d)))
dm_pval <- 2 * (1 - pnorm(abs(dm_stat)))
cat(sprintf("Diebold-Mariano test: DM = %.3f, p-value = %.4f\n", dm_stat, dm_pval))
cat(ifelse(dm_pval < 0.05,
    "Reject H0: forecast accuracy differs significantly.\n\n",
    "Fail to reject H0: no significant difference in accuracy.\n\n"))

# Save forecast comparison results
forecast_results <- data.frame(
  Model         = c("Phillips Curve", "AR(4) Benchmark"),
  MSFE          = round(c(msfe_pc, msfe_ar), 4),
  Relative_MSFE = round(c(relative_msfe, 1), 4),
  DM_stat       = c(round(dm_stat, 3), NA),
  DM_pval       = c(round(dm_pval, 4), NA)
)
write.csv(forecast_results, "outputs/tables/forecast_comparison.csv", row.names = FALSE)

# Save forecast errors over time as a plot
forecast_dates <- df$date[forecast_rows]
p_errors <- ggplot(
  data.frame(
    date  = rep(forecast_dates, 2),
    error = c(errors_pc, errors_ar),
    model = rep(c("Phillips Curve", "AR(4)"), each = length(forecast_rows))
  ),
  aes(x = date, y = error, color = model)
) +
  geom_line(alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(title = "Out-of-Sample Forecast Errors (h = 12 months)",
       x = NULL, y = "Forecast Error", color = "Model") +
  theme_minimal()

ggsave("outputs/figures/forecast_errors.png", p_errors, width = 10, height = 5)
cat("Forecast error plot saved.\n\n")


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
# FIX: p3 was never saved in the original; now saved alongside the other figures
p3 <- ggplot(data.frame(r = resid_pc), aes(sample = r)) +
  stat_qq() + stat_qq_line(color = "red") +
  labs(title = "Q-Q Plot of Phillips Curve Residuals",
       x = "Theoretical Quantiles", y = "Sample Quantiles") +
  theme_minimal()

ggsave("outputs/figures/residual_qqplot.png", p3, width = 8, height = 5)

# --- Breusch-Godfrey test for serial correlation in residuals ---
cat("\nBreusch-Godfrey test for serial correlation (12 lags):\n")
bg_test <- bgtest(pc_ols, order = 12)
print(bg_test)

# --- Breusch-Pagan test for heteroskedasticity ---
cat("\nBreusch-Pagan test for heteroskedasticity:\n")
bp_test <- bptest(pc_ols)
print(bp_test)

# Save residual diagnostic test summary
resid_tests <- data.frame(
  Test      = c("Jarque-Bera", "KS normality", "Breusch-Godfrey (12 lags)",
                "Breusch-Pagan"),
  Statistic = c(round(jb_resid$statistic, 3), round(ks_resid$statistic, 4),
                round(bg_test$statistic, 3),   round(bp_test$statistic, 3)),
  p_value   = c(round(jb_resid$p.value, 6), round(ks_resid$p.value, 4),
                round(bg_test$p.value, 4),   round(bp_test$p.value, 4))
)
write.csv(resid_tests, "outputs/tables/residual_diagnostics.csv", row.names = FALSE)
cat("\nResidual diagnostic results saved.\n\n")


# ==============================================================================
# PART 7: STRUCTURAL STABILITY (Chow test / QLR)
# Topics: Week 7 (F-test), regression diagnostics
# S&W find evidence of parameter instability in the Phillips curve
# ==============================================================================

cat("=== PART 7: Structural Stability ===\n\n")

# Simple Chow test: split sample at a candidate break date
# Common choices: 1984:01 (Great Moderation), 1990:01, etc.
break_date <- as.Date("1984-01-01")

pre  <- pc_data[pc_data$date <  break_date, ]
post <- pc_data[pc_data$date >= break_date, ]

fit_pre  <- lm(pc_formula, data = pre)
fit_post <- lm(pc_formula, data = post)

# Chow F-statistic
SSR_full <- sum(residuals(pc_ols)^2)
SSR_pre  <- sum(residuals(fit_pre)^2)
SSR_post <- sum(residuals(fit_post)^2)
# FIX: renamed k -> k_pc to avoid collision with k_bins defined in Part 1
k_pc   <- length(coef(pc_ols))
n_full <- nrow(pc_data)

chow_F <- ((SSR_full - SSR_pre - SSR_post) / k_pc) /
          ((SSR_pre + SSR_post) / (n_full - 2 * k_pc))
chow_p <- 1 - pf(chow_F, k_pc, n_full - 2 * k_pc)

cat(sprintf("Chow test (break at %s):\n", break_date))
cat(sprintf("  F = %.3f, df1 = %d, df2 = %d, p-value = %.4f\n",
            chow_F, k_pc, n_full - 2 * k_pc, chow_p))
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
# FIX: removed unused post_errors_pc / post_errors_ar intermediaries from
# the original; compute sub-period MSFEs directly from the stored error vectors.
idx_post <- which(forecast_rows %in% post2000_rows)

if (length(idx_post) > 10) {
  msfe_pc_post <- mean(errors_pc[idx_post]^2, na.rm = TRUE)
  msfe_ar_post <- mean(errors_ar[idx_post]^2, na.rm = TRUE)
  cat(sprintf("Post-2000 Phillips Curve MSFE:   %.4f\n", msfe_pc_post))
  cat(sprintf("Post-2000 AR(4) Benchmark MSFE:  %.4f\n", msfe_ar_post))
  cat(sprintf("Post-2000 Relative MSFE (PC/AR): %.4f\n\n", msfe_pc_post / msfe_ar_post))
  cat("This lets you discuss why the Phillips curve's forecasting power may\n")
  cat("have weakened: anchored inflation expectations, flattening of the\n")
  cat("wage-price relationship, and structural changes in the labor market.\n\n")
}

cat("=== Analysis complete. ===\n")
cat("Figures saved to: outputs/figures/\n")
cat("Tables  saved to: outputs/tables/\n")
