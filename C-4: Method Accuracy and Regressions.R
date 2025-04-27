# ============================================================================
#  Appendix C:  Beta‑Metric Adjustments, Summary Tables & Overall Regressions
# -----------------------------------------------------------------------------
#  Commands have been reordered for logical flow; no new functions were added.
#  The script:
#     1. Loads master_data.csv.
#     2. Defines helper `calculate_metrics()` (already present) followed by
#        evaluation blocks for RMSE tables, regression grids, SD stats, and a
#        24‑month average‑beta plot.
# ============================================================================

# ----------------------------- 1.  Setup ------------------------------------
library(readr)
library(dplyr)
library(lubridate)
library(tidyr)
library(tidyverse)
library(ggplot2)
library(gridExtra)

setwd("C:/Users/ikolesnikov25/OneDrive - Claremont McKenna College/Desktop/Thesis")

data  <- read_csv("master_data.csv")
gamma <- 1.0

# --------------------------- 2.  Helpers ------------------------------------
# (calculate_metrics and run_regressions are unchanged)

# ---- 2a.  Metric calculator -------------------------------------------------
calculate_metrics <- function(data, gamma, period, suffix) {
  data %>%
    mutate(
      !!paste0("skew_", suffix, "_oic") := !!sym(paste0("skew_", suffix)) +
        gamma * (!!sym(paste0("stock_kurt_", period, "m_past")) - 3) *
        !!sym(paste0("stock_vol_", period, "m_past")),
      !!paste0("skew_", suffix, "_oic_m") := !!sym(paste0("skew_", suffix, "_m")) +
        gamma * (!!sym(paste0("market_kurt_", period, "m_past")) - 3) *
        !!sym(paste0("market_vol_", period, "m_past")),
      
      !!paste0("vol_", suffix, "_oic") := !!sym(paste0("vol_", suffix)) /
        (1 - gamma * !!sym(paste0("stock_skew_", period, "m_past")) *
           !!sym(paste0("stock_vol_", period, "m_past"))),
      !!paste0("vol_", suffix, "_oic_m") := !!sym(paste0("vol_", suffix, "_m")) /
        (1 - gamma * !!sym(paste0("market_skew_", period, "m_past")) *
           !!sym(paste0("market_vol_", period, "m_past"))),
      
      !!paste0("iv_", period, "m") := approx(
        x = c(!!sym(paste0("lower_yte_", suffix)),
              !!sym(paste0("higher_yte_", suffix))),
        y = c(!!sym(paste0("lower_Iv_", suffix)),
              !!sym(paste0("higher_Iv_", suffix))),
        xout = as.numeric(suffix))$y,
      !!paste0("iv_", period, "m_m") := approx(
        x = c(!!sym(paste0("lower_yte_", suffix, "_m")),
              !!sym(paste0("higher_yte_", suffix, "_m"))),
        y = c(!!sym(paste0("lower_Iv_", suffix, "_m")),
              !!sym(paste0("higher_Iv_", suffix, "_m"))),
        xout = as.numeric(suffix))$y,
      
      !!paste0("corr_", suffix, "_oic") :=
        (!!sym(paste0("skew_", suffix, "_oic")) /
           !!sym(paste0("skew_", suffix, "_oic_m")))^(1/3),
      
      beta      = !!sym(paste0("beta_", period, "m_future")),
      beta_h    = !!sym(paste0("beta_", period, "m_past")),
      beta_oi   = (!!sym(paste0("skew_", suffix)) / !!sym(paste0("skew_", suffix, "_m")))^(1/3) *
        !!sym(paste0("vol_", suffix)) / !!sym(paste0("vol_", suffix, "_m")),
      beta_c    = !!sym(paste0("corr_", period, "m_past")) *
        !!sym(paste0("vol_", suffix)) / !!sym(paste0("vol_", suffix, "_m")),
      beta_atm  = !!sym(paste0("corr_", period, "m_past")) *
        !!sym(paste0("iv_", period, "m")) / !!sym(paste0("iv_", period, "m_m")),
      beta_oic  = (!!sym(paste0("skew_", suffix, "_oic")) /
                     !!sym(paste0("skew_", suffix, "_oic_m")))^(1/3) *
        !!sym(paste0("vol_", suffix, "_oic")) /
        !!sym(paste0("vol_", suffix, "_oic_m")),
      beta_oic_2= !!sym(paste0("corr_", suffix, "_oic")) *
        !!sym(paste0("iv_", period, "m")) / !!sym(paste0("iv_", period, "m_m"))
    ) %>%
    filter(!!sym(paste0("num_obs_higher_", suffix))  > 4,
           !!sym(paste0("num_obs_lower_",  suffix))  > 4,
           !!sym(paste0("num_obs_higher_", suffix, "_m")) > 4,
           !!sym(paste0("num_obs_lower_",  suffix, "_m")) > 4)
}

# ---- 2b.  Regression‑table helper ------------------------------------------
run_regressions <- function(df, suffix) {   # unchanged content inside
  model_formulas <- list(
    "beta_h + beta_oi"                    = beta ~ beta_h + beta_oi,
    "beta_h + beta_c"                     = beta ~ beta_h + beta_c,
    "beta_c + beta_oi"                    = beta ~ beta_c + beta_oi,
    paste0("beta_oi + num_obs_higher_", suffix) %>% as.formula(),
    "beta_h"                              = beta ~ beta_h,
    "beta_oi"                             = beta ~ beta_oi,
    "beta_c"                              = beta ~ beta_c)
  
  results <- map_dfr(names(model_formulas), function(nm) {
    fit <- lm(model_formulas[[nm]], data = df)
    co   <- summary(fit)$coefficients
    preds <- setdiff(rownames(co), "(Intercept)")
    tibble(
      regression      = nm,
      r_squared       = summary(fit)$r.squared,
      adj_r_squared   = summary(fit)$adj.r.squared,
      intercept_coef  = co["(Intercept)", "Estimate"],
      intercept_t     = co["(Intercept)", "t value"],
      intercept_se    = co["(Intercept)", "Std. Error"],
      x1              = preds[1],
      x1_coef         = co[preds[1], "Estimate"],
      x1_t            = co[preds[1], "t value"],
      x1_se           = co[preds[1], "Std. Error"],
      x2              = preds[2] %||% NA,
      x2_coef         = co[preds[2], "Estimate"] %||% NA,
      x2_t            = co[preds[2], "t value"]   %||% NA,
      x2_se           = co[preds[2], "Std. Error"] %||% NA)
  }) %>% arrange(desc(adj_r_squared))
}

# --------------------------- 3.  RMSE Tables --------------------------------
beta_methods <- c("beta_h", "beta_oi", "beta_c", "beta_atm", "beta_oic", "beta_oic_2")
periods <- list(c(6, "0.5"), c(12, "1"), c(24, "2"))

walk(periods, function(p) {
  period_num <- p[[1]]; suf <- p[[2]]
  calc <- calculate_metrics(data, gamma, period_num, suf)
  tbl  <- calc %>%
    pivot_longer(all_of(beta_methods), names_to = "Method", values_to = "Estimate") %>%
    summarise(across(beta:Estimate, ~ .x, .names = "{col}"), .by = Method) %>%
    mutate(RMSE           = sqrt(mean((beta - Estimate)^2, na.rm = TRUE)),
           Avg_Error      = mean(beta - Estimate, na.rm = TRUE),
           Avg_Abs_Error  = mean(abs(beta - Estimate), na.rm = TRUE),
           Pearson_Corr   = cor(beta, Estimate, use = "complete.obs"),
           Spearman_Corr  = cor(beta, Estimate, method = "spearman", use = "complete.obs")) %>%
    select(Method, RMSE, Avg_Error, Avg_Abs_Error, Pearson_Corr, Spearman_Corr) %>%
    arrange(RMSE)
  write_csv(tbl, paste0(period_num, "m_results_overall.csv"))
})

# --------------------------- 4.  Overall Regressions ------------------------
all_reg <- map_dfr(periods, function(p) {
  dat <- calculate_metrics(data, gamma, p[[1]], p[[2]])
  run_regressions(dat, p[[2]]) %>% mutate(period = p[[1]])
})
write_csv(all_reg, "overall_regression_results.csv")

# --------------------------- 5.  SD & Beta Plot -----------------------------
# 24‑month overall beta SD ----------------------------------------------------
calc_24 <- calculate_metrics(data, gamma, 24, "2")
cat("SD of 24‑month beta: ", sd(calc_24$beta, na.rm = TRUE), "\n")

avg_beta <- calc_24 %>% group_by(date) %>% summarise(avg_beta = mean(beta, na.rm = TRUE))

ggplot(avg_beta, aes(date, avg_beta)) +
  geom_line() +
  labs(title = "Average 24‑Month Beta Over Time", x = "Date", y = "Average Beta") +
  theme_minimal()
# ===========================================================================
