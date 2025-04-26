# ============================================================================
#  Appendix C: CAPM Error Term Skewness (R)
# -----------------------------------------------------------------------------
#  Objective
#  ---------
#  •  Draw 1 000 random (ticker, end‑date) pairs across the S&P 500 universe
#     (2005‑01‑01 .. today).
#  •  For each pair, run an OLS regression of daily stock returns on market
#     returns over the previous 126 trading days (≈ 6 calendar months):
#           r_i,t  = α_i  + β_i r_m,t  + ε_i,t.
#  •  Store in two outputs:
#       1.  `half_year_regression_results.csv`   – one row per regression with
#            β̂,  s.e.(β̂),  R² and descriptive moments of r_m, r_i, ε.
#       2.  `half_year_regression_residuals.csv` – long panel of *all* residuals
#            (ε_i,t) to support the skewness histogram in Figure 5.
#  •  Additionally, run D’Agostino skewness tests on each residual window and
#     tabulate counts by significance sign (negative, positive, none).
#
#  Usage
#  -----
#  Source this script after the price‑pull libraries are installed.  Runtime is
#  ≈ 2–3 minutes with a decent network connection (Yahoo API) and modern CPU.
# ============================================================================

# ----------------------------- 1.  Setup ------------------------------------
library(tidyquant)   # price retrieval
library(dplyr)
library(lubridate)
library(purrr)
library(broom)       # tidy()/glance()
library(moments)     # agostino.test()
library(readr)

setwd("/Users/ivankolesnikov/Desktop/Spring 2025/Thesis")  # adjust as needed

# --------------------------- 2.  Universe -----------------------------------
companies <- read_csv("companies.csv", show_col_types = FALSE)
tickers   <- unique(companies$company_name)    # 236 S&P 500 names

start_date <- "2005-01-01"
end_date   <- Sys.Date()
window     <- 126                              # ½‑year trading days (~6m)

# ------------------------- 3.  Fetch Prices ---------------------------------
stock_px  <- tq_get(tickers, from = start_date, to = end_date) %>%
  select(symbol, date, adjusted) %>% arrange(symbol, date)
market_px <- tq_get("^GSPC",  from = start_date, to = end_date) %>%
  select(date, adjusted) %>% rename(market_adj = adjusted)

# ---------------------- 4.  Daily Return Panel ------------------------------
panel <- stock_px %>%
  left_join(market_px, by = "date") %>%
  group_by(symbol) %>% arrange(date, .by_group = TRUE) %>%
  mutate(r_i = adjusted   / lag(adjusted)   - 1,
         r_m = market_adj / lag(market_adj) - 1,
         obs_idx = row_number()) %>%
  ungroup() %>% drop_na(r_i, r_m)

# ------------------- 5.  Candidate End‑Dates --------------------------------
candidates <- panel %>% filter(obs_idx > window) %>%
  select(symbol, date, obs_idx)
set.seed(42)
sampled <- sample_n(candidates, 1000)

# -------------------- 6.  Regression Helper ---------------------------------
run_one <- function(sym, d, idx) {
  dat <- panel %>% filter(symbol == sym,
                          between(obs_idx, idx - window + 1, idx)) %>%
    select(date, r_i, r_m)
  
  fit <- lm(r_i ~ r_m, data = dat)
  
  # --- 6.1  Summary row ------------------------------------------------------
  coef_row <- tidy(fit) %>% filter(term == "r_m")
  g_info   <- glance(fit)
  
  summary_row <- tibble(
    date        = d,
    ticker      = sym,
    beta        = coef_row$estimate,
    se_beta     = coef_row$std.error,
    r2          = g_info$r.squared,
    skew_mkt    = skewness(dat$r_m),
    sd_mkt      = sd(dat$r_m),
    skew_asset  = skewness(dat$r_i),
    sd_asset    = sd(dat$r_i),
    skew_resid  = skewness(residuals(fit)),
    sd_resid    = sd(residuals(fit))
  )
  
  # --- 6.2  Long residual table ---------------------------------------------
  resid_tbl <- tibble(
    ticker    = sym,
    end_date  = d,
    obs_date  = dat$date,
    residual  = residuals(fit)
  )
  
  list(summary = summary_row, residuals = resid_tbl)
}

# ---------------------- 7.  Execute 1 000 Loops -----------------------------
out <- pmap(list(sampled$symbol, sampled$date, sampled$obs_idx), run_one)
results       <- map_dfr(out, "summary")
all_residuals <- map_dfr(out, "residuals")

# ---------------------- 8.  Save Artifacts ----------------------------------
write_csv(results,       "half_year_regression_results.csv")
write_csv(all_residuals, "half_year_regression_residuals.csv")

# ---------------------- 9.  Skewness Summary --------------------------------
agostino_results <- all_residuals %>%
  group_by(ticker, end_date) %>%
  summarise(
    n_obs       = n(),
    skew_sample = agostino.test(residual)$statistic["skew"],
    z_value     = agostino.test(residual)$statistic["z"],
    p_value     = agostino.test(residual)$p.value,
    .groups     = "drop")

skew_summary <- agostino_results %>%
  mutate(skew_class = case_when(
    p_value < 0.05 & skew_sample < 0 ~ "Negative skew (p < .05)",
    p_value < 0.05 & skew_sample > 0 ~ "Positive skew (p < .05)",
    TRUE                              ~ "Not significant")) %>%
  count(skew_class, name = "n_windows")

print(skew_summary)
# ===========================================================================
