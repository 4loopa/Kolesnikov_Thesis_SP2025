# ============================================================================
#  Appendix C:  Historical Rolling Metrics & Beta Construction (R)
# -----------------------------------------------------------------------------
#  This script pulls daily total‑return data via tidyquant for every stock in
#  `companies.csv` plus the S&P 500 index, then computes rolling volatility,
#  correlation and CAPM beta over 6‑, 12‑, and 24‑month windows.  It also stores
#  higher moments (skewness, kurtosis) and shifts the past estimates forward so
#  that “future” metrics line up with the option‑implied moments extracted in
#  Appendix B‑X.
#
#  Key outputs per stock‑date:
#      •  stock_vol_*_past   / *_future   — annualised σ of stock returns
#      •  market_vol_*_*                    — annualised σ of S&P 500
#      •  corr_*_*                         — Pearson correlation( r_i , r_m )
#      •  beta_*_*                         — β_i = ρ × σ_i / σ_m
#      •  skew_*_*, kurt_*_*               — sample skewness / kurtosis
#
#  The resulting tibble contains ~1.4 million rows (2005‑present) and is saved
#  as `rolling_metrics.csv`, ready for the forecasting exercises in Section 5.
# ============================================================================

# ----------------------------- 1.  Setup ------------------------------------
library(tidyquant)   # price retrieval from Yahoo! Finance
library(dplyr)       # data wrangling
library(lubridate)   # dates
library(zoo)         # rollapply
library(tictoc)      # timing
library(readr)       # CSV I/O
library(tidyverse)   # misc helpers
library(moments)     # skewness / kurtosis

setwd("/Users/ivankolesnikov/Desktop/Spring 2025/Thesis")  # adjust as needed

# --------------------------- 2.  Universe -----------------------------------
companies <- read_csv("companies.csv")
tickers   <- unique(companies$company_name)   # 236 perennial S&P 500 firms

# --------------------------- 3.  Rolling helper -----------------------------
# Returns an 8‑column vector for each rolling window idx ----------------------
compute_rolling_metrics <- function(stock_ret, market_ret, window) {
  rollapply(
    seq_along(stock_ret),                # index vector
    width = window,                      # 126 / 252 / 504 trading days
    FUN = function(idx) {
      s_vol <- sd(stock_ret[idx],  na.rm = TRUE)
      m_vol <- sd(market_ret[idx], na.rm = TRUE)
      rho   <- cor(stock_ret[idx], market_ret[idx], use = "complete.obs")
      beta  <- rho * (s_vol / m_vol)
      c(s_vol, m_vol, rho, beta,
        skewness(stock_ret[idx],  na.rm = TRUE),
        skewness(market_ret[idx], na.rm = TRUE),
        kurtosis(stock_ret[idx],  na.rm = TRUE),
        kurtosis(market_ret[idx], na.rm = TRUE))
    },
    by.column = FALSE, fill = NA, align = "right")
}

# ------------------------ 4.  Master function -------------------------------
get_stock_metrics <- function(tickers,
                              start_date   = "2005-01-01",
                              market_index = "^GSPC") {
  # -- 4.1  Pull prices -------------------------------------------------------
  stock_px  <- tq_get(tickers,      from = start_date) %>%
    select(symbol, date, adjusted) %>% arrange(symbol, date)
  market_px <- tq_get(market_index, from = start_date) %>%
    select(date, adjusted) %>% rename(market_adj = adjusted)
  
  # -- 4.2  Merge & daily log‑returns ----------------------------------------
  stock_rt <- stock_px %>%
    left_join(market_px, by = "date") %>%
    group_by(symbol) %>%
    mutate(stock_ret  = adjusted      / lag(adjusted)      - 1,
           market_ret = market_adj   / lag(market_adj)   - 1) %>%
    drop_na(stock_ret, market_ret) %>% ungroup()
  
  # -- 4.3  Rolling metrics ---------------------------------------------------
  stock_rt <- stock_rt %>% group_by(symbol) %>%
    mutate(
      # past windows ----------------------------------------------------------
      roll_6  = compute_rolling_metrics(stock_ret, market_ret, 126),
      roll_12 = compute_rolling_metrics(stock_ret, market_ret, 252),
      roll_24 = compute_rolling_metrics(stock_ret, market_ret, 504),
      
      # unpack 6‑month --------------------------------------------------------
      stock_vol_6m_past   = roll_6[,1],  market_vol_6m_past   = roll_6[,2],
      corr_6m_past        = roll_6[,3],  beta_6m_past         = roll_6[,4],
      stock_skew_6m_past  = roll_6[,5],  market_skew_6m_past  = roll_6[,6],
      stock_kurt_6m_past  = roll_6[,7],  market_kurt_6m_past  = roll_6[,8],
      
      # unpack 12‑month -------------------------------------------------------
      stock_vol_12m_past  = roll_12[,1], market_vol_12m_past  = roll_12[,2],
      corr_12m_past       = roll_12[,3], beta_12m_past        = roll_12[,4],
      stock_skew_12m_past = roll_12[,5], market_skew_12m_past = roll_12[,6],
      stock_kurt_12m_past = roll_12[,7], market_kurt_12m_past = roll_12[,8],
      
      # unpack 24‑month -------------------------------------------------------
      stock_vol_24m_past  = roll_24[,1], market_vol_24m_past  = roll_24[,2],
      corr_24m_past       = roll_24[,3], beta_24m_past        = roll_24[,4],
      stock_skew_24m_past = roll_24[,5], market_skew_24m_past = roll_24[,6],
      stock_kurt_24m_past = roll_24[,7], market_kurt_24m_past = roll_24[,8],
      
      # future windows (lead by same length) ---------------------------------
      stock_vol_6m_future   = lead(stock_vol_6m_past,   126),
      market_vol_6m_future  = lead(market_vol_6m_past,  126),
      corr_6m_future        = lead(corr_6m_past,        126),
      beta_6m_future        = lead(beta_6m_past,        126),
      stock_skew_6m_future  = lead(stock_skew_6m_past,  126),
      market_skew_6m_future = lead(market_skew_6m_past, 126),
      stock_kurt_6m_future  = lead(stock_kurt_6m_past,  126),
      market_kurt_6m_future = lead(market_kurt_6m_past, 126),
      
      stock_vol_12m_future   = lead(stock_vol_12m_past,   252),
      market_vol_12m_future  = lead(market_vol_12m_past,  252),
      corr_12m_future        = lead(corr_12m_past,        252),
      beta_12m_future        = lead(beta_12m_past,        252),
      stock_skew_12m_future  = lead(stock_skew_12m_past,  252),
      market_skew_12m_future = lead(market_skew_12m_past, 252),
      stock_kurt_12m_future  = lead(stock_kurt_12m_past,  252),
      market_kurt_12m_future = lead(market_kurt_12m_past, 252),
      
      stock_vol_24m_future   = lead(stock_vol_24m_past,   504),
      market_vol_24m_future  = lead(market_vol_24m_past,  504),
      corr_24m_future        = lead(corr_24m_past,        504),
      beta_24m_future        = lead(beta_24m_past,        504),
      stock_skew_24m_future  = lead(stock_skew_24m_past,  504),
      market_skew_24m_future = lead(market_skew_24m_past, 504),
      stock_kurt_24m_future  = lead(stock_kurt_24m_past,  504),
      market_kurt_24m_future = lead(market_kurt_24m_past, 504)
    ) %>%
    select(-roll_6, -roll_12, -roll_24) %>%
    ungroup()
  
  stock_rt %>% rename(ticker = symbol)
}

# ----------------------------- 5.  Execute ----------------------------------
tic("Rolling metrics")
metrics_results <- get_stock_metrics(tickers)
toc()

write_csv(metrics_results, "rolling_metrics.csv")
# ============================================================================
