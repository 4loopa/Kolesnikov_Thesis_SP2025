# ============================================================================
#  Appendix C:  S&P 500 Skewness Diagnostics (R)
# -----------------------------------------------------------------------------
#  Purpose
#  -------
#  This short script documents the normality checks reported in Section 6.2 of
#  the thesis.  Using the longest‑available S&P 500 (^GSPC) time series from
#  Yahoo Finance, it constructs:
#      •  180‑day independent returns sampled every 6 months (Jan‑1 & Jul‑1)
#        — ‘Full sample’ : 1929‑07‑01 .. today
#        — ‘Modern era’  : 2005‑07‑01 .. today
#      •  Daily log‑returns over the same horizon for comparison.
#
#  Each return series is subjected to the D’Agostino skewness test
#  (package **moments**).  P‑values reported in the thesis originate from this
#  exact code.
#
#  Notes
#  -----
#  •  Yahoo’s ^GSPC history starts 1928‑10‑01, so setting `from = "1900‑01‑01"`
#     simply requests the full available range.
#  •  A 180‑day log‑return is approximated by `adjusted / lag(adjusted, 124) ‑ 1`.
#     Because the index is traded almost every calendar day, 124 trading days ≈
#     180 calendar days.
#  •  The semi‑annual subsample (Jan 1 / Jul 1) yields ~41 returns between
#     2005‑07‑01 and 2025‑01‑01, matching the count quoted in the thesis.
# ============================================================================

# ----------------------------- 1.  Setup ------------------------------------
library(tidyquant)   # Yahoo Finance interface
library(readr)       # CSV I/O (not used but loaded for completeness)
library(tidyverse)   # dplyr + ggplot2 helpers
library(moments)     # agostino.test()

setwd("/Users/ivankolesnikov/Desktop/Spring 2025/Thesis")  # adjust as needed

# --------------------------- 2.  Fetch Prices -------------------------------
sp500_prices <- tq_get(
  "^GSPC",                # S&P 500 composite index
  from = "1900-01-01",    # earliest possible; Yahoo truncates to 1928‑10‑01
  to   = Sys.Date())

# --------------------- 3.  180‑day Independent Returns ----------------------
# Approximate 180‑day total return (not log‑return) ---------------------------
mkt <- sp500_prices %>%
  select(date, adjusted) %>%
  mutate(ret = adjusted / lag(adjusted, 124) - 1) %>%   # 124 ≈ 180 trading days
  drop_na()

# --- 3.1  Full-sample semi‑annual observation grid --------------------------
mkt_iid <- mkt %>%
  filter(month(date) %in% c(1, 7),       # Jan or Jul
         day(date)   %in% c(1, 7))       # 1st or 7th (cushion for holidays)

cat("D'Agostino test – full sample (1929‑present):\n")
print(agostino.test(mkt_iid$ret))

# --- 3.2  Modern‑era subsample (2005‑07‑01 .. today) ------------------------
mkt_iid_modern <- mkt_iid %>%
  filter(date >= as.Date("2005-07-01"))

cat("\nD'Agostino test – modern era (since 2005‑07‑01):\n")
print(agostino.test(mkt_iid_modern$ret))

# --------------------- 4.  Daily Returns Normality Check --------------------
mkt_daily <- sp500_prices %>%
  select(date, adjusted) %>%
  mutate(ret = adjusted / lag(adjusted) - 1) %>%
  drop_na()

cat("\nD'Agostino test – daily returns (1928‑present):\n")
print(agostino.test(mkt_daily$ret))

# ---------------------------------------------------------------------------
# End of script – outputs printed to console; capture with sink() if needed.
# ===========================================================================
