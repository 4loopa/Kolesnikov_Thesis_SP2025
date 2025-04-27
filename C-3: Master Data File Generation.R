# ============================================================================
#  Appendix C:  Master Dataset Construction – Option‑Implied + Historical
# -----------------------------------------------------------------------------
#  This script merges the option‑implied moment files (final_results_YYYY.csv)
#  with the rolling historical metrics (rolling_metrics.csv) to create the
#  canonical `master_data.csv` used throughout Sections 4–6 of the thesis.
#
#  Outputs
#  --------
#  •  `master_data.csv` – panel keyed by (ticker, date) containing:
#        – risk‑neutral quad / cubic / vol / skew at 0.5, 1, 2 years
#        – interpolated ATM IV and diagnostics (for stock and S&P 500 index)
#        – rolling historical vol, correlation, beta, skew, kurtosis (6/12/24 m)
#  •  The wide format suffixes `_0.5`, `_1`, `_2` follow the option horizon;
#    market fields carry `_m` to denote index counterparts.
# ============================================================================

# ----------------------------- 1.  Setup ------------------------------------
library(readr)       # fast I/O
library(dplyr)       # data wrangling
library(lubridate)   # dates
library(tidyr)       # pivot_wider()
library(purrr)       # map_dfr()

setwd("C:/Users/ikolesnikov25/OneDrive - Claremont McKenna College/Desktop/Thesis")

# --------------------------- 2.  Load Sources -------------------------------
# — historical rolling metrics ----------------------------------------------
betas <- read_csv("rolling_metrics.csv")

# — option‑implied results, 2007‑2024 ----------------------------------------
years <- 2007:2024
opt_files <- paste0("final_results_", years, ".csv")

combined_df <- map_dfr(opt_files, read_csv, .id = "source") %>%
  select(-source, -file_source)          # drop helpers from earlier pipeline

# ---------------------- 3.  Make Option Table Wide --------------------------
combined_wide <- combined_df %>%
  pivot_wider(
    id_cols   = c(ticker, date),
    names_from  = yte_target,                     # 0.5 / 1 / 2 years
    values_from = -c(ticker, date, yte_target),   # every other column
    names_sep   = "_")

# ---------------------- 4.  Add S&P 500 Market Columns ----------------------
market  <- combined_wide %>% filter(ticker == "SPX")
merged  <- combined_wide %>% left_join(market, by = "date", suffix = c("", "_m"))

# ---------------------- 5.  Merge Historical Metrics ------------------------
master <- merged %>% left_join(betas, by = c("date", "ticker"))

write_csv(master, "master_data.csv")
# ===========================================================================
