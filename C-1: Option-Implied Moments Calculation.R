# ============================================================================
#  Appendix B‑X:  Option‑Implied Moments Extraction Script (R)
# -----------------------------------------------------------------------------
#  This script reproduces the entire preprocessing pipeline used in the thesis
#  “Forecasting Equity Betas Using Option‑Implied Moments.”  Executing it on the
#  Orats Near‑EOD dataset will
#     1.  Clean raw option quotes (2007‑2024) for the 236 perennial S&P 500 firms;
#     2.  Fit natural‑cubic‑spline volatility smiles at each maturity;
#     3.  Re‑price a dense strike grid via Black–Scholes to obtain synthetic
#         option prices;
#     4.  Compute option‑implied variance, skewness, and the risk‑neutral mean
#         return using the Bakshi–Kapadia–Madan (2003) integrals;
#     5.  Interpolate these moments to the target horizons (0.5, 1, and 2 years);
#     6.  Export a per‑year CSV of all moments plus diagnostics.
#
#  The code is heavily commented to serve as a self‑contained reference for
#  replication.  Paths assume a Windows workstation; adjust `setwd()` and
#  `folder_path` as needed.  Run time on the full dataset (~500 GB) is roughly
#  15–20 hours on a modern 8‑core CPU with 32 GB RAM.
# ============================================================================

# ----------------------------- 1.  Setup ------------------------------------

# Load required libraries ------------------------------------------------------
library(readr)      # fast CSV I/O
library(readxl)     # Excel import (not used here but loaded for completeness)
library(dplyr)      # data manipulation verbs
library(lubridate)  # date handling
library(splines)    # natural cubic splines via bs()
library(tidyverse)  # ggplot2, purrr, etc.
library(ggplot2)

# Working directory -----------------------------------------------------------
# NOTE: Change to the root of your Orats directory before running.
setwd("C:/Users/ikolesnikov25/Desktop/Thesis")

# --------------------------- 2.  Global Params ------------------------------

years        <- 2017:2024           # Years to process
yte_targets  <- c(0.5, 1, 2)        # Target maturities (in years)

# Universe of tickers ---------------------------------------------------------
companies <- read_csv("companies.csv")
tickers   <- as.list(companies$company_name)  # 236 perennial S&P 500 members

# Single moneyness grid (1 % – 300 %) ----------------------------------------
moneyness_grid <- seq(0.01, 3, length.out = 1000)

# --------------------------- 3.  Helper Functions ---------------------------

# Black–Scholes price ----------------------------------------------------------
# `type`   : "call" | "put"
# `S`      : underlying price
# `K`      : strike price
# `r`      : continuously‑compounded risk‑free rate
# `tao`    : time‑to‑expiry in years
# `sigma`  : implied volatility (σ)
european_option_price <- function(type, S, K, r, tao, sigma) {
  d1 <- (log(S / K) + (r + 0.5 * sigma^2) * tao) / (sigma * sqrt(tao))
  d2 <- d1 - sigma * sqrt(tao)
  if (type == "call") {
    S * pnorm(d1) - K * exp(-r * tao) * pnorm(d2)
  } else if (type == "put") {
    K * exp(-r * tao) * pnorm(-d2) - S * pnorm(-d1)
  } else {
    stop("Invalid option type. Use 'call' or 'put'.")
  }
}

# Fit natural‑cubic‑spline IV smile -------------------------------------------
fit_volatility_spline <- function(moneyness, vols, df = 3) {
  # Returns an lm() object mapping moneyness → implied vol.
  lm(vols ~ splines::bs(moneyness, df = df),
     data = data.frame(moneyness = moneyness, vols = vols))
}

# Predict IV on grid with flat extrapolation beyond observed strikes ----------
predict_with_extrapolation <- function(spline_model, m_grid,
                                       orig_moneyness, orig_vols) {
  lower <- min(orig_moneyness)
  upper <- max(orig_moneyness)
  iv_hat <- predict(spline_model, newdata = data.frame(moneyness = m_grid))
  # Flat extrapolation outside fitted region
  iv_hat[m_grid < lower] <- orig_vols[which.min(orig_moneyness)]
  iv_hat[m_grid > upper] <- orig_vols[which.max(orig_moneyness)]
  iv_hat
}

# Construct full IV curve by stitching put‑ and call‑based smiles -------------
type_based_prediction <- function(call_spline, put_spline, call_data, put_data) {
  data.frame(
    moneyness = moneyness_grid,
    iv = ifelse(moneyness_grid <= 1,
                predict_with_extrapolation(put_spline,  moneyness_grid,
                                           put_data$moneyness,  put_data$pMidIv),
                predict_with_extrapolation(call_spline, moneyness_grid,
                                           call_data$moneyness, call_data$cMidIv))
  )
}

# Linear interpolation of option‑implied moments -----------------------------
interpolate_moments <- function(lower, higher, lower_yte, higher_yte, target_yte) {
  w <- (target_yte - lower_yte) / (higher_yte - lower_yte)
  list(
    quad  = lower$quad  + w * (higher$quad  - lower$quad ),
    cubic = lower$cubic + w * (higher$cubic - lower$cubic),
    vol   = lower$vol   + w * (higher$vol   - lower$vol  ),
    skew  = lower$skew  + w * (higher$skew  - lower$skew ),
    Er_R  = lower$Er_R  + w * (higher$Er_R  - lower$Er_R )
  )
}

# Compute option‑implied moments (Bakshi–Kapadia–Madan integrals) ------------
compute_moments <- function(grid, spot, r, yte) {
  # Split grid into call and put regions (K≥S and K≤S) ------------------------
  call_grid <- grid %>% filter(K >= spot)
  put_grid  <- grid %>% filter(K <= spot)
  
  # --- 3.1 Quadratic (variance swap) payoff ----------------------------------
  call_q <- sum(((2 * (1 - log(call_grid$K / spot))) / (call_grid$K^2)) *
                  call_grid$cPrice) * mean(diff(call_grid$K))
  put_q  <- sum(((2 * (1 + log(spot / put_grid$K))) / (put_grid$K^2)) *
                  put_grid$pPrice) * mean(diff(put_grid$K))
  quad <- call_q + put_q
  
  # --- 3.2 Cubic payoff (skewness contract) ----------------------------------
  cubic_call <- sum(((6 * log(call_grid$K / spot) - 3 * (log(call_grid$K / spot))^2) /
                       (call_grid$K^2)) * call_grid$cPrice) * mean(diff(call_grid$K))
  cubic_put  <- sum(((6 * log(spot / put_grid$K) + 3 * (log(spot / put_grid$K))^2) /
                       (put_grid$K^2)) * put_grid$pPrice) * mean(diff(put_grid$K))
  cubic <- cubic_call - cubic_put
  
  # --- 3.3 Risk‑neutral moments ---------------------------------------------
  Er_R <- exp(r * yte) - 1 - 0.5 * exp(r * yte) * quad - (1/6) * exp(r * yte) * cubic
  var  <- exp(r * yte) * quad - Er_R^2
  skew <- (exp(r * yte) * cubic - 3 * Er_R * exp(r * yte) * quad + 2 * Er_R^3) /
    var^(3/2)
  
  list(quad = quad, cubic = cubic, vol = sqrt(var), skew = skew, Er_R = Er_R)
}

# Master wrapper for a single ticker‑date‑YTE combination ---------------------
calculation <- function(data, ticker, yte_target) {
  stock <- data |> filter(ticker == !!ticker)
  # Identify nearest maturities around target ---------------------------------
  lower_yte  <- max(stock$yte[stock$yte <= yte_target], na.rm = TRUE)
  higher_yte <- min(stock$yte[stock$yte >= yte_target], na.rm = TRUE)
  if (is.infinite(lower_yte) || is.infinite(higher_yte)) {
    skipped_tickers <<- rbind(skipped_tickers,
                              data.frame(ticker = ticker,
                                         date   = if (nrow(stock) > 0)
                                           stock$trade_date[1] else NA))
    return(NULL)
  }
  
  # Split data by maturity ----------------------------------------------------
  lower_stock  <- filter(stock, yte == lower_yte)
  higher_stock <- filter(stock, yte == higher_yte)
  
  # Fit IV smiles -------------------------------------------------------------
  lower_call_spline  <- fit_volatility_spline(lower_stock$moneyness,
                                              lower_stock$cMidIv)
  lower_put_spline   <- fit_volatility_spline(lower_stock$moneyness,
                                              lower_stock$pMidIv)
  higher_call_spline <- fit_volatility_spline(higher_stock$moneyness,
                                              higher_stock$cMidIv)
  higher_put_spline  <- fit_volatility_spline(higher_stock$moneyness,
                                              higher_stock$pMidIv)
  
  # Spot prices and rates -----------------------------------------------------
  price_lower  <- lower_stock$stkPx[1]
  r_lower      <- lower_stock$iRate[1]
  price_higher <- higher_stock$stkPx[1]
  r_higher     <- higher_stock$iRate[1]
  
  # Build full grids of IV + synthetic prices ---------------------------------
  lower_grid <- type_based_prediction(lower_call_spline, lower_put_spline,
                                      lower_stock, lower_stock) |>
    mutate(K      = moneyness * price_lower,
           cPrice = european_option_price("call", price_lower, K, r_lower,
                                          lower_yte, iv),
           pPrice = european_option_price("put",  price_lower, K, r_lower,
                                          lower_yte, iv))
  
  higher_grid <- type_based_prediction(higher_call_spline, higher_put_spline,
                                       higher_stock, higher_stock) |>
    mutate(K      = moneyness * price_higher,
           cPrice = european_option_price("call", price_higher, K, r_higher,
                                          higher_yte, iv),
           pPrice = european_option_price("put",  price_higher, K, r_higher,
                                          higher_yte, iv))
  
  # Compute lower / higher moments -------------------------------------------
  lower_moments  <- compute_moments(lower_grid,  price_lower,  r_lower,  lower_yte)
  higher_moments <- compute_moments(higher_grid, price_higher, r_higher, higher_yte)
  
  # Interpolate to target YTE --------------------------------------------------
  target <- interpolate_moments(lower_moments, higher_moments,
                                lower_yte, higher_yte, yte_target)
  
  tibble(
    ticker      = ticker,
    date        = stock$trade_date[1],
    stkPx       = stock$stkPx[1],
    quad        = target$quad,
    cubic       = target$cubic,
    vol         = target$vol,
    skew        = target$skew,
    eqr         = target$Er_R,
    r_lower     = r_lower,
    r_higher    = r_higher,
    lower_yte   = lower_yte,
    higher_yte  = higher_yte,
    # Diagnostics -------------------------------------------------------------
    num_obs_higher = nrow(higher_stock),
    num_obs_lower  = nrow(lower_stock),
    lower_Iv   = mean(predict(lower_call_spline,  data.frame(moneyness = 1)),
                      predict(lower_put_spline,   data.frame(moneyness = 1))),
    higher_Iv  = mean(predict(higher_call_spline, data.frame(moneyness = 1)),
                      predict(higher_put_spline,  data.frame(moneyness = 1))),
    yte_target = yte_target
  )
}

# --------------------------- 4.  Main Loop ----------------------------------

skipped_tickers <- data.frame(ticker = character(), date = as.Date(character()))

for (year in years) {
  all_results <- list()
  folder_path <- paste0("C:/Users/ikolesnikov25/Desktop/Thesis/", year)
  csv_files   <- list.files(folder_path, pattern = "\\.csv$", full.names = TRUE)
  
  for (file in csv_files) {
    message("Processing: ", file)
    
    # --- 4.1  Load and pre‑filter -------------------------------------------
    data <- read_csv(file, show_col_types = FALSE) |>
      select(ticker, stkPx, trade_date, yte, strike, cMidIv, pMidIv, iRate,
             cValue, pValue, cVolu, pVolu, cAskPx, cBidPx, pAskPx, pBidPx) |>
      filter(stkPx - strike <= cValue,            # no‑arbitrage: intrinsic ≤ price
             strike - stkPx <= pValue) |>
      mutate(trade_date      = as.Date(trade_date, format = "%m/%d/%Y"),
             moneyness       = strike / stkPx,
             call_spread_perc = (cAskPx - cBidPx) / cAskPx,
             put_spread_perc  = (pAskPx - pBidPx) / pAskPx) |>
      filter(abs(cMidIv - pMidIv) / ((cMidIv + pMidIv) / 2) <= 0.2)  # IV parity
    
    # --- 4.2  Compute moments at each target YTE ----------------------------
    for (yte in yte_targets) {
      res <- map_dfr(companies$company_name,           # for each ticker
                     ~ calculation(data, .x, yte))
      all_results[[paste(file, yte, sep = "_")]] <- res
    }
  }
  
  # --- 4.3  Save annual results --------------------------------------------
  final_yr <- bind_rows(all_results, .id = "file_source")
  assign(paste0("final_results_", year), final_yr)
  write_csv(final_yr, paste0("final_results_", year, ".csv"))
}

# ---------------------------------------------------------------------------
# End of script
# ===========================================================================
