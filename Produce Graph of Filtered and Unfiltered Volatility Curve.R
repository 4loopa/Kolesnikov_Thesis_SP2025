# ============================================================================
#  Appendix C:  Implied-Volatility Smile Visualization (Filtered vs Unfiltered)
# -----------------------------------------------------------------------------
#  This script produces Figure 1 and Figure 2 in Appendix B.10, contrasting raw
#  and filtered implied-volatility smiles for three representative S&P 500 names
#  (Apple, Broadcom, Eli Lilly) on 4 Jan 2024.
#
#  •  Section 1 loads the ORATS daily file, applies basic no‑arbitrage filters,
#    then plots the *unfiltered* call-IV smiles with a cubic‑spline overlay.
#  •  Section 2 re-applies the stricter parity filter |cIV−pIV|/avg ≤ 20 %, then
#    re‑plots the *filtered* smiles.
#
#  Fonts: Times New Roman via **extrafont** (install once with font_import()).
# ============================================================================

# ----------------------------- 1.  Setup ------------------------------------
library(readr)
library(dplyr)
library(lubridate)
library(splines)    # ns()
library(ggplot2)
library(extrafont)  # Times New Roman embedding

setwd("C:/Users/ikolesnikov25/OneDrive - Claremont McKenna College/Desktop/Thesis")

file_040124 <- "ORATS_SMV_strikes_20240104.csv"   # Jan‑04‑2024 snapshot

# Helper to pick nearest‑below YTE -------------------------------------------
nearest_yte <- function(df, tgt = 0.5) {
  y <- max(df$yte[df$yte <= tgt], na.rm = TRUE)
  filter(df, yte == y)
}

# Target tickers --------------------------------------------------------------
sel <- c("AAPL", "AVGO", "LLY")

# ----------------------- 2.  Unfiltered Smile -------------------------------
raw <- read_csv(file_040124, show_col_types = FALSE) %>%
  select(ticker, stkPx, trade_date, yte, strike, cMidIv, pMidIv, iRate,
         cValue, pValue, cVolu, pVolu, cAskPx, cBidPx, pAskPx, pBidPx) %>%
  filter((stkPx - strike) <= cValue,       # intrinsic ≤ option price
         (strike - stkPx) <= pValue) %>%
  mutate(trade_date       = as.Date(trade_date, "%m/%d/%Y"),
         moneyness        = strike / stkPx,
         call_spread_perc = (cAskPx - cBidPx) / cAskPx,
         put_spread_perc  = (pAskPx - pBidPx) / pAskPx) %>%
  filter(ticker %in% sel) %>%
  group_by(ticker) %>% group_modify(~ nearest_yte(.x)) %>% ungroup()

# Facet labels ----------------------------------------------------------------
raw$Company <- recode(raw$ticker, AAPL = "Apple", AVGO = "Broadcom", LLY = "Eli Lilly")

p1 <- ggplot(raw, aes(moneyness, cMidIv)) +
  geom_point(alpha = 0.7, size = 2, colour = "#2c3e50") +
  geom_smooth(method = "lm", formula = y ~ ns(x, df = 3), se = FALSE,
              colour = "#e74c3c", size = 1.5) +
  facet_wrap(~ Company, ncol = 3) +
  labs(title = "Implied Volatility vs. Moneyness (Unfiltered, 0.5 Y)",
       x = "Moneyness (K/S)", y = "Call Mid IV") +
  theme_minimal(base_size = 14, base_family = "Times New Roman") +
  theme(plot.title = element_text(size = 18, hjust = 0.5, colour = "#34495e"),
        strip.background = element_rect(fill = "#ecf0f1", colour = "#bdc3c7"))

# ----------------------- 3.  Filtered Smile ---------------------------------
flt <- read_csv(file_040124, show_col_types = FALSE) %>%
  select(ticker, stkPx, trade_date, yte, strike, cMidIv, pMidIv, iRate,
         cValue, pValue, cVolu, pVolu, cAskPx, cBidPx, pAskPx, pBidPx) %>%
  filter((stkPx - strike) <= cValue,
         (strike - stkPx) <= pValue) %>%
  mutate(trade_date = as.Date(trade_date, "%m/%d/%Y"),
         moneyness  = strike / stkPx) %>%
  filter(abs(cMidIv - pMidIv) / ((cMidIv + pMidIv) / 2) <= 0.2,  # IV parity
         ticker %in% sel) %>%
  group_by(ticker) %>% group_modify(~ nearest_yte(.x)) %>% ungroup()

flt$Company <- recode(flt$ticker, AAPL = "Apple", AVGO = "Broadcom", LLY = "Eli Lilly")

p2 <- ggplot(flt, aes(moneyness, cMidIv)) +
  geom_point(alpha = 0.7, size = 2, colour = "#2c3e50") +
  geom_smooth(method = "lm", formula = y ~ ns(x, df = 3), se = FALSE,
              colour = "#e74c3c", size = 1.5) +
  facet_wrap(~ Company, ncol = 3) +
  labs(title = "Implied Volatility vs. Moneyness (Filtered, 0.5 Y)",
       x = "Moneyness (K/S)", y = "Call Mid IV") +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(size = 18, hjust = 0.5, colour = "#34495e"),
        strip.background = element_rect(fill = "#ecf0f1", colour = "#bdc3c7"))

# ----------------------- 4.  Display or Save --------------------------------
# gridExtra::grid.arrange(p1, p2, nrow = 2)  # shows both in RStudio
# ggsave("iv_smile_unfiltered.png", p1, width = 9, height = 3)
# ggsave("iv_smile_filtered.png",   p2, width = 9, height = 3)
# ===========================================================================
