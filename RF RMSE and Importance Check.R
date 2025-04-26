# ============================================================================
#  Appendix C‑3:  Random‑Forest Beta Forecasting Workflow (R)
# -----------------------------------------------------------------------------
#  This script builds and evaluates random‑forest models for 6‑, 12‑, and
#  24‑month realized betas.  Commands have been reordered for clarity; no new
#  functionality has been added.
# ============================================================================

# ----------------------------- 1.  Setup ------------------------------------
# Core data wrangling / modelling libs
library(readr)        # CSV I/O
library(dplyr)        # data manipulation
library(tidyverse)    # ggplot2, purrr, etc.
library(tidymodels)   # recipes, rsample, yardstick, workflows
library(randomForest) # RF engine
library(tictoc)       # timing

setwd("/Users/ivankolesnikov/Desktop/Spring 2025/Thesis")  # adjust as needed

# --------------------------- 2.  Load & Join --------------------------------
master <- read_csv("master_data.csv") %>%
  select(-...1)                             # drop stray index col

companies <- read_csv("companies.csv") %>%
  rename(ticker = company_name,
         sector = `GICS Sector`) %>%
  select(ticker, sector)

master <- master %>% left_join(companies, by = "ticker")

# ----------------------- 3.  Feature Engineering ----------------------------
master <- master %>%
  select(-ticker) %>%
  mutate(
    # option‑implied (OI) betas ----------------------------------------------
    beta_oi_6m  = (skew_0.5_m / skew_0.5)^(1/3) * vol_0.5 / vol_0.5_m,
    beta_oi_12m = (skew_1_m   / skew_1)  ^(1/3) * vol_1   / vol_1_m,
    beta_oi_24m = (skew_2_m   / skew_2)  ^(1/3) * vol_2   / vol_2_m,
    # combined (C) betas ------------------------------------------------------
    beta_c_6m  = corr_6m_past  * vol_0.5 / vol_0.5_m,
    beta_c_12m = corr_12m_past * vol_1   / vol_1_m,
    beta_c_24m = corr_24m_past * vol_2   / vol_2_m) %>%
  drop_na()

# ----------------------- 4.  6‑Month Random Forest --------------------------
# — data split ----------------------------------------------------------------
six_m <- master %>%
  rename(beta_realized = beta_6m_future) %>%
  select(-ends_with("future"))

set.seed(47)
split_6  <- initial_split(six_m)
train_6  <- training(split_6)
test_6   <- testing(split_6)

# — preprocessing recipe ------------------------------------------------------
rec_6 <- recipe(beta_realized ~ ., data = train_6) %>%
  step_mutate(sector = factor(sector))

# — model spec ----------------------------------------------------------------
rf_spec <- rand_forest(mtry = 7, trees = 100) %>%
  set_engine("randomForest", oob.error = TRUE) %>%
  set_mode("regression")

wf_6 <- workflow() %>% add_model(rf_spec) %>% add_recipe(rec_6)

tic("6m RF fit")
fit_6 <- wf_6 %>% fit(data = train_6)
toc()

# — predictions & RMSE --------------------------------------------------------
train_pred_6 <- predict(fit_6, train_6) %>% bind_cols(train_6 %>% select(beta_realized))
test_pred_6 <- predict(fit_6,  test_6) %>% bind_cols( test_6 %>% select(beta_realized))

rmse_train_6 <- rmse(train_pred_6, truth = beta_realized, estimate = .pred)
rmse_test_6  <- rmse( test_pred_6, truth = beta_realized, estimate = .pred)

print(rmse_train_6)
print(rmse_test_6)

# — variable importance -------------------------------------------------------
imp_6 <- importance(extract_fit_engine(fit_6)) %>%
  as_tibble(rownames = "variable") %>%
  arrange(desc(IncNodePurity))
write_csv(imp_6, "Importance_scores_6m.csv")

# ----------------------- 5.  Helper for other horizons -----------------------
run_rf_beta <- function(data, beta_future_col, seed = 47, mtry = 7, trees = 100) {
  set.seed(seed)
  df <- data %>% rename(beta_realized = !!sym(beta_future_col)) %>%
    select(-ends_with("future"))
  split  <- initial_split(df)
  train  <- training(split)
  test   <- testing(split)
  
  rec <- recipe(beta_realized ~ ., data = train) %>% step_mutate(sector = factor(sector))
  mod <- rand_forest(mtry = mtry, trees = trees) %>%
    set_engine("randomForest", oob.error = TRUE) %>%
    set_mode("regression")
  fit <- workflow() %>% add_model(mod) %>% add_recipe(rec) %>% fit(train)
  
  train_rmse <- rmse(bind_cols(predict(fit, train), train %>% select(beta_realized)),
                     truth = beta_realized, estimate = .pred)
  test_rmse  <- rmse(bind_cols(predict(fit, test),  test  %>% select(beta_realized)),
                     truth = beta_realized, estimate = .pred)
  list(rmse_train = train_rmse, rmse_test = test_rmse,
       importance = importance(extract_fit_engine(fit)) %>%
         as_tibble(rownames = "variable") %>% arrange(desc(IncNodePurity)))
}

# ----------------------- 6.  12m & 24m Full Models --------------------------
res_12m <- run_rf_beta(master, "beta_12m_future")
res_24m <- run_rf_beta(master, "beta_24m_future")

# ----------------------- 7.  Limited Feature Sets ---------------------------
master_lim_oi <- master %>% select(sector,
                                   beta_6m_future, beta_12m_future, beta_24m_future,
                                   beta_6m_past, beta_12m_past, beta_24m_past,
                                   beta_oi_6m, beta_oi_12m, beta_oi_24m,
                                   beta_c_6m,  beta_c_12m,  beta_c_24m)

master_lim_hist <- master %>% select(sector,
                                     beta_6m_future, beta_12m_future, beta_24m_future,
                                     beta_6m_past,   beta_12m_past,   beta_24m_past)

lim_6m_oi   <- run_rf_beta(master_lim_oi,   "beta_6m_future")
lim_6m_hist <- run_rf_beta(master_lim_hist, "beta_6m_future")
lim_12m_oi  <- run_rf_beta(master_lim_oi,   "beta_12m_future")
lim_12m_hist<- run_rf_beta(master_lim_hist, "beta_12m_future")
lim_24m_oi  <- run_rf_beta(master_lim_oi,   "beta_24m_future")
lim_24m_hist<- run_rf_beta(master_lim_hist, "beta_24m_future")

# Print headline RMSEs --------------------------------------------------------
print(list(full_12m = res_12m[1:2], full_24m = res_24m[1:2],
           lim6_oi  = lim_6m_oi[1:2],  lim6_hist  = lim_6m_hist[1:2],
           lim12_oi = lim_12m_oi[1:2], lim12_hist = lim_12m_hist[1:2],
           lim24_oi = lim_24m_oi[1:2], lim24_hist = lim_24m_hist[1:2]))
# ===========================================================================
