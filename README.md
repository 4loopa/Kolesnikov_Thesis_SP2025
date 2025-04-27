# Forecasting Betas with Option-Implied Moments
Companion repository for my senior thesis **“Forecasting Equity Betas Using Option-Implied Moments.”**  
The R scripts reproduce every result: extracting risk-neutral variance/skew from ORATS Near-EOD option quotes, computing rolling historical betas, building a master panel, and benchmarking forward-looking estimators (OLS & random forest).

---

## Quick start

```bash
# clone & enter
git clone https://github.com/<your-handle>/forecasting-betas-option-moments.git
cd forecasting-betas-option-moments

# install R packages (renv or plain install.packages)
Rscript -e "install.packages(c(
  'tidyverse','tidymodels','tidyquant','splines',
  'randomForest','moments','extrafont','gridExtra'))"

# populate data/
#   ORATS/                    # yearly folders of daily CSVs (≈500 GB unzipped)
#   companies.csv             # ticker → GICS sector lookup
#   ORATS_SMV_strikes_20240104.csv   # for IV-smile figures

# run the full pipeline step-by-step
Rscript 'C-1 Option-Implied Moments Calculation.R'     # 15-20 h
Rscript 'C-2 Calculate Historical Moments.R'           # 1 h
Rscript 'C-3 Master Data File Generation.R'            # seconds
Rscript 'C-4 Method Accuracy and Regressions.R'        # 10 min
Rscript 'C-5 Produce Graph of Filtered and Unfiltered IV.R'
Rscript 'C-6 RF RMSE and Importance Check.R'           # 5 min
Rscript 'C-7 Market Returns Skewness Calculations.R'
Rscript 'C-8 Residual Skewness.R'
