# Description: EU-SILC Indicators
# Content: This file creates EU-SILC indicators on a sample data. It is expected to be adapted to the needs of each country
#			For analysis over 1 year, adapt the codes from line 13
#			For analysis over multiple years (panel), see example from line 99 with a panel data sample


# Setup --------------------------------------------------------------------

library(haven)
library(dplyr)
library(purrr)
library(tidyr)
library(stringr)
library(readr)
library(Hmisc)
library(openxlsx)
library(reldist)
library(survey)

# Define paths to data ---------------------------------------------------
data_raw <- "./data/raw/silk-sample" # path to the SLIC data. Uses simulated data for now.
data_aux <- "./data/temp"            # Path to where aggregated data is saved

# General guidance ----------------------------------------------------------------

# 1. Ensure that sub national unit variables are properly identifies (if applies)
# 2. Combine multiple years of data into one indicating the year variable

# Typical data structure --------------------------------------------

# db040  - sub national id. must uniquely identify your geographic level
# hhid   - household identifier (Possibly DB030 variable)
# hx090  - Equivalised disposable income
# hx040  - Household size
# hx050  - Equivalised household size
# rb050  - personal cross-sectional weigh
# rx070  - At risk of poverty or social exclusion
# psu    - First stage randomization (possibly DB060)
# strata - Stratum of the second stage randomization (possibly DB062)
# code   - Country code

# Loading data ----------------------------------------------------------------

df <- 
  file.path(data_raw, "dataexample2021.dta") |> read_dta() |> 
  bind_rows(file.path(data_raw, "dataexample2022.dta") |> read_dta()) |> 
  bind_rows(file.path(data_raw, "dataexample2023.dta") |> read_dta()) |> 
  
  # Clean any Stata-related metadata
  zap_labels() |> 
  
  # Rename native variable into the relavant names here, specifically:
  # hhid, subcode, psu, strata, and others if needed. Use syntax:
  # rename(NEW_WAR_NAME = OLD_VAR_NAME). All names are case-sensitive
  rename(subcode = db040)

# Aggregating a household-level data
df_hh <-
  df |> 
  group_by(code, year, hhid) |> 
  summarise(
    across(c(hx090, hx040, hx050, rb050, rx070, psu, strata, subcode),
           ~ first(.)), 
    .groups = "drop"
  ) |> 
    
  # We must group by country/year in order to compute all statistics on the annual basis.
  group_by(code, year) |> 
  mutate(
    popw = rb050 * hx040,                           # population weights 
    inc_pc_euro = (hx090 * hx050) / hx040,          # HH disposable income per capita
    inc_adeq_euro = hx090,                          # HH disposable income per adult equivalent
    
    arop = 0.6 * Hmisc::wtd.quantile(
      x = inc_adeq_euro, 
      weights = popw, 
      probs = 0.5, 
      na.rm = TRUE),
    
    arope = ifelse(is.na(rx070), NA, as.numeric(rx070 != 0)),
    
    fgt0_arop = as.integer(inc_adeq_euro < arop),
    fgt1_arop = ifelse(fgt0_arop, 1 - (inc_adeq_euro / arop), 0),
    fgt2_arop = ifelse(fgt0_arop, (1 - (inc_adeq_euro / arop)) ^ 2, 0)
  ) |> 
  select(-rb050, -hx040, -hx090, -hx050) |> 
  group_by(code, year, subcode) |> 
  mutate(
    b40_inc_adeq_euro =
      ifelse(
        inc_adeq_euro <= Hmisc::wtd.quantile(x = inc_adeq_euro, weights = popw,
                                             probs = 0.4, na.rm = TRUE),
        inc_adeq_euro,
        NA_real_
      ),

    b20_inc_adeq_euro =
      ifelse(
        inc_adeq_euro <= Hmisc::wtd.quantile(x = inc_adeq_euro, weights = popw,
                                             probs = 0.2, na.rm = TRUE),
        inc_adeq_euro,
        NA_real_
      ),

  ) |> 
  ungroup() 

# Survey-base aggregates that include design -----------------

survey_design <- svydesign(
  id = ~ psu,
  weights = ~ popw,
  strata = ~ strata,
  nest = TRUE,
  data = df_hh
)

# Variables list for which we compute means and variances ------------------
# Note: Make sure check and report the number of missing observations in these variables
stats_cols <- c("inc_adeq_euro", 
                "fgt0_arop", "fgt1_arop", "fgt2_arop", "arope", 
                "b40_inc_adeq_euro", "b20_inc_adeq_euro")

# Regular population weights-weighted statistics
df_wtd <- 
  df_hh |> 
  group_by(code, year, subcode) |> 
  summarise(
    n_obs = n(), # Number of observations per region
    sum_popw = sum(popw),  # Sum of the population weights
    gini_inc_adeq_euro = reldist::gini(x = inc_adeq_euro, weights = popw),
    across(any_of(stats_cols), 
           ~ Hmisc::wtd.mean(x = ., weights = popw, na.rm = TRUE), 
           .names = "mean_{.col}"),
    across(c(inc_adeq_euro), 
           ~ Hmisc::wtd.quantile(x = ., weights = popw, probs = 0.5, na.rm = TRUE), 
           .names = "median_{.col}"), 
    .groups = "drop" 
  )

# Survey-design corrected variances
df_svywtd <-
  stats_cols |>
  map( ~ {
    var_focus <- .x
    svyby(
      formula = paste("~", var_focus) |> as.formula(),
      by = ~ code + year + subcode,
      design = survey_design,
      FUN = svymean,
      deff = FALSE,
      na.rm = TRUE,
      vartype = "var"
    ) |>
      as_tibble() |>
      select(-any_of(var_focus)) |>
      rename_with( ~ str_c("var_", var_focus), .cols = "var")
  }) |>
  reduce(left_join, by = join_by(code, year, subcode)) |>
  arrange(year, subcode) 

# Compiling all data ------------------------------------------------------
df_out <- df_wtd |> left_join(df_svywtd, by = join_by(code, year, subcode))

# Saving data -----------------------------------------------------------
write_rds(df_out, file.path(data_aux, "eu-silc.rds"), compress = "gz")
write_csv(df_out, file.path(data_aux, "eu-silc.csv"))

# Save as Excel (.xlsx), if needed
write.xlsx(df_out,
           file.path(data_aux, "EU_SILC_Indicators.xlsx"),
           sheetName = "R",
           overwrite = F)
