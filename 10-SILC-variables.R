# Description: EU-SILC Indicators
# Content: This file creates EU-SILC indicators on a sample data. It is expected to be adapted to the needs of each country
#			For analysis over 1 year, adapt the codes from line 13
#			For analysis over multiple years (panel), see example from line 99 with a panel data sample


# Setup --------------------------------------------------------------------

library(haven)
library(dplyr)
library(purrr)
library(tidyr)
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
    across(c(hx090, hx040, hx050, rb050, rx070, psu, strata, subcode), ~ first(.))
  ) |> 
  ungroup() |> 
    
  # We must group by country/year in order to compute all statistics on the annual basis.
  group_by(code, year) |> 
  mutate(
    popw = rb050 * hx040,                           # population weights 
    
    inc_pc_euro = (hx090 * hx050) / hx040,          # HH disposable income per capita
    
    inc_adeq_euro = hx090,                          # HH disposable income per adult equivalent
    
    inc_adeq_euro_40 = 
      ifelse(inc_adeq_euro <= Hmisc::wtd.quantile(inc_adeq_euro, popw, .4, na.rm = TRUE), 
             inc_adeq_euro, 
             NA_real_),
    
    inc_adeq_euro_20 = 
      ifelse(inc_adeq_euro <= Hmisc::wtd.quantile(inc_adeq_euro, popw, .2, na.rm = TRUE), 
             inc_adeq_euro, 
             NA_real_),
    
    # inc_adeq_euro_median = Hmisc::wtd.quantile(inc_adeq_euro, popw, probs = 0.5, na.rm = TRUE),
    # # inc_adeq_euro_mean = Hmisc::wtd.mean(inc_adeq_euro, popw, na.rm = TRUE),
    
    arop = 0.6 * Hmisc::wtd.quantile(inc_adeq_euro, popw, probs = 0.5, na.rm = TRUE),
    # arope = ifelse(is.na(rx070), NA, as.numeric(rx070 != 0)),
    
    fgt0_arop = as.integer(inc_adeq_euro < arop),
    fgt1_arop = ifelse(fgt0_arop, 1 - (inc_adeq_euro / arop), 0),
    fgt2_arop = ifelse(fgt0_arop, (1 - (inc_adeq_euro / arop)) ^ 2, 0)
    
  ) |> 
  select(-rb050, -hx040, -hx090, -hx050) |> 
  ungroup() |> 
  glimpse()

# Survey-base aggregates that include design -----------------

survey_design <- svydesign(
  id = ~ psu,
  weights = ~ popw,
  strata = ~ strata,
  nest = TRUE,
  data = df_hh
)

stats_cols <- c("inc_adeq_euro", "inc_pc_euro", "fgt0_arop", "fgt1_arop", "fgt2_arop", 
                "inc_adeq_euro_40", "inc_adeq_euro_20")

# Regular population weights-weighted statistics
df_wtd <- 
  df_hh |> 
  group_by(code, year, subcode) |> 
  summarise(
    gini_adeq = reldist::gini(inc_adeq_euro, weights = popw),
    across(any_of(stats_cols), 
           ~ Hmisc::wtd.mean(., popw, na.rm = TRUE), 
           .names = "wtd.mean_{.col}"),
    
    across(any_of(stats_cols), 
           ~ Hmisc::wtd.quantile(., popw, probs = 0.5, na.rm = TRUE), 
           .names = "wtd.median_{.col}"),
    
    across(any_of(stats_cols),
           ~ Hmisc::wtd.var(., popw, na.rm = TRUE), 
           .names = "wtd.var_{.col}")
  ) |> 
  ungroup() |> 
  glimpse()

# Survey-design corrected means
df_svywtd.mean <-
  svyby(
    formula = paste("~", paste(stats_cols, collapse = "+")) |> as.formula(),
    by = ~ code + year + subcode,
    design = survey_design,
    FUN = svymean,
    deff = TRUE,
    na.rm = TRUE
  ) |>
  as_tibble() |> 
  rename_with(~ paste0("svywtd.mean_", .), .cols = contains(stats_cols)) |> 
  glimpse()

# Survey-design corrected medians
df_svywtd.median <-
  svyby(
    formula = paste("~", paste(stats_cols, collapse = "+")) |> as.formula(),
    by = ~ code + year + subcode,
    design = survey_design,
    FUN = svyquantile,
    deff = TRUE,
    quantiles = 0.5,
    na.rm = TRUE
  )  |>
  as_tibble() |> 
  rename_with(~ paste0("svywtd.median_", .), .cols = contains(stats_cols)) |> 
  glimpse()

# Survey-design corrected variance
df_svywtd.var <-
  stats_cols |> 
  map(~{
    one_var <- .x
    new_se_var <- paste0("se.", one_var)
    svyby(
      formula = paste("~", .x) |> as.formula(),
      by = ~ code + year + subcode,
      design = survey_design,
      FUN = svyvar,
      deff = TRUE,
      na.rm = TRUE
    )  |>
      as_tibble() |> 
      rename_with(~ new_se_var, .cols = contains("se")) |> 
      rename_with(~ paste0("svywtd.var_", .), .cols = contains(stats_cols)) 
  }) |> 
  reduce(left_join, join_by(code, year, subcode))

# Compiling all data ------------------------------------------------------
df_out <- 
  df_wtd |> 
  left_join(df_svywtd.mean, by = join_by(code, year, subcode))|> 
  left_join(df_svywtd.median, by = join_by(code, year, subcode))|> 
  left_join(df_svywtd.var, by = join_by(code, year, subcode)) |> 
  glimpse()

# Saving data -----------------------------------------------------------
write_rds(df_out, file.path(data_aux, "eu-silc.rds"), compress = "gz")
write_csv(df_out, file.path(data_aux, "eu-silc.csv"))

# Save as Excel (.xlsx)
write.xlsx(df_out, file.path(data_aux, "eu-silc.xlsx"),
           sheetName = "EU-SILC Indicators",
           overwrite = TRUE)
