#===========================================================================================================================
# Project: EU
# Description: EU-SILC Indicators for EL in 2021
# Content: This file creates EU-SILC indicators on a sample data. It is expected to be adapted to the needs of each country
#			For analysis over 1 year, adapt the codes from line 13
#			For analysis over multiple years (panel), see example from line 99 with a panel data sample
#==============================================================================================================================
rm(list = ls())

# ========================================================
# PART 1: SINGLE YEAR (NO LOOP)
# ========================================================

# Load necessary packages

library(haven)
library(dplyr)
library(Hmisc)
library(DescTools)
library(data.table)
library(openxlsx)
library(reldist)
library(matrixStats)
library(survey)

# Define paths
data <- "./data/raw/silk-sample"
output <- "./data/temp"

# Load data
df <- read_dta(file.path(data, "dataexample.dta"))

# NB: Important notes before running the analysis:

# 1. Ensure subnational unit codes (db040) are harmonized and numeric.

# 2. The variable "db040" must uniquely identify your geographic level (NUTS2 or NUTS3), depending on your specific analytical requirements.

# 3. Example of converting "db040" from string to numeric (if necessary):
df <- df %>%
  mutate(subcode = if (is.character(db040))
    as.numeric(as.factor(db040))
    else
      db040)

# 4. For comparisons with international poverty lines, convert income data into euros PPP using appropriate conversion factors.

# Grouping and calculation
df <- df %>%
  group_by(code, year, hhid) %>%
  summarise(
    hx090 = first(hx090),
    hx040 = first(hx040),
    hx050 = first(hx050),
    rb050 = first(rb050),
    rx070 = first(rx070),
    psu = first(psu),
    strata = first(strata),
    subcode = first(subcode)
  ) %>%
  ungroup() %>%
  mutate(
    popw = rb050 * hx040,
    # population weights for the group function
    inc_pc_euro = (hx090 * hx050) / hx040,
    # total household disposable income per capita
    inc_adeq_euro = hx090  # total household disposable euros per AE
  )

# Survey design for variance estimates
survey_design <- svydesign(
  id = ~ psu,
  weights = ~ popw,
  strata = ~ strata,
  nest = TRUE,
  data = df
)

# Relative poverty & FGT indices

median_inc <- weightedMedian(df$inc_adeq_euro, w = df$popw, na.rm = TRUE)
arop <- median_inc * 0.6

df <- df %>%
  mutate(
    thres_b20_inc_adeq_euro = median_inc * 0.2,
    thres_b40_inc_adeq_euro = median_inc * 0.4,
    arope = ifelse(is.na(rx070), NA, as.numeric(rx070 != 0)),
    fgt0_arop = as.numeric(inc_adeq_euro < arop),
    fgt1_arop = ifelse(inc_adeq_euro < arop, 1 - (inc_adeq_euro / arop), 0),
    fgt2_arop = ifelse(inc_adeq_euro < arop, (1 - (inc_adeq_euro / arop))^2, 0)
  )


# Measuring the Gini by code-subcode
df <- df %>%
  mutate(subcode = ifelse(is.na(subcode), 0, subcode)) %>%
  group_by(subcode) %>%
  mutate(gini_adeq = gini(inc_adeq_euro, weights = popw)) %>%
  ungroup()

#Measuring required statistics by subcode

df <- df %>%
  group_by(subcode, year) %>%
  mutate(
    mean_inc_adeq_euro = weighted.mean(inc_adeq_euro, popw, na.rm = TRUE),
    p50_inc_adeq_euro = weightedMedian(inc_adeq_euro, w = popw, na.rm = TRUE)
  ) %>%
  ungroup()


# Making the deciles
weighted_decile <- function(x, w) {
  breaks <- unique(Hmisc::wtd.quantile(x, w, seq(0, 1, 0.1), na.rm = TRUE))
  cut(x,
      c(-Inf, breaks[-c(1, length(breaks))], Inf),
      labels = FALSE,
      include.lowest = TRUE)
}

df <- df %>%
  group_by(code, subcode, year) %>%
  mutate(dec_inc_adeq_euro = weighted_decile(inc_adeq_euro, popw)) %>%
  ungroup()

df <- df %>%
  group_by(subcode, year) %>%
  mutate(
    mean_b40_inc_adeq_euro = weighted.mean(inc_adeq_euro[dec_inc_adeq_euro <= 4], popw[dec_inc_adeq_euro <= 4], na.rm = TRUE),
    mean_b20_inc_adeq_euro = weighted.mean(inc_adeq_euro[dec_inc_adeq_euro <= 2], popw[dec_inc_adeq_euro <= 2], na.rm = TRUE)
  ) %>%
  ungroup()

# The variances of the variables
#df <- df %>%
# mutate(
#  obs = 1,
#   var_inc_adeq_euro = inc_adeq_euro,
#  var_fgt0_arop = fgt0_arop,
#  var_fgt1_arop = fgt1_arop,
#   var_fgt2_arop = fgt2_arop
# )

# Estimating the variables with svy design
df$obs <- 1
df$var_fgt0_arop <- NA
df$var_inc_adeq_euro <- NA

# Repeat the survey design to incorporate new variables

survey_design <- svydesign(
  id = ~ psu,
  weights = ~ popw,
  strata = ~ strata,
  nest = TRUE,
  data = df
)

for (sub in unique(df$subcode)) {
  sub_design <- subset(survey_design, subcode == sub)
  var_fgt0 <- as.numeric(attr(svymean( ~ fgt0_arop, sub_design, na.rm = TRUE), "var"))
  var_inc  <- as.numeric(attr(svymean( ~ inc_adeq_euro, sub_design, na.rm = TRUE), "var"))
  df$var_fgt0_arop[df$subcode == sub] <- var_fgt0
  df$var_inc_adeq_euro[df$subcode == sub] <- var_inc
}


# Weighted statistics

df_final <- df %>%
  group_by(code, subcode, year) %>%
  summarise(
    across(c(starts_with("fgt"), "arope"), ~ weighted.mean(.x, popw, na.rm = TRUE)),
    across(
      c(
        starts_with("thres"),
        starts_with("p50"),
        starts_with("mean_inc"),
        "gini_adeq",
        starts_with("var")
      ),
      ~ first(na.omit(.x))
    ),
    across(starts_with("mean_b"), ~ max(.x, na.rm = TRUE)),
    popw = sum(popw, na.rm = TRUE),
    obs = sum(obs, na.rm = TRUE),
    .groups = 'drop'
  )

df_final <- df_final %>%
  mutate(across(c(
    starts_with("fgt"),
    starts_with("gini"),
    starts_with("var_fgt")
  ), ~ round(.x, 4))) %>%
  select(
    code,
    subcode,
    year,
    popw,
    obs,
    fgt0_arop,
    var_fgt0_arop,
    mean_inc_adeq_euro,
    var_inc_adeq_euro,
    p50_inc_adeq_euro,
    everything()
  )


# Export final dataset as RDS and excel

saveRDS(df_final, file.path(output, "EU_SILC_Indicators.rds"))
write.xlsx(
  df_final,
  file.path(output, "EU_SILC_Indicators.xlsx"),
  sheetName = "EU-SILC Indicators",
  overwrite = TRUE
)




# ========================================================
# PART 2: MULTIPLE YEARS (WITH LOOP)
# ========================================================
df_combined <- NULL
for (year in 2021:2023) {
  message("Processing year: ", year)
  file_name <- paste0("dataexample", year, ".dta")
  full <- file.path(data, file_name)
  
  if (!file.exists(full)) {
    message("File not found: ", full)
    next
  }
  
  df <- read_dta(full)
  
  df <- df %>%
    mutate(subcode = if (is.character(db040))
      as.numeric(as.factor(db040))
      else
        db040)
  
  df <- df %>%
    group_by(code, year, hhid) %>%
    summarise(
      hx090 = first(hx090),
      hx040 = first(hx040),
      hx050 = first(hx050),
      rb050 = first(rb050),
      rx070 = first(rx070),
      psu = first(psu),
      strata = first(strata),
      subcode = first(subcode)
    ) %>%
    ungroup() %>%
    mutate(
      popw = rb050 * hx040,
      # population weights for the group function
      inc_pc_euro = (hx090 * hx050) / hx040,
      # total household disposable income per capita
      inc_adeq_euro = hx090  # total household disposable euros per AE
    )
  
  # Survey design for variance estimates
  survey_design <- svydesign(
    id = ~ psu,
    weights = ~ popw,
    strata = ~ strata,
    nest = TRUE,
    data = df
  )
  
  # Relative poverty & FGT indices
  
  median_inc <- weightedMedian(df$inc_adeq_euro, w = df$popw, na.rm = TRUE)
  arop <- median_inc * 0.6
  
  df <- df %>%
    mutate(
      thres_b20_inc_adeq_euro = median_inc * 0.2,
      thres_b40_inc_adeq_euro = median_inc * 0.4,
      arope = ifelse(is.na(rx070), NA, as.numeric(rx070 != 0)),
      fgt0_arop = as.numeric(inc_adeq_euro < arop),
      fgt1_arop = ifelse(inc_adeq_euro < arop, 1 - (inc_adeq_euro / arop), 0),
      fgt2_arop = ifelse(inc_adeq_euro < arop, (1 - (
        inc_adeq_euro / arop
      ))^2, 0)
    )
  
  
  # Measuring the Gini by code-subcode
  df <- df %>%
    mutate(subcode = ifelse(is.na(subcode), 0, subcode)) %>%
    group_by(subcode) %>%
    mutate(gini_adeq = gini(inc_adeq_euro, weights = popw)) %>%
    ungroup()
  
  #Measuring required statistics by subcode
  
  df <- df %>%
    group_by(subcode, year) %>%
    mutate(
      mean_inc_adeq_euro = weighted.mean(inc_adeq_euro, popw, na.rm = TRUE),
      p50_inc_adeq_euro = weightedMedian(inc_adeq_euro, w = popw, na.rm = TRUE)
    ) %>%
    ungroup()
  
  
  # Making the deciles
  weighted_decile <- function(x, w) {
    breaks <- unique(Hmisc::wtd.quantile(x, w, seq(0, 1, 0.1), na.rm = TRUE))
    cut(x,
        c(-Inf, breaks[-c(1, length(breaks))], Inf),
        labels = FALSE,
        include.lowest = TRUE)
  }
  
  df <- df %>%
    group_by(code, subcode, year) %>%
    mutate(dec_inc_adeq_euro = weighted_decile(inc_adeq_euro, popw)) %>%
    ungroup()
  
  df <- df %>%
    group_by(subcode, year) %>%
    mutate(
      mean_b40_inc_adeq_euro = weighted.mean(inc_adeq_euro[dec_inc_adeq_euro <= 4], popw[dec_inc_adeq_euro <= 4], na.rm = TRUE),
      mean_b20_inc_adeq_euro = weighted.mean(inc_adeq_euro[dec_inc_adeq_euro <= 2], popw[dec_inc_adeq_euro <= 2], na.rm = TRUE)
    ) %>%
    ungroup()
  
  # The variances of the variables
  #df <- df %>%
  # mutate(
  #  obs = 1,
  #   var_inc_adeq_euro = inc_adeq_euro,
  #  var_fgt0_arop = fgt0_arop,
  #  var_fgt1_arop = fgt1_arop,
  #   var_fgt2_arop = fgt2_arop
  # )
  
  # Estimating the variables with svy design
  df$obs <- 1
  df$var_fgt0_arop <- NA
  df$var_inc_adeq_euro <- NA
  
  # Repeat the survey design to incorporate new variables
  
  survey_design <- svydesign(
    id = ~ psu,
    weights = ~ popw,
    strata = ~ strata,
    nest = TRUE,
    data = df
  )
  
  for (sub in unique(df$subcode)) {
    sub_design <- subset(survey_design, subcode == sub)
    var_fgt0 <- as.numeric(attr(svymean( ~ fgt0_arop, sub_design, na.rm = TRUE), "var"))
    var_inc  <- as.numeric(attr(svymean( ~ inc_adeq_euro, sub_design, na.rm = TRUE), "var"))
    df$var_fgt0_arop[df$subcode == sub] <- var_fgt0
    df$var_inc_adeq_euro[df$subcode == sub] <- var_inc
  }
  
  
  # Weighted statistics
  
  df_final <- df %>%
    group_by(code, subcode, year) %>%
    summarise(
      across(
        c(starts_with("fgt"), "arope"),
        ~ weighted.mean(.x, popw, na.rm = TRUE)
      ),
      across(
        c(
          starts_with("thres"),
          starts_with("p50"),
          starts_with("mean_inc"),
          "gini_adeq",
          starts_with("var")
        ),
        ~ first(na.omit(.x))
      ),
      across(starts_with("mean_b"), ~ max(.x, na.rm = TRUE)),
      popw = sum(popw, na.rm = TRUE),
      obs = sum(obs, na.rm = TRUE),
      .groups = 'drop'
    )
  
  df_final <- df_final %>%
    mutate(across(c(
      starts_with("fgt"),
      starts_with("gini"),
      starts_with("var_fgt")
    ), ~ round(.x, 4))) %>%
    select(
      code,
      subcode,
      year,
      popw,
      obs,
      fgt0_arop,
      var_fgt0_arop,
      mean_inc_adeq_euro,
      var_inc_adeq_euro,
      p50_inc_adeq_euro,
      everything()
    )
  
  
  df_combined <- bind_rows(df_combined, df_final)
  message("Processed and appended data for year ", year)
  
}
# Save as RDS (recommended)
saveRDS(df_combined,
        file.path(output, "EU_SILC_Indicators_allyears.rds"))

# Save as Excel (.xlsx)
write.xlsx(
  df_combined,
  file.path(output, "EU_SILC_Indicators_allyears.xlsx"),
  sheetName = "EU-SILC Indicators",
  overwrite = TRUE
)

rm(list = ls())
