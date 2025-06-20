## Libraries and file paths

pak_list <- c("sf", "data.table", "tidyverse", "car", 
              "sae", "survey", "spdep", "knitr", "MASS", "caret",
              "caret", "purrr", "pins", "gt",  "scales",
              "viridis", "emdi", "rlang", "matrixcalc")

sapply(pak_list,
       library,
       character.only = TRUE)


# Loading locally-developed
list.files("R", pattern="*.R$", full.names=TRUE, ignore.case=TRUE) |> 
  walk(source)

# Raw data root
root_raw <- "./data/raw"
root_temp <- "./data/temp"
root_clean <- "./data/clean-example"

# Data-storage boards
bd_raw <- root_raw |> file.path("api") |> board_folder(versioned = T)
bd_aux <- root_temp |> board_folder(versioned = T)
bd_clean <- root_clean |> board_folder(versioned = T)
bd_out <- board_folder("data/clean-example", versioned = TRUE)

## Data

survey_dt <- 
  bd_clean |>
  pin_read("pov_direct")


rhs_dt <- 
  bd_clean |>
  pin_read("sae_data")


shp_dt <- 
  bd_clean |>
  pin_read("geometries")

## Specify the variable names

area_vars <- c("prov", "provlab") ### both variables are at the same level. If the levels vary, you would need to combine both variables for effect use

cluster_var <- "ea_id"

weight_var <- "weight"

year_var <- "year"

outcome_var <- "income"

povline_var <- "povline"

candidate_vars <- colnames(rhs_dt)[!colnames(rhs_dt) %in% c(area_vars, year_var)]

## For the univariate model, we only use a single year; here we use a loop to run multiple years one after the other
# results data frame
pov_fh <- data.frame(matrix(ncol = 8, nrow = 0))
colnames(pov_fh) <- c(area_vars[1], "Direct", "Direct_MSE", "Direct_CV", "FH", "FH_MSE", "FH_CV", "year")

singleyears <- as.character(unique(survey_dt$year))
for (i in singleyears){
  survey_dt_tmp <- survey_dt |>
    filter(year == i)
  
  rhs_dt_tmp <- rhs_dt |>
    filter(year == i)
  
  ## calculate sample size for each province
  sampsize_dt <- 
    survey_dt_tmp |>
    group_by(!!!syms(area_vars[1])) |>
    summarize(N = n())
  
  ## computation of direct estimates and their variances (the poverty line is already included within the data)
  ## creating the poverty indicator
  survey_dt_tmp <- 
    survey_dt_tmp |>
    mutate(pov_indicator = ifelse(income < povline, 1, 0))
  
  
  ### creating a survey object
  design_obj <- survey::svydesign(ids = eval(expr(~!!sym(cluster_var))), 
                                  weights = eval(expr(~!!sym(weight_var))), 
                                  data = survey_dt_tmp)
  
  var_dt <- survey::svyby(~pov_indicator, by=eval(expr(~!!sym(area_vars[1]))), design = design_obj, FUN = survey::svymean)
  
  
  direct_dt <- 
    var_dt |>
    rename(direct_povrate = "pov_indicator") |>
    rename(SD = "se") |>
    mutate(vardir = SD^2) |>
    mutate(CV = SD / direct_povrate) |>
    merge(sampsize_dt, 
          by = area_vars[[1]]) |>
    mutate(var_SRS = direct_povrate * (1 - direct_povrate) / N) |>
    mutate(deff = vardir / var_SRS) |>
    mutate(n_eff = N/deff)
  
  ## set zero variance to OOS
  direct_dt <- direct_dt[complete.cases(direct_dt), ]
  
  
  ## Variance smoothing
  var_smooth <- varsmoothie_king(domain = direct_dt[[area_vars[1]]],
                                 direct_var = direct_dt$vardir,
                                 sampsize = direct_dt$N)
  
  direct_dt <- var_smooth |> merge(direct_dt, by.x = "Domain",
                                   by.y = area_vars[1])
  
  # Replace the variances that are zero with their smoothed counterparts
  direct_dt <- 
    direct_dt |>
    mutate(across(
      starts_with("v_"),
      ~ if_else(abs(.x) <= 1e-4, get(paste0("vsv", str_remove(cur_column(), "^v"))), .x),
      .names = "{.col}"
    ))
  
  # Combine data sets
  fh_dt <- merge(direct_dt, rhs_dt_tmp,
                 by.x = "Domain", by.y = area_vars[[1]],
                 all = TRUE)
  
  # Remove NAs in auxiliary variables
  rowsNAcovariates <- rowSums(sapply(fh_dt[,..candidate_vars], is.na))
  fh_dt <- fh_dt[rowsNAcovariates == 0, ]
  
  # Variable selection
  fh_step <- step_wrapper_fh(dt = fh_dt,
                             xvars = candidate_vars,
                             y = "direct_povrate",
                             cor_thresh = 0.7,
                             criteria = "BIC",
                             vardir = "vardir", 
                             transformation = "arcsin", 
                             eff_smpsize = "n_eff")
  
  # Run model
  set.seed(123)
  fh_model <- fh(fixed = formula(fh_step$fixed),
                 vardir = "vardir", 
                 combined_data = fh_dt, 
                 domains = "Domain",
                 method = "ml", 
                 transformation = "arcsin", 
                 backtransformation = "bc",
                 eff_smpsize = "n_eff", 
                 MSE = TRUE, 
                 mse_type = "boot", B = c(200, 0))
  
  # Save results
  pov_fh_tmp <- as.data.frame(estimators(fh_model, MSE = TRUE, CV = TRUE))

  pov_fh_tmp <- pov_fh_tmp |>
    rename(!!area_vars[1] := "Domain") |>
    mutate(year = i)
  
  pov_fh <- rbind(pov_fh, pov_fh_tmp)
  
}


bd_out |>
  pin_write(x = pov_fh,
            name = "pov_fh_multipleyears",
            type = "rds")
write.csv(pov_fh, "data/clean-example/pov_fh_multipleyears.csv")
