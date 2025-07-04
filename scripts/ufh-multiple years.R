## Libraries and file paths

pacman::p_load(
  sf, 
  data.table, 
  msae,
  sae, 
  pins,
  emdi, 
  scales, 
  patchwork,
  tictoc,
  tidyverse
)

# Loading locally-developed
list.files("R", pattern="*.R$", full.names=TRUE, ignore.case=TRUE) |> 
  walk(source)

# Raw data root
root_raw <- "./data/raw"
root_temp <- "./data/temp"
root_clean <- "./data/clean"

# Data-storage boards
bd_raw <- root_raw |> file.path("api") |> board_folder(versioned = T)
bd_aux <- root_temp |> board_folder(versioned = T)
bd_clean <- root_clean |> board_folder(versioned = T)
bd_out <- "output/res-40-fh" |> board_folder(versioned = T)


# Data loading -------------------------------------------------------------


## Variable names ---------------------------------------------------------
area_var <- c("id") 
year_var <- "year"
directpov <- "pov"
directvar <- "vardir"
year_set <- c(2019:2023)

## Direct -----------------------------------------------------------------
pov_direct <- bd_clean |> pin_read("pov_direct") |> 
  rename(id = subcode) |> mutate(estimate = "Direct")

## RHS variables ---------------------------------------------------------
rhs_dta <- bd_clean |> pin_read("rhs_dta")
rhs_meta <- bd_clean |> pin_read("rhs_meta")

### Log transformation ---------------------------------------------------
var_to_log <- 
  rhs_dta |>
  pivot_longer(contains("x"), names_to = "var", cols_vary = "slowest") |>
  group_by(year, var) |> 
  summarise(min = min(value), max = max(value)) |> 
  filter(min * 10 < max & max > 10 | max > 100) |> pull(var) |> unique()

rhs_dta1 <- 
  rhs_dta |>
  group_by(year) |>
  mutate(across(all_of(var_to_log), ~ log(.))) |> 
  ungroup() |> 
  select(-x15)

# UFH in a loop ---------------------------------------------------------------
pov_fh <- tibble()

for (method in c("reml", "ml")) {
for (criteria in c("AIC", "BIC")) {
for (type_poverty in unique(pov_direct$type)) {
for (i in as.character(year_set)){
  survey_dt_tmp <- pov_direct |> filter(year == i, type == type_poverty)
  rhs_dt_tmp <- rhs_dta1  |> filter(year == as.character(as.numeric(i)))
  
  dta_ufh <-
    rhs_dt_tmp |> 
    select(-year) |> 
    left_join(survey_dt_tmp, by = join_by(!!area_var)) |> 
    select(-where(~any(is.na(.))))
  
  candidate_vars <- dta_ufh |> select(contains("x")) |> names()
  # Variable selection
  fh_step <- step_wrapper_fh(dt = dta_ufh,
                             xvars = candidate_vars,
                             y = "pov",
                             cor_thresh = 0.7,
                             criteria = criteria,
                             vardir = "vardir", 
                             transformation = "arcsin", 
                             eff_smpsize = "n_eff",
                             method = "ml")
  
  # Run model
  set.seed(123)
  fh_model <- fh(fixed = formula(fh_step$fixed),
                 vardir = "vardir", 
                 combined_data = as.data.frame(dta_ufh), 
                 domains = "id",
                 method = method, 
                 transformation = "arcsin", 
                 backtransformation = "bc",
                 eff_smpsize = "n_eff", 
                 MSE = TRUE, 
                 mse_type = "boot", B = c(200, 0))
  
  bd_aux |> pin_write(
    fh_model,
    name = paste0("fh_model__", type_poverty, "__", criteria, "__", i),
    type = "rds"
  )
  
  # Save results
  # pov_fh_tmp <- 
    as_tibble(estimators(fh_model, MSE = TRUE, CV = TRUE)) |> 
    select(!!area_var := Domain, pov = FH, MSE = FH_MSE, CV = FH_CV) |>
    mutate(
      year = i,
      estimate = "UFH",
      type = type_poverty,
      criteria = criteria,
      method = method,
      re = fh_model$model$random_effects[,1],
      n_var = length(coef(fh_model)) - 1,
      r2_adj = fh_model$model$model_select$AdjR2,
      r2_fh = fh_model$model$model_select$FH_R2
    )
  
  pov_fh <- bind_rows(pov_fh, pov_fh_tmp)
  
}
}
}
}

bd_clean |>
  pin_write(x = pov_fh,
            name = "pov_fh_multipleyears",
            type = "rds")


# Selected Diagnosis ------------------------------------------------------


#### Post estimation analysis
fh_model_2019 <- bd_out |> pin_read("fh_model__arop__2019")
fh_model_2020 <- bd_out |> pin_read("fh_model__arop__2020")
fh_model_2021 <- bd_out |> pin_read("fh_model__arop__2021")
fh_model_2022 <- bd_out |> pin_read("fh_model__arop__2022")
fh_model_2023 <- bd_out |> pin_read("fh_model__arop__2023")

for (i in singleyears){
  pdf(paste0("Plots-", i ,".pdf"))
 # summary(eval(parse(text = paste0("fh_model_",i))))
  plot(eval(parse(text = paste0("fh_model_",i))))
 # compare_plot(val(parse(text = paste0("fh_model_",i))), MSE = TRUE, CV = TRUE)
 # compare(eval(paste0("fh_model_",i)))
  dev.off()
}

for (i in singleyears){
  pdf(paste0("Compare-plots-", i ,".pdf"))
  # summary(eval(parse(text = paste0("fh_model_",i))))
  #plot(eval(parse(text = paste0("fh_model_",i))))
   compare_plot(eval(parse(text = paste0("fh_model_",i))), MSE = TRUE, CV = TRUE)
  # compare(eval(paste0("fh_model_",i)))
  dev.off()
}


summary(fh_model_2019)
plot(fh_model_2019)
compare_plot(fh_model_2019, MSE = TRUE, CV = TRUE)
compare(fh_model_2019)

summary(fh_model_2020)
plot(fh_model_2020)
compare_plot(fh_model_2020, MSE = TRUE, CV = TRUE)
compare(fh_model_2020)

summary(fh_model_2021)
plot(fh_model_2021)
compare_plot(fh_model_2021, MSE = TRUE, CV = TRUE)
compare(fh_model_2021)

summary(fh_model_2022)
plot(fh_model_2022)
compare_plot(fh_model_2022, MSE = TRUE, CV = TRUE)
compare(fh_model_2022)

summary(fh_model_2023)
plot(fh_model_2023)
compare_plot(fh_model_2023, MSE = TRUE, CV = TRUE)
compare(fh_model_2023)
