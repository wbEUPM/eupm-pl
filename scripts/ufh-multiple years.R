## Libraries and file paths

pacman::p_load(
  sf, 
  data.table, 
  car, 
  msae,
  sae, 
  survey, 
  spdep,
  knitr, 
  MASS, 
  caret,
  purrr,
  pins,
  rlang,
  emdi,
  gt,
  lmtest, 
  scales,
  viridis, 
  patchwork,
  dplyr, 
  tictoc,
  gt,
  matrixcalc,
  tidyverse,
  modelsummary,
  flextable,
  correlation,
  datawizard
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


## Data

# Direct
pov_direct <- bd_clean |> pin_read("pov_direct") |> 
  rename(id = subcode) |> filter(type == "arop") |> 
  mutate(estimate = "Direct")


## Specify the variable names

area_var <- c("id") 
year_var <- "year"
directpov <- "pov"
directvar <- "vardir"
year_set <- c(2019:2023)

# Auxiliary variables
rhs_dta_0 <- 
  suppressMessages(file.path(root_raw, "data-other", "Auxiliary (potential) variables.xlsx") |>
                     readxl::read_excel()) 

new_names <- 
  rhs_dta_0 |> slice(1) |> t() |> as.data.frame() |>
  rownames_to_column() |> as_tibble() |> 
  mutate(rowname = ifelse(str_detect(rowname, "\\.{2,}"), NA_character_, rowname)) |> 
  fill(rowname) |> 
  mutate(rowname = str_to_lower(rowname),
         rowname = ifelse(!is.na(V1), str_c(rowname,"__", V1), rowname)) |> 
  pull(rowname)

rhs_dta <- 
  rhs_dta_0 |> 
  slice(-1) |> 
  (\(x){colnames(x) <- new_names; x})() |> 
  select(-no., -subregion) |> 
  pivot_longer(
    cols = matches("\\d{4}"), 
    values_transform = \(x) suppressWarnings(as.numeric(x))
  ) |> 
  filter(!is.na(value)) |> 
  separate(name, into = c("var", "year")) |> 
  mutate(year = as.numeric(year)) |> 
  pivot_wider(names_from = "var", values_from = "value")


rhs_meta <-
  file.path(root_raw, "data-other", "Auxiliary (potential) variables.xlsx") |>
  readxl::read_excel(sheet = 2) |> 
  select(var = Variable,
         category = `Category of variable`,
         name = `Variable (in English)`) |> 
  mutate(var = str_to_lower(var))

rhs_dta0 <- 
  rhs_dta |> 
  mutate(
    x101 = x4 / x3, # Share of Working Age pop in the region,
    x102 = x6 / (x3 / x42),
    x103 = x7 / (x3 / x42),
    x104 = x8 / (x3 / x42),
    x105 = x9 / (x3/10000),
    x106 = x12 / (x3 / x42),
    x118 = x12 / (x3 / x42),
    x107 = x14 / (x3 / x42),
    x108 = x23 / (x3/10000),
    x109 = x41 /  (x3/10000),
    x110 = x40 /  (x3/x42),
    x111 = x22 /  (x3/x42),
    x112 = x17 /  (x12),
    x113 = x18 /  (x12),
    x114 = x19 /  (x12),
    x115 = x13 /  (x12),
    x116 = x1 /  (x3 / 1000),
    x117 = x2 /  (x3 / 1000),
    x119 = x29 /  (x4 / 1000),
    x120 = x30 /  (x4 / 1000)
  ) |> 
  select(#-x1, -x2, -x3, 
    -x6, -x7, -x8, -x12, -x13,
    -x17, -x18, -x19, -x31, -x33)

rhs_meta_extra <- 
  tribble(
    ~var, ~name,
    "x101", "Workage pop share",
    "x102", "Sahre fam. poverty ben.",
    "x103", "Sahre fam. homelessness ben.",
    "x104", "Sahre fam. unemployment ben.",
    "x105", "Market places per 10000 person",
    "x106", "Dwelling appartment per family",
    "x107", "Dwelling alowance per family",
    "x108", "Pharmacies per 10000 person",
    "x109", "Care facil. resid. per 10000 person",
    "x110", "Childrent 0-3 per family",
    "x111", "Share of families with Large Fam Card",
    "x112", "Water, share dwel.",
    "x113", "Heating, share dwel.",
    "x114", "Gas, share dwel.",
    "x115", "Rooms per apt.",
    "x116", "Birth rate per 1000",
    "x117", "Deth rate per 1000",
    "x118", "Dwelling rooms per family",
    "x119", "Immigrants per 1000 working pop",
    "x120", "Emmigrants per 1000 working pop"
  )

rhs_meta_0 <- rhs_meta |> bind_rows(rhs_meta_extra)

# Descriptive statistics
rhs_dta0 |>
  pivot_longer(contains("x"), names_to = "var", cols_vary = "slowest") |>
  mutate(var_order = str_extract(var, "\\d{1,4}") |> as.numeric()) |>
  arrange(var_order) |>
  left_join(rhs_meta_0, by = join_by(var)) |>
  mutate(name = str_c(var, " ", name) |> as_factor()) |>
  mutate(year = as.character(year)) |>
  (\(x) bind_rows(x |> mutate(year = "All years"), x))() |>
  filter(year == "All years") |>
  mutate(year = as_factor(year)) |>
  datasummary(
    formula = name  ~ value * (Mean + SD + Min + P50 + Max),
    output = "flextable",
    data = _
  ) |>
  FitFlextableToPage()

var_to_log <- 
  rhs_dta0 |>
  pivot_longer(contains("x"), names_to = "var", cols_vary = "slowest") |>
  group_by(year, var) |> 
  summarise(min = min(value), max = max(value)) |> 
  filter(min * 10 < max & max > 10 | max > 100) |> pull(var) |> unique()

rhs_dta1 <- 
  rhs_dta0 |>
  group_by(year) |>
  mutate(across(all_of(var_to_log), ~ log(.))) |> 
  ungroup()

# Descriptive statistics after transformation is:
rhs_dta1 |> mutate(type = "log") |> 
  bind_rows(rhs_dta0 |> mutate(type = "level")) |> 
  pivot_longer(contains("x"), names_to = "var", cols_vary = "slowest") |>
  mutate(var_order = str_extract(var, "\\d{1,4}") |> as.numeric()) |>
  arrange(var_order) |>
  left_join(rhs_meta_0, by = join_by(var)) |>
  mutate(name = str_c(var, " ", name) |> as_factor()) |>
  mutate(year = as.character(year)) |>
  (\(x) bind_rows(x |> mutate(year = "All years"), x))() |>
  filter(year == "All years") |>
  mutate(year = as_factor(year)) |> #glimpse()
  datasummary(
    formula = name  ~ value * (Mean ) * type,
    output = "flextable",
    data = _
  ) |>
  FitFlextableToPage()

# Correlatoin between poverty rates and key variables

stars.pval <- function(x){
  stars <- c("***", "**", "*", "")
  var <- c(0, 0.01, 0.05, 0.10, 1)
  i <- findInterval(x, var, left.open = T, rightmost.closed = T)
  stars[i]
}

dta_cor <-
  rhs_dta1 |> mutate(type = "log") |> 
  bind_rows(rhs_dta0 |> mutate(type = "level")) |> 
  mutate(year = year + 1) |> 
  left_join(pov_direct |> filter(type == "arop") |> select(id, year, pov)) |> 
  pivot_longer(contains("x"), names_to = "var", cols_vary = "slowest") |>
  filter(!is.na(pov)) |> 
  mutate(var_order = str_extract(var, "\\d{1,4}") |> as.numeric()) |>
  arrange(var_order) |>
  left_join(rhs_meta_0, by = join_by(var)) |>
  mutate(name = str_c(var, " ", name) |> as_factor()) |>
  mutate(year = as.character(year)) |>
  (\(x) bind_rows(x |> mutate(year = "All years"), x))() |>
  filter(year == "All years") |>
  mutate(year = as_factor(year)) |> 
  # filter(str_detect(var, "x10")) |> 
  select(id, year, type, pov, value, name) |> 
  pivot_wider(names_from = name, values_from = value) |> 
  datawizard::data_group(type) |>
  correlation::correlation() |> 
  as_tibble() |> 
  filter(Parameter1 == "pov") |> 
  mutate(Parameter2 = as_factor(Parameter2) |> fct_reorder(abs(r), .desc = T)) |> 
  mutate(stats = str_c(number(r, 0.001), "", stars.pval(p))) |> 
  select(Group, Variable = Parameter2, stats) |>
  pivot_wider(names_from = Group, values_from = stats) |> 
  arrange(Variable)

dta_cor |> 
  flextable() |>
  FitFlextableToPage()




## For the univariate model, we only use a single year; here we use a loop to run multiple years one after the other
# results data frame
pov_fh <- data.frame(matrix(ncol = 8, nrow = 0))
colnames(pov_fh) <- c(area_var, "Direct", "Direct_MSE", "Direct_CV", "FH", "FH_MSE", "FH_CV", "year")

singleyears <- as.character(year_set)
for (i in singleyears){
  survey_dt_tmp <- pov_direct |>
    filter(year == i)
  
  rhs_dt_tmp <- rhs_dta1  |>
    filter(year == as.character(as.numeric(i) - 1))
  
  dta_ufh <-
    rhs_dt_tmp |> 
    left_join(survey_dt_tmp, by = join_by(!!area_var)) |> 
    select(-where(~any(is.na(.))))
  
  candidate_vars <- dta_ufh |> select(contains("x")) |> names()
  
  # Variable selection
  fh_step <- step_wrapper_fh(dt = dta_ufh,
                             xvars = candidate_vars,
                             y = "pov",
                             cor_thresh = 0.7,
                             criteria = "BIC",
                             vardir = "vardir", 
                             transformation = "arcsin", 
                             eff_smpsize = "n_eff")
  
  # Run model
  set.seed(123)
  fh_model <- fh(fixed = formula(fh_step$fixed),
                 vardir = "vardir", 
                 combined_data = as.data.frame(dta_ufh), 
                 domains = "id",
                 method = "ml", 
                 transformation = "arcsin", 
                 backtransformation = "bc",
                 eff_smpsize = "n_eff", 
                 MSE = TRUE, 
                 mse_type = "boot", B = c(200, 0))
  
  # Save results
  pov_fh_tmp <- as.data.frame(estimators(fh_model, MSE = TRUE, CV = TRUE))

  pov_fh_tmp <- pov_fh_tmp |>
    rename(!!area_var[1] := "Domain") |>
    mutate(year = i)
  
  pov_fh <- rbind(pov_fh, pov_fh_tmp)
  
}


bd_out |>
  pin_write(x = pov_fh,
            name = "pov_fh_multipleyears",
            type = "rds")
write.csv(pov_fh, "data/clean-example/pov_fh_multipleyears.csv")
