################################################################################
############# PREPARE POVERTY RATES ESTIMATES FOR THE SAE_DATA #################
################################################################################

pacman::p_load("sf", "tidyverse", "pins", "sae")

#### lets include read in and write up data for pins folders

# first declare pins folder
bd_out <- board_folder("data/clean-example", versioned = TRUE)

income_dt <- readRDS("data/raw/example-data/incomedata_survey.RDS")

rhs_dt <- readRDS("data/raw/example-data/simadmin.RDS")

#### lets simulate enumeration areas into the data
set.seed(123)  # For reproducibility

# Simulate 10 EAs per province, assigned randomly
income_dt <- 
  income_dt %>%
  group_by(prov) %>%
  mutate(ea_num = sample(1:10, n(), replace = TRUE)) %>%
  mutate(
    ea_id = paste0(
      sprintf("%02d", prov),    # pad prov to 2 digits if necessary
      sprintf("%02d", ea_num)   # pad EA number to 2 digits
    )
  ) %>%
  dplyr::select(-ea_num) %>%
  ungroup() |>
  mutate(ea_id = as.integer(ea_id))


### convert both datasets from wide to long to mirror what might be coming out 
### of poland
income_dt <- 
income_dt |>
  pivot_longer(cols = c(starts_with("income"), starts_with("povline")),
               names_to = c(".value", "year"),
               names_pattern = "(income|povline)(\\d{4})") |>
  mutate(year = as.integer(year)) |>
  arrange(prov, ea_id, year)


years <- unique(income_dt$year)


rhs_dt <- 
  lapply(years, 
         function(x){
           
           zz <- 
           rhs_dt |>
             mutate(year = x)
           
           return(zz)
           
         }) |>
  bind_rows()


### save the data to pin board

bd_out |>
  pin_write(x = income_dt,
            name = "pov_direct",
            type = "rds")

write.csv(income_dt, "data/clean-example/pov_direct.csv")

bd_out |>
  pin_write(x = rhs_dt,
            name = "sae_data",
            type = "rds")

write.csv(rhs_dt, "data/clean-example/sae_data.csv")


















