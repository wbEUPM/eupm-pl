
### Main R Script Template
# This scrip can help executing all analysis file by file.
# It uses parameters to activate certain components of the analysis. 

# Set a working directory to the folder with analysis
# Do not change here anything if you use RStudio project (recommended)
setwd(".")

# Restore the environment and load relevant libaray for documents developing
library(renv)
renv::restore()

library(quarto)

# Step 1. Preparing geometries -----------------------------------------
quarto::quarto_render(
  "05-geometries.qmd",
  execute_params = list(
    eval_all = TRUE # Set to TRUE to run all calculations.
  )
)

# Step 3. Computing SILC-based variables ---------------------------------
source("10-SILC-variables.R")

# Step 4. Re-harvesting and compiling right-hand-side variables --------------
quarto::quarto_render(
  "20-other-variables.qmd",
  execute_params = list(
    eval_all = TRUE,# Set to TRUE to run all calculations.
    harvest = FALSE # Setting this to TRUE will start data re-harvesting process, 
    # It is very slow! The best is to supervise it manually.
  )
)

# Step 5. Compiling data into a single database ready-for analysis --------------
quarto::quarto_render(
  "30-data-mege.qmd",
  execute_params = list(
    eval_all = TRUE, # Set to TRUE to run all calculations.
  )
)
