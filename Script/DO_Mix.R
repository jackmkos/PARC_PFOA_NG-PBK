# --------------------------------------------------------------------------- #
# SCRIPT FOR RUNNING THE MODEL WITH VARIOUS PFAS IN SERIES 
# By: Jack Koster 
# Date: 22-05-2025
# --------------------------------------------------------------------------- #

# Packages
library(here)
library(tidyverse)
library(deSolve)
library(PKNCA)
library(pracma)
library(showtext)

# PFAS Loop ####
PFAS_list <- c("PFOA", "PFOS", "PFNA", "PFHxS")

for (compound in PFAS_list) {
  PFAS <- compound
  Mix <- TRUE #Flag for toggling rm() in DO_PFOA_PBK.R
  message("Running PBK model for ", PFAS)
  source(here("Script", "DO_PFOA_PBK.R"), local = FALSE)
  rm(All_PFAS)
}

# Results ####

## Recursively list all RESULTS_*_PBK.RData files in Output
RESULTS_total <- list.files(
  path = here("Output"),
  pattern = "^RESULTS_.*_PBK\\.RData$",
  recursive = TRUE,
  full.names = TRUE
)

## Loop through each file, load into its own RESULTS_PFAS object
for (file_path in RESULTS_total) {
  filename <- basename(file_path)
  pfas <- sub("^RESULTS_(.*)_PBK\\.RData$", "\\1", filename)
  tmp_env <- new.env()
  load(file_path, envir = tmp_env)
  assign(paste0("RESULTS_", pfas), tmp_env$RESULTS)
}

## Load manually
# temp_env <- new.env()
# load("RESULTS_PFOA_PBK.RData", envir = temp_env)
# RESULTS_PFOA <- temp_env$RESULTS

## Results extraction
PFOA.c <- RESULTS_PFOA$ANALYSED_data[[1]]$HalfLife
PFOS.c <- RESULTS_PFOS$ANALYSED_data[[1]]$HalfLife
PFNA.c <- RESULTS_PFNA$ANALYSED_data[[1]]$HalfLife
PFHxS.c <- RESULTS_PFHxS$ANALYSED_data[[1]]$HalfLife

# For DO_PFOA_PNK.R:
if (!exists("Mix")) { # check for DO_Mix.R
  rm(list = ls())
}
  
# Choose PFAS compound ("PFOA", "PFOS", "PFNA", or "PFHxS")
if (!exists("Mix")) # Override for running All_PFAS
PFAS = "PFOA"
