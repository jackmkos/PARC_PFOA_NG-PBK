  # --------------------------------------------------------------------------- #
  # SCRIPT FOR RUNNING THE PBK MODEL & OBTAINING RESULTS
  # By: Chrysanthi Pachoulide
  # Date: 14-04-2025
  # --------------------------------------------------------------------------- #
  
  rm(list=ls()) # to clear out the global environment
  
  # Packages
  library(here)
  library(tidyverse)
  library(deSolve)
  library(PKNCA)
  library(pracma)
  
  # Set exposure type, choose between "Oral", "Dermal", "Inhalation"
  exposure_type = "Oral"
  
  # Set output storage directory
  OUTPUT <- here("Output", format(Sys.Date(), "%Y-%m-%d"), format(Sys.time(), "%H-%M-%S"), exposure_type)
  dir.create(OUTPUT, recursive = TRUE)
  
  # Load input files
  Physio.c <- read_csv(here("Input", "PhysioVariables.csv"))
  Tissue.c <- read_csv(here("Input", "TissueComposition.csv"))
  
  source(here("Script", "RUN_and_OUTPUT.R"))
  
  RESULTS <- RUNandOUT(exposure_type = "Oral", # type of exposure
                       expCONC = 0.048, # ug/kg/day
                       Tinput = 1, # for repeated exposure or so, default = 1
                       tinterval = 1, # for repeated exposure or so, default = 1
                       expAGE = 70, # age at exposure, needed to select physiological parameters
                       expSTOP = 1, # time in days after which the exposure stopped 
                       Tstart = 0, # days, start of the simulation
                       Tstop = 2*356, # days, stop of the simulation
                       Dt = 1/10 # days, iteration steps (decrease/increase depending on run time)
                       )

  save(RESULTS, file = here(OUTPUT, "RESULTS_PFOA_PBK.RData"))