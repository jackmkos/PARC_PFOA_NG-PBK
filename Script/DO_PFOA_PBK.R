  # --------------------------------------------------------------------------- #
  # SCRIPT FOR RUNNING THE PBK MODEL & OBTAINING RESULTS
  # By: Chrysanthi Pachoulide
  # Date: 14-04-2025
  # --------------------------------------------------------------------------- #
  
  rm(list=ls()) 
  
  # Packages
  library(here)
  library(tidyverse)
  library(deSolve)
  library(PKNCA)
  library(pracma)
  
  # Set exposure type, choose between "Oral", "Dermal", "Oral_Dermal" (if exposure is both Oral and Dermal), Inhalation"
  exposure_type = "Oral"
  
  # Set output storage directory
  OUTPUT <- here("Output", format(Sys.Date(), "%Y-%m-%d"), format(Sys.time(), "%H-%M-%S"), exposure_type)
  dir.create(OUTPUT, recursive = TRUE)
  
  # Load input files
  Physio.c <- read_csv(here("Input", "PhysioVariables.csv"))
  Tissue.c <- read_csv(here("Input", "TissueComposition.csv"))
  
  source(here("Script", "RUN_and_OUTPUT.R"))
  
  
  
  RESULTS <- RUNandOUT(

    # Exposure-relevant information
    exposure_type = exposure_type, # type of exposure
    expCONC = 0.024, # ug/kg/day, concentration
    expCONC_Oral = NA, # ug/kg/day, concentration, to be used only when both oral and dermal are used
    expCONC_Dermal = NA, # ug/kg/day, concentration, to be used only when both oral and dermal are used
    Tinput = 1, # for repeated exposure or so, default = 1
    tinterval = 1, # for repeated exposure or so, default = 1
    expSTOP = 10*356, # time in days after which the exposure stopped

    # Subject-relevant information
    expAGE = 20, # years old, age at exposure, if not provided then age argument is not used, physiology is based on BW
    expBW = NA, # kg, if not provided then the BW of the corresponding age and sex is taken; if both BW and Age are not given then a default BW = 70 is taken; if BW is higher than the BW from the lifestage equations then the actual BW overwrites the calculated one
    sex = "M", # sex, either "F" or "M", if none then default is "M"

    # Simulation relevant information
    Tstart = 0, # days, start of the simulation
    Tstop = 30*356, # days, stop of the simulation
    Dt = 1 # days, iteration steps (decrease/increase depending on run time)
                       )

  save(RESULTS, file = here(OUTPUT, "RESULTS_PFOA_PBK.RData"))
  
 