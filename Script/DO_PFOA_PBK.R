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
  
  # !!! FOR DERMAL SEE BELOW NOTES
  # Comment Chrysa: the below section can be used to calculate the dermal concentration based on cosmetic product type and PFOA concentrations in the cosmetic product
  # q = 123.20         # mg/kg/day amount of cosmetic product applied per day, SCCS 2021 table 3 (https://health.ec.europa.eu/document/download/89af1a70-a2b1-44da-a868-e7d80a8e736c_en?filename=sccs_o_250.pdf)
  # fret = 1           # fraction of the cosmetic product retained on the skin, SCCS 2021 table 3 (https://health.ec.europa.eu/document/download/89af1a70-a2b1-44da-a868-e7d80a8e736c_en?filename=sccs_o_250.pdf)
  # Aproduct = 1       # ug/mg amount of PFOA in the cosmetic product
  # fbac = 0.78        # fraction bioaccessible (fraction of PFOA that leaves the cosmetic product and is accessible for absorption) Namazkar et al 2024 DOI: 10.1039/D3EM00461A 
  # CDermal = q*fret*Aproduct*fbac # ug/kd/day
  
  RESULTS <- RUNandOUT(exposure_type = exposure_type, # type of exposure
                       expCONC = 0.048, # ug/kg/day
                       Tinput = 1, # for repeated exposure or so, default = 1
                       tinterval = 1, # for repeated exposure or so, default = 1
                       expAGE = 70, # age at exposure, needed to select physiological parameters
                       expSTOP = 20*356, # time in days after which the exposure stopped 
                       Tstart = 0, # days, start of the simulation
                       Tstop = 50*356, # days, stop of the simulation
                       Dt = 1 # days, iteration steps (decrease/increase depending on run time)
                       )

  save(RESULTS, file = here(OUTPUT, "RESULTS_PFOA_PBK.RData"))