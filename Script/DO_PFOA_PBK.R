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
  library(showtext)
  font_add(family = "Garamond", regular = "GARA.TTF")
  showtext_auto()
  
  # Set output storage directory
  OUTPUT <- here("Output", format(Sys.Date(), "%Y-%m-%d"), format(Sys.time(), "%H-%M-%S"))
  dir.create(OUTPUT, recursive = TRUE)
  
  # Choose if physiology should change with age ("Yes" to include physiological changes due to age and "No" to assume the same physiology over time)
  Lifestage = "Yes" 
  
  # Choose to include population or individual exposure ("Yes" to include population based and "No" to only run the model for one person)
  Population = "Yes" 
  
  # Load files
  Physio.c <- read_csv(here("Input", "PhysioVariables.csv"))
  Tissue.c <- read_csv(here("Input", "TissueComposition.csv"))
  source(here("Script", "RUN_and_OUTPUT.R"))
  

  # Load input ----
  
  if(Population == "Yes"){
    
    INPUT_dummy <- read.csv(here("Input", "INPUT_dummy.csv")) 
    Input <- INPUT_dummy %>% filter(Study == "dummy") # Choose the study you're interested in
    
    
    if(Lifestage == "Yes" && any(is.na(INPUT_dummy$expAGE))){
      
      warning("Removing ", sum(is.na(INPUT_dummy$expAGE)), " samples as they had NA(s) as expAGE; exposure AGE is needed to run the lifestage model")
      Input <- filter (INPUT_dummy, !is.na(expAGE))
      Input <- Input  %>% filter(Study == "dummy")
    } 
    
    nPeople <- as.numeric(nrow(Input)) # number of people
    Pop.RESULTS <- list(
      CALC_Parameters = vector("list", nPeople), # list of length nPeople
      OUT_RAW_data = vector("list", nPeople),
      ANALYSED_data = vector("list", nPeople),
      OUT_Plots = vector("list", nPeople) # could be removed if it's too heavy for R
    )
    
    
    } else {
    
    # Set exposure type, choose between "Oral", "Dermal", "Oral_Dermal" (if exposure is both Oral and Dermal), Inhalation"
    exposure_type = "Oral"
    
    # Add input information
    
    # Exposure-relevant information
    # Current input is that of the Abraham study; Abraham et al. 2024 https://doi.org/10.1016/j.envint.2024.109047 
    exposure_type = exposure_type # type of exposure
    exp_Oral = 0.048 # ug/kg/day, (for Abraham: 3.96/BW of 82Kg),to be used only when both oral and dermal are used
    exp_Dermal = 0 # ug/kg/day, to be used only when both oral and dermal are used
    exp = exp_Oral + exp_Dermal # ug/kg/day 
    Tinput = 1 # for repeated exposure or so default = 1
    tinterval = 1 # for repeated exposure or so default = 1
    expSTOP = 1 # time in days after which the exposure stopped
    
    # Subject-relevant information
    expAGE = 67 # years, old age at exposure if not provided then age argument is not used physiology is based on BW
    expBW = 82 # kg, if not provided then the BW of the corresponding age and sex is taken; if both BW and Age are not given then a default BW = 70 is taken; if BW is higher than the BW from the lifestage equations then the actual BW overwrites the calculated one
    sex = "M" # sex either "F" or "M" if none then default is "M"
    
    # Simulation relevant information
    Tstart = 0 # days, start of the simulation
    Tstop = 450 # days, stop of the simulation
    Dt = 1/10 # days, iteration steps (decrease/increase depending on run time)
    
    
    # List for storing results
    RESULTS <- list(
      CALC_Parameters = vector("list", 1), 
      OUT_RAW_data = vector("list", 1),
      ANALYSED_data = vector("list", 1),
      OUT_Plots = vector("list", 1) 
    )
    
    if(Lifestage == "Yes"){
      
      expBW <- NA # forcing BW to be NA as BW is taken based on age in the lifestage model
      if(is.na(expAGE)){
        stop("expAGE required to run the lifestage option")
      }
      
    }
    
  }

 
  MODEL_OUTPUT <- if(Population == "Yes"){ # Population exposure
    
    if(Lifestage == "Yes"){ # Population and lifestage
      
      for (i in 1:nPeople) {
        
        # Run the model per person
        Pop.MODEL_OUTPUT <- RUNandOUT_lifestage(
          
          exposure_type = as.character(Input[i, "exposure_type"]),            # Exposure type, choose between "Oral", "Dermal", "Oral_Dermal" (if exposure is both Oral and Dermal), Inhalation"
          exp = as.numeric(Input[i, "exp"]),                           # ug/kg/day,(total exposure)
          exp_Oral = ifelse(is.na(Input[i, "exp_Oral"]), 0,
                                 as.numeric(Input[i, "exp_Oral"])),        # ug/kg/day, to be used only when both oral and dermal are used
          exp_Dermal = ifelse(is.na(Input[i, "exp_Dermal"]), 0,
                                   as.numeric(Input[i, "exp_Dermal"])),    # ug/kg/day, to be used only when both oral and dermal are used
          Tinput = ifelse(is.na(Input[i, "Tinput"]), 1,
                           as.numeric(Input[i, "Tinput"])),          # for repeated exposure or so default = 1
          tinterval = ifelse(is.na(Input[i, "tinterval"]), 1,
                              as.numeric(Input[i, "tinterval"])),    # for repeated exposure or so default = 1
          expSTOP = as.numeric(Input[i, "expSTOP"]),                # days, time after which the exposure stopped
          
          # Subject-relevant information
          expAGE = ifelse(is.na(Input[i, "expAGE"]), stop("Error as expAGE required to run the lifestage model"),
                           as.numeric(Input[i, "expAGE"])),     # years old, age at exposure
          expBW = if(!is.na(Input[i, "expBW"])){
            warning("expBW overwritten as NA as it should be automatically calculated based on age in the lifestage model")      # kg
            NA
          } else { Input[i, "expBW"] },
          sex = ifelse(is.na(Input[i, "sex"]), "M",
                        as.character(Input[i, "sex"])),         # sex, either "F" or "M" if none then default is "M"
          
          # Simulation relevant information
          Tstart = as.numeric(Input[i, "Tstart"]),          # days, start of the simulation
          Tstop = as.numeric(Input[i, "Tstop"]),      # days, stop of the simulation
          Dt = as.numeric(Input[i, "Dt"])               # days, iteration steps (decrease/increase depending on run time)
          
        )
        
        # Collect population results together
        Pop.RESULTS$CALC_Parameters[[i]] <- Pop.MODEL_OUTPUT$CALC_Parameters
        Pop.RESULTS$OUT_RAW_data[[i]] <- Pop.MODEL_OUTPUT$OUT_RAW_data
        Pop.RESULTS$ANALYSED_data[[i]] <- Pop.MODEL_OUTPUT$ANALYSED_data
        Pop.RESULTS$OUT_Plots[[i]] <- Pop.MODEL_OUTPUT$OUT_Plots # could be removed if it's too heavy for R
        
        
        
      }
      
    } else { # Population but no lifestage
      
      # Run the model for the population
      for (i in 1:nPeople) {
        
        # Run the model per person
        Pop.MODEL_OUTPUT <- RUNandOUT(
          
          exposure_type = as.character(Input[i, "exposure_type"]),            # Exposure type, choose between "Oral", "Dermal", "Oral_Dermal" (if exposure is both Oral and Dermal), Inhalation"
          exp = as.numeric(Input[i, "exp"]),                           # ug/kg/day, (total exposure)
          exp_Oral = ifelse(is.na(Input[i, "exp_Oral"]), 0,
                                 as.numeric(Input[i, "exp_Oral"])),        # ug/kg/day, to be used only when both oral and dermal are used
          exp_Dermal = ifelse(is.na(Input[i, "exp_Dermal"]), 0,
                                   as.numeric(Input[i, "exp_Dermal"])),    # ug/kg/day, to be used only when both oral and dermal are used
          Tinput = ifelse(is.na(Input[i, "Tinput"]), 1,
                           as.numeric(Input[i, "Tinput"])),          # for repeated exposure or so default = 1
          tinterval = ifelse(is.na(Input[i, "tinterval"]), 1,
                              as.numeric(Input[i, "tinterval"])),    # for repeated exposure or so default = 1
          expSTOP = as.numeric(Input[i, "expSTOP"]),                # time in days after which the exposure stopped
          
          # Subject-relevant information
          expAGE = ifelse(is.na(Input[i, "expAGE"]), NA,
                           as.numeric(Input[i, "expAGE"])),     # years old age at exposure if not provided then age argument is not used physiology is based on BW
          expBW = ifelse(is.na(Input[i, "expBW"]), NA,
                          as.numeric(Input[i, "expBW"])),       # kg, if not provided then the BW of the corresponding age and sex is taken; if both BW and Age are NA, then a default BW = 70 is taken; if BW is higher than the BW from the lifestage equations then the actual BW overwrites the calculated one
          sex = ifelse(is.na(Input[i, "sex"]), "M",
                        as.character(Input[i, "sex"])),         # sex, either "F" or "M" if none then default is "M"
          
          # Simulation relevant information
          Tstart = as.numeric(Input[i, "Tstart"]),          # days, start of the simulation
          Tstop = as.numeric(Input[i, "Tstop"]),      # days, stop of the simulation
          Dt = as.numeric(Input[i, "Dt"])               # days, iteration steps (decrease/increase depending on run time)
          
        )
        
        # Collect population results together
        Pop.RESULTS$CALC_Parameters[[i]] = Pop.MODEL_OUTPUT$CALC_Parameters
        Pop.RESULTS$OUT_RAW_data[[i]] = Pop.MODEL_OUTPUT$OUT_RAW_data
        Pop.RESULTS$ANALYSED_data[[i]] = Pop.MODEL_OUTPUT$ANALYSED_data
        Pop.RESULTS$OUT_Plots[[i]] = Pop.MODEL_OUTPUT$OUT_Plots # could be removed if it's too heavy for R
        
      }
      
    }
    
    
  } else { # Individual
    
    if(Lifestage == "Yes"){ # Individual and lifestage
      
      RUNandOUT_lifestage(
        
        # Exposure-relevant information
        exposure_type = exposure_type, 
        exp = exp, 
        exp_Oral = exp_Oral, 
        exp_Dermal = exp_Dermal, 
        Tinput = Tinput, 
        tinterval = tinterval, 
        expSTOP = expSTOP, 
        
        # Subject-relevant information
        expAGE = expAGE, 
        expBW = expBW, 
        sex = sex, 
        
        # Simulation relevant information
        Tstart = Tstart, 
        Tstop = Tstop, 
        Dt = Dt 
        
        )
      
    } else { # Individual, no lifestage
      
      RUNandOUT(
        
        # Exposure-relevant information
        exposure_type = exposure_type, 
        exp = exp, 
        exp_Oral = exp_Oral, 
        exp_Dermal = exp_Dermal, 
        Tinput = Tinput, 
        tinterval = tinterval, 
        expSTOP = expSTOP, 
        
        # Subject-relevant information
        expAGE = expAGE, 
        expBW = expBW, 
        sex = sex, 
        
        # Simulation relevant information
        Tstart = Tstart, 
        Tstop = Tstop, 
        Dt = Dt 
        
      )
      
    }
    
  }
  
  
  # Save Results ----

  if(Population == "Yes"){
    RESULTS <- Pop.RESULTS
    } else {
    RESULTS <- MODEL_OUTPUT
    }
  
  save(RESULTS, file = here(OUTPUT, "RESULTS_PFOA_PBK.RData"))
  
 