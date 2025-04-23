# --------------------------------------------------------------------------- #
# SCRIPT FOR RUNNING THE PBK MODEL & OBTAINING RESULTS FOR A POPULATION
# By: Chrysanthi Pachoulide
# Date: 23-04-2025
# --------------------------------------------------------------------------- #

rm(list=ls()) 

# Packages
library(here)
library(tidyverse)
library(deSolve)
library(PKNCA)
library(pracma)


# Load PBK files
Physio.c <- read_csv(here("Input", "PhysioVariables.csv"))
Tissue.c <- read_csv(here("Input", "TissueComposition.csv"))
source(here("Script", "RUN_and_OUTPUT.R"))


# Load population 
INPUT_dummy <- read.csv(here("Input", "INPUT_dummy.csv")) 
nPeople <- as.numeric(nrow(INPUT_dummy)) # number of people


# List for storing results
Pop.RESULTS <- list(
  CALC_Parameters = vector("list", nPeople), # list of length nPeople
  OUT_RAW_data = vector("list", nPeople),
  ANALYSED_data = vector("list", nPeople),
  OUT_Plots = vector("list", nPeople) # could be removed if it's too heavy for R
)


# Set output storage directory
OUTPUT <- here("Output", format(Sys.Date(), "%Y-%m-%d"), format(Sys.time(), "%H-%M-%S"))
dir.create(OUTPUT, recursive = TRUE)


# Run the model for the population
for (i in 1:nPeople) {

  # Run the model per person
  RESULTS <- RUNandOUT(
    
    exposure_type <- as.character(INPUT_dummy[i, "exposure_type"]),            # Exposure type, choose between "Oral", "Dermal", "Oral_Dermal" (if exposure is both Oral and Dermal), Inhalation"
    expCONC = as.numeric(INPUT_dummy[i, "expCONC"]),                           # ug/kg/day concentration (total exposure concentration)
    expCONC_Oral <- ifelse(is.na(INPUT_dummy[i, "expCONC_Oral"]), 0,
                           as.numeric(INPUT_dummy[i, "expCONC_Oral"])),        # ug/kg/day concentration to be used only when both oral and dermal are used
    expCONC_Dermal <- ifelse(is.na(INPUT_dummy[i, "expCONC_Dermal"]), 0,
                             as.numeric(INPUT_dummy[i, "expCONC_Dermal"])),    # ug/kg/day concentration to be used only when both oral and dermal are used
    Tinput <- ifelse(is.na(INPUT_dummy[i, "Tinput"]), 1,
                     as.numeric(INPUT_dummy[i, "Tinput"])),          # for repeated exposure or so default = 1
    tinterval <- ifelse(is.na(INPUT_dummy[i, "tinterval"]), 1,
                        as.numeric(INPUT_dummy[i, "tinterval"])),    # for repeated exposure or so default = 1
    expSTOP <- as.numeric(INPUT_dummy[i, "expSTOP"]),                # time in days after which the exposure stopped
    
    # Subject-relevant information
    expAGE <- ifelse(is.na(INPUT_dummy[i, "expAGE"]), NA,
                     as.numeric(INPUT_dummy[i, "expAGE"])),     # years old age at exposure if not provided then age argument is not used physiology is based on BW
    expBW <- ifelse(is.na(INPUT_dummy[i, "expBW"]), 70,
                    as.numeric(INPUT_dummy[i, "expBW"])),       # kg if not provided then the BW of the corresponding age and sex is taken; if both BW and Age are not given then a default BW = 70 is taken; if BW is higher than the BW from the lifestage equations then the actual BW overwrites the calculated one
    sex <- ifelse(is.na(INPUT_dummy[i, "sex"]), "M",
                  as.character(INPUT_dummy[i, "sex"])),         # sex either "F" or "M" if none then default is "M"
    
    # Simulation relevant information
    Tstart <- 0,          # days, start of the simulation
    Tstop <- 50*456,      # days, stop of the simulation
    Dt <- 1               # days, iteration steps (decrease/increase depending on run time)
    
    )

  # Collect population results together
  Pop.RESULTS$CALC_Parameters[[i]] <- RESULTS$CALC_Parameters
  Pop.RESULTS$OUT_RAW_data[[i]] <- RESULTS$OUT_RAW_data
  Pop.RESULTS$ANALYSED_data[[i]] <- RESULTS$ANALYSED_data
  Pop.RESULTS$OUT_Plots[[i]] <- RESULTS$OUT_Plots # could be removed if it's too heavy for R


  }







