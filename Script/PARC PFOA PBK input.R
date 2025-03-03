# --------------------------------------------------------------------------- #
# PBK MODEL FOR PFOA, TO BE USED TOGETHER WITH THE LATEST HBM DATA
# Input file
# CP, 10-11-2024
# --------------------------------------------------------------------------- #

rm(list=ls()) # to clear out the global environment


# Set working directory

HOME = "C:/Users/pacho003/OneDrive - Wageningen University & Research/C Channel/R/PARC_PFAS_PBPKmodel/Codes_PFAS_models/Chrysanthi_Pachoulide_2024_PFOA_oral_dermal_blood"
# HOME = "/home/westerj"
setwd(HOME)

# Set input storage directory
INPUT = file.path("Input", Sys.Date())
dir.create(INPUT, recursive = TRUE)
setwd(INPUT)


# Load packages

library(lubridate)
library(ggplot2)
library(deSolve)
library(writexl)
library(ggpubr)
library(tidyverse)


# For plotting
library(showtext)
font_add(family = "Garamond", regular = "GARA.TTF")
showtext_auto()


# PFOA SPECIFIC KINETIC PARAMETERS ####
# ------------------------------------------------------ #

## Calculate partition coefficients based on tissue matrix distribution experimental data ####
Physio.data <- read_csv("C:/Users/pacho003/OneDrive - Wageningen University & Research/C Channel/R/PARC_PFAS_PBPKmodel/Codes_PFAS_models/Chrysanthi_Pachoulide_2024_PFOA_oral_dermal_blood/Input/TissueComposition.csv")

# PFOA Matrix/Water Distribution Coefficients

logML <- 3.52  # 3.52 ± 0.08 from Ebert A., Allendorf F. et al (2020), https://dx.doi.org/10.1021/acs.est.0c00175  (Liposomes composed of POPC (1-palmitoyl-2-oleoyl-glycero-3-phosphocholine))
logSP <- 1.61 # 1.61 ± 0.15 from Allendorf, F., Goss, K.-U. and Ulrich, N. (2021), https://doi.org/10.1002/etc.4954  (Structural proteins from chicken breast fillet (actin & myosin 60-95%), Recovery 95%)
logALB <- 4.33 # 4.33 ± 0.05 from Allendorf, F., Goss, K.-U. and Ulrich, N. (2021), https://doi.org/10.1002/etc.4954 (BSA (fatty acid free) Molar ratio compound to BSA < 0.1 72-96h, Recovery 94%))
logSL <- -1.37 # -1.37 ± 0.01 from Allendorf, F., Goss, K.-U. and Ulrich, N. (2021), https://doi.org/10.1002/etc.4954 (Olive oil with a high fraction of unsaturated fatty acids, Recovery 95%)
logFABP <- 4.3 # calculated by Allendorf, F., Goss, K.-U. and Ulrich, N. (2021), https://doi.org/10.1002/etc.4954

k_ML <- 10^logML       # neutral phospholipids:water partition coefficient
k_SP <- 10^logSP       # structural protein:water partition coefficient
k_ALB <- 10^logALB     # albumin:water partition coefficient
k_SL <- 10^logSL       # structural lipids:water partition coefficient
k_FABP <- 10^logFABP   # fatty acid-binding protein:water partition coefficient


Physio.dat <- Physio.data %>% filter(!Tissue %in% c("Comment", "NamingInUtsey")) %>% 
  select(-Comment) %>% 
  mutate(across(-Tissue, as.numeric))

Physio.Tissues <- Physio.dat %>% filter(!Tissue %in% c("Blood", "Plasma")) 

Kp <- (Physio.Tissues$f_W +
         (k_ML*Physio.Tissues$f_ML) +
         (k_SL*Physio.Tissues$f_SL) +
         (k_SP*Physio.Tissues$f_SP) +
         (k_ALB*Physio.Tissues$f_ALB) +
         (k_FABP*Physio.Tissues$f_FABP)
) #* fup Removing fup from here, Kps will be recalculated directly in the model to enable sensitivity analysis

Names <- Physio.Tissues$Tissue %>% substr(1,2)
Names <- paste("Kp",Names ,sep="")
Kp.df <- as.data.frame(list(Tissue = Names)) %>%
  mutate(Value = Kp) %>% 
  pivot_wider(names_from = Tissue, values_from = Value)


PFOA.params_df <- tibble(
  MW = 414.07, # PFOA molecular weight g/mol
  fup = 0.061/100, # fFraction unbound in plasma, unitless (fraction) Fischer et al. 2024 https://doi.org/10.1021/acs.est.3c07415; OLD fup = 0.02 # fup fraction of PFOA in plasma
  
  # These are the values from Trine's model
  # Comment Chrysa 10-11-2024: the initial partition coefficients already accound for the fraction unbound as we assume the tissues/plasma were in equilibrium in the animal's body and only the fraction ubound was able to move between plasma and tissues
  # Comment Chrysa 10-11-2024: dividing with fraction unbound to be able to see the influence of the fraction unbound in the sensitivity analysis
  PF = 0.04/fup,  # Plasma/fat partition coefficient; Rat tissue data (Kudo et al. 2007
  PG = 0.05/fup,  # Plasma/gut partition coefficient; Rat tissue data (Kudo et al. 2007 %>%
  PK = 1.05/fup,  # Plasma/kidney partition coefficient; Rat tissue data (Kudo et al. 2007
  PL = 2.2/fup,  # Plasma/liver partition coefficient; Rat tissue data (Kudo et al. 2007
  PSk = 0.1/fup,  # Plasma/skin partition coefficient; Rat tissue data (Kudo et al. 2007
  PR = 0.12/fup  # Plasma/rest of the body partition coefficient; Rat tissue data (Kudo et al. 2007
  )

# Adding the calculated partition coefficients
PFOA.params_df <- cbind(PFOA.params_df, Kp.df)

# Adding the rest of parameters
PFOA.params_df <- PFOA.params_df %>% mutate(
  # Dermal absorption
  
  # Comment Chrysa 21-10-2024: AbsPFOA is used in the Dermaldose input so we're already correcting before using the Papp?! not sure I agree with this
  fBAc_dust = 69.6/100, # 69.6±5 # Fraction bio accessible from dust; Ragnarsdottir et al.: https://doi.org/10.1016/j.envres.2023.117093; fraction bioaccessible = fraction of the compound released from the matrix (cosmetic formulation dust etc and is available to be absorbed from the epidermis
  fBAc_cosm = 97.8/100, # 97.8±19.0 !! this is an estimated value, based on [PFOA]_cosmetics and the LOQ # Fraction bio accessible from cosmetics; Namazkar et al.: DOI: 10.1039/d3em00461a
  AbsPFOA = 0.016, # from Trine; 0.00048 # Changed to the absorption measured by Abraham and Monien 2022 of 1.6% of applied dose from sunscreen.
  Papp = 3.82 * 10^-3, # cm/h ref: https://doi.org/10.1016/j.envint.2024.108772
  
  # Kidney clearance
  
  # Comment Chrysa 05-11-2024:
  ## Assuming that the kidney PFAS concentrations never reaches Km concentrations therefore transforming Vmax and Km to a Clearance; in the paper what they call transporter efficiency : Louisse Pedroni et al. 2024 https://doi.org/10.1016/j.tox.2024.153961
  ## OAT1 and OAT3 are determining transport between blood and proximal tubule
  ## OAT4 is determining transport between proximal tubule and filtrate
  ## All transporters are bi-directional therefore the equation is written assuming that direction is determined based on the equilibrium between the two compartments
  ## Introducing the "affinity constants: kAap and kAbl to compensate for the fact that the transporters have an affinity to one side.
  kAbl = 0.01, # affinity constant basolateral this is about OAT1 and OAT3 which have affinity to uptake (movement from plasma to cells; this is fitted value for now; kAbl = 0.01 is driving the equilibrium towards uptake into the proximal tubule cells
  kAap = 0.01, # affinity constant apical this is about OAT4 which has affinity to re-abs (movement from filtrate to cells; this is fitted value for now; kAap = 0.01 is driving the equilibrium towards re-absorption into the proximal tubule cells
  
  # Comment Chrysa 12-11-2024: I'm not sure about these scaling factors as the in vitro clearance is expressed in ul/min/mg protein and not ul/min/HEK cells
  # Comment Chrysa 16-01-2025: conversion mistake in the below equations, I'm correcting from / 10^-6 *60*24 * 10^-6 to / 10^-6 *60*24 * 10^6
  # Comment Chrysa 23-01-2024: Conversion mistake again.. I'm correcting from: ex. CL_OAT1 = 19/ 10^-6 *60*24 * 10^6 to  CL_OAT1 = 19*10^-6*60*24*10^6
  CL_OAT1 = 19 * 10^-6 *60*24 * 10^6, # L/d/kg protein; initial ul/min/mg protein -> * 10^-6 L/min/mg protein, *60*24 L/d/mg protein, *10^6 L/d/Kg protein
  CL_OAT3 = 17 * 10^-6 *60*24 * 10^6, # L/d/kg protein; initial ul/min/mg protein -> * 10^-6 L/min/mg protein, *60*24 L/d/mg protein, *10^6 L/d/Kg protein
  CL_OAT4 = 96 * 10^-6 *60*24 * 10^6, # L/d/kg protein; initial ul/min/mg protein -> * 10^-6 L/min/mg protein, *60*24 L/d/mg protein, *10^6 L/d/Kg protein
  
  PTCPGK = 9.94* 10^7 * 10^3, # proximal tubule cells/kg kidney cortex; initial 99.4 million PTC/g kidney https://doi.org/10.1021/acs.molpharmaceut.4c00504
  
  InVivo_OAT1 = 5.3, # 5.3 ± 1.9 OAT1 /mg membrane protein in the human kidney cortex http://dx.doi.org/10.1124/dmd.116.072066; alternative : 4.3 ± 0.3 pmol OAT1 /mg membrane protein in the human kidney cortex https://doi.org/10.1124/dmd.121.000367;
  InVitro_OAT1 = 26.6, # 26.6 ± 3.4 pmol/mg membrane protein: OAT1 expression in HEK293-OAT1 cells https://doi.org/10.1124/dmd.121.000367;
  InVivo_OAT3 = 3.5, # 3.5 ± 1.6 OAT3/mg membrane protein in the human kidney cortex http://dx.doi.org/10.1124/dmd.116.072066; alternative : 2.7 ± 0.1 pmol OAT3 /mg membrane protein in the human kidney cortex https://doi.org/10.1124/dmd.121.000367;
  InVitro_OAT3 = 7.3, # 7.3 ± 0.5 pmol/mg membrane protein: OAT3 expression in HEK293-OAT1 cells https://doi.org/10.1124/dmd.121.000367;
  InVivo_OAT4 = 0.52, # 0.52 ± 0.23 pmol OAT4 /mg membrane protein in the human kidney cortex http://dx.doi.org/10.1124/dmd.116.072066;
  InVitro_OAT4 = (InVitro_OAT1 + InVitro_OAT3)/2, # OAT4 expression in HEK293-OAT4 cells not found so assuming average of OAT1 and OAT3 expression
  
  REF_OAT1 = InVivo_OAT1/InVitro_OAT1, # relative expression factor OAT1: expression in the human kidneys /expression in the cells
  REF_OAT3 = InVivo_OAT3/InVitro_OAT3, # relative expression factor OAT3: expression in the human kidneys /expression in the cells
  REF_OAT4 = 1, # assuming 1 as Invitro expression not found InVivo_OAT4/InVitro_OAT4 # relative expression factor OAT4: expression in the human kidneys /expression in the cells
  
  CLurinec = 0.000044,  # L/d/kg; 0.044 mL/d/kg taken from Fujii et al 2015
  VmaxOAT4c = 4.5*MW/1000*60*24, # ug/d/mg protein; 45 nmol/min/mg protein *MW/1000*60*24 = ug/d/mg protein ref: Louisse et al. 2023 https://doi.org/10.1007/s00204-022-03428-6
  KmOAT4c = 47*MW, # ug/L; 47 uM*MW = ug/L ref: Louisse et al. 2023
  
  
  # Hepatic clearance
  # These are the values from Trine's model
  CLbiliaryc = 0.00262, # L/d/kg ; 2.62 +/- 3.6 mL/d/kg from Fujii et al 2015 DOI: 10.1539/joh.14-0136-OA
  CLfaecesc = 0.000052 # L/d/kg ; 0.052 +/- 0.05 mL/d/kg clearance in faeces taken from Fujii et al 2015 DOI: 10.1539/joh.14-0136-OA
  )

write.csv(PFOA.params_df, "PFOAParams.csv", row.names = FALSE)


# LIFETIME EQUATIONS ####
# ------------------------------------------------------ #

lifeTSTOP = 80 # duration of lifetime (0 - 80 years old)  of simulation
TSTART = 0
TSTOP = 365*lifeTSTOP # years in days
DT = 1
TIME = seq(TSTART,TSTOP,by=DT)


# Creating a dataframe for all the variables
Variables_df = as.data.frame(list(TIME = TIME)) #df column 1 = simulation time, every step is 1 day
Variables_df = Variables_df %>%
  mutate(age = TIME/365) # add column 2 = age in days

Variables_df = Variables_df %>%
  
## Fractional Volumes ####
# Body weight
# BW_M_Ratier_2024 & BW_F_Ratier_2024 = Equation extracted from supplemental material from Ratier et al., 2024
mutate(BW_M_Ratier_2024 = if_else(age <19.00093277, 74.16235828-(2*(74.16235828-57.19957758)/(exp(0.63466182*(age-13.31018000))+exp(0.05457656*(age-13.31018000)))),
                                  -0.01129273*age^2 + 1.11817056*age + 56.74397436)) %>%
  mutate(BW_F_Ratier_2024 = if_else(age <17.9374115, 62.95490567-(2*(62.95490567-49.36574299)/(exp(0.84039606*(age-11.56691488))+exp(0.06710088*(age-11.56691488)))),
                                    -0.01258006*age^2 + 1.25029379*age + 44.4459234)) %>%
  mutate(BDW_M_Ratier_2024 = 74.16235828-(2*(74.16235828-57.19957758)/(exp(0.63466182*(age-13.31018000))+exp(0.05457656*(age-13.31018000))))) %>%
  mutate(BDW_F_Ratier_2024 = 62.95490567-(2*(62.95490567-49.36574299)/(exp(0.84039606*(age-11.56691488))+exp(0.06710088*(age-11.56691488)))))

# Blood/Plasma/Hematocrit
# Using Ratier et al. (2024) model, Fraction of arterial plasma, calculated from Filser 2000 p.43
Fr_art_plasma = 0.0178 / (0.0178 + 0.0533) #fraction of arterial blood (corrected for plasma)

# Hematocrit - male                                                              # From Supp mat of Brochot et al. 2019
Param1_M = 33.455469
Param2_M = 53.206039
Param3_M = 8.277945
Param4_M = 40.492556
Param5_M = 46.899695

b1_M = (Param4_M - Param1_M -(Param2_M - Param1_M) * exp(-Param3_M))/5
a1_M = Param4_M - 6*b1_M
b2_M = (Param5_M - Param4_M)/5
a2_M = Param4_M - 15*b2_M

# Hematocrit - non pragnant female
Param1_F = 32.617402
Param2_F = 53.188459
Param3_F = 7.699418
Param4_F = 37.531463
Param5_F = 40.055284

b1_F = (Param4_F - Param1_F - (Param2_F-Param1_F)*exp(-Param3_F))/2
a1_F = Param4_F - 3*b1_F
b2_F = (Param5_F - Param4_F)/7
a2_F = Param5_F - 10*b2_F

Variables_df = Variables_df %>%
  select(TIME,age,BW_M_Ratier_2024,BW_F_Ratier_2024,BDW_M_Ratier_2024,BDW_F_Ratier_2024) %>%
  rename(BW_M = BW_M_Ratier_2024) %>% #could actually be ignored as we only use BDW and not BW
  rename(BW_F = BW_F_Ratier_2024) %>% #could actually be ignored as we only use BDW and not BW
  rename(BDW_M = BDW_M_Ratier_2024) %>%
  rename(BDW_F = BDW_F_Ratier_2024) %>%
  
  # Adrenal; compartment [1] in Ratier 2024 (not used in our model, but needed for calculation of adipose tissue)
  mutate(V_adrenalFraction_M = 0.0002 + (0.00171 - 0.0002)*exp(-2.02*age)) %>%
  mutate(V_adrenalFraction_F = 0.0002 + (0.00171 - 0.0002)*exp(-2.02*age)) %>%
  
  # Bone; compartment [2] in Ratier 2024 (not used in our model, but needed for calculation of adipose tissue)
  mutate(V_boneFraction_M = (0.313 + (0.506 - 0.313)*exp(-0.0907*age))*0.095) %>%
  mutate(V_bonenonperfusedFraction_M = 0.095 - V_boneFraction_M) %>%
  mutate(V_boneFraction_F = (0.298 + (0.505 - 0.298)*exp(-0.0792*age))*0.085) %>%
  mutate(V_bonenonperfusedFraction_F = 0.085 - V_boneFraction_F) %>%
  
  # Brain; compartment [3] in Ratier 2024 (not used in our model, but needed for calculation of adipose tissue)
  mutate(V_brainFraction_M = (1.450 + (0.353 - 1.450) * exp (-0.440*age))/BDW_M) %>%
  mutate(V_brainFraction_F = (1.300 + (0.347 - 1.300) * exp (-0.573*age))/BDW_F) %>%
  
  # Breast; compartment [4] in Ratier 2024 (not used in our model, but needed for calculation of adipose tissue)
  mutate(V_breastFraction_M = 3.42E-4*1/(1 + exp(-1.42*age + 20.1))) %>%
  mutate(V_breastFraction_F = 0.00833/(1 + exp(-1.92*age+ 28.6))) %>%
  
  # Heart; compartment [5] in Ratier 2024 (not used in our model, but needed for calculation of adipose tissue)
  mutate(V_heartFraction_M = 0.0045) %>%
  mutate(V_heartFraction_F = 0.0042) %>%
  
  # Marrow; compartment [6] in Ratier 2024 (not used in our model, but needed for calculation of adipose tissue)
  mutate(V_marrowFraction_M = 0.05 + (0.0138 - 0.05)*exp(-0.112*age)) %>%
  mutate(V_marrowFraction_F = 0.045 + (0.0138 - 0.045)*exp(-0.136*age)) %>%
  
  # Muscle; compartment [7] in Ratier 2024 (not used in our model, but needed for calculation of adipose tissue)
  mutate(MuscleAtrophy_M = if_else(age < 24.3, 1,
                                   (-0.0001264*age^2 + 0.006131*age + 0.926))) %>%
  mutate(V_muscleFraction_M = (0.3973 + (0.201 - 0.3973)*exp(-0.141*age)) * MuscleAtrophy_M) %>%
  mutate(MuscleAtrophy_F = if_else(age < 25.90709, 1,
                                   (-0.0001264*age^2 + 0.006131*age + 0.926))) %>%
  mutate(V_muscleFraction_F = (0.2917 + (0.207 - 0.2917)*exp(-0.339*age)) * MuscleAtrophy_F) %>%
  
  # Reproductive organs; compartment [8] in Ratier 2024 (not used in our model, but needed for calculation of adipose tissue)
  mutate(V_reproFraction_M = if_else(age < 20.01, -1.5156E-07*age^3 + 9.3351E-06*age^2 - 1.1177E-04*age + 4.7966E-04,
                                     0.0008)) %>%
  mutate(V_reproFraction_F = if_else(age < 1, -1.064E-3*age + 1.338E-3,
                                     if_else(age < 20, 2.6380E-7*age^3 - 1.7943E-6*age^2 - 5.6465E-6*age + 2.8105E-4,
                                             0.001552))) %>%
  
  # Pancreas; compartment [9] in Ratier 2024 (not used in our model, but needed for calculation of adipose tissue)
  mutate(V_pancreasFraction_M = 0.00192) %>%
  mutate(V_pancreasFraction_F = 0.002) %>%
  
  # Skin; compartment [10] in Ratier 2024
  mutate(V_skinFraction_M = if_else(age < 20.01, -1.1706E-05*age^3 + 5.4130E-04*age^2 - 6.1966E-03*age + 4.6231E-02,
                                    0.0452)) %>%
  mutate(V_skinFraction_F = if_else(age < 19.45, -7.8882E-06*age^3 + 4.0224E-04*age^2 - 5.2146E-03*age + 4.5605E-02,
                                    0.0383)) %>%
  
  # Spleen; compartment [11] in Ratier 2024 (not used in our model, but needed for calculation of adipose tissue)
  mutate(V_spleenFraction_M = 0.0021) %>%
  mutate(V_spleenFraction_F = 0.0022) %>%
  
  # Thyroid; compartment [12] in Ratier 2024 (not used in our model, but needed for calculation of adipose tissue)
  mutate(V_thyroidFraction_M = 0.000274) %>%
  mutate(V_thyroidFraction_F = 0.0003) %>%
  
  # Urinary tract (bladder, ureters, urethra); compartment [13] in Ratier 2024
  mutate(V_urinarytractFraction_M = 0.00104) %>%
  mutate(V_urinarytractFraction_F = 0.0010) %>%
  
  # Kidney; compartment [14] in Ratier 2024
  mutate(V_kidneyFraction_M = 0.0042 + (0.00767 - 0.0042)*exp(-0.206*age)) %>%
  mutate(V_kidneyFraction_F = 0.0046 + (0.0071 - 0.0046)*exp(-0.221*age)) %>%
  
  # Lungs; compartment [15] in Ratier 2024 (not used in our model, but needed for calculation of adipose tissue)
  # Comment Chrysa on 18-10-2024: Shouldn't the lungs take 100% of the blood flow?
  mutate(V_lungFraction_M = 0.0068) %>%
  mutate(V_lungFraction_F = 0.0070) %>%
  
  # Gut; compartment [16] in Ratier 2024
  mutate(V_gutFraction_M = if_else(age < 16, -0.000082562*age^2 + 0.0013523*age + 0.01293,
                                   0.0140)) %>%
  mutate(V_gutFraction_F = if_else(age < 14.453301, -7.42E-5*age^2 + 1.28E-3*age + 1.30E-2,
                                   0.0160)) %>%
  
  # Stomach; compartment [17] in Ratier 2024
  mutate(V_stomachFraction_M = 0.0021) %>%
  mutate(V_stomachFraction_F = 0.0023) %>%
  
  # Liver; compartments [18] (liver) and [21] (liver artery) in Ratier 2024
  mutate(V_liverFraction_M = 0.0247 + (0.0409 - 0.0247)*exp(-0.218*age)) %>%
  mutate(V_liverFraction_F = 0.0233 + (0.038 - 0.0233)*exp(-0.122*age)) %>%
  
  # Plasma volume; compartment [22] in Ratier 2024 -> USED TO BE BLOOD volume, as it's corrected for hematocrit then it's plasma
  mutate(V_plasmaFraction_M = if_else(age < 1, (-0.0273*age + 0.0771),
                                      0.0761 + (0.0289 - 0.0761)*exp(-0.592*age))) %>%
  mutate(V_plasmaFraction_F = if_else(age < 1, (-0.0273*age + 0.0771),
                                      if_else(age < 14.019723, 3.28E-5*age^3 - 1.21E-3*age^2 + 1.24E-2*age + 3.86E-2,
                                              0.065))) %>%
  
  # Adipose tissue
  mutate(V_adiposeFraction_M = 0.96 - V_adrenalFraction_M - V_boneFraction_M - V_bonenonperfusedFraction_M - V_brainFraction_M - V_breastFraction_M +
           - V_heartFraction_M - V_marrowFraction_M - V_muscleFraction_M - V_reproFraction_M - V_pancreasFraction_M +
           - V_skinFraction_M - V_spleenFraction_M - V_thyroidFraction_M - V_urinarytractFraction_M - V_kidneyFraction_M +
           - V_lungFraction_M - V_gutFraction_M - V_stomachFraction_M - V_liverFraction_M - V_plasmaFraction_M) %>%
  mutate(AdiposeMass_M = if_else(age < 19.00093277, 0,
                                 (-0.01129273*age^2 + 1.11817056*age + 56.74397436) - (74.16235828-(2*(74.16235828-57.19957758)/(exp(0.63466182*(age-13.31018000))+exp(0.05457656*(age-13.31018000))))))) %>% # age and not BW dependent
  mutate(V_adiposeFraction_F = 0.96 - V_adrenalFraction_F - V_boneFraction_F - V_bonenonperfusedFraction_F - V_brainFraction_F - V_breastFraction_F +
           - V_heartFraction_F - V_marrowFraction_F - V_muscleFraction_F - V_reproFraction_F - V_pancreasFraction_F +
           - V_skinFraction_F - V_spleenFraction_F - V_thyroidFraction_F - V_urinarytractFraction_F - V_kidneyFraction_F +
           - V_lungFraction_F - V_gutFraction_F - V_stomachFraction_F - V_liverFraction_F - V_plasmaFraction_F) %>%
  mutate(AdiposeMass_F = if_else(age < 17.9374115, 0,
                                 if_else(((-0.01258006*age^2 + 1.25029379*age + 44.4459234) - (62.95490567-(2*(62.95490567-49.36574299)/(exp(0.84039606*(age-11.56691488))+exp(0.06710088*(age-11.56691488)))))) < 0, 0,
                                         (-0.01258006*age^2 + 1.25029379*age + 44.4459234) - (62.95490567-(2*(62.95490567-49.36574299)/(exp(0.84039606*(age-11.56691488))+exp(0.06710088*(age-11.56691488)))))))) # what does this mean?!!


## Fractional Blood Flows ####

Variables_df = Variables_df %>%
  
  # Hematocrit
  mutate(Hct_ven_M = if_else(age < 1, (Param1_M +(Param2_M-Param1_M)*exp(-Param3_M*age))*0.01,
                             if_else(age < 6, (a1_M + b1_M*age)*0.01,
                                     if_else(age < 15, Param4_M*0.01,
                                             if_else(age < 20, (a2_M + b2_M*age)*0.01,
                                                     Param5_M*0.01))))) %>%
  mutate(Hct_M = Hct_ven_M*0.91) %>%
  mutate(Hct_ven_F = if_else(age < 1, (Param1_F +(Param2_F-Param1_F)*exp(-Param3_F*age))*0.01,
                             if_else(age < 3, (a1_F + b1_F*age)*0.01,
                                     if_else(age < 10, (a2_F + b2_F*age)*0.01,
                                             Param5_F*0.01)))) %>%
  mutate(Hct_F = Hct_ven_F*0.91) %>%
  
  # Cardiac output (plasma; L/min*60*24 = L/d)
  mutate(CardOut_M = if_else(age < 33.37, (6.642 + (0.6 - 6.642)*exp(-0.1323*age))*(1-Hct_M)*60*24,
                             (-0.000895*age^2 + 0.0607*age + 5.54)*(1-Hct_M)*60*24)) %>%
  mutate(CardOut_F = if_else(age < 16.027, (7.734 + (0.6 - 7.734)*exp(-0.09747*age))*(1-Hct_F)*60*24,
                             (0.000473*age^2 - 0.0782*age + 7.37)*(1-Hct_F)*60*24)) %>%
  
  # Adrenal; compartment [1] in Ratier 2024 (not used in our model, but needed for calculation of adipose tissue)
  mutate(Q_adrenalFraction_M = (V_adrenalFraction_M/0.0002)*0.003) %>% # sc_F[0-17] = (sc_V[i]  / sc_V_adult[i])  * sc_F_adult[i];
  mutate(Q_adrenalFraction_F = (V_adrenalFraction_F/0.0002)*0.003) %>% # sc_F[0-17] = (sc_V[i]  / sc_V_adult[i])  * sc_F_adult[i];
  
  # Bone; compartment [2] in Ratier 2024 (not used in our model, but needed for calculation of adipose tissue)
  mutate(Q_boneFraction_M = (V_boneFraction_M/(0.095*0.32))*0.021) %>% # sc_F[0-17] = (sc_V[i]  / sc_V_adult[i])  * sc_F_adult[i];
  mutate(Q_boneFraction_F = (V_boneFraction_F/(0.085*0.298))*0.021) %>% # sc_F[0-17] = (sc_V[i]  / sc_V_adult[i])  * sc_F_adult[i];
  
  # Brain; compartment [3] in Ratier 2024 (not used in our model, but needed for calculation of adipose tissue)
  mutate(Q_brainFraction_M = (V_brainFraction_M/0.01986)*0.124) %>% # sc_F[0-17] = (sc_V[i]  / sc_V_adult[i])  * sc_F_adult[i];
  mutate(Q_brainFraction_F = (V_brainFraction_F/0.0217)*0.124) %>% # sc_F[0-17] = (sc_V[i]  / sc_V_adult[i])  * sc_F_adult[i];
  
  # Breast; compartment [4] in Ratier 2024 (not used in our model, but needed for calculation of adipose tissue)
  mutate(Q_breastFraction_M = (V_breastFraction_M/0.00035)*0.0002) %>% # sc_F[0-17] = (sc_V[i]  / sc_V_adult[i])  * sc_F_adult[i];
  mutate(Q_breastFraction_F = (V_breastFraction_F/0.0083)*0.004) %>% # sc_F[0-17] = (sc_V[i]  / sc_V_adult[i])  * sc_F_adult[i];
  
  # Heart; compartment [5] in Ratier 2024 (not used in our model, but needed for calculation of adipose tissue)
  mutate(Q_heartFraction_M = (V_heartFraction_M/0.0045)*0.041) %>% # sc_F[0-17] = (sc_V[i]  / sc_V_adult[i])  * sc_F_adult[i];
  mutate(Q_heartFraction_F = (V_heartFraction_F/0.0042)*0.051) %>% # sc_F[0-17] = (sc_V[i]  / sc_V_adult[i])  * sc_F_adult[i];
  
  # Marrow; compartment [6] in Ratier 2024 (not used in our model, but needed for calculation of adipose tissue)
  mutate(Q_marrowFraction_M = (V_marrowFraction_M/0.050)*0.031) %>% # sc_F[0-17] = (sc_V[i]  / sc_V_adult[i])  * sc_F_adult[i];
  mutate(Q_marrowFraction_F = (V_marrowFraction_F/0.045)*0.031) %>% # sc_F[0-17] = (sc_V[i]  / sc_V_adult[i])  * sc_F_adult[i];
  
  # Muscle; compartment [7] in Ratier 2024 (not used in our model, but needed for calculation of adipose tissue)
  mutate(Q_muscleFraction_M = (V_muscleFraction_M/0.3973)*0.175) %>% # sc_F[0-17] = (sc_V[i]  / sc_V_adult[i])  * sc_F_adult[i];
  mutate(Q_muscleFraction_F = (V_muscleFraction_F/0.2917)*0.124) %>% # sc_F[0-17] = (sc_V[i]  / sc_V_adult[i])  * sc_F_adult[i];
  
  # Reproductive organs; compartment [8] in Ratier 2024 (not used in our model, but needed for calculation of adipose tissue)
  mutate(Q_reproFraction_M = (V_reproFraction_M/0.0008)*0.001) %>% # sc_F[0-17] = (sc_V[i]  / sc_V_adult[i])  * sc_F_adult[i];
  mutate(Q_reproFraction_F = (V_reproFraction_F/0.0016)*0.004) %>% # sc_F[0-17] = (sc_V[i]  / sc_V_adult[i])  * sc_F_adult[i];
  
  # Pancreas; compartment [9] in Ratier 2024 (not used in our model, but needed for calculation of adipose tissue)
  mutate(Q_pancreasFraction_M = (V_pancreasFraction_M/0.00192)*0.01) %>% # sc_F[0-17] = (sc_V[i]  / sc_V_adult[i])  * sc_F_adult[i];
  mutate(Q_pancreasFraction_F = (V_pancreasFraction_F/0.002)*0.01) %>% # sc_F[0-17] = (sc_V[i]  / sc_V_adult[i])  * sc_F_adult[i];
  
  # Skin; compartment [10] in Ratier 2024
  mutate(Q_skinFraction_M = (V_skinFraction_M/0.0452)*0.052) %>% # sc_F[0-17] = (sc_V[i]  / sc_V_adult[i])  * sc_F_adult[i];
  mutate(Q_skinFraction_F = (V_skinFraction_F/0.0383)*0.051) %>% # sc_F[0-17] = (sc_V[i]  / sc_V_adult[i])  * sc_F_adult[i];
  
  # Spleen; compartment [11] in Ratier 2024 (not used in our model, but needed for calculation of adipose tissue)
  mutate(Q_spleenFraction_M = (V_spleenFraction_M/0.0021)*0.031) %>% # sc_F[0-17] = (sc_V[i]  / sc_V_adult[i])  * sc_F_adult[i];
  mutate(Q_spleenFraction_F = (V_spleenFraction_F/0.0022)*0.031) %>% # sc_F[0-17] = (sc_V[i]  / sc_V_adult[i])  * sc_F_adult[i];
  
  # Thyroid; compartment [12] in Ratier 2024 (not used in our model, but needed for calculation of adipose tissue)
  mutate(Q_thyroidFraction_M = (V_thyroidFraction_M/0.000274)*0.015) %>% # sc_F[0-17] = (sc_V[i]  / sc_V_adult[i])  * sc_F_adult[i];
  mutate(Q_thyroidFraction_F = (V_thyroidFraction_F/0.0003)*0.015) %>% # sc_F[0-17] = (sc_V[i]  / sc_V_adult[i])  * sc_F_adult[i];
  
  # Urinary tract (bladder, ureters, urethra); compartment [13] in Ratier 2024
  mutate(Q_urinarytractFraction_M = (V_urinarytractFraction_M/0.00104)*0.001) %>% # sc_F[0-17] = (sc_V[i]  / sc_V_adult[i])  * sc_F_adult[i];
  mutate(Q_urinarytractFraction_F = (V_urinarytractFraction_F/0.0010)*0.001) %>% # sc_F[0-17] = (sc_V[i]  / sc_V_adult[i])  * sc_F_adult[i];
  
  # Kidney; compartment [14] in Ratier 2024
  mutate(Q_kidneyFraction_M = (V_kidneyFraction_M/0.0042)*0.196) %>% # sc_F[0-17] = (sc_V[i]  / sc_V_adult[i])  * sc_F_adult[i];
  mutate(Q_kidneyFraction_F = (V_kidneyFraction_F/0.0046)*0.175) %>% # sc_F[0-17] = (sc_V[i]  / sc_V_adult[i])  * sc_F_adult[i];
  
  # Lungs; compartment [15] in Ratier 2024 (not used in our model, but needed for calculation of adipose tissue)
  # Comment Chrysa on 18-10-2024: Shouldn't the lungs take 100% of the blood flow?
  mutate(Q_lungFraction_M = (V_lungFraction_M/0.0068)*0.026) %>% # sc_F[0-17] = (sc_V[i]  / sc_V_adult[i])  * sc_F_adult[i];
  mutate(Q_lungFraction_F = (V_lungFraction_F/0.0070)*0.026) %>% # sc_F[0-17] = (sc_V[i]  / sc_V_adult[i])  * sc_F_adult[i];
  
  # Gut; compartment [16] in Ratier 2024
  mutate(Q_gutFraction_M = (V_gutFraction_M/0.0140)*0.144) %>% # sc_F[0-17] = (sc_V[i]  / sc_V_adult[i])  * sc_F_adult[i];
  mutate(Q_gutFraction_F = (V_gutFraction_F/0.0160)*0.165) %>% # sc_F[0-17] = (sc_V[i]  / sc_V_adult[i])  * sc_F_adult[i];
  
  # Stomach; compartment [17] in Ratier 2024
  mutate(Q_stomachFraction_M = (V_stomachFraction_M/0.0021)*0.01) %>% # sc_F[0-17] = (sc_V[i]  / sc_V_adult[i])  * sc_F_adult[i];
  mutate(Q_stomachFraction_F = (V_stomachFraction_F/0.0023)*0.01) %>% # sc_F[0-17] = (sc_V[i]  / sc_V_adult[i])  * sc_F_adult[i];
  
  # Liver artery; compartments [18] (liver) and [21] (liver artery) in Ratier 2024
  mutate(Q_liverFraction_M = (V_liverFraction_M/0.0247)*0.065) %>% # sc_F[21] = (sc_V[18] / sc_V_adult[18]) * sc_F_adult[21];
  mutate(Q_liverFraction_F = (V_liverFraction_F/0.0233)*0.065) %>% # sc_F[21] = (sc_V[18] / sc_V_adult[18]) * sc_F_adult[21];
  
  # Adipose tissue
  mutate(Q_adiposeFraction_M = (V_adiposeFraction_M/0.20)*0.052) %>% # sc_F[0-17] = (sc_V[i]  / sc_V_adult[i])  * sc_F_adult[i];
  mutate(Q_adiposeFraction_F = (V_adiposeFraction_F/0.3167)*0.087) # sc_F[0-17] = (sc_V[i]  / sc_V_adult[i])  * sc_F_adult[i];
write.csv(Variables_df, "PhysioVariables.csv", row.names = FALSE)


# ## Check mass balance volumes and flows ####
# Variables_M_df = Variables_df %>%
#   select(.,ends_with("_M")) %>%
#   mutate(BloodFlowSum = rowSums(select(., starts_with("Q_")))) %>%
#   mutate(VolumesSum = rowSums(select(., starts_with("V_")))) %>%
#   mutate(TIME = TIME) %>%
#   mutate(age = TIME/365)
# 
# Variables_F_df = Variables_df %>%
#   select(.,ends_with("_F")) %>%
#   mutate(BloodFlowSum = rowSums(select(., starts_with("Q_")))) %>%
#   mutate(VolumesSum = rowSums(select(., starts_with("V_")))) %>%
#   mutate(TIME = TIME) %>%
#   mutate(age = TIME/365)
# 
# PLOT_VolumeTotal =
#   ggplot() +
#   geom_path(data = Variables_M_df, aes(age, VolumesSum), colour = "lavenderblush4") +
#   geom_path(data = Variables_F_df, aes(age, VolumesSum), colour = "purple")
# PLOT_VolumeTotal # Is 1
# 
# PLOT_BloodFlowTotal =
#   ggplot()+
#   geom_path(data = Variables_M_df, aes(age, BloodFlowSum), colour = "lavenderblush4") +
#   geom_path(data = Variables_F_df, aes(age, BloodFlowSum), colour = "purple")
# PLOT_BloodFlowTotal # Not 1; shouldn't it be 1?
# 
# PLOT_BWChanges =
#   ggplot()+
#   geom_path(data = Variables_M_df, aes(age, BDW_M, color = "Male")) +
#   geom_path(data = Variables_F_df, aes(age, BDW_F, color = "Female")) +
#   scale_color_manual(values = c("Male" = "mediumpurple",
#                                 "Female" = "mediumpurple4"),
#                      name = "") +
#   ylab("Body Weight (Kg)") +
#   xlab("Age (years)") +
#   theme_CP()
# PLOT_BWChanges
# 
# PLOT_CardiacOutChanges =
#   ggplot()+
#   geom_path(data = Variables_M_df, aes(age, CardOut_M, color = "Male")) +
#   geom_path(data = Variables_F_df, aes(age, CardOut_F, color = "Female")) +
#   scale_color_manual(values = c("Male" = "mediumpurple",
#                                 "Female" = "mediumpurple4"),
#                      name = "") +
#   ylab("Cardiac Output (L/h)") +
#   xlab("Age (years)") +
#   theme_CP()
# PLOT_CardiacOutChanges
# 
# ## Actual Volumes and Blood Flows ####
# MaleVariables_df = Variables_df %>%
#   
#   # Final organ volumes
#   select(ends_with('_M')) %>%
#   select(!c(V_boneFraction_M, V_bonenonperfusedFraction_M, V_adiposeFraction_M, V_plasmaFraction_M)) %>%
#   mutate(across(starts_with('V_'), ~ . * BDW_M)) %>% # Final Organ Volume = Fractional Organ Volume * Body Weight (age specific)
#   mutate(V_boneFraction_M = Variables_df$V_boneFraction_M*BDW_M/2, # 2 is the bone density
#          V_bonenonperfusedFraction_M = Variables_df$V_bonenonperfusedFraction_M*BDW_M/2, # 2 is the bone density
#          V_adiposeFraction_M = (Variables_df$V_adiposeFraction_M*BDW_M/0.9) + (AdiposeMass_M/0.9), # 0.9 is the bone density
#          V_plasmaFraction_M = Variables_df$V_plasmaFraction_M*BDW_M*(1-Variables_df$Hct_M)) %>%  # (1-Hematocrite) corrects for the plasma volume (if not is total blood)
#   mutate(TotalVolume_M = rowSums(select(., starts_with("V_")))) %>%
#   
#   # Final blood flow to organs
#   mutate(TotalBloodFlow_M = Variables_M_df$BloodFlowSum) %>%
#   mutate(across(starts_with("Q_"), ~ . /TotalBloodFlow_M * CardOut_M)) %>% # Final Blood Flow = Fractional Blood Flow/Total Fractional Blood Flow * Cardiac Output (age specific)
#   mutate(TotalBloodFlow_M = rowSums(select(., starts_with("Q_")))) %>%
#   select(BDW_M, TotalVolume_M, Hct_M, CardOut_M, TotalBloodFlow_M, (matches("(Q_|V_)"))) %>%
#   mutate(TIME = Variables_df$TIME,
#          age = Variables_df$age)
# 
# colnames(MaleVariables_df) = str_remove(colnames(MaleVariables_df), "Fraction")
# write.csv(MaleVariables_df, "MaleVariables.csv", row.names = FALSE)
# 
# FemaleVariables_df = Variables_df %>%
#   
#   # Final organ volumes
#   select(ends_with('_F')) %>%
#   select(!c(V_boneFraction_F, V_bonenonperfusedFraction_F, V_adiposeFraction_F, V_plasmaFraction_F)) %>%
#   mutate(across(starts_with('V_'), ~ . * BDW_F)) %>% # Final Organ Volume = Fractional Organ Volume * Body Weight (age specific)
#   mutate(V_boneFraction_F = Variables_df$V_boneFraction_F*BDW_F/2, # 2 is the bone density
#          V_bonenonperfusedFraction_F = Variables_df$V_bonenonperfusedFraction_F*BDW_F/2, # 2 is the bone density
#          V_adiposeFraction_F = (Variables_df$V_adiposeFraction_F*BDW_F/0.9) + (AdiposeMass_F/0.9), # 0.9 is the bone density
#          V_plasmaFraction_F = Variables_df$V_plasmaFraction_F*BDW_F*(1-Variables_df$Hct_F)) %>%  # (1-Hematocrite) corrects for the plasma volume (if not is total blood)
#   mutate(TotalVolume_F = rowSums(select(., starts_with("V_")))) %>%
# 
#   # Final blood flow to organs
#   mutate(TotalBloodFlow_F = Variables_F_df$BloodFlowSum) %>%
#   mutate(across(starts_with("Q_"), ~ . /TotalBloodFlow_F * CardOut_F)) %>% # Final Blood Flow = Fractional Blood Flow/Total Fractional Blood Flow * Cardiac Output (age specific)
#   mutate(TotalBloodFlow_F = rowSums(select(., starts_with("Q_")))) %>%
#   select(BDW_F, TotalVolume_F, Hct_F, CardOut_F, TotalBloodFlow_F, (matches("(Q_|V_)"))) %>%
#   mutate(TIME = Variables_df$TIME,
#          age = Variables_df$age)
# 
# colnames(FemaleVariables_df) = str_remove(colnames(FemaleVariables_df), "Fraction")
# write.csv(FemaleVariables_df, "FemaleVariables.csv", row.names = FALSE)
# 
### Notes ####

# Difference from original Aude's code:
# Here we have one compartment for plasma, while in Aude's paper it was divided in arterial and venous, see notes below
# ## Plasma volume; compartment [22] in Ratier 2024 -> USED TO BE BLOOD volume, as it's corrected for hematocrit then it's plasma
# mutate(Vart_M = Vplasma_M*Fr_art_plasma*(1-Hct_M)) %>% #arterial blood (plasma)
# mutate(Vven_M = Vplasma_M*(1-Hct_M) - Vart_M) %>% #venous blood (plasma)#
# mutate(Vart_F = Vplasma_F*Fr_art_plasma*(1-Hct_F)) %>%
# mutate(Vven_F = Vplasma_F*(1-Hct_F) - Vart_F) %>%

### Plots ####
# PLOT_VolumeChanges =
#   ggplot()+
#   geom_path(data = FemaleVariables_df, aes(age, V_liver_F, color = "Liver")) +
#   geom_path(data = FemaleVariables_df, aes(age, V_kidney_F, color = "Kidney")) +
#   # geom_path(data = FemaleVariables_df, aes(age, V_adipose_F, color = "Adipose")) +
#   geom_path(data = FemaleVariables_df, aes(age, V_gut_F, color = "Gut")) +
#   geom_path(data = FemaleVariables_df, aes(age, V_skin_F, color = "Skin")) +
#   scale_color_manual(values = c("Liver" = "goldenrod2",
#                                 "Kidney" = "royalblue2",
#                                 #"Adipose" = "seagreen2",
#                                 "Gut" = "turquoise2",
#                                 "Skin" = "deeppink2"),
#                      name = "") +
#   ylab("Organ Volumes (Kg)") +
#   xlab("Age (years)") +
#   theme_CP()
# PLOT_VolumeChanges
# 
# PLOT_FlowChanges =
#   ggplot()+
#   geom_path(data = FemaleVariables_df, aes(age, Q_liver_F, color = "Liver")) +
#   geom_path(data = FemaleVariables_df, aes(age, Q_kidney_F, color = "Kidney")) +
#   # geom_path(data = FemaleVariables_df, aes(age, Q_adipose_F, color = "Adipose")) +
#   geom_path(data = FemaleVariables_df, aes(age, Q_gut_F, color = "Gut")) +
#   geom_path(data = FemaleVariables_df, aes(age, Q_skin_F, color = "Skin")) +
#   scale_color_manual(values = c("Liver" = "goldenrod2",
#                                 "Kidney" = "royalblue2",
#                                 #"Adipose" = "seagreen2",
#                                 "Gut" = "turquoise2",
#                                 "Skin" = "deeppink2"),
#                      name = "") +
#   ylab("Organ Blood Flows (L/h)") +
#   xlab("Age (years)") +
#   theme_CP()
# PLOT_FlowChanges



















































# PHYSIOLOGICAL, PHYSICOCHEMICAL and BIOKINETIC PARAMETERS ####
# ------------------------------------------------------ #
# To directly call .csv files instead of re-running the life-stage code below everytime
PhysioVariables_M_df = read_csv("C:/Users/pacho003/OneDrive - Wageningen University & Research/C Channel/R/PARC_PFAS_PBPKmodel/Codes_PFAS_models/Chrysanthi_Pachoulide_2024_PFOA_oral_dermal_blood/Input/2024-11-14/PhysioVariablesMale.csv") %>% as.data.frame()
# PhysioVariables_F_df = read_csv("C:/Users/pacho003/OneDrive - Wageningen University & Research/C Channel/R/PARC_PFAS_PBPKmodel/Codes_PFAS_models/Chrysanthi_Pachoulide_2024_PFOA_oral_dermal_blood/Input/2024-11-14/PhysioVariablesFemale.csv") %>% as.data.frame()

# PFOA
MW = 414.07

# Fraction unbound in plasma

fup = 0.061/100 # Fischer et al. 2024 https://doi.org/10.1021/acs.est.3c07415; OLD fup = 0.02 # fup fraction of PFOA in plasma


# Partition coefficients

# These are the values from Trine's model
# Note Chrysa 24-20-2024: plasma/tissue partition coefficients to be changed using the Allendorf paper:  https://doi.org/10.1002/etc.4954
PF = 0.04  # Plasma/fat partition coefficient; Rat tissue data (Kudo et al. 2007
PG = 0.05  # Plasma/gut partition coefficient; Rat tissue data (Kudo et al. 2007 %>%
PK = 1.05  # Plasma/kidney partition coefficient; Rat tissue data (Kudo et al. 2007
PL = 2.2  # Plasma/liver partition coefficient; Rat tissue data (Kudo et al. 2007
PSk = 0.1  # Plasma/skin partition coefficient; Rat tissue data (Kudo et al. 2007
PR = 0.12  # Plasma/rest of the body partition coefficient; Rat tissue data (Kudo et al. 2007

# Calculated partition coefficients
Kp.df <- read_csv("C:/Users/pacho003/OneDrive - Wageningen University & Research/C Channel/R/PARC_PFAS_PBPKmodel/Codes_PFAS_models/Chrysanthi_Pachoulide_2024_PFOA_oral_dermal_blood/Input/PFOAPBK.PartCoefs.csv")
  
KpAd = Kp.df$KpAd
KpGu = Kp.df$KpGu
KpKi = Kp.df$KpKi
KpLi = Kp.df$KpLi
KpSk = Kp.df$KpSk
KpRe = Kp.df$KpRe

# Dermal absorption

# Comment Chrysa 21-10-2024: AbsPFOA is used in the Dermaldose input so we're already correcting before using the Papp?! not sure I agree with this
fBAc_dust = 69.6/100 # 69.6±5 # Fraction bio accessible from dust; Ragnarsdottir et al.: https://doi.org/10.1016/j.envres.2023.117093; fraction bioaccessible = fraction of the compound released from the matrix (cosmetic formulation dust etc and is available to be absorbed from the epidermis
fBAc_cosm = 97.8/100 # 97.8±19.0 !! this is an estimated value, based on [PFOA]_cosmetics and the LOQ # Fraction bio accessible from cosmetics; Namazkar et al.: DOI: 10.1039/d3em00461a
AbsPFOA = 0.016 # from Trine; 0.00048 # Changed to the absorption measured by Abraham and Monien 2022 of 1.6% of applied dose from sunscreen.
Papp = 3.82 * 10^-3 # cm/h ref: https://doi.org/10.1016/j.envint.2024.108772


# Kidney clearance

# Comment Chrysa 05-11-2024:
## Assuming that the kidney PFAS concentrations never reaches Km concentrations therefore transforming Vmax and Km to a Clearance; in the paper what they call transporter efficiency : Louisse Pedroni et al. 2024 https://doi.org/10.1016/j.tox.2024.153961
## OAT1 and OAT3 are determining transport between blood and proximal tubule
## OAT4 is determining transport between proximal tubule and filtrate
## All transporters are bi-directional therefore the equation is written assuming that direction is determined based on the equilibrium between the two compartments
## Introducing the "affinity constants: kAap and kAbl to compensate for the fact that the transporters have an affinity to one side.
kAbl = 0.01 # affinity constant basolateral this is about OAT1 and OAT3 which have affinity to uptake (movement from plasma to cells; this is fitted value for now; kAbl = 0.01 is driving the equilibrium towards uptake into the proximal tubule cells
kAap = 0.01 # affinity constant apical this is about OAT4 which has affinity to re-abs (movement from filtrate to cells; this is fitted value for now; kAap = 0.01 is driving the equilibrium towards re-absorption into the proximal tubule cells

# Comment Chrysa 12-11-2024: I'm not sure about these scaling factors as the in vitro clearance is expressed in ul/min/mg protein and not ul/min/HEK cells
CL_OAT1 = 19/ 10^-6 *60*24 * 10^-6 # L/d/kg protein; initial ul/min/mg protein
CL_OAT3 = 17/ 10^-6 *60*24 * 10^-6 # L/d/kg protein; initial ul/min/mg protein
CL_OAT4 = 96/ 10^-6 *60*24 * 10^-6 # L/d/kg protein; initial ul/min/mg protein

PTCPGK = 9.94* 10^7 * 10^3 # proximal tubule cells/kg kidney cortex; initial 99.4 million PTC/g kidney https://doi.org/10.1021/acs.molpharmaceut.4c00504

InVivo_OAT1 = 5.3 # 5.3 ± 1.9 OAT1 /mg membrane protein in the human kidney cortex http://dx.doi.org/10.1124/dmd.116.072066; alternative : 4.3 ± 0.3 pmol OAT1 /mg membrane protein in the human kidney cortex https://doi.org/10.1124/dmd.121.000367;
InVitro_OAT1 = 26.6 # 26.6 ± 3.4 pmol/mg membrane protein: OAT1 expression in HEK293-OAT1 cells https://doi.org/10.1124/dmd.121.000367;
InVivo_OAT3 = 3.5 # 3.5 ± 1.6 OAT3/mg membrane protein in the human kidney cortex http://dx.doi.org/10.1124/dmd.116.072066; alternative : 2.7 ± 0.1 pmol OAT3 /mg membrane protein in the human kidney cortex https://doi.org/10.1124/dmd.121.000367;
InVitro_OAT3 = 7.3 # 7.3 ± 0.5 pmol/mg membrane protein: OAT3 expression in HEK293-OAT1 cells https://doi.org/10.1124/dmd.121.000367;
InVivo_OAT4 = 0.52 # 0.52 ± 0.23 pmol OAT4 /mg membrane protein in the human kidney cortex http://dx.doi.org/10.1124/dmd.116.072066;
# InVitro_OAT4 = (InVitro_OAT1 + InVitro_OAT3)/2 # OAT4 expression in HEK293-OAT4 cells not found so assuming average of OAT1 and OAT3 expression


REF_OAT1 = InVivo_OAT1/InVitro_OAT1 # relative expression factor OAT1: expression in the human kidneys /expression in the cells
REF_OAT3 = InVivo_OAT3/InVitro_OAT3 # relative expression factor OAT3: expression in the human kidneys /expression in the cells
REF_OAT4 = 1 # assuming 1 as Invitro expression not found InVivo_OAT4/InVitro_OAT4 # relative expression factor OAT4: expression in the human kidneys /expression in the cells

# CLurinec = 0.000044  # L/d/kg; 0.044 mL/d/kg taken from Fujii et al 2015
# Vmaxc = 4.5*MW/1000*60*24 # ug/d/mg protein; 45 nmol/min/mg protein *MW/1000*60*24 = ug/d/mg protein ref: Louisse et al. 2023 https://doi.org/10.1007/s00204-022-03428-6
# Km = 47*MW # ug/L; 47 uM*MW = ug/L ref: Louisse et al. 2023


# Hepatic clearance
# These are the values from Trine's model
CLbiliaryc = 0.00262 # L/d/kg ; 2.62 +/- 3.6 mL/d/kg from Fujii et al 2015 DOI: 10.1539/joh.14-0136-OA
CLfaecesc = 0.000052 # L/d/kg ; 0.052 +/- 0.05 mL/d/kg clearance in faeces taken from Fujii et al 2015 DOI: 10.1539/joh.14-0136-OA

# Male
Final_variables_M_df <- PhysioVariables_M_df %>%

  # Independent variables
  mutate(
    AbsPFOA = AbsPFOA, # fraction absorbed, from Trine
    kAbl = kAbl, # affinity constant basolateral this is about OAT1 and OAT3 which have affinity to uptake (movement from plasma to cells; this is fitted value for now; kAbl = 0.01 is driving the equilibrium towards uptake into the proximal tubule cells
    kAap = kAap, # affinity constant apical this is about OAT4 which has affinity to re-abs (movement from filtrate to cells; this is fitted value for now; kAap = 0.01 is driving the equilibrium towards re-absorption into the proximal tubule cells
    fup = fup, # Unbound fraction in plasmal Fischer et al. 2024 https://doi.org/10.1021/acs.est.3c07415;
    PL = PL,  # Plasma/liver partition coefficient; Rat tissue data (Kudo et al. 2007)
    PF = PF, # Plasma/fat partition coefficient; Rat tissue data (Kudo et al. 2007)
    PK = PK,  # Plasma/kidney partition coefficient; Rat tissue data (Kudo et al. 2007)
    PSk = PSk,  # Plasma/skin partition coefficient; Rat tissue data (Kudo et al. 2007)
    PR = PR,  # Plasma/rest of the body partition coefficient; Rat tissue data (Kudo et al. 2007)
    PG = PG,  # Plasma/gut partition coefficient; Rat tissue data (Kudo et al. 2007)
    KpAd = KpAd, # Calculated
    KpGu = KpGu, # Calculated
    KpKi = KpKi, # Calculated
    KpLi = KpLi, # Calculated
    KpSk = KpSk, # Calculated
    KpRe = KpRe, # Calculated
    Papp = Papp, # Calculated
    CL_OAT1 = CL_OAT1, 
    REF_OAT1 = REF_OAT1,
    
  ) %>%
  
  # Age-dependent variables
  mutate(
    SA_hands_M = SA_hands_M, #Surface area of hands
  )
  # # Dermal absorption
  # mutate(CLdermalabs_M = ((Papp*SA_hands_M)/1000)*24) %>% # (L/d); cm/h*cm^2 = mL/h /1000 = L/h * 24 = L/d
  # 
  # # Kidney clearance
  # mutate(# CL_PltPT_M = ((CL_OAT1*REF_OAT1) + (CL_OAT3*REF_OAT3)) * PTC_kidneyTissue_M, # L/d plasma to proximal tubule clearance
  #        CL_PltPT_M = ((CL_OAT1*REF_OAT1) + (CL_OAT3*REF_OAT3)) * 0.17 * 0.7 * V_kidney_M,# L/d scaling for protein content instead of cell content; 17% of kidney is protein and 70% of kidney is cortex [ICRP 89], assuming that all the kidney protein is found in the cortex; this is an overestimation though; double ref for 17% protein Ruark 2020: DOI: https://doi.org/10.1016/B978-0-12-818596-4.00006-0
  #        # CL_FiltPT_M = (CL_OAT4*REF_OAT4) * PTC_kidneyTissue_M, #L/d filtrate to proximal tubule clearance
  #        CL_FiltPT_M = CL_OAT4 * REF_OAT4 * 0.17 * 0.7 * V_kidney_M # L/d scaling for protein content instead of cell content; 17% of kidney is protein and 70% of kidney is cortex [ICRP 89], assuming that all the kidney protein is found in the cortex; this is an overestimation though; double ref for 17% protein Ruark 2020: DOI: https://doi.org/10.1016/B978-0-12-818596-4.00006-0
  #        # Trine's values
  #        # CLurine_M = CLurinec*BDW_M^(-0.25) #from Husoy; L/d clearance urine
  # ) %>%
  # 
  # # Biliary clearance
  # mutate(CLbiliary_M = CLbiliaryc*(BDW_M^0.1)) %>% #from Husoy; L/d biliary clearance rate
  # 
  # # Fecal clearance
  # mutate(CLfecal_M = CLfaecesc*(BDW_M^0.001)) #from Husoy; L/d faeces clearance rate

write.csv(Final_variables_M_df, "Final_variables_M.csv", row.names = FALSE)

# Female
Final_variables_F_df <- PhysioVariables_F_df %>%
  
  # Independent variables
  mutate(
    AbsPFOA = AbsPFOA, # fraction absorbed, from Trine
    kAbl = kAbl, # affinity constant basolateral this is about OAT1 and OAT3 which have affinity to uptake (movement from plasma to cells; this is fitted value for now; kAbl = 0.01 is driving the equilibrium towards uptake into the proximal tubule cells
    kAap = kAap, # affinity constant apical this is about OAT4 which has affinity to re-abs (movement from filtrate to cells; this is fitted value for now; kAap = 0.01 is driving the equilibrium towards re-absorption into the proximal tubule cells
    fup = fup, # Unbound fraction in plasmal Fischer et al. 2024 https://doi.org/10.1021/acs.est.3c07415;
    PL = PL,  # Plasma/liver partition coefficient; Rat tissue data (Kudo et al. 2007)
    PF = PF, # Plasma/fat partition coefficient; Rat tissue data (Kudo et al. 2007)
    PK = PK,  # Plasma/kidney partition coefficient; Rat tissue data (Kudo et al. 2007)
    PSk = PSk,  # Plasma/skin partition coefficient; Rat tissue data (Kudo et al. 2007)
    PR = PR,  # Plasma/rest of the body partition coefficient; Rat tissue data (Kudo et al. 2007)
    PG = PG,  # Plasma/gut partition coefficient; Rat tissue data (Kudo et al. 2007)
    KpAd = KpAd, # Calculated
    KpGu = KpGu, # Calculated
    KpKi = KpKi, # Calculated
    KpLi = KpLi, # Calculated
    KpSk = KpSk, # Calculated
    KpRe = KpRe
  ) %>%
  
  # Dermal absorption
  mutate(CLdermalabs_F = ((Papp*SA_hands_F)/1000)*24) %>% # (L/d); cm/h*cm^2 = mL/h /1000 = L/h * 24 = L/d
  
  # Kidney clearance
  mutate(# CL_PltPT_F = ((CL_OAT1*REF_OAT1) + (CL_OAT3*REF_OAT3)) * PTC_kidneyTissue_F, # L/d plasma to proximal tubule clearance
    CL_PltPT_F = ((CL_OAT1*REF_OAT1) + (CL_OAT3*REF_OAT3)) * 0.17 * 0.7 * V_kidney_F,# L/d scaling for protein content instead of cell content; 17% of kidney is protein and 70% of kidney is cortex [ICRP 89], assuming that all the kidney protein is found in the cortex; this is an overestimation though; double ref for 17% protein Ruark 2020: DOI: https://doi.org/10.1016/B978-0-12-818596-4.00006-0
    # CL_FiltPT_F = (CL_OAT4*REF_OAT4) * PTC_kidneyTissue_F, #L/d filtrate to proximal tubule clearance
    CL_FiltPT_F = CL_OAT4 * REF_OAT4 * 0.17 * 0.7 * V_kidney_F # L/d scaling for protein content instead of cell content; 17% of kidney is protein and 70% of kidney is cortex [ICRP 89], assuming that all the kidney protein is found in the cortex; this is an overestimation though; double ref for 17% protein Ruark 2020: DOI: https://doi.org/10.1016/B978-0-12-818596-4.00006-0
    # Trine's values
    # CLurine_F = CLurinec*BDW_F^(-0.25) #from Husoy; L/d clearance urine
  ) %>%
  
  # Biliary clearance
  mutate(CLbiliary_F = CLbiliaryc*(BDW_F^0.1)) %>% #from Husoy; L/d biliary clearance rate
  
  # Fecal clearance
  mutate(CLfecal_F = CLfaecesc*(BDW_F^0.001)) #from Husoy; L/d faeces clearance rate

write.csv(Final_variables_F_df, "Final_variables_F.csv", row.names = FALSE)

# #Plot differences in partition coefficients
# 
# KpCompare.df <- data.frame(
#   Tissue = c("Adipose", "Gut", "Kidney", "Liver", "Skin", "Rest"),
#   Method = c("InVivoRat", "InVivoRat", "InVivoRat", "InVivoRat", "InVivoRat", "InVivoRat", "Calculated", "Calculated", "Calculated", "Calculated", "Calculated", "Calculated"),
#   Coefficient = c(PF, PG, PK, PL, PSk, PR, KpAd, KpGu, KpKi, KpLi, KpSk, KpRe)
# )
# 
# KpCompare.plot <- KpCompare.df %>% 
#   ggplot(aes(x = Tissue, y = Coefficient, fill = Method)) +
#   geom_bar(stat = "identity", position = "dodge") +
#   labs(
#     x = "",
#     y = "Partition Coefficient") +
#   theme_minimal() +
#   theme(
#     axis.text.x = element_text(angle = 45, hjust = 1)
#   ) +
#   scale_fill_brewer(palette = "Accent")
# KpCompare.plot



# PFAS MODEL SPECIFIC PHYSIOLOGICAL DATA ####
# Call .csv files instead of re-running the life-stage code below everytime
MaleVariables_df <- read_csv("C:/Users/pacho003/OneDrive - Wageningen University & Research/C Channel/R/PARC_PFAS_PBPKmodel/Codes_PFAS_models/Chrysanthi_Pachoulide_2024_PFOA_oral_dermal_blood/Input/2024-11-12/MaleVariables.csv") %>% as.data.frame()
FemaleVariables_df <- read_csv("C:/Users/pacho003/OneDrive - Wageningen University & Research/C Channel/R/PARC_PFAS_PBPKmodel/Codes_PFAS_models/Chrysanthi_Pachoulide_2024_PFOA_oral_dermal_blood/Input/2024-11-12/FemaleVariables.csv") %>% as.data.frame()

# Model compartments: Skin (Barrier and Skin -to become Plasma and Tissue-), Kidney (Plasma, Tissue, Filtrate), Gut, Liver, Plasma, Adipose

# Male
PhysioVariables_M_df = MaleVariables_df %>%
  
  # Rest compartment: lumps all organs that are not specified in the model
  mutate(V_rest_M = TotalVolume_M - rowSums(
    select(., starts_with("V_")) %>%
      select(contains(c("skin", "kidney", "gut", "liver", "plasma", "adipose")))
    )) %>%
  mutate(Q_rest_M = TotalBloodFlow_M - rowSums(
    select(., starts_with("Q_")) %>%
      select(contains(c("skin", "kidney", "gut", "liver", "adipose")))
  )) %>%

  # Comment Chrysa 12-11-2024: Adding hepatic artery blood flow
  mutate(Q_hepatic_M = rowSums( # hepatic artery = liver blood flow + portal blood flow
    select(., starts_with("Q_")) %>%
      select(contains(c("liver", "gut"))) # Q spleen, stomach, gut and pancreas is the portal artery blood flow to the liver; ignoring spleen, stomach and pancreas as they are not explicit compartments of the model
  )) %>%

  # Kidney compartment is permeability limited => divided in Plasma, Tissue, Filtrate and Urine compartments
  # Comment Chrysa 12-11-2024: I think the blood volume be corrected with Hct also to call it plasma right?!
  mutate(V_kidneyPlasma_M = V_kidney_M*0.36, # Volume fraction of blood in the kidneys 0.36+-0.01 [Brown 1997 table 30]
         V_kidneyTissue_M = V_kidney_M*0.64,
         V_filtrate_M = V_kidney_M*0.05, # Kidney filtrate compartment, corresponds to the volume of the collecting system in [ICRP 89 page 149] http://www.icrp.org/publication.asp?id=ICRP%20Publication%2089
         QUr_M = 0.022*BDW_M, # L/d, Urine flow rate to the bladder 22 mL/kg BW/d [ICRP 89 page 161]
         GFR_M = 0.18*Q_kidney_M, # L/d, Glomerular filtration rate 18% of total renal plasma flow [ICRP 89 page 159] http://www.icrp.org/publication.asp?id=ICRP%20Publication%2089
         # Needed for scaling clearances
         # Comment Chrysa 12-11-2024: these scaling factors assume that the cellularity of the kidneys is constant through the lifestages
         KW_cortex_M = 0.7*V_kidney_M, # Kg kidney cortexes, only scaling to kidney cortex volume as proximal tubule cells are in the cortex; 70% of the total kidney volume according to ICRP89; PT are in the cortex https://doi.org/10.1021/acs.molpharmaceut.4c00504; alternatively we could have 68% of kidney weight https://doi.org/10.1124/dmd.117.075242
         PTCPGK = 99.4 / 1000, # proximal tubule cells/kg kidney cortex initial PTC/g kidney https://doi.org/10.1021/acs.molpharmaceut.4c0050
         PTC_kidneyTissue_M = PTCPGK*KW_cortex_M #actual number of cells in the kidney
         ) %>%

  # Skin compartment
  mutate(# Needed for scaling Papp
         # Comment Chrysa 12-11-2024: According to ICRP 89 page 64 skin surface area depends on both height and BDW. Consider changing SA_skin_M parameter in the future..
         SA_skin_M = 9.1*((BDW_M*1000)^0.666), # Used to be called SkbTarea_M; Total surface area of the skin; taken from Trine's model
         SA_hands_M = 0.05*SA_skin_M, # Used to be called fSkbarea_M; 0.05 used to be called skin_fraction; Surface area of the hands; hands are 5% https://www.epa.gov/sites/default/files/2015-09/documents/efh-chapter07.pdf
         Skb_thickness_M = 83.1/1000, # Used to be called Skbthickness_M (cm) ref: DOI: 10.1080/00015550310015419; in Husoy this was 0.1
         V_skinBarrier_M = (SA_hands_M*Skb_thickness_M)/1000, # (L); Skin barrier volume; as previously coded by Trine

         # Permeability limited skin compartment => divided in Plasma and Tissue compartments; for future use
         V_skinPlasma_M = V_skin_M*0.08, # Not used yet; to be used when updating skin. Volume fraction of blood in skin 0.08 [Brown 1997 table 30]
         V_skinTissue_M = V_skin_M*0.92 # Not used yet; to be used when updating skin.
                  ) %>%

  select(TIME, age,
         BDW_M, TotalVolume_M, Hct_M, CardOut_M, TotalBloodFlow_M,
         QUr_M, GFR_M, PTC_kidneyTissue_M, SA_hands_M,
         contains(c("skin", "kidney", "filtrate", "gut", "liver", "hepatic", "plasma", "adipose", "rest")))
write.csv(PhysioVariables_M_df, "PhysioVariablesMale.csv", row.names = FALSE)


# Female
PhysioVariables_F_df = FemaleVariables_df %>%

  # Rest compartment: lumps all organs that are not specified in the model
  mutate(V_rest_F = TotalVolume_F - rowSums(
    select(., starts_with("V_")) %>%
      select(contains(c("skin", "kidney", "gut", "liver", "plasma", "adipose")))
  )) %>%
  mutate(Q_rest_F = TotalBloodFlow_F - rowSums(
    select(., starts_with("Q_")) %>%
      select(contains(c("skin", "kidney", "gut", "liver", "adipose")))
  )) %>%

  # Comment Chrysa 12-11-2024: Adding hepatic artery blood flow
  mutate(Q_hepatic_F = rowSums( # hepatic artery = liver blood flow + portal blood flow
    select(., starts_with("Q_")) %>%
      select(contains(c("liver", "gut"))) # Q spleen, stomach, gut and pancreas is the portal artery blood flow to the liver; ignoring spleen, stomach and pancreas as they are not explicit compartments of the model
  )) %>%

  # Kidney compartment is permeability limited => divided in Plasma, Tissue, Filtrate and Urine compartments
  # Comment Chrysa 12-11-2024: I think the blood volume be corrected with Hct also to call it plasma right?!
  mutate(V_kidneyPlasma_F = V_kidney_F*0.36, # Volume fraction of blood in the kidneys 0.36+-0.01 [Brown 1997 table 30]
         V_kidneyTissue_F = V_kidney_F*0.64,
         V_filtrate_F = V_kidney_F*0.05, # Kidney filtrate compartment, corresponds to the volume of the collecting system in [ICRP 89 page 149] http://www.icrp.org/publication.asp?id=ICRP%20Publication%2089
         QUr_F = 0.022*BDW_F, # L/d, Urine flow rate to the bladder 22 mL/kg BW/d [ICRP 89 page 161]
         GFR_F = 0.18*Q_kidney_F, # L/d, Glomerular filtration rate 18% of total renal plasma flow [ICRP 89 page 159] http://www.icrp.org/publication.asp?id=ICRP%20Publication%2089
         # Needed for scaling clearances
         # Comment Chrysa 12-11-2024: these scaling factors assume that the cellularity of the kidneys is constant through the lifestages
         KW_cortex_F = 0.7*V_kidney_F, # Kg kidney cortexes, only scaling to kidney cortex volume as proximal tubule cells are in the cortex; 70% of the total kidney volume according to ICRP89; PT are in the cortex https://doi.org/10.1021/acs.molpharmaceut.4c00504; alternatively we could have 68% of kidney weight https://doi.org/10.1124/dmd.117.075242
         PTCPGK = 99.4 / 1000, # proximal tubule cells/kg kidney cortex initial PTC/g kidney https://doi.org/10.1021/acs.molpharmaceut.4c0050
         PTC_kidneyTissue_F = PTCPGK*KW_cortex_F #actual number of cells in the kidney
  ) %>%

  # Skin compartment
  mutate(# Needed for scaling Papp
    # Comment Chrysa 12-11-2024: According to ICRP 89 page 64 skin surface area depends on both height and BDW. Consider changing SA_skin_F parameter in the future..
    SA_skin_F = 9.1*((BDW_F*1000)^0.666), # Used to be called SkbTarea_F; Total surface area of the skin; taken from Trine's model
    SA_hands_F = 0.05*SA_skin_F, # Used to be called fSkbarea_F; 0.05 used to be called skin_fraction; Surface area of the hands; hands are 5% https://www.epa.gov/sites/default/files/2015-09/documents/efh-chapter07.pdf
    Skb_thickness_F = 83.1/1000, # Used to be called Skbthickness_F (cm) ref: DOI: 10.1080/00015550310015419; in Husoy this was 0.1
    V_skinBarrier_F = (SA_hands_F*Skb_thickness_F)/1000, # (L); Skin barrier volume; as previously coded by Trine

    # Permeability limited skin compartment => divided in Plasma and Tissue compartments; for future use
    V_skinPlasma_F = V_skin_F*0.08, # Not used yet; to be used when updating skin. Volume fraction of blood in skin 0.08 [Brown 1997 table 30]
    V_skinTissue_F = V_skin_F*0.92 # Not used yet; to be used when updating skin.
  ) %>%

  select(TIME, age,
         BDW_F, TotalVolume_F, Hct_F, CardOut_F, TotalBloodFlow_F,
         QUr_F, GFR_F, PTC_kidneyTissue_F, SA_hands_F,
         contains(c("skin", "kidney", "filtrate", "gut", "liver", "hepatic", "plasma", "adipose", "rest")))
write.csv(PhysioVariables_F_df, "PhysioVariablesFemale.csv", row.names = FALSE)




# Male
Final_variables_M_df <- PhysioVariables_M_df %>%
  
  # Independent variables
  mutate(
    AbsPFOA = AbsPFOA, # fraction absorbed, from Trine
    kAbl = kAbl, # affinity constant basolateral this is about OAT1 and OAT3 which have affinity to uptake (movement from plasma to cells; this is fitted value for now; kAbl = 0.01 is driving the equilibrium towards uptake into the proximal tubule cells
    kAap = kAap, # affinity constant apical this is about OAT4 which has affinity to re-abs (movement from filtrate to cells; this is fitted value for now; kAap = 0.01 is driving the equilibrium towards re-absorption into the proximal tubule cells
    fup = fup, # Unbound fraction in plasmal Fischer et al. 2024 https://doi.org/10.1021/acs.est.3c07415;
    PL = PL,  # Plasma/liver partition coefficient; Rat tissue data (Kudo et al. 2007)
    PF = PF, # Plasma/fat partition coefficient; Rat tissue data (Kudo et al. 2007)
    PK = PK,  # Plasma/kidney partition coefficient; Rat tissue data (Kudo et al. 2007)
    PSk = PSk,  # Plasma/skin partition coefficient; Rat tissue data (Kudo et al. 2007)
    PR = PR,  # Plasma/rest of the body partition coefficient; Rat tissue data (Kudo et al. 2007)
    PG = PG,  # Plasma/gut partition coefficient; Rat tissue data (Kudo et al. 2007)
    KpAd = KpAd, # Calculated
    KpGu = KpGu, # Calculated
    KpKi = KpKi, # Calculated
    KpLi = KpLi, # Calculated
    KpSk = KpSk, # Calculated
    KpRe = KpRe, # Calculated
    Papp = Papp, # Calculated
    CL_OAT1 = CL_OAT1, 
    REF_OAT1 = REF_OAT1,
    
  ) %>%
  
  # Age-dependent variables
  mutate(
    SA_hands_M = SA_hands_M, #Surface area of hands
  )
# # Dermal absorption
# mutate(CLdermalabs_M = ((Papp*SA_hands_M)/1000)*24) %>% # (L/d); cm/h*cm^2 = mL/h /1000 = L/h * 24 = L/d
# 
# # Kidney clearance
# mutate(# CL_PltPT_M = ((CL_OAT1*REF_OAT1) + (CL_OAT3*REF_OAT3)) * PTC_kidneyTissue_M, # L/d plasma to proximal tubule clearance
#        CL_PltPT_M = ((CL_OAT1*REF_OAT1) + (CL_OAT3*REF_OAT3)) * 0.17 * 0.7 * V_kidney_M,# L/d scaling for protein content instead of cell content; 17% of kidney is protein and 70% of kidney is cortex [ICRP 89], assuming that all the kidney protein is found in the cortex; this is an overestimation though; double ref for 17% protein Ruark 2020: DOI: https://doi.org/10.1016/B978-0-12-818596-4.00006-0
#        # CL_FiltPT_M = (CL_OAT4*REF_OAT4) * PTC_kidneyTissue_M, #L/d filtrate to proximal tubule clearance
#        CL_FiltPT_M = CL_OAT4 * REF_OAT4 * 0.17 * 0.7 * V_kidney_M # L/d scaling for protein content instead of cell content; 17% of kidney is protein and 70% of kidney is cortex [ICRP 89], assuming that all the kidney protein is found in the cortex; this is an overestimation though; double ref for 17% protein Ruark 2020: DOI: https://doi.org/10.1016/B978-0-12-818596-4.00006-0
#        # Trine's values
#        # CLurine_M = CLurinec*BDW_M^(-0.25) #from Husoy; L/d clearance urine
# ) %>%
# 
# # Biliary clearance
# mutate(CLbiliary_M = CLbiliaryc*(BDW_M^0.1)) %>% #from Husoy; L/d biliary clearance rate
# 
# # Fecal clearance
# mutate(CLfecal_M = CLfaecesc*(BDW_M^0.001)) #from Husoy; L/d faeces clearance rate

write.csv(Final_variables_M_df, "Final_variables_M.csv", row.names = FALSE)

# Female
Final_variables_F_df <- PhysioVariables_F_df %>%
  
  # Independent variables
  mutate(
    AbsPFOA = AbsPFOA, # fraction absorbed, from Trine
    kAbl = kAbl, # affinity constant basolateral this is about OAT1 and OAT3 which have affinity to uptake (movement from plasma to cells; this is fitted value for now; kAbl = 0.01 is driving the equilibrium towards uptake into the proximal tubule cells
    kAap = kAap, # affinity constant apical this is about OAT4 which has affinity to re-abs (movement from filtrate to cells; this is fitted value for now; kAap = 0.01 is driving the equilibrium towards re-absorption into the proximal tubule cells
    fup = fup, # Unbound fraction in plasmal Fischer et al. 2024 https://doi.org/10.1021/acs.est.3c07415;
    PL = PL,  # Plasma/liver partition coefficient; Rat tissue data (Kudo et al. 2007)
    PF = PF, # Plasma/fat partition coefficient; Rat tissue data (Kudo et al. 2007)
    PK = PK,  # Plasma/kidney partition coefficient; Rat tissue data (Kudo et al. 2007)
    PSk = PSk,  # Plasma/skin partition coefficient; Rat tissue data (Kudo et al. 2007)
    PR = PR,  # Plasma/rest of the body partition coefficient; Rat tissue data (Kudo et al. 2007)
    PG = PG,  # Plasma/gut partition coefficient; Rat tissue data (Kudo et al. 2007)
    KpAd = KpAd, # Calculated
    KpGu = KpGu, # Calculated
    KpKi = KpKi, # Calculated
    KpLi = KpLi, # Calculated
    KpSk = KpSk, # Calculated
    KpRe = KpRe
  ) %>%
  
  # Dermal absorption
  mutate(CLdermalabs_F = ((Papp*SA_hands_F)/1000)*24) %>% # (L/d); cm/h*cm^2 = mL/h /1000 = L/h * 24 = L/d
  
  # Kidney clearance
  mutate(# CL_PltPT_F = ((CL_OAT1*REF_OAT1) + (CL_OAT3*REF_OAT3)) * PTC_kidneyTissue_F, # L/d plasma to proximal tubule clearance
    CL_PltPT_F = ((CL_OAT1*REF_OAT1) + (CL_OAT3*REF_OAT3)) * 0.17 * 0.7 * V_kidney_F,# L/d scaling for protein content instead of cell content; 17% of kidney is protein and 70% of kidney is cortex [ICRP 89], assuming that all the kidney protein is found in the cortex; this is an overestimation though; double ref for 17% protein Ruark 2020: DOI: https://doi.org/10.1016/B978-0-12-818596-4.00006-0
    # CL_FiltPT_F = (CL_OAT4*REF_OAT4) * PTC_kidneyTissue_F, #L/d filtrate to proximal tubule clearance
    CL_FiltPT_F = CL_OAT4 * REF_OAT4 * 0.17 * 0.7 * V_kidney_F # L/d scaling for protein content instead of cell content; 17% of kidney is protein and 70% of kidney is cortex [ICRP 89], assuming that all the kidney protein is found in the cortex; this is an overestimation though; double ref for 17% protein Ruark 2020: DOI: https://doi.org/10.1016/B978-0-12-818596-4.00006-0
    # Trine's values
    # CLurine_F = CLurinec*BDW_F^(-0.25) #from Husoy; L/d clearance urine
  ) %>%
  
  # Biliary clearance
  mutate(CLbiliary_F = CLbiliaryc*(BDW_F^0.1)) %>% #from Husoy; L/d biliary clearance rate
  
  # Fecal clearance
  mutate(CLfecal_F = CLfaecesc*(BDW_F^0.001)) #from Husoy; L/d faeces clearance rate

write.csv(Final_variables_F_df, "Final_variables_F.csv", row.names = FALSE)


