# --------------------------------------------------------------------------- #
# SCRIPT FOR CALCULATING PHYSIOLOGICAL CONSTANTS 
# Note: Actually not needed to run the model as exported csv file is used as input
# By: Chrysanthi Pachoulide, Joost Westernhout, Jack Koster initial ref. A. Ratier et al. 2024 doi: 10.1016/j.envint.2024.108621
# Date: 14-04-2025
# --------------------------------------------------------------------------- #

rm(list=ls()) # to clear out the global environment

# Packages
library(here)
library(tidyverse)
library(showtext)
font_add(family = "Garamond", regular = "GARA.TTF")
showtext_auto()

# Set storage directory
INPUT <- here("Input")


# LIFETIME EQUATIONS ####
# ------------------------------------------------------ #

lifeTSTOP = 80 # duration of lifetime (0 - 80 years old)  of simulation
TSTART = 0
TSTOP = 365*lifeTSTOP # years in days
DT = 1
TIME = seq(TSTART,TSTOP,by=DT)


# Creating a dataframe for all the variables
Variables_df = as.data.frame(list(TIME = TIME)) # df column 1 = simulation time, every step is 1 day
Variables_df = Variables_df %>%
  mutate(age = TIME/365) # add column 2 = age in days

Variables_df = Variables_df %>%
  
## Base Physiology ####

# Body weight (kg), reference: Deepika et al. 2021, https://doi.org/10.1016/j.envres.2021.111287
# mutate(BW_M=3.382e+00+
#          2.866e+00*age+
#          (1.694e-01)*age^2-(1.169e-02)*age^3+
#          (2.577e-04)*age^4-(2.484e-06)*age^5+
#          (8.891e-09)*age^6,
#        BW_F=2.354+4.050*age +
#          -(3.240e-02)*age^2 +
#          -(3.057e-03)*age^3 +
#          (9.353e-05)*age^4 +
#          -(1.022e-06)*age^5 +
#          (3.918e-09)*age^6  ) %>% 
# BW_M_Ratier_2024 & BW_F_Ratier_2024 = Equation extracted from supplemental material from Ratier et al., 2024
mutate(BW_M_Ratier_2024 = if_else(age <19.00093277, 74.16235828-(2*(74.16235828-57.19957758)/(exp(0.63466182*(age-13.31018000))+exp(0.05457656*(age-13.31018000)))),
                                  -0.01129273*age^2 + 1.11817056*age + 56.74397436)) %>%
  mutate(BW_F_Ratier_2024 = if_else(age <17.9374115, 62.95490567-(2*(62.95490567-49.36574299)/(exp(0.84039606*(age-11.56691488))+exp(0.06710088*(age-11.56691488)))),
                                    -0.01258006*age^2 + 1.25029379*age + 44.4459234)) %>%
  mutate(BDW_M_Ratier_2024 = 74.16235828-(2*(74.16235828-57.19957758)/(exp(0.63466182*(age-13.31018000))+exp(0.05457656*(age-13.31018000))))) %>% 
  mutate(BDW_F_Ratier_2024 = 62.95490567-(2*(62.95490567-49.36574299)/(exp(0.84039606*(age-11.56691488))+exp(0.06710088*(age-11.56691488))))) 

# Body height (cm2), reference: Deepika et al. 2021, https://doi.org/10.1016/j.envres.2021.111287
Variables_df = Variables_df %>% 
  mutate(BH_M = (5.869e+01)+(1.265e+01)*age-(4.665e-01)*age^2+(7.198e-03)*age^3-(3.224e-05)*age^4-(2.512e-07)*age^5+(2.071e-09)*age^6,
         BH_F = (5.373e+01)+(1.296e+01)*age-(5.506e-01)*age^2+(1.113e-02)*age^3-(1.106e-04)*age^4+(4.697e-07)*age^5-(4.416e-10)*age^6) 

# Body Surface Area (m2), reference: Gastellu et al. 2024, 10.1016/j.envres.2024.120393 (supplementary file Physio_equations_detailed.xlsx, eq. from Pendse et al. 2020)
Variables_df = Variables_df %>%  
  mutate(BSA_M = exp(-3.75 + 0.42*log(BH_M)+0.52*log(BW_M_Ratier_2024)),
         BSA_F =  exp(-3.75 + 0.42*log(BH_F)+0.52*log(BW_F_Ratier_2024)))  


  
# Blood/Plasma/Hematocrit/Creatinine
# Using Ratier et al. (2024) model, Fraction of arterial plasma, calculated from Filser 2000 p.43
Fr_art_blood = 0.0178 / (0.0178 + 0.0533) #fraction of arterial blood 

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

# Default serum creatinine values
SCr_M = 0.90
SCr_F = 0.70


## Fractional Volumes ####

Variables_df = Variables_df %>%
  # select(TIME,age,BW_M_Ratier_2024,BW_F_Ratier_2024,BDW_M_Ratier_2024,BDW_F_Ratier_2024) %>%
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
  mutate(V_brainFraction_M = (1.450 + (0.353 - 1.450) * exp (-0.440*age))/BW_M) %>%
  mutate(V_brainFraction_F = (1.300 + (0.347 - 1.300) * exp (-0.573*age))/BW_F) %>%
  
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
  
  # Blood volume; compartment [22] in Ratier 2024 
  mutate(V_bloodFraction_M = if_else(age < 1, (-0.0273*age + 0.0771),
                                      0.0761 + (0.0289 - 0.0761)*exp(-0.592*age))) %>%
  mutate(V_bloodFraction_F = if_else(age < 1, (-0.0273*age + 0.0771),
                                      if_else(age < 14.019723, 3.28E-5*age^3 - 1.21E-3*age^2 + 1.24E-2*age + 3.86E-2,
                                              0.065))) %>%
  
  # Adipose tissue
  mutate(V_adiposeFraction_M = 0.96 - V_adrenalFraction_M - V_boneFraction_M - V_bonenonperfusedFraction_M - V_brainFraction_M - V_breastFraction_M +
           - V_heartFraction_M - V_marrowFraction_M - V_muscleFraction_M - V_reproFraction_M - V_pancreasFraction_M +
           - V_skinFraction_M - V_spleenFraction_M - V_thyroidFraction_M - V_urinarytractFraction_M - V_kidneyFraction_M +
           - V_lungFraction_M - V_gutFraction_M - V_stomachFraction_M - V_liverFraction_M - V_bloodFraction_M) %>%
  mutate(AdiposeMass_M = if_else(age < 19.00093277, 0,
                                 (-0.01129273*age^2 + 1.11817056*age + 56.74397436) - (74.16235828-(2*(74.16235828-57.19957758)/(exp(0.63466182*(age-13.31018000))+exp(0.05457656*(age-13.31018000))))))) %>% # age and not BW dependent
  mutate(V_adiposeFraction_F = 0.96 - V_adrenalFraction_F - V_boneFraction_F - V_bonenonperfusedFraction_F - V_brainFraction_F - V_breastFraction_F +
           - V_heartFraction_F - V_marrowFraction_F - V_muscleFraction_F - V_reproFraction_F - V_pancreasFraction_F +
           - V_skinFraction_F - V_spleenFraction_F - V_thyroidFraction_F - V_urinarytractFraction_F - V_kidneyFraction_F +
           - V_lungFraction_F - V_gutFraction_F - V_stomachFraction_F - V_liverFraction_F - V_bloodFraction_F) %>%
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
# write.csv(Variables_df, "PhysioVariables.csv", row.names = FALSE)


## Glomerular Filtration Rate (L/day) : 

## # Baseline neonatal GFR should be 20.0 mL/min: Qi = (20*0.9)/107.3 = 0.1678 mg/dL (Smeets 2022, https://doi.org/10.1681/ASN.2021101326)
Variables_df <- Variables_df %>% 
  # Calculate median SCr value for age/sex
  mutate(
    SCrQ_M = if_else(age < 18, 0.1678 + ((0.90  - 0.1678) / 18) * age,
                     SCr_M),
    SCrQ_F = if_else(age < 18, 0.1305 + ((0.70  - 0.1305) / 18) * age,
                     SCr_F)) %>% 
  # Baseline GFR for males and females (in L/day)
  # (mL/min/1.73m^2 -> L/day)  # scale to actual BSA: SA_B*1e-4 / 1.73
  mutate( 
    GFR_M_base = (107.3 * 1.44*(BSA_M)/1.73) / (SCr_M/SCrQ_M), # Pottel et al. 2016 doi: 10.1093/ndt/gfv454
    GFR_F_base = (107.3 * 1.44*(BSA_F)/1.73) / (SCr_F/SCrQ_F)  # Pottel et al. 2016 doi: 10.1093/ndt/gfv454
  ) %>% 
  # Exponential decline after age 40
  mutate(
    GFR_M = if_else(age <= 40, GFR_M_base,
                    GFR_M_base * 0.988^(age - 40)),
    GFR_F = if_else(age <= 40, GFR_F_base,
                    GFR_F_base * 0.988^(age - 40))
  ) %>% 
  # mutate(
  #   GFR_M = Q_GFRc_M*0.6944444444*1.73,
  #   GFR_F = Q_GFRc_F*0.6944444444*1.73
  # )
  # Serum Albumin per age: (g_HSA /L_serum, or mg/ml);
  # Derived from Weaving et al 2016; DOI: 10.1177/0004563215593561, Figure 1
  mutate( 
    SAlb_M = 3.9330e+01 + (5.8156e-01*age) - (1.9072e-02*(age^2)) +
      (2.3220e-04*(age^3)) - (1.0313e-06*(age^4)), 
    SAlb_F = 4.0540e+01 + (6.0605e-01*age) - (3.5316e-02*(age^2)) +
      7.9681e-04*(age^3) - (7.8062e-06*(age^4)) + (2.7227e-08*(age^5))
  )

# write.csv(Variables_df, here("Input", "PhysioVariables.csv"), row.names = FALSE)


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
#   geom_path(data = Variables_M_df, aes(age, BW_M, color = "Male")) +
#   geom_path(data = Variables_F_df, aes(age, BW_F, color = "Female")) +
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

## Plots ####

CP_theme <- theme_minimal() +
  theme(
    axis.title = element_text(size = 44, face = "bold"),
    axis.text = element_text(size = 40),
    plot.title = element_text(hjust = 0, size = 48, face = "bold"),
    plot.margin = margin(0.3, 0.3, 0.3, 0.3, "cm"), 
    axis.line = element_line(size = 0.1),
    panel.grid = element_line(size = 0.3),
    plot.caption = element_text(hjust = 0, size = 40, face= "italic"),
    plot.title.position = "plot", 
    plot.caption.position = "plot",
    legend.position = "none"
  )  

Flows <- Variables_df %>%
  mutate(FGFR_M = Q_kidneyFraction_M*CardOut_M*0.18,
         FGFR_F = Q_kidneyFraction_F*CardOut_F*0.18) %>%
  mutate(QK_M = Q_kidneyFraction_M*CardOut_M,
         QK_F = Q_kidneyFraction_M*CardOut_F) %>%
  mutate(QA_M = Q_adiposeFraction_M*CardOut_M,
         QA_F = Q_adiposeFraction_F*CardOut_F,
         QL_M = Q_liverFraction_M*CardOut_M,
         QL_F = Q_liverFraction_F*CardOut_F,
         QG_M = Q_gutFraction_M*CardOut_M,
         QG_F = Q_gutFraction_F*CardOut_F)

Flows %>%
  ggplot()+
  geom_path(aes(age, QL_M, color = "Liver", linetype = "Male")) +
  geom_path(aes(age, QL_F, color = "Liver", linetype = "Female")) +
  geom_path(aes(age, QK_M, color = "Kidney", linetype = "Male")) +
  geom_path(aes(age, QK_F, color = "Kidney", linetype = "Female")) +
  geom_path(aes(age, QA_M, color = "Adipose", linetype = "Male")) +
  geom_path(aes(age, QA_F, color = "Adipose", linetype = "Female")) +
  geom_path(aes(age, QG_M, color = "Gut", linetype = "Male")) +
  geom_path(aes(age, QG_F, color = "Gut", linetype = "Female")) +
  scale_color_manual(values = c("Liver" = "goldenrod2",
                                "Kidney" = "royalblue2",
                                "Adipose" = "seagreen2",
                                "Gut" = "turquoise2",
                                "Skin" = "deeppink2"),
                     name = "") +
  scale_linetype_manual(values = c("Male" = "dashed",
                                   "Female" = "solid"),
                        name = "") +
  CP_theme +
  ylab("Blood flow (L/d)") +
  xlab("Age (years)")
  # theme(axis.title = element_text(size = 15),
  #       axis.text = element_text(size = 14),
  #       legend.position = "bottom",
  #       legend.text = element_text(size = 14))
ggsave(filename = here("OrganFlows.png"),
       dpi = 600,
       width = 10,
       height = 8,
       units = "cm")


ggplot() +
  # geom_path(data = Variables_df, aes(age, GFR_M), linetype = "dashed", color = "blue") +
  # geom_path(data = Variables_df, aes(age, GFR_F), linetype = "solid", color = "blue") + 
  # geom_path(data = Flows, aes(age, FGFR_M), linetype = "dashed", color = "grey") + 
  # geom_path(data = Flows, aes(age, FGFR_F), linetype = "solid", color = "grey") + 
  geom_path(data = Variables_df, aes(age, GFR_M, linetype = "Male", color = "FAS_equation")) +
  geom_path(data = Variables_df, aes(age, GFR_F, linetype = "Female", color = "FAS_equation")) +
  geom_path(data = Flows, aes(age, FGFR_M, linetype = "Male", color = "Fraction_filtrated")) +
  geom_path(data = Flows, aes(age, FGFR_F, linetype = "Female", color = "Fraction_filtrated")) +
  # Manually set colors and line types
  scale_color_manual(
    name = "Method",
    values = c("FAS_equation" = "blue", "Fraction_filtrated" = "grey"),
    labels = c("FAS_equation", "Fraction_filtrated")
  ) +
  scale_linetype_manual(
    name = "Sex",
    values = c("Male" = "dashed", "Female" = "solid"),
    labels = c("Male", "Female")
  ) +
  # Combine legends
  guides(
    color = guide_legend(title = "Method"),
    linetype = guide_legend(title = "Sex")
  ) +
  CP_theme +
  ylab("GFR (ml/min)") +
  xlab("Age (years)") +
  theme(legend.position = c(1, 1),       
        legend.justification = c(1, 1), 
        legend.text = element_text(size = 40),
        legend.title = element_blank(),
        legend.key.size = unit(0.5, "lines"),
        legend.spacing.y = unit(0, "cm"),
        legend.margin = margin(0, 0, 0, 0))

ggsave(filename = here("GFR.png"),
       dpi = 600,
       width = 9,
       height = 9,
       units = "cm")

ggplot() +
  geom_path(data = Variables_df, aes(age, BW_M, linetype = "Male")) +
  geom_path(data = Variables_df, aes(age, BW_F, linetype = "Female")) +
  scale_linetype_manual(values = c("Male" = "dashed",
                                   "Female" = "solid"),
                        name = "") +
  CP_theme +
  ylab("BW (kg)") +
  xlab("Age (years)")+
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 14),
        legend.position = "bottom",
        legend.text = element_text(size = 14))
ggsave(filename = here("BW.png"),
       dpi = 600,
       width = 10,
       height = 8,
       units = "cm")

ggplot() +
  geom_path(data = Variables_df, aes(age, CardOut_M, linetype = "Male")) +
  geom_path(data = Variables_df, aes(age, CardOut_F, linetype = "Female")) +
  scale_linetype_manual(values = c("Male" = "dashed",
                                   "Female" = "solid"),
                        name = "") +
  CP_theme+
  ylab("QC (L/d)") +
  xlab("Age (years)")
  # theme(axis.title = element_text(size = 15),
  #       axis.text = element_text(size = 14),
  #       legend.position = "bottom",
  #       legend.text = element_text(size = 14))
ggsave(filename = here("QC.png"),
       dpi = 600,
       width = 10,
       height = 8,
       units = "cm")

Volumes <-Variables_df %>%
  mutate(VK_M = V_kidneyFraction_M*BW_M,
         VK_F = V_kidneyFraction_M*BW_F,
         VA_M = V_adiposeFraction_M*BW_M,
         VA_F = V_adiposeFraction_F*BW_F,
         VL_M = V_liverFraction_M*BW_M,
         VL_F = V_liverFraction_F*BW_F,
         VG_M = V_gutFraction_M*BW_M,
         VG_F = V_gutFraction_F*BW_F)

Volumes %>%
  ggplot()+
  geom_path(aes(age, VL_M, color = "Liver", linetype = "Male")) +
  geom_path(aes(age, VL_F, color = "Liver", linetype = "Female")) +
  geom_path(aes(age, VK_M, color = "Kidney", linetype = "Male")) +
  geom_path(aes(age, VK_F, color = "Kidney", linetype = "Female")) +
  geom_path(aes(age, VA_M, color = "Adipose", linetype = "Male")) +
  geom_path(aes(age, VA_F, color = "Adipose", linetype = "Female")) +
  geom_path(aes(age, VG_M, color = "Gut", linetype = "Male")) +
  geom_path(aes(age, VG_F, color = "Gut", linetype = "Female")) +
  scale_color_manual(values = c("Liver" = "goldenrod2",
                                "Kidney" = "royalblue2",
                                "Adipose" = "seagreen2",
                                "Gut" = "turquoise2",
                                "Skin" = "deeppink2"),
                     name = "") +
    scale_linetype_manual(values = c("Male" = "dashed",
                                   "Female" = "solid"),
                        name = "") +
  CP_theme +
  ylab("Weight (Kg)") +
  xlab("Age (years)")
  # theme(axis.title = element_text(size = 15),
  #       axis.text = element_text(size = 14),
  #       legend.position = "bottom",
  #       legend.text = element_text(size = 14))
ggsave(filename = here("OrganVolumes.png"),
       dpi = 600,
       width = 10,
       height = 8,
       units = "cm")
# 
# CalcPhysioParams <- cbind(Flows, Volumes)
# write.csv(CalcPhysioParams, here("Input", "CalculatedPhysiologicalParams.csv"), row.names = FALSE)
