# *************************************************************************** #
# PBK MODEL FOR PFOA, TO BE USED TOGETHER WITH THE LATEST HBM DATA
# The model is an update from the Husoy 2023 model
# Lifestage equations were taken from Ratierau and updated by JW 
# Code updates by CP, AN and JW
# *************************************************************************** #

rm(list=ls()) # to clear out the global environment

# Set working directory

HOME <- "C:/Users/pacho003/OneDrive - Wageningen University & Research/C Channel/R/PARC_PFAS_PBPKmodel/Codes_PFAS_models/Chrysanthi_Pachoulide_2024_PFOA_oral_dermal_blood"
# HOME <- "/home/westerj"
setwd(HOME)


# Set output storage directory
 
workingtime <- gsub(":", "-", Sys.time())
OUTPUT <- file.path("Output", Sys.Date(), workingtime)
dir.create(OUTPUT, recursive = TRUE)
setwd(OUTPUT)


# Load packages

library(lubridate)
library(ggplot2)
library(deSolve)
library(writexl)
library(ggpubr)
library(tidyverse)


# PBK MODEL PARAMETERS ####
# ---------------------------------------------------------------------------- #


## EXPOSURE SCENARIO ####
## ------------------------------------------------------ #

exposure_stop <- 50
sim_stop <- 80

Oraldose <- 0.000187 # ug/kg/day [EFSA 2020; page 143]
Dermconc <- 0.000542 # ug/kg/day; mean of #as.numeric(SumExpPFOA_LB_val[i,14])
skin_fraction <- 0.05 # fraction of the total skin surface exposed; hands are 5% https://www.epa.gov/sites/default/files/2015-09/documents/efh-chapter07.pdf
# Note Chrysa 18-10-2024: the mean % of total surface area that is hands is 5.2 (in adult male 21+ years) and 4.8 (in adult female), so that would be 5 +/- 0.2 %
# Note Chrysa 24-20-2024: assumption that the main skin area exposed is the hands

# For future use; to enable input for dermal exposure from different cosmetics/ aggregate exposure scenarios etc
# Eproduct = # mg/kg bw/d Effective dermal exposure depending to a product category, this takes into account the rinsing off, dilution etc; https://health.ec.europa.eu/system/files/2022-08/sccs_o_250.pdf
# Dermconc = # Compound concentration in the specific product; https://health.ec.europa.eu/system/files/2022-08/sccs_o_250.pdf
# Edermal = Eproduct * Dermconc # Dermal exposure, based on EFSA; https://health.ec.europa.eu/system/files/2022-08/sccs_o_250.pdf

# From Ragnarsdottir et al.
# DED = (Dermconc * BSA * DAS * Fa * IEF)/BW # daily exposure dose (ng/kg bw/day)
# Dermconc # (ng/g) concentration in dust
# BSA # body surface area exposed (cm^2)
# DAS # (mg/cm^2) dust adhered to skin
# Fa # fraction absorbed by the skin; 50% of fBAc
# IEF # indoor exposure fraction



## CHEMICAL SPECIFIC PARAMETERS ####
## ------------------------------------------------------ #

# PFOA 
MW = 414.07 #PFOA MW
fup <- 0.061 # from Fischer et al. 2024 https://doi.org/10.1021/acs.est.3c07415; OLD fup = 0.02 # fup fraction of PFOA in plasma

# Note Chrysa 24-20-2024: plasma/tissue partition coefficients to be changed, using the Allendorf paper:  https://doi.org/10.1002/etc.4954
PL <- 2.2  # Plasma/liver partition coefficient; Rat tissue data (Kudo et al., 2007)
PF <- 0.04  # Plasma/fat partition coefficient; Rat tissue data (Kudo et al., 2007)
PK <- 1.05  # Plasma/kidney partition coefficient; Rat tissue data (Kudo et al., 2007)
PSk <- 0.1  # Plasma/skin partition coefficient; Rat tissue data (Kudo et al., 2007)
PR <- 0.12  # Plasma/rest of the body partition coefficient; Rat tissue data (Kudo et al., 2007)
PG <- 0.05  # Plasma/gut partition coefficient; Rat tissue data (Kudo et al., 2007)

# Comment Chrysa 21-10-2024: AbsPFOA is used in the Dermaldose input, so we're already correcting before using the Papp?! not sure I agree with this
#fBAc <- # Fraction bio accessible; fraction of the compound released from the matrix (cosmetic formulation, dust etc) and is available to be absorbed from the epidermis
AbsPFOA <- 0.016 # 0.00048 # Changed to the absorption measured by Abraham and Monien 2022, of 1.6% of applied dose from sunscreen. 
Papp = 3.82 * 10^-3 # cm/h ref: https://doi.org/10.1016/j.envint.2024.108772

# Renal clearance parameters
# # CLurinec = 0.000044  # L/d/kg; 0.044 mL/d/kg taken from Fujii et al 2015 
# Vmaxc = 4.5*MW/1000*60*24 # ug/d/mg protein; 45 nmol/min/mg protein *MW/1000*60*24 = ug/d/mg protein ref: Louisse et al. 2023 https://doi.org/10.1007/s00204-022-03428-6
# Km = 47*MW # ug/L; 47 uM*MW = ug/L ref: Louisse et al. 2023 https://doi.org/10.1007/s00204-022-03428-
# Comment Chrysa 05-11-2024: 
## Assuming that the kidney PFAS concentrations never reaches Km concentrations, therefore transforming Vmax and Km to a Clearance; in the paper what they call transporter efficiency : Louisse, Pedroni et al. 2024 https://doi.org/10.1016/j.tox.2024.153961)
## OAT1 and OAT3 are determining transport between blood and proximal tubule
## OAT4 is determining transport between proximal tubule and filtrate
## All transporters are bi-directional, therefore the equation is written assuming that direction is determined based on the equilibrium between the two compartments
## Introducing the "affinity constants: kAap and kAbl, to compensate for the fact that the transporters have an affinity to one side.
kAbl <- 0.01 # affinity constant basolateral, this is about OAT1 and OAT3 which have affinity to uptake (movement from plasma to cells); this is fitted value for now; kAbl = 0.01 is driving the equilibrium towards uptake into the proximal tubule cells
kAap <- 0.01 # affinity constant apical, this is about OAT4 which has affinity to re-abs (movement from filtrate to cells); this is fitted value for now; kAap = 0.01 is driving the equilibrium towards re-absorption into the proximal tubule cells
CL_OAT1 <- 19 * 10^-6 *60*24 * 10^-6 #L/d/kg protein protein; initial ul/min/mg protein
CL_OAT3 <- 17* 10^-6 *60*24 * 10^-6 #L/d/kg protein; initial ul/min/mg protein
CL_OAT4 <- 96* 10^-6 *60*24 * 10^-6 #L/d/kg protein; initial ul/min/mg protein
PTCPGK <- 9.94* 10^7 * 10^3 # proximal tubule cells/kg kidney cortexl initial PTC/g kidney https://doi.org/10.1021/acs.molpharmaceut.4c00504
REF_OAT1 <- 4.3/26.6 # relative expression factor: expression in the human kidneys /expression in the cells (4.3 ± 0.3 pmol/mg membrane protein in the human kidney cortex and 26.6 ± 3.4 pmol/mg membrane protein: OAT1 expression in HEK293-OAT1 cells  https://doi.org/10.1124/dmd.121.000367); alternative:  5.33 ± 1.88 pmol/mg protein in the human kidney cortex http://dx.doi.org/10.1124/dmd.116.072066
REF_OAT3 <- 2.7/7.3 # relative expression factor: expression in the human kidneys /expression in the cells (2.7 ± 0.1 pmol/mg membrane protein in the human kidney cortex and 7.3 ± 0.5 pmol/mg membrane protein: OAT3 expression in HEK293-OAT3 cells  https://doi.org/10.1124/dmd.121.000367);  alternative: 3.50 ± 1.55 pmol/mg membrane protein in the human kidney cortex http://dx.doi.org/10.1124/dmd.116.072066
REF_OAT4 <- 0.52/16 # relative expression factor: expression in the human kidneys /expression in the cells  (0.52 ± 0.23 pmol/mg protein in the human kidney cortex http://dx.doi.org/10.1124/dmd.116.072066; OAT4 expression in HEK293-OAT4 cells not found therefore mean of OAT1 and OAT3 used) for alternative input https://doi.org/10.1002/cpt.2396) 
Papp_kidney <- 1.46 * 10^-6 #cm/s https://doi.org/10.1016/j.chemosphere.2024.142390
PT_Sarea <- 0.0176*2 #cm^2; multiplying by 2 as we have 2 kidneys DOI 10.1007/978-1-4614-3785-7, (14 mm long × 40 µm thick) => 1.76 mm^2

# Hepatic clearance parameters
CLbiliaryc = 0.00262 # L/d/kg ; 2.62 +/- 3.6 mL/d/kg from Fujii et al 2015 DOI: 10.1539/joh.14-0136-OA
CLfaecesc = 0.000052 # L/d/kg ; 0.052 +/- 0.05 mL/d/kg clearance in faeces taken from Fujii et al 2015 DOI: 10.1539/joh.14-0136-OA


## PHYSIOLOGICAL PARAMETERS ####
## ------------------------------------------------------ #

# Kidney scaling factors
# Kcells = 6E7	# number of cells per gram kidney from DOI 10.1007/978-1-4614-8229-1_7, chapter 7; they say it should be included in the sensitivity analysis
# Kprotein = 2.0e-9	# gram protein/proximal tubule cell
# SFOAT4 = 15 # for now it's 1; could be 1 pmole/g tissue from https://doi.org/10.1124/dmd.119.086579, but this is used for scaling between animals; as they mention in the article: "the data obtained from the absolute peptide approach should not be considered as absolute molar protein abundance data because complete trypsin digestion may not be confirmed"



### Lifetime equations ####

#### Duration of lifetime (0 - 80 years old) ----
TSTART <- 0
TSTOP <- 365*sim_stop # years in days
DT <- 1
TIME <- seq(TSTART,TSTOP,by=DT)

Variables_df <- as.data.frame(list(TIME = TIME)) #df column 1 = simulation time, every step is 1 day
Variables_df <- Variables_df %>%
  mutate(age = TIME/365) # add column 2 = age in days

##### Body weight ----
Variables_df <- Variables_df %>%
  # BW_M_Ratier_2024 & BW_F_Ratier_2024 = Equation extracted from supplemental material from Ratier et al., 2024
  mutate(BW_M_Ratier_2024 = if_else(age <19.00093277, 74.16235828-(2*(74.16235828-57.19957758)/(exp(0.63466182*(age-13.31018000))+exp(0.05457656*(age-13.31018000)))),
                                    -0.01129273*age^2 + 1.11817056*age + 56.74397436)) %>%
  mutate(BW_F_Ratier_2024 = if_else(age <17.9374115, 62.95490567-(2*(62.95490567-49.36574299)/(exp(0.84039606*(age-11.56691488))+exp(0.06710088*(age-11.56691488)))),
                                    -0.01258006*age^2 + 1.25029379*age + 44.4459234)) %>%
  mutate(BDW_M_Ratier_2024 = 74.16235828-(2*(74.16235828-57.19957758)/(exp(0.63466182*(age-13.31018000))+exp(0.05457656*(age-13.31018000))))) %>%
  mutate(BDW_F_Ratier_2024 = 62.95490567-(2*(62.95490567-49.36574299)/(exp(0.84039606*(age-11.56691488))+exp(0.06710088*(age-11.56691488)))))

##### Plots ----
# Plot BW changes over time
# Comment Chrysa on 18-10-2024: MassBalance issue: should we then use the BDW term that is also changing during adulthood or are we OK with the massbalance issue in the total body volume?
# PlotBDW <- Variables_df %>% select(c(age, BDW_M_Ratier_2024, BDW_F_Ratier_2024)) %>%
#   rename(Male = BDW_M_Ratier_2024, Female = BDW_F_Ratier_2024) %>%
#   pivot_longer(names_to = "Gender", values_to = "BDW", Male:Female)
# 
# PlotBW <- Variables_df %>% select(c(age, BW_M_Ratier_2024, BW_F_Ratier_2024)) %>%
#   rename(Male = BW_M_Ratier_2024, Female = BW_F_Ratier_2024) %>%
#   pivot_longer(names_to = "Gender", values_to = "BW", Male:Female)
# ggplot() +
#   geom_path(data = PlotBDW, aes(age, BDW, colour = Gender, linetype = "BDW")) +
#   geom_path(data = PlotBW, aes(age, BW, colour = Gender, linetype = "BW")) +
#   labs(x = "Age", 
#        y = "Values",
#        colour = "Gender", 
#        linetype = "Bodyweight method") +
#   scale_linetype_manual(values = c("BDW" = "dashed", "BW" = "solid"))
#ggsave("BWovertime.png", dpi = 300)

##### Blood/Plasma/Hematocrit ----
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

##### Storage data frame ----                                            
Variables_df <- Variables_df %>%
  select(TIME,age,BW_M_Ratier_2024,BW_F_Ratier_2024,BDW_M_Ratier_2024,BDW_F_Ratier_2024) %>%
  rename(BW_M = BW_M_Ratier_2024) %>% #could actually be ignored as we only use BDW and not BW
  rename(BW_F = BW_F_Ratier_2024) %>% #could actually be ignored as we only use BDW and not BW
  rename(BDW_M = BDW_M_Ratier_2024) %>%
  rename(BDW_F = BDW_F_Ratier_2024) %>%
  
  mutate(Oraldose_M = Oraldose*BW_M) %>%
  mutate(Oraldose_F = Oraldose*BW_F) %>%
  mutate(Oraldose_M = if_else(age <= exposure_stop,Oraldose_M,0)) %>%
  mutate(Oraldose_F = if_else(age <= exposure_stop,Oraldose_F,0)) %>%
  
  mutate(Dermaldose_M = AbsPFOA*Dermconc*BW_M) %>% # dermal concentration is expressed as ug/kg BW/day
  mutate(Dermaldose_F = AbsPFOA*Dermconc*BW_F) %>% # dermal concentration is expressed as ug/kg BW/day
  mutate(Dermaldose_M = if_else(age <= exposure_stop,Dermaldose_M,0)) %>%
  mutate(Dermaldose_F = if_else(age <= exposure_stop,Dermaldose_F,0)) %>%
  
  ## Hematocrit
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
  
  
  ## Cardiac output (plasma; L/min*60*24 = L/d)
  mutate(CardOut_M = if_else(age < 33.37, (6.642 + (0.6 - 6.642)*exp(-0.1323*age))*(1-Hct_M)*60*24,
                             (-0.000895*age^2 + 0.0607*age + 5.54)*(1-Hct_M)*60*24)) %>%
  mutate(CardOut_F = if_else(age < 16.027, (7.734 + (0.6 - 7.734)*exp(-0.09747*age))*(1-Hct_F)*60*24,
                             (0.000473*age^2 - 0.0782*age + 7.37)*(1-Hct_F)*60*24)) %>%
  
  ## Plasma volume; compartment [22] in Ratier 2024 -> USED TO BE BLOOD volume, as it's corrected for hematocrit then it's plasma
  mutate(VplasmaFraction_M = if_else(age < 1, (-0.0273*age + 0.0771),
                                     0.0761 + (0.0289 - 0.0761)*exp(-0.592*age))) %>%
  mutate(Vplasma_M = VplasmaFraction_M*BDW_M) %>%
  mutate(Vart_M = Vplasma_M*Fr_art_plasma*(1-Hct_M)) %>% #arterial blood (plasma)
  mutate(Vven_M = Vplasma_M*(1-Hct_M) - Vart_M) %>% #venous blood (plasma)
  mutate(VplasmaFraction_F = if_else(age < 1, (-0.0273*age + 0.0771),
                                     if_else(age < 14.019723, 3.28E-5*age^3 - 1.21E-3*age^2 + 1.24E-2*age + 3.86E-2,
                                             0.065))) %>%
  mutate(Vplasma_F = VplasmaFraction_F*BDW_F) %>%
  mutate(Vart_F = Vplasma_F*Fr_art_plasma*(1-Hct_F)) %>%
  mutate(Vven_F = Vplasma_F*(1-Hct_F) - Vart_F) %>%
  
  ## Liver; compartments [18] (liver) and [21] (liver artery) in Ratier 2024
  mutate(VliverFraction_M = 0.0247 + (0.0409 - 0.0247)*exp(-0.218*age)) %>%
  mutate(Vliver_M = VliverFraction_M*BDW_M) %>%
  mutate(QliverFraction_M = (VliverFraction_M/0.0247)*0.065) %>% # sc_F[21] = (sc_V[18] / sc_V_adult[18]) * sc_F_adult[21];
  mutate(VliverFraction_F = 0.0233 + (0.038 - 0.0233)*exp(-0.122*age)) %>%
  mutate(Vliver_F = VliverFraction_F*BDW_F) %>%
  mutate(QliverFraction_F = (VliverFraction_F/0.0233)*0.065) %>% # sc_F[21] = (sc_V[18] / sc_V_adult[18]) * sc_F_adult[21];
  
  ## Stomach; compartment [17] in Ratier 2024
  mutate(VstomachFraction_M = 0.0021) %>%
  mutate(Vstomach_M = VstomachFraction_M*BDW_M) %>%
  mutate(QstomachFraction_M = (VstomachFraction_M/0.0021)*0.01) %>% # sc_F[0-17] = (sc_V[i]  / sc_V_adult[i])  * sc_F_adult[i];
  mutate(VstomachFraction_F = 0.0023) %>%
  mutate(Vstomach_F = VstomachFraction_F*BDW_F) %>%
  mutate(QstomachFraction_F = (VstomachFraction_F/0.0023)*0.01) %>% # sc_F[0-17] = (sc_V[i]  / sc_V_adult[i])  * sc_F_adult[i];
  
  ## Gut; compartment [16] in Ratier 2024
  mutate(VgutFraction_M = if_else(age < 16, -0.000082562*age^2 + 0.0013523*age + 0.01293,
                                  0.0140)) %>%
  mutate(Vgut_M = VgutFraction_M*BDW_M) %>%
  mutate(QgutFraction_M = (VgutFraction_M/0.0140)*0.144) %>% # sc_F[0-17] = (sc_V[i]  / sc_V_adult[i])  * sc_F_adult[i];
  mutate(VgutFraction_F = if_else(age < 14.453301, -7.42E-5*age^2 + 1.28E-3*age + 1.30E-2,
                                  0.0160)) %>%
  mutate(Vgut_F = VgutFraction_F*BDW_F) %>%
  mutate(QgutFraction_F = (VgutFraction_F/0.0160)*0.165) %>% # sc_F[0-17] = (sc_V[i]  / sc_V_adult[i])  * sc_F_adult[i];
  
  ## Kidney; compartment [14] in Ratier 2024
  mutate(VkidneyFraction_M = 0.0042 + (0.00767 - 0.0042)*exp(-0.206*age)) %>%
  mutate(Vkidney_M = VkidneyFraction_M*BDW_M) %>%
  mutate(QkidneyFraction_M = (VkidneyFraction_M/0.0042)*0.196) %>% # sc_F[0-17] = (sc_V[i]  / sc_V_adult[i])  * sc_F_adult[i];
  mutate(VkidneyFraction_F = 0.0046 + (0.0071 - 0.0046)*exp(-0.221*age)) %>%
  mutate(Vkidney_F = VkidneyFraction_F*BDW_F) %>%
  mutate(QkidneyFraction_F = (VkidneyFraction_F/0.0046)*0.175) %>% # sc_F[0-17] = (sc_V[i]  / sc_V_adult[i])  * sc_F_adult[i];
  
  ## New Chrysa 05-11-2014
  ## Kidney blood and Kidney tissue volumes; based on Brown1997, Table30
  mutate(VkidneyBlood_M = Vkidney_M*0.36) %>% # 0.36+-0.01 volume fraction of blood in the kidneys
  mutate(VkidneyBlood_F = Vkidney_F*0.36) %>% # 0.36+-0.01 volume fraction of blood in the kidneys
  mutate(VkidneyTissue_M = Vkidney_M*0.64) %>% # 1-0.36
  mutate(VkidneyTissue_F = Vkidney_F*0.64) %>% # 1-0.36
  mutate(VFil_M = Vkidney_M*0.05) %>% #corresponds to the volume of the collecting system in ICRP89
  mutate(VFil_F = Vkidney_F*0.05) %>% #corresponds to the volume of the collecting system in ICRP89
  
  ## Urinary tract (bladder, ureters, urethra); compartment [13] in Ratier 2024
  mutate(VurinarytractFraction_M = 0.00104) %>%
  mutate(Vurinarytract_M = VurinarytractFraction_M*BDW_M) %>%
  mutate(QurinarytractFraction_M = (VurinarytractFraction_M/0.00104)*0.001) %>% # sc_F[0-17] = (sc_V[i]  / sc_V_adult[i])  * sc_F_adult[i];
  mutate(VurinarytractFraction_F = 0.0010) %>%
  mutate(Vurinarytract_F = VurinarytractFraction_F*BDW_F) %>%
  mutate(QurinarytractFraction_F = (VurinarytractFraction_F/0.0010)*0.001) %>% # sc_F[0-17] = (sc_V[i]  / sc_V_adult[i])  * sc_F_adult[i];
  
  ## Skin; compartment [10] in Ratier 2024
  mutate(VskinFraction_M = if_else(age < 20.01, -1.1706E-05*age^3 + 5.4130E-04*age^2 - 6.1966E-03*age + 4.6231E-02,
                                   0.0452)) %>%
  mutate(Vskin_M = VskinFraction_M*BDW_M) %>%
  mutate(QskinFraction_M = (VskinFraction_M/0.0452)*0.052) %>% # sc_F[0-17] = (sc_V[i]  / sc_V_adult[i])  * sc_F_adult[i];
  mutate(VskinFraction_F = if_else(age < 19.45, -7.8882E-06*age^3 + 4.0224E-04*age^2 - 5.2146E-03*age + 4.5605E-02,
                                   0.0383)) %>%
  mutate(Vskin_F = VskinFraction_F*BDW_F) %>%
  mutate(QskinFraction_F = (VskinFraction_F/0.0383)*0.051) %>% # sc_F[0-17] = (sc_V[i]  / sc_V_adult[i])  * sc_F_adult[i];
  
  ## Adrenal; compartment [1] in Ratier 2024 (not used in our model, but needed for calculation of adipose tissue)
  mutate(VadrenalFraction_M = 0.0002 + (0.00171 - 0.0002)*exp(-2.02*age)) %>%
  mutate(Vadrenal_M = VadrenalFraction_M*BDW_M) %>%
  mutate(QadrenalFraction_M = (VadrenalFraction_M/0.0002)*0.003) %>% # sc_F[0-17] = (sc_V[i]  / sc_V_adult[i])  * sc_F_adult[i];
  mutate(VadrenalFraction_F = 0.0002 + (0.00171 - 0.0002)*exp(-2.02*age)) %>%
  mutate(Vadrenal_F = VadrenalFraction_F*BDW_F) %>%
  mutate(QadrenalFraction_F = (VadrenalFraction_F/0.0002)*0.003) %>% # sc_F[0-17] = (sc_V[i]  / sc_V_adult[i])  * sc_F_adult[i];
  
  ## Bone; compartment [2] in Ratier 2024 (not used in our model, but needed for calculation of adipose tissue)
  mutate(VboneFraction_M = (0.313 + (0.506 - 0.313)*exp(-0.0907*age))*0.095) %>%
  mutate(VbonenonperfusedFraction_M = 0.095 - VboneFraction_M) %>%
  mutate(Vbone_M = VboneFraction_M*BDW_M/2) %>% # 2 is bone density
  mutate(Vbonenonperfused_M = VbonenonperfusedFraction_M*BDW_M/2) %>% # 2 is bone density
  mutate(QboneFraction_M = (VboneFraction_M/(0.095*0.32))*0.021) %>% # sc_F[0-17] = (sc_V[i]  / sc_V_adult[i])  * sc_F_adult[i];
  mutate(VboneFraction_F = (0.298 + (0.505 - 0.298)*exp(-0.0792*age))*0.085) %>%
  mutate(VbonenonperfusedFraction_F = 0.085 - VboneFraction_F) %>%
  mutate(Vbone_F = VboneFraction_F*BDW_F/2) %>% # 2 is bone density
  mutate(Vbonenonperfused_F = VbonenonperfusedFraction_F*BDW_F/2) %>% # 2 is bone density
  mutate(QboneFraction_F = (VboneFraction_F/(0.085*0.298))*0.021) %>% # sc_F[0-17] = (sc_V[i]  / sc_V_adult[i])  * sc_F_adult[i];
  
  ## Brain; compartment [3] in Ratier 2024 (not used in our model, but needed for calculation of adipose tissue)
  mutate(VbrainFraction_M = (1.450 + (0.353 - 1.450) * exp (-0.440*age))/BDW_M) %>%
  mutate(Vbrain_M = VbrainFraction_M*BDW_M) %>%
  mutate(QbrainFraction_M = (VbrainFraction_M/0.01986)*0.124) %>% # sc_F[0-17] = (sc_V[i]  / sc_V_adult[i])  * sc_F_adult[i];
  mutate(VbrainFraction_F = (1.300 + (0.347 - 1.300) * exp (-0.573*age))/BDW_F) %>%
  mutate(Vbrain_F = VbrainFraction_F*BDW_F) %>%
  mutate(QbrainFraction_F = (VbrainFraction_F/0.0217)*0.124) %>% # sc_F[0-17] = (sc_V[i]  / sc_V_adult[i])  * sc_F_adult[i];
  
  ## Breast; compartment [4] in Ratier 2024 (not used in our model, but needed for calculation of adipose tissue)
  mutate(VbreastFraction_M = 3.42E-4*1/(1 + exp(-1.42*age + 20.1))) %>%
  mutate(Vbreast_M = VbreastFraction_M*BDW_M) %>%
  mutate(QbreastFraction_M = (VbreastFraction_M/0.00035)*0.0002) %>% # sc_F[0-17] = (sc_V[i]  / sc_V_adult[i])  * sc_F_adult[i];
  mutate(VbreastFraction_F = 0.00833/(1 + exp(-1.92*age+ 28.6))) %>%
  mutate(Vbreast_F = VbreastFraction_F*BDW_F) %>%
  mutate(QbreastFraction_F = (VbreastFraction_F/0.0083)*0.004) %>% # sc_F[0-17] = (sc_V[i]  / sc_V_adult[i])  * sc_F_adult[i];
  
  ## Heart; compartment [5] in Ratier 2024 (not used in our model, but needed for calculation of adipose tissue)
  mutate(VheartFraction_M = 0.0045) %>%
  mutate(Vheart_M = VheartFraction_M*BDW_M) %>%
  mutate(QheartFraction_M = (VheartFraction_M/0.0045)*0.041) %>% # sc_F[0-17] = (sc_V[i]  / sc_V_adult[i])  * sc_F_adult[i];
  mutate(VheartFraction_F = 0.0042) %>%
  mutate(Vheart_F = VheartFraction_F*BDW_F) %>%
  mutate(QheartFraction_F = (VheartFraction_F/0.0042)*0.051) %>% # sc_F[0-17] = (sc_V[i]  / sc_V_adult[i])  * sc_F_adult[i];
  
  ## Marrow; compartment [6] in Ratier 2024 (not used in our model, but needed for calculation of adipose tissue)
  mutate(VmarrowFraction_M = 0.05 + (0.0138 - 0.05)*exp(-0.112*age)) %>%
  mutate(Vmarrow_M = VmarrowFraction_M*BDW_M) %>%
  mutate(QmarrowFraction_M = (VmarrowFraction_M/0.050)*0.031) %>% # sc_F[0-17] = (sc_V[i]  / sc_V_adult[i])  * sc_F_adult[i];
  mutate(VmarrowFraction_F = 0.045 + (0.0138 - 0.045)*exp(-0.136*age)) %>%
  mutate(Vmarrow_F = VmarrowFraction_F*BDW_F) %>%
  mutate(QmarrowFraction_F = (VmarrowFraction_F/0.045)*0.031) %>% # sc_F[0-17] = (sc_V[i]  / sc_V_adult[i])  * sc_F_adult[i];
  
  ## Muscle; compartment [7] in Ratier 2024 (not used in our model, but needed for calculation of adipose tissue)
  mutate(MuscleAtrophy_M = if_else(age < 24.3, 1,
                                   (-0.0001264*age^2 + 0.006131*age + 0.926))) %>%
  mutate(VmuscleFraction_M = (0.3973 + (0.201 - 0.3973)*exp(-0.141*age)) * MuscleAtrophy_M) %>%
  mutate(Vmuscle_M = VmuscleFraction_M*BDW_M) %>%
  mutate(QmuscleFraction_M = (VmuscleFraction_M/0.3973)*0.175) %>% # sc_F[0-17] = (sc_V[i]  / sc_V_adult[i])  * sc_F_adult[i];
  mutate(MuscleAtrophy_F = if_else(age < 25.90709, 1,
                                   (-0.0001264*age^2 + 0.006131*age + 0.926))) %>%
  mutate(VmuscleFraction_F = (0.2917 + (0.207 - 0.2917)*exp(-0.339*age)) * MuscleAtrophy_F) %>%
  mutate(Vmuscle_F = VmuscleFraction_F*BDW_F) %>%
  mutate(QmuscleFraction_F = (VmuscleFraction_F/0.2917)*0.124) %>% # sc_F[0-17] = (sc_V[i]  / sc_V_adult[i])  * sc_F_adult[i];
  
  ## Reproductive organs; compartment [8] in Ratier 2024 (not used in our model, but needed for calculation of adipose tissue)
  mutate(VreproFraction_M = if_else(age < 20.01, -1.5156E-07*age^3 + 9.3351E-06*age^2 - 1.1177E-04*age + 4.7966E-04,
                                    0.0008)) %>%
  mutate(Vrepro_M = VreproFraction_M*BDW_M) %>%
  mutate(QreproFraction_M = (VreproFraction_M/0.0008)*0.001) %>% # sc_F[0-17] = (sc_V[i]  / sc_V_adult[i])  * sc_F_adult[i];
  mutate(VreproFraction_F = if_else(age < 1, -1.064E-3*age + 1.338E-3,
                                    if_else(age < 20, 2.6380E-7*age^3 - 1.7943E-6*age^2 - 5.6465E-6*age + 2.8105E-4,
                                            0.001552))) %>%
  mutate(Vrepro_F = VreproFraction_F*BDW_F) %>%
  mutate(QreproFraction_F = (VreproFraction_F/0.0016)*0.004) %>% # sc_F[0-17] = (sc_V[i]  / sc_V_adult[i])  * sc_F_adult[i];
  
  ## Pancreas; compartment [9] in Ratier 2024 (not used in our model, but needed for calculation of adipose tissue)
  mutate(VpancreasFraction_M = 0.00192) %>%
  mutate(Vpancreas_M = VpancreasFraction_M*BDW_M) %>%
  mutate(QpancreasFraction_M = (VpancreasFraction_M/0.00192)*0.01) %>% # sc_F[0-17] = (sc_V[i]  / sc_V_adult[i])  * sc_F_adult[i];
  mutate(VpancreasFraction_F = 0.002) %>%
  mutate(Vpancreas_F = VpancreasFraction_F*BDW_F) %>%
  mutate(QpancreasFraction_F = (VpancreasFraction_F/0.002)*0.01) %>% # sc_F[0-17] = (sc_V[i]  / sc_V_adult[i])  * sc_F_adult[i];
  
  ## Spleen; compartment [11] in Ratier 2024 (not used in our model, but needed for calculation of adipose tissue)
  mutate(VspleenFraction_M = 0.0021) %>%
  mutate(Vspleen_M = VspleenFraction_M*BDW_M) %>%
  mutate(QspleenFraction_M = (VspleenFraction_M/0.0021)*0.031) %>% # sc_F[0-17] = (sc_V[i]  / sc_V_adult[i])  * sc_F_adult[i];
  mutate(VspleenFraction_F = 0.0022) %>%
  mutate(Vspleen_F = VspleenFraction_F*BDW_F) %>%
  mutate(QspleenFraction_F = (VspleenFraction_F/0.0022)*0.031) %>% # sc_F[0-17] = (sc_V[i]  / sc_V_adult[i])  * sc_F_adult[i];
  
  ## Thyroid; compartment [12] in Ratier 2024 (not used in our model, but needed for calculation of adipose tissue)
  mutate(VthyroidFraction_M = 0.000274) %>%
  mutate(Vthyroid_M = VthyroidFraction_M*BDW_M) %>%
  mutate(QthyroidFraction_M = (VthyroidFraction_M/0.000274)*0.015) %>% # sc_F[0-17] = (sc_V[i]  / sc_V_adult[i])  * sc_F_adult[i];
  mutate(VthyroidFraction_F = 0.0003) %>%
  mutate(Vthyroid_F = VthyroidFraction_F*BDW_F) %>%
  mutate(QthyroidFraction_F = (VthyroidFraction_F/0.0003)*0.015) %>% # sc_F[0-17] = (sc_V[i]  / sc_V_adult[i])  * sc_F_adult[i];
  
  ## Lungs; compartment [15] in Ratier 2024 (not used in our model, but needed for calculation of adipose tissue)
  # Comment Chrysa on 18-10-2024: Shouldn't the lungs take 100% of the blood flow?
  mutate(VlungFraction_M = 0.0068) %>%
  mutate(Vlung_M = VlungFraction_M*BDW_M) %>%
  mutate(QlungFraction_M = (VlungFraction_M/0.0068)*0.026) %>% # sc_F[0-17] = (sc_V[i]  / sc_V_adult[i])  * sc_F_adult[i];
  mutate(VlungFraction_F = 0.0070) %>%
  mutate(Vlung_F = VlungFraction_F*BDW_F) %>%
  mutate(QlungFraction_F = (VlungFraction_F/0.0070)*0.026) %>% # sc_F[0-17] = (sc_V[i]  / sc_V_adult[i])  * sc_F_adult[i];
  
  ## Adipose tissue 
  mutate(VadiposeFraction_M = 0.96 - VadrenalFraction_M - VboneFraction_M - VbonenonperfusedFraction_M - VbrainFraction_M - VbreastFraction_M - 
           VheartFraction_M - VmarrowFraction_M - VmuscleFraction_M - VreproFraction_M - VpancreasFraction_M -
           VskinFraction_M - VspleenFraction_M - VthyroidFraction_M - VurinarytractFraction_M - VkidneyFraction_M -
           VlungFraction_M - VgutFraction_M - VstomachFraction_M - VliverFraction_M - VplasmaFraction_M) %>%
  mutate(AdiposeMass_M = if_else(age < 19.00093277, 0,
                                 (-0.01129273*age^2 + 1.11817056*age + 56.74397436)-(74.16235828-(2*(74.16235828-57.19957758)/(exp(0.63466182*(age-13.31018000))+exp(0.05457656*(age-13.31018000))))))) %>%
  mutate(Vadipose_M = (AdiposeMass_M/0.9) + VadiposeFraction_M*BDW_M/0.9) %>% # 0.9 is adipose tissue density
  mutate(QadiposeFraction_M = (VadiposeFraction_M/0.20)*0.052) %>% # sc_F[0-17] = (sc_V[i]  / sc_V_adult[i])  * sc_F_adult[i];
  mutate(VadiposeFraction_F = 0.96 - VadrenalFraction_F - VboneFraction_F - VbonenonperfusedFraction_F - VbrainFraction_F - VbreastFraction_F - 
           VheartFraction_F - VmarrowFraction_F - VmuscleFraction_F - VreproFraction_F - VpancreasFraction_F -
           VskinFraction_F - VspleenFraction_F - VthyroidFraction_F - VurinarytractFraction_F - VkidneyFraction_F -
           VlungFraction_F - VgutFraction_F - VstomachFraction_F - VliverFraction_F - VplasmaFraction_F) %>%
  mutate(AdiposeMass_F = if_else(age < 17.9374115, 0,
                                 if_else(((-0.01258006*age^2 + 1.25029379*age + 44.4459234) - (62.95490567-(2*(62.95490567-49.36574299)/(exp(0.84039606*(age-11.56691488))+exp(0.06710088*(age-11.56691488)))))) < 0, 0,
                                         (-0.01258006*age^2 + 1.25029379*age + 44.4459234) - (62.95490567-(2*(62.95490567-49.36574299)/(exp(0.84039606*(age-11.56691488))+exp(0.06710088*(age-11.56691488)))))))) %>%
  mutate(Vadipose_F = (AdiposeMass_F/0.9) + VadiposeFraction_F*BDW_F/0.9) %>% # 0.9 is adipose tissue density
  mutate(QadiposeFraction_F = (VadiposeFraction_F/0.3167)*0.087) %>% # sc_F[0-17] = (sc_V[i]  / sc_V_adult[i])  * sc_F_adult[i];
  
  ## Total volume (sum of all organs) 
  # Comment Chrysa: Why is this changing over time and not staying constant?
  mutate(VtotalFraction_M = VplasmaFraction_M + VadrenalFraction_M + VboneFraction_M + VbrainFraction_M + VbreastFraction_M + 
           VheartFraction_M + VmarrowFraction_M + VmuscleFraction_M + VreproFraction_M + VpancreasFraction_M +
           VskinFraction_M + VspleenFraction_M + VthyroidFraction_M + VurinarytractFraction_M + VkidneyFraction_M +
           VlungFraction_M + VgutFraction_M + VstomachFraction_M + VliverFraction_M + VadiposeFraction_M) %>%
  mutate(VtotalFraction_F = VplasmaFraction_F + VadrenalFraction_F + VboneFraction_F + VbrainFraction_F + VbreastFraction_F + 
           VheartFraction_F + VmarrowFraction_F + VmuscleFraction_F + VreproFraction_F + VpancreasFraction_F +
           VskinFraction_F + VspleenFraction_F + VthyroidFraction_F + VurinarytractFraction_F + VkidneyFraction_F +
           VlungFraction_F + VgutFraction_F + VstomachFraction_F + VliverFraction_F + VadiposeFraction_F) %>%
  
  mutate(TestV_M = BDW_M - VtotalFraction_M) %>% 
  mutate(TestV_F = BDW_F - VtotalFraction_F) %>% 
  
  mutate(Vtotal_M = Vplasma_M + Vadrenal_M + Vbone_M + Vbrain_M + Vbreast_M + 
           Vheart_M + Vmarrow_M + Vmuscle_M + Vrepro_M + Vpancreas_M +
           Vskin_M + Vspleen_M + Vthyroid_M + Vurinarytract_M + Vkidney_M +
           Vlung_M + Vgut_M + Vstomach_M + Vliver_M + Vadipose_M) %>%
  mutate(Vtotal_F = Vplasma_F + Vadrenal_F + Vbone_F + Vbrain_F + Vbreast_F + 
           Vheart_F + Vmarrow_F + Vmuscle_F + Vrepro_F + Vpancreas_F +
           Vskin_F + Vspleen_F + Vthyroid_F + Vurinarytract_F + Vkidney_F +
           Vlung_F + Vgut_F + Vstomach_F + Vliver_F + Vadipose_F) %>%
  
  ## Flow rates
  mutate(QtotalFraction_M = QadrenalFraction_M + QboneFraction_M + QbrainFraction_M + QbreastFraction_M +
           QheartFraction_M + QmarrowFraction_M + QmuscleFraction_M + QreproFraction_M + QpancreasFraction_M +
           QskinFraction_M + QspleenFraction_M + QthyroidFraction_M + QurinarytractFraction_M + QkidneyFraction_M +
           QlungFraction_M + QgutFraction_M + QstomachFraction_M + QliverFraction_M + QadiposeFraction_M) %>%
  mutate(QtotalFraction_F = QadrenalFraction_F + QboneFraction_F + QbrainFraction_F + QbreastFraction_F +
           QheartFraction_F + QmarrowFraction_F + QmuscleFraction_F + QreproFraction_F + QpancreasFraction_F +
           QskinFraction_F + QspleenFraction_F + QthyroidFraction_F + QurinarytractFraction_F + QkidneyFraction_F +
           QlungFraction_F + QgutFraction_F + QstomachFraction_F + QliverFraction_F + QadiposeFraction_F) %>%
  
  mutate(TestQF_M = CardOut_M - QtotalFraction_M) %>% 
  mutate(TestQF_F = CardOut_F - QtotalFraction_F) %>% 
  
  mutate(Qliver_M = (QliverFraction_M/QtotalFraction_M)*CardOut_M) %>% # Note the inclusion of the free fraction: F[21] = (sc_F[21] / (sum_sc_F / 1)) * CardOut * Free;
  # Qliver is what was used to be Qhepatic
  mutate(Qstomach_M = (QstomachFraction_M/QtotalFraction_M)*CardOut_M) %>% # Note the inclusion of the free fraction: F[0-17] = (sc_F[i]  / (sum_sc_F / 1)) * CardOut * Free;
  mutate(Qgut_M = (QgutFraction_M/QtotalFraction_M)*CardOut_M) %>% # Note the inclusion of the free fraction: F[0-17] = (sc_F[i]  / (sum_sc_F / 1)) * CardOut * Free;
  mutate(Qkidney_M = (QkidneyFraction_M/QtotalFraction_M)*CardOut_M) %>% # Note the inclusion of the free fraction: F[0-17] = (sc_F[i]  / (sum_sc_F / 1)) * CardOut * Free;
  mutate(Qurinarytract_M = (QurinarytractFraction_M/QtotalFraction_M)*CardOut_M) %>% # Note the inclusion of the free fraction: F[0-17] = (sc_F[i]  / (sum_sc_F / 1)) * CardOut * Free;
  mutate(Qskin_M = (QskinFraction_M/QtotalFraction_M)*CardOut_M) %>% # Note the inclusion of the free fraction: F[0-17] = (sc_F[i]  / (sum_sc_F / 1)) * CardOut * Free;
  mutate(Qadrenal_M = (QadrenalFraction_M/QtotalFraction_M)*CardOut_M) %>% # Note the inclusion of the free fraction: F[0-17] = (sc_F[i]  / (sum_sc_F / 1)) * CardOut * Free;
  mutate(Qbone_M = (QboneFraction_M/QtotalFraction_M)*CardOut_M) %>% # Note the inclusion of the free fraction: F[0-17] = (sc_F[i]  / (sum_sc_F / 1)) * CardOut * Free;
  mutate(Qbrain_M = (QbrainFraction_M/QtotalFraction_M)*CardOut_M) %>% # Note the inclusion of the free fraction: F[0-17] = (sc_F[i]  / (sum_sc_F / 1)) * CardOut * Free;
  mutate(Qbreast_M = (QbreastFraction_M/QtotalFraction_M)*CardOut_M) %>% # Note the inclusion of the free fraction: F[0-17] = (sc_F[i]  / (sum_sc_F / 1)) * CardOut * Free;
  mutate(Qheart_M = (QheartFraction_M/QtotalFraction_M)*CardOut_M) %>% # Note the inclusion of the free fraction: F[0-17] = (sc_F[i]  / (sum_sc_F / 1)) * CardOut * Free;
  mutate(Qmarrow_M = (QmarrowFraction_M/QtotalFraction_M)*CardOut_M) %>% # Note the inclusion of the free fraction: F[0-17] = (sc_F[i]  / (sum_sc_F / 1)) * CardOut * Free;
  mutate(Qmuscle_M = (QmuscleFraction_M/QtotalFraction_M)*CardOut_M) %>% # Note the inclusion of the free fraction: F[0-17] = (sc_F[i]  / (sum_sc_F / 1)) * CardOut * Free;
  mutate(Qrepro_M = (QreproFraction_M/QtotalFraction_M)*CardOut_M) %>% # Note the inclusion of the free fraction: F[0-17] = (sc_F[i]  / (sum_sc_F / 1)) * CardOut * Free;
  mutate(Qpancreas_M = (QpancreasFraction_M/QtotalFraction_M)*CardOut_M) %>% # Note the inclusion of the free fraction: F[0-17] = (sc_F[i]  / (sum_sc_F / 1)) * CardOut * Free;
  mutate(Qspleen_M = (QspleenFraction_M/QtotalFraction_M)*CardOut_M) %>% # Note the inclusion of the free fraction: F[0-17] = (sc_F[i]  / (sum_sc_F / 1)) * CardOut * Free;
  mutate(Qthyroid_M = (QthyroidFraction_M/QtotalFraction_M)*CardOut_M) %>% # Note the inclusion of the free fraction: F[0-17] = (sc_F[i]  / (sum_sc_F / 1)) * CardOut * Free;
  mutate(Qlung_M = (QlungFraction_M/QtotalFraction_M)*CardOut_M) %>% # Note the inclusion of the free fraction: F[0-17] = (sc_F[i]  / (sum_sc_F / 1)) * CardOut * Free;
  mutate(Qadipose_M = (QadiposeFraction_M/QtotalFraction_M)*CardOut_M) %>% # Note the inclusion of the free fraction: F[0-17] = (sc_F[i]  / (sum_sc_F / 1)) * CardOut * Free;
  
  mutate(Qliver_F = (QliverFraction_F/QtotalFraction_F)*CardOut_F) %>% # Note the inclusion of the free fraction: F[21] = (sc_F[21] / (sum_sc_F / 1)) * CardOut * Free;
  mutate(Qstomach_F = (QstomachFraction_F/QtotalFraction_F)*CardOut_F) %>% # Note the inclusion of the free fraction: F[0-17] = (sc_F[i]  / (sum_sc_F / 1)) * CardOut * Free;
  mutate(Qgut_F = (QgutFraction_F/QtotalFraction_F)*CardOut_F) %>% # Note the inclusion of the free fraction: F[0-17] = (sc_F[i]  / (sum_sc_F / 1)) * CardOut * Free;
  mutate(Qkidney_F = (QkidneyFraction_F/QtotalFraction_F)*CardOut_F) %>% # Note the inclusion of the free fraction: F[0-17] = (sc_F[i]  / (sum_sc_F / 1)) * CardOut * Free;
  mutate(Qurinarytract_F = (QurinarytractFraction_F/QtotalFraction_F)*CardOut_F) %>% # Note the inclusion of the free fraction: F[0-17] = (sc_F[i]  / (sum_sc_F / 1)) * CardOut * Free;
  mutate(Qskin_F = (QskinFraction_F/QtotalFraction_F)*CardOut_F) %>% # Note the inclusion of the free fraction: F[0-17] = (sc_F[i]  / (sum_sc_F / 1)) * CardOut * Free;
  mutate(Qadrenal_F = (QadrenalFraction_F/QtotalFraction_F)*CardOut_F) %>% # Note the inclusion of the free fraction: F[0-17] = (sc_F[i]  / (sum_sc_F / 1)) * CardOut * Free;
  mutate(Qbone_F = (QboneFraction_F/QtotalFraction_F)*CardOut_F) %>% # Note the inclusion of the free fraction: F[0-17] = (sc_F[i]  / (sum_sc_F / 1)) * CardOut * Free;
  mutate(Qbrain_F = (QbrainFraction_F/QtotalFraction_F)*CardOut_F) %>% # Note the inclusion of the free fraction: F[0-17] = (sc_F[i]  / (sum_sc_F / 1)) * CardOut * Free;
  mutate(Qbreast_F = (QbreastFraction_F/QtotalFraction_F)*CardOut_F) %>% # Note the inclusion of the free fraction: F[0-17] = (sc_F[i]  / (sum_sc_F / 1)) * CardOut * Free;
  mutate(Qheart_F = (QheartFraction_F/QtotalFraction_F)*CardOut_F) %>% # Note the inclusion of the free fraction: F[0-17] = (sc_F[i]  / (sum_sc_F / 1)) * CardOut * Free;
  mutate(Qmarrow_F = (QmarrowFraction_F/QtotalFraction_F)*CardOut_F) %>% # Note the inclusion of the free fraction: F[0-17] = (sc_F[i]  / (sum_sc_F / 1)) * CardOut * Free;
  mutate(Qmuscle_F = (QmuscleFraction_F/QtotalFraction_F)*CardOut_F) %>% # Note the inclusion of the free fraction: F[0-17] = (sc_F[i]  / (sum_sc_F / 1)) * CardOut * Free;
  mutate(Qrepro_F = (QreproFraction_F/QtotalFraction_F)*CardOut_F) %>% # Note the inclusion of the free fraction: F[0-17] = (sc_F[i]  / (sum_sc_F / 1)) * CardOut * Free;
  mutate(Qpancreas_F = (QpancreasFraction_F/QtotalFraction_F)*CardOut_F) %>% # Note the inclusion of the free fraction: F[0-17] = (sc_F[i]  / (sum_sc_F / 1)) * CardOut * Free;
  mutate(Qspleen_F = (QspleenFraction_F/QtotalFraction_F)*CardOut_F) %>% # Note the inclusion of the free fraction: F[0-17] = (sc_F[i]  / (sum_sc_F / 1)) * CardOut * Free;
  mutate(Qthyroid_F = (QthyroidFraction_F/QtotalFraction_F)*CardOut_F) %>% # Note the inclusion of the free fraction: F[0-17] = (sc_F[i]  / (sum_sc_F / 1)) * CardOut * Free;
  mutate(Qlung_F = (QlungFraction_F/QtotalFraction_F)*CardOut_F) %>% # Note the inclusion of the free fraction: F[0-17] = (sc_F[i]  / (sum_sc_F / 1)) * CardOut * Free;
  mutate(Qadipose_F = (QadiposeFraction_F/QtotalFraction_F)*CardOut_F) %>% # Note the inclusion of the free fraction: F[0-17] = (sc_F[i]  / (sum_sc_F / 1)) * CardOut * Free;
  
  ## Total flow
  mutate(Qtotal_M = Qadrenal_M + Qbone_M + Qbrain_M + Qbreast_M + 
           Qheart_M + Qmarrow_M + Qmuscle_M + Qrepro_M + Qpancreas_M +
           Qskin_M + Qspleen_M + Qthyroid_M + Qurinarytract_M + Qkidney_M +
           Qlung_M + Qgut_M + Qstomach_M + Qliver_M + Qadipose_M) %>%
  mutate(Qtotal_F = Qadrenal_F + Qbone_F + Qbrain_F + Qbreast_F + 
           Qheart_F + Qmarrow_F + Qmuscle_F + Qrepro_F + Qpancreas_F +
           Qskin_F + Qspleen_F + Qthyroid_F + Qurinarytract_F + Qkidney_F +
           Qlung_F + Qgut_F + Qstomach_F + Qliver_F + Qadipose_F) %>% 
  
  
##### CHANGED ----
  
  ## Skin barrier or epidermis (stratum corneum & viable epidermis)
  mutate(SkbTarea_M = 9.1*((BW_M*1000)^0.666)) %>%  #Skin barrier (epidermis = stratum corneum and viable epidermis) area; assumed to be the same as the total area of the skin (cm^2) from Husoy
  mutate(SkbTarea_F = 9.1*((BW_F*1000)^0.666)) %>%  #Skin barrier (epidermis = stratum corneum and viable epidermis) area; assumed to be the same as the total area of the skin (cm^2) from Husoy
  mutate(fSkbarea_M = skin_fraction*SkbTarea_M) %>% # 1070 (cm^2); exposed skin area = surface area of the hands [https://www.epa.gov/sites/default/files/2015-09/documents/efh-chapter07.pdf]
  mutate(fSkbarea_F = skin_fraction*SkbTarea_F) %>% # 1070 (cm^2); exposed skin area = surface area of the hands [https://www.epa.gov/sites/default/files/2015-09/documents/efh-chapter07.pdf]
  # checked the reference, 0.972m^2 is the 95th percentile of the lower extremities; I think that was a mistake when copying the information as the hand data is just above. For hands the 95th percentile is 0.131 (mean value is 0.107m^2)
  mutate(Skbthickness_M = 83.1/1000) %>% # (cm) ref: DOI: 10.1080/00015550310015419; in Husoy this was 0.1
  mutate(Skbthickness_F = 83.1/1000) %>% # (cm) ref: DOI: 10.1080/00015550310015419; in Husoy this was 0.1
  mutate(VSkb_M = (fSkbarea_M*Skbthickness_M)/1000) %>%  #(L); Skin barrier volume; as previously coded by Trine
  mutate(VSkb_F = (fSkbarea_F*Skbthickness_F)/1000) %>%  #(L); Skin barrier volume; as previously coded by Trine 
  #mutate(QSkb_M = Qskin_M * (fSkbarea_M/SkbTarea_M)) %>% #(L/h); NOT USED IN THE CODE Plasma flow to the skin barrier volume; as previously coded by Trine
  #mutate(QSkb_F = Qskin_F * (fSkbarea_F/SkbTarea_F)) %>% #(L/h); NOT USED IN THE CODE Plasma flow to the skin barrier volume; as previously coded by Trine
  mutate(CLdermalabs_M = ((Papp*fSkbarea_M)/1000)*24) %>% # (L/d) ; cm/h*cm^2 = mL/h /1000 = L/h * 24 = L/d
  mutate(CLdermalabs_F = ((Papp*fSkbarea_F)/1000)*24) %>% # (L/d) ; cm/h*cm^2 = mL/h /1000 = L/h * 24 = L/d
  
  ## This is the portal vein input to the liver, all excluding Qgut which is added individually as we have a gut compartment
  mutate(Qhepatic_M = Qliver_M + Qspleen_M + Qstomach_M + Qpancreas_M) %>% 
  mutate(Qhepatic_F = Qliver_F + Qspleen_F + Qstomach_F + Qpancreas_F) %>% 
  
  ## Rest (Rest = Total - Organs included in the model)
  # mutate(VrestFraction_M = VtotalFraction_M - VskinFraction_M - VurinarytractFraction_M - VkidneyFraction_M - VgutFraction_M - VliverFraction_M - VplasmaFraction_M - VadiposeFraction_M) %>%
  # mutate(Vrest_M = VrestFraction_M*BDW_M) %>%
  # mutate(QrestFraction_M = Qtotal_M - Qskin_M - Qurinarytract_M - Qkidney_M - Qhepatic_M - Qgut_M - Qadipose_M) %>% #Qhepatic includes Qliver, Qstomach and Qpancreas
  # mutate(Qrest_M = QrestFraction_M/QtotalFraction_M*CardOut_M) %>% 
  # mutate(VrestFraction_F = VtotalFraction_F - VskinFraction_F - VurinarytractFraction_F - VkidneyFraction_F - VgutFraction_F - VliverFraction_F - VplasmaFraction_F - VadiposeFraction_F) %>%
  # mutate(Vrest_F = VrestFraction_F*BDW_F) %>%
  # mutate(QrestFraction_F = Qtotal_F - Qskin_F - Qurinarytract_F - Qkidney_F - Qhepatic_F - Qgut_F - Qadipose_F) %>% #Qhepatic includes Qliver, Qstomach and Qpancreas
  # mutate(Qrest_F = QrestFraction_F/QtotalFraction_F*CardOut_F) %>% 
  mutate(Vrest_M = Vtotal_M - Vskin_M - Vkidney_M - Vgut_M - Vliver_M - Vplasma_M - Vadipose_M) %>%
  mutate(Qrest_M = Qtotal_M - Qskin_M - Qkidney_M - Qhepatic_M - Qgut_M - Qadipose_M) %>%
  mutate(Vrest_F = Vtotal_F - Vskin_F - Vkidney_F - Vgut_F - Vliver_F - Vplasma_F - Vadipose_F) %>%
  mutate(Qrest_F = Qtotal_F - Qskin_F - Qkidney_F - Qhepatic_F - Qgut_F - Qadipose_F) %>%
  
  ## MassBalance Flow #better be 0
  # !!!!! Issue with the mass balance of the flows
  mutate(Qmass_balance_M = CardOut_M - (Qskin_M + Qkidney_M + Qgut_M + Qhepatic_M + Qadipose_M + Qrest_M)) %>% #Qhepatic includes Qliver, Qstomach and Qpancreas
  mutate(Qmass_balance_F = CardOut_F - (Qskin_F + Qkidney_F + Qgut_F + Qhepatic_F + Qadipose_F + Qrest_F)) %>% #Qhepatic includes Qliver, Qstomach and Qpancreas
  
  ## MassBalance Volumes #better be 0
  # !!!!! Issue with the mass balance of the flows
  mutate(Vmass_balance_M = BW_M - (Vplasma_M + Vskin_M + Vkidney_M + Vgut_M + Vliver_M + Vadipose_M + Vrest_M)) %>% 
  mutate(Vmass_balance_F = BW_F - (Vplasma_F + Vskin_F + Vkidney_F + Vgut_F + Vliver_F + Vadipose_F + Vrest_F)) %>% 
  
  ## Biliary clearance
  # mutate(CLbiliaryFraction_M = 0.000109) %>% #from Husoy; 2.62 biliary clearance L/h/kg calculated from the biliary clearance of 2.62 ml/day taken from Fujii et al 2015 
  # mutate(CLbiliary_M = CLbiliaryFraction_M*BDW_M^VliverFraction_M) %>% #from Husoy; L/h biliary clearance rate, corrected for the fraction of liver to the total BDW (as done by Husoy)
  # mutate(CLbiliaryFraction_F = 0.000109) %>% #from Husoy; 2.62/24/1000 biliary clearance L/h/kg calculated from the biliary clearance of 2.62 ml/day taken from Fujii et al 2015 
  # mutate(CLbiliary_F = CLbiliaryFraction_M*BDW_F^VliverFraction_F) %>% #from Husoy; L/h biliary clearance rate, corrected for the fraction of liver to the total BDW (as done by Husoy)
  # Comment Chrysa 18-10-2024: Shouldn't we be correcting for the fraction of liver changing during life?
  mutate(CLbiliary_M = CLbiliaryc*(BW_M^0.1)) %>% #from Husoy; L/d biliary clearance rate, corrected for the fraction of liver to the total BDW (as done by Husoy)
  mutate(CLbiliary_F = CLbiliaryc*(BW_F^0.1)) %>% #from Husoy; L/d biliary clearance rate, corrected for the fraction of liver to the total BDW (as done by Husoy)
  
  ## Fecal clearance
  # mutate(CLfecalFraction_M = 0.00262) %>% #from Husoy; 0.052/24/1000 feces clearance L/h/kg clearance in feces taken from Fujii et al 2015, calculated from 0.052 ml/day/kg
  # mutate(CLfecal_M = CLfecalFraction_M*BDW_M^VgutFraction_M) %>% # L/h faeces clearance, BDW adjusted to the volume of GI tract as done by Husoy
  # mutate(CLfecalFraction_F = 0.000052) %>% #from Husoy; 0.052/24/1000 feces clearance L/h/kg clearance in feces taken from Fujii et al 2015, calculated from 0.052 ml/day/kg
  # mutate(CLfecal_F = CLfecalFraction_F*BDW_F^VgutFraction_F) %>% # L/h faeces clearance, BDW adjusted to the volume of GI tract as done by Husoy
  mutate(CLfecal_M = CLfaecesc*(BW_M^0.001)) %>% #from Husoy; L/d faeces clearance, BDW adjusted to the volume of GI tract as done by Husoy
  mutate(CLfecal_F = CLfaecesc*(BW_F^0.001)) %>% #from Husoy; L/d faeces clearance, BDW adjusted to the volume of GI tract as done by Husoy
  
  ## Urinary clearance 
  mutate(QUr_M = 0.022*BW_M) %>% # 22 mL/kg BW/d -> L/d urine flow rate to the bladder [ICRP 89 page 161] http://www.icrp.org/publication.asp?id=ICRP%20Publication%2089
  mutate(QUr_F = 0.022*BW_F) %>% # 22 mL/kg BW/d -> L/d urine flow rate to the bladder [ICRP 89 page 161] http://www.icrp.org/publication.asp?id=ICRP%20Publication%2089 
  mutate(GFR_M = 0.18*Qkidney_M) %>% #L/d ; as it's 18% of total renal plasma flow [ICRP 89] http://www.icrp.org/publication.asp?id=ICRP%20Publication%2089
  mutate(GFR_F = 0.18*Qkidney_F) %>% # L/d; as it's 18% of total renal plasma flow [ICRP 89]
  # mutate(CLurine_M = CLurinec*BW_M^(-0.25)) %>% # clearance urine (L/d) # NOT USED
  # mutate(CLurine_F = CLurinec*BW_F^(-0.25)) %>% # clearance urine (L/d) # NOT USED
  # mutate(MPT_M = Vkidney_M*1000*Kcells*Kprotein) %>%	# mass proximal tubule cells in gram based on BW
  # mutate(MPT_F = Vkidney_F*1000*Kcells*Kprotein) %>%	# mass proximal tubule cells in gram based on BW
  # mutate(Vmax_M = Vmaxc * MPT_M * SFOAT4) %>% #ug/d
  # mutate(Vmax_F = Vmaxc * MPT_F * SFOAT4) %>% #ug/d
  mutate(KW_cortex = 0.7*Vkidney_M) %>% # Kg kidney cortexes, only scaling to kidney cortex volume as proximal tubule cells are in the cortex; 70% of the total kidney volume according to ICRP89; PT are in the cortex https://doi.org/10.1021/acs.molpharmaceut.4c00504; alternatively we could have 68% of kidney weight https://doi.org/10.1124/dmd.117.075242 
  mutate(CL_PltPT = ((CL_OAT1*REF_OAT1) + (CL_OAT3*REF_OAT3)) * PTCPGK * KW_cortex) %>% #L/d plasma to proximal tubule clearance
  mutate(CL_FiltPT = (CL_OAT4*REF_OAT4) * PTCPGK * KW_cortex) %>%  #L/d filtrate to proximal tubule clearance 
  mutate(CL_passive <- ((Papp_kidney*PT_Sarea)/1000*60*24)) # (L/D) ; 


# view(Variables_df)

##### Gender-specific data frames ----  
# Variables_M <- Variables_df %>%
#   select(TIME, age, BW_M, Hct_M, CardOut_M, Vplasma_M, Vart_M, Vven_M, Vliver_M, Vstomach_M, Vgut_M, Vkidney_M, Vurinarytract_M,
#          Vskin_M, Vadipose_M, Vadrenal_M, Vbone_M, Vbonenonperfused_M, Vbrain_M, Vbreast_M, Vheart_M, Vmarrow_M, Vmuscle_M,
#          Vrepro_M, Vpancreas_M, Vspleen_M, Vthyroid_M, Vlung_M,
#          Qliver_M, Qstomach_M, Qgut_M, Qkidney_M, Qurinarytract_M,
#          Qskin_M, Qadipose_M, Qadrenal_M, Qbone_M, Qbrain_M, Qbreast_M, Qheart_M, Qmarrow_M, Qmuscle_M,
#          Qrepro_M, Qpancreas_M, Qspleen_M, Qthyroid_M, Qlung_M) %>%
#   mutate(Vtotal_M = Vplasma_M + Vliver_M + Vstomach_M + Vgut_M + Vkidney_M + Vurinarytract_M +
#            Vskin_M + Vadipose_M + Vadrenal_M + Vbone_M + Vbonenonperfused_M + Vbrain_M + Vbreast_M + Vheart_M + Vmarrow_M + Vmuscle_M +
#            Vrepro_M + Vpancreas_M + Vspleen_M + Vthyroid_M + Vlung_M) %>%
#   mutate(Fraction_M = Vtotal_M/BW_M) %>%
#   mutate(Qtotal_M = Qliver_M + Qstomach_M + Qgut_M + Qkidney_M + Qurinarytract_M +
#            Qskin_M + Qadipose_M + Qadrenal_M + Qbone_M + Qbrain_M + Qbreast_M + Qheart_M + Qmarrow_M + Qmuscle_M +
#            Qrepro_M + Qpancreas_M + Qspleen_M + Qthyroid_M + Qlung_M)
# 
# Variables_fractions_M <- Variables_df %>%
#   select(TIME, age, VadrenalFraction_M, VboneFraction_M, VbonenonperfusedFraction_M, VbrainFraction_M, VbreastFraction_M, 
#          VheartFraction_M, VmarrowFraction_M, VmuscleFraction_M, VreproFraction_M, VpancreasFraction_M,
#          VskinFraction_M, VspleenFraction_M, VthyroidFraction_M, VurinarytractFraction_M, VkidneyFraction_M,
#          VlungFraction_M, VgutFraction_M, VstomachFraction_M, VliverFraction_M, VplasmaFraction_M, VadiposeFraction_M,
#          QadrenalFraction_M, QboneFraction_M, QbrainFraction_M, QbreastFraction_M, 
#          QheartFraction_M, QmarrowFraction_M, QmuscleFraction_M, QreproFraction_M, QpancreasFraction_M,
#          QskinFraction_M, QspleenFraction_M, QthyroidFraction_M, QurinarytractFraction_M, QkidneyFraction_M,
#          QlungFraction_M, QgutFraction_M, QstomachFraction_M, QliverFraction_M, QadiposeFraction_M) %>%
#   mutate(VtotalFraction_M = VadrenalFraction_M + VboneFraction_M + VbonenonperfusedFraction_M + VbrainFraction_M + VbreastFraction_M + 
#            VheartFraction_M + VmarrowFraction_M + VmuscleFraction_M + VreproFraction_M + VpancreasFraction_M +
#            VskinFraction_M + VspleenFraction_M + VthyroidFraction_M + VurinarytractFraction_M + VkidneyFraction_M +
#            VlungFraction_M + VgutFraction_M + VstomachFraction_M + VliverFraction_M + VplasmaFraction_M + VadiposeFraction_M) %>%
#   mutate(QtotalFraction_M = QadrenalFraction_M + QboneFraction_M + QbrainFraction_M + QbreastFraction_M + 
#            QheartFraction_M + QmarrowFraction_M + QmuscleFraction_M + QreproFraction_M + QpancreasFraction_M +
#            QskinFraction_M + QspleenFraction_M + QthyroidFraction_M + QurinarytractFraction_M + QkidneyFraction_M +
#            QlungFraction_M + QgutFraction_M + QstomachFraction_M + QliverFraction_M + QadiposeFraction_M)
# 
# Variables_F <- Variables_df %>%
#   select(TIME, age, BW_F, Hct_F, CardOut_F, Vplasma_F, Vart_F, Vven_F, Vliver_F, Vstomach_F, Vgut_F, Vkidney_F, Vurinarytract_F,
#          Vskin_F, Vadipose_F, Vadrenal_F, Vbone_F, Vbonenonperfused_F, Vbrain_F, Vbreast_F, Vheart_F, Vmarrow_F, Vmuscle_F,
#          Vrepro_F, Vpancreas_F, Vspleen_F, Vthyroid_F, Vlung_F,
#          Qliver_F, Qstomach_F, Qgut_F, Qkidney_F, Qurinarytract_F,
#          Qskin_F, Qadipose_F, Qadrenal_F, Qbone_F, Qbrain_F, Qbreast_F, Qheart_F, Qmarrow_F, Qmuscle_F,
#          Qrepro_F, Qpancreas_F, Qspleen_F, Qthyroid_F, Qlung_F) %>%
#   mutate(Vtotal_F = Vplasma_F + Vliver_F + Vstomach_F + Vgut_F + Vkidney_F + Vurinarytract_F +
#            Vskin_F + Vadipose_F + Vadrenal_F + Vbone_F + Vbonenonperfused_F + Vbrain_F + Vbreast_F + Vheart_F + Vmarrow_F + Vmuscle_F +
#            Vrepro_F + Vpancreas_F + Vspleen_F + Vthyroid_F + Vlung_F) %>%
#   mutate(Fraction_F = Vtotal_F/BW_F) %>%
#   mutate(Qtotal_F = Qliver_F + Qstomach_F + Qgut_F + Qkidney_F + Qurinarytract_F +
#            Qskin_F + Qadipose_F + Qadrenal_F + Qbone_F + Qbrain_F + Qbreast_F + Qheart_F + Qmarrow_F + Qmuscle_F +
#            Qrepro_F + Qpancreas_F + Qspleen_F + Qthyroid_F + Qlung_F)
# 
# Variables_fractions_F <- Variables_df %>%
#   select(TIME, age, VadrenalFraction_F, VboneFraction_F, VbonenonperfusedFraction_F, VbrainFraction_F, VbreastFraction_F, 
#          VheartFraction_F, VmarrowFraction_F, VmuscleFraction_F, VreproFraction_F, VpancreasFraction_F,
#          VskinFraction_F, VspleenFraction_F, VthyroidFraction_F, VurinarytractFraction_F, VkidneyFraction_F,
#          VlungFraction_F, VgutFraction_F, VstomachFraction_F, VliverFraction_F, VplasmaFraction_F, VadiposeFraction_F,
#          QadrenalFraction_F, QboneFraction_F, QbrainFraction_F, QbreastFraction_F, 
#          QheartFraction_F, QmarrowFraction_F, QmuscleFraction_F, QreproFraction_F, QpancreasFraction_F,
#          QskinFraction_F, QspleenFraction_F, QthyroidFraction_F, QurinarytractFraction_F, QkidneyFraction_F,
#          QlungFraction_F, QgutFraction_F, QstomachFraction_F, QliverFraction_F, QadiposeFraction_F) %>%
#   mutate(VtotalFraction_F = VadrenalFraction_F + VboneFraction_F + VbonenonperfusedFraction_F + VbrainFraction_F + VbreastFraction_F + 
#            VheartFraction_F + VmarrowFraction_F + VmuscleFraction_F + VreproFraction_F + VpancreasFraction_F +
#            VskinFraction_F + VspleenFraction_F + VthyroidFraction_F + VurinarytractFraction_F + VkidneyFraction_F +
#            VlungFraction_F + VgutFraction_F + VstomachFraction_F + VliverFraction_F + VplasmaFraction_F + VadiposeFraction_F) %>%
#   mutate(QtotalFraction_F = QadrenalFraction_F + QboneFraction_F + QbrainFraction_F + QbreastFraction_F + 
#            QheartFraction_F + QmarrowFraction_F + QmuscleFraction_F + QreproFraction_F + QpancreasFraction_F +
#            QskinFraction_F + QspleenFraction_F + QthyroidFraction_F + QurinarytractFraction_F + QkidneyFraction_F +
#            QlungFraction_F + QgutFraction_F + QstomachFraction_F + QliverFraction_F + QadiposeFraction_F)
# 
# Figure male organ volumes
# Figure_M <- ggplot() +
#   geom_line(data = Variables_M, aes(x=age, y=BW_M, color="01_BW")) +
#   geom_line(data = Variables_M, aes(x=age, y=Vplasma_M, color="02_Vplasma")) +
#   geom_line(data = Variables_M, aes(x=age, y=Vliver_M, color="03_Vliver")) +
#   geom_line(data = Variables_M, aes(x=age, y=Vstomach_M, color="04_Vstomach")) +
#   geom_line(data = Variables_M, aes(x=age, y=Vgut_M, color="05_Vgut")) +
#   geom_line(data = Variables_M, aes(x=age, y=Vkidney_M, color="06_Vkidney")) +
#   geom_line(data = Variables_M, aes(x=age, y=Vskin_M, color="07_Vskin")) +
#   geom_line(data = Variables_M, aes(x=age, y=Vadipose_M, color="08_Vadipose")) +
#   geom_line(data = Variables_M, aes(x=age, y=Vbrain_M, color="09_Vbrain")) +
#   geom_line(data = Variables_M, aes(x=age, y=Vmuscle_M, color="10_Vmuscle")) +
#   geom_line(data = Variables_M, aes(x=age, y=Vtotal_M, color="11_Vtotal")) +
# 
#   scale_colour_manual(name='',
#                       values=c('01_BW'='black',
#                                '02_Vplasma'='red',
#                                '03_Vliver'='darkred',
#                                '04_Vstomach'='indianred',
#                                '05_Vgut'='lightpink',
#                                '06_Vkidney'='indianred4',
#                                '07_Vskin'='darksalmon',
#                                '08_Vadipose'='khaki',
#                                '09_Vbrain'='mistyrose',
#                                '10_Vmuscle'='rosybrown1',
#                                '11_Vtotal'='grey'),
#                       labels=c('01_BW'='Body weight',
#                                '02_Vplasma'='Plasma',
#                                '03_Vliver'='Liver',
#                                '04_Vstomach'='Stomach',
#                                '05_Vgut'='Gut',
#                                '06_Vkidney'='Kidney',
#                                '07_Vskin'='Skin',
#                                '08_Vadipose'='Adipose tissue',
#                                '09_Vbrain'='Brain',
#                                '10_Vmuscle'='Muscle tissue',
#                                '11_Vtotal'='Sum of all tissues')) +
# 
#   theme_bw() +
#   labs(title="Organ weight changes over time - males",x="\nAge (y)", y="Weight (kg)\n") +
#   theme(plot.title = element_text(hjust = 0.5))
# 
# Figure_M
# 
# Figure male organ fractions 
# Figure_fractions_M <- ggplot() +
#   geom_line(data = Variables_fractions_M, aes(x=age, y=VplasmaFraction_M, color="02_Vplasma")) +
#   geom_line(data = Variables_fractions_M, aes(x=age, y=VliverFraction_M, color="03_Vliver")) +
#   geom_line(data = Variables_fractions_M, aes(x=age, y=VstomachFraction_M, color="04_Vstomach")) +
#   geom_line(data = Variables_fractions_M, aes(x=age, y=VgutFraction_M, color="05_Vgut")) +
#   geom_line(data = Variables_fractions_M, aes(x=age, y=VkidneyFraction_M, color="06_Vkidney")) +
#   geom_line(data = Variables_fractions_M, aes(x=age, y=VskinFraction_M, color="07_Vskin")) +
#   geom_line(data = Variables_fractions_M, aes(x=age, y=VadiposeFraction_M, color="08_Vadipose")) +
#   geom_line(data = Variables_fractions_M, aes(x=age, y=VbrainFraction_M, color="09_Vbrain")) +
#   geom_line(data = Variables_fractions_M, aes(x=age, y=VmuscleFraction_M, color="10_Vmuscle")) +
#   geom_line(data = Variables_fractions_M, aes(x=age, y=VtotalFraction_M, color="11_Vtotal")) +
#   geom_line(data = Variables_M, aes(x=age, y=Fraction_M, color="12_Fraction")) +
# 
#   scale_colour_manual(name='',
#                       values=c('02_Vplasma'='red',
#                                '03_Vliver'='darkred',
#                                '04_Vstomach'='indianred',
#                                '05_Vgut'='lightpink',
#                                '06_Vkidney'='indianred4',
#                                '07_Vskin'='darksalmon',
#                                '08_Vadipose'='khaki',
#                                '09_Vbrain'='mistyrose',
#                                '10_Vmuscle'='rosybrown1',
#                                '11_Vtotal'='grey',
#                                '12_Fraction'='black'),
#                       labels=c('02_Vplasma'='Plasma',
#                                '03_Vliver'='Liver',
#                                '04_Vstomach'='Stomach',
#                                '05_Vgut'='Gut',
#                                '06_Vkidney'='Kidney',
#                                '07_Vskin'='Skin',
#                                '08_Vadipose'='Adipose tissue',
#                                '09_Vbrain'='Brain',
#                                '10_Vmuscle'='Muscle tissue',
#                                '11_Vtotal'='Sum of all tissues',
#                                '12_Fraction'='Sum tissues / BW')) +
# 
#   theme_bw() +
#   labs(title="Organ fractional volume changes over time - males",x="\nAge (y)", y="Fraction of body weight\n") +
#   theme(plot.title = element_text(hjust = 0.5))
# 
# Figure_fractions_M
# 
# Figure male blood flows 
# Figure_Q_M <- ggplot() +
#   geom_line(data = Variables_M, aes(x=age, y=CardOut_M, color="01_CardOut")) +
#   geom_line(data = Variables_M, aes(x=age, y=Qliver_M, color="03_Qliver")) +
#   geom_line(data = Variables_M, aes(x=age, y=Qstomach_M, color="04_Qstomach")) +
#   geom_line(data = Variables_M, aes(x=age, y=Qgut_M, color="05_Qgut")) +
#   geom_line(data = Variables_M, aes(x=age, y=Qkidney_M, color="06_Qkidney")) +
#   geom_line(data = Variables_M, aes(x=age, y=Qskin_M, color="07_Qskin")) +
#   geom_line(data = Variables_M, aes(x=age, y=Qadipose_M, color="08_Qadipose")) +
#   geom_line(data = Variables_M, aes(x=age, y=Qbrain_M, color="09_Qbrain")) +
#   geom_line(data = Variables_M, aes(x=age, y=Qmuscle_M, color="10_Qmuscle")) +
#   geom_line(data = Variables_M, aes(x=age, y=Qtotal_M, color="11_Qtotal")) +
# 
#   scale_colour_manual(name='',
#                       values=c('01_CardOut'='black',
#                                '03_Qliver'='darkred',
#                                '04_Qstomach'='indianred',
#                                '05_Qgut'='lightpink',
#                                '06_Qkidney'='indianred4',
#                                '07_Qskin'='darksalmon',
#                                '08_Qadipose'='khaki',
#                                '09_Qbrain'='mistyrose',
#                                '10_Qmuscle'='rosybrown1',
#                                '11_Qtotal'='grey'),
#                       labels=c('01_CardOut'='Cardiac output',
#                                '03_Qliver'='Liver',
#                                '04_Qstomach'='Stomach',
#                                '05_Qgut'='Gut',
#                                '06_Qkidney'='Kidney',
#                                '07_Qskin'='Skin',
#                                '08_Qadipose'='Adipose tissue',
#                                '09_Qbrain'='Brain',
#                                '10_Qmuscle'='Muscle tissue',
#                                '11_Qtotal'='Sum of all tissues')) +
# 
#   theme_bw() +
#   labs(title="Plasma flow changes over time - males",x="\nAge (y)", y="Plasma flow (L/d)\n") +
#   theme(plot.title = element_text(hjust = 0.5))
# 
# Figure_Q_M
# 
# Figure female organ volumes 
# Figure_F <- ggplot() +
#   geom_line(data = Variables_F, aes(x=age, y=BW_F, color="01_BW")) +
#   geom_line(data = Variables_F, aes(x=age, y=Vplasma_F, color="02_Vplasma")) +
#   geom_line(data = Variables_F, aes(x=age, y=Vliver_F, color="03_Vliver")) +
#   geom_line(data = Variables_F, aes(x=age, y=Vstomach_F, color="04_Vstomach")) +
#   geom_line(data = Variables_F, aes(x=age, y=Vgut_F, color="05_Vgut")) +
#   geom_line(data = Variables_F, aes(x=age, y=Vkidney_F, color="06_Vkidney")) +
#   geom_line(data = Variables_F, aes(x=age, y=Vskin_F, color="07_Vskin")) +
#   geom_line(data = Variables_F, aes(x=age, y=Vadipose_F, color="08_Vadipose")) +
#   geom_line(data = Variables_F, aes(x=age, y=Vbrain_F, color="09_Vbrain")) +
#   geom_line(data = Variables_F, aes(x=age, y=Vmuscle_F, color="10_Vmuscle")) +
#   geom_line(data = Variables_F, aes(x=age, y=Vtotal_F, color="11_Vtotal")) +
#   
#   scale_colour_manual(name='',
#                       values=c('01_BW'='black',
#                                '02_Vplasma'='red',
#                                '03_Vliver'='darkred',
#                                '04_Vstomach'='indianred',
#                                '05_Vgut'='lightpink',
#                                '06_Vkidney'='indianred4',
#                                '07_Vskin'='darksalmon',
#                                '08_Vadipose'='khaki',
#                                '09_Vbrain'='mistyrose',
#                                '10_Vmuscle'='rosybrown1',
#                                '11_Vtotal'='grey'),
#                       labels=c('01_BW'='Body weight',
#                                '02_Vplasma'='Plasma',
#                                '03_Vliver'='Liver',
#                                '04_Vstomach'='Stomach',
#                                '05_Vgut'='Gut',
#                                '06_Vkidney'='Kidney',
#                                '07_Vskin'='Skin',
#                                '08_Vadipose'='Adipose tissue',
#                                '09_Vbrain'='Brain',
#                                '10_Vmuscle'='Muscle tissue',
#                                '11_Vtotal'='Sum of all tissues')) +
#   
#   theme_bw() +
#   labs(title="Organ weight changes over time - females",x="\nAge (y)", y="Weight (kg)\n") +
#   theme(plot.title = element_text(hjust = 0.5))
# 
# Figure_F
# 
# Figure female organ fractions 
# Figure_fractions_F <- ggplot() +
#   geom_line(data = Variables_fractions_F, aes(x=age, y=VplasmaFraction_F, color="02_Vplasma")) +
#   geom_line(data = Variables_fractions_F, aes(x=age, y=VliverFraction_F, color="03_Vliver")) +
#   geom_line(data = Variables_fractions_F, aes(x=age, y=VstomachFraction_F, color="04_Vstomach")) +
#   geom_line(data = Variables_fractions_F, aes(x=age, y=VgutFraction_F, color="05_Vgut")) +
#   geom_line(data = Variables_fractions_F, aes(x=age, y=VkidneyFraction_F, color="06_Vkidney")) +
#   geom_line(data = Variables_fractions_F, aes(x=age, y=VskinFraction_F, color="07_Vskin")) +
#   geom_line(data = Variables_fractions_F, aes(x=age, y=VadiposeFraction_F, color="08_Vadipose")) +
#   geom_line(data = Variables_fractions_F, aes(x=age, y=VbrainFraction_F, color="09_Vbrain")) +
#   geom_line(data = Variables_fractions_F, aes(x=age, y=VmuscleFraction_F, color="10_Vmuscle")) +
#   geom_line(data = Variables_fractions_F, aes(x=age, y=VtotalFraction_F, color="11_Vtotal")) +
#   geom_line(data = Variables_F, aes(x=age, y=Fraction_F, color="12_Fraction")) +
#   
#   scale_colour_manual(name='',
#                       values=c('02_Vplasma'='red',
#                                '03_Vliver'='darkred',
#                                '04_Vstomach'='indianred',
#                                '05_Vgut'='lightpink',
#                                '06_Vkidney'='indianred4',
#                                '07_Vskin'='darksalmon',
#                                '08_Vadipose'='khaki',
#                                '09_Vbrain'='mistyrose',
#                                '10_Vmuscle'='rosybrown1',
#                                '11_Vtotal'='grey',
#                                '12_Fraction'='black'),
#                       labels=c('02_Vplasma'='Plasma',
#                                '03_Vliver'='Liver',
#                                '04_Vstomach'='Stomach',
#                                '05_Vgut'='Gut',
#                                '06_Vkidney'='Kidney',
#                                '07_Vskin'='Skin',
#                                '08_Vadipose'='Adipose tissue',
#                                '09_Vbrain'='Brain',
#                                '10_Vmuscle'='Muscle tissue',
#                                '11_Vtotal'='Sum of all tissues',
#                                '12_Fraction'='Sum tissues / BW')) +
#   
#   theme_bw() +
#   labs(title="Organ fractional volume changes over time - females",x="\nAge (y)", y="Fraction of body weight\n") +
#   theme(plot.title = element_text(hjust = 0.5))
# 
# Figure_fractions_F
# 
# Figure female blood flows 
# Figure_Q_F <- ggplot() +
#   geom_line(data = Variables_F, aes(x=age, y=CardOut_F, color="01_CardOut")) +
#   geom_line(data = Variables_F, aes(x=age, y=Qliver_F, color="03_Qliver")) +
#   geom_line(data = Variables_F, aes(x=age, y=Qstomach_F, color="04_Qstomach")) +
#   geom_line(data = Variables_F, aes(x=age, y=Qgut_F, color="05_Qgut")) +
#   geom_line(data = Variables_F, aes(x=age, y=Qkidney_F, color="06_Qkidney")) +
#   geom_line(data = Variables_F, aes(x=age, y=Qskin_F, color="07_Qskin")) +
#   geom_line(data = Variables_F, aes(x=age, y=Qadipose_F, color="08_Qadipose")) +
#   geom_line(data = Variables_F, aes(x=age, y=Qbrain_F, color="09_Qbrain")) +
#   geom_line(data = Variables_F, aes(x=age, y=Qmuscle_F, color="10_Qmuscle")) +
#   geom_line(data = Variables_F, aes(x=age, y=Qtotal_F, color="11_Qtotal")) +
#   
#   scale_colour_manual(name='',
#                       values=c('01_CardOut'='black',
#                                '03_Qliver'='darkred',
#                                '04_Qstomach'='indianred',
#                                '05_Qgut'='lightpink',
#                                '06_Qkidney'='indianred4',
#                                '07_Qskin'='darksalmon',
#                                '08_Qadipose'='khaki',
#                                '09_Qbrain'='mistyrose',
#                                '10_Qmuscle'='rosybrown1',
#                                '11_Qtotal'='grey'),
#                       labels=c('01_CardOut'='Cardiac output',
#                                '03_Qliver'='Liver',
#                                '04_Qstomach'='Stomach',
#                                '05_Qgut'='Gut',
#                                '06_Qkidney'='Kidney',
#                                '07_Qskin'='Skin',
#                                '08_Qadipose'='Adipose tissue',
#                                '09_Qbrain'='Brain',
#                                '10_Qmuscle'='Muscle tissue',
#                                '11_Qtotal'='Sum of all tissues')) +
#   
#   theme_bw() +
#   labs(title="Plasma flow changes over time - females",x="\nAge (y)", y="Plasma flow (L/d)\n") +
#   theme(plot.title = element_text(hjust = 0.5))
# 
# Figure_Q_F
# 
# Plot mass balance organ flow  
# 
# # !!!!! THERE SEEMS TO BE AN ISSUE !!!!!
# FlowMassBalance <- Variables_df %>% select(age, Qmass_balance_M, Qmass_balance_F) %>% 
#   pivot_longer(names_to = "Gender", values_to = "MB", Qmass_balance_M:Qmass_balance_F) %>% 
#   ggplot(aes(age, MB)) +
#   geom_path() +
#   facet_wrap(~Gender)
# FlowMassBalance
# 
# # Plot mass balance organ volumes
# # !!!!! THERE SEEMS TO BE AN ISSUE !!!!!
# VolumeMassBalance <- Variables_df %>% select(age, Vmass_balance_M, Vmass_balance_F) %>% 
#   pivot_longer(names_to = "Gender", values_to = "MB", Vmass_balance_M:Vmass_balance_F) %>% 
#   ggplot(aes(age, MB)) +
#   geom_path() +
#   facet_wrap(~Gender)
# VolumeMassBalance
# 
# # Plot mass balance organ volumes
# # !!!!! THERE SEEMS TO BE AN ISSUE !!!!!
# VolumeTest <- Variables_df %>% select(age, VtotalFraction_M, VtotalFraction_F) %>% 
#   pivot_longer(names_to = "Gender", values_to = "Volume", VtotalFraction_M:VtotalFraction_F) %>% 
#   ggplot(aes(age, Volume, color = Gender)) +
#   geom_path() 
# VolumeTest




## Connecting the dose/flow/volume parameter to a timepoint  ####

### Male ----
### Doses
varOraldose_M <- approxfun(Variables_df$TIME, Variables_df$Oraldose_M, rule = 2)
varOraldose_F <- approxfun(Variables_df$TIME, Variables_df$Oraldose_F, rule = 2)
varDermaldose_M <- approxfun(Variables_df$TIME, Variables_df$Dermaldose_M, rule = 2)
varDermaldose_F <- approxfun(Variables_df$TIME, Variables_df$Dermaldose_F, rule = 2)

### Volumes - Male
varBW_M <- approxfun(Variables_df$TIME, Variables_df$BW_M, rule = 2)
varVplasma_M <- approxfun(Variables_df$TIME, Variables_df$Vplasma_M, rule = 2)
varVliver_M <- approxfun(Variables_df$TIME, Variables_df$Vliver_M, rule = 2)
varVstomach_M <- approxfun(Variables_df$TIME, Variables_df$Vstomach_M, rule = 2)
varVgut_M <- approxfun(Variables_df$TIME, Variables_df$Vgut_M, rule = 2)
varVkidney_M <- approxfun(Variables_df$TIME, Variables_df$Vkidney_M, rule = 2)
varVskin_M <- approxfun(Variables_df$TIME, Variables_df$Vskin_M, rule = 2)
varVadipose_M <- approxfun(Variables_df$TIME, Variables_df$Vadipose_M, rule = 2)
varVbrain_M <- approxfun(Variables_df$TIME, Variables_df$Vbrain_M, rule = 2)
varVmuscle_M <- approxfun(Variables_df$TIME, Variables_df$Vmuscle_M, rule = 2)
varVtotal_M <- approxfun(Variables_df$TIME, Variables_df$Vtotal_M, rule = 2)
varVurinarytract_M <- approxfun(Variables_df$TIME, Variables_df$Vurinarytract_M, rule = 2)
# varVart_M <- approxfun(Variables_df$TIME, Variables_df$Vart_M, rule = 2) #we have a central plasma compartment 
# varVven_M <- approxfun(Variables_df$TIME, Variables_df$Vven_M, rule = 2) #we have a central plasma compartment
varVadrenal_M <- approxfun(Variables_df$TIME, Variables_df$Vadrenal_M, rule = 2)
varVbone_M <- approxfun(Variables_df$TIME, Variables_df$Vbone_M, rule = 2)
varVbonenonperfused_M <- approxfun(Variables_df$TIME, Variables_df$Vbonenonperfused_M, rule = 2)
varVbreast_M <- approxfun(Variables_df$TIME, Variables_df$Vbreast_M, rule = 2)
varVheart_M <- approxfun(Variables_df$TIME, Variables_df$Vheart_M, rule = 2)
varVmarrow_M <- approxfun(Variables_df$TIME, Variables_df$Vmarrow_M, rule = 2)
varVrepro_M <- approxfun(Variables_df$TIME, Variables_df$Vrepro_M, rule = 2)
varVpancreas_M <- approxfun(Variables_df$TIME, Variables_df$Vpancreas_M, rule = 2)
varVspleen_M <- approxfun(Variables_df$TIME, Variables_df$Vspleen_M, rule = 2)
varVthyroid_M <- approxfun(Variables_df$TIME, Variables_df$Vthyroid_M, rule = 2)
varVlung_M <- approxfun(Variables_df$TIME, Variables_df$Vlung_M, rule = 2)

### CHANGED
varVrest_M <- approxfun(Variables_df$TIME, Variables_df$Vrest_M, rule = 2)
varVSkb_M <- approxfun(Variables_df$TIME, Variables_df$VSkb_M, rule = 2)
varSkbTarea_M <- approxfun(Variables_df$TIME, Variables_df$SkbTarea_M, rule = 2)
varVkidneyBlood_M <- approxfun(Variables_df$TIME, Variables_df$VkidneyBlood_M, rule = 2)
varVkidneyTissue_M <- approxfun(Variables_df$TIME, Variables_df$VkidneyTissue_M, rule = 2)
varVFil_M <- approxfun(Variables_df$TIME, Variables_df$VFil_M, rule = 2)


## Blood flows - Male
varCardOut_M <- approxfun(Variables_df$TIME, Variables_df$CardOut_M, rule = 2)
varQliver_M <- approxfun(Variables_df$TIME, Variables_df$Qliver_M, rule = 2)
varQstomach_M <- approxfun(Variables_df$TIME, Variables_df$Qstomach_M, rule = 2)
varQgut_M <- approxfun(Variables_df$TIME, Variables_df$Qgut_M, rule = 2)
varQkidney_M <- approxfun(Variables_df$TIME, Variables_df$Qkidney_M, rule = 2)
varQskin_M <- approxfun(Variables_df$TIME, Variables_df$Qskin_M, rule = 2)
varQadipose_M <- approxfun(Variables_df$TIME, Variables_df$Qadipose_M, rule = 2)
varQbrain_M <- approxfun(Variables_df$TIME, Variables_df$Qbrain_M, rule = 2)
varQmuscle_M <- approxfun(Variables_df$TIME, Variables_df$Qmuscle_M, rule = 2)
varQtotal_M <- approxfun(Variables_df$TIME, Variables_df$Qtotal_M, rule = 2)
varQurinarytract_M <- approxfun(Variables_df$TIME, Variables_df$Qurinarytract_M, rule = 2)
varQadrenal_M <- approxfun(Variables_df$TIME, Variables_df$Qadrenal_M, rule = 2)
varQbone_M <- approxfun(Variables_df$TIME, Variables_df$Qbone_M, rule = 2)
varQbreast_M <- approxfun(Variables_df$TIME, Variables_df$Qbreast_M, rule = 2)
varQheart_M <- approxfun(Variables_df$TIME, Variables_df$Qheart_M, rule = 2)
varQmarrow_M <- approxfun(Variables_df$TIME, Variables_df$Qmarrow_M, rule = 2)
varQrepro_M <- approxfun(Variables_df$TIME, Variables_df$Qrepro_M, rule = 2)
varQpancreas_M <- approxfun(Variables_df$TIME, Variables_df$Qpancreas_M, rule = 2)
varQspleen_M <- approxfun(Variables_df$TIME, Variables_df$Qspleen_M, rule = 2)
varQthyroid_M <- approxfun(Variables_df$TIME, Variables_df$Qthyroid_M, rule = 2)
varQlung_M <- approxfun(Variables_df$TIME, Variables_df$Qlung_M, rule = 2)

### CHANGED
varQhepatic_M <- approxfun(Variables_df$TIME, Variables_df$Qhepatic_M, rule = 2)
varQrest_M <- approxfun(Variables_df$TIME, Variables_df$Qrest_M, rule = 2)

## Clearances - Male
varCLdermalabs_M <- approxfun(Variables_df$TIME, Variables_df$CLdermalabs_M, rule = 2)
varVmax_M <- approxfun(Variables_df$TIME, Variables_df$Vmax_M, rule = 2)
varCLfecal_M <- approxfun(Variables_df$TIME, Variables_df$CLfecal_M, rule = 2)
varCLbiliary_M <- approxfun(Variables_df$TIME, Variables_df$CLbiliary_M, rule = 2)
varQUr_M <- approxfun(Variables_df$TIME, Variables_df$QUr_M, rule = 2)
#varkUr_M <- approxfun(Variables_df$TIME, Variables_df$kUr_M, rule = 2)
varGFR_M <- approxfun(Variables_df$TIME, Variables_df$GFR_M)  #L/d ; as it's 18% of total renal plasma flow [ICRP 89] http://www.icrp.org/publication.asp?id=ICRP%20Publication%2089
varKW_cortex_M <- approxfun(Variables_df$TIME, Variables_df$KW_cortex)  # g kidney cortexes, only scaling to kidney cortex volume as proximal tubule cells are in the cortex; 70% of the total kidney volume according to ICRP89; PT are in the cortex https://doi.org/10.1021/acs.molpharmaceut.4c00504; alternatively we could have 68% of kidney weight https://doi.org/10.1124/dmd.117.075242 
varCL_PltPT_M <- approxfun(Variables_df$TIME, Variables_df$CL_PltPT)  # to discuss
varCL_FiltPT_M <- approxfun(Variables_df$TIME, Variables_df$CL_FiltPT) #proximal tubule to filtrate


### Female ----

## Volumes - Female
varBW_F <- approxfun(Variables_df$TIME, Variables_df$BW_F, rule = 2)
varVplasma_F <- approxfun(Variables_df$TIME, Variables_df$Vplasma_F, rule = 2)
varVliver_F <- approxfun(Variables_df$TIME, Variables_df$Vliver_F, rule = 2)
varVstomach_F <- approxfun(Variables_df$TIME, Variables_df$Vstomach_F, rule = 2)
varVgut_F <- approxfun(Variables_df$TIME, Variables_df$Vgut_F, rule = 2)
varVkidney_F <- approxfun(Variables_df$TIME, Variables_df$Vkidney_F, rule = 2)
varVskin_F <- approxfun(Variables_df$TIME, Variables_df$Vskin_F, rule = 2)
varVadipose_F <- approxfun(Variables_df$TIME, Variables_df$Vadipose_F, rule = 2)
varVbrain_F <- approxfun(Variables_df$TIME, Variables_df$Vbrain_F, rule = 2)
varVmuscle_F <- approxfun(Variables_df$TIME, Variables_df$Vmuscle_F, rule = 2)
varVtotal_F <- approxfun(Variables_df$TIME, Variables_df$Vtotal_F, rule = 2)
varVurinarytract_F <- approxfun(Variables_df$TIME, Variables_df$Vurinarytract_F, rule = 2)
# varVart_F <- approxfun(Variables_df$TIME, Variables_df$Vart_F, rule = 2) #we have a central plasma compartment 
# varVven_F <- approxfun(Variables_df$TIME, Variables_df$Vven_F, rule = 2) #we have a central plasma compartment
varVadrenal_F <- approxfun(Variables_df$TIME, Variables_df$Vadrenal_F, rule = 2)
varVbone_F <- approxfun(Variables_df$TIME, Variables_df$Vbone_F, rule = 2)
varVbonenonperfused_F <- approxfun(Variables_df$TIME, Variables_df$Vbonenonperfused_F, rule = 2)
varVbreast_F <- approxfun(Variables_df$TIME, Variables_df$Vbreast_F, rule = 2)
varVheart_F <- approxfun(Variables_df$TIME, Variables_df$Vheart_F, rule = 2)
varVmarrow_F <- approxfun(Variables_df$TIME, Variables_df$Vmarrow_F, rule = 2)
varVrepro_F <- approxfun(Variables_df$TIME, Variables_df$Vrepro_F, rule = 2)
varVpancreas_F <- approxfun(Variables_df$TIME, Variables_df$Vpancreas_F, rule = 2)
varVspleen_F <- approxfun(Variables_df$TIME, Variables_df$Vspleen_F, rule = 2)
varVthyroid_F <- approxfun(Variables_df$TIME, Variables_df$Vthyroid_F, rule = 2)
varVlung_F <- approxfun(Variables_df$TIME, Variables_df$Vlung_F, rule = 2)
### CHANGED
varVrest_F <- approxfun(Variables_df$TIME, Variables_df$Vrest_F, rule = 2)
varVSkb_F <- approxfun(Variables_df$TIME, Variables_df$VSkb_F, rule = 2)
varSkbTarea_F <- approxfun(Variables_df$TIME, Variables_df$SkbTarea_F, rule = 2)
varVkidneyBlood_F <- approxfun(Variables_df$TIME, Variables_df$VkidneyBlood_F, rule = 2)
varVkidneyTissue_F <- approxfun(Variables_df$TIME, Variables_df$VkidneyTissue_F, rule = 2)
varVFil_F <- approxfun(Variables_df$TIME, Variables_df$VFil_F, rule = 2)

## Blood flows - Female
varCardOut_F <- approxfun(Variables_df$TIME, Variables_df$CardOut_F, rule = 2)
varQliver_F <- approxfun(Variables_df$TIME, Variables_df$Qliver_F, rule = 2)
varQstomach_F <- approxfun(Variables_df$TIME, Variables_df$Qstomach_F, rule = 2)
varQgut_F <- approxfun(Variables_df$TIME, Variables_df$Qgut_F, rule = 2)
varQkidney_F <- approxfun(Variables_df$TIME, Variables_df$Qkidney_F, rule = 2)
varQskin_F <- approxfun(Variables_df$TIME, Variables_df$Qskin_F, rule = 2)
varQadipose_F <- approxfun(Variables_df$TIME, Variables_df$Qadipose_F, rule = 2)
varQbrain_F <- approxfun(Variables_df$TIME, Variables_df$Qbrain_F, rule = 2)
varQmuscle_F <- approxfun(Variables_df$TIME, Variables_df$Qmuscle_F, rule = 2)
varQtotal_F <- approxfun(Variables_df$TIME, Variables_df$Qtotal_F, rule = 2)
varQurinarytract_F <- approxfun(Variables_df$TIME, Variables_df$Qurinarytract_F, rule = 2)
varQadrenal_F <- approxfun(Variables_df$TIME, Variables_df$Qadrenal_F, rule = 2)
varQbone_F <- approxfun(Variables_df$TIME, Variables_df$Qbone_F, rule = 2)
varQbreast_F <- approxfun(Variables_df$TIME, Variables_df$Qbreast_F, rule = 2)
varQheart_F <- approxfun(Variables_df$TIME, Variables_df$Qheart_F, rule = 2)
varQmarrow_F <- approxfun(Variables_df$TIME, Variables_df$Qmarrow_F, rule = 2)
varQrepro_F <- approxfun(Variables_df$TIME, Variables_df$Qrepro_F, rule = 2)
varQpancreas_F <- approxfun(Variables_df$TIME, Variables_df$Qpancreas_F, rule = 2)
varQspleen_F <- approxfun(Variables_df$TIME, Variables_df$Qspleen_F, rule = 2)
varQthyroid_F <- approxfun(Variables_df$TIME, Variables_df$Qthyroid_F, rule = 2)
varQlung_F <- approxfun(Variables_df$TIME, Variables_df$Qlung_F, rule = 2)
### CHANGED
varQhepatic_F <- approxfun(Variables_df$TIME, Variables_df$Qhepatic_F, rule = 2)
varQrest_F <- approxfun(Variables_df$TIME, Variables_df$Qrest_F, rule = 2)

## Clearances - Female
varCLdermalabs_F <- approxfun(Variables_df$TIME, Variables_df$CLdermalabs_F, rule = 2)
varVmax_F <- approxfun(Variables_df$TIME, Variables_df$Vmax_F, rule = 2)
varCLfecal_F <- approxfun(Variables_df$TIME, Variables_df$CLfecal_F, rule = 2)
varCLbiliary_F <- approxfun(Variables_df$TIME, Variables_df$CLbiliary_F, rule = 2)
varQUr_F <- approxfun(Variables_df$TIME, Variables_df$QUr_F, rule = 2)
#varkUr_F <- approxfun(Variables_df$TIME, Variables_df$kUr_F, rule = 2)
varGFR_F <- approxfun(Variables_df$TIME, Variables_df$GFR_F)  #L/d ; as it's 18% of total renal plasma flow [ICRP 89] http://www.icrp.org/publication.asp?id=ICRP%20Publication%2089
varKW_cortex_F <- approxfun(Variables_df$TIME, Variables_df$KW_cortex)  # this is male for now; need to change KW_cortex to Vkidney_F for Female !!! g kidney cortexes, only scaling to kidney cortex volume as proximal tubule cells are in the cortex; 70% of the total kidney volume according to ICRP89; PT are in the cortex https://doi.org/10.1021/acs.molpharmaceut.4c00504; alternatively we could have 68% of kidney weight https://doi.org/10.1124/dmd.117.075242 
varCL_PltPT_F <- approxfun(Variables_df$TIME, Variables_df$CL_PltPT)  # to discuss
varCL_FiltPT_F <- approxfun(Variables_df$TIME, Variables_df$CL_FiltPT) #proximal tubule to filtrate
varCL_passive <- approxfun(Variables_df$TIME, Variables_df$CL_passive) #proximal tubule to filtrate


# ---------------------------------------------------------------------------- #


## FINAL PARAMS ----
parms <- c(fup, PG, PL, PF, PK, PSk, PR, kAap, kAbl) #Km

  
# PBK MODEL ####
# ---------------------------------------------------------------------------- #

### Male ####
PBPKmodPFOA_M <- function(t, state, parameters){
  with(as.list(c(state, parameters)), {
    
    ### Time dependent variables ----
    Oraldose <- varOraldose_M(t)
    Dermaldose <- varDermaldose_M(t)
    
    # Clearances
    CLdermalabs <- varCLdermalabs_M(t)
    # Vmax <- varVmax_M(t)
    CLfaeces <- varCLfecal_M(t)
    CLbiliary <- varCLbiliary_M(t)
    QUr <- varQUr_M(t) #Urine secretion flow rate NEEDS TO BE CHANGED WHEN WE HAVE THE CORRECT VALUE. Comment Chrysa 05-11-2024: I believe this value is ok, it's from ICRP89
    #kUr <- varkUr_M(t)
    CL_PltPT <- varCL_PltPT_M(t)
    CL_FiltPT <- varCL_FiltPT_M(t)
    CL_passive <- varCL_passive(t)
    
    # Flow rates 
    QCP <- varCardOut_M(t) #cardiac output plasma
    QG <- varQgut_M(t) #gut (=intestine only)
    QH <- varQhepatic_M(t) #portal vein except gut
    QL <- varQliver_M(t) #vena cava
    QF <- varQadipose_M(t) #adipose (used to be called QF)
    QK <- varQkidney_M(t) #kidney
    # used to be QFil <- varQurinarytract_M(t) #filtrate #Comment Chrysa 05-11-2024 this is GFR now
    GFR <- varGFR_M(t) #filtrate
    QSk <- varQskin_M(t) #skin
    QR <- varQrest_M(t) #rest of the body
    #Qp <- varQp_M(t) # needs to be incorporated when inhalation exposure is included
    
    
    # Volumes
    VP <- varVplasma_M(t) #plasma
    VG <- varVgut_M(t) #gut
    VL <- varVliver_M(t) #liver
    VF <- varVadipose_M(t) #adipose (used to be called VF)
    VK <- varVkidney_M(t) #kidney
    # VKT <- varVkidneyTissue_M(t)
    # VKB <- varVkidneyBlood_M(t)
    VFil <- varVFil_M(t) #filtrate
    VSk <- varVskin_M(t) #skin
    VSkb <- varVSkb_M(t) #skin
    VR <- varVrest_M(t) #rest of the body
    VUr <- 1 # QUr*Timeframe #NEED TO DEFINE THE TIMEFRAME
    
    
    
    ### Concentrations ----
    
    # Organ concentrations (ug/L); these are TOTAL concentrations
    CP <- AP/VP  # concentration in plasma (ug/L)
    CG <- AG/VG  # concentration of PFOA in gut (ug/L)
    CL <- AL/VL  # concentration of PFOA in liver (ug/L)
    CF <- AF/VF  # concentration of PFOA in adipose (ug/L)
    CFil <- AFil/VFil # concentration of PFOA in filtrate compartment
    CK <- AK/VK  # concentration of PFOA in kidney (ug/L)
    # NEW
    CUr <- AUr/VUr # concentration of PFOA in urine (ug/L)
    CSk <- ASk/VSk  # concentration of PFOA in skin (ug/L)
    CSkb <- ASkb/VSkb # concentration of PFOA in skin barrier (ug/L)
    CR <- AR/VR  # concentration of PFOA in rest of the body (ug/L)
    # CKB <- AKB/VKB # concentration in kidney blood
    # CKT <- AKT/VKT # concentration in kidney tissue
    
    # Venous concentrations (ug/L); these are the concentrations leaving the organs
    CVG <- CG/PG # concentration of PFOA leaving gut (ug/L)
    CVL <- CL/PL  # concentration of PFOA leaving liver (ug/L)
    CVF <- CF/PF # concentration of PFOA leaving adipose (ug/L)
    CVK <- CK/PK  # concentration of PFOA leaving the kidney (ug/L)
    CVSk <- CSk/PSk  # Concentration of PFOA leaving the skin (ug/L)
    # NEW
    CVR <- CR/PR  # concentration of PFOA leaving the rest of the body (ug/L)
    # CVKB <- CKB/PK # concentration of PFOA leaving kidney blood (ug/L)
    # 
    
    
    ### Differential equations ----
    # Fat compartment
    dAF <- QF*(CP-CVF) # (ug/h)
    
    # Rest compartment
    dAR <- QR*(CP-CVR) # (ug/h)
    
    #### NEW/CHANGED---------------
    
    # Gut compartment: plasma to gut then to liver. Biliary clearance to gut, based on free concentration in the liver; faecal clearance based on total gut concentration; not only the free fraction as partitioning is not needed
    dAG <- Oraldose + QG*CP - QG*CVG + CLbiliary*CL*fup - CLfaeces*CG #(ug/h)
    
    # Excretion fecal: cumulative
    dAEx_feces <- CLfaeces*CG #ug/h, used to be CVG but I think this is wrong as it's not the tissue plasma partition that defines fecal excretion
    
    # Liver compartment
    dAL <- QH*CP + QG*CVG - (QH+QG)*CVL - CLbiliary*CL*fup #input from the hepatic artery
    # Rate of PFOA amount change in the liver (ug/h)
    
    # Comment Chrysa 04-11-2024: see updated kidney compartment below
    #### OLD Kidney compartment ----
    # ### Kidney: this is basically the kidney blood compartment
    # ### I think we're missing active secretion
    # dAK <- QK*(CP-CVK) + (Vmax*CFil)/(Km+CFil) - GFR*CK*fup
    # #                    re-absorption     ultrafiltration
    # 
    # ### Filtrate compartment: rate of formation of the filtrate(=urine) in the lumen; this is the urinary tract from Aude's lifestage equations
    # dAFil <- GFR*CK*fup - (Vmax*CFil)/(Km+CFil) - QUr*CFil
    # #        ultrafiltration  REABSORPTION  urine flow rate to the bladder
    # 
    # ### Remove this compartment
    # ### Urine: urine in the bladder
    # ### ??? why do we need this compartment? I would try having a simple compartment first; i.e not having an excretion compartment
    # dAUr <- QUr*CFil - kUr*AUr #(ug/h)
    # #       urineflow  urine excretion
    # 
    # ### Excretion urinary: cumulative PFOA concentration in the urine
    # dAEx_urine <- kUr*AUr #(ug/h)
    
    #### Updated Kidney compartment ----
    # Kidney Blood
    # dAKB <- QK*(CP-CVKB) - GFR*fup*CKB #- CL_PltPT*fup* (CKB - (CKT*kAbl)) #- CL_passive*fup*(CKB - CKT) #  
    
    dAK <- QK*(CP-CVK) - GFR*fup*CK + CL_FiltPT*fup*CFil #- (CVKB*kAap)) #- CL_PltPT*fup* (CKB - (CKT*kAbl)) #- CL_passive*fup*(CKB - CKT) #  
    
    
    # Filtrate
    dAFil <- GFR*fup*CK - CL_FiltPT*fup*CFil - QUr*CFil  #- (CKT*kAap)) - QUr*CFil #
    
    # Urine
    dAUr <- QUr*CFil
    
    # # Kidney Tissue
    # dAKT <-  CL_FiltPT*fup*(CFil- (CKT*kAap)) #- (CKT*kAbl)), + CL_passive*fup*(CKB - CKT) CL_PltPT*fup*CKB
    # 
    # Skin compartment
    # Skin barrier
    dASkb <- Dermaldose - CLdermalabs*CSkb #(ug/h)
    # Skin tissue
    dASk <- CLdermalabs*CSkb + QSk*(CP - CVSk)
    
    # Plasma compartment
    dAP <- - (QSk + QG + QH + QF + QK + QR)*CP +
      QSk*CVSk + (QH+QG)*CVL + QF*CVF + QK*CVK + QR*CVR #(ug/h)
    
    Atot <- AP + ASkb + ASk + AG + AL + AF + AK  + AFil + AUr + AR + AEx_feces # + AEx_urine + AKT
    dInput <- Oraldose + Dermaldose
    MB = Input - Atot
    
    list(c(dAP, dASkb, dASk, dAG, dAL, dAF, dAK, dAFil, dAUr, dAR, dAEx_feces, dInput), #dAEx_urine dAKT
         c(CP=CP, 
           CSkb=CSkb, CSk=CSk,
           CG=CG, CVG=CVG, 
           CL=CL, CVL=CVL, 
           CF=CF, CVF=CVF,
           CK=CK, #CKT=CKT, 
           CFil=CFil, #CVKB=CVKB,
           CR=CR, CVR=CVR,
           Atot=Atot,MB=MB))
    
  }
  )
}

### Initials ####

A_init <- c(AP=0, 
            ASkb=0, 
            ASk=0, 
            AG=0, 
            AL=0, 
            AF=0, 
            AK=0,
            # AKT=0,
            AFil=0, 
            AUr=0, 
            AR=0, 
            AEx_feces=0, 
            #AEx_urine=0,
            Input=0)


### Solving the model ####
output_PFOA <- lsoda(A_init, TIME, PBPKmodPFOA_M, parms)
output.PFOA.df <- as.data.frame(output_PFOA)
output.df <- Variables_df %>%
  rename(time = TIME) %>%
  left_join(output.PFOA.df)


# RESULTS ####
# ---------------------------------------------------------------------------- #

### Figures ####
# Plot_PFOA_doses <- ggplot()+
#   geom_line(data = Variables_df, aes(x = age, y = Oraldose_M),col="red")+
#   geom_line(data = Variables_df, aes(x = age, y = Dermaldose_M),col="blue")+
#   CP_theme()+
#   theme(axis.text.x = element_text(size = 7),axis.text.y = element_text(size = 7), axis.title = element_text(size = 8))+
#   #scale_colour_hue()+
#   ylab("Dose (ug)")
# 
# Plot_PFOA_doses


# PFOA in plasma of one individual
Plot_PFOA_Plasma <- ggplot()+
  geom_path(data = output.df, aes(x = age, y = CP))+
  theme(axis.text.x = element_text(size = 7),axis.text.y = element_text(size = 7), axis.title = element_text(size = 8))+
  ylab("Plasma (ng/ml)")

Plot_PFOA_Plasma
ggsave("PlasmaConcentration.1.png", dpi = 300)


# PFOA in Kidney of one individual
Plot_PFOA_Kidney <- ggplot()+
  geom_path(data = output.df, aes(x = age, y = CKB, color="Kidney blood"))+
  geom_path(data = output.df, aes(x = age, y = CKT, color="Kidney tissue"))+
  geom_path(data = output.df, aes(x = age, y = CFil, color="Filtrate"))+
  theme(axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7), 
        axis.title = element_text(size = 8)) +
  ylab("Kidney Compartments (ng/ml)") +
  scale_color_manual(values = c("Kidney blood" = "lightblue", 
                                "Kidney tissue" = "blue", 
                                "Filtrate" = "darkblue"))

Plot_PFOA_Kidney
ggsave("KidneyConcentration.2.png", dpi = 300)
