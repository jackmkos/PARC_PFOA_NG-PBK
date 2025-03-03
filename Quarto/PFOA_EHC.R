# --------------------------------------------------------------------------- #
# PBK MODEL FOR PFOA, TO BE USED TOGETHER WITH THE LATEST HBM DATA
# Model File
# CP, 27-01-2025
# --------------------------------------------------------------------------- #

rm(list=ls()) # to clear out the global environment

# Set working directory
HOME <- "C:/Users/pacho003/OneDrive - Wageningen University & Research/CP_L_R/PARC_PFAS_PBPKmodel/Codes_PFAS_models/Chrysanthi_Pachoulide_2024_PFOA_mechanistic"
# HOME <- "/home/westerj"
setwd(HOME)


# Set output storage directory
workingtime <- gsub(":", "-", Sys.time())
OUTPUT <- file.path("Output", Sys.Date(), workingtime)
dir.create(OUTPUT, recursive = TRUE)
setwd(OUTPUT)


# Load packages

library(ggplot2)
library(deSolve)
library(tidyverse)
library(PKNCA)
library(pracma)
library(readxl)

# For plotting
library(showtext)
font_add(family = "Garamond", regular = "GARA.TTF")
showtext_auto()
theme_CP <- function() {
  theme_bw()+
    theme(
      text = element_text(size = 25, lineheight = unit(0.5, "lines")), # lineheight is adjusting the space between lines
      axis.title = element_text(size = 20),
      axis.text = element_text(size = 20),
      axis.line = element_line(colour = "grey"),
      axis.ticks = element_line(colour = "grey"),
      # plot.margin = margin(0.2, 0.2, 0.2, 0.2, "cm"),
      panel.border = element_blank(), 
      panel.background = element_rect((fill = "white")),
      panel.grid = element_line(linewidth = 0.1,5, colour = "grey"), 
      strip.background = element_blank(),
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.box.margin = margin(0, 0, 0, 0, "cm"),
      legend.key.width = unit(0.2, "cm"),  
      legend.key.height = unit(0.2, "cm"),
      legend.text = element_text(size = 15)
    )
}


# INPUT ####
# ------------------------------------------------------ #

# Input variables
# Calculated variables from: "PARC PFOA PBK input: file ../Script/PARC PFOA PBK input.R
# Final_variables_M_df <- read_csv("C:/Users/pacho003/OneDrive - Wageningen University & Research/CP_L_R/PARC_PFAS_PBPKmodel/Codes_PFAS_models/Chrysanthi_Pachoulide_2024_PFOA_mechanistic/Input/2024-11-18/Final_variables_M.csv") %>% as.data.frame()
# Final_variables_M_df <- read_csv("C:/Users/pacho003/OneDrive - Wageningen University & Research/CP_L_R/PARC_PFAS_PBPKmodel/Codes_PFAS_models/Chrysanthi_Pachoulide_2024_PFOA_mechanistic/Input/2024-11-18/Final_variables_MVolunteer.csv")
Physio_params <- read_csv("C:/Users/pacho003/OneDrive - Wageningen University & Research/CP_L_R/PARC_PFAS_PBPKmodel/Codes_PFAS_models/Chrysanthi_Pachoulide_2024_PFOA_mechanistic/Input/PhysioVariables.csv")
PFOA_params <- read_csv("C:/Users/pacho003/OneDrive - Wageningen University & Research/CP_L_R/PARC_PFAS_PBPKmodel/Codes_PFAS_models/Chrysanthi_Pachoulide_2024_PFOA_mechanistic/Input/PFOAParams.csv")


# EXPOSURE SCENARIO ####
# ------------------------------------------------------ #

# exposure_stop <-  50*365      # days
sim_stop <- 450                 # days Abraham study follow-up time  

TSTART <- 0
TSTOP <- sim_stop               # days
DT <- 1/10                      # days
TIME <- seq(TSTART,TSTOP,by=DT)


# Oralconc <- 0.000187 # ug/kg/day [EFSA 2020; page 143] 3.96ug
# Dermconc <- 0.000542 # ug/kg/day; mean of #as.numeric(SumExpPFOA_LB_val[i,14])


Oraldose <- 3.96 # ug https://doi.org/10.1016/j.envint.2024.109047 # Oralconc*BDW # ug/day


# PBK MODEL PARAMETERS ####
# ---------------------------------------------------------------------------- #

## Constants ####

### Physiological  ####
Physio_params <- Physio_params %>% 
  filter(age == 60) %>%  # 60 years old adult, male, as in the Abraham study
  select(ends_with('_M')) %>% 
  mutate(BloodFlowSum = rowSums(select(., starts_with("Q_")))) %>% # 0.9935; total blood flow as the sum of the fractional blood flows of all organs on which we have data
  mutate(VolumesSum = rowSums(select(., starts_with("V_")))) # 0.96; total volume as the sum of the fractional organ volumes of all organs on which we have data

BW <- 82 #Physio_params$BDW_M 
QC <- Physio_params$CardOut_M # This is corrected for hematocrit already so it's plasma


# Fractional organ volumes
VAc <- Physio_params$V_adiposeFraction_M   
MassAc <- Physio_params$AdiposeMass_M      
VGc <- Physio_params$V_gutFraction_M
VKc <- Physio_params$V_kidneyFraction_M
VLc <- Physio_params$V_liverFraction_M
VPc <- Physio_params$V_plasmaFraction_M
Hct <- 46.7/100 #Physio_params$Hct_M
VSkc <- Physio_params$V_skinFraction_M
VTotc <- Physio_params$VolumesSum 

# Fractional organ blood flows
QAc <- Physio_params$Q_adiposeFraction_M/Physio_params$BloodFlowSum
QGc <- Physio_params$Q_gutFraction_M/Physio_params$BloodFlowSum 
QKc <- Physio_params$Q_kidneyFraction_M/Physio_params$BloodFlowSum 
QLc <- Physio_params$Q_liverFraction_M/Physio_params$BloodFlowSum 
QSkc <- Physio_params$Q_skinFraction_M/Physio_params$BloodFlowSum 


#### Mechanistic Kidney Model ####

# Fractional volumes recalculated from Pletz et al. 2020 https://doi.org/10.1016/j.comtox.2021.100172
fVGlom <- 0.0400  #Fractional volume of glomeruli to total kidney volume
fVGlomB <- 0.6815 #Fractional volume of glomerular blood to the total glomeruli volume
fVGlomL <- 0.3185 #Fractional volume of glomerular space(lumen) to the total glomeruli volume

fVPT <- 0.3581  #Fractional volume of proximal tubule to total kidney volume
fVPTB <- 0.0766 #Fractional volume of proximal tubule blood to proximal tubule volume
fVPTC <- 0.5110 #Fractional volume of proximal tubule cell to proximal tubule volume
fVPTL <- 0.4124 #Fractional volume of proximal tubule lumen to proximal tubule volume

fVRK <- 1-fVPT-fVGlom # should be 0.6018,Fractional volume of rest of kidney to total kidney volume
fVRKB <- 0.2019 #Fractional volume of rest of kidney blood to rest of kidney volume
fVRKC <- 0.5154 #Fractional volume of rest of kidney cell to rest of kidney volume
fVRKL <- 0.2826 #Fractional volume of rest of kidney lumen to rest of kidney volume



### Chemical Specific ####

MW <- PFOA_params$MW 

fup <- PFOA_params$fup

#### Partition coefficients ####
# Calculate Plasma/Rest of the body partition coefficient
# KpRe = partition coefficient of each of the lumped organs * fractional volume of the respective organ / sum of the fractional volume of all these organs
KpRe <- (PFOA_params$KpBr*Physio_params$V_brainFraction_M +
           PFOA_params$KpHe*Physio_params$V_heartFraction_M + 
           PFOA_params$KpLu*Physio_params$V_lungFraction_M +
           PFOA_params$KpMu*Physio_params$V_muscleFraction_M +
           PFOA_params$KpSp*Physio_params$V_spleenFraction_M +
           PFOA_params$KpGo*Physio_params$V_reproFraction_M +
           PFOA_params$KpBo*Physio_params$V_boneFraction_M) / (
             Physio_params$V_brainFraction_M +
               Physio_params$V_heartFraction_M +
               Physio_params$V_lungFraction_M +
               Physio_params$V_muscleFraction_M +
               Physio_params$V_spleenFraction_M +
               Physio_params$V_reproFraction_M +
               Physio_params$V_boneFraction_M)  
PFOA_params$KpRe <- KpRe

# Choose between calculated partition coefficients or initial ones (rat)
# In PFOA_params, PF, PG, PK, PL, PSK and PR are the initial PCs from Kudo 2007 (rat). These were recalculated to total tissue (/fup) to enable differentiating the fup for the sensitivity analysis
PAc <- PFOA_params$KpAd #PF # Adipose
PGc <- PFOA_params$KpGu #PG # Gut
PKc <- PFOA_params$KpKi #PK # Kidney
PLc <- PFOA_params$KpLi #PL # Liver
PSkc <- PFOA_params$KpSk #PSk # Skin
PRc <- PFOA_params$KpRe #PR #Rest


#### Biliary and Fecal Clearances ####
CLbiliaryc <- PFOA_params$CLbiliaryc #L/d/kg Fujii et al 2015 
CLfecalc <- PFOA_params$CLfaecesc # L/d/kg Fujii et al 2015 


## Final Parameters ####

### Physiological  ####

# Body volumes (L)
VA <- VAc * BW/0.9 + MassAc/0.9
VG <- VGc * BW
VK <- VKc * BW
VL <- VLc * BW
VP <- VPc * (1-Hct) * BW
VSk <- VSkc * BW
VR <- (VTotc*BW) - (VA + VG + VK + VL + VP + VSk)

VFil <- VK * 0.05 # Kidney filtrate compartment, corresponds to the volume of the collecting system in [ICRP 89 page 149] http://www.icrp.org/publication.asp?id=ICRP%20Publication%2089

# Body flows (L/d)
QA <- QAc * QC
QG <- QGc * QC
QK <- QKc * QC
QL <- QLc * QC
QSk <- QSkc * QC
QR <- QC - (QA + QG + QK + QL + QSk)

QUr = 0.022 * BW # L/d, Urine flow rate to the bladder 22 mL/kg BW/d [ICRP 89 page 161]
GFR = 113.7*60*24/1000 #0.18 * QK # L/d, Glomerular filtration rate 18% of total renal plasma flow [ICRP 89 page 159] http://www.icrp.org/publication.asp?id=ICRP%20Publication%2089
QT = 43.2*60*24/1000 # L/d; Tubular flow rate at the end of the PT 43.2 ml/min TFR from Scotcher et al. 2016 https://doi.org/10.1016/j.ejps.2016.03.018 and Pletz



#### Actual volumes of mechanistic kidney ####
VGlom <- fVGlom*VK  # Volume of glomeruli 
VGlomB <- fVGlomB*VGlom # Volume of glomerular blood 
VGlomL <- fVGlomL*VGlom # Volume of glomerular space(lumen)

VPT <- fVPT*VK  # Volume of proximal tubule 
VPTB <- fVPTB*VPT # Volume of proximal tubule 
VPTC <- fVPTC*VPT # Volume of proximal tubule 
VPTL <- fVPTL*VPT # Volume of proximal tubule 

VRK <- fVRK*VK # Volume of rest of kidney 
VRKB <- fVRKB*VRK # Volume of rest of kidney blood 
VRKC <- fVRKC*VRK # Volume of rest of kidney cell 
VRKL <- fVRKL*VRK # Volume of rest of kidney lumen 
RK_balance <- VRK - VRKB - VRKC - VRKL
RK_balance <- round(RK_balance, 4)

MB_kidney <- VK - sum(VGlomB, VGlomL, VPTB, VPTC, VPTL, VRKB, VRKC, VRKL)
MB_kidney <- round(MB_kidney, 4)

##### In the model currently ####
VGlomP <- VGlomB*Hct # Volume of glomerular blood, corrected to plasma
VGlomL <- VGlomL # Volume of glomerular space(lumen)
VPTP <- VPTB*Hct # Glom + prox + rest
VPTC <- VPTC #+ VRKC #prox + rest
VPTL <- VPTL #+ VRKL #glom + prox + rest
VRKP <- (VRKB*Hct)  
VRKC <- VRKC
VRKL <- VRKL

VPTP <- VPTP + VGlomP 
VPTC <- VPTC 

VPT <- VPTP + VPTC
VPTL <- VPTL + VGlomL 

VRKL <- VRKL
VRK <- VK - VPT - VPTL - VRKL


### Compound specific ####

#### Partition coefficients ####
# Here correcting for fraction unbound (as it was not incorporated in the input calculating file)
PA <- PAc * fup # Adipose
PG <- PGc * fup # Gut
PK <- PKc * fup # Kidney
PL <- PLc * fup # Liver
PSk <- PSkc * fup # Skin
PR <- PRc * fup #Rest

#### Biliary and Fecal Clearances ####
CLbiliary <- CLbiliaryc * BW #L/d  
CLfecal <- CLfecalc * BW # L/d 

kco = 2.016 # /d  0.0014 #/min
#kt = 50.4 # /d 0.035/min
#fecal excretion: what is in - what leaves (kt*agut) - what gets absorbed (ka*agut) + (EHR*Vmax.bile/(Kmbile +fuli*Cli))*fuli*Cli*Vbile


#### Renal clearance ####
# SF <- 0.17 * 10e6 # 17% of kidney is protein [ICRP 89], 10e6 is scaling from mg protein to kg protein, double ref for 17% protein Ruark 2020: DOI: https://doi.org/10.1016/B978-0-12-818596-4.00006-0
SF <- 10.9e-7 * 99.4e6 * 1e3 * (fVPT*VK) #6.54 mgprotein/HEK293cell (ref: Han and Ni, 2004; Ho et al., 2004) * PTCPGK cells/g kidney * 1e3 as Vkidney is in Kg (could be 99.4e6 or 60e6 see below comment ref: Neuhoff et al., 2013), equation from: https://doi.org/10.1016/j.comtox.2021.100172 
# Comment regarding PTCPGK: from Tang et al. 2024 https://doi.org/10.1021/acs.molpharmaceut.4c00504  a value of 60 million PTCPGKis commonly used but the observed value as high as 209 million PTCPGK has been reported. In this study, avalue of 99.4 million PTCPGK was applied based on the mostrecent meta-analysis.45 
# These scaling factors are not too far off. Scaling based on protein gives 1.7e6, scaling for cells gives 1.08346e6

#### Renal Clearance ####
# In vitro clearance
Vmax_OAT1 = 3.5*MW*1e-3*60*24*SF # ug/d/mg protein (nmol -> ug, min -> d) (scaled from nmol/min/mg protein); Louisse et al. 2024 doi.org/10.1016/j.tox.2024.153961
Vmax_OAT3 = 1.5*MW*1e-3*60*24*SF # ug/d/mg protein (nmol -> ug, min -> d) (scaled from nmol/min/mg protein); Louisse et al. 2024 doi.org/10.1016/j.tox.2024.153961
Vmax_OAT4 = 4.5*MW*1e-3*60*24*SF # ug/d/mg protein (nmol -> ug, min -> d) (scaled from nmol/min/mg protein); Louisse et al. 2024 doi.org/10.1016/j.tox.2024.153961

Km_OAT1 = 185*MW # ug/L (scaled from uM); Louisse et al. 2024 doi.org/10.1016/j.tox.2024.153961
Km_OAT3 = 90*MW # ug/L (scaled from uM); Louisse et al. 2024 doi.org/10.1016/j.tox.2024.153961
Km_OAT4 = 47*MW # ug/L (scaled from uM); Louisse et al. 2024 doi.org/10.1016/j.tox.2024.153961

CL_OAT1 = 19*1e-6*60*24 # L/d/mg protein; initial ul/min/mg protein, Louisse et al. 2024
CL_OAT3 = 17*1e-6*60*24 # L/d/mg protein; initial ul/min/mg protein, Louisse et al. 2024
CL_OAT4 = 96*1e-6*60*24 # L/d/mg protein; initial ul/min/mg protein,Louisse et al. 2024

# Relative expression factor
REF_OAT4 <- 1 #PFOA_params$REF_OAT4 # is equal to 1; as we don't have data on the in vitro expression of OAT4
REF_OAT1 <- 1 #PFOA_params$REF_OAT1 #
REF_OAT3 <- 1 #PFOA_params$REF_OAT3 #

# CL_FiltPT <- CL_OAT4 * REF_OAT4 * 0.17 * 0.7 * VK # L/d scaling for protein content instead of cell content; 17% of kidney is protein and 70% of kidney is cortex [ICRP 89], assuming that all the kidney protein is found in the cortex; this is an overestimation though; double ref for 17% protein Ruark 2020: DOI: https://doi.org/10.1016/B978-0-12-818596-4.00006-0
CL_FiltPT <- CL_OAT4 * REF_OAT4 * SF  # L/d, corrected to proximal tubule volume only
CL_PltPT <- ((CL_OAT1 * REF_OAT1) + (CL_OAT3 * REF_OAT3)) * SF  #L/d, corrected to proximal tubule volume only

##### Passive diffusion ####

pKa <- 1.886 # pKa of PFOA, average of experimental values from https://pfas-1.itrcweb.org, Table 4-1 excel file, pKa sheet

# pH at the different sections
pH_KC <- 7    # pH of intracellular tissue water R&R 2004 DOI 10.1002/jps.20322
pH_P <- 7.4    # pH of plasma
pH_PTL <- 7    # pH of proximal tubule urine Huang and Isoherranen 2018 doi:10.1002/psp4.12321
pH_RKL <- 6.75 # calculated average pH of the urine in the rest of the kidney, data Huang and Isoherranen 2018 doi:10.1002/psp4.12321

# fraction unionised
f.union_KC <- 1/(1 + 10^(pH_KC - pKa)) #in kidney cell (both proximal and rest of the kidney)
f.union_P <- 1/(1 + 10^(pH_P - pKa))
f.union_PTL <- 1/(1 + 10^(pH_PTL - pKa))
f.union_RKL <- 1/(1 + 10^(pH_RKL - pKa))

# surface area of the different sections
SA_PTL <- 611*30 # dm3 = L, correcting for the microvilli, as done by Huang and Isoherranen 2018 doi:10.1002/psp4.12321
SA_PT <- 611          # dm3 = L
SA_RK <- 280.5       # dm3 = L, sum of the surface area of all other compartments, SA_RK <- 61+61+125+6.7+6.7+6.7+6.7+6.7, from Huang and Isoherranen 2018 doi:10.1002/psp4.12321, Table 1

# input data
Papp <- 1.46 *1e-6 #cm/s in vitro permeability at apical compartment pH 7.4, PFAS were added to the donor wells and transport buffer containing 0.4% BSA was added to the receiver wells
f.union_exp <- f.union_P #as exp.pH=7.4 = plasma pH
Pint <- Papp/f.union_exp #cm/s intrinsic permeability, corrected for fraction unionised in the experiment

# Effective passive diffusion between either tubule and cell or cell and blood
# In Huang and Isoherranen 2018 doi:10.1002/psp4.12321 they assume the same Peff (what they call CL_PD) for apical and basolateral sides except for the proximal tubule where apical side has 30 fold higher TSA than basolateral side, due to the presence of microvilli
CLdif_PTLtPTC <- (Pint*SA_PTL*f.union_PTL/1000)*60*60*24 # L/d; cm/s = L/s /1000 = L/d *60*60*24; Proximal tubule lumen to proximal tubule cell
CLdif_PTCtPTP <- (Pint*SA_PT*f.union_KC/1000)*60*60*24   # L/d; cm/s = L/s /1000 = L/d *60*60*24; Proximal tubule cell to proximal tubule plasma
CLdif_PTCtPTL <- (Pint*SA_PT*f.union_KC/1000)*60*60*24   # L/d; cm/s = L/s /1000 = L/d *60*60*24; Proximal tubule cell to proximal tubule lumen
CLdif_PTPtPTC <- (Pint*SA_PT*f.union_P/1000)*60*60*24    # L/d; cm/s = L/s /1000 = L/d *60*60*24; Proximal tubule plasma to proximal tubule cell


##### Calculating fraction unbound in tissues ####

# To calculate fraction unbound in tissue based on the Poulin and Theil equation
# R is the ratio of average albumin and lipoprotein in tissue / plasma

# Albumin concentrations
# From Akihiro Tojo and Satoshi Kinugasa 2012 doi:10.1155/2012/481520
# Calb_GlomL <- 22.9 #ug/ml bowmans capsule
# Calb_PT <- 14.4 #ug/ml proximal tubule 
# Calb_DT <- 1.3 #ug/ml distal tubule
# Calb_Ur <- 0.7 #ug/ml urine
# In the same paper: the proximal tubule reabsorbes 71% of albumin, while LoH and DT 23% and the collecting duct 3%
Calb_P <- 37.0      #mg/ml plasma
Calb_PTL <- 14.4e-3 #mg/ml proximal tubule 
Calb_RKL <- 2.88e-3 #mg/ml rest of kidney tubule, as the PT is filtering 71%, then the concentration of albumin leaving the PT is 6.641ul/ml, in the urine it's 0.7ul/ml, therefore I'm doing the average here
Calb_exp <- 1e-10   #albumin, or protein concentration not reported in the experiments, therefore assuming a very low number

# Albumin ratio
R_T <- 0.5 #albumin and lipoprotein ratio between the tissue interstitial fluid and plasma
R_PTL <- Calb_PTL/Calb_P 
R_RKL <- Calb_RKL/Calb_P
R_exp <- Calb_exp/Calb_P

# Fraction unbound
fuT <- 1/(1 + ((1 - fup)/fup) * R_T)     # Poulin and Theil, 2009, below Table 6
fuPTL <- 1/(1 + ((1 - fup)/fup) * R_PTL)
fuRKL <- 1/(1 + ((1 - fup)/fup) * R_RKL)
fuexp <- 1/(1 + ((1 - fup)/fup) * R_exp) # is actually 1

# Enterohepatic circulation ####
VmaxBSEPc = 7.1 # umole/min/mg BSEP, average active transport of bile acids, assuming that the maximum velocity of PFOA transport by BSEP corresponds to that of bile acids 
KmBSEP = 16.4 # uM, average affinity constant of bile acids to BSEP, following the above assumption

aBSEP_all = 0.839 # pmoles/10^6 hepatocytes, BSEP protein abundance De Bruijn et al. 2024
MWBSEP_all = 140000 # g/mole, MW BSEP De Bruijn et al. 2024
Hep_all = 99 # 10^6 hepatocytes per g liver, hepatocellularity De Bruijn et al. 2024

SF_BSEP = aBSEP_all*MWBSEP_all*Hep_all*1e-9*1e3*VL; # mg BSEP/entire liver; scaling factor for BSEP mediated hepatic efflux for GCA and GCDC De Bruijn et al. 2024
VmaxBSEP = VmaxBSEPc*SF_BSEP*60*24*MW;   # ug/d   maximum   speed   for   BSEP-mediated GCDCA efflux (calculated) De Bruijn et al. 2024


VGL_up = 0.014/3*BW
VGL_low = 0.014/3*BW
VGL_col = 0.014/3*BW
ktj = 2.76*24 # /d
kti = 2.76*24 # /d
kco = 2.76*24 # /d
VmaxASBT = 
KmASBT =




# Model parameters ####
parms <- unlist(c(data.frame(VP,
                             VG, 
                             VL, 
                             PG,
                             PL, 
                             CLbiliary, 
                             fuT, 
                             kco, 
                             CLfecal,
                             VmaxBSEP,
                             KmBSEP, 
                             VGL_up, 
                             VGL_low, 
                             VGL_col, 
                             ka, 
                             ktj,
                             kti,
                             kco,
                             VmaxASBT,
                             KmASBT,
                             fuexp
                             )))

parms


# PBK MODEL ####
# ---------------------------------------------------------------------------- #

PBPKmodPFOA_M <- function(t, state, parameters){
  with(as.list(c(state, parameters)), {
    
    ## Dose ----
    # Oraldose <- if_else(t <= exposure_stop, Oralconc * BW, 0)
    # Dermaldose <- if_else(t <= exposure_stop, Dermconc * BW, 0) #+ AbsPFOA
    
    ## Concentrations ----
    
    # Organ concentrations (ug/L); these are TOTAL concentrations
    CP <- AP/VP  # Concentration in plasma (ug/L)
    
    CG <- AG/VG  # Concentration in gut (ug/L)
    CGL_up <- AGL_up/VGL_up 
    CGL_low <- AGL_low/VGL_low
    CGL_col <- AGL_col/VGL_col
    
    CL <- AL/VL  # Concentration in liver (ug/L)
    
    CVG <- CG/PG # Concentration leaving gut (ug/L)
    CVL <- CL/PL  # Concentration leaving liver (ug/L)
   
    
    ## Differential equations ----
    dD = - D
   
    # Gut compartment: 
    dAGL_up <- D - CLdif_PTLtPTC*CGL_up - ktj*AGL_up + (VmaxBSEP/(KmBSEP + (CL*fuT)))*CL*fuT 
    dAGL_low <- ktj*AGL_up - CLdif_PTLtPTC*CGL_low - kti*AGL_low - (VmaxASBT/(KmASBT + (CGL_low*fuexp)))*CGL_low*fuexp
    dAGL_col <- kti*AGL_low - CLdif_PTLtPTC*CGL_low - kco*AG_col
    
    dAG <- QG*(CP-CVG) +  CLdif_PTLtPTC*CGL_up + CLdif_PTLtPTC*CGL_low + CLdif_PTLtPTC*CGL_low + (VmaxASBT/(KmASBT + (CGL_low*fuexp)))*CGL_low*fuexp #+ Oraldose - kco*AG + CLbiliary*CL*fuT 
    
    # Excretion fecal: cumulative
    dAEx_feces <- kco*AG_col #CLfecal*CG #kco*AG
    
    # Liver compartment
    dAL <- QL*CP + QG*CVG - (QL+QG)*CVL - (VmaxBSEP/(KmBSEP + (CL*fuT)))*CL*fuT #CLbiliary*CL*fuT 
    
    
    # Plasma compartment
    dAP <- - (QG + QL)*CP +
      (QL+QG)*CVL
    
    
    # Mass Balance
    Atot <- D + AP + AG + AL + AEx_feces + AGL_up + AGL_low + AGL_col
    
    # dInput <- Oraldose + Dermaldose
    # MB = dInput - Atot
    MB = Oraldose - Atot + 1
    
    # End
    
    list(c(dD,
           dAP, 
           dAGL_up, 
           dAGL_low, 
           dAGL_col,
           dAG, 
           dAL, 
           dAEx_feces
    ), 
    c(CP = CP, 
      CGL_up = CGL_up, 
      CGL_low = CGL_low, 
      CGL_col = CGL_col,
      CG = CG, CVG = CVG, 
      CL = CL, CVL = CVL,
      Atot = Atot,MB = MB)
    )
  })
}

## Initials ####

A_init <- c(D = Oraldose,
            AP = 0,
            AG = 0, 
            AGL_up = 0, 
            AGL_low = 0, 
            AGL_col = 0, 
            AL = 0, 
            AEx_feces = 0
)


## Solving the model ####
output_PFOA <- lsoda(y = A_init, 
                     times = TIME, 
                     func = PBPKmodPFOA_M, 
                     parms = parms)
output.PFOA.df <- as.data.frame(output_PFOA) 


# RESULTS ####
# ---------------------------------------------------------------------------- #

output.PFOA.df <- output.PFOA.df %>% 
  # mutate(time = time/365) %>% 
  rename(Days = time)

# MASSBALANCE
MB.df <- output.PFOA.df %>% select(Days, Atot, MB)
MB.df$MB <- round(MB.df$MB, 3)
MB.df$ERROR <- (Oraldose - MB.df$Atot) / MB.df$Atot * 100
MB.df$ERROR <- round(MB.df$ERROR, 3)
MB_plot <- ggplot(data = MB.df)+
  geom_line(aes(x = Days, y = ERROR, color = "ERROR")) +
  geom_line(aes(x = Days, y = MB, color = "MB")) +
  scale_color_manual(values = c("ERROR" = "blue", "MB" = "black")) +
  # ylim(0,1) +
  theme_minimal() +
  ylab("MB / ERROR")
MB_plot

# PFOA in plasma of one individual
Plot_PFOA_Plasma <- ggplot()+
  geom_path(data = output.PFOA.df, aes(x = Days, y = CP))+
  theme_CP()+
  # theme(axis.text.x = element_text(size = 7),axis.text.y = element_text(size = 7), axis.title = element_text(size = 8))+
  ylab("Plasma (ng/ml)")
Plot_PFOA_Plasma

ggsave("PlasmaConcentration.png", dpi = 300)


# PFOA in Liver and Gut
Plot_PFOA_LiverGut <- output.PFOA.df %>%
  select(Days, CP, CL, CG) %>% 
  filter(Days > 0.05) %>%
  rename(Total_Plasma = CP,
         Liver = CL,
         Gut = CG
  ) %>%
  pivot_longer(names_to = "Organ", values_to = "Concentration", Total_Plasma:Gut) %>%
  ggplot()+
  geom_path(aes(x = Days, y = Concentration, color = Organ), linewidth = 0.5) +
  facet_wrap(~ Organ)+
  theme_CP() +
  theme(legend.position = "null")+
  #       panel.background = element_rect(fill = "white"),
  #       panel.border = element_blank()
  #       ) +
  ylab("Concentration (ng/ml)") +
  xlab("Time (d)")
Plot_PFOA_LiverGut
ggsave("Plot_PFOA_LiverGut.png", dpi = 300)


# # Calculate AUC and Half life ####
# # ---------------------------------------------------------------------------- #
# 
# AUC <- trapz(output_PFOA[ , "time"], output_PFOA[ , "CP"]) #ug*day/L
# 
# # Calculate predicted half-life
# time <- output_PFOA[ , "time"] #in days
# conc <- output_PFOA[ , "CP"] #(ug/L) or (ng/ml)
# Cmax <- max(conc)
# Tmax <- time[which.max(conc)]
# tlast <- max(time[conc > 0])
# 
# half_life <- pk.calc.half.life(
#   conc,
#   time,
#   Tmax,
#   tlast
# )
# 
# HalfLife <- half_life$half.life/365 #halflife in years
# 
# # Experimental
# ExpData <- read_csv("C:/Users/pacho003/OneDrive - Wageningen University & Research/CP_L_R/PARC_PFAS_PBPKmodel/Codes_PFAS_models/Chrysanthi_Pachoulide_2024_PFOA_mechanistic/Input/HalfLifes.csv")
# 
# Exp_HalfLifes <- ExpData %>%
#   filter(species == "human",
#          chemical == "pfoa",
#          parameter== "HalfLife") %>%
#   select(value_average) %>%
#   rename(HalfLife = value_average)
# Exp_HalfLifes$HalfLife <- as.numeric(Exp_HalfLifes$HalfLife) #years
# 
# Predicted.df <- data.frame(
#   HalfLife = HalfLife,
#   Origin = "Predicted",
#   value = 1)
# Experimental.df <- data.frame(
#   HalfLife = Exp_HalfLifes$HalfLife,
#   Origin = "Experimental",
#   value = 1)
# 
# HalfLifes <- rbind(Predicted.df, Experimental.df)
# 
# Plot_HalfLifes <- HalfLifes %>%
#   ggplot()+
#   geom_violin(data = Experimental.df, aes(value, HalfLife),
#               color = "transparent",
#               fill = "grey")+
#   geom_point(data = Predicted.df, aes(value, HalfLife),
#              color = "slateblue3", size = 5, shape = 18) +
#   ylab("Half life (years)") +
#   theme_CP() +
#   theme(axis.text.x=element_blank(),
#         axis.ticks.x=element_blank(),
#         axis.title.x = element_blank()
#   )
# Plot_HalfLifes
# ggsave("ExpVsSimHalfLife.png", dpi = 300)
# 
# 
# Experimental Vs Simulated Plasma concentrations ####
ExpPlasma <- read_excel("C:/Users/pacho003/OneDrive - Wageningen University & Research/CP_L_R/PARC_PFAS_PBPKmodel/Codes_PFAS_models/Chrysanthi_Pachoulide_2024_PFOA_mechanistic/Input/Experimental.Plasma.PFOA.xlsx", 
                        col_types = c("numeric", "numeric"))


ExpPlasma <- ExpPlasma %>% 
  # mutate(Years = Time_days/365) %>% 
  # select(!Time_days) %>% 
  rename(Days = Time_days) %>% 
  rename(CP = MPFOA_µg_per_L) %>% # ug/L or ng/ml
  mutate(CP = CP - 0.130) %>%  #substracting the pre-existing level of 0.130ug/L from their previous study, as also done in the ref. article: https://doi.org/10.1016/j.envint.2024.109047 (table 3)
  filter(Days <=TSTOP)


Plot_Plasma <- ggplot()+
  geom_path(data = output.PFOA.df, aes(x = Days, y = CP), color = "aquamarine", linewidth = 1.5)+
  geom_point(data = ExpPlasma, aes(x = Days, y = CP), color = "black")+
  theme_CP()+
  # scale_color_manual(values = c("Simulated" = "aquamarine",
  #                               "Experimental" = "black")) +
  # # theme(axis.text.x = element_text(size = 7),axis.text.y = element_text(size = 7), axis.title = element_text(size = 8))+
  ylab("Plasma (ng/ml)")
Plot_Plasma
ggsave("PlasmaExpVsPredicted.png", dpi = 300)

# 
# print(AUC)
# print(HalfLife)