  # --------------------------------------------------------------------------- #
  # PBK MODEL FOR PFOA, TO BE USED TOGETHER WITH THE LATEST HBM DATA
  # Model File
  # CP, 27-01-2025
  # --------------------------------------------------------------------------- #
  
  rm(list=ls()) # to clear out the global environment
  
  # Set working directory
  HOME <- "C:/Users/pacho003/OneDrive - Wageningen University & Research/CP_L_R/PARC_PFOA_mechanistic"
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
  library(sensitivity)
  
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
  
  Physio_params <- read_csv("C:/Users/pacho003/OneDrive - Wageningen University & Research/CP_L_R/PARC_PFOA_mechanistic/Input/PhysioVariables.csv")
  PFOA_params <- read_csv("C:/Users/pacho003/OneDrive - Wageningen University & Research/CP_L_R/PARC_PFOA_mechanistic/Input/PFOAParams.csv")
  
  
  # EXPOSURE SCENARIO ####
  # ------------------------------------------------------ #
  
  # exposure_stop <-  50*365      # days
  sim_stop <- 450                 # days Abraham study follow-up time  
  
  TSTART <- 0
  TSTOP <- sim_stop               # days
  DT <- 1/10                    # days
  TIME <- seq(TSTART,TSTOP,by=DT)
  
  
  Oraldose <- 3.96 # ug https://doi.org/10.1016/j.envint.2024.109047
  # Oralconc <- 0.000187 # ug/kg/day [EFSA 2020, page 143] 3.96ug
  # Dermconc <- 0.000542 # ug/kg/day, mean of #as.numeric(SumExpPFOA_LB_val[i,14])
  
  
  # PBK MODEL PARAMETERS ####
  # ---------------------------------------------------------------------------- #
  
  ## Physiological  ####
  Physio_params <- Physio_params %>% 
    filter(age == 60) %>%  # 60 years old adult, male, as in the Abraham study
    select(ends_with('_M')) %>% 
    mutate(BloodFlowSum = rowSums(select(., starts_with("Q_")))) %>% # 0.9935, total blood flow as the sum of the fractional blood flows of all organs on which we have data
    mutate(VolumesSum = rowSums(select(., starts_with("V_")))) # 0.96, total volume as the sum of the fractional organ volumes of all organs on which we have data
  
  BW <- 82                         # L (kg), The body weight of the subject in Abraham et al. 2024 https://doi.org/10.1016/j.envint.2024.109047 Physio_params$BDW_M 
  QC <- Physio_params$CardOut_M    # L/d, This is corrected for hematocrit already so it's plasma
  Hct <- 46.7/100                  # The hematocrit of the subject in Abraham et al. 2024 https://doi.org/10.1016/j.envint.2024.109047Physio_params$Hct_M
  
  
  ### Organ volumes -------------------------
  
  # Kidney
  VK <- Physio_params$V_kidneyFraction_M * BW             # L, Volume of kidney
  
  # Kidney is divided in to proximal tubule and rest of kidney. Each compartment is divided in to tissue and lumen (or filtrate). 
  # This separation is necessary because of the active reabsorption happening in the proximal tubule, but also due to the differences in composition and pH of primary and terminal urine.
  # Fractional volumes were recalculated from Pletz et al. 2020 https://doi.org/10.1016/j.comtox.2021.100172
  
  VGlomc <- 0.0400*VK                # L, Total volume of glomeruli 
  VGlomPc <- 0.6815*VGlomc*Hct       # L, Volume of glomerular plasma 
  VGlomLc <- 0.3185*VGlomc           # L, Volume of glomerular space(lumen)
  
  VPTc <- 0.3581*VK                  # L, Total volume of proximal tubule 
  VPTPc <- 0.0766*VPTc*Hct           # L, Volume of proximal tubule plasma 
  VPTCc <- 0.5110*VPTc               # L, Volume of proximal tubule cell
  VPTLc <- 0.4124*VPTc               # L, Volume of proximal tubule lumen
  
  VPT <- VGlomPc + VPTPc + VPTCc     # L, Volume of proximal tubule tissue, used in the model (proximal tubule tissue = volume of proximal tubule and glomerular blood + volume of cells in the proximal tubule)
  VPTL <- VGlomLc + VPTLc            # L, Volume of proximal tubule lumen, used in the model (proximal tubule lumen = volume of proximal tubule lumen and glomerular space)
  
  VRKc <- VK - VGlomc - VPTc         # L, Total volume of rest of kidney 
  # VRKP <- 0.2019*VRK                # L, Volume of rest of kidney blood 
  # VRKCc <- 0.5154*VRK               # L, Volume of rest of kidney cell 
  VRKL <- 0.2826*VRKc                # L, Volume of rest of kidney lumen, used in the model 
  VRK <- VK - VPT - VPTL - VRKL      # L, Volume of rest of kidney tissue, used in the model
  
  # Surface area of the different sections 
  SA_PTL <- 611*30                   # Proximal tubule, lumen side, correcting for the microvilli, as done by Huang and Isoherranen 2018 doi:10.1002/psp4.12321
  SA_PT <- 611                       # Proximal tubule, cell side
  SA_RK <- 280.5                     # Rest of kidney, sum of the surface area of all other compartments, SA_RK <- 61+61+125+6.7+6.7+6.7+6.7+6.7, from Huang and Isoherranen 2018 doi:10.1002/psp4.12321, Table 1
  
  
  
  # Gut
  VG <- Physio_params$V_gutFraction_M * BW
  
  # Gut lumen is taken as a compartment outside the intestine (weight of gut lumen is outside of the bodyweight)
  L = 280                           # cm, (adult of 70kg) Willmann2004 doi: 10.1021/jm030999b
  R = (1.75+1)/2                    # cm, (adult of 70kg mean value) Willmann2004 doi: 10.1021/jm030999b
  VGL = pi*L*(R^2)/1000             # L, cm3/1000 #Value is the same as Punt et al. 2021 https://dx.doi.org/10.1021/acs.chemrestox.0c00307
  
  # L_up = 124                        # cm, (adult of 70kg) (duodenum 20, jejunum 104 cm) Willmann2004 doi: 10.1021/jm030999b
  # R_up = 1.75                       # cm
  # VGL_up = pi*L_up*(R_up^2)/1000    # L, cm3/1000
  # 
  # L_low = 156                       # cm, (adult of 70kg) ileum
  # R_low = 1                         # cm, Willmann2004 doi: 10.1021/jm030999b
  # VGL_low = pi*L_low*(R_low^2)/1000 # L, cm3/1000
  # 
  # L_col = 1.5/100                   # cm, (1.5m/100) Willmann2004 doi: 10.1021/jm030999b
  # R_col = 3.5                       # cm, Willmann2004 doi: 10.1021/jm030999b
  # VGL_col = pi*L_col*(R_col^2)/1000 # L, cm3/1000

  SA_SI = 2*pi*R*L*25               # cm2, amplification factor of 25 for the microvilli in the intestinal lumen, Willmann2004 doi: 10.1021/jm030999b
  
  
  # Liver
  VL <- Physio_params$V_liverFraction_M * BW             # L, Volume of liver
  
  # # Liver is divided in intracellular and extracellular compartments. 
  # # Extracellular compartment combines both the vascular and interstitial space
  # # Fractional volume of intracellular space was taken from Utsey et al. 2020 https://doi.org/10.1124/dmd.120.090498, https://github.com/metrumresearchgroup/PBPK_PC/blob/master/data/unified_tissue_comp.csv
  VL_icc = 0.573                    # Fractional volume of intracellular space in the liver
  VL_ecc = 1 - 0.573                # Fractional volume of extracellular space in the liver
  
  VL_ic = VL_icc*VL                 # L, Volume of liver intracellular compartment, used in the model
  VL_ec = VL_ecc*VL                 # L, Volume of liver extracellular compartment, used in the model
  
  
  # Skin
  VSk <- Physio_params$V_skinFraction_M * BW              # L, Volume of skin
  
  
  # Other compartments
  VA <- Physio_params$V_adiposeFraction_M * BW/0.9 + Physio_params$AdiposeMass_M/0.9
  VP <- Physio_params$V_plasmaFraction_M * (1-Hct) * BW   # L, Volume of plasma
  VTot <- Physio_params$VolumesSum * BW                   # L, Total body volume (used for mass balance)
  VR <- VTot - (VA + VG + VK + VL + VP + VSk)             # L, Volume of the lumped rest compartment
  
  
  ### Organ blood flows -------------------------
  
  QK <- Physio_params$Q_kidneyFraction_M/Physio_params$BloodFlowSum * QC  # L/d, Kidney
  QG <- Physio_params$Q_gutFraction_M/Physio_params$BloodFlowSum * QC     # L/d, Gut
  QL <- Physio_params$Q_liverFraction_M/Physio_params$BloodFlowSum * QC   # L/d, Liver
  QSk <- Physio_params$Q_skinFraction_M/Physio_params$BloodFlowSum * QC   # L/d, Skin 
  QA <- Physio_params$Q_adiposeFraction_M/Physio_params$BloodFlowSum * QC # L/d, Rest
  QR <- QC - (QA + QG + QK + QL + QSk)
  
  
  ### Other physiological flows and constants -------------------------
  
  QUr = 0.022 * BW            # L/d, Urine flow rate to the bladder 22 mL/kg BW/d [ICRP 89 page 161]
  GFR = 113.7*60*24/1000      # The measured GFR of the subject in Abraham et al. 2024 https://doi.org/10.1016/j.envint.2024.109047Physio_params$Hct_M  or 0.18 * QK # L/d 18% of total renal plasma flow [ICRP 89 page 159] http://www.icrp.org/publication.asp?id=ICRP%20Publication%2089
  QT = 43.2*60*24/1000        # L/d, Tubular flow rate at the end of the proximal tubule, 43.2 ml/min TFR from Scotcher et al. 2016 https://doi.org/10.1016/j.ejps.2016.03.018 and Pletz
  
  # tge = 2*24                  #/d (24*/h) or 30min (10-60min), gastric emptying rate, fasted Willmann2004 doi: 10.1021/jm030999b, Punt et al 2021 https://dx.doi.org/10.1021/acs.chemrestox.0c00307
  tsi = 0.3*24                #/d (24*/h) or 4h (2-6 h), intestinal transit time, Willmann2004 doi: 10.1021/jm030999b (#transit time duodenum 14min, jejunum 71min, ileum 114min B. Agoram et al. 2001), Punt et al 2021 https://dx.doi.org/10.1021/acs.chemrestox.0c00307
  # tsi_up = 0.7*24             #/d (24*/h) 1.417 h, transit time from duodenum and jejunoum ,14+17 min (duodenum 14min, jejunum 71min, from B. Agoram et al. 2001), Punt et al 2021 https://dx.doi.org/10.1021/acs.chemrestox.0c00307
  # tsi_low = 1.9*24            #/d (24*/h) 1.9 h, transit time of ileum ,114 min from B. Agoram et al. 2001, Punt et al 2021 https://dx.doi.org/10.1021/acs.chemrestox.0c00307
  tco = 0.09*24               #/d (24*/h) or 7h willmann2004 doi: 10.1021/jm030999b
  
  
  ### Physiological pHs in different matrices -------------------------
  
  pH_KC <- 7      # intracellular tissue water, R&R 2004 DOI 10.1002/jps.20322
  pH_P <- 7.4     # plasma
  pH_PTL <- 7     # proximal tubule urine, Huang and Isoherranen 2018 doi:10.1002/psp4.12321
  pH_RKL <- 6.75  # calculated average pH of the urine in the rest of the kidney, data from Huang and Isoherranen 2018 doi:10.1002/psp4.12321
  pH_L <- 7       # intestinal gut lumen, average
  pH_up <- 6.5    # fasted (6.0 duodenum, 7.0 jejunum) willmann2004 doi: 10.1021/jm030999b
  pH_low <- 7.5   # fasted willmann2004 doi: 10.1021/jm030999b
  
  
  ### Albumin concentrations in different matrices -------------------------
  
  # From Akihiro Tojo and Satoshi Kinugasa 2012 doi:10.1155/2012/481520
  # In the same paper: the proximal tubule reabsorbes 71% of albumin, while LoH and DT 23% and the collecting duct 3%
  Calb_P <- 37.0      # mg/ml plasma
  Calb_PTL <- 14.4e-3 # mg/ml proximal tubule 
  Calb_RKL <- 2.88e-3 # mg/ml rest of kidney tubule, as the PT is filtering 71%, then the concentration of albumin leaving the PT is 6.641ul/ml, in the urine it's 0.7ul/ml, therefore I'm doing the average here
  Calb_exp <- 1e-10   # albumin, or protein concentration not reported in the experiments, therefore assuming a very low number
  
  
  
  ## Chemical Specific ####
  
  MW <- PFOA_params$MW 
  fup <- PFOA_params$fup        # fraction unbound in plasma from Fisher et al.2024
  pKa <- 1.886                  # pKa of PFOA, average of experimental values from https://pfas-1.itrcweb.org, Table 4-1 excel file, pKa sheet
  
  
  ### Partition coefficients -------------------------
  
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
  # Correcting for fraction unbound (as it was not incorporated in the input calculating file)
  PK <- PFOA_params$KpKi * fup  #PK   # Kidney
  PG <- PFOA_params$KpGu * fup  #PG   # Gut
  PL <- PFOA_params$KpLi * fup  #PL   # Liver
  PSk <- PFOA_params$KpSk * fup #PSk  # Skin
  PA <- PFOA_params$KpAd * fup  #PF   # Adipose
  PR <- PFOA_params$KpRe * fup  #PR   # Rest
  
  
  ### Fraction unionised in different matrices -------------------------
  
  f.union_KC <- 1/(1 + 10^(pH_KC - pKa))     # Kidney cell (both proximal and rest of the kidney)
  f.union_P <- 1/(1 + 10^(pH_P - pKa))       # Plasma
  f.union_PTL <- 1/(1 + 10^(pH_PTL - pKa))   # Proximal tubule lumen
  f.union_RKL <- 1/(1 + 10^(pH_RKL - pKa))   # Rest of kidney lumen
  f.union_exp <- f.union_P                   # As exp.pH=7.4 = plasma pH
  f.union_GT <- 1/(1 + 10^(pH_L - pKa))      # Gut lumen
  
  
  ### Fraction unbound in different matrices -------------------------
  
  # Based on the Poulin and Theil 2009, below Table 6

  # Albumin ratio
  R_T <- 0.5 # albumin and lipoprotein ratio between the tissue interstitial fluid and plasma
  R_PTL <- Calb_PTL/Calb_P 
  R_RKL <- Calb_RKL/Calb_P
  R_exp <- Calb_exp/Calb_P
  R_L_ec <- 0.086 # albumin ratio liver, Utsey et al. 2020 https://doi.org/10.1124/dmd.120.090498, https://github.com/metrumresearchgroup/PBPK_PC/blob/master/data/unified_tissue_comp.csv
  
  # Fraction unbound
  fuT <- 1/(1 + ((1 - fup)/fup) * R_T)       # Tissue
  fuPTL <- 1/(1 + ((1 - fup)/fup) * R_PTL)   # Proximal tubule lumen
  fuRKL <- 1/(1 + ((1 - fup)/fup) * R_RKL)   # Rest of kidney lumen
  fuexp <- 1/(1 + ((1 - fup)/fup) * R_exp)   # Experiment (is actually 1)
  fuL_ec <- 1/(1 + ((1 - fup)/fup) * R_L_ec) # Liver extracellular space
  # Note Chrysa: need to check if I find the amount of albumin in liver interstitial space, as here the fuL_ec is 10 times higher than the fup, but extracellular space is actually mainly the albumin from the vascular space. 
  
  
  ### Renal Clearance -------------------------
  
  ## Active transport
  
  # Input data, in vitro clearance
  # Vmax_OAT1c = 3.5                 # nmol/min/mg protein, Louisse et al. 2024 doi.org/10.1016/j.tox.2024.153961
  # Vmax_OAT3c = 1.5                 # nmol/min/mg protein, Louisse et al. 2024 doi.org/10.1016/j.tox.2024.153961
  Vmax_OAT4c = 4.5                 # nmol/min/mg protein, Louisse et al. 2024 doi.org/10.1016/j.tox.2024.153961
  
  # Km_OAT1 = 185*MW                 # ug/L, scaled from uM, Louisse et al. 2024 doi.org/10.1016/j.tox.2024.153961
  # Km_OAT3 = 90*MW                  # ug/L, scaled from uM, Louisse et al. 2024 doi.org/10.1016/j.tox.2024.153961
  Km_OAT4 = 47*MW                  # ug/L, scaled from uM, Louisse et al. 2024 doi.org/10.1016/j.tox.2024.153961
  
  # CL_OAT1 = 19*1e-6*60*24          # L/d/mg protein, initial ul/min/mg protein, Louisse et al. 2024
  # CL_OAT3 = 17*1e-6*60*24          # L/d/mg protein, initial ul/min/mg protein, Louisse et al. 2024
  CL_OAT4 = 96*1e-6*60*24          # L/d/mg protein, initial ul/min/mg protein,Louisse et al. 2024
  
  # Relative expression factor
  # REF_OAT1 <- 1 #PFOA_params$REF_OAT1 #
  # REF_OAT3 <- 1 #PFOA_params$REF_OAT3 #
  REF_OAT4 <- 1 #PFOA_params$REF_OAT4 # is equal to 1, as we don't have data on the in vitro expression of OAT4
  
  
  # SF_OAT <- 0.17 * 10e6 # 17% of kidney is protein [ICRP 89], 10e6 is scaling from mg protein to kg protein, double ref for 17% protein Ruark 2020: DOI: https://doi.org/10.1016/B978-0-12-818596-4.00006-0
  SF_OAT <- 10.9e-7 * 99.4e6 * 1e3 * VPT #6.54 mgprotein/HEK293cell (ref: Han and Ni, 2004, Ho et al., 2004) * PTCPGK cells/g kidney * 1e3 as Vkidney is in Kg (could be 99.4e6 or 60e6 see below comment ref: Neuhoff et al., 2013), equation from: https://doi.org/10.1016/j.comtox.2021.100172 
  # Comment regarding PTCPGK: from Tang et al. 2024 https://doi.org/10.1021/acs.molpharmaceut.4c00504  a value of 60 million PTCPGKis commonly used but the observed value as high as 209 million PTCPGK has been reported. In this study, avalue of 99.4 million PTCPGK was applied based on the mostrecent meta-analysis.45 
  
  
  # Scaled clearances for active transport
  # CL_FiltPT <- CL_OAT4 * REF_OAT4 * SF_OAT                            # L/d, reabsorption
  # CL_PltPT <- ((CL_OAT1 * REF_OAT1) + (CL_OAT3 * REF_OAT3)) * SF_OAT  # L/d, excretion
  
  # Vmax_OAT1 = Vmax_OAT1c*MW*1e-3*60*24*SF_OAT # ug/d (nmol -> ug, min -> d)
  # Vmax_OAT3 = Vmax_OAT3c*MW*1e-3*60*24*SF_OAT # ug/d (nmol -> ug, min -> d)
  Vmax_OAT4 = Vmax_OAT4c*MW*1e-3*60*24*SF_OAT # ug/d (nmol -> ug, min -> d)
  
  
  ## Passive permeability
  
  # Input data, in vitro apparent permeability
  # Papp_PT <- 1.46*1e-6                       # cm/s, In vitro permeability at apical compartment pH 7.4, PFAS were added to the donor wells and transport buffer containing 0.4% BSA was added to the receiver wells
  # Pint_PT <- Papp_PT/f.union_exp             # cm/s, Intrinsic permeability, corrected for fraction unionised in the experiment
  
  # Final clearance for passive permeability (what is called effective passive diffusion)
  # Equations from  Huang and Isoherranen 2018, what it's called effective passive diffusion doi:10.1002/psp4.12321
  # In Huang and Isoherranen, it is assumed the same Peff (what they call CL_PD) for apical and basolateral sides except for the proximal tubule where apical side has 30 fold higher TSA than basolateral side, due to the presence of microvilli
  # CLdif_PTLtPTC <- (Pint_PT*SA_PTL*f.union_PTL/1000)*60*60*24 # L/d, Proximal tubule lumen to proximal tubule cell (calculations: cm/s = L/s /1000 = L/d *60*60*24)
  # CLdif_PTCtPTP <- (Pint_PT*SA_PT*f.union_KC/1000)*60*60*24   # L/d, Proximal tubule cell to proximal tubule plasma (calculations: cm/s = L/s /1000 = L/d *60*60*24)
  # CLdif_PTCtPTL <- (Pint_PT*SA_PT*f.union_KC/1000)*60*60*24   # L/d, Proximal tubule cell to proximal tubule lumen (calculations: cm/s = L/s /1000 = L/d *60*60*24)
  # CLdif_PTPtPTC <- (Pint_PT*SA_PT*f.union_P/1000)*60*60*24    # L/d, Proximal tubule plasma to proximal tubule cell (calculations: cm/s = L/s /1000 = L/d *60*60*24), 
  
  
  ### Uptake from gut -------------------------
  
  # Input data, in vitro clearance
  Papp_SI = 7.31*1e-6                         # cm/s, 7.31 ± 0.43, Janssen et al. 2024
  Pint_SI <- Papp_SI/f.union_exp # cm/s, Intrinsic permeability, corrected for fraction unionised in the experiment
  
  # The ka way
  # log_Peff = 0.4926*log10(Papp_SI) - 0.1454   # Equation for passively and actively absorbed compounds, Hou, Zhang et al. 2004
  # ka = 10^(log_Peff)*2/R * 3600 * 24          # /d, (24*/h), Yu and Amidon 1999
  
  # The Clearance way
  CL_GLtG <- (Pint_SI*SA_SI*f.union_GT/1000)*60*60*24 # L/d, Gut lumen to gut cell (calculations: cm/s = L/s /1000 = L/d *60*60*24) 
  
  
  ### Uptake to the liver -------------------------
  
  # Input data, in vitro clearance
  Vmax_OATP1B1c = 2.305 * 1e-6                # umol/min/mg protein, 2.305± 0.295 pmol/min/mg protein [@lin2023]
  Km_OATP1B1 = 52.65 * MW                     # ug/L, 52.65 ± 23.28 uM [@lin2023]
  
  Vmax_OATP1B3c = 2.694 * 1e-6                # umol/min/mg protein, 2.694± 0.470 pmol/min/mg protein [@lin2023]
  Km_OATP1B3 = 91.61 * MW                     # ug/L 91.61 ± 47.70 uM [@lin2023]
  

  # Relative expression factors
  OATP1B1_vitro = 0.120                       # [@lin2023, tables6, ref23]
  OATP1B1_vivo = 2.000                        # pmol/mg membrane protein [@lin2023, tables7, ref23]
  REF_OATP1B1 = OATP1B1_vivo/OATP1B1_vitro
  
  OATP1B3_vitro = 0.719                       # [@lin2023, tables6, ref23]
  OATP1B3_vivo = 1.000                        # pmol/mg membrane protein [@lin2023, tables7, ref23]
  REF_OATP1B3 = OATP1B3_vivo/OATP1B3_vitro
  
  # Scaled clearances
  Vmax_OATP1B1 = Vmax_OATP1B1c*REF_OATP1B1
  Vmax_OATP1B3 = Vmax_OATP1B3c*REF_OATP1B3
  
  
  ### Biliary clearance -------------------------
  
  # Initial value, from in vivo study
  # CLbiliary <- PFOA_params$CLbiliaryc * BW       # L/d, Biliary clearance, Fujii et al 2015 
  
  # Input data, in vitro clearance
  VmaxBSEPc <- 7.1                               # umol/min/mg BSEP, Average active transport of bile acids, assuming that the maximum velocity of PFOA transport by BSEP corresponds to that of bile acids 
  KmBSEP <- 16.4 * MW                            # ug/L, uM, Average affinity constant of bile acids to BSEP, following the above assumption
  
  # Scaling factor
  SF_BSEP <- 0.839*140000 * 99 * 1e-9 * 1e3 * VL # Scaling factor for BSEP mediated hepatic efflux for GCA and GCDC De Bruijn et al. 2024 (calculation: SF_BSEP = aBSEP_all*MWBSEP_all * Hep_all * 1e-9 * 1e3 * VL )
  
  # Scaled clearance
  VmaxBSEP <- VmaxBSEPc*SF_BSEP*60*24*MW         # ug/d
  
  
  ### Fecal clearance -------------------------
  
  # CLfecal <- PFOA_params$CLfaecesc * BW          # L/d, Fecal clearance, Fujii et al 2015 
  

  
  
  parms <- unlist(c(data.frame(BW, 
                               VPT,
                               VPTL,
                               VRK,
                               VRKL,
                               VGL,
                               VL_ec,
                               VL_ic,
                               VSk,
                               VG, 
                               VL, 
                               VA, 
                               VP, 
                               VR, 
                               QK, 
                               QSk,
                               QG,  
                               QL, 
                               QA, 
                               QR, 
                               QUr,
                               GFR,
                               QT,
                               tco,
                               PA, 
                               PG, 
                               PK, 
                               PL, 
                               PSk,
                               PR, 
                               fup,
                               fuT,
                               fuPTL,
                               fuL_ec,
                               Vmax_OAT4, 
                               Km_OAT4,
                               VmaxBSEP,
                               KmBSEP, 
                               CL_GLtG,
                               Vmax_OATP1B1,
                               Km_OATP1B1,
                               Vmax_OATP1B3,
                               Km_OATP1B3
  )))
  
  parms
  
  
  # PBK MODEL ####
  # ---------------------------------------------------------------------------- #
  
  PBPKmodPFOA_M <- function(t, state, parameters){
    with(as.list(c(state, parameters)), {
      
      ## Dose -------------------------
      
      # Oraldose <- if_else(t <= exposure_stop, Oralconc * BW, 0)
      # Dermaldose <- if_else(t <= exposure_stop, Dermconc * BW, 0) #+ AbsPFOA
      
      
      ## Concentrations -------------------------
      
      CPT <- APT/VPT               # ug/L, Kidney, proximal tubule
      CPTL <- APTL/VPTL            # ug/L, Kidney, proximal tubule lumen
      CRK <- ARK/VRK               # ug/L, Rest of kidney
      CRKL <- ARKL/VRKL            # ug/L, Rest of kidney lumen
      CVRK <- CRK/PK               # ug/L, Rest of kidney venous
      
      CSk <- ASk/VSk               # ug/L, Skin
      CVSk <- CSk/PSk              # ug/L, Skin venous 
      
      CGL <- AGL/VGL               # ug/L, Gut lumen
      # CGL_up <- AGL_up/VGL_up      # ug/L, Gut lumen, duodenum + jejunum
      # CGL_low <- AGL_low/VGL_low   # ug/L, Gut lumen, ileum
      # CGL_col <- AGL_col/VGL_col   # ug/L, Gut lumen, colon
      CG <- AG/VG                  # ug/L, Gut 
      CVG <- CG/PG                 # ug/L, Gut venous
      
      CL_ec <- AL_ec/VL_ec         # ug/L, Liver extracellular
      CL_ic <- AL_ic/VL_ic         # ug/L, Liver intracellular
      CVL_ec <- CL_ec/PL           # ug/L, Liver extracellular venous
      
      CA <- AA/VA                  # ug/L, Adipose
      CVA <- CA/PA                 # ug/L, Adipose venous
      
      CR <- AR/VR                  # ug/L, Rest
      CVR <- CR/PR                 # ug/L, Rest venous
      
      CP <- AP/VP                  # ug/L, Plasma
     
      
      ## Differential equations -------------------------
      dD = - D                     # ug/d, Dose input
      
      
      dAPT <- QK*(CP - CPT) + 
        + (Vmax_OAT4/(Km_OAT4+(CPTL*fuPTL)))*CPTL*fuPTL #+   
        # + CLdif_PTLtPTC*CPTL*fuPTL +
        # - CLdif_PTCtPTL*CPT*fuT                              # ug/d, Proximal tubule tissue 
      
      dAPTL <- + fup*GFR*CP - QT*CPTL +
        - (Vmax_OAT4/(Km_OAT4+(CPTL*fuPTL)))*CPTL*fuPTL #+
        # - CLdif_PTLtPTC*CPTL*fuPTL +
        # + CLdif_PTCtPTL*CPT*fuT #QT*CPTL QUr*CPTL            # ug/d, Proximal tubule lumen    
      
      dARK <- QK*(CPT - CVRK)                                  # ug/d, Rest of kidney
      dARKL <- QT*CPTL - QUr*CRKL                              # ug/d, Rest of kidney lumen
      dAUr <- QUr*CRKL #CRKL                                   # ug/d, Urine
      
      
      dASk <- QSk*(CP-CVSk) #+ Dermaldose                      # ug/d, Skin
      
      
      dAGL <- + D - tco*AGL - CL_GLtG*CGL + 
        + (VmaxBSEP/(KmBSEP + (CL_ic*fuT)))*CL_ic*fuT          # ug/d, Gut lumen 
      
      # dAGL_up <- 0 #D - CL_GLtG*CGL_up + (VmaxBSEP/(KmBSEP + (CL_ic*fuT)))*CL_ic*fuT - tsi_up*AGL_up
      # dAGL_low <- 0 #tsi_up*AGL_up - CL_GLtG*CGL_low - tsi_low*AGL_low
      # dAGL_col <- 0 #tsi_low*AGL_low - CL_GLtG*CGL_col - tco*AGL_col
      
      dAG <- QG*(CP - CVG) + CL_GLtG*CGL                       # ug/d, Gut
      
      dAFe <-  tco*AGL #tco*AGL_col                            # ug/d, Feces
      
      
      dAL_ec <- QG*CVG + QL*CP - (QG+QL)*CVL_ec +
        - (Vmax_OATP1B1/(Km_OATP1B1 + (CL_ec*fup)))*CL_ec*fuL_ec +
        - (Vmax_OATP1B3/(Km_OATP1B3 + (CL_ec*fup)))*CL_ec*fuL_ec            # ug/d, Liver extracellular space (vascular + interstitial space)
      
      dAL_ic <- (Vmax_OATP1B1/(Km_OATP1B1 + (CL_ec*fup)))*CL_ec*fuL_ec +
        + (Vmax_OATP1B3/(Km_OATP1B3 + (CL_ec*fup)))*CL_ec*fuL_ec +
        - (VmaxBSEP/(KmBSEP + (CL_ic*fuT)))*CL_ic*fuT                       # ug/d, Liver intracellular space
      
      
      dAA <- QA*(CP-CVA)                                      # ug/d, Adipose
      
      dAR <- QR*(CP-CVR)                                      # ug/d, Rest
      
      dAP <- - (QSk + QG + QL + QA + QR)*CP +
        + QSk*CVSk + (QL+QG)*CVL_ec + QA*CVA + QR*CVR +     
        - QK*CP + (QK*CVRK) +
        - fup*GFR*CP                                          # ug/d, Plasma
      
      
      # Mass Balance
      Atot <- D + APT + APTL + ARK + ARKL + AUr +
        ASk +
        AG + AGL + #AGL_up + AGL_low + AGL_col + AFe + 
        AL_ec + AL_ic +
        AA + AR +
        AP 
      
      # dInput <- Oraldose + Dermaldose
      MB = Oraldose - Atot + 1
      
      # End
      
      list(c(dD,
             dAPT,
             dAPTL,
             dARK,
             dARKL,
             dAUr,
             dASk, 
             dAGL,
             dAG, 
             dAFe,
             dAL_ec,
             dAL_ic,
             dAA,
             dAR, 
             dAP
             ), 
           c(CPT = CPT,
             CPTL = CPTL,
             CRK = CRK,
             CRKL = CRKL,
             CSk = CSk, 
             CGL = CGL,
             CG = CG, CVG = CVG, 
             CL_ec = CL_ec,
             CL_ic = CL_ic,
             CVL_ec = CVL_ec,
             CA = CA, CVA = CVA,
             CR = CR, CVR = CVR,
             CP = CP, 
             Atot = Atot, 
             MB = MB
             )
           )
    })
  }
  
  ## Initials ####
  
  A_init <- c(D = Oraldose,
              APT = 0,
              APTL = 0,
              ARK = 0,
              ARKL = 0,
              AUr = 0,
              ASk = 0, 
              AGL = 0,
              AG = 0,
              AFe = 0,
              AL_ec = 0,
              AL_ic = 0,
              AA = 0, 
              AR = 0, 
              AP = 0
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
  
  ## Mass Balance ###
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
  
  
  ## Concentration over time plots ####
  
  Plot_PFOA_Plasma <- ggplot()+
    geom_path(data = output.PFOA.df, aes(x = Days, y = CP))+
    theme_CP()+
    ylab("Plasma (ng/ml)")
  Plot_PFOA_Plasma
  ggsave("PlasmaConcentration.png", dpi = 300)
  
  Plot_PFOA_All <- output.PFOA.df %>% 
    mutate(CK = CPTL+CRKL+CPT+CRK,
           CL = CL_ec + CL_ic
           ) %>% 
    select(Days, CK, CSk, CL, CG, CA, CR, CP) %>% 
    rename(Kidney = CK, Skin = CSk, Liver = CL, Gut = CG, Adipose = CA, Rest = CR, Plasma = CP) %>% 
    pivot_longer(names_to = "Organ", values_to = "Concentration", Kidney:Plasma) %>% 
    ggplot()+
    geom_path(aes(x = Days, y = Concentration, color = Organ)) +
    facet_wrap(~ Organ)+
    theme_CP()+
    theme(legend.position = "none")+
    ylab("Organ Concentration (ng/ml)")
  Plot_PFOA_All
  ggsave("OrganConcentrations.png", dpi = 300)
  
  
  ## AUC and Half life ####

  AUC <- trapz(output_PFOA[ , "time"], output_PFOA[ , "CP"])   # ug*day/L
  
  # Calculate predicted half-life
  time <- output_PFOA[ , "time"]                               # days
  conc <- output_PFOA[ , "CP"]                                 # ug/L or ng/ml
  Cmax <- max(conc)
  Tmax <- time[which.max(conc)]
  tlast <- max(time[conc > 0])
  
  half_life <- pk.calc.half.life(
      conc,
      time,
      Tmax,
      tlast
    )
  
  HalfLife <- half_life$half.life/365                          # half-life in years
  
  # Experimental
  ExpData <- read_csv("C:/Users/pacho003/OneDrive - Wageningen University & Research/CP_L_R/PARC_PFOA_mechanistic/Input/HalfLifes.csv")
  
  Exp_HalfLifes <- ExpData %>%
    filter(species == "human",
           chemical == "pfoa",
           parameter== "HalfLife") %>%
    select(value_average) %>%
    rename(HalfLife = value_average)
  Exp_HalfLifes$HalfLife <- as.numeric(Exp_HalfLifes$HalfLife) # years
  
  Predicted.df <- data.frame(
    HalfLife = HalfLife,
    Origin = "Predicted",
    value = 1)
  Experimental.df <- data.frame(
    HalfLife = Exp_HalfLifes$HalfLife,
    Origin = "Experimental",
    value = 1)
  
  HalfLifes <- rbind(Predicted.df, Experimental.df)
  
  Plot_HalfLifes <- HalfLifes %>%
    ggplot()+
    geom_violin(data = Experimental.df, aes(value, HalfLife),
                color = "transparent",
                fill = "grey")+
    geom_point(data = Predicted.df, aes(value, HalfLife),
               color = "slateblue3", size = 5, shape = 18) +
    ylab("Half life (years)") +
    theme_CP() +
    theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x = element_blank()
        )
  Plot_HalfLifes
  ggsave("ExpVsSimHalfLife.png", dpi = 300)
  
  
  ## Experimental Vs Simulated ####
  ExpPlasma <- read_excel("C:/Users/pacho003/OneDrive - Wageningen University & Research/CP_L_R/PARC_PFOA_mechanistic/Input/Experimental.Plasma.PFOA.xlsx", 
                          col_types = c("numeric", "numeric"))
  
  
  ExpPlasma <- ExpPlasma %>% 
    rename(Days = Time_days) %>% 
    rename(CP = MPFOA_µg_per_L) %>%    # ug/L or ng/ml
    mutate(CP = CP - 0.130) %>%        # substracting the pre-existing level of 0.130ug/L from their previous study, as also done in the ref. article: https://doi.org/10.1016/j.envint.2024.109047 (table 3)
    filter(Days <=TSTOP)
  
  
  Plot_Plasma <- ggplot()+
    geom_path(data = output.PFOA.df, aes(x = Days, y = CP), color = "aquamarine", linewidth = 1.5)+
    geom_point(data = ExpPlasma, aes(x = Days, y = CP), color = "black")+
    theme_CP()+
    ylab("Plasma (ng/ml)")
  Plot_Plasma
  ggsave("PlasmaExpVsPredicted.png", dpi = 300)
  
  
  print(AUC)
  print(HalfLife)