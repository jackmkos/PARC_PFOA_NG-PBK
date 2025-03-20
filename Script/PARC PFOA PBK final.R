  # --------------------------------------------------------------------------- #
  # PBK MODEL FOR PFOA, TO BE USED TOGETHER WITH THE LATEST HBM DATA
  # Model File
  # CP, 20-03-2025
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
  
  EXP_STOP <- 50 #*365      # days, duration of the exposure; 1, Abraham et al. 2024 https://doi.org/10.1016/j.envint.2024.109047
  SIM_STOP <- 450 # 5000 #*365 #450  # days, duration of the simulation; 450 days follow-up period from Abraham et al. 2024 https://doi.org/10.1016/j.envint.2024.109047  
  
  TSTART <- 0
  TSTOP <- SIM_STOP                # days
  DT <- 1/10 #1/10                    # days
  TIME <- seq(TSTART,TSTOP,by=DT)
  
  
  ## Dermal exposure ##
  # Comment Chrysa: the below section can be used to calculate the dermal concentration based on cosmetic product type and PFOA concentrations in the cosmetic product
  # q = 123.20         # mg/kg/day amount of cosmetic product applied per day, SCCS 2021 table 3 (https://health.ec.europa.eu/document/download/89af1a70-a2b1-44da-a868-e7d80a8e736c_en?filename=sccs_o_250.pdf)
  # fret = 1           # fraction of the cosmetic product retained on the skin, SCCS 2021 table 3 (https://health.ec.europa.eu/document/download/89af1a70-a2b1-44da-a868-e7d80a8e736c_en?filename=sccs_o_250.pdf)
  # Aproduct = 1       # ug/mg amount of PFOA in the cosmetic product
  # fbac = 0.78        # fraction bioaccessible (fraction of PFOA that leaves the cosmetic product and is accessible for absorption) Namazkar et al 2024 DOI: 10.1039/D3EM00461A 
  # CDermal = q*fret*Aproduct*fbac # ug/kd/day
  CDermal = 0# 0.000542   # ug/kg/day, mean of #as.numeric(SumExpPFOA_LB_val[i,14])
  
  
  
  ## Oral exposure ##
  COral = 0.048   #0.000187 #0.048       # ug/kg/day, calculated back from Abraham et al. 2024 https://doi.org/10.1016/j.envint.2024.109047 (3.96/BW of 82Kg); 0.000187 # ug/kg/day [EFSA 2020, page 143] 3.96ug
  DOral = COral*82 #*DoseOn         # ug, PFOA oral dose 
  
  
  Tinput = 1          # day, duration of dose day
  tinterval = 1       # day, dose interval 
  
  
  # PBK MODEL PARAMETERS ####
  # ---------------------------------------------------------------------------- #
  
  ## Physiological ####
  Physio_params <- Physio_params %>% 
    filter(age == 60) %>%  # 60 years old adult, male, as in the Abraham study
    select(ends_with('_M')) %>% 
    mutate(BloodFlowSum = rowSums(select(., starts_with("Q_")))) %>% # 0.9935, total blood flow as the sum of the fractional blood flows of all organs on which we have data
    mutate(VolumesSum = rowSums(select(., starts_with("V_")))) # 0.96, total volume as the sum of the fractional organ volumes of all organs on which we have data
  
  BW <- 82                         # L (kg), The body weight of the subject in Abraham et al. 2024 https://doi.org/10.1016/j.envint.2024.109047 Physio_params$BDW_M 
  QC <- Physio_params$CardOut_M    # L/d, This is corrected for hematocrit already so it's plasma
  Hct <- 46.7/100                  # The hematocrit of the subject in Abraham et al. 2024 https://doi.org/10.1016/j.envint.2024.109047 Physio_params$Hct_M
  
  
  ### Organ volumes -------------------------
  
  # Skin
  VSkc <- Physio_params$V_skinFraction_M  
  SA_SkB = 15670                   # cm2, total body surface area except head, SCCS 2021 table 4 (https://health.ec.europa.eu/document/download/89af1a70-a2b1-44da-a868-e7d80a8e736c_en?filename=sccs_o_250.pdf)
  H_SkB = 83.1                     # cm, average thickness of the skin barrier, 83.7 +- 16.6 J.Sandby-Moller et al. 2003, Table II DOI: 10.1080/00015550310015419
  VSkB = SA_SkB*H_SkB * 1e-3       # L, volume of the skin barrier 
  # # Notes Chrysa: Trine calculated the surface area based on the BodyWeight with this formula: SA_SkB = 9.1*(BW*1000)^0.666  # cm2, total body area of the skin (Husoy)
  # # Trine's value provides the surface total surface area, but according to the SCCS, if we take the body lotion,then we should remove the surface area of the head from our calculations
  # # The general equation is SA = 0.02350*H^0.4226*W^0.51456, ref. https://www.rivm.nl/bibliotheek/rapporten/090013003.pdf
  
  # Stomach
  VStc <- Physio_params$V_stomachFraction_M
  SA_St <- 525                     # cm2, ICRP 2006, table 7.3, page 100
  
  # Intestine
  VIc <- Physio_params$V_gutFraction_M
  
  # Intestinal lumen is taken as a compartment outside the intestine (weight of Intestinal lumen is outside of the bodyweight)
  L = 280                           # cm, (adult of 70kg) Willmann2004 doi: 10.1021/jm030999b
  R = (1.75+1)/2                    # cm, (adult of 70kg mean value) Willmann2004 doi: 10.1021/jm030999b
  VIL = pi*L*(R^2)/1000             # L, cm3/1000 #Value is the same as Punt et al. 2021 https://dx.doi.org/10.1021/acs.chemrestox.0c00307
  VILc = VIL/70                     # deriving the constant as from the paper they used the default BW of 70kg
  
  SA_SI = 2*pi*R*L*25               # cm2, amplification factor of 25 for the microvilli in the intestinal lumen, Willmann2004 doi: 10.1021/jm030999b
  
  
  # Liver
  VLc <- Physio_params$V_liverFraction_M           
  
  # # Liver is divided in intracellular and extracellular compartments. 
  # # Extracellular compartment combines both the vascular and interstitial space
  # # Fractional volume of intracellular space was taken from Utsey et al. 2020 https://doi.org/10.1124/dmd.120.090498, https://github.com/metrumresearchgroup/PBPK_PC/blob/master/data/unified_tissue_comp.csv
  VL_icc = 0.573                   # Fractional volume of intracellular space in the liver
  VL_ecc = 1 - 0.573               # Fractional volume of extracellular space in the liver
  
  
  # Kidney
  VKc <- Physio_params$V_kidneyFraction_M
  
  # Kidney is divided in to proximal tubule and rest of kidney. Each compartment is divided in to tissue and lumen (or filtrate). 
  # This separation is necessary because of the active reabsorption happening in the proximal tubule, but also due to the differences in composition and pH of primary and terminal urine.
  # Fractional volumes were recalculated from Pletz et al. 2020 https://doi.org/10.1016/j.comtox.2021.100172
  
  VPTc <- 0.398                   # Fractional volume proximal tubule 
  VPTBc <- 0.137                  # Fractional volume proximal tubule blood
  VPTCc <- 0.460                  # Fractional volume proximal tubule cell
  VPTTc <- VPTBc*(1-Hct) + VPTCc  # Fractional volume proximal tubule tissue
  VPTLc <- 0.403                  # Fractional volume proximal tubule lumen
  
  VRKLc <- 0.283                   # Fractional volume rest of kidney lumen
  
  
  # Adipose
  VAc <- Physio_params$V_adiposeFraction_M
  MAc <- Physio_params$AdiposeMass_M
  
  # Lungs
  VLuc <- Physio_params$V_lungFraction_M
  
  # Plasma
  VPc <- Physio_params$V_plasmaFraction_M
  VAPc <- 0.39*VPc                  # Arterial plasma
  VVPc <- 0.61*VPc                  # Venous plasma
  
  # Total body volume
  VTotc <- Physio_params$VolumesSum                  
  
  ### Organ blood flows -------------------------
  QSkc <- Physio_params$Q_skinFraction_M    
  QTotc <- Physio_params$BloodFlowSum
  QStc <- Physio_params$Q_stomachFraction_M
  QIc <- Physio_params$Q_gutFraction_M   
  QLc <- Physio_params$Q_liverFraction_M 
  QKc <- Physio_params$Q_kidneyFraction_M  
  QAc <- Physio_params$Q_adiposeFraction_M 
  
  
  ### Other physiological flows and constants -------------------------
  
  tge = 0.77*24               # /d, Gastric emptying time for solids (as PFAS were administered in a muffin) this is about 75-80min in men and 100-110 min in women [ICRP 2006, page 84] (for liquids the gastric transit time is 10-60minutes willmann2004 doi: 10.1021/jm030999b)
  tco = 0.14*24               # /d, Bowel residence/transit time in the colon, (24*/h) or 7h willmann2004 doi: 10.1021/jm030999b
  
  QUrc = 0.022                # L/d, Urine flow rate to the bladder 22 mL/kg BW/d [ICRP 89 page 161]
  GFRc = 0.18                 # L/d 18% of total renal plasma flow [ICRP 89 page 159] http://www.icrp.org/publication.asp?id=ICRP%20Publication%2089
  QT = 43.2*60*24/1000        # L/d, Tubular flow rate at the end of the proximal tubule, 43.2 ml/min TFR from Scotcher et al. 2016 https://doi.org/10.1016/j.ejps.2016.03.018 and Pletz
  
  
  ### Physiological pHs in different matrices -------------------------
  
  pH_P <- 7.4     # plasma
  pH_IL <- 7       # intestinal lumen, average
  pH_S <- 2 # stomach, fasted (between 1.5 and 2.5)
  
  ### Albumin concentrations in different matrices -------------------------
  
  # From Akihiro Tojo and Satoshi Kinugasa 2012 doi:10.1155/2012/481520
  # In the same paper: the proximal tubule reabsorbes 71% of albumin, while LoH and DT 23% and the collecting duct 3%
  Calb_P <- 37.0      # mg/ml plasma
  Calb_PTL <- 14.4e-3 # mg/ml proximal tubule 
  Calb_RKL <- 2.88e-3 # mg/ml rest of kidney tubule, as the PT is filtering 71%, then the concentration of albumin leaving the PT is 6.641ul/ml, in the urine it's 0.7ul/ml, therefore I'm doing the average here
  Calb_exp <- 1e-10   # albumin, or protein concentration not reported in the experiments, therefore assuming a very low number
  
  # Based on the Poulin and Theil 2009, below Table 6
  # Albumin ratio
  R_T <- 0.5 # albumin and lipoprotein ratio between the tissue interstitial fluid and plasma
  R_PTL <- Calb_PTL/Calb_P 
  R_L_ec <- 0.086 # albumin ratio liver, Utsey et al. 2020 https://doi.org/10.1124/dmd.120.090498, https://github.com/metrumresearchgroup/PBPK_PC/blob/master/data/unified_tissue_comp.csv
  
  
  ## Chemical Specific ####
  
  MW <- PFOA_params$MW 
  fup <- PFOA_params$fup        # fraction unbound in plasma from Fisher et al.2024
  pKa <- 1.886                  # pKa of PFOA, average of experimental values from https://pfas-1.itrcweb.org, Table 4-1 excel file, pKa sheet
  
  
  ### Partition coefficients -------------------------
  
  # Calculate Plasma/Rest of the body partition coefficient
  # KpRe = partition coefficient of each of the lumped organs * fractional volume of the respective organ / sum of the fractional volume of all these organs
  KpRe <- (PFOA_params$KpBr*Physio_params$V_brainFraction_M +
             PFOA_params$KpHe*Physio_params$V_heartFraction_M + 
             PFOA_params$KpMu*Physio_params$V_muscleFraction_M +
             PFOA_params$KpSp*Physio_params$V_spleenFraction_M +
             PFOA_params$KpGo*Physio_params$V_reproFraction_M +
             PFOA_params$KpBo*Physio_params$V_boneFraction_M) / (
               Physio_params$V_brainFraction_M +
                 Physio_params$V_heartFraction_M +
                 Physio_params$V_muscleFraction_M +
                 Physio_params$V_spleenFraction_M +
                 Physio_params$V_reproFraction_M +
                 Physio_params$V_boneFraction_M)  
  PFOA_params$KpRe <- KpRe
  
  # Choose between calculated partition coefficients or initial ones (rat)
  # In PFOA_params, PF, PG, PK, PL, PSK and PR are the initial PCs from Kudo 2007 (rat). These were recalculated to total tissue (/fup) to enable differentiating the fup for the sensitivity analysis
  # Correcting for fraction unbound (as it was not incorporated in the input calculating file)
  PSkc <- PFOA_params$KpSk #PSk  # Skin
  PStc <- PFOA_params$KpSt       # Stomach
  PIc <- PFOA_params$KpIn  #PG   # Intestine
  PLc <- PFOA_params$KpLi  #PL   # Liver
  PKc <- PFOA_params$KpKi  #PK   # Kidney
  PAc <- PFOA_params$KpAd  #PF   # Adipose
  PRc <- PFOA_params$KpRe  #PR   # Rest
  PLuc <- PFOA_params$KpLu       # Lungs
  

  ### Uptake from skin -------------------------
  Papp_SkB = 3.82*1e-3 * 60*60                 # cm/s, Ragnarsdottir et al. 2024 https://doi.org/10.1016/j.envint.2024.108772 (calculations 3.82*1e-3 cm/h -> *60*60 cm/s)
  tlag = 1/6.21 * 24                           # /d, lag time for dermal absorption, Ragnarsdottir et al. 2024 https://doi.org/10.1016/j.envint.2024.108772 (calculations 6.21h *24d)
  
  ### Uptake from the gastro-intestinal duct -------------------------
  
  # Input data, in vitro clearance
  Papp_SI = 7.31*1e-6                         # cm/s, 7.31 ± 0.43, Janssen et al. 2024
  
  ### Uptake to the liver -------------------------
  
  # Input data, in vitro clearance
  Vmax_OATP1B1c = 2.305 * 1e-6                # umol/min/mg protein, 2.305± 0.295 pmol/min/mg protein [@lin2023]
  Km_OATP1B1c = 52.65                         # ug/L, 52.65 ± 23.28 uM [@lin2023]
  
  Vmax_OATP1B3c = 2.694 * 1e-6                # umol/min/mg protein, 2.694± 0.470 pmol/min/mg protein [@lin2023]
  Km_OATP1B3c = 91.6                          # ug/L, 91.61 ± 47.70 uM [@lin2023]
  
  
  # Relative expression factors
  OATP1B1_vitro = 0.120                       # [@lin2023, tables6, ref23]
  OATP1B1_vivo = 2.000                        # pmol/mg membrane protein [@lin2023, tables7, ref23]
  REF_OATP1B1 = OATP1B1_vivo/OATP1B1_vitro
  SF_OATP1B1 = REF_OATP1B1
  
  OATP1B3_vitro = 0.719                       # [@lin2023, tables6, ref23]
  OATP1B3_vivo = 1.000                        # pmol/mg membrane protein [@lin2023, tables7, ref23]
  REF_OATP1B3 = OATP1B3_vivo/OATP1B3_vitro
  SF_OATP1B3 = REF_OATP1B3
  
  ### Biliary clearance -------------------------
  
  # Input data, in vitro clearance
  VmaxBSEPc <- 7.1                            # umol/min/mg BSEP, Average active transport of bile acids, assuming that the maximum velocity of PFOA transport by BSEP corresponds to that of bile acids 
  KmBSEPc <- 16.4                             # ug/L, uM, Average affinity constant of bile acids to BSEP, following the above assumption
  
  # Scaling factor
  SF_BSEP <- 0.839*140000 * 99 * 1e-9 * 1e3   # Scaling factor for BSEP mediated hepatic efflux for GCA and GCDC De Bruijn et al. 2024 (calculation: SF_BSEP = aBSEP_all*MWBSEP_all * Hep_all * 1e-9 * 1e3 * VL )
  
  
  ### Renal Clearance -------------------------
  
  ## Active transport
  
  Vmax_OAT4c = 4.5                 # nmol/min/mg protein, Louisse et al. 2024 doi.org/10.1016/j.tox.2024.153961
  Km_OAT4c = 47                    # ug/L, scaled from uM, Louisse et al. 2024 doi.org/10.1016/j.tox.2024.153961
  REF_OAT <- 0.56                  # average of REF_OAT values of OAT1 and OAT3, check excel file for detailed information
  SF_OAT <- 0.17 * REF_OAT * 1e6   # 0.17 is the mg of protein per gram kidney, ICRP 89, Ruark et al 2020
  
  
  ### Final parameter constants -------------------------
  
  parm.c <- unlist(c(data.frame(BW,
                                QC,
                                Hct,
                                VSkB,
                                VSkc,
                                VStc,
                                SA_St,
                                VIc,
                                VILc,
                                SA_SI,
                                VLc,
                                VL_icc,
                                VL_ecc,
                                VKc,
                                VPTc, 
                                VPTTc, 
                                VPTLc, 
                                VRKLc,
                                VAc,
                                MAc,
                                VLuc,
                                VAPc,
                                VVPc,
                                QSkc, 
                                QStc,
                                QIc,
                                QLc,
                                QKc,
                                QAc,
                                QUrc,
                                GFRc,
                                QT,
                                tge,
                                tco,
                                R_T,
                                R_PTL,
                                R_L_ec,
                                fup,
                                pKa,
                                PSkc,
                                PStc,
                                PIc,
                                PLc,
                                PKc,
                                PAc,
                                PRc,
                                PLuc,
                                Papp_SkB,
                                SA_SkB,
                                Papp_SI,
                                Vmax_OATP1B1c,
                                Km_OATP1B1c,
                                Vmax_OATP1B3c,
                                Km_OATP1B3c,
                                REF_OATP1B1,
                                REF_OATP1B3,
                                VmaxBSEPc,
                                KmBSEPc,
                                SF_BSEP,
                                Vmax_OAT4c,
                                Km_OAT4c,
                                SF_OAT, 
                                EXP_STOP,
                                Tinput,
                                tinterval,
                                CDermal,
                                COral
  )))
  
  write.csv(parm.c, "GSA_parms.csv")
  
  
  
  # PBK MODEL ####
  # ---------------------------------------------------------------------------- #
  
  PBK.model <- function(t, state, parameters){
    with(as.list(c(state, parameters)), {
      
      ### Physiological ----
      
      VSkB <- VSkB                    # L, Volume of skin barrier
      VSk <- VSkc * BW                # L, Volume of skin
      
      VSt <- VStc * BW                # L, Volume of stomach
      
      VIL <- VILc * BW                # L, Volume of intestinal lumen
      VI <- VIc * BW                  # L, Volume of intestine
      
      VL <- VLc * BW                  # L, Volume of liver
      VL_ic <- VL_icc*VL              # L, Volume liver intracellular                 
      VL_ec <- VL_ecc*VL              # L, Volume liver extracellular
      
      VK <- VKc*BW                    # L, Volume of kidney
      VPT <- VPTc*VK                  # L, Volume of proximal tubule
      VPTT <- VPT*VPTTc               # L, Volume of proximal tubule tissue
      VPTL <- VPT*VPTLc               # L, Volume of proximal tubule lumen
      VRKL <- (VK - VPT)*VRKLc        # L, Volume of rest of kidney lumen
      VRKT <- VK - VPTT - VPTL - VRKL # L, Volume of rest of kidney tissue, used in the model
      
      VA <- VAc * BW/0.9 + MAc/0.9    # L, Volume of adipose
      
      VLu <- VLuc * BW                # L, Volume of lungs
      
      VAP <- VAPc * (1-Hct) * BW      # L, Volume of arterial plasma
      VVP <- VVPc * (1-Hct) * BW      # L, Volume of venous plasma
      
      VTotc <- 0.96
      VTot <- VTotc * BW              # L, Total body volume (used for mass balance)
      
      VR <- VTot - (VSk + VI + VL + VK + VA + VLu + VAP + VVP)             # L, Volume of the lumped rest compartment
      
      QTotc <- 0.988
      QSk <- QSkc/QTotc * QC                # L/d, Skin 
      QSt <- QStc/QTotc * QC                # L/d, Stomach
      QI <- QIc/QTotc * QC                  # L/d, Intestinal
      QL <- QLc/QTotc * QC                  # L/d, Liver
      QK <- QKc/QTotc * QC                  # L/d, Kidney
      QA <- QAc/QTotc * QC                  # L/d, Adipose 
      QR <- QC - (QA + QI + QK + QL + QSk)  # L/d, Rest
      
      QUr <- QUrc * BW              # L/d, Urine flow rate to the bladder 22 mL/kg BW/d [ICRP 89 page 161]
      GFR <- GFRc * QK              # L/d 18% of total renal plasma flow [ICRP 89 page 159] http://www.icrp.org/publication.asp?id=ICRP%20Publication%2089
      QT <- QT                      # L/d, Proximal tubule fluid flow
      
      tge <- tge                    # /d, Gastric emptying time
      tco <- tco                    # /d, Bowel residence time in the colon
      
      
      ### Physicochemical ----
      MW <- 414.07
      
      PSk <- PSkc * fup  #PSk  # Skin
      PSt <- PStc * fup        # Stomach
      PI <- PIc * fup    #PI   # Intestinal
      PL <- PLc * fup    #PL   # Liver
      PK <- PKc * fup    #PK   # Kidney
      PA <- PAc * fup    #PF   # Adipose
      PR <- PRc * fup    #PR   # Rest
      PLu <- PLuc * fup        # Lungs
      
      # Fraction unionised
      pH_P <- 7.4     # plasma
      pH_IL <- 7      # intestinal Intestinal lumen, average
      pH_S <- 2       # stomach, fasted (between 1.5 and 2.5)
      
      f.union_exp <- 1/(1 + 10^(pH_P - pKa))     # Is the same as plasma as pH in the experiment is 7.4
      f.union_IL <- 1/(1 + 10^(pH_IL - pKa))      # Intestinal lumen
      f.union_S <- 1/(1 + 10^(pH_S - pKa))      # Intestinal lumen
      
      # Fraction unbound
      fuT <- 1/(1 + ((1 - fup)/fup) * R_T)       # Tissue
      fuPTL <- 1/(1 + ((1 - fup)/fup) * R_PTL)   # Proximal tubule lumen
      fuL_ec <- 1/(1 + ((1 - fup)/fup) * R_L_ec) # Liver extracellular space
      # Note Chrysa: need to check if I find the amount of albumin in liver interstitial space, as here the fuL_ec is 10 times higher than the fup, but extracellular space is actually mainly the albumin from the vascular space. 
      
      ### Kinetic ----
      
      # Skin uptake
      CL_SkBtSk <- Papp_SkB*SA_SkB*1e-3*60*60*24           # L/d, Skin barrier to skin (calculations: cm/s -> L/s * 1e-3 -> L/d *60*60*24)
      
      # Gastro-intestinal uptake
      Pint_SI <- Papp_SI/f.union_exp                       # cm/s, Intrinsic permeability, corrected for fraction unionised in the experiment
      CL_GL <- (Pint_SI*SA_SI*f.union_IL*1e-3)*60*60*24    # L/d, Intestinal lumen to intestinal tissue (calculations: cm/s -> L/s /1000 -> L/d *60*60*24) 
      CL_St <- (Pint_SI*SA_St*f.union_S*1e-3)*60*60*24     # L/d, Stomach uptake to the portal vein
      
      # Liver uptake
      Vmax_OATP1B1 <- Vmax_OATP1B1c*REF_OATP1B1            # ug/d
      Km_OATP1B1 <- Km_OATP1B1c*MW                         # ug/L, 52.65 ± 23.28 uM [@lin2023]
      Vmax_OATP1B3 <- Vmax_OATP1B3c*REF_OATP1B3            # ug/d
      Km_OATP1B3 <- Km_OATP1B3c*MW                         # ug/L, 91.61 ± 47.70 uM [@lin2023]
      
      # Biliary excretion
      VmaxBSEP <- VmaxBSEPc*SF_BSEP*60*24*MW               # ug/d
      KmBSEP <- KmBSEPc*MW                                 # ug/L, uM, Average affinity constant of bile acids to BSEP, following the above assumption
      
      # Renal clearance
      Vmax_OAT4 = Vmax_OAT4c*MW*1e-3*60*24*SF_OAT*VPTT     # ug/d (nmol -> ug, min -> d)
      Km_OAT4 = Km_OAT4c*MW                                # ug/L, scaled from uM, Louisse et al. 2024 doi.org/10.1016/j.tox.2024.153961
      
      
      ## Dose -------------------------
      
      # if(t<EXP_STOP){DoseOn=1} else{DoseOn=0}
      # 
      # ## Dermal exposure ##
      # DDermal = CDermal*BW*DoseOn     # ug, PFOA dermal dose (what is available to be absorbed from the skin barrier)
      # DermalD <- DDermal/Tinput*(t %% tinterval<Tinput)
      # 
      # ## Oral exposure ##
      # DOral = COral*BW*DoseOn         # ug, PFOA oral dose
      # OralD = DOral #/Tinput*(t %% tinterval<Tinput)
      # 
      ## Concentrations -------------------------
      
      CSkB <- ASkB/VSkB            # ug/L, Skin barrier
      CSk <- ASk/VSk               # ug/L, Skin
      CVSk <- CSk/PSk              # ug/L, Skin venous 
      
      CSt <- ASt/VSt               # ug/L, Stomach
      CVSt <- CSt/PSt              # ug/L, Stomach venous
      
      CIL <- AIL/VIL               # ug/L, Intestinal lumen
      CI <- AI/VI                  # ug/L, Intestinal 
      CVI <- CI/PI                 # ug/L, Intestinal venous
      
      CL_ec <- AL_ec/VL_ec         # ug/L, Liver extracellular
      CL_ic <- AL_ic/VL_ic         # ug/L, Liver intracellular
      CVL_ec <- CL_ec/PL           # ug/L, Liver extracellular venous
      
      CPTT <- APTT/VPTT            # ug/L, Kidney, proximal tubule tissue
      CPTL <- APTL/VPTL            # ug/L, Kidney, proximal tubule lumen
      CRKT <- ARKT/VRKT            # ug/L, Rest of kidney tissue
      CRKL <- ARKL/VRKL            # ug/L, Rest of kidney lumen
      CVRKT <- CRKT/PK             # ug/L, Rest of kidney venous
      
      CA <- AA/VA                  # ug/L, Adipose
      CVA <- CA/PA                 # ug/L, Adipose venous
      
      CR <- AR/VR                  # ug/L, Rest
      CVR <- CR/PR                 # ug/L, Rest venous
      
      CLu <- ALu/VLu               # ug/L, Lungs
      CVLu <- CLu/PLu              # ug/L, Lung venous
      
      CAP <- AAP/VAP               # ug/L, Arterial plasma
      CVP <- AVP/VVP               # ug/L, Venous plasma
      
      
      ## Differential equations -------------------------
      
      dDD = 0 #DermalD - DD                             # ug/d, Dermal dose input
      dOD = - OD # OralD - OD                           # ug/d, Oral dose input
      
      
      dASkB <- DD - CL_SkBtSk*CSkB                                # ug/d, Skin barrier
      
      dASk <- + CL_SkBtSk*CSkB + QSk*(CAP-CVSk)                   # ug/d, Skin
      
      
      dASt <- + OD - CL_St*CVSt*QSt + QSt*CAP - tge*ASt       # ug/d, Stomach 
      
      
      dAIL <- + tge*ASt - tco*AIL - CL_GL*CIL + #+ tge*ASt
        + (VmaxBSEP/(KmBSEP + (CL_ic*fuT)))*CL_ic*fuT        # ug/d, Intestinal lumen
      
      dAI <- QI*(CAP - CVI) + CL_GL*CIL                      # ug/d, Intestinal
      
      dAFe <-  tco*AIL                                       # ug/d, Feces
      
      
      dAL_ec <- + QI*CVI + QL*CAP - (QI+QL+QSt)*CVL_ec + CL_St*CVSt*QSt +
        - (Vmax_OATP1B1/(Km_OATP1B1 + (CL_ec*fup)))*CL_ec*fuL_ec +
        - (Vmax_OATP1B3/(Km_OATP1B3 + (CL_ec*fup)))*CL_ec*fuL_ec            # ug/d, Liver extracellular space (vascular + interstitial space)
      
      dAL_ic <- (Vmax_OATP1B1/(Km_OATP1B1 + (CL_ec*fup)))*CL_ec*fuL_ec +
        + (Vmax_OATP1B3/(Km_OATP1B3 + (CL_ec*fup)))*CL_ec*fuL_ec +
        - (VmaxBSEP/(KmBSEP + (CL_ic*fuT)))*CL_ic*fuT                       # ug/d, Liver intracellular space
      
      
      dAPTT <- QK*(CAP - CPTT) + 
        + (Vmax_OAT4/(Km_OAT4+(CPTL*fuPTL)))*CPTL*fuPTL        # ug/d, Proximal tubule tissue 
      
      dAPTL <- + fup*GFR*CAP - QT*CPTL +
        - (Vmax_OAT4/(Km_OAT4+(CPTL*fuPTL)))*CPTL*fuPTL        # ug/d, Proximal tubule lumen    
      
      dARKT <- QK*(CPTT - CVRKT)                               # ug/d, Rest of kidney
      
      dARKL <- QT*CPTL - QUr*CRKL                              # ug/d, Rest of kidney lumen
      
      dAUr <- QUr*CRKL                                         # ug/d, Urine
      
      
      dAA <- QA*(CAP-CVA)                                      # ug/d, Adipose
      
      
      dAR <- QR*(CAP-CVR)                                      # ug/d, Rest
      
      dALu <- QC*(CVP - CVLu)
      
      dAAP <- - (QSk + QI + QL + QA + QR + QK + QSt)*CAP + QC*CVLu - fup*GFR*CAP    #+ QSt       # ug/d, Arterial Plasma
      dAVP <- + QSk*CVSk + (QL+QI+QSt)*CVL_ec + QK*CVRKT + QA*CVA + QR*CVR - QC*CVP      # ug/d, Venous Plasma
      
      # Mass Balance
      Atot <- DD + OD +
        ASkB + ASk +
        ASt +
        AIL + AI + AFe +   
        AL_ec + AL_ic +
        APTT + APTL + ARKT + ARKL + AUr +
        AA + 
        AR +
        ALu +
        AAP + AVP  
      
      dAin <- 0 #OralD + DermalD # to be used if repeated exposure
      # MB <- Ain - Atot + 1    # to be used if repeated exposure
      MB <- DOral - Atot + 1
      
      # End
      
      list(c(dDD,
             dOD,
             dASkB,
             dASk, 
             dASt,
             dAIL,
             dAI, 
             dAFe,
             dAL_ec,
             dAL_ic,
             dAPTT,
             dAPTL,
             dARKT,
             dARKL,
             dAUr,
             dAA,
             dAR, 
             dALu,
             dAAP, 
             dAVP,
             dAin
      ), 
      c(CSkB = CSkB, 
        CSk = CSk, 
        CSt = CSt,
        CIL = CIL,
        CI = CI, CVI = CVI, 
        CL_ec = CL_ec,
        CL_ic = CL_ic,
        CVL_ec = CVL_ec,
        CPTT = CPTT,
        CPTL = CPTL,
        CRKT = CRKT,
        CRKL = CRKL,
        CA = CA, CVA = CVA,
        CR = CR, CVR = CVR,
        CLu = CLu, CVLu = CVLu,
        CAP = CAP, CVP = CVP,
        CP = CVP,
        Atot = Atot, 
        MB = MB
      )
      )
    })
  }
  
  A_init <- c(DD = 0, OD = DOral, #0
              ASkB = 0, ASk = 0,
              ASt = 0,
              AIL = 0, AI = 0, AFe = 0,   
              AL_ec = 0, AL_ic = 0,
              APTT = 0, APTL = 0, ARKT = 0, ARKL = 0, AUr = 0,
              AA = 0, 
              AR = 0,
              ALu = 0,
              AAP = 0, AVP = 0,
              Ain = 0)
  
  output_PFOA <- lsoda(y = A_init, 
                       times = TIME, 
                       func = PBK.model, 
                       parms = parm.c)
  output.PFOA.df <- as.data.frame(output_PFOA) 
  

  
  # RESULTS ####
  # ---------------------------------------------------------------------------- #
  
  output.PFOA.df <- output.PFOA.df %>% 
    # mutate(time = time/365) %>% 
    rename(Days = time)
  
  ## Mass Balance ###
  MB.df <- output.PFOA.df %>% select(Days, Atot, MB)
  MB.df$MB <- round(MB.df$MB, 3)
  MB.df$ERROR <- (DOral - MB.df$Atot) / MB.df$Atot * 100
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
    mutate(CK = CPTT + CPTL + CRKT + CRKL,
           CL = CL_ec + CL_ic
           ) %>% 
    select(Days, CK, CSk, CL, CI, CSt, CA, CR, CP) %>% 
    rename(Kidney = CK, Skin = CSk, Liver = CL, Intestine = CI, Stomach = CSt, Adipose = CA, Rest = CR, Plasma = CP) %>% 
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