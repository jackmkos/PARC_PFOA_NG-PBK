# --------------------------------------------------------------------------- #
# SCRIPT FOR CALCULATING BASE PARAMETERS
# By: Chrysanthi Pachoulide
# Date: 14-04-2025
# --------------------------------------------------------------------------- #

# !! max BW in Physio_params is 76 for Female

# COMMON PARAMETERS ####
BASE_PARAMS <- function(expAGE = NULL, expBW = NULL, sex = NULL) { 
  
  ## Function conditions ####
  if (is.na(expAGE) && is.na(expBW)) {
    expBW <- 70  # Default body weight
  } 
  
  if (!is.na(expAGE) && expAGE > 80) {
    warning("Age above 80 which is the max age in the lifestage model. Physiological parameters assumed as those of an 80 year old")
    expAGE <- 80
  } 
  
  if (is.na(sex)) {
    sex <- "M"  # Default sex is male...
  }
  # 25ml plasma per cycle (29.2 days per cycle until; 12.5 cycles per year ); Ruark
  
  # Menstruation, reference Ruark et al 2016 http://dx.doi.org/10.1016/j.envint.2016.11.030
  if (sex == "F" && expAGE >= 12 && expAGE < 45){  
    CL_menses <-  25e-3/29.2 # L/day
  } else if (sex == "F" && expAGE >= 45 && expAGE <= 51){
    CL_menses <- 27.1e-3/35.1 # L/day
  } else {
    CL_menses <- 0 # L/day
  }
  
  # Serum albumin (g_HSA /L_serum, or mg/ml); Weaving et al 2016; DOI: 10.1177/0004563215593561
  if (expAGE <= 5) {
    SAlb <- 43  
  } else if (expAGE <= 15 && sex == "M") {
    SAlb <- 44
  } else if (expAGE <= 25 && sex == "M") {
    SAlb <- 46
  } else if (expAGE <= 35 && sex == "M") {  
    SAlb <- 45
  } else if (expAGE <= 50 && sex == "M") {
    SAlb <- 44
  } else if (expAGE <= 60 && sex == "M") {
    SAlb <- 43
  } else if (expAGE <= 75 && sex == "M") {
    SAlb <- 42
  } else if (expAGE <= 80 && sex == "M") {
    SAlb <- 41
  } else if (expAGE <= 90 && sex == "M") {
    SAlb <- 40
  } else if (expAGE <= 10 && sex == "F") {
    SAlb <- 45
  } else if (expAGE <= 20 && sex == "F") {
    SAlb <- 44
  } else if (expAGE <= 40 && sex == "F") {
    SAlb <- 43
  } else if (expAGE <= 75 && sex == "F") {
    SAlb <- 42
  } else if (expAGE <= 85 && sex == "F") {
    SAlb <- 41
  } else if (expAGE <= 90 && sex == "F") {
    SAlb <- 40
  } else if (expAGE <= 94) { 
    SAlb <- 39
  } else { 
    SAlb <- 38
  } 
  
  ## Input
  Physio_params <- Physio.c
  
  suffix <- ifelse(sex == "F", "_F", "_M")
  
  if (!is.na(expAGE)) { # filter based on age
    Physio_params <- Physio_params %>% 
      mutate(age_diff = abs(age - expAGE)) %>% 
      filter(age_diff == min(age_diff)) %>% 
      slice(1)
  } else if (!is.na(expBW)) { # filter based on BW if age is not provided
    bw_col <- paste0("BW", suffix)
    Physio_params <- Physio_params %>%
      mutate(bw_diff = abs(!!sym(bw_col) - expBW)) %>% 
      filter(bw_diff == min(bw_diff)) %>% 
      slice(1)
  }
  
  Physio_params <- Physio_params %>% select(ends_with(suffix))
  
  
  ## Physiological Parameters ####
  Physio_params <- Physio_params %>% 
    mutate(BloodFlowSum = rowSums(select(., starts_with("Q_")))) %>% # 0.9935, total blood flow as the sum of the fractional blood flows of all organs on which we have data
    mutate(VolumesSum = rowSums(select(., starts_with("V_")))) # 0.96, total volume as the sum of the fractional organ volumes of all organs on which we have data
  
  BW <- if (!is.na(expBW)) {
    expBW # BW is the actual BW at the time of measurement
  } else {
    Physio_params[[paste0("BW", suffix)]]}
  # BW <- Physio_params[[paste0("BW", suffix)]]   # Body weight calculated only from the lifestage equations
  
  QC <- Physio_params[[paste0("CardOut", suffix)]]              # L/d, This is corrected for hematocrit already so it's plasma
  Hct <- Physio_params[[paste0("Hct", suffix)]]                 # Hematocrit
  
  BH <- Physio_params[[paste0("BH", suffix)]]                   # Height (cm)
  BSA <- exp(-3.75 + 0.42*log(BH)+0.52*log(BW))*1e4             # Body Surface area (cm2)
  
  GFR <- Physio_params[[paste0("GFR", suffix)]]                 # Age/creatinine dependent glomerular filtration rate (L/day)
  
  ### Organ volumes -------------------------
  
  
  # Intestine
  VIc <- Physio_params[[paste0("V_gutFraction", suffix)]] 
  
  # Intestinal lumen is taken as a compartment outside the intestine (weight of Intestinal lumen is outside of the bodyweight)
  L = 280                           # cm, (adult of 70kg) Willmann2004 doi: 10.1021/jm030999b
  R = 1.75                          # cm, (adult of 70kg mean value) Willmann2004 doi: 10.1021/jm030999b
  VIL = pi*L*(R^2)/1000             # L, cm3/1000 #Value is the same as Punt et al. 2021 https://dx.doi.org/10.1021/acs.chemrestox.0c00307
  VILc = VIL/70                     # deriving the constant by dividing by BW, (adult of 70kg) Willmann2004 doi: 10.1021/jm030999b
  
  SA_SIc = 2*pi*R*L*25/70           # cm2, amplification factor of 25 for the microvilli in the intestinal lumen, Willmann2004 doi: 10.1021/jm030999b
  # SA_SIc = 70.1 * 1e3/70            # cm2 converted to a constant by dividing by BW, (adult of 70kg) Willmann2004 doi: 10.1021/jm030999b
  
  # Liver
  VLc <- Physio_params[[paste0("V_liverFraction", suffix)]]  # cm2, (adult of 70kg) Willmann2004 doi: 10.1021/jm030999b        
  
  # # Liver is divided in intracellular and extracellular compartments. 
  # # Extracellular compartment combines both the vascular and interstitial space
  # # Fractional volume of intracellular space was taken from Utsey et al. 2020 https://doi.org/10.1124/dmd.120.090498, https://github.com/metrumresearchgroup/PBPK_PC/blob/master/data/unified_tissue_comp.csv
  VL_icc = 0.573                   # Fractional volume of intracellular space in the liver
  VL_ecc = 1 - 0.573               # Fractional volume of extracellular space in the liver
  
  
  # Kidney
  VKc <- Physio_params[[paste0("V_kidneyFraction", suffix)]]
  
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
  VAc <- Physio_params[[paste0("V_adiposeFraction", suffix)]]
  MAc <- Physio_params[[paste0("AdiposeMass", suffix)]]
  
  
  # Plasma
  VPc <- Physio_params[[paste0("V_bloodFraction", suffix)]]
  
  
  # Total body volume
  VTotc <- Physio_params$VolumesSum                  
  
  ### Organ blood flows -------------------------
  QTotc <- Physio_params$BloodFlowSum
  QIc <- Physio_params[[paste0("Q_gutFraction", suffix)]]
  QLc <- Physio_params[[paste0("Q_liverFraction", suffix)]]
  QKc <- Physio_params[[paste0("Q_kidneyFraction", suffix)]] 
  QAc <- Physio_params[[paste0("Q_adiposeFraction", suffix)]]
  
  
  ### Other physiological flows and constants -------------------------
  
  tco = 0.14*24               # /d, Bowel residence/transit time in the colon, (24*/h) or 7h willmann2004 doi: 10.1021/jm030999b
  
  QUrc = 0.022                # L/d, Urine flow rate to the bladder 22 mL/kg BW/d [ICRP 89 page 161]
  GFRc = 0.18                 # L/d 18% of total renal plasma flow [ICRP 89 page 159] http://www.icrp.org/publication.asp?id=ICRP%20Publication%2089
  QTc = 0.36 # QT normalised to GFR, based on reference values of GFR 120 and QT of 43.2 from Scotcher et al. 2016 https://doi.org/10.1016/j.ejps.2016.03.018 and Pletz
  
  ### Physiological pHs in different matrices -------------------------
  
  pH_P <- 7.4     # plasma
  pH_IL <- 7      # intestinal lumen, average
  
  ### Albumin concentrations in different matrices -------------------------

  # From Akihiro Tojo and Satoshi Kinugasa 2012 doi:10.1155/2012/481520
  # In the same paper: the proximal tubule reabsorbes 71% of albumin, while LoH and DT 23% and the collecting duct 3%
  Calb_P <- SAlb      # g/L plasma (see above for values)
  Calb_PTL <- 14.4e-3 # g/L proximal tubule 
  Calb_IL <- 0.007e-3    # g/L intestinal lumen
  Calb_L_ec <- 0.0012e-3 # g/L Liver extracellular 

  # Based on the Poulin and Theil 2009, below Table 6
  # Albumin ratio
  R_PTL <- Calb_P/Calb_PTL # plasma to proximal tubule lumen albumin ratio
  R_L_ec <- Calb_P/Calb_L_ec # plasma to liver albumin ratio
  R_IL <- Calb_P/Calb_IL # plasma to intestinal lumen albumin ratio
  
  ## Chemical Specific ####
  
  # Distribution coefficients to different matrices 
  logML <- 3.52  # 3.52 ± 0.08 from Ebert A., Allendorf F. et al (2020), https://dx.doi.org/10.1021/acs.est.0c00175  (Liposomes composed of POPC (1-palmitoyl-2-oleoyl-glycero-3-phosphocholine))
  logSP <- 1.61  # 1.61 ± 0.15 from Allendorf, F., Goss, K.-U. and Ulrich, N. (2021), https://doi.org/10.1002/etc.4954  (Structural proteins from chicken breast fillet (actin & myosin 60-95%), Recovery 95%)
  logALB <- 4.48 # Fischer et al 2024 logALB <- 4.48±0.04 https://pubs.acs.org/doi/pdf/10.1021/acs.est.3c07415?ref=article_openPDF; 4.33 # 4.33 ± 0.05 from Allendorf, F., Goss, K.-U. and Ulrich, N. (2021), https://doi.org/10.1002/etc.4954 (BSA (fatty acid free) Molar ratio compound to BSA < 0.1 72-96h, Recovery 94%))
  logSL <- -1.37 # -1.37 ± 0.01 from Allendorf, F., Goss, K.-U. and Ulrich, N. (2021), https://doi.org/10.1002/etc.4954 (Olive oil with a high fraction of unsaturated fatty acids, Recovery 95%)
  logFABP <- 4.3 # calculated by Allendorf, F., Goss, K.-U. and Ulrich, N. (2021), https://doi.org/10.1002/etc.4954
  
  k_ML <- 10^logML       # neutral phospholipids:water partition coefficient
  k_SP <- 10^logSP       # structural protein:water partition coefficient
  k_ALB <- 10^logALB     # albumin:water partition coefficient
  k_SL <- 10^logSL       # structural lipids:water partition coefficient
  k_FABP <- 10^logFABP   # fatty acid-binding protein:water partition coefficient
  
  
  ### Fraction unbound in plasma -------------------------
  # Equation for calculating fraction unbound in plasma was taken from https://doi.org/10.1021/acs.est.3c07415
  
  DAlb <- 1.36 #kg/L density of albumin
  D_ALB <- k_ALB*DAlb # L_w/L_hsa; distribution coefficient
  FV_Alb <- SAlb*10^-3 # kg/kg of albumin in serum

  DGlob <- 1.20 #kg/L density of globulin
  k_GLOB <- 10^2.16 # Lw/Kg_globulin; globulin/water partition coefficient 
  D_Glob <- 10^2.16*DGlob # L_w/L_hsa; distribution coefficient
  FV_Glob <- 0.031 # kg/kg of globulin in serum
  
  FV_W <- 1 - FV_Alb - FV_Glob # serum mass balance
  
  fup <- 1/(1 + (D_ALB*(FV_Alb/FV_W)) + (D_Glob*(FV_Glob/FV_W)))
  
  ### Partition coefficients -------------------------
  # Code for calculating the partition coefficients is a slightly modified version from Utsey et al. 2020 https://github.com/metrumresearchgroup/PBPK_PC/blob/master/script/CalcKp_Schmitt.R

  Physio.data <- Tissue.c %>% 
    filter(Species == "Human") %>%
    filter(!Tissue %in% c("Comment", "NamingInUtsey")) %>%
    select(- c(Comment, Species)) %>% 
    mutate(across(-Tissue, as.numeric))
  
  Physio.Tissues <- Physio.data %>% filter(!Tissue %in% c("Blood", "Plasma")) 
  
  Kp <- (Physio.Tissues$f_W +
           (k_ML*Physio.Tissues$f_ML) +
           (k_SL*Physio.Tissues$f_SL) +
           (k_SP*Physio.Tissues$f_SP) +
           (k_ALB*Physio.Tissues$f_ALB) +
           (k_FABP*Physio.Tissues$f_FABP)
  ) #* fup Removing fup from here, Kps will be recalculated directly in the model for correct sensitivity analysis
  
  Names <- Physio.Tissues$Tissue %>% substr(1,2)
  Names <- paste("Kp",Names ,sep="")
  Kp.df <- as.data.frame(list(Tissue = Names)) %>%
    mutate(Value = Kp) %>% 
    pivot_wider(names_from = Tissue, values_from = Value)
  
  
  # Calculate Plasma/Rest of the body partition coefficient
  # KpRe = partition coefficient of each of the lumped organs * fractional volume of the respective organ / sum of the fractional volume of all these organs
  KpRe <- (Kp.df$KpSk * Physio_params[[paste0("V_skinFraction", suffix)]] +
             Kp.df$KpLu * Physio_params[[paste0("V_lungFraction", suffix)]] +
             Kp.df$KpSt * Physio_params[[paste0("V_stomachFraction", suffix)]] +
             Kp.df$KpBr * Physio_params[[paste0("V_brainFraction", suffix)]] +
             Kp.df$KpHe * Physio_params[[paste0("V_heartFraction", suffix)]] + 
             Kp.df$KpMu * Physio_params[[paste0("V_muscleFraction", suffix)]] +
             Kp.df$KpSp * Physio_params[[paste0("V_spleenFraction", suffix)]] +
             Kp.df$KpGo * Physio_params[[paste0("V_reproFraction", suffix)]] +
             Kp.df$KpBo * Physio_params[[paste0("V_boneFraction", suffix)]]) / (
               Physio_params[[paste0("V_skinFraction", suffix)]] +
               Physio_params[[paste0("V_lungFraction", suffix)]] +
                 Physio_params[[paste0("V_stomachFraction", suffix)]] +
                 Physio_params[[paste0("V_brainFraction", suffix)]] +
                 Physio_params[[paste0("V_heartFraction", suffix)]] +
                 Physio_params[[paste0("V_muscleFraction", suffix)]] +
                 Physio_params[[paste0("V_spleenFraction", suffix)]] +
                 Physio_params[[paste0("V_reproFraction", suffix)]] +
                 Physio_params[[paste0("V_boneFraction", suffix)]])  
  Kp.df$KpRe <- KpRe
  
  PIc <- Kp.df$KpIn   # Intestine
  PLc <- Kp.df$KpLi   # Liver
  PKc <- Kp.df$KpKi   # Kidney
  PAc <- Kp.df$KpAd   # Adipose
  PRc <- Kp.df$KpRe   # Rest
  
  ### Uptake from the gastro-intestinal duct -------------------------
  
  # Input data, in vitro clearance
  Papp_SI = 7.31 * 1e-6                       # cm/s, 7.31 ± 0.43, Janssen et al. 2024 http://dx.doi.org/10.1007/s00204-024-03851-x, epiIntestinal
  
  Vmax_OATP2B1c <- 0.001493                    # umol/min/mg protein, 1493 ± 543 pmol/min/mg protein; Lin et al. 2023 http://dx.doi.org/10.1021/acs.est.2c05642 (pmol -> umol)
  Km_OATP2B1c <- 148.68                        # umol/L, 148.68 ± 132.89 uM, Lin et al. 2023
  
  # IVIVE scaling factors
  REF_OATP2B1 = 1.6/0.812                     # relative expression factor OATP2B1, calculated from: OATP2B1_vivo/OATP2B1_vitro (pmol/mg membrane protein / pmol/mg membrane protein); Lin et al. 2023 http://dx.doi.org/10.1021/acs.est.2c05642 tables7 and s6, ref23
  mgOATP.mgI = 0.13                           # mg of protein per mg intestine; Utsey et al. 2020 http://dx.doi.org/10.1124/dmd.120.090498
  SF_OATP2B1 = mgOATP.mgI * REF_OATP2B1 * 1e6       # (mg liver * 1e6 -> kg liver)
  
  
  ### Uptake to the liver -------------------------
  
  # Input data, in vitro clearance
  Vmax_OATP1B1c = 0.002305                     # umol/min/mg protein, 2305 ± 295 pmol/min/mg protein, Lin et al. 2023 (pmol -> umol)
  Km_OATP1B1c = 52.65                          # umol/L, 52.65 ± 23.28 uM, Lin et al. 2023
  
  Vmax_OATP1B3c = 0.002694                     # umol/min/mg protein, 2694 ± 470 pmol/min/mg protein, Lin et al. 2023 (pmol -> umol)
  Km_OATP1B3c = 91.61                           # umol/L, 91.61 ± 47.70 uM, Lin et al. 2023
  
  # IVIVE scaling factors
  REF_OATP1B1 = 2/0.120                             # relative expression factor OATP1B1, calculated from: OATP1B1_vivo / OATP1B1_vitro  (pmol/mg membrane protein/pmol/mg membrane protein); Lin et al. 2023 http://dx.doi.org/10.1021/acs.est.2c05642 tables7 and s6, ref23
  REF_OATP1B3 = 1/0.719                             # relative expression factor OATP1B3, calculated from: OATP1B3_vivo / OATP1B3_vitro  (pmol/mg membrane protein/pmol/mg membrane protein); Lin et al. 2023 http://dx.doi.org/10.1021/acs.est.2c05642 tables7 and s6, ref23
  mgOATP.mgL = 0.18                                 # mg of protein per mg liver; Utsey et al. 2020 http://dx.doi.org/10.1124/dmd.120.090498
  SF_OATP1B1 = mgOATP.mgL * REF_OATP1B1 * 1e6       # (mg liver * 1e6 -> kg liver)
  SF_OATP1B3 = mgOATP.mgL * REF_OATP1B3 * 1e6       # (mg liver * 1e6 -> kg liver)
  
  
  ### Biliary clearance -------------------------
  # Input data, in vitro clearance
  Vmax_BSEPc <- 7                            # umol/min/mg BSEP, Average active transport of bile acids, assuming that the maximum velocity of PFOA transport by BSEP corresponds to that of bile acids ; de Bruijn et al. 2024 https://doi.org/10.14573/altex.2302011 
  Km_BSEPc <- 16.4                           # umol/L, uM, Average affinity constant of bile acids to BSEP, following the above assumption ; de Bruijn et al. 2024 https://doi.org/10.14573/altex.2302011
  
  # IVIVE 
  mgBSEP.HC <-  0.839 * 140000 * 1e-9         # mg BSEP/1e6 hepatocytes (0.839 pmole BSEP/1e6 hepatocytes (amound of BSEP per hepatocyte) * 140000 g/mole (MW BSEP) -> pg * 1e-9 -> mg) de Bruijn et al. 2024 https://doi.org/10.14573/altex.2302011]
  HC.GL <- 99                                 # 1e6 hepatocytes in the liver ; de Bruijn et al. 2024 https://doi.org/10.14573/altex.2302011
  SF_BSEP <- mgBSEP.HC * HC.GL * 1e3          # Scaling factor for BSEP mediated hepatic efflux for GCA and GCDC ; de Bruijn et al. 2024 https://doi.org/10.14573/altex.2302011
  
  
  ### Renal Clearance -------------------------
  
  ## Active transport
  
  Vmax_OAT4c = 4.5 *1e-3           # umol/min/mg protein, Louisse et al. 2024 doi.org/10.1016/j.tox.2024.153961 (nmol -> umol)
  Km_OAT4c = 47                    # ug/L, scaled from uM, Louisse et al. 2024 doi.org/10.1016/j.tox.2024.153961
  
  # IVIVE 
  REF_OAT4 <- 1.9/0.063                 # relative expression factor OAT,  calculated from: OAT4_vivo / OAT4_vitro  (pmol/mg membrane protein/ pmol/mg membrane protein); in vivo: Al-Madjoub et al 2021 (https://doi.org/10.1002/cpt.2396), in vitro:; Lin et al. 2023 http://dx.doi.org/10.1021/acs.est.2c05642 tables5
  mgOAT.mgK <- 0.17                     # mg of protein per mg kidney; Utsey et al. 2020 http://dx.doi.org/10.1124/dmd.120.090498
  SF_OAT <- mgOAT.mgK * REF_OAT4 * 1e6   # (mg kidney * 1e6 -> kg kidney)
  
  
  ### Final parameter constants -------------------------
  
  parm.c <- list(
   BW = BW,
   QC = QC,
   Hct = Hct,
   VIc = VIc,
   VILc = VILc,
   SA_SIc = SA_SIc,
   VLc = VLc,
   VL_icc = VL_icc,
   VL_ecc = VL_ecc,
   VKc = VKc,
   VPTc = VPTc, 
   VPTTc = VPTTc, 
   VPTLc = VPTLc, 
   VRKLc = VRKLc,
   VAc = VAc,
   MAc = MAc,
   VPc = VPc,
   QIc = QIc,
   QLc = QLc,
   QKc = QKc,
   QAc = QAc,
   QUrc = QUrc,
   GFR = GFR,
   QTc = QTc,
   tco = tco,
   R_PTL = R_PTL,
   R_L_ec = R_L_ec,
   fup = fup,
   PIc = PIc,
   PLc = PLc,
   PKc = PKc,
   PAc = PAc,
   PRc = PRc,
   Papp_SI = Papp_SI,
   Vmax_OATP2B1c = Vmax_OATP2B1c,
   Km_OATP2B1c = Km_OATP2B1c,
   SF_OATP2B1 = SF_OATP2B1,
   Vmax_OATP1B1c = Vmax_OATP1B1c,
   Km_OATP1B1c = Km_OATP1B1c,
   Vmax_OATP1B3c = Vmax_OATP1B3c,
   Km_OATP1B3c = Km_OATP1B3c,
   SF_OATP1B1 = SF_OATP1B1,
   SF_OATP1B3 = SF_OATP1B3,
   Vmax_BSEPc = Vmax_BSEPc,
   Km_BSEPc = Km_BSEPc,
   SF_BSEP = SF_BSEP,
   Vmax_OAT4c = Vmax_OAT4c,
   Km_OAT4c = Km_OAT4c,
   SF_OAT = SF_OAT,
   CL_menses = CL_menses
  )
  
  return(parm.c)
}

# PARAMETERS FOR DERMAL EXPOSURE ####
DERMAL_PARAMS <- function(expAGE = NULL, expBW = NULL, sex = "M", base.parm.c) { 
  
  ## Function conditions ####
  if (is.na(expAGE) && is.na(expBW)) {
    expBW <- 70  # Default body weight
  } 
  
  if (!is.na(expAGE) && expAGE > 80) {
    warning("Age above 80 which is the max age in the lifestage model. Physiological parameters assumed as those of an 80 year old")
    expAGE <- 80
  } 
  
  if (is.na(sex)) {
    sex <- "M"  # Default sex to male
  }
  
  ## Input
  Physio_params <- Physio.c
  
  suffix <- ifelse(sex == "F", "_F", "_M")
  
  
  if (!is.na(expAGE)) { # filters based on Age
    Physio_params <- Physio_params %>% 
      mutate(age_diff = abs(age - expAGE)) %>%  
      filter(age_diff == min(age_diff)) %>% 
      slice(1)
  } else if (!is.na(expBW)) { # filters based on BW if Age is not provided
    bw_col <- paste0("BW", suffix)
    Physio_params <- Physio_params %>%
      mutate(bw_diff = abs(!!sym(bw_col) - expBW)) %>% 
      filter(bw_diff == min(bw_diff)) %>% 
      slice(1)
  }
  
  Physio_params <- Physio_params %>% select(ends_with(suffix))
  
  
  BW <- if (!is.na(expBW)) {
    expBW
  } else {
    Physio_params[[paste0("BW", suffix)]]}
  # BW <- Physio_params[[paste0("BW", suffix)]]   # Body weight calculated only from the lifestage equations
  
  BH <- Physio_params[[paste0("BH", suffix)]]                   # Height (cm)
  BSA <- exp(-3.75 + 0.42*log(BH)+0.52*log(BW))*1e4             # Body Surface area (cm2), reference: Gastellu et al. 2024, 10.1016/j.envres.2024.120393 (supplementary file Physio_equations_detailed.xlsx, eq. from Pendse et al. 2020)
  # Trine calculated the surface area based on the BodyWeight with this formula: SA_SkB = 9.1*(BW*1000)^0.666  # cm2, total body area of the skin (Husoy)
  # Another option: BSA = 0.02350*BH^0.4226*BW^0.51456, ref. https://www.rivm.nl/bibliotheek/rapporten/090013003.pdf

  # Physiological Parameters ####
  Physio_params <- Physio_params %>% 
    mutate(BloodFlowSum = rowSums(select(., starts_with("Q_")))) %>% # 0.9935, total blood flow as the sum of the fractional blood flows of all organs on which we have data
    mutate(VolumesSum = rowSums(select(., starts_with("V_")))) # 0.96, total volume as the sum of the fractional organ volumes of all organs on which we have data
  
  
  # Skin
  VSkc <- Physio_params[[paste0("V_skinFraction", suffix)]]  # fractional skin volume
  SA_SkB <- BSA                       # cm2, total body surface area
  # SA_SkB <- 15670                   # cm2, total body surface area except head, SCCS 2021 table 4 (https://health.ec.europa.eu/document/download/89af1a70-a2b1-44da-a868-e7d80a8e736c_en?filename=sccs_o_250.pdf)

  H_SkB <- 83.1                     # cm, average thickness of the skin barrier, 83.7 +- 16.6 J.Sandby-Moller et al. 2003, Table II DOI: 10.1080/00015550310015419
  VSkB <- SA_SkB*H_SkB * 1e-3       # L, volume of the skin barrier 
  
  QSkc <- Physio_params[[paste0("Q_skinFraction", suffix)]]  # fractional skin blood flow
  
  # Partition coefficients -------------------------
  
  # Code for calculating the partition coefficients is a slidely modified version from Utsey et al. 2020 https://github.com/metrumresearchgroup/PBPK_PC/blob/master/script/CalcKp_Schmitt.R
  
  logML <- 3.52  # 3.52 ± 0.08 from Ebert A., Allendorf F. et al (2020), https://dx.doi.org/10.1021/acs.est.0c00175  (Liposomes composed of POPC (1-palmitoyl-2-oleoyl-glycero-3-phosphocholine))
  logSP <- 1.61  # 1.61 ± 0.15 from Allendorf, F., Goss, K.-U. and Ulrich, N. (2021), https://doi.org/10.1002/etc.4954  (Structural proteins from chicken breast fillet (actin & myosin 60-95%), Recovery 95%)
  logALB <- 4.33 # 4.33 ± 0.05 from Allendorf, F., Goss, K.-U. and Ulrich, N. (2021), https://doi.org/10.1002/etc.4954 (BSA (fatty acid free) Molar ratio compound to BSA < 0.1 72-96h, Recovery 94%))
  logSL <- -1.37 # -1.37 ± 0.01 from Allendorf, F., Goss, K.-U. and Ulrich, N. (2021), https://doi.org/10.1002/etc.4954 (Olive oil with a high fraction of unsaturated fatty acids, Recovery 95%)
  logFABP <- 4.3 # calculated by Allendorf, F., Goss, K.-U. and Ulrich, N. (2021), https://doi.org/10.1002/etc.4954
  
  k_ML <- 10^logML       # neutral phospholipids:water partition coefficient
  k_SP <- 10^logSP       # structural protein:water partition coefficient
  k_ALB <- 10^logALB     # albumin:water partition coefficient
  k_SL <- 10^logSL       # structural lipids:water partition coefficient
  k_FABP <- 10^logFABP   # fatty acid-binding protein:water partition coefficient
  
  Physio.data <- Tissue.c %>% 
    filter(Species == "Human") %>%
    filter(!Tissue %in% c("Comment", "NamingInUtsey")) %>%
    select(- c(Comment, Species)) %>% 
    mutate(across(-Tissue, as.numeric))
  
  Physio.Tissues <- Physio.data %>% filter(!Tissue %in% c("Blood", "Plasma")) 
  
  Kp <- (Physio.Tissues$f_W +
           (k_ML*Physio.Tissues$f_ML) +
           (k_SL*Physio.Tissues$f_SL) +
           (k_SP*Physio.Tissues$f_SP) +
           (k_ALB*Physio.Tissues$f_ALB) +
           (k_FABP*Physio.Tissues$f_FABP)
  ) #* fup Removing fup from here, Kps will be recalculated directly in the model for correct sensitivity analysis
  
  Names <- Physio.Tissues$Tissue %>% substr(1,2)
  Names <- paste("Kp",Names ,sep="")
  Kp.df <- as.data.frame(list(Tissue = Names)) %>%
    mutate(Value = Kp) %>% 
    pivot_wider(names_from = Tissue, values_from = Value)
  
  
  # Calculate Plasma/Rest of the body partition coefficient
  # KpRe = partition coefficient of each of the lumped organs * fractional volume of the respective organ / sum of the fractional volume of all these organs
  KpRe <- (Kp.df$KpLu * Physio_params[[paste0("V_lungFraction", suffix)]] +
             Kp.df$KpSt * Physio_params[[paste0("V_stomachFraction", suffix)]] +
             Kp.df$KpBr * Physio_params[[paste0("V_brainFraction", suffix)]] +
             Kp.df$KpHe * Physio_params[[paste0("V_heartFraction", suffix)]] + 
             Kp.df$KpMu * Physio_params[[paste0("V_muscleFraction", suffix)]] +
             Kp.df$KpSp * Physio_params[[paste0("V_spleenFraction", suffix)]] +
             Kp.df$KpGo * Physio_params[[paste0("V_reproFraction", suffix)]] +
             Kp.df$KpBo * Physio_params[[paste0("V_boneFraction", suffix)]]) / (
                 Physio_params[[paste0("V_lungFraction", suffix)]] +
                 Physio_params[[paste0("V_stomachFraction", suffix)]] +
                 Physio_params[[paste0("V_brainFraction", suffix)]] +
                 Physio_params[[paste0("V_heartFraction", suffix)]] +
                 Physio_params[[paste0("V_muscleFraction", suffix)]] +
                 Physio_params[[paste0("V_spleenFraction", suffix)]] +
                 Physio_params[[paste0("V_reproFraction", suffix)]] +
                 Physio_params[[paste0("V_boneFraction", suffix)]]) 
  Kp.df$KpRe <- KpRe
  
  # Choose between calculated partition coefficients or initial ones (rat)
  # In Kp.df, PF, PI, PK, PL, PSK and PR are the initial PCs from Kudo 2007 (rat). These were recalculated to total tissue (/fup) to enable differentiating the fup for the sensitivity analysis
  # Correcting for fraction unbound (as it was not incorporated in the input calculating file)
  PSkc <- Kp.df$KpSk #PSk  # Skin
  PRc <- Kp.df$KpRe  #PR   # Rest
  
  base.parm.c["PRc"] <- PRc
  
  
  # Uptake from skin -------------------------
  Papp_SkB = 3.82*1e-3 * 60*60                 # cm/s, Ragnarsdottir et al. 2024 https://doi.org/10.1016/j.envint.2024.108772 (calculations 3.82*1e-3 cm/h -> *60*60 cm/s)
  tlag = 1/6.21 * 24                           # /d, lag time for dermal absorption, Ragnarsdottir et al. 2024 https://doi.org/10.1016/j.envint.2024.108772 (calculations 6.21h *24d)
  
  # Final parameter constants -------------------------
  
  parm.c <- c(
    base.parm.c, 
    VSkc = VSkc, 
    SA_SkB = SA_SkB,
    VSkB = VSkB,
    QSkc = QSkc,
    PSkc = PSkc,
    Papp_SkB = Papp_SkB,
    tlag = tlag)
  
  
  # Function out -------------------------
  return(parm.c)
}

# PARAMETERS FOR INHALATION EXPOSURE ####
INHALATION_PARAMS <- function(expAGE = NULL, expBW = NULL, sex = "M", base.parm.c) { 
  
  
  ## Function conditions ####
  if (is.na(expAGE) && is.na(expBW)) {
    expBW <- 70  # Default body weight
  } 
  
  if (!is.na(expAGE) && expAGE > 80) {
    warning("Age above 80 which is the max age in the lifestage model. Physiological parameters assumed as those of an 80 year old")
    expAGE <- 80
  } 
  
  if (is.na(sex)) {
    sex <- "M"  # Default sex to male
  }
  
  ## Input
  Physio_params <- Physio.c
  
  suffix <- ifelse(sex == "F", "_F", "_M")
  
  
  if (!is.na(expAGE)) { # filters based on Age
    Physio_params <- Physio_params %>% 
      mutate(age_diff = abs(age - expAGE)) %>%  
      filter(age_diff == min(age_diff)) %>% 
      slice(1)
  } else if (!is.na(expBW)) { # filters based on BW if Age is not provided
    bw_col <- paste0("BW", suffix)
    Physio_params <- Physio_params %>%
      mutate(bw_diff = abs(!!sym(bw_col) - expBW)) %>% 
      filter(bw_diff == min(bw_diff)) %>% 
      slice(1)
  }
  
  Physio_params <- Physio_params %>% select(ends_with(suffix))
  

  # Physiological Parameters ####
  Physio_params <- Physio_params %>% 
    mutate(BloodFlowSum = rowSums(select(., starts_with("Q_")))) %>% # 0.9935, total blood flow as the sum of the fractional blood flows of all organs on which we have data
    mutate(VolumesSum = rowSums(select(., starts_with("V_")))) # 0.96, total volume as the sum of the fractional organ volumes of all organs on which we have data
  
  # Lungs
  VLuc <- Physio_params[[paste0("V_lungFraction", suffix)]]
  
  # Plasma
  VPc <- Physio_params[[paste0("V_plasmaFraction", suffix)]]
  VAPc <- 0.39*VPc                  # Arterial plasma
  VVPc <- 0.61*VPc                  # Venous plasma  
  
  # Partition coefficients -------------------------
  
  # Code for calculating the partition coefficients is a slidely modified version from Utsey et al. 2020 https://github.com/metrumresearchgroup/PBPK_PC/blob/master/script/CalcKp_Schmitt.R
  
  logML <- 3.52  # 3.52 ± 0.08 from Ebert A., Allendorf F. et al (2020), https://dx.doi.org/10.1021/acs.est.0c00175  (Liposomes composed of POPC (1-palmitoyl-2-oleoyl-glycero-3-phosphocholine))
  logSP <- 1.61  # 1.61 ± 0.15 from Allendorf, F., Goss, K.-U. and Ulrich, N. (2021), https://doi.org/10.1002/etc.4954  (Structural proteins from chicken breast fillet (actin & myosin 60-95%), Recovery 95%)
  logALB <- 4.33 # 4.33 ± 0.05 from Allendorf, F., Goss, K.-U. and Ulrich, N. (2021), https://doi.org/10.1002/etc.4954 (BSA (fatty acid free) Molar ratio compound to BSA < 0.1 72-96h, Recovery 94%))
  logSL <- -1.37 # -1.37 ± 0.01 from Allendorf, F., Goss, K.-U. and Ulrich, N. (2021), https://doi.org/10.1002/etc.4954 (Olive oil with a high fraction of unsaturated fatty acids, Recovery 95%)
  logFABP <- 4.3 # calculated by Allendorf, F., Goss, K.-U. and Ulrich, N. (2021), https://doi.org/10.1002/etc.4954
  
  k_ML <- 10^logML       # neutral phospholipids:water partition coefficient
  k_SP <- 10^logSP       # structural protein:water partition coefficient
  k_ALB <- 10^logALB     # albumin:water partition coefficient
  k_SL <- 10^logSL       # structural lipids:water partition coefficient
  k_FABP <- 10^logFABP   # fatty acid-binding protein:water partition coefficient
  
  Physio.data <- Tissue.c %>% 
    filter(Species == "Human") %>%
    filter(!Tissue %in% c("Comment", "NamingInUtsey")) %>%
    select(- c(Comment, Species)) %>% 
    mutate(across(-Tissue, as.numeric))
  
  Physio.Tissues <- Physio.data %>% filter(!Tissue %in% c("Blood", "Plasma")) 
  
  Kp <- (Physio.Tissues$f_W +
           (k_ML*Physio.Tissues$f_ML) +
           (k_SL*Physio.Tissues$f_SL) +
           (k_SP*Physio.Tissues$f_SP) +
           (k_ALB*Physio.Tissues$f_ALB) +
           (k_FABP*Physio.Tissues$f_FABP)
  ) #* fup Removing fup from here, Kps will be recalculated directly in the model for correct sensitivity analysis
  
  Names <- Physio.Tissues$Tissue %>% substr(1,2)
  Names <- paste("Kp",Names ,sep="")
  Kp.df <- as.data.frame(list(Tissue = Names)) %>%
    mutate(Value = Kp) %>% 
    pivot_wider(names_from = Tissue, values_from = Value)
  
  
  # Calculate Plasma/Rest of the body partition coefficient
  # KpRe = partition coefficient of each of the lumped organs * fractional volume of the respective organ / sum of the fractional volume of all these organs
  KpRe <- (Kp.df$KpSk * Physio_params[[paste0("V_skinFraction", suffix)]] +
             Kp.df$KpSt * Physio_params[[paste0("V_stomachFraction", suffix)]] +
             Kp.df$KpBr * Physio_params[[paste0("V_brainFraction", suffix)]] +
             Kp.df$KpHe * Physio_params[[paste0("V_heartFraction", suffix)]] + 
             Kp.df$KpMu * Physio_params[[paste0("V_muscleFraction", suffix)]] +
             Kp.df$KpSp * Physio_params[[paste0("V_spleenFraction", suffix)]] +
             Kp.df$KpGo * Physio_params[[paste0("V_reproFraction", suffix)]] +
             Kp.df$KpBo * Physio_params[[paste0("V_boneFraction", suffix)]]) / (
               Physio_params[[paste0("V_skinFraction", suffix)]] +
                 Physio_params[[paste0("V_stomachFraction", suffix)]] +
                 Physio_params[[paste0("V_brainFraction", suffix)]] +
                 Physio_params[[paste0("V_heartFraction", suffix)]] +
                 Physio_params[[paste0("V_muscleFraction", suffix)]] +
                 Physio_params[[paste0("V_spleenFraction", suffix)]] +
                 Physio_params[[paste0("V_reproFraction", suffix)]] +
                 Physio_params[[paste0("V_boneFraction", suffix)]])  
  Kp.df$KpRe <- KpRe
  
  # Choose between calculated partition coefficients or initial ones (rat)
  # Correcting for fraction unbound (as it was not incorporated in the input calculating file)
  PRc <- Kp.df$KpRe   # Rest
  PLuc <- Kp.df$KpLu  # Lungs
  
  base.parm.c["PRc"] <- PRc
  
  # Final parameter constants -------------------------
  
  parm.c <- c(
    base.parm.c, 
    VLuc = VLuc,
    VAPc = VAPc, 
    VVPc = VVPc,
    PLuc = PLuc)
  
  # Function out -------------------------
  return(parm.c)
}
