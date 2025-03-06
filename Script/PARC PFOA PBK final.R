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
  # Final_variables_M_df <- read_csv("C:/Users/pacho003/OneDrive - Wageningen University & Research/CP_L_R/PARC_PFOA_mechanistic/Input/2024-11-18/Final_variables_M.csv") %>% as.data.frame()
  # Final_variables_M_df <- read_csv("C:/Users/pacho003/OneDrive - Wageningen University & Research/CP_L_R/PARC_PFOA_mechanistic/Input/2024-11-18/Final_variables_MVolunteer.csv")
  Physio_params <- read_csv("C:/Users/pacho003/OneDrive - Wageningen University & Research/CP_L_R/PARC_PFOA_mechanistic/Input/PhysioVariables.csv")
  PFOA_params <- read_csv("C:/Users/pacho003/OneDrive - Wageningen University & Research/CP_L_R/PARC_PFOA_mechanistic/Input/PFOAParams.csv")
  
  
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
  VmaxBSEPc = 7.1 # umol/min/mg BSEP, average active transport of bile acids, assuming that the maximum velocity of PFOA transport by BSEP corresponds to that of bile acids 
  KmBSEP = 16.4 # uM, average affinity constant of bile acids to BSEP, following the above assumption
  
  aBSEP_all = 0.839 # pmoles/10^6 hepatocytes, BSEP protein abundance De Bruijn et al. 2024
  MWBSEP_all = 140000 # g/mole, MW BSEP De Bruijn et al. 2024
  Hep_all = 99 # 10^6 hepatocytes per g liver, hepatocellularity De Bruijn et al. 2024
  
  SF_BSEP = aBSEP_all*MWBSEP_all*Hep_all*1e-9*1e3*VL; # mg BSEP/entire liver; scaling factor for BSEP mediated hepatic efflux for GCA and GCDC De Bruijn et al. 2024
  VmaxBSEP = VmaxBSEPc*SF_BSEP*60*24*MW;   # ug/d   maximum   speed   for   BSEP-mediated GCDCA efflux (calculated) De Bruijn et al. 2024
  
  L = 280 #cm (adult of 70kg) willmann2004 doi: 10.1021/jm030999b
  R = (1.75+1)/2 #cm (adult of 70kg mean value) willmann2004 doi: 10.1021/jm030999b
  VGL = pi*L*(R^2)/1000 #L cm3/1000 #value is the same as Ans Punt 
  
  L_up = 124 #cm (adult of 70kg) duodenum 20, jejunum 104 cm) willmann2004 doi: 10.1021/jm030999b
  R_up = 1.75 #cm !amplification factor of 3 regarding the folds found in the lumen part surface area expansion, then added amplification factor of 25 for the microvilli willmann2004 doi: 10.1021/jm030999b
  VGL_up = pi*L_up*(R_up^2)/1000 #L cm3/1000  
  
  
  L_low = 156 #cm (adult of 70kg) (ileum) !amplification factor 1.5 as it drops to 1 in the distal part of the ileum willmann2004 doi: 10.1021/jm030999b
  R_low = 1 #cm willmann2004 doi: 10.1021/jm030999b
  VGL_low = pi*L_low*(R_low^2)/1000 #L cm3/1000 #value is the same as Ans Punt 
  
  
  L_col = 1.5 #m willmann2004 doi: 10.1021/jm030999b
  R_col = 3.5 #cm willmann2004 doi: 10.1021/jm030999b
  VGL_col = pi*L_col*(R_col^2)/1000 #L cm3/1000 #value is the same as Ans Punt 
  
  pH_L = 7
  pH_up = 6.5 # fasted (6.0 duodenum, 7.0 jejunum) willmann2004 doi: 10.1021/jm030999b
  pH_low = 7.5 # fasted willmann2004 doi: 10.1021/jm030999b
  
  
  tge = 2*24 #/d (24*/h) or 30min (10-60min) gastric emptying time fasted willmann2004 doi: 10.1021/jm030999b, Punt et al 2021 https://dx.doi.org/10.1021/acs.chemrestox.0c00307
  tsi = 0.3*24 #/d (24*/h) or 4h (2-6 h) willmann2004 doi: 10.1021/jm030999b (#transit time duodenum 14min, jejunum 71min, ileum 114min B. Agoram et al. 2001), Punt et al 2021 https://dx.doi.org/10.1021/acs.chemrestox.0c00307
  tsi_up = 0.7*24  #/d (24*/h) 1.417 h  ;14+17 min (#transit time duodenum 14min, jejunum 71min, ileum 114min B. Agoram et al. 2001), Punt et al 2021 https://dx.doi.org/10.1021/acs.chemrestox.0c00307
  tsi_low = 1.9*24 #/d (24*/h) 1.9 h;114 min (#transit time duodenum 14min, jejunum 71min, ileum 114min B. Agoram et al. 2001), Punt et al 2021 https://dx.doi.org/10.1021/acs.chemrestox.0c00307
  tco = 0.09*24 #/d (24*/h) or 7h willmann2004 doi: 10.1021/jm030999b
  
  # Passive permeability lumen to gut tissue
  # The ka way
  Papp_SI = 7.31 *1e-6 #cm/s 7.31 ± 0.43; Janssen et al. 2024
  log_Peff = 0.4926*log10(Papp_SI) - 0.1454 #equation for passively and actively absorbed compounds; Hou, Zhang et al. 2004
  ka = 10^(log_Peff)*2/R * 3600 * 24 #/d (24*/h); Yu and Amidon 1999
  
  # The Clearance way
  SA_SI = 2*pi*R*L*25 #cm2 !amplification factor of 3 regarding the folds found in the lumen part surface area expansion, then added amplification factor of 25 for the microvilli, !amplification factor 1.5 as it drops to 1 in the distal part of the ileum willmann2004 doi: 10.1021/jm030999b
  f.union_GT <- 1/(1 + 10^(pH_L - pKa))
  Pint_SI <- Papp_SI/f.union_exp #cm/s intrinsic permeability, corrected for fraction unionised in the experiment
  CL_GLtG <- (Pint_SI*SA_SI*f.union_GT/1000)*60*60*24 # L/d; cm/s = L/s /1000 = L/d *60*60*24; Proximal tubule lumen to proximal tubule cell
  
  VL_ic = 0.573
  VL_ec = 1 - 0.573
  
  Vmax_OATP1B1 = 2.305 * 1e-6 #umol/min/mg protein; 2.305± 0.295 pmol/min/mg protein [@lin2023]
  Km_OATP1B1 = 52.65 #52.65 ± 23.28 uM [@lin2023]
  OATP1B1_vitro = 0.120 #[@lin2023, tables6, ref23]
  OATP1B1_vivo = 2.000 #pmol/mg membrane protein [@lin2023, tables7, ref23]
  REF_OATP1B1 = OATP1B1_vivo/OATP1B1_vitro
  
  Vmax_OATP1B3 = 2.694 * 1e-6 #umol/min/mg protein; 2.694± 0.470 pmol/min/mg protein [@lin2023]
  Km_OATP1B3 = 91.61 #91.61 ± 47.70 uM [@lin2023]
  OATP1B3_vitro = 0.719 #[@lin2023, tables6, ref23]
  OATP1B3_vivo = 1.000 #pmol/mg membrane protein [@lin2023, tables7, ref23]
  REF_OATP1B3 = OATP1B3_vivo/OATP1B3_vitro
  
  Vmax_OATP2B1 = 1.493 * 1e-6 #umol/min/mg protein; 1.493± 0.543 pmol/min/mg protein [@lin2023]
  Km_OATP2B1 = 148.68 #148.68± 132.89 uM [@lin2023]
  OATP2B1_vitro = 0.812 #pmol/mg membrane protein
  OATP2B1_vivo = 1.600 #pmol/mg membrane protein [@lin2023, tables7, ref23]
  REF_OATP2B1 = OATP2B1_vivo/OATP2B1_vitro
  
  
  
  parms <- unlist(c(data.frame(BW, #1
                               QC, #2
                               VA, #3
                               VG, #4
                               VK, #5
                               VL, #6
                               VP, #7
                               VSk, #8
                               VR, #9
                               VFil, #10
                               QA, #11
                               QG, #12 
                               QK, #13
                               QL, #14
                               QSk, #15
                               QR, #16
                               QUr, #17
                               GFR, #18
                               PA, #19
                               PG, #20
                               PK, #21
                               PL, #22
                               PSk, #23
                               PR, #24
                               CL_FiltPT, #25
                               CLbiliary, #26
                               CLfecal, #27
                               fup, #28
                               CL_OAT4, #46 
                               REF_OAT4,
                               VPTL,
                               VRKL,
                               PK,
                               QK,
                               fup,
                               fuc = fup,
                               GFR = GFR,
                               CL_PltPT,
                               CL_FiltPT,
                               CLdif_PTLtPTC,
                               CLdif_PTCtPTP,
                               CLdif_PTCtPTL,
                               CLdif_PTPtPTC,
                               QT,
                               QUr,
                               fuT,
                               fuPTL,
                               fuRKL,
                               fuexp, 
                               Vmax_OAT4, 
                               Vmax_OAT3,
                               Vmax_OAT1,
                               Km_OAT4,
                               Km_OAT3,
                               Km_OAT1,
                               VPT,
                               VRK,
                               VGL,
                               VGL_up,
                               VGL_low,
                               VGL_col,
                               VL_ec,
                               VL_ic,
                               VmaxBSEP,
                               KmBSEP, 
                               CL_GLtG,
                               ka,
                               tsi_up,
                               tsi_low,
                               tco,
                               fup,
                               Vmax_OATP1B1,
                               Km_OATP1B1,
                               REF_OATP1B1,
                               Vmax_OATP1B3,
                               Km_OATP1B3,
                               REF_OATP1B3,
                               Vmax_OATP2B1,
                               Km_OATP2B1,
                               REF_OATP2B1
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
      
      CPT <- APT/VPT
      CPTL <- APTL/VPTL
      CRK <- ARK/VRK
      CRKL <- ARKL/VRKL
      
      CSk <- ASk/VSk  # Concentration in skin (ug/L)
      CVSk <- CSk/PSk  # Concentration leaving skin (ug/L) 
      
      # CG <- AG/VG  # Concentration in gut (ug/L)
      # CVG <- CG/PG # Concentration leaving gut (ug/L)
      # 
      # CL <- AL/VL  # Concentration in liver (ug/L)
      # CVL <- CL/PL  # Concentration leaving liver (ug/L)
      # 
      CGL <- AGL/VGL # Concentration in gut lumen (ug/L)
      CGL_up <- AGL_up/VGL_up
      CGL_low <- AGL_low/VGL_low
      CGL_col <- AGL_col/VGL_col
      
      CG <- AG/VG  # Concentration in gut (ug/L)
      CVG <- CG/PG # Concentration leaving gut (ug/L)
      
      
      CL <- AL/VL  # Concentration in liver (ug/L)
      CL_ec <- AL_ec/VL_ec
      CL_ic <- AL_ic/VL_ic
      
      CVG <- CG/PG # Concentration leaving gut (ug/L)
      CVL <- CL/PL  # Concentration leaving liver (ug/L)
      CVL_ec <- CL_ec/PL
      
      
      
      CA <- AA/VA  # Concentration in adipose (ug/L)
      CVA <- CA/PA # Concentration leaving adipose (ug/L)
      
      CR <- AR/VR  # Concentration in rest the body (ug/L)
      CVR <- CR/PR  # Concentration leaving rest of the body (ug/L)
      
      CP <- AP/VP  # Concentration in plasma (ug/L)
     
      
      ## Differential equations ----
      dD = - D
      
      # Proximal tubule tissue
      dAPT <- QK*(CP - CPT) + 
        + (Vmax_OAT4/(Km_OAT4+(CPTL*fuPTL)))*CPTL*fuPTL #+ 
        # + CLdif_PTLtPTC*CPTL*fuPTL +
        # - CLdif_PTCtPTL*CPT*fuT
      
      # Proximal tubule lumen
      dAPTL <- + fup*GFR*CP - QT*CPTL +
        - (Vmax_OAT4/(Km_OAT4+(CPTL*fuPTL)))*CPTL*fuPTL #+
        # - CLdif_PTLtPTC*CPTL*fuPTL +
        # + CLdif_PTCtPTL*CPT*fuT #QT*CPTL QUr*CPTL
      
      # Rest of kidney
      dARK <- QK*(CPT - CRK/PK) 
      
      # Rest of kidney lumen
      dARKL <- QT*CPTL - QUr*CRKL
      
      # Urine
      dAUr <- QUr*CRKL #CRKL
      
      
      
      # Skin compartment
      dASk <- QSk*(CP-CVSk) #+ Dermaldose 
      
      
      # Gut compartment: 
      dAGL <- + D - tco*AGL - CL_GLtG*CGL + (VmaxBSEP/(KmBSEP + (CL_ic*fuT)))*CL_ic*fuT # CL_GLtG*CGL, CLfecal*CGL, ka*AGL
      dAGL_up <- 0 #D - CL_GLtG*CGL_up + (VmaxBSEP/(KmBSEP + (CL_ic*fuT)))*CL_ic*fuT - tsi_up*AGL_up
      dAGL_low <- 0 #tsi_up*AGL_up - CL_GLtG*CGL_low - tsi_low*AGL_low
      dAGL_col <- 0 #tsi_low*AGL_low - CL_GLtG*CGL_col - tco*AGL_col
      
      dAG <- QG*(CP - CVG) + CL_GLtG*CGL #CL_GLtG*CGL_up + CL_GLtG*CGL_low + CL_GLtG*CGL_col #  CL_GLtG*CGL ka*AGL
      # dAG <- D + QG*CP - QG*CVG + CLbiliary*CL*fuT - CLfecal*CG # CLfecal*CG #+ Oraldose kco*AG
      
      # Excretion fecal: cumulative
      dAEx_feces <-  tco*AGL #_col #tco*AGL_col #tco*AGL # #CLfecal*CGL
      # dAEx_feces <-  CLfecal*CG #kco*AG #
      
      # Liver compartment
      # dAL <- QL*CP + QG*CVG - (QL+QG)*CVL - CLbiliary*CL*fuT 
      # dAL <-  0 #QG*CVG + QL*CP - (QG+QL)*CVL - (VmaxBSEP/(KmBSEP + (CL*fuT)))*CL*fuT ##CLbiliary*CL*fuT 
      dAL <- 0 #QG*CVG + QL*CP - (QG+QL)*CVL - (VmaxBSEP/(KmBSEP + (CL_ic*fuT)))*CL*fuT
      dAL_ec <- QG*CVG + QL*CP - (QG+QL)*CVL_ec +
        - (Vmax_OATP1B1/(Km_OATP1B1 + (CL_ec*fup)))*CL_ec*fup*REF_OATP1B1 +
        - (Vmax_OATP1B3/(Km_OATP1B3 + (CL_ec*fup)))*CL_ec*fup*REF_OATP1B3 #+
        # - (Vmax_OATP2B1/(Km_OATP2B1 + (CL_ec*fup)))*CL_ec*fup*REF_OATP1B1
      
      dAL_ic <- (Vmax_OATP1B1/(Km_OATP1B1 + (CL_ec*fup)))*CL_ec*fup*REF_OATP1B1 +
        + (Vmax_OATP1B3/(Km_OATP1B3 + (CL_ec*fup)))*CL_ec*fup*REF_OATP1B3 +
        # + (Vmax_OATP2B1/(Km_OATP2B1 + (CL_ec*fup)))*CL_ec*fup*REF_OATP1B1 #+
        - (VmaxBSEP/(KmBSEP + (CL_ic*fuT)))*CL_ic*fuT ##CLbiliary*CL*fuT 
      
      # Adipose compartment
      dAA <- QA*(CP-CVA) # (ug/d)
      
      # Rest compartment
      dAR <- QR*(CP-CVR) # (ug/d)
      
      # Plasma compartment
      dAP <- - (QSk + QG + QL + QA + QR)*CP +
        QSk*CVSk + (QL+QG)*CVL_ec + QA*CVA + QR*CVR + #(QL+QG)*CVL
        - QK*CP + (QK*CRK/PK) - fup*GFR*CP
      
      
      # Mass Balance
      Atot <- D + APT + APTL + ARK + ARKL + AUr +
        ASk +
        AG + AGL + AGL_up + AGL_low + AGL_col + AEx_feces + 
        AL + AL_ec + AL_ic +
        AA + AR +
        AP 
      
      # dInput <- Oraldose + Dermaldose
      # MB = dInput - Atot
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
             dAGL_up, 
             dAGL_low, 
             dAGL_col,
             dAG, 
             dAEx_feces,
             dAL,
             dAL_ec,
             dAL_ic,
             # dAG, 
             # dAEx_feces,
             # dAL, 
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
             CGL_up = CGL_up, 
             CGL_low = CGL_low, 
             CGL_col = CGL_col,
             CG = CG, CVG = CVG, 
             CL = CL, CVL = CVL,
             CL_ec = CL_ec,
             CL_ic = CL_ic,
             CVL_ec = CVL_ec,
             # CG = CG, CVG = CVG, 
             # CL = CL, CVL = CVL, 
             CA = CA, CVA = CVA,
             CR = CR, CVR = CVR,
             CP = CP, 
             Atot = Atot,MB = MB
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
              AGL_up = 0, 
              AGL_low = 0, 
              AGL_col = 0,
              AG = 0,
              AEx_feces = 0,
              AL = 0, 
              AL_ec = 0,
              AL_ic = 0,
              # AG = 0, 
              # AEx_feces = 0,
              # AL = 0, 
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
  
  
  # # PFOA in Kidney of one individual
  # Plot_PFOA_NoPTC <- output.PFOA.df %>%
  #   select(Days, CP, CPTL, CRKL, CPT, CRK) %>% # CGlomP, CGlomL,
  #   filter(Days > 0.05) %>%
  #   rename(Total_Plasma = CP,
  #          ProximalT_Lumen = CPTL,
  #          RestK_Lumen = CRKL,
  #          Proximal = CPT,
  #          Rest = CRK
  #   ) %>%
  #   pivot_longer(names_to = "Organ", values_to = "Concentration", Total_Plasma:Rest) %>%
  #   filter(Organ != "ProximalT_Cell") %>%
  #   ggplot()+
  #   geom_path(aes(x = Days, y = Concentration, color = Organ), linewidth = 0.5) +
  #   facet_wrap(~ Organ)+
  #   theme_CP() +
  #   theme(legend.position = "null")+
  #   #       panel.background = element_rect(fill = "white"),
  #   #       panel.border = element_blank()
  #   #       ) +
  #   ylab("Concentration (ng/ml)") +
  #   xlab("Time (d)")
  # Plot_PFOA_NoPTC
  # ggsave("Plot_PFOA_NoPTC.png", dpi = 300)
  # 
  # # PFOA proximal tubule cell
  # Plot_PFOA_PT <- output.PFOA.df %>%
  #   select(Days, CPT) %>%
  #   rename(Concentration = CPT) %>%
  #   mutate(Organ = "Proximal") %>%
  #   ggplot()+
  #   geom_path(aes(x = Days, y = Concentration), linewidth = 0.5) +
  #   theme_CP() +
  #   labs(title = "Proximal Tubule Tissue") +
  #   ylab("Concentration (ng/ml)") +
  #   xlab("Time (d)")
  # Plot_PFOA_PT
  ggsave("Plot_PFOA_PT.png", dpi = 300)
  
  # PFOA in all organs of one individual
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
  
  
  # Calculate AUC and Half life ####
  # ---------------------------------------------------------------------------- #
  
  AUC <- trapz(output_PFOA[ , "time"], output_PFOA[ , "CP"]) #ug*day/L
  
  # Calculate predicted half-life
  time <- output_PFOA[ , "time"] #in days
  conc <- output_PFOA[ , "CP"] #(ug/L) or (ng/ml)
  Cmax <- max(conc)
  Tmax <- time[which.max(conc)]
  tlast <- max(time[conc > 0])
  
  half_life <- pk.calc.half.life(
      conc,
      time,
      Tmax,
      tlast
    )
  
  HalfLife <- half_life$half.life/365 #halflife in years
  
  # Experimental
  ExpData <- read_csv("C:/Users/pacho003/OneDrive - Wageningen University & Research/CP_L_R/PARC_PFOA_mechanistic/Input/HalfLifes.csv")
  
  Exp_HalfLifes <- ExpData %>%
    filter(species == "human",
           chemical == "pfoa",
           parameter== "HalfLife") %>%
    select(value_average) %>%
    rename(HalfLife = value_average)
  Exp_HalfLifes$HalfLife <- as.numeric(Exp_HalfLifes$HalfLife) #years
  
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
  
  
  # Experimental Vs Simulated Plasma concentrations ####
  ExpPlasma <- read_excel("C:/Users/pacho003/OneDrive - Wageningen University & Research/CP_L_R/PARC_PFOA_mechanistic/Input/Experimental.Plasma.PFOA.xlsx", 
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
  
  
  
  
  print(AUC)
  print(HalfLife)