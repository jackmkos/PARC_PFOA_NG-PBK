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
  VRKP <- VK - VPTP - VPTC - VPTL - VRKL
  VRK <- VK - VPT - VPTL - VRKL
  
  VKL <- VPTL + VRKL 
  VKT <- VK - VPTL - VRKL
  
  
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
                               VPTP,
                               VPTC,
                               VPTL,
                               VRKP,
                               VRKC,
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
                               VKT,
                               VKL
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
      CL <- AL/VL  # Concentration in liver (ug/L)
      CA <- AA/VA  # Concentration in adipose (ug/L)
      CR <- AR/VR  # Concentration in rest the body (ug/L)
      
      CSk <- ASk/VSk  # Concentration in skin (ug/L)
      
      # CK <- AK/VK  # Concentration in kidney (ug/L)
      # # CKP <- AKP/VKP # Concentration in kidney plasma (ug/L) 
      # # CKT <- AKT/VKT # Concentration in kidney tissue (ug/l) 
      # CFil <- AFil/VFil # Concentration in filtrate compartment
      
      CPT <- APT/VPT
      CPTP <- APTP/VPTP
      CPTC <- APTC/VPTC
      CPTL <- APTL/VPTL
      CRK <- ARK/VRK
      CRKP <- ARKP/VRKP
      CRKL <- ARKL/VRKL
      CKT <- AKT/VKT
      CKL <- AKL/VKL
      
      # Organ venous concentrations (ug/L)
      CVG <- CG/PG # Concentration leaving gut (ug/L)
      CVL <- CL/PL  # Concentration leaving liver (ug/L)
      CVA <- CA/PA # Concentration leaving adipose (ug/L)
      CVSk <- CSk/PSk  # Concentration leaving skin (ug/L) 
      CVR <- CR/PR  # Concentration leaving rest of the body (ug/L)
      
      # CVK <- CK/PK # Concentration leaving the kidney (ug/L)
      # CVKP <- CKP/PK # Concentration leaving kidney tissue (ug/l)
      
      
      ## Differential equations ----
      
      # Adipose compartment
      dAA <- QA*(CP-CVA) # (ug/d)
      
      # Rest compartment
      dAR <- QR*(CP-CVR) # (ug/d)
      
      # Gut compartment: 
      dAG <- + QG*CP - QG*CVG + CLbiliary*CL*fuT - CLfecal*CG #+ Oraldose
      
      # Excretion fecal: cumulative
      dAEx_feces <- CLfecal*CG 
      
      # Liver compartment
      dAL <- QL*CP + QG*CVG - (QL+QG)*CVL - CLbiliary*CL*fuT 
      
      # Proximal tubule plasma
      dAPTP <- 0 #QK*(CP - CPTP) - CL_PltPT*CPTP*fup + CLdif_PTCtPTP*CPTC*fuT - CLdif_PTPtPTC*CPTP*fup
      #                         plasTcell           cellTplas                plasTcell
      
      # Proximal tubule cell
      dAPTC <- 0 #CL_PltPT*CPTP*fup + CL_FiltPT*CPTL*fuPTL + CLdif_PTLtPTC*CPTL*fuPTL - CLdif_PTCtPTP*CPTC*fuT - CLdif_PTCtPTL*CPTC*fuT + CLdif_PTPtPTC*CPTP*fup
      #        plasTcell           filTcell               filTcell                   cellTplas                cellTfil                 plasTcell
      
      # Proximal tubule tissue
      dAPT <- QK*(CP - CPT) + (Vmax_OAT4/(Km_OAT4+(CPTL*fuPTL)))*CPTL*fuPTL + CLdif_PTLtPTC*CPTL*fuPTL - CLdif_PTCtPTL*CPT*fuT
      
      #dAPT <- QK*(CP - CPT) + CL_FiltPT*CPTL*fuPTL + CLdif_PTLtPTC*CPTL*fuPTL - CLdif_PTCtPTL*CPT*fuT
      #                       filTcell               filTcell                   cellTfil
      
      # Proximal tubule lumen
      dAPTL <- + fup*GFR*CP - QT*CPTL - (Vmax_OAT4/(Km_OAT4+(CPTL*fuPTL)))*CPTL*fuPTL - CLdif_PTLtPTC*CPTL*fuPTL + CLdif_PTCtPTL*CPT*fuT #QT*CPTL QUr*CPTL
      #dAPTL <- + fup*GFR*CP - QT*CPTL - CL_FiltPT*CPTL*fuPTL - CLdif_PTLtPTC*CPTL*fuPTL + CLdif_PTCtPTL*CPT*fuT #QT*CPTL QUr*CPTL
      #          GFR            flow      filTcell               filTcell                   cellTfill
      
      # Rest of kidney
      dARK <- QK*(CPT - CRK/PK) 
      
      
      # Rest of kidney plasma (rest of kidney cell is ignored)
      dARKP <- 0 #QK*(CPTP - CRKP) 
      
      # Rest of kidney lumen
      dARKL <- QT*CPTL - QUr*CRKL
      
      # Kidney tissue
      dAKT <- 0 #QK*(CP - CKT/PK) + CL_FiltPT*CKL*fuPTL + CLdif_PTLtPTC*CKL*fuPTL - CLdif_PTCtPTL*CKT*fuT
      
      # Kidney Lumen
      dAKL <- 0 #+ fup*GFR*CP - QUr*CKL - CL_FiltPT*CKL*fuPTL - CLdif_PTLtPTC*CKL*fuPTL + CLdif_PTCtPTL*CKT*fuT #QT*CPTL QUr*CPTL
      
      # Urine
      dAUr <- QUr*CRKL #QUr*CPTL
      
      # Skin compartment
      dASk <- QSk*(CP-CVSk) #+ Dermaldose 
      
      # Plasma compartment
      dAP <- - (QSk + QG + QL + QA + QR)*CP +
        QSk*CVSk + (QL+QG)*CVL + QA*CVA + QR*CVR +
        - QK*CP + (QK*CRK/PK) - fup*GFR*CP
      
      
      # Mass Balance
      Atot <- AP + ASk + AG + AL + AA + AR + AEx_feces +
        + APTP + APTC + APTL + ARKP + ARKL + AUr +
        + APT + ARK + AKT + AKL
        
      # dInput <- Oraldose + Dermaldose
      # MB = dInput - Atot
      MB = Oraldose - Atot + 1
      
      # End
      
      list(c(dAP, 
             dASk, 
             dAG, 
             dAL, 
             dAA,
             # dAK,
             # # dAKT,
             # # dAKP,
             # dAFil,
             # dAUr, 
             dAR, 
             dAEx_feces,
             dAPTP,
             dAPTC,
             dAPTL,
             dARKP,
             dARKL,
             dAUr,
             dAPT,
             dARK,
             dAKT,
             dAKL
             #,
             # dInput
             ), 
           c(CP = CP, 
             CSk = CSk, 
             CG = CG, CVG = CVG, 
             CL = CL, CVL = CVL, 
             CA = CA, CVA = CVA,
             # CK = CK,
             # # CKP = CKP, CKT = CKT, 
             # CFil = CFil, 
             # CVK = CVK,
             # CVKP = CVKP, 
             CPTP = CPTP,
             CPTC = CPTC,
             CPTL = CPTL,
             CRKP = CRKP,
             CRKL = CRKL,
             CPT = CPT,
             CRK = CRK,
             CKT = CKT,
             CKL = CKL,
             CR = CR, CVR = CVR,
             Atot = Atot,MB = MB)
           )
    })
  }
  
  ## Initials ####
  
  A_init <- c(AP = 0, 
              # ASkb = 0, 
              ASk = 0, 
              AG = Oraldose, 
              AL = 0, 
              AA = 0, 
              # AK = 0,
              # # AKP = 0,
              # # AKT = 0,
              # AFil = 0,
              # AUr = 0, 
              AR = 0, 
              AEx_feces = 0,
              APTP = 0,
              APTC = 0,
              APTL = 0,
              ARKP = 0,
              ARKL = 0,
              AUr = 0,
              APT = 0,
              ARK = 0,
              AKT = 0,
              AKL = 0
              #Input = 0
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
  
  
  # PFOA in Kidney of one individual
  Plot_PFOA_NoPTC <- output.PFOA.df %>%
    select(Days, CP, CPTP, CPTC, CPTL, CRKP, CRKL, CPT, CRK) %>% # CGlomP, CGlomL,
    filter(Days > 0.05) %>%
    rename(Total_Plasma = CP,
           # Glomerular_Plasma = CGlomP,
           # Glomerular_Lumen = CGlomL,
           ProximalT_Plasma = CPTP,
           ProximalT_Cell = CPTC,
           ProximalT_Lumen = CPTL,
           RestK_Plasma = CRKP,
           RestK_Lumen = CRKL,
           Proximal = CPT,
           Rest = CRK
    ) %>%
    pivot_longer(names_to = "Organ", values_to = "Concentration", Total_Plasma:Rest) %>%
    filter(Organ != "ProximalT_Cell") %>%
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
  Plot_PFOA_NoPTC
  ggsave("Plot_PFOA_NoPTC.png", dpi = 300)
  
  # PFOA proximal tubule cell
  # Plot_PFOA_PTCell <- output.PFOA.df %>%
  #   select(Days, CPTC) %>%
  #   rename(Concentration = CPTC) %>%
  #   mutate(Organ = "ProximalT_Cell") %>% 
  #   ggplot()+
  #   geom_path(aes(x = Days, y = Concentration), linewidth = 0.5) +
  #   theme_CP() +
  #   labs(title = "Proximal Tubule Cell") +
  #   ylab("Concentration (ng/ml)") +
  #   xlab("Time (d)")
  # Plot_PFOA_PTCell
  # ggsave("Plot_PFOA_PTCell.png", dpi = 300)
  
  # PFOA in all organs of one individual
  Plot_PFOA_All <- output.PFOA.df %>% 
    mutate(CK = CPTP+CPTC+CPTL+CRKP+CRKL+CPT+CRK) %>% 
    select(Days, CK, CSk, CL, CG, CA, CR, CP) %>% 
    rename(Kidney = CK, Skin = CSk, Liver = CL, Gut = CG, Adipose = CA, Rest = CR, Plasma = CP) %>% 
    pivot_longer(names_to = "Organ", values_to = "Concentration", Kidney:Plasma) %>% 
    ggplot()+
    geom_path(aes(x = Days, y = Concentration, color = Organ)) +
    facet_wrap(~ Organ)+
    theme_CP()+
    theme(legend.position = "none")+
    # theme(axis.text.x = element_text(size = 7),
    #       axis.text.y = element_text(size = 7), 
    #       axis.title = element_text(size = 8)) +
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
  ExpData <- read_csv("C:/Users/pacho003/OneDrive - Wageningen University & Research/CP_L_R/PARC_PFAS_PBPKmodel/Codes_PFAS_models/Chrysanthi_Pachoulide_2024_PFOA_mechanistic/Input/HalfLifes.csv")
  
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
  
  
  
  
  print(AUC)
  print(HalfLife)