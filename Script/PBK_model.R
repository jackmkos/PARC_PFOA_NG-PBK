# --------------------------------------------------------------------------- #
# SCRIPT FOR PBK MODEL FOR ORAL EXPOSURE
# By: Chrysanthi Pachoulide
# Date: 14-04-2025
# --------------------------------------------------------------------------- #

# ORAL ####
ORAL_PBK_RUN <- function(y, parms, times){ # Input for ode
  
  PBK.model <- function(t, state, parameters){
    with(as.list(c(state, parameters)), {
      
      ### Physiological ----
      
      VIL <- VILc * BW                # L, Volume of intestinal lumen
      VI <- VIc * BW                  # L, Volume of intestine
      SA_SI <- SA_SIc * BW
      
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
      
      VP <- VPc * (1-Hct) * BW                  # L, Volume of plasma
      
      VTotc <- 0.96
      VTot <- VTotc * BW              # L, Total body volume (used for mass balance)
      
      VR <- VTot - (VI + VL + VK + VA + VP)             # L, Volume of the lumped rest compartment
      
      QTotc <- 0.988
      QI <- QIc/QTotc * QC                  # L/d, Intestinal
      QL <- QLc/QTotc * QC                  # L/d, Liver
      QK <- QKc/QTotc * QC                  # L/d, Kidney
      QA <- QAc/QTotc * QC                  # L/d, Adipose 
      QR <- QC - (QA + QI + QK + QL)  # L/d, Rest
      
      QUr <- QUrc * BW              # L/d, Urine flow rate to the bladder 22 mL/kg BW/d [ICRP 89 page 161]
      GFR <- GFRc * QK              # L/d 18% of total renal plasma flow [ICRP 89 page 159] http://www.icrp.org/publication.asp?id=ICRP%20Publication%2089
      # GFR <- GFR
      QT <- QT                      # L/d, Proximal tubule fluid flow
      
      tco <- tco                    # /d, Bowel residence times in the colon
      
      
      ### Physicochemical ----
      MW <- 414.07
      
      PI <- PIc * fup     # Intestinal
      PL <- PLc * fup     # Liver
      PK <- PKc * fup     # Kidney
      PA <- PAc * fup     # Adipose
      PR <- PRc * fup     # Rest
      
      # Fraction unionised
      pH_P <- 7.4     # plasma
      pH_IL <- 7      # intestinal Intestinal lumen, average
      
      f.union_p <- 1/(1 + 10^(pH_P - pKa))       # Plasma
      f.union_exp <- 1/(1 + 10^(pH_P - pKa))     # Is the same as plasma as pH in the experiment is 7.4
      f.union_IL <- 1/(1 + 10^(pH_IL - pKa))     # Intestinal lumen
      
      # Fraction unbound, calculated based on Poulin and Haddad, 2018 https://doi.org/10.1016/j.xphs.2018.03.012
      # Equation was adapted to not account for fraction unionised
      # OAT and OATP transporters transport the ionised compound, given that the ratio of fraction ionised at plasma to cellular pH is 1, this can be ignored (fraction unionised of PFOA is 0.9999923 at pH 7.4 and 0.9999963 at pH 7)
      fu_PTL <- R_PTL*fup/(1 + ((R_PTL-1)*fup)) # Proximal tubule lumen
      fu_Lic <- R_L_ec*fup/(1 + ((R_L_ec-1)*fup)) # Liver intracellular space
      
      
      ### Kinetic ----
      
      # Gastro-intestinal uptake
      Pint_SI <- Papp_SI/f.union_exp                       # cm/s, Intrinsic permeability, corrected for fraction unionised in the experiment
      CL_IL <- (Pint_SI*SA_SI*f.union_IL*1e-3)*60*60*24    # L/d, Intestinal lumen to intestinal tissue (calculations: cm/s -> L/s /1000 -> L/d *60*60*24) 
      
      # Liver uptake
      Vmax_OATP1B1 <- Vmax_OATP1B1c*MW*60*24*SF_OATP1B1*VL_ec             # ug/d
      Km_OATP1B1 <- Km_OATP1B1c*MW                                        # ug/L (uM -> ug/L)
      Vmax_OATP1B3 <- Vmax_OATP1B3c*MW*60*24*SF_OATP1B3*VL_ec             # ug/d
      Km_OATP1B3 <- Km_OATP1B3c*MW                                        # ug/L (uM -> ug/L)
      
      # Biliary excretion
      VmaxBSEP <- VmaxBSEPc*MW*60*24*SF_BSEP*VL_ic         # ug/d
      KmBSEP <- KmBSEPc*MW                                 # ug/L (uM -> ug/L)
      
      # Renal clearance
      Vmax_OAT4 = Vmax_OAT4c*MW*60*24*SF_OAT*VPT           # ug/d (umol -> ug, min -> d)
      Km_OAT4 = Km_OAT4c*MW                                # ug/L (uM -> ug/L)
      
      
      ## Dose -------------------------
      
      if(t<expSTOP){DoseOn=1} else{DoseOn=0}
      
      ## Oral exposure ##
      DOral = expOral*BW*DoseOn         # ug, PFOA oral dose
      OralD = DOral/Tinput*(t %% tinterval<Tinput)
      
      ## Concentrations -------------------------
      
      CIL <- AIL/VIL               # ug/L, Intestinal lumen
      CI <- AI/VI                  # ug/L, Intestine 
      CVI <- CI/PI                 # ug/L, Intestine venous
      
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
      
      CP <- AP/VP                  # ug/L, Plasma
      
      
      ## Differential equations -------------------------
      
      dOD = OralD - OD             # ug/d, Oral dose input 
      
      dAIL <- + OD - tco*AIL - CL_IL*CIL + 
        + (VmaxBSEP/(KmBSEP + (CL_ic*fu_Lic)))*CL_ic*fu_Lic          # ug/d, Intestine lumen
      
      dAI <- QI*(CP - CVI) + CL_IL*CIL                      # ug/d, Intestinal
      
      dAFe <-  tco*AIL                                       # ug/d, Feces
      
      
      dAL_ec <- + QI*CVI + QL*CP - (QI+QL)*CVL_ec + 
        - (Vmax_OATP1B1/(Km_OATP1B1 + (CL_ec*fup)))*CL_ec*fup +
        - (Vmax_OATP1B3/(Km_OATP1B3 + (CL_ec*fup)))*CL_ec*fup            # ug/d, Liver extracellular space (vascular + interstitial space)
      
      dAL_ic <- (Vmax_OATP1B1/(Km_OATP1B1 + (CL_ec*fup)))*CL_ec*fup +
        + (Vmax_OATP1B3/(Km_OATP1B3 + (CL_ec*fup)))*CL_ec*fup +
        - (VmaxBSEP/(KmBSEP + (CL_ic*fu_Lic)))*CL_ic*fu_Lic             # ug/d, Liver intracellular space
      
      
      dAPTT <- QK*(CP - CPTT) + 
        + (Vmax_OAT4/(Km_OAT4+(CPTL*fu_PTL)))*CPTL*fu_PTL      # ug/d, Proximal tubule tissue 
      
      dAPTL <- + fup*GFR*CP - QT*CPTL +
        - (Vmax_OAT4/(Km_OAT4+(CPTL*fu_PTL)))*CPTL*fu_PTL      # ug/d, Proximal tubule lumen    
      
      dARKT <- QK*(CPTT - CVRKT)                               # ug/d, Rest of kidney
      
      dARKL <- QT*CPTL - QUr*CRKL                              # ug/d, Rest of kidney lumen
      
      dAUr <- QUr*CRKL                                         # ug/d, Urine
      
      
      dAA <- QA*(CP-CVA)                                      # ug/d, Adipose
      
      
      dAR <- QR*(CP-CVR)                                      # ug/d, Rest
      
      dAP <- - (QI + QL + QA + QR + QK)*CP - fup*GFR*CP +     # ug/d, Arterial Plasma
        + (QL+QI)*CVL_ec + QK*CVRKT + QA*CVA + QR*CVR    # ug/d, Venous Plasma
      
      # Mass Balance
      Atot <- OD +
        AIL + AI + AFe + 
        AL_ec + AL_ic +
        APTT + APTL + ARKT + ARKL + AUr +
        AA + 
        AR +
        AP
      
      dAin <- OralD # to be used if repeated exposure
      MB <- Ain - Atot + 1    # to be used if repeated exposure
      # MB <- DOral - Atot + 1
      
      # End
      
      list(c(dOD,
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
             dAP, 
             dAin
      ), 
      c(CIL = CIL,
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
        CP = CP,
        Atot = Atot, 
        MB = MB
      )
      )
    })
  }
  
  output_PBK <-lsoda(y = y, 
                      times = times, 
                      func = PBK.model,
                      parms = parms)
  
  return(as.data.frame(output_PBK))
}
  
# DERMAL ####
DERMAL_PBK_RUN <- function(y, parms, times){ # Input for ode
  
  
  PBK.model <- function(t, state, parameters){
    with(as.list(c(state, parameters)), {
      
      ### Physiological ----
      
      VSkB <- VSkB                    # L, Volume of skin barrier
      VSk <- VSkc * BW                # L, Volume of skin
      
      VIL <- VILc * BW                # L, Volume of intestinal lumen
      VI <- VIc * BW                  # L, Volume of intestine
      SA_SI <- SA_SIc * BW
      
      
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
      
      VP <- VPc * (1-Hct) * BW                  # L, Volume of plasma
      
      VTotc <- 0.96
      VTot <- VTotc * BW              # L, Total body volume (used for mass balance)
      
      VR <- VTot - (VSk + VI + VL + VK + VA + VP)             # L, Volume of the lumped rest compartment
      
      QTotc <- 0.988
      QSk <- QSkc/QTotc * QC                # L/d, Skin 
      QI <- QIc/QTotc * QC                  # L/d, Intestinal
      QL <- QLc/QTotc * QC                  # L/d, Liver
      QK <- QKc/QTotc * QC                  # L/d, Kidney
      QA <- QAc/QTotc * QC                  # L/d, Adipose 
      QR <- QC - (QA + QI + QK + QL + QSk)  # L/d, Rest
      
      QUr <- QUrc * BW              # L/d, Urine flow rate to the bladder 22 mL/kg BW/d [ICRP 89 page 161]
      # GFR <- GFRc * QK              # L/d 18% of total renal plasma flow [ICRP 89 page 159] http://www.icrp.org/publication.asp?id=ICRP%20Publication%2089
      GFR <- GFR
      QT <- QT                      # L/d, Proximal tubule fluid flow
      
      tco <- tco                    # /d, Bowel residence times in the colon
      
      
      ### Physicochemical ----
      MW <- 414.07
      
      PSk <- PSkc * fup  # Skin
      PI <- PIc * fup    # Intestinal
      PL <- PLc * fup    # Liver
      PK <- PKc * fup    # Kidney
      PA <- PAc * fup    # Adipose
      PR <- PRc * fup    # Rest
      
      # Fraction unionised
      pH_P <- 7.4     # plasma
      pH_IL <- 7      # intestinal Intestinal lumen, average
      
      f.union_p <- 1/(1 + 10^(pH_P - pKa))       # Plasma
      f.union_exp <- 1/(1 + 10^(pH_P - pKa))     # Is the same as plasma as pH in the experiment is 7.4
      f.union_IL <- 1/(1 + 10^(pH_IL - pKa))      # Intestinal lumen
      
      # Fraction unbound, calculated based on Poulin and Haddad, 2018 https://doi.org/10.1016/j.xphs.2018.03.012
      # Equation was adapted to not account for fraction unionised
      # OAT and OATP transporters transport the ionised compound, given that the ratio of fraction ionised at plasma to cellular pH is 1, this can be ignored (fraction unionised of PFOA is 0.9999923 at pH 7.4 and 0.9999963 at pH 7)
      fu_PTL <- R_PTL*fup/(1 + ((R_PTL-1)*fup)) # Proximal tubule lumen
      fu_Lic <- R_L_ec*fup/(1 + ((R_L_ec-1)*fup)) # Liver intracellular space
      
      
      
      ### Kinetic ----
      
      # Skin uptake
      CL_SkBtSk <- Papp_SkB*SA_SkB*1e-3*60*60*24           # L/d, Skin barrier to skin (calculations: cm/s -> L/s * 1e-3 -> L/d *60*60*24)
      
      # Gastro-intestinal uptake
      Pint_SI <- Papp_SI/f.union_exp                       # cm/s, Intrinsic permeability, corrected for fraction unionised in the experiment
      CL_IL <- (Pint_SI*SA_SI*f.union_IL*1e-3)*60*60*24    # L/d, Intestinal lumen to intestinal tissue (calculations: cm/s -> L/s /1000 -> L/d *60*60*24) 
      
      # Liver uptake
      Vmax_OATP1B1 <- Vmax_OATP1B1c*MW*60*24*SF_OATP1B1*VL_ec             # ug/d
      Km_OATP1B1 <- Km_OATP1B1c*MW                                        # ug/L (uM -> ug/L)
      Vmax_OATP1B3 <- Vmax_OATP1B3c*MW*60*24*SF_OATP1B3*VL_ec             # ug/d
      Km_OATP1B3 <- Km_OATP1B3c*MW                                        # ug/L (uM -> ug/L)
      
      # Biliary excretion
      VmaxBSEP <- VmaxBSEPc*MW*60*24*SF_BSEP*VL_ic         # ug/d
      KmBSEP <- KmBSEPc*MW                                 # ug/L (uM -> ug/L)
      
      # Renal clearance
      Vmax_OAT4 = Vmax_OAT4c*MW*60*24*SF_OAT*VPT           # ug/d (umol -> ug, min -> d)
      Km_OAT4 = Km_OAT4c*MW                                # ug/L, scaled from uM, Louisse et al. 2024 doi.org/10.1016/j.tox.2024.153961
      
      
      ## Dose -------------------------
      
      if(t<expSTOP){DoseOn=1} else{DoseOn=0}
      
      ## Dermal exposure ##
      DDermal = expDermal*BW*DoseOn         # ug, PFOA dermal dose
      DermalD = DDermal/Tinput*(t %% tinterval<Tinput)
      
      ## Concentrations -------------------------
      
      CSkB <- ASkB/VSkB            # ug/L, Skin barrier
      CSk <- ASk/VSk               # ug/L, Skin
      CVSk <- CSk/PSk              # ug/L, Skin venous 
      
      CIL <- AIL/VIL               # ug/L, Intestinal lumen
      CI <- AI/VI                  # ug/L, Intestine 
      CVI <- CI/PI                 # ug/L, Intestine venous
      
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
      
      CP <- AP/VP                  # ug/L, Plasma
      
      
      ## Differential equations -------------------------
      
      dDD = DermalD - DD                             # ug/d, Dermal dose input
      
      dASkB <- DD - CL_SkBtSk*CSkB                                # ug/d, Skin barrier
      
      dASk <- + CL_SkBtSk*CSkB + QSk*(CP-CVSk)                   # ug/d, Skin
      
      
      dAIL <- - tco*AIL - CL_IL*CIL + 
        + (VmaxBSEP/(KmBSEP + (CL_ic*fu_Lic)))*CL_ic*fu_Lic        # ug/d, Intestinal lumen
      
      dAI <- QI*(CP - CVI) + CL_IL*CIL                      # ug/d, Intestine
      
      dAFe <-  tco*AIL                                       # ug/d, Feces
      
      
      dAL_ec <- + QI*CVI + QL*CP - (QI+QL)*CVL_ec + 
        - (Vmax_OATP1B1/(Km_OATP1B1 + (CL_ec*fup)))*CL_ec*fup +
        - (Vmax_OATP1B3/(Km_OATP1B3 + (CL_ec*fup)))*CL_ec*fup            # ug/d, Liver extracellular space (vascular + interstitial space)
      
      dAL_ic <- (Vmax_OATP1B1/(Km_OATP1B1 + (CL_ec*fup)))*CL_ec*fup +
        + (Vmax_OATP1B3/(Km_OATP1B3 + (CL_ec*fup)))*CL_ec*fup +
        - (VmaxBSEP/(KmBSEP + (CL_ic*fu_Lic)))*CL_ic*fu_Lic                       # ug/d, Liver intracellular space
      
      
      dAPTT <- QK*(CP - CPTT) + 
        + (Vmax_OAT4/(Km_OAT4+(CPTL*fu_PTL)))*CPTL*fu_PTL        # ug/d, Proximal tubule tissue 
      
      dAPTL <- + fup*GFR*CP - QT*CPTL +
        - (Vmax_OAT4/(Km_OAT4+(CPTL*fu_PTL)))*CPTL*fu_PTL        # ug/d, Proximal tubule lumen    
      
      dARKT <- QK*(CPTT - CVRKT)                               # ug/d, Rest of kidney
      
      dARKL <- QT*CPTL - QUr*CRKL                              # ug/d, Rest of kidney lumen
      
      dAUr <- QUr*CRKL                                         # ug/d, Urine
      
      
      dAA <- QA*(CP-CVA)                                      # ug/d, Adipose
      
      
      dAR <- QR*(CP-CVR)                                      # ug/d, Rest
      
      dAP <- - (QSk + QI + QL + QA + QR + QK)*CP - fup*GFR*CP +     # ug/d, Arterial Plasma
        + (QL+QI)*CVL_ec + QSk*CVSk + QK*CVRKT + QA*CVA + QR*CVR    # ug/d, Venous Plasma
      
      # Mass Balance
      Atot <- DD + 
        ASk +
        AIL + AI + AFe + 
        AL_ec + AL_ic +
        APTT + APTL + ARKT + ARKL + AUr +
        AA + 
        AR +
        AP  
      
      dAin <- DermalD # to be used if repeated exposure
      MB <- Ain - Atot + 1    # to be used if repeated exposure
      # MB <- DOral - Atot + 1
      
      # End
      
      list(c(dDD,
             dASkB,
             dASk, 
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
             dAP, 
             dAin
      ), 
      c(CSkB = CSkB, 
        CSk = CSk, 
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
        CP = CP,
        Atot = Atot, 
        MB = MB
      )
      )
    })
  }
  
  output_PBK <-lsoda(y = y, 
                      times = times, 
                      func = PBK.model,
                      parms = parms)
  
  return(as.data.frame(output_PBK))
}


# ORAL & DERMAL ####
ORAL_DERMAL_PBK_RUN <- function(y, parms, times){ # Input for ode
  
  
  PBK.model <- function(t, state, parameters){
    with(as.list(c(state, parameters)), {
      
      ### Physiological ----
      
      VSkB <- VSkB                    # L, Volume of skin barrier
      VSk <- VSkc * BW                # L, Volume of skin
      
      VIL <- VILc * BW                # L, Volume of intestinal lumen
      VI <- VIc * BW                  # L, Volume of intestine
      SA_SI <- SA_SIc * BW
      
      
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
      
      VP <- VPc * (1-Hct) * BW                  # L, Volume of plasma
      
      VTotc <- 0.96
      VTot <- VTotc * BW              # L, Total body volume (used for mass balance)
      
      VR <- VTot - (VSk + VI + VL + VK + VA + VP)             # L, Volume of the lumped rest compartment
      
      QTotc <- 0.988
      QSk <- QSkc/QTotc * QC                # L/d, Skin 
      QI <- QIc/QTotc * QC                  # L/d, Intestinal
      QL <- QLc/QTotc * QC                  # L/d, Liver
      QK <- QKc/QTotc * QC                  # L/d, Kidney
      QA <- QAc/QTotc * QC                  # L/d, Adipose 
      QR <- QC - (QA + QI + QK + QL + QSk)  # L/d, Rest
      
      QUr <- QUrc * BW              # L/d, Urine flow rate to the bladder 22 mL/kg BW/d [ICRP 89 page 161]
      GFR <- GFRc * QK              # L/d 18% of total renal plasma flow [ICRP 89 page 159] http://www.icrp.org/publication.asp?id=ICRP%20Publication%2089
      # GFR <- GFR
      QT <- QT                      # L/d, Proximal tubule fluid flow
      
      tco <- tco                    # /d, Bowel residence times in the colon
      
      
      ### Physicochemical ----
      MW <- 414.07
      
      PSk <- PSkc * fup  # Skin
      PI <- PIc * fup    # Intestinal
      PL <- PLc * fup    # Liver
      PK <- PKc * fup    # Kidney
      PA <- PAc * fup    # Adipose
      PR <- PRc * fup    # Rest
      
      # Fraction unionised
      pH_P <- 7.4     # plasma
      pH_IL <- 7      # intestinal Intestinal lumen, average
      
      f.union_p <- 1/(1 + 10^(pH_P - pKa))       # Plasma
      f.union_exp <- 1/(1 + 10^(pH_P - pKa))     # Is the same as plasma as pH in the experiment is 7.4
      f.union_IL <- 1/(1 + 10^(pH_IL - pKa))      # Intestinal lumen
      
      # Fraction unbound, calculated based on Poulin and Haddad, 2018 https://doi.org/10.1016/j.xphs.2018.03.012
      # Equation was adapted to not account for fraction unionised
      # OAT and OATP transporters transport the ionised compound, given that the ratio of fraction ionised at plasma to cellular pH is 1, this can be ignored (fraction unionised of PFOA is 0.9999923 at pH 7.4 and 0.9999963 at pH 7)
      fu_PTL <- R_PTL*fup/(1 + ((R_PTL-1)*fup)) # Proximal tubule lumen
      fu_Lic <- R_L_ec*fup/(1 + ((R_L_ec-1)*fup)) # Liver intracellular space
      
      
      
      ### Kinetic ----
      
      # Skin uptake
      CL_SkBtSk <- Papp_SkB*SA_SkB*1e-3*60*60*24           # L/d, Skin barrier to skin (calculations: cm/s -> L/s * 1e-3 -> L/d *60*60*24)
      
      # Gastro-intestinal uptake
      Pint_SI <- Papp_SI/f.union_exp                       # cm/s, Intrinsic permeability, corrected for fraction unionised in the experiment
      CL_IL <- (Pint_SI*SA_SI*f.union_IL*1e-3)*60*60*24    # L/d, Intestinal lumen to intestinal tissue (calculations: cm/s -> L/s /1000 -> L/d *60*60*24) 
      
      # Liver uptake
      Vmax_OATP1B1 <- Vmax_OATP1B1c*MW*60*24*SF_OATP1B1*VL_ec             # ug/d
      Km_OATP1B1 <- Km_OATP1B1c*MW                                        # ug/L (uM -> ug/L)
      Vmax_OATP1B3 <- Vmax_OATP1B3c*MW*60*24*SF_OATP1B3*VL_ec             # ug/d
      Km_OATP1B3 <- Km_OATP1B3c*MW                                        # ug/L (uM -> ug/L)
      
      # Biliary excretion
      VmaxBSEP <- VmaxBSEPc*MW*60*24*SF_BSEP*VL_ic         # ug/d
      KmBSEP <- KmBSEPc*MW                                 # ug/L (uM -> ug/L)
      
      # Renal clearance
      Vmax_OAT4 = Vmax_OAT4c*MW*60*24*SF_OAT*VPT           # ug/d (umol -> ug, min -> d)
      Km_OAT4 = Km_OAT4c*MW                                # ug/L, scaled from uM, Louisse et al. 2024 doi.org/10.1016/j.tox.2024.153961
      
      
      ## Dose -------------------------
      
      if(t<expSTOP){DoseOn=1} else{DoseOn=0}
      
      ## Oral exposure ##
      DOral = expOral*BW*DoseOn         # ug, PFOA oral dose
      OralD = DOral/Tinput*(t %% tinterval<Tinput)
      
      ## Dermal exposure ##
      DDermal = expDermal*BW*DoseOn         # ug, PFOA dermal dose
      DermalD = DDermal/Tinput*(t %% tinterval<Tinput)
      
      ## Concentrations -------------------------
      
      CSkB <- ASkB/VSkB            # ug/L, Skin barrier
      CSk <- ASk/VSk               # ug/L, Skin
      CVSk <- CSk/PSk              # ug/L, Skin venous 
      
      CIL <- AIL/VIL               # ug/L, Intestinal lumen
      CI <- AI/VI                  # ug/L, Intestine 
      CVI <- CI/PI                 # ug/L, Intestine venous
      
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
      
      CP <- AP/VP                  # ug/L, Plasma
      
      
      ## Differential equations -------------------------
      
      dOD = OralD - OD             # ug/d, Oral dose input 
      dDD = DermalD - DD                             # ug/d, Dermal dose input
      
      
      dASkB <- DD - CL_SkBtSk*CSkB                                # ug/d, Skin barrier
      
      dASk <- + CL_SkBtSk*CSkB + QSk*(CP-CVSk)                   # ug/d, Skin
      
      
      dAIL <- + OD - tco*AIL - CL_IL*CIL + 
        + (VmaxBSEP/(KmBSEP + (CL_ic*fu_Lic)))*CL_ic*fu_Lic          # ug/d, Intestine lumen
      
      dAI <- QI*(CP - CVI) + CL_IL*CIL                      # ug/d, Intestine
      
      dAFe <-  tco*AIL                                       # ug/d, Feces
      
      
      dAL_ec <- + QI*CVI + QL*CP - (QI+QL)*CVL_ec + 
        - (Vmax_OATP1B1/(Km_OATP1B1 + (CL_ec*fup)))*CL_ec*fup +
        - (Vmax_OATP1B3/(Km_OATP1B3 + (CL_ec*fup)))*CL_ec*fup            # ug/d, Liver extracellular space (vascular + interstitial space)
      
      dAL_ic <- (Vmax_OATP1B1/(Km_OATP1B1 + (CL_ec*fup)))*CL_ec*fup +
        + (Vmax_OATP1B3/(Km_OATP1B3 + (CL_ec*fup)))*CL_ec*fup +
        - (VmaxBSEP/(KmBSEP + (CL_ic*fu_Lic)))*CL_ic*fu_Lic                       # ug/d, Liver intracellular space
      
      
      dAPTT <- QK*(CP - CPTT) + 
        + (Vmax_OAT4/(Km_OAT4+(CPTL*fu_PTL)))*CPTL*fu_PTL        # ug/d, Proximal tubule tissue 
      
      dAPTL <- + fup*GFR*CP - QT*CPTL +
        - (Vmax_OAT4/(Km_OAT4+(CPTL*fu_PTL)))*CPTL*fu_PTL        # ug/d, Proximal tubule lumen    
      
      dARKT <- QK*(CPTT - CVRKT)                               # ug/d, Rest of kidney
      
      dARKL <- QT*CPTL - QUr*CRKL                              # ug/d, Rest of kidney lumen
      
      dAUr <- QUr*CRKL                                         # ug/d, Urine
      
      
      dAA <- QA*(CP-CVA)                                      # ug/d, Adipose
      
      
      dAR <- QR*(CP-CVR)                                      # ug/d, Rest
      
      dAP <- - (QSk + QI + QL + QA + QR + QK)*CP - fup*GFR*CP +     # ug/d, Arterial Plasma
        + (QL+QI)*CVL_ec + QSk*CVSk + QK*CVRKT + QA*CVA + QR*CVR    # ug/d, Venous Plasma
      
      # Mass Balance
      Atot <- OD + DD + 
        ASk +
        AIL + AI + AFe + 
        AL_ec + AL_ic +
        APTT + APTL + ARKT + ARKL + AUr +
        AA + 
        AR +
        AP  
      
      dAin <- OralD + DermalD # to be used if repeated exposure
      MB <- Ain - Atot + 1    # to be used if repeated exposure
      # End
      
      list(c(dOD,
             dDD,
             dASkB,
             dASk, 
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
             dAP, 
             dAin
      ), 
      c(CSkB = CSkB, 
        CSk = CSk, 
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
        CP = CP,
        Atot = Atot, 
        MB = MB
      )
      )
    })
  }
  
  output_PBK <-lsoda(y = y, 
                      times = times, 
                      func = PBK.model,
                      parms = parms)
  
  return(as.data.frame(output_PBK))
}


# INHALATION ####
INHALATION_PBK_RUN <- function(y, parms, times){ # Input for ode
  
  
  PBK.model <- function(t, state, parameters){
    with(as.list(c(state, parameters)), {
      
      ### Physiological ----
      
      VIL <- VILc * BW                # L, Volume of intestinal lumen
      VI <- VIc * BW                  # L, Volume of intestine
      SA_SI <- SA_SIc * BW
      
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
      
      Fr_art_blood = 0.0178 / (0.0178 + 0.0533) #fraction of arterial blood (corrected for plasma), ref. A. Ratier et al. 2024 doi: 10.1016/j.envint.2024.108621
      VP <- VPc * (1-Hct) * BW                  # L, Volume of plasma
      VAP <- VP*Fr_art_blood         # L, Volume of arterial plasma
      VVP <- VP - VAP                 # L, Volume of venous plasma
      
      VTotc <- 0.96
      VTot <- VTotc * BW              # L, Total body volume (used for mass balance)
      
      VR <- VTot - (VI + VL + VK + VA + VLu + VAP + VVP)             # L, Volume of the lumped rest compartment
      
      QTotc <- 0.988
      QI <- QIc/QTotc * QC                  # L/d, Intestinal
      QL <- QLc/QTotc * QC                  # L/d, Liver
      QK <- QKc/QTotc * QC                  # L/d, Kidney
      QA <- QAc/QTotc * QC                  # L/d, Adipose 
      QLu <- QC                             # L/d, Lungs
      QR <- QC - (QA + QI + QK + QL)        # L/d, Rest
      
      QUr <- QUrc * BW              # L/d, Urine flow rate to the bladder 22 mL/kg BW/d [ICRP 89 page 161]
      # GFR <- GFRc * QK              # L/d 18% of total renal plasma flow [ICRP 89 page 159] http://www.icrp.org/publication.asp?id=ICRP%20Publication%2089
      GFR <- GFR
      QT <- QT                      # L/d, Proximal tubule fluid flow
      
      tco <- tco                    # /d, Bowel residence times in the colon
      
      
      ### Physicochemical ----
      MW <- 414.07
      
      PI <- PIc * fup     # Intestinal
      PL <- PLc * fup     # Liver
      PK <- PKc * fup     # Kidney
      PA <- PAc * fup     # Adipose
      PLu <- PLuc * fup   # Lungs
      PR <- PRc * fup     # Rest
      
      # Fraction unionised
      pH_P <- 7.4     # plasma
      pH_IL <- 7      # intestinal Intestinal lumen, average
      
      f.union_p <- 1/(1 + 10^(pH_P - pKa))       # Plasma
      f.union_exp <- 1/(1 + 10^(pH_P - pKa))     # Is the same as plasma as pH in the experiment is 7.4
      f.union_IL <- 1/(1 + 10^(pH_IL - pKa))     # Intestinal lumen
      
      # Fraction unbound, calculated based on Poulin and Haddad, 2018 https://doi.org/10.1016/j.xphs.2018.03.012
      # Equation was adapted to not account for fraction unionised
      # OAT and OATP transporters transport the ionised compound, given that the ratio of fraction ionised at plasma to cellular pH is 1, this can be ignored (fraction unionised of PFOA is 0.9999923 at pH 7.4 and 0.9999963 at pH 7)
      fu_PTL <- R_PTL*fup/(1 + ((R_PTL-1)*fup)) # Proximal tubule lumen
      fu_Lic <- R_L_ec*fup/(1 + ((R_L_ec-1)*fup)) # Liver intracellular space
      
      
      ### Kinetic ----
      
      # Gastro-intestinal uptake
      Pint_SI <- Papp_SI/f.union_exp                       # cm/s, Intrinsic permeability, corrected for fraction unionised in the experiment
      CL_IL <- (Pint_SI*SA_SI*f.union_IL*1e-3)*60*60*24    # L/d, Intestinal lumen to intestinal tissue (calculations: cm/s -> L/s /1000 -> L/d *60*60*24) 
      
      # Liver uptake
      Vmax_OATP1B1 <- Vmax_OATP1B1c*MW*60*24*SF_OATP1B1*VL_ec             # ug/d
      Km_OATP1B1 <- Km_OATP1B1c*MW                                        # ug/L (uM -> ug/L)
      Vmax_OATP1B3 <- Vmax_OATP1B3c*MW*60*24*SF_OATP1B3*VL_ec             # ug/d
      Km_OATP1B3 <- Km_OATP1B3c*MW                                        # ug/L (uM -> ug/L)
      
      # Biliary excretion
      VmaxBSEP <- VmaxBSEPc*MW*60*24*SF_BSEP*VL_ic         # ug/d
      KmBSEP <- KmBSEPc*MW                                 # ug/L (uM -> ug/L)
      
      # Renal clearance
      Vmax_OAT4 = Vmax_OAT4c*MW*60*24*SF_OAT*VPT           # ug/d (umol -> ug, min -> d)
      Km_OAT4 = Km_OAT4c*MW                                # ug/L (uM -> ug/L)
      
      
      ## Dose -------------------------
      
      if(t<expSTOP){DoseOn=1} else{DoseOn=0}
      
      ## Inhalation exposure ##
      DLung = expLung*BW*DoseOn         # ug, PFOA inhalation dose
      LungD = DLung/Tinput*(t %% tinterval<Tinput)
      
      ## Concentrations -------------------------
      
      CIL <- AIL/VIL               # ug/L, Intestinal lumen
      CI <- AI/VI                  # ug/L, Intestine 
      CVI <- CI/PI                 # ug/L, Intestine venous
      
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
      CVLu <- CLu/PLu              # ug/L, Lungs venous
      
      CAP <- AAP/VAP               # ug/L, Arterial plasma
      CVP <- AVP/VVP               # ug/L, Venous plasma
      
      
      ## Differential equations -------------------------
      
      dLuD = LungD - LuD             # ug/d, Inhaled dose input 
      
      dALu <- LuD + QC*(CVP - CVLu)                            # ug/d, Lungs
      
      
      dAIL <- - tco*AIL - CL_IL*CIL + 
        + (VmaxBSEP/(KmBSEP + (CL_ic*fu_Lic)))*CL_ic*fu_Lic          # ug/d, Intestine lumen
      
      dAI <- QI*(CAP - CVI) + CL_IL*CIL                      # ug/d, Intestinal
      
      dAFe <-  tco*AIL                                       # ug/d, Feces
      
      
      dAL_ec <- + QI*CVI + QL*CAP - (QI+QL)*CVL_ec + 
        - (Vmax_OATP1B1/(Km_OATP1B1 + (CL_ec*fup)))*CL_ec*fup +
        - (Vmax_OATP1B3/(Km_OATP1B3 + (CL_ec*fup)))*CL_ec*fup            # ug/d, Liver extracellular space (vascular + interstitial space)
      
      dAL_ic <- (Vmax_OATP1B1/(Km_OATP1B1 + (CL_ec*fup)))*CL_ec*fup +
        + (Vmax_OATP1B3/(Km_OATP1B3 + (CL_ec*fup)))*CL_ec*fup +
        - (VmaxBSEP/(KmBSEP + (CL_ic*fu_Lic)))*CL_ic*fu_Lic             # ug/d, Liver intracellular space
      
      
      dAPTT <- QK*(CAP - CPTT) + 
        + (Vmax_OAT4/(Km_OAT4+(CPTL*fu_PTL)))*CPTL*fu_PTL      # ug/d, Proximal tubule tissue 
      
      dAPTL <- + fup*GFR*CAP - QT*CPTL +
        - (Vmax_OAT4/(Km_OAT4+(CPTL*fu_PTL)))*CPTL*fu_PTL      # ug/d, Proximal tubule lumen    
      
      dARKT <- QK*(CPTT - CVRKT)                               # ug/d, Rest of kidney
      
      dARKL <- QT*CPTL - QUr*CRKL                              # ug/d, Rest of kidney lumen
      
      dAUr <- QUr*CRKL                                         # ug/d, Urine
      
      
      dAA <- QA*(CAP-CVA)                                      # ug/d, Adipose
      
      
      dAR <- QR*(CAP-CVR)                                      # ug/d, Rest
      
      
      dAAP <- - (QI + QL + QA + QR + QK)*CAP + QC*CVLu - fup*GFR*CAP    #+ QSt       # ug/d, Arterial Plasma
      dAVP <- (QL+QI)*CVL_ec + QK*CVRKT + QA*CVA + QR*CVR - QC*CVP      # ug/d, Venous Plasma
      
      
      # Mass Balance
      Atot <- LuD +
        ALu +
        AIL + AI + AFe + 
        AL_ec + AL_ic +
        APTT + APTL + ARKT + ARKL + AUr +
        AA + 
        AR +
        AAP + AVP
      
      dAin <- LungD 
      MB <- Ain - Atot + 1    
      
      
      # End
      
      list(c(dLuD,
             dALu, 
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
             dAAP,
             dAVP,
             dAin
      ), 
      c(CLu = CLu,
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
        CP = CAP + CVP,
        Atot = Atot, 
        MB = MB
      )
      )
    })
  }
  
  output_PBK <-lsoda(y = y, 
                      times = times, 
                      func = PBK.model,
                      parms = parms)
  
  return(as.data.frame(output_PBK))
}


