  # --------------------------------------------------------------------------- #
  # PBK MODEL FOR PFOA, TO BE USED TOGETHER WITH THE LATEST HBM DATA
  # Core PBK model, minimalistic
  # CP, 29-03-2025
  # --------------------------------------------------------------------------- #
  
  rm(list=ls()) # to clear out the global environment
  
  # # Set output storage directory
  # OUTPUT <- here("Output", format(Sys.Date(), "%Y-%m-%d"), format(Sys.time(), "%H-%M-%S"))
  # dir.create(OUTPUT, recursive = TRUE)
  # 
  
  # Load packages
  library(here)
  library(ggplot2)
  library(deSolve)
  library(tidyverse)
  library(sensitivity)
  library(ggrepel)
  library(scales)
  library(patchwork)
  library(pracma)
  library(readxl)
  library(purrr)

  # For plotting
  library(showtext)
  font_add(family = "Garamond", regular = "GARA.TTF")
  showtext_auto()
  CP_theme <- theme_minimal() +
    theme(
      axis.text = element_text(size = 56),
      axis.title = element_text(size = 62),
      plot.title = element_text(size = 65, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 56, hjust = 0.5),
      legend.position = "none",
      plot.margin = margin(0.3, 0.3, 0.3, 0.3, "cm") 
    )
  
  # The PBK model ####
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
      GFR <- GFR
      QT <- QTc*GFR                      # L/d, Proximal tubule fluid flow
      
      tco <- tco                    # /d, Bowel residence times in the colon
      
      
      ### Physicochemical ----
      MW <- 414.07
      
      PI <- PIc * fup     # Intestinal
      PL <- PLc * fup     # Liver
      PK <- PKc * fup     # Kidney
      PA <- PAc * fup     # Adipose
      PR <- PRc * fup     # Rest
      
      
      # Fraction unbound, calculated based on Poulin and Haddad, 2018 https://doi.org/10.1016/j.xphs.2018.03.012
      # Equation was adapted to not account for fraction unionised
      # OAT and OATP transporters transport the ionised compound, given that the ratio of fraction ionised at plasma to cellular pH is 1, this can be ignored (fraction unionised of PFOA is 0.9999923 at pH 7.4 and 0.9999963 at pH 7)
      fu_PTL <- R_PTL*fup/(1 + ((R_PTL-1)*fup)) # Proximal tubule lumen
      fu_Lic <- R_L_ec*fup/(1 + ((R_L_ec-1)*fup)) # Liver intracellular space
      
      ### Kinetic ----
      
      # Gastro-intestinal uptake
      CL_IL <- (Papp_SI*SA_SI*1e-3)*60*60*24    # L/d, Intestinal lumen to intestinal tissue (calculations: cm/s -> L/s /1000 -> L/d *60*60*24)
      
      VmaxOATP2B1 <- Vmax_OATP2B1c*MW*60*24*SF_OATP2B1*VIL # ug/d
      Km_OATP2B1 <- Km_OATP2B1c*MW                         # ug/L (uM -> ug/L)
      
      # Liver uptake
      Vmax_OATP1B1 <- Vmax_OATP1B1c*MW*60*24*SF_OATP1B1*VL_ec             # ug/d
      Km_OATP1B1 <- Km_OATP1B1c*MW                                        # ug/L (uM -> ug/L)
      Vmax_OATP1B3 <- Vmax_OATP1B3c*MW*60*24*SF_OATP1B3*VL_ec             # ug/d
      Km_OATP1B3 <- Km_OATP1B3c*MW                                        # ug/L (uM -> ug/L)
      
      # Biliary excretion
      VmaxBSEP <- Vmax_BSEPc*MW*60*24*SF_BSEP*VL_ic         # ug/d
      Km_BSEP <- Km_BSEPc*MW                                 # ug/L (uM -> ug/L)
      
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
        - (VmaxOATP2B1/(Km_OATP2B1 + CIL))*CIL +
        + (VmaxBSEP/(Km_BSEP + (CL_ic*fu_Lic)))*CL_ic*fu_Lic          # ug/d, Intestine lumen
      
      dAI <- QI*(CP - CVI) + CL_IL*CIL + 
        + (VmaxOATP2B1/(Km_OATP2B1 + CIL))*CIL                     # ug/d, Intestinal
      
      dAFe <-  tco*AIL                                       # ug/d, Feces
      
      
      dAL_ec <- + QI*CVI + QL*CP - (QI+QL)*CVL_ec + 
        - (Vmax_OATP1B1/(Km_OATP1B1 + (CL_ec*fup)))*CL_ec*fup +
        - (Vmax_OATP1B3/(Km_OATP1B3 + (CL_ec*fup)))*CL_ec*fup            # ug/d, Liver extracellular space (vascular + interstitial space)
      
      dAL_ic <- (Vmax_OATP1B1/(Km_OATP1B1 + (CL_ec*fup)))*CL_ec*fup +
        + (Vmax_OATP1B3/(Km_OATP1B3 + (CL_ec*fup)))*CL_ec*fup +
        - (VmaxBSEP/(Km_BSEP + (CL_ic*fu_Lic)))*CL_ic*fu_Lic             # ug/d, Liver intracellular space
      
      
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
  
  
  # # SENSITIVITY FUNCTION ####
  # # ---------------------------------------------------------------------------- #
  
  ## eFAST ####

  GSA_model <- function(parm.GSA, parm.c) {
    
    # Final parameters per test
    names(parm.GSA) <- GSA.parms
    parm.final <- parm.c
    parm.final[names(parm.GSA)] <- parm.GSA
    
    # (Optional) Print parameters for debugging
    print(parm.final[names(parm.GSA)])
    
    A_init = c(OD = 0,
               AIL = 0,
               AI = 0, 
               AFe = 0,
               AL_ec = 0,
               AL_ic = 0,
               APTT = 0,
               APTL = 0,
               ARKT = 0,
               ARKL = 0,
               AUr = 0,
               AA = 0,
               AR = 0, 
               AP = 0, 
               Ain = 0)
    
    PBK.out <- lsoda(y = A_init,
                     times = seq(0, 20*365, by=10),
                     func = PBK.model,
                     parms = parm.final)
    
    CP <- PBK.out[,"CP"]
    CPTT <- PBK.out[,"CPTT"]
    CL_ic <- PBK.out[,"CL_ic"]
    CL_ec <- PBK.out[,"CL_ec"]
    
    return(c(CP=CP, CPTT = CPTT, CL_ic = CL_ic, CL_ec = CL_ec)) 
  }
  
  # SET eFAST PARAMETERS ####
  # ---------------------------------------------------------------------------- #
  
  # Parameter upper and lower bounds

  Parameters <- read_excel(here("Input", "ParameterseFAST.F.xlsx")) #"Input/Parameters.GSA.xlsx" 
  # Parameters <- read_excel(here("Output/efast1000/Parameters.eFAST.F.xlsx")) #"Input/Parameters.GSA.xlsx" 

  P <- Parameters %>% select(Parameter, Distribution, Value, Std, CV, eFAST) 
  # P <- Parameters %>% select(Abbreviation, Distribution, Initial_value, Binf, Bsup, sdlog, zscore) 
  
  parm.c <- as.list(setNames(Parameters$Value, Parameters$Parameter))
  # parm.c <- as.list(setNames(Parameters$Initial_value, Parameters$Abbreviation))

  # SET eFAST EXPERIMENT DESIGN ####
  # ---------------------------------------------------------------------------- #

  GSA.P <- P %>% filter(eFAST == "Y") %>% 
    mutate(
      Mean = as.numeric(Value),
      Std = as.numeric(Std),
      Binf = Value*0.9, 
      Bsup = Value*1.1,
      meanlog = if_else(Distribution == "LogNormal",
                       log(Mean / sqrt(1 + (Std/Mean)^2)), 
                       NA),
      sdlog = if_else(Distribution == "LogNormal",
                      sqrt(log(1 + (Std/Mean)^2)),
                      NA),
      zscore = c(3,NA,NA,5,3.5,3,3,3,NA,NA,NA,NA,4,4,3,NA,NA,4,4),
      q = ifelse(Distribution == "LogNormal", "qlnorm", "qunif"))
  
  eFAST.factors <- nrow(GSA.P)
  
  GSA.parms <- GSA.P$Parameter
  # GSA.parms <- P$Abbreviation
  
  # # P$Mean <- as.numeric(P$Mean)
  # P$sdlog <- as.numeric(P$sdlog)

  # Define q: q needs to be a list of character strings, giving the names of the quantile functions
  # q <- ifelse(P$Distribution == "LogNormal", "qlnorm", "qunif") 
  
  # Define q.arg: q.arg needs to be a list of lists
  GSA.P$q.arg <- lapply(1:nrow(GSA.P), function(i) { # apply the function to eachone of the rows of P
    if (GSA.P$q[i] == "qunif") { #if the distribution in uniform
      list(min = GSA.P$Binf[i], max = GSA.P$Bsup[i]) #create a list containing min (Bing of row i) and max (Bsup of row i) per row
    } else {
      list(meanlog = GSA.P$meanlog[i], sdlog = GSA.P$sdlog[i]) #if not then create a list containing mean and sdlog
    }
  })
  
  q <- GSA.P$q
  names(q) <- GSA.P$Parameter
  
  q.arg <- GSA.P$q.arg
  names(q.arg) <- GSA.P$Parameter
  # eFAST test function
  set.seed(1234)
  
  length(GSA.P$q) == length(GSA.P$q.arg)

  Experience <- fast99(
    model = NULL, #PBK_4_GSA,
    factors = GSA.parms, # These are the parameters that will be varying in the model
    n = 10000,   #1000         # Integer giving the sample size, i.e. the length of the discretization of the s-space
    q = q,
    q.arg = q.arg
  )

  INITIALdesign <- Experience$X
  write.csv(INITIALdesign, "INITIALeFAST.ExpDesign.csv", row.names = FALSE)
  
  pdf("INITIALDistr_ParFull.pdf")
  par(mfrow=c(4,3))
  for( i in 1:ncol(INITIALdesign)){ hist(INITIALdesign[,i], breaks=100, col="purple", main=colnames(INITIALdesign)[i] ) }
  dev.off()
  save(Experience, file = "INITIALExperienceFull.RData")
  
  # P$zscore <- as.numeric(P$zscore)
  
  for (i in 1:nrow(GSA.P)) {
    if (GSA.P$Distribution[i] != "Uniform" && GSA.P$Parameter[i] != "fup" ) {
      
      par_name <- GSA.parms[i]
     
      min_val <- exp(GSA.P$meanlog[i] - GSA.P$zscore[i] * GSA.P$sdlog[i])
      max_val <- exp(GSA.P$meanlog[i] + GSA.P$zscore[i] * GSA.P$sdlog[i])
      
      Experience$X[, par_name] <- pmax(Experience$X[, par_name], min_val)
      Experience$X[, par_name] <- pmin(Experience$X[, par_name], max_val)
      
    } 
    # else if (GSA.P$Parameter[i] == "fup") {
    # 
    #   log_fup <- log(Experience$X[, "fup"] + 1e-10)
    #   log_fup <- pmin(log_fup, quantile(log_fup, 0.99))
    #   Experience$X[,"fup"] <- exp(log_fup)
    # 
    # }
  }

  design <- Experience$X
  write.csv(design, "eFAST.ExpDesign.csv", row.names = FALSE)
  
  pdf("Distr_ParFull.pdf")
  par(mfrow=c(4,3))
  for( i in 1:ncol(design)){ hist(design[,i], breaks=100, col="purple", main=colnames(design)[i] ) }
  dev.off()
  save(Experience, file = "ExperienceFull.RData")
  
  
  # RUN eFAST DESIGN ####
  # ---------------------------------------------------------------------------- #
  results <- apply(design, 1, function(row) {
    GSA_model(row, parm.c = parm.c)  
  })
  
  y <- results
  y <- as.data.frame(y) 
  save(y, file = "eFAST.y.RData")

  # CALCULATE eFAST INDICES ####
  # ---------------------------------------------------------------------------- #
  
  #All outputs
  sim.results.eFAST <- as.matrix(y)
  
  for (i in (1:length(sim.results.eFAST[,1])))
  {
    tell(Experience, sim.results.eFAST[i,])
  }
  
  plot(Experience)
  
  Variance <- Experience$V # Total variance
  names(Variance) <- GSA.parms
  
  Done <- Experience$D1
  names(Done) <- GSA.parms
  
  Dt <- Experience$Dt
  names(Dt) <- GSA.parms
  
  first_order <- Done / Variance  # D1 is the estimated Variance of the Conditional Expectation (VCE) with respect to each factor, normalized by total variance
  names(first_order) <- GSA.parms
  
  total_order <- Dt / Experience$V  # Dt is the estimated VCE with respect to each factor complementary set of factors ("all but Xi"), normalized by total variance
  names(total_order) <- GSA.parms
  
  lowry_data <- data.frame(
    Parameter = GSA.parms,
    Main.Effect = first_order,
    Interaction = total_order) %>% 
    mutate(
      Main.Effect = Main.Effect/sum(Main.Effect), # to normalise
      Interaction = Interaction/sum(Interaction)  # to normalise
    ) 
  
  write.csv(lowry_data, file = here("eFASTresults.csv"), row.names = FALSE)
  
  # CP
  CP.sim.results.eFAST <- sim.results.eFAST[str_detect(rownames(sim.results.eFAST), "CP") &
                                              !str_detect(rownames(sim.results.eFAST), "TT"), ]
  
  for (i in (1:length(CP.sim.results.eFAST[,1])))
  {
    tell(Experience, CP.sim.results.eFAST[i,])
  }
  
  plot(Experience)
  
  Variance <- Experience$V # Total variance
  names(Variance) <- GSA.parms
  
  Done <- Experience$D1
  names(Done) <- GSA.parms
  
  Dt <- Experience$Dt
  names(Dt) <- GSA.parms
  
  first_order <- Done / Variance  # D1 is the estimated Variance of the Conditional Expectation (VCE) with respect to each factor, normalized by total variance
  names(first_order) <- GSA.parms
  
  total_order <- Dt / Experience$V  # Dt is the estimated VCE with respect to each factor complementary set of factors ("all but Xi"), normalized by total variance
  names(total_order) <- GSA.parms
  
  lowry_data <- data.frame(
    Parameter = GSA.parms,
    Main.Effect = first_order,
    Interaction = total_order) %>% 
    mutate(
      Main.Effect = Main.Effect/sum(Main.Effect), # to normalise
      Interaction = Interaction/sum(Interaction)  # to normalise
    ) 
  
  write.csv(lowry_data, file = here("CP.eFASTresults.csv"), row.names = FALSE)
  
  #CPTT
  CPTT.sim.results.eFAST <- sim.results.eFAST[str_detect(rownames(sim.results.eFAST), "CPTT"), ]
  
  for (i in (1:length(CPTT.sim.results.eFAST[,1])))
  {
    tell(Experience, CPTT.sim.results.eFAST[i,])
  }
  
  plot(Experience)
  
  Variance <- Experience$V # Total variance
  names(Variance) <- GSA.parms
  
  Done <- Experience$D1
  names(Done) <- GSA.parms
  
  Dt <- Experience$Dt
  names(Dt) <- GSA.parms
  
  first_order <- Done / Variance  # D1 is the estimated Variance of the Conditional Expectation (VCE) with respect to each factor, normalized by total variance
  names(first_order) <- GSA.parms
  
  total_order <- Dt / Experience$V  # Dt is the estimated VCE with respect to each factor complementary set of factors ("all but Xi"), normalized by total variance
  names(total_order) <- GSA.parms
  
  lowry_data <- data.frame(
    Parameter = GSA.parms,
    Main.Effect = first_order,
    Interaction = total_order) %>% 
    mutate(
      Main.Effect = Main.Effect/sum(Main.Effect), # to normalise
      Interaction = Interaction/sum(Interaction)  # to normalise
    ) 
  
  write.csv(lowry_data, file = here("CPTT.eFASTresults.csv"), row.names = FALSE)
  
  #CPTL
  CPTL.sim.results.eFAST <- sim.results.eFAST[str_detect(rownames(sim.results.eFAST), "CPTL"), ]
  
  for (i in (1:length(CPTL.sim.results.eFAST[,1])))
  {
    tell(Experience, CPTL.sim.results.eFAST[i,])
  }
  
  plot(Experience)
  
  Variance <- Experience$V # Total variance
  names(Variance) <- GSA.parms
  
  Done <- Experience$D1
  names(Done) <- GSA.parms
  
  Dt <- Experience$Dt
  names(Dt) <- GSA.parms
  
  first_order <- Done / Variance  # D1 is the estimated Variance of the Conditional Expectation (VCE) with respect to each factor, normalized by total variance
  names(first_order) <- GSA.parms
  
  total_order <- Dt / Experience$V  # Dt is the estimated VCE with respect to each factor complementary set of factors ("all but Xi"), normalized by total variance
  names(total_order) <- GSA.parms
  
  lowry_data <- data.frame(
    Parameter = GSA.parms,
    Main.Effect = first_order,
    Interaction = total_order) %>% 
    mutate(
      Main.Effect = Main.Effect/sum(Main.Effect), # to normalise
      Interaction = Interaction/sum(Interaction)  # to normalise
    ) 
  
  write.csv(lowry_data, file = here("CPTL.eFASTresults.csv"), row.names = FALSE)
  
  #CL_ic
  CL_ic.sim.results.eFAST <- sim.results.eFAST[str_detect(rownames(sim.results.eFAST), "CL_ic"), ]

  for (i in (1:length(CL_ic.sim.results.eFAST[,1])))
  {
    tell(Experience, CL_ic.sim.results.eFAST[i,])
  }
  
  plot(Experience)
  
  Variance <- Experience$V # Total variance
  names(Variance) <- GSA.parms
  
  Done <- Experience$D1
  names(Done) <- GSA.parms
  
  Dt <- Experience$Dt
  names(Dt) <- GSA.parms
  
  first_order <- Done / Variance  # D1 is the estimated Variance of the Conditional Expectation (VCE) with respect to each factor, normalized by total variance
  names(first_order) <- GSA.parms
  
  total_order <- Dt / Experience$V  # Dt is the estimated VCE with respect to each factor complementary set of factors ("all but Xi"), normalized by total variance
  names(total_order) <- GSA.parms
  
  lowry_data <- data.frame(
    Parameter = GSA.parms,
    Main.Effect = first_order,
    Interaction = total_order) %>% 
    mutate(
      Main.Effect = Main.Effect/sum(Main.Effect), # to normalise
      Interaction = Interaction/sum(Interaction)  # to normalise
    ) 
  
  write.csv(lowry_data, file = here("CL_ic.eFASTresults.csv"), row.names = FALSE)
  
  
  
  
  
  
  # lowry_data <- read.csv(here("eFASTresults.csv"))
  
  ordered_data <- lowry_data %>%
    arrange(desc(Main.Effect)) %>%
    mutate(
      Parameter = factor(Parameter, levels = Parameter), 
      Total.effect = Main.Effect + Interaction,
      Cumulative.main = cumsum(Main.Effect))
  long_data <- ordered_data %>% pivot_longer(c(Main.Effect, Interaction))
  long_data$name <- factor(long_data$name, levels = c("Main.Effect", "Interaction"))
  
  ordered_data$Cumulative.total = cumsum(ordered_data$Total.effect)
  
  lowry_plot <- 
    ggplot(long_data) +
    geom_col(aes(x = Parameter, y = value, fill = name),
             position = position_stack(reverse = TRUE),
             width = 0.8,
             alpha = 0.8) +
    geom_ribbon(
      data = ordered_data,
      aes(x = as.numeric(Parameter),
          ymin = lag(Cumulative.main, default = 0),
          ymax = pmin(Cumulative.total,1)),
      fill = "grey50", alpha = 0.3, color = "grey50", linewidth = 0.2
    ) +
    scale_fill_manual(
      values = c("Main.Effect" = "darkblue", "Interaction" = "blueviolet"),
      labels = c("Main Effect", "Interaction")
    ) +
    scale_y_continuous(
      limits = c(0, 1),
      expand = expansion(mult = c(0, 0.05)),
      labels = percent_format(),
      name = "Sensitivity Index"
    ) +
    labs(
      x = "Parameter",
      title = "eFAST Sensitivity Analysis"
    ) +
    CP_theme +
    theme(
      legend.position = "top",
      legend.title = element_blank(),
      legend.text = element_text(size = 42),
      legend.key.size = unit(0.5, "lines"),
      legend.spacing.y = unit(0, "cm"),
      legend.margin = margin(0, 0, 0, 0),
      axis.text.x = element_text(angle = 45, hjust = 1,),
      plot.margin = margin(0, 0, 0, 0.15, "cm")
    )
  
  lowry_plot
  ggsave(
    filename = here("lowry_plot.png"),
    plot = lowry_plot,
    dpi = 1000, 
    width = 12, height = 9, units = "cm"
  )
 