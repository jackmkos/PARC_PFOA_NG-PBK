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

  Parameters <- read_excel(here("Input/Parameters.GSA.xlsx")) 
  
  P <- Parameters %>% select(Abbreviation, Distribution, Initial_value, Binf, Bsup, Mean, sdlog, zscore) %>%
    filter(Abbreviation %in% c("expOral",
                               "fup", 
                               "SF_OAT", 
                               "R_PTL", 
                               "Km_OAT4c", 
                               "QT",
                               "BW",
                               "QC",
                               "SF_OATP1B1",
                               "KmBSEPc",
                               "GFRc",
                               "Vmax_OAT4c",
                               "SF_BSEP",
                               "Km_OATP1B1c",
                               "VKc",
                               "QKc",
                               "VILc",
                               "VPTc",
                               "VmaxBSEPc",
                               "PLc",
                               "Vmax_OATP1B1c",
                               "VL_icc",
                               "Hct",
                               "SF_OATP1B3",
                               "R_L_ec")) # "SA_SI", "QKc", "SF_OATP1B1", "Vmax_OAT4c", "Km_OATP1B1c", "tco", "Papp_SI", "VPc", "VKc", "VPTc", "Hct", "PAc", "Vmax_OATP1B1c", "Km_OATP1B3c"))

  parm.c <- as.list(setNames(Parameters$Initial_value, Parameters$Abbreviation))

  # SET eFAST EXPERIMENT DESIGN ####
  # ---------------------------------------------------------------------------- #
  
  GSA.parms <- P$Abbreviation

  eFAST.factors <- nrow(P)

  # Define q: q needs to be a list of character strings, giving the names of the quantile functions
  q <- ifelse(P$Distribution == "LogNormal", "qlnorm", "qunif") # if distribution is lognormal then q should be qlnorm, if not then qunif

  # Define q.arg: q.arg needs to be a list of lists
  q.arg <- lapply(1:nrow(P), function(i) { # apply the function to eachone of the rows of P
    if (P$Distribution[i] == "Uniform") { #if the distribution in uniform
      list(min = P$Binf[i], max = P$Bsup[i]) #create a list containing min (Bing of row i) and max (Bsup of row i) per row
    } else {
      list(meanlog = log(P$Mean[i]), sdlog = P$sdlog[i]) #if not then create a list containing mean and sdlog
    }
  })

  # q.arg <- lapply(1:nrow(P), function(i) {
  #   if (P$Distribution[i] == "Uniform") {
  #     list(min = P$Binf[i], max = P$Bsup[i])
  #   } else {
  #     list(meanlog = log(P$Mean[i]), sdlog = P$sdlog[i])
  #   }
  # })
  # names(q.arg) <- GSA.parms
  # 
  # 
  #   # Perform eFAST test
  # set.seed(1234)
  # 
  # length(q) == length(q.arg)
  # 
  # eFAST <- fast99(
  #   model = NULL, #PBK_4_GSA,
  #   factors = GSA.parms, # These are the parameters that will be varying in the model
  #   n = 1000,            # Integer giving the sample size, i.e. the length of the discretization of the s-space
  #   q = q,
  #   q.arg = q.arg
  # )
  # 
  # dim(eFAST$X) # dataframe of: n factors(number of parameters) * n n(1000) number of observations, of n factors(number of parameters) number of variables
  # INITIALdesign <- eFAST$X
  # write.csv(INITIALdesign, "INITIALeFAST.ExpDesign.csv", row.names = FALSE)
  # 
  # pdf("INITIALDistr_ParFull.pdf")
  # par(mfrow=c(4,3))
  # for( i in 1:ncol(INITIALdesign)){ hist(INITIALdesign[,i], breaks=100, col="purple", main=colnames(INITIALdesign)[i] ) }
  # dev.off()
  # save(eFAST, file = "INITIALExperienceFull.RData")
  # 
  # for (i in 1:nrow(P)) {
  #   if (P$Distribution[i] != "Uniform") {
  #     par_name <- GSA.parms[i]
  #     min_val <- exp(log(P$Mean[i]) - P$zscore[i] * P$sdlog[i])
  #     max_val <- exp(log(P$Mean[i]) + P$zscore[i] * P$sdlog[i])
  # 
  #     eFAST$X[, par_name] <- pmax(eFAST$X[, par_name], min_val)
  #     eFAST$X[, par_name] <- pmin(eFAST$X[, par_name], max_val)
  #   }
  # }
  # 
  # design <- eFAST$X
  # write.csv(design, "eFAST.ExpDesign.csv", row.names = FALSE)
  # 
  # pdf("Distr_ParFull.pdf")
  # par(mfrow=c(4,3))
  # for( i in 1:ncol(design)){ hist(design[,i], breaks=100, col="purple", main=colnames(design)[i] ) }
  # dev.off()
  # save(eFAST, file = "ExperienceFull.RData")
  # 
  # GSA_model <- function(parm.GSA, parm.c) {
  #   
  #   # Final parameters per test
  #   names(parm.GSA) <- GSA.parms  
  #   parm.final <- parm.c
  #   parm.final[names(parm.GSA)] <- parm.GSA
  #   
  #   # (Optional) Print parameters for debugging
  #   print(parm.final[names(parm.GSA)])
  #   
  #   A_init <- c(OD = 0,
  #               ASk = 0,
  #               AIL = 0, AI = 0, AFe = 0,
  #               AL_ec = 0, AL_ic = 0,
  #               APTT = 0, APTL = 0, ARKT = 0, ARKL = 0, AUr = 0,
  #               AA = 0,
  #               AR = 0,
  #               AP = 0,
  #               Ain = 0)
  #   
  #   out <- lsoda(y = A_init,
  #                times = seq(0, 10, 1),
  #                func = PBK.model,
  #                parms = parm.final)
  #   
  #   # Return selected outputs
  #   CP <- out[,"CP"]
  #   AUC <- sum(diff(out[,"time"]) * (out[-nrow(out),"CP"] + out[-1,"CP"]))/2
  #   
  #   return(c(CP = CP, AUC = AUC)) #CP = CP
  # }
  # 
  # 
  # results <- apply(design, 1, function(row) { 
  #   GSA_model(row, parm.c = parm.c)  # Run model row-wise from design (i.e each row is an eFAST test). parm.GSA = row
  # })
  # 
  # AUC <- results["AUC",]
  # CP <- results[1:11,]
  # eFAST$y <- CP
  # tell(eFAST)
  # plot(eFAST)
 