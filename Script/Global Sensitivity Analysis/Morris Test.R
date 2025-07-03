  # --------------------------------------------------------------------------- #
  # Script for running the Morris Test
  # CP, 29-03-2025
  # --------------------------------------------------------------------------- #
  
  rm(list=ls()) # to clear out the global environment
  
  # Load packages
  library(here)
  library(ggplot2)
  library(deSolve)
  library(tidyverse)
  library(sensitivity)
  library(ggrepel)
  library(scales)
  library(patchwork)
  library(readxl)
  
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
  
  # Set output storage directory
  # OUTPUT <- here("Output", format(Sys.Date(), "%Y-%m-%d"), format(Sys.time(), "%H-%M-%S"))
  # dir.create(OUTPUT, recursive = TRUE)
  
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
  
  
  # SENSITIVITY FUNCTION ####
  # ---------------------------------------------------------------------------- #
  
  ## Morris test ####
  SENSI_model <- function(parm.c){
    
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
                     times = seq(0, 20*365,by=10),
                     func =  ORAL_PBK.model,
                     parms = parm.c)
    
    CP <- PBK.out[,"CP"]
    CPTT <- PBK.out[,"CPTT"]
    CL_ic <- PBK.out[,"CL_ic"]
    CL_ec <- PBK.out[,"CL_ec"]
    
    return(c(CP=CP, CPTT = CPTT, CL_ic = CL_ic, CL_ec = CL_ec))
    
  }
  
  # SET MORRIS PARAMETERS ####
  # ---------------------------------------------------------------------------- #
  
  # Parameter upper and lower bounds
  
  ParametersMorris <- read_excel(here("Input/Parameters.GSA.xlsx")) 

  P <- ParametersMorris %>% select(Abbreviation, Initial_value, Binf, Bsup) %>% 
    tibble::column_to_rownames("Abbreviation")
  
  binf <- P$Binf
  bsup <- P$Bsup
  
  Morris.factors <- nrow(P)  
  
  # Experiment design
  Morris.design <- list(type = "oat", levels = 6, grid.jump = 3)
  
  # Repetitions
  Morris.r <- 50
  
  # Perform Morris test
  Morris <- morris(model =  NULL, 
                   factors = Morris.factors, # number of parameters
                   r = Morris.r,  # number of repetitions
                   design = Morris.design,  
                   binf = binf,  
                   bsup = bsup, 
                   scale = TRUE)
  
  # Save Morris test design
  design <- Morris$X
  colnames(design) <- row.names(P)
  
  write.csv(design, "MorrisExpDesign.csv", row.names = FALSE)
  
  pdf("Distr_ParFull.pdf")
  par(mfrow=c(4,3))
  for( i in 1:ncol(design)){ hist(design[,i], breaks=100, col="purple", main=colnames(design)[i] ) }
  dev.off()
  save(Morris, file = "ExperienceFull.RData")
  
  
  # RUN MORRIS DESIGN ####
  # ---------------------------------------------------------------------------- #
  
  # Runs r (n of repetitions) * (param+1) simulations * (n) model outputs * (n) model outputs per time
  results <- apply(design, 1,  SENSI_model) 
  
  y <- results
  save.y <- as.data.frame(y) #%>% mutate(Tout = "times = seq(0, 5*365,by=1/10)")
  save(save.y, file = "y.RData")
  
  load(here("y.RData"))
  
  
  # CALCULATE MORRIS INDICES ####
  # ---------------------------------------------------------------------------- #
  
  sim.results.morris <- as.matrix(y)
  
  # CB.sim.results.morris <- sim.results.morris[str_detect(rownames(sim.results.morris), "CB") &
  #                                               !str_detect(rownames(sim.results.morris), "aZEL"), ]
  # 
  # CBaZEL.sim.results.morris <- sim.results.morris[str_detect(rownames(sim.results.morris), "CB") &
  #                                                   str_detect(rownames(sim.results.morris), "aZEL"), ]
  # 
  # CLGLU.sim.results.morris <- sim.results.morris[str_detect(rownames(sim.results.morris), "CL") &
  #                                                  str_detect(rownames(sim.results.morris), "GLU"), ]
  # 
  # CL.sim.results.morris <- sim.results.morris[str_detect(rownames(sim.results.morris), "CL") &
  #                                               !str_detect(rownames(sim.results.morris), "GLU"), ]
  # 
  
  
  
  mu = mu.star = sigma = NULL
  
  # Tell morris y for each model output
  for (i in (1:length(sim.results.morris[,1])))
  {
    tell(Experience, sim.results.morris[i,])
    mu = rbind( mu, apply(Experience$ee, 2, mean, na.rm = TRUE))
    mu.star = rbind( mu.star, apply(Experience$ee, 2, function(x) mean(abs(x), na.rm = TRUE)))
    sigma = rbind( sigma, apply(Experience$ee, 2, sd, na.rm = TRUE))
    
  }
  plot(Experience)
  
  # Save mu.star, mu, sigma, scaled indices and global indices
  rownames( mu.star ) = rownames( mu ) = rownames( sigma ) = row.names(sim.results.morris)
  
  mu.star_scale = mu.star/apply(mu.star, 1, mean, na.rm=TRUE) 
  sigma_scale = sigma/apply(sigma, 1, mean, na.rm=TRUE)
  
  # Calculate global indices
  IGlobal_mu.star = apply(mu.star_scale, 2, mean, na.rm=TRUE)
  IGlobal_sigma = apply(sigma_scale, 2, mean, na.rm = TRUE )
  
  Morris_IGlobal.df <- data.frame(
    IGlobal_mu.star = apply(mu.star_scale, 2, mean, na.rm=TRUE),
    IGlobal_sigma = apply(sigma_scale, 2, mean, na.rm = TRUE )
  )
  
  IGlobal <- data.frame(IGlobal_sigma, IGlobal_mu.star) 
  IGlobal <- IGlobal %>% mutate(parm_names = rownames(IGlobal))
  
  Morris.indices <- list(
    mu = mu,
    mu_star = mu.star,
    sigma = sigma,
    mu_star_scale = mu.star_scale,
    sigma_scale = sigma_scale,
    IGlobal = IGlobal
  )
  
  save(Morris.indices, file = here("Morris.indices.RData"))
  
  
  # PLOT MORRIS RESULTS ####
  # ---------------------------------------------------------------------------- #
  
  ggplot(IGlobal) +
    aes(x = IGlobal_mu.star, y = IGlobal_sigma, label = rownames(IGlobal)) +
    # xlim (min,2.5) + #to zoom in/out, determine max value
    # ylim (min,2.5) + #to zoom in/out, determine max value
    geom_point()
  
  
  morris_plot <- IGlobal %>% 
    # filter(mu.star_scale >= 0.01) %>%
    ggplot(aes(x = IGlobal_mu.star, y = IGlobal_sigma, label = parm_names)) +
    geom_point(alpha = 0.5, color = "darkblue") +
    geom_text_repel(
      size = 7,
      angle = 30,
      max.overlaps = Inf,
      segment.color = "darkblue",
      segment.size = 0.1,
      segment.alpha = 0.5,
      segment.linetype = "dotted",
      direction = "x",
      ylim = c(7, 7),
      point.padding = 0.5,
      box.padding = 0.2,
      set.seed(123)
    ) +
    coord_cartesian(clip = "off") +
    labs(
      title = "Morris Sensitivity Analysis",
      x = "μ* (Mean Effect)",
      y = "σ (Standard Deviation)"
      # subtitle = "Parameters with relative μ* ≥ 5%"
    ) +
    # theme_minimal() +
    CP_theme+
    theme(plot.margin = margin(0, 0.15, 0, 0, "cm"))
  morris_plot
  ggsave(
    filename = here("morrisGI_plot.png"),
    plot = morris_plot,
    dpi = 1000, 
    width = 12, height = 9, units = "cm"
  )
  
  #Zoom out
  morris_plot <- IGlobal %>% 
    filter(IGlobal_mu.star >= 0.15) %>%
    ggplot(aes(x = IGlobal_mu.star, y = IGlobal_sigma, label = parm_names)) +
    geom_point(alpha = 0.5, color = "darkblue") +
    geom_text_repel(
      size = 9,
      angle = 30,
      max.overlaps = Inf,
      segment.color = "darkblue",
      segment.size = 0.1,
      segment.alpha = 0.5,
      segment.linetype = "dotted",
      direction = "x",
      ylim = c(7, 7),
      point.padding = 0.4,
      box.padding = 0.1,
      set.seed(123)
    ) +
    coord_cartesian(clip = "off") +
    labs(
      title = "Morris Sensitivity Analysis",
      x = "μ* (Mean Effect)",
      y = "σ (Standard Deviation)",
      subtitle = "Parameters with μ* ≥ 0.15"
    ) +
    CP_theme+
    # theme_minimal() +
    theme(plot.margin = margin(0, 0.15, 0, 0, "cm"))
  morris_plot
  ggsave(
    filename = here("morrisZoomOut_plot.png"),
    plot = morris_plot,
    dpi = 1000, 
    width = 12, height = 9, units = "cm"
  )
  
  #Zoom in
  morris_plot <- IGlobal %>% 
    filter(IGlobal_mu.star <=0.15) %>%
    ggplot(aes(x = IGlobal_mu.star, y = IGlobal_sigma, label = parm_names)) +
    geom_point(alpha = 0.5, color = "darkblue") +
    geom_text_repel(
      size = 9,
      angle = 30,
      max.overlaps = Inf,
      segment.color = "darkblue",
      segment.size = 0.1,
      segment.alpha = 0.5,
      segment.linetype = "dotted",
      direction = "x",
      ylim = c(0.035, 0.035),
      point.padding = 0.4,
      box.padding = 0.3,
      set.seed(123)
    ) +
    coord_cartesian(clip = "off") +
    labs(
      title = "Morris Sensitivity Analysis",
      x = "μ* (Mean Effect)",
      y = "σ (Standard Deviation)",
      subtitle = "Parameters with μ* <= 0.15"
    ) +
    CP_theme+
    # theme_minimal() +
    theme(plot.margin = margin(0, 0.15, 0, 0, "cm"))
  morris_plot
  ggsave(
    filename = here("morrisZoomIn_plot.png"),
    plot = morris_plot,
    dpi = 1000, 
    width = 12, height = 9, units = "cm"
  )
  # D <- ggplotly(Plot)
  