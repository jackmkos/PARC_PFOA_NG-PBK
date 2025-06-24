# --------------------------------------------------------------------------- #
# SCRIPT FOR CALCULATING MODEL RESULTS
# By: Chrysanthi Pachoulide
# Date: 14-04-2025
# --------------------------------------------------------------------------- #

# Without physiological changes with lifestage (depending on age) ####
RUNandOUT <- function(exposure_type, 
                      exp_Oral = NULL, exp_Dermal = NULL, exp, 
                      Tinput, tinterval, expSTOP, 
                      expAGE = NULL, expBW = NULL, sex = "M", 
                      Tstart, Tstop, Dt) {
  
  
  # Load R files 
  source(here("Script", "CALC_Parameters.R"), local = TRUE)
  source(here("Script", "PBK_model.R"), local = TRUE)  
  
  ## Add state ####
  Astate <- switch (exposure_type,
                    "Oral" = c(OD = 0,
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
                               Ain = 0),
                    "Dermal" = c(DD = 0,
                                 ASkB = 0,
                                 ASk = 0, 
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
                                 Ain = 0),
                    "Oral_Dermal" = c(OD = 0,
                                      DD = 0,
                                      ASkB = 0,
                                      ASk = 0, 
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
                                      Ain = 0),
                    "Inhalation" = c(LuD = 0,
                                     ALu = 0, 
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
                                     AAP = 0,
                                     AVP = 0,
                                     Ain = 0)
  )
  
  
  # Calculate parameters
  base.parm.c = BASE_PARAMS(expAGE, expBW, sex) # common parameters
  parm.c <- switch(exposure_type,
                   "Oral" = c(base.parm.c, 
                              expSTOP = expSTOP, 
                             expOral = exp, 
                              Tinput = Tinput, 
                              tinterval = tinterval),  
                   "Dermal" = c(DERMAL_PARAMS(expAGE, expBW, sex, base.parm.c), 
                                list(expSTOP = expSTOP, 
                                   expDermal = exp, 
                                     Tinput = Tinput, 
                                     tinterval = tinterval)),
                   "Oral_Dermal" = c(DERMAL_PARAMS(expAGE, expBW, sex = "M", base.parm.c), 
                                     list(expSTOP = expSTOP, 
                                         expOral = exp_Oral,
                                        expDermal = exp_Dermal, 
                                          Tinput = Tinput, 
                                          tinterval = tinterval)),
                   "Inhalation" = c(INHALATION_PARAMS(expAGE, expBW, sex = "M", base.parm.c), 
                                    list(expSTOP = expSTOP, 
                                      expLung = exp, 
                                         Tinput = Tinput, 
                                         tinterval = tinterval))
  )
  # write.csv(parm.c, file = here(OUTPUT, "ModelParameters.csv"), row.names = FALSE)
  
  
  
  ## Run the PBK model ####
  PBK_OUTPUT <- switch(exposure_type,
                        "Oral" = ORAL_PBK_RUN(y = Astate, 
                                              parms = parm.c, 
                                              times = seq(Tstart,Tstop,by=Dt)),
                        "Dermal" = DERMAL_PBK_RUN(y = Astate, 
                                                  parms = parm.c, 
                                                  times = seq(Tstart,Tstop,by=Dt)),
                        "Oral_Dermal" = ORAL_DERMAL_PBK_RUN(y = Astate, 
                                                            parms = parm.c, 
                                                            times = seq(Tstart, Tstop, by=Dt)),
                        "Inhalation" = INHALATION_PBK_RUN(y = Astate, 
                                                          parms = parm.c, 
                                                          times = seq(Tstart,Tstop,by=Dt))
  )

  return(list(
    CALC_Parameters = parm.c,
    PBK_OUTPUT = PBK_OUTPUT
  ))
  
  
}


# With physiological changes with lifestage (depending on age) ####
RUNandOUT_lifestage <- function(exposure_type, 
                      exp_Oral = NULL, exp_Dermal = NULL, exp, 
                      Tinput, tinterval, expSTOP, 
                      expAGE = NULL, expBW = NULL, sex = "M", 
                      Tstart, Tstop, Dt) {
  
  
  source(here("Script", "CALC_Parameters.R"), local = TRUE)
  source(here("Script", "PBK_model.R"), local = TRUE)  
  
  
  ## Add state ####
  Astate <- switch (exposure_type,
                    "Oral" = c(OD = 0,
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
                               Ain = 0),
                    "Dermal" = c(DD = 0,
                                 ASkB = 0,
                                 ASk = 0, 
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
                                 Ain = 0),
                    "Oral_Dermal" = c(OD = 0,
                                      DD = 0,
                                      ASkB = 0,
                                      ASk = 0, 
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
                                      Ain = 0),
                    "Inhalation" = c(LuD = 0,
                                     ALu = 0, 
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
                                     AAP = 0,
                                     AVP = 0,
                                     Ain = 0)
  )
  
  # Initial input for lifestage changes
  fullPBK_OUTPUT <- NULL
  newTstart <- Tstart
  newAstate <- Astate
  newAGE <- expAGE
  newBW <- expBW
  
  while (newTstart < Tstop) {
    
    yearlyTstop <- min(newTstart + 365, Tstop) 
    
    
    ## Calculare parameters, per newAGE ####
    
    base.parm.c = BASE_PARAMS(newAGE, newBW, sex) # common parameters
    parm.c <- switch(exposure_type,
                     "Oral" = c(base.parm.c, 
                                expSTOP = expSTOP, 
                               expOral = exp, 
                                Tinput = Tinput, 
                                tinterval = tinterval),  
                     "Dermal" = c(DERMAL_PARAMS(newAGE, newBW, sex, base.parm.c), 
                                  list(expSTOP = expSTOP, 
                                     expDermal = exp, 
                                       Tinput = Tinput, 
                                       tinterval = tinterval)),
                     "Oral_Dermal" = c(DERMAL_PARAMS(newAGE, newBW, sex, base.parm.c), 
                                       list(expSTOP = expSTOP, 
                                           expOral = exp_Oral,
                                          expDermal = exp_Dermal, 
                                            Tinput = Tinput, 
                                            tinterval = tinterval)),
                     "Inhalation" = c(INHALATION_PARAMS(newAGE, newBW, sex, base.parm.c), 
                                      list(expSTOP = expSTOP, 
                                        expLung = exp, 
                                           Tinput = Tinput, 
                                           tinterval = tinterval))
    )
    # write.csv(parm.c, file = here(OUTPUT, "ModelParameters.csv"), row.names = FALSE)
    
    
    
    ## Run the PBK model, per current age ####
    yearlyPBK_OUTPUT <- switch(exposure_type,
                               "Oral" = ORAL_PBK_RUN(y = newAstate, 
                                                     parms = parm.c, 
                                                     times = seq(newTstart, yearlyTstop, by=Dt)),
                               "Dermal" = DERMAL_PBK_RUN(y = newAstate, 
                                                         parms = parm.c, 
                                                         times = seq(newTstart, yearlyTstop, by=Dt)),
                               "Oral_Dermal" = ORAL_DERMAL_PBK_RUN(y = newAstate, 
                                                                   parms = parm.c, 
                                                                   times = seq(newTstart, yearlyTstop, by=Dt)),
                               "Inhalation" = INHALATION_PBK_RUN(y = newAstate, 
                                                                 parms = parm.c, 
                                                                 times = seq(newTstart, yearlyTstop, by=Dt))
    )
    
    if(is.null(fullPBK_OUTPUT)){
      fullPBK_OUTPUT <- yearlyPBK_OUTPUT # this is for the first iteration
    } else {
      fullPBK_OUTPUT <- rbind(fullPBK_OUTPUT, yearlyPBK_OUTPUT[-1,]) # binding all iteration results 
    }
    
    
    # Update input for next iteration
    newAstate <- switch(exposure_type, # Last state values
                        "Oral" = unlist(yearlyPBK_OUTPUT[nrow(yearlyPBK_OUTPUT), 2:16]),  
                        "Dermal" = unlist(yearlyPBK_OUTPUT[nrow(yearlyPBK_OUTPUT), 2:18]),
                        "Oral_Dermal" = unlist(yearlyPBK_OUTPUT[nrow(yearlyPBK_OUTPUT), 2:19]),
                        "Inhalation" = unlist(yearlyPBK_OUTPUT[nrow(yearlyPBK_OUTPUT), 2:18])
    )    
    newTstart <- yearlyTstop
    newAGE <- expAGE + (newTstart/365)
    newBW <- NA
    
  }

  PBK_OUTPUT <- fullPBK_OUTPUT
  
  return(list(
    CALC_Parameters = parm.c,
    PBK_OUTPUT = PBK_OUTPUT # time is in days
  ))
  
  
}

# Post-run analysis and plots for 1 individual ####

POST.Run <- function(exposure_type, 
                     exp_Oral = NULL, exp_Dermal = NULL, exp, 
                     Tinput, tinterval, expSTOP, 
                     expAGE = NULL, expBW = NULL, sex = "M", 
                     Tstart, Tstop, Dt,
                     RawData){
  
  ## Input ####
  PBK_OUTPUT <- RawData$PBK_OUTPUT
  Lit.HalfLifes <- read_csv(here("Input", "HalfLifes.csv"))
  
  ## Plot Mass Balance/Error ####
  MB.df <- PBK_OUTPUT %>% select(time, Ain, Atot, MB)  # change days to years if needed
  MB.df$MB <- round(MB.df$MB, 10) # rounding significance points
  MB.df$ERROR <- (MB.df$Ain - MB.df$Atot) / MB.df$Atot * 100
  MB.df$ERROR <- round(MB.df$ERROR, 10) # rounding significance points
  MB_plot <- ggplot(data = MB.df)+
    geom_line(aes(x = time, y = ERROR, color = "ERROR")) +
    geom_line(aes(x = time, y = MB, color = "MB")) +
    scale_color_manual(values = c("ERROR" = "blue", "MB" = "black"), 
                       name = NULL) +
    # ylim(0,1) +
    CP_theme +
    # theme(axis.text = element_text(size = 10),
    #       axis.title = element_text(size = 12)
    # )+
    ylab("MB / ERROR")
  MB_plot
  
  ## Plot organ concentrations ####
  # Create df with organ concentrations
  C_organs.df <- switch (exposure_type,
                         "Oral" = PBK_OUTPUT %>% 
                           transmute(
                             time = time,
                             CI = CI + CIL,
                             CL = CL_ec + CL_ic,
                             CK = CPTT + CPTL + CRKT + CRKL,
                             CA, CR, CP  
                           ) %>% 
                           rename("Intestine" = CI, 
                                  "Liver" = CL, 
                                  "Kidney" = CK, 
                                  "Adipose" = CA, 
                                  "Rest" = CR, 
                                  "Plasma" = CP) %>% 
                           pivot_longer(cols = Intestine:Plasma, names_to = "Organ", values_to = "Concentration"),
                         "Dermal" = PBK_OUTPUT %>% 
                           transmute(
                             time = time,
                             CI = CI + CIL,
                             CL = CL_ec + CL_ic,
                             CK = CPTT + CPTL + CRKT + CRKL,
                             CSk, CA, CR, CP  
                           ) %>% 
                           rename("Intestine" = CI, 
                                  "Liver" = CL, 
                                  "Kidney" = CK, 
                                  "Adipose" = CA, 
                                  "Rest" = CR,
                                  "Skin" = CSk,
                                  "Plasma" = CP) %>% 
                           pivot_longer(cols = Intestine:Plasma, names_to = "Organ", values_to = "Concentration"), 
                         "Oral_Dermal" = PBK_OUTPUT %>% 
                           transmute(
                             time = time,
                             CI = CI + CIL,
                             CL = CL_ec + CL_ic,
                             CK = CPTT + CPTL + CRKT + CRKL,
                             CSk, CA, CR, CP  
                           ) %>% 
                           rename("Intestine" = CI, 
                                  "Liver" = CL, 
                                  "Kidney" = CK, 
                                  "Adipose" = CA, 
                                  "Rest" = CR,
                                  "Skin" = CSk,
                                  "Plasma" = CP) %>% 
                           pivot_longer(cols = Intestine:Plasma, names_to = "Organ", values_to = "Concentration"), 
                         "Inhalation" = PBK_OUTPUT %>% 
                           transmute(
                             time = time,
                             CI = CI + CIL,
                             CL = CL_ec + CL_ic,
                             CK = CPTT + CPTL + CRKT + CRKL,
                             CLu, CA, CR, CP  
                           ) %>% 
                           rename("Intestine" = CI, 
                                  "Liver" = CL, 
                                  "Kidney" = CK, 
                                  "Adipose" = CA, 
                                  "Rest" = CR, 
                                  "Lungs" = CLu,
                                  "Plasma" = CP) %>% 
                           pivot_longer(cols = Intestine:Plasma, names_to = "Organ", values_to = "Concentration")
  ) %>% mutate(time = time/365) #transforming time in years
  
  
  # Plot organ concentrations
  Plot_C_organs <- C_organs.df %>% 
    ggplot(aes(time, Concentration)) +
    geom_path(linewidth = 0.5) +
    facet_wrap(~Organ) +
    labs(title = "PFOA organ concentrations",
         x = "Time (years)", # check that time is indeed in days and not years
         y = "Concentration (ng/ml)") +
    CP_theme #+
  # theme(
  #   axis.text = element_text(size = 10),
  #   axis.title = element_text(size = 12)
  # )
  Plot_C_organs
  
  ## Calculate AUC and Half life ####
  AUC <- trapz(PBK_OUTPUT[ , "time"], PBK_OUTPUT[ , "CP"])  # ug*day/L
  print(AUC)
  
  time <- PBK_OUTPUT[ , "time"] # days
  conc <- PBK_OUTPUT[ , "CP"] # ug/L or ng/ml
  Cmax <- max(conc)
  Tmax <- time[which.max(conc)]
  tlast <- max(time[conc > 0])
  half_life <- pk.calc.half.life(
    conc,
    time,
    Tmax,
    tlast
  )
  HalfLife <- half_life$half.life/365  # half-life in years
  print(HalfLife)
  
  # Evaluate against HBM data ####
  HL_literature <- Lit.HalfLifes %>%
    filter(species == "human", chemical == "pfoa", parameter == "HalfLife") %>%
    select(c(value_average, n, sex)) %>%
    rename(HL_observed = value_average) %>%
    mutate(HL_observed = as.numeric(HL_observed),  # years
           n = as.numeric(n)) %>% 
    filter(sex %in% c("F", "M")) %>% 
    mutate(value = case_when(sex == "F" ~ 0.75,  
                             sex == "M" ~ 1.5),
           Origin = "Observed") 
  
  
  HL_predicted <- data.frame(
    HL_predicted = HalfLife, 
    sex = sex,
    value = ifelse(sex == "F", 0.75, 1.5),
    Origin = "Predicted")
  
  HL_violin_plot <- 
    ggplot() +
    geom_violin(data = HL_literature, 
                aes(x = 1.5, y = HL_observed), 
                fill = "grey89", color = NA, width = 0.5, trim = FALSE) +
    geom_point(data = HL_literature,
               aes(x = 1.5, y = HL_observed, size = n),
               color = "grey70", alpha = 0.5, position = position_jitter(width = 0.05)) +
    geom_point(data = HL_predicted,
               aes(x = 1.5, y = HL_predicted),
               shape = 18, size = 5,  alpha = 0.9, position = position_jitter(width = 0.01)) +
    
    scale_size_continuous(range = c(1, 5)) +
    
    labs(title = "Predicted over observed half lives",
         x = "", y = "Half life (years)") +
    guides(size = guide_legend(title = "HBM sample size")) +
    CP_theme +
    theme(
      axis.text.x = element_blank())
  HL_violin_plot
  
  message("Post run analysis finished")
  
  return(list(
    ANALYSED_data = data.frame(exposure_type = paste(exposure_type),
                               exp = paste(round(exp, digits = 10),"ug/kg/day", sep = "_"),
                               expAGE = paste(expAGE, "years", sep = "_"),
                               Tstart = paste(round(Tstart/365, digits = 10), "years", sep = "_"),
                               expSTOP = paste(round(expSTOP/365, digits = 10), "years", sep = "_"),
                               Tstop = paste(round(Tstop/365, digits = 10), "years", sep = "_"),
                               expBW = paste(round(expBW, digits = 10), "kg", sep = "_"),
                               sex = paste(sex),
                               AUC = paste(round(AUC, digits = 10), "ug*day/L", sep = "_"),
                               HalfLife = paste(round(HalfLife, digits = 10), "years", sep = "_")),
    OUT_Plots = list(MB_plot, Plot_C_organs, HL_violin_plot) # could be removed if it's too heavy for R
  ))
}


# Post-run analysis and plots per subject ####

Pers.POST.Run <- function(exposure_type, 
                     exp_Oral = NULL, exp_Dermal = NULL, exp, 
                     Tinput, tinterval, expSTOP, 
                     expAGE = NULL, expBW = NULL, sex = "M", 
                     Tstart, Tstop, Dt,
                     RawData){
  
  ## Input ####
  PBK_OUTPUT <- RawData #$PBK_OUTPUT
  Lit.HalfLifes <- read_csv(here("Input", "HalfLifes.csv"))
  
  ## Plot Mass Balance/Error ####
  MB.df <- PBK_OUTPUT %>% select(time, Ain, Atot, MB)  # change days to years if needed
  MB.df$MB <- round(MB.df$MB, 10) # rounding significance points
  MB.df$ERROR <- (MB.df$Ain - MB.df$Atot) / MB.df$Atot * 100
  MB.df$ERROR <- round(MB.df$ERROR, 10) # rounding significance points
  MB_plot <- ggplot(data = MB.df)+
    geom_line(aes(x = time, y = ERROR, color = "ERROR")) +
    geom_line(aes(x = time, y = MB, color = "MB")) +
    scale_color_manual(values = c("ERROR" = "blue", "MB" = "black"), 
                       name = NULL) +
    # ylim(0,1) +
    CP_theme +
    # theme(axis.text = element_text(size = 10),
    #       axis.title = element_text(size = 12)
    # )+
    ylab("MB / ERROR")
  MB_plot
  
  ## Plot organ concentrations ####
  # Create df with organ concentrations
  C_organs.df <- switch (exposure_type,
                         "Oral" = PBK_OUTPUT %>% 
                           transmute(
                             time = time,
                             CI = CI + CIL,
                             CL = CL_ec + CL_ic,
                             CK = CPTT + CPTL + CRKT + CRKL,
                             CA, CR, CP  
                           ) %>% 
                           rename("Intestine" = CI, 
                                  "Liver" = CL, 
                                  "Kidney" = CK, 
                                  "Adipose" = CA, 
                                  "Rest" = CR, 
                                  "Plasma" = CP) %>% 
                           pivot_longer(cols = Intestine:Plasma, names_to = "Organ", values_to = "Concentration"),
                         "Dermal" = PBK_OUTPUT %>% 
                           transmute(
                             time = time,
                             CI = CI + CIL,
                             CL = CL_ec + CL_ic,
                             CK = CPTT + CPTL + CRKT + CRKL,
                             CSk, CA, CR, CP  
                           ) %>% 
                           rename("Intestine" = CI, 
                                  "Liver" = CL, 
                                  "Kidney" = CK, 
                                  "Adipose" = CA, 
                                  "Rest" = CR,
                                  "Skin" = CSk,
                                  "Plasma" = CP) %>% 
                           pivot_longer(cols = Intestine:Plasma, names_to = "Organ", values_to = "Concentration"), 
                         "Oral_Dermal" = PBK_OUTPUT %>% 
                           transmute(
                             time = time,
                             CI = CI + CIL,
                             CL = CL_ec + CL_ic,
                             CK = CPTT + CPTL + CRKT + CRKL,
                             CSk, CA, CR, CP  
                           ) %>% 
                           rename("Intestine" = CI, 
                                  "Liver" = CL, 
                                  "Kidney" = CK, 
                                  "Adipose" = CA, 
                                  "Rest" = CR,
                                  "Skin" = CSk,
                                  "Plasma" = CP) %>% 
                           pivot_longer(cols = Intestine:Plasma, names_to = "Organ", values_to = "Concentration"), 
                         "Inhalation" = PBK_OUTPUT %>% 
                           transmute(
                             time = time,
                             CI = CI + CIL,
                             CL = CL_ec + CL_ic,
                             CK = CPTT + CPTL + CRKT + CRKL,
                             CLu, CA, CR, CP  
                           ) %>% 
                           rename("Intestine" = CI, 
                                  "Liver" = CL, 
                                  "Kidney" = CK, 
                                  "Adipose" = CA, 
                                  "Rest" = CR, 
                                  "Lungs" = CLu,
                                  "Plasma" = CP) %>% 
                           pivot_longer(cols = Intestine:Plasma, names_to = "Organ", values_to = "Concentration")
  ) %>% mutate(time = time/365) #transforming time in years
  
  
  # Plot organ concentrations
  Plot_C_organs <- C_organs.df %>% 
    ggplot(aes(time, Concentration)) +
    geom_path(linewidth = 0.5) +
    facet_wrap(~Organ) +
    labs(title = "PFOA organ concentrations",
         x = "Time (years)", # check that time is indeed in days and not years
         y = "Concentration (ng/ml)") +
    CP_theme #+
    # theme(
    #   axis.text = element_text(size = 10),
    #   axis.title = element_text(size = 12)
    # )
  Plot_C_organs

  ## Calculate AUC and Half life ####
  AUC <- trapz(PBK_OUTPUT[ , "time"], PBK_OUTPUT[ , "CP"])  # ug*day/L
  print(AUC)
  
  time <- PBK_OUTPUT[ , "time"] # days
  conc <- PBK_OUTPUT[ , "CP"] # ug/L or ng/ml
  Cmax <- max(conc)
  Tmax <- time[which.max(conc)]
  tlast <- max(time[conc > 0])
  half_life <- pk.calc.half.life(
    conc,
    time,
    Tmax,
    tlast
  )
  HalfLife <- half_life$half.life/365  # half-life in years
  print(HalfLife)
  
  # Evaluate against HBM data ####
  HL_literature <- Lit.HalfLifes %>%
    filter(species == "human", chemical == "pfoa", parameter == "HalfLife") %>%
    select(c(value_average, n, sex)) %>%
    rename(HL_observed = value_average) %>%
    mutate(HL_observed = as.numeric(HL_observed),  # years
           n = as.numeric(n)) %>% 
    filter(sex %in% c("F", "M")) %>% 
    mutate(value = case_when(sex == "F" ~ 0.75,  
                             sex == "M" ~ 1.5),
           Origin = "Observed") 
  
  
  HL_predicted <- data.frame(
    HL_predicted = HalfLife, 
    sex = sex,
    value = ifelse(sex == "F", 0.75, 1.5),
    Origin = "Predicted")
  
  HL_violin_plot <- 
    ggplot() +
    geom_violin(data = HL_literature, 
                aes(x = 1.5, y = HL_observed), 
                fill = "grey89", color = NA, width = 0.5, trim = FALSE) +
    geom_point(data = HL_literature,
               aes(x = 1.5, y = HL_observed, size = n),
               color = "grey70", alpha = 0.5, position = position_jitter(width = 0.05)) +
    geom_point(data = HL_predicted,
               aes(x = 1.5, y = HL_predicted),
               shape = 18, size = 5,  alpha = 0.9, position = position_jitter(width = 0.01)) +
    
    scale_size_continuous(range = c(1, 5)) +
    
    labs(title = "Predicted over observed half lives",
         x = "", y = "Half life (years)") +
    guides(size = guide_legend(title = "HBM sample size")) +
    CP_theme +
    theme(
          axis.text.x = element_blank())
  HL_violin_plot
  
  message("Post run analysis finished")
  
  return(list(
    ANALYSED_data = data.frame(exposure_type = paste(exposure_type),
                               exp = paste(round(exp, digits = 10),"ug/kg/day", sep = "_"),
                               expAGE = paste(expAGE, "years", sep = "_"),
                               Tstart = paste(round(Tstart/365, digits = 10), "years", sep = "_"),
                               expSTOP = paste(round(expSTOP/365, digits = 10), "years", sep = "_"),
                               Tstop = paste(round(Tstop/365, digits = 10), "years", sep = "_"),
                               expBW = paste(round(expBW, digits = 10), "kg", sep = "_"),
                               sex = paste(sex),
                               AUC = paste(round(AUC, digits = 10), "ug*day/L", sep = "_"),
                               HalfLife = paste(round(HalfLife, digits = 10), "years", sep = "_")),
    OUT_Plots = list(MB_plot, Plot_C_organs, HL_violin_plot) # could be removed if it's too heavy for R
                               ))
}


# Post-run analysis and plots per population ####

Pop.POST.Run <- function(Input,
                         RawData,
                         Pers.POST.RUN_RESULTS){
  
  
  ## Input ####
  
  PBK_OUT <- RawData$PBK_OUTPUT # PBK model results per subject
  ANALYSED_data <- lapply(Pers.POST.RUN_RESULTS, function(x) x$ANALYSED_data)
  Lit.HalfLifes <- read_csv(here("Input", "HalfLifes.csv"))
  Parameters <- RawData$CALC_Parameters
  BW <- sapply(Parameters, function(x) x$BW)
  GFR <- sapply(Parameters, function(x) x$GFR)
  QKc <- sapply(Parameters, function(x) x$QKc)
  QC <- sapply(Parameters, function(x) x$QC)
  GFR.ff <- 0.18*QKc*QC
  
  
  PredictedObserved <- data.frame(
    Idcode = Input$Idcode,
    exp = Input$exp,
    exp_Oral = Input$exp_Oral,
    exp_Dermal = Input$exp_Dermal,
    expSTOP = Input$expSTOP, # this is the time of the stop of exposure in days
    CP_observed = Input$CP_measured, # measured plasma PFOA concentration
    samplingT = Input$samplingT, # time at which the plasma concentration was measured
    HL_observed = Input$HL_observed, # half life reported in the HBM study
    sex = Input$sex,
    expAGE = Input$expAGE,
    BW = BW,
    GFR = GFR,
    GFR.ff = GFR.ff
  ) 
  PredictedObserved <- PredictedObserved %>% 
    mutate(log_exp = log10(exp),
           log_exp_Oral = log10(exp_Oral),
           log_exp_Dermal = log10(exp_Dermal))
  
  # Add predicted concentration at the sampling time 
  PredictedObserved <-  PredictedObserved %>% 
    mutate(
      CP_predicted = map2_dbl( 
        PBK_OUT,samplingT, ~ {
          idx <- which.min(abs(.x$time - .y))
          .x$CP[idx]
        }
      )
    ) 
  
  # Add predicted half life
  PredictedObserved <- PredictedObserved %>% mutate(
    Idcode = seq_along(ANALYSED_data),
    HalfLife = sapply(ANALYSED_data, function(x) x$HalfLife)) %>%
    separate(col = HalfLife, into = c("HL_predicted", "unit"), sep = "_") %>%
    mutate(HL_predicted = as.numeric(HL_predicted), 
           HL_predicted = round(HL_predicted,1))
  
  
  ## Plots ####
  
  hist(PredictedObserved$exp)
  hist(PredictedObserved$exp_Oral)
  hist(PredictedObserved$exp_Dermal)
  hist(PredictedObserved$log_exp)
  hist(PredictedObserved$log_exp_Oral)
  hist(PredictedObserved$log_exp_Dermal)
  
  hist(PredictedObserved$CP_observed)
  hist(PredictedObserved$HL_observed)
  
  
  ### Exposure estimate vs predicted plasma concentration ####
  Plot_exp_vs_P_CP <- PredictedObserved %>% 
    ggplot(aes(exp, CP_predicted)) +
    geom_point(color = "black", size = 0.5) +
    CP_theme +
    labs(title="Exposure vs predicted plasma concentration",
         x="\n Exposure (\u03BCg/kg bw/dayL)", 
         y="Predicted concentration (\u03BCg/L)\n") 
  Plot_exp_vs_P_CP
  
  ### Exposure age  vs predicted plasma concentration ####
  Plot_age_vs_P_CP <- PredictedObserved %>% 
    ggplot(aes(expAGE, CP_predicted)) +
    geom_point(color = "black", size = 0.5) +
    CP_theme +
    labs(title="Age at the start of exposure vs predicted plasma concentration",
         x="\n Exposure age (years)", 
         y="Predicted concentration (\u03BCg/L)\n") 
  Plot_age_vs_P_CP
  
  ### GFR  vs predicted plasma concentration ####
  Plot_GFR_vs_P_CP <- PredictedObserved %>% 
    ggplot() +
    geom_point(aes(GFR, CP_predicted), color = "black", size = 0.5) +
    geom_point(aes(GFR.ff, CP_predicted), color = "blue", size = 0.5) +
    CP_theme +
    labs(title="Glomerular filtration rate vs predicted plasma concentration",
         x="\n GFR (L/d)", 
         y="Predicted concentration (\u03BCg/L)\n") 
  Plot_GFR_vs_P_CP
  
  ### Exposure estimate vs predicted plasma concentration ####
  Plot_exp_vs_P_HL <- PredictedObserved %>% 
    ggplot(aes(exp, HL_predicted)) +
    geom_point(color = "black", size = 0.5) +
    CP_theme +
    labs(title="Exposure vs predicted half life",
         x="\n Exposure (\u03BCg/kg bw/dayL)", 
         y="Predicted half life (years)") 
  Plot_exp_vs_P_HL
  
  ### Exposure age vs predicted half life ####
  Plot_age_vs_P_HL <- PredictedObserved %>% 
    ggplot(aes(expAGE, HL_predicted)) +
    geom_point(color = "black", size = 0.5) +
    CP_theme +
    labs(title="Age at the start of exposure vs predicted half life",
         x="\n Exposure age (years)", 
         y="Predicted half life (years)") 
  Plot_age_vs_P_HL
  
  ### GFRvs predicted plasma concentration ####
  Plot_GFR_vs_P_HL <- PredictedObserved %>% 
    ggplot() +
    geom_point(aes(GFR, HL_predicted), color = "black", size = 0.5) +
    geom_point(aes(GFR.ff, HL_predicted), color = "blue", size = 0.5) +
    CP_theme +
    labs(title="Glomerular filtration rate vs predicted half life",
         x="\n GFR (L/d)", 
         y="Predicted half life (years)") 
  Plot_GFR_vs_P_HL
  
  
  ### Observed vs Predicted plasma concentrations ####
  CP_Regression <- lm(CP_predicted~CP_observed, data=PredictedObserved) #lm(y~x) (y is the dependent variable)
  summary(CP_Regression)
  
  CP_Regression_plot <- CP_Regression %>%
    ggplot(aes(log(CP_predicted), log(CP_observed))) +
    # geom_smooth(method='lm', color = "black", se = TRUE) +
    geom_abline(intercept = 0, slope = 1, linetype = "solid", linewidth = 0.3, color = "grey50") +  
    geom_abline(intercept = log(2), slope = 1, linetype = "dashed", linewidth = 0.2, color = "grey50") + #2 fold
    geom_abline(intercept = log(0.5), slope = 1, linetype = "dashed", linewidth = 0.2, color = "grey50") + #2 fold
    geom_abline(intercept = log(3), slope = 1, linetype = "dotted", linewidth = 0.1, color = "grey50") +  #3 fold
    geom_abline(intercept = log(1/3), slope = 1, linetype = "dotted", linewidth = 0.1, color = "grey50") + #3 fold
    geom_point(color = "black", size = 0.2) +
    CP_theme +
    labs(title="Observed vs predicted plasma concentrations",
         x="\n log10 (Predicted concentration) (\u03BCg/L)", 
         y=" log10 (Observed concentration) (\u03BCg/L)\n") 
  CP_Regression_plot
  
  
  ### Observed vs Predicted half lives ####
  HL_Regression <- lm(HL_predicted~HL_observed, data=PredictedObserved)
  summary(HL_Regression)
  
  HL_Regression_plot <- HL_Regression %>%
    ggplot(aes(log10(HL_predicted), log10(HL_observed))) +
    # geom_smooth(method='lm', color = "black", se = TRUE) +
    geom_abline(intercept = 0, slope = 1, linetype = "solid", linewidth = 0.3, color = "grey50") +  
    geom_abline(intercept = log10(2), slope = 1, linetype = "dashed", linewidth = 0.2, color = "grey50") +  #2 fold
    geom_abline(intercept = log10(0.5), slope = 1, linetype = "dashed", linewidth = 0.2, color = "grey50") + #2 fold
    geom_abline(intercept = log10(3), slope = 1, linetype = "dotted", linewidth = 0.1, color = "grey50") +  #3 fold
    geom_abline(intercept = log10(1/3), slope = 1, linetype = "dotted", linewidth = 0.1, color = "grey50") + #3 fold
    geom_point(color = "black", size = 0.2) +
    CP_theme +
    labs(title="Observed vs predicted half lives",
         x = "log10 (Predicted Half life) (years)", 
         y = "log10 (Observed Half life) (years)")
  HL_Regression_plot
  
  
  ### Plot predicted half lives over those reported in literature ####
  HL_literature <- Lit.HalfLifes %>%
    filter(species == "human", chemical == "pfoa", parameter == "HalfLife") %>%
    select(c(value_average, n, sex)) %>%
    rename(HL_observed = value_average) %>%
    mutate(HL_observed = as.numeric(HL_observed),  # years
           n = as.numeric(n)) %>% 
    filter(sex %in% c("F", "M")) %>% 
    mutate(value = case_when(sex == "F" ~ 0.75,  
                             sex == "M" ~ 1.5),
           Origin = "Observed") 
  
  
  HL_predicted <- data.frame(
    HL_predicted = PredictedObserved$HL_predicted, 
    sex = PredictedObserved$sex) %>% 
    mutate(value = case_when(sex == "F" ~ 0.75,  
                             sex == "M" ~ 1.5),
           Origin = "Predicted")
  
  HL_violin_plot <- 
    ggplot() +
    geom_violin(data = filter(HL_literature, sex == "M"), 
                aes(x = 1.5, y = HL_observed), 
                fill = "grey89", color = NA, width = 0.5, trim = FALSE) +
    geom_point(data = filter(HL_literature, sex == "M"),
               aes(x = 1.5, y = HL_observed, size = n),
               color = "grey70", alpha = 0.5, position = position_jitter(width = 0.05)) +
    geom_point(data = filter(HL_predicted, sex == "M"),
               aes(x = 1.5, y = HL_predicted),
               shape = 18, size = 5,  alpha = 0.9, position = position_jitter(width = 0.01)) +
    
    geom_violin(data = filter(HL_literature, sex == "F"), 
                aes(x = 0.75, y = HL_observed), 
                fill = "grey89", color = NA, width = 0.5, trim = FALSE) +
    geom_point(data = filter(HL_literature, sex == "F"),
               aes(x = 0.75, y = HL_observed, size = n),
               color = "grey70", alpha = 0.5, position = position_jitter(width = 0.05)) +
    geom_point(data = filter(HL_predicted, sex == "F"),
               aes(x = 0.75, y = HL_predicted),
               shape = 18, size = 5,  alpha = 0.9, position = position_jitter(width = 0.01)) +
    scale_x_continuous(breaks = c(0.75, 1.5),       
                       labels = c("Female", "Male")) + 
    scale_size_continuous(range = c(1, 5)) +
    
    labs(title = "Predicted half lifes plotted over the range of observed",
         x = "", y = "Half life (years)") +
    guides(size = guide_legend(title = "HBM sample size")) +
    CP_theme +
    theme(legend.position = "right")
  HL_violin_plot
  
  message("Post run analysis finished")
  
  return(list(Plot_exp_vs_P_CP,
              Plot_age_vs_P_CP,
              Plot_GFR_vs_P_CP,
              Plot_exp_vs_P_HL,
              Plot_age_vs_P_HL,
              Plot_GFR_vs_P_HL,
              CP_Regression_plot,
              HL_Regression_plot,
              HL_violin_plot))
  

  }