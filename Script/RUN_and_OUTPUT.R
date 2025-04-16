# --------------------------------------------------------------------------- #
# SCRIPT FOR CALCULATING MODEL RESULTS
# By: Chrysanthi Pachoulide
# Date: 14-04-2025
# --------------------------------------------------------------------------- #


RUNandOUT <- function(exposure_type, 
                      expCONC, expCONC_Oral = NULL, expCONC_Dermal = NULL, 
                      Tinput, tinterval, expSTOP, 
                      expAGE = NULL, expBW = NULL, sex = "M", 
                      Tstart, Tstop, Dt) {
  
  
  # Load R file to calculate parameters
  source(here("Script", "CALC_Parameters.R"), local = TRUE)
  
  
  # Calculate parameters
  base.parm.c = BASE_PARAMS(expAGE, expBW, sex) # common parameters
  parm.c <- switch(exposure_type,
                   "Oral" = c(base.parm.c, 
                                 expSTOP = expSTOP, 
                                 COral = expCONC, 
                                 Tinput = Tinput, 
                                 tinterval = tinterval),  
                   "Dermal" = c(DERMAL_PARAMS(expAGE, expBW, sex, base.parm.c), 
                                list(expSTOP = expSTOP, 
                                     CDermal = expCONC, 
                                     Tinput = Tinput, 
                                     tinterval = tinterval)),
                   "Oral_Dermal" = c(DERMAL_PARAMS(expAGE, expBW, sex = "M", base.parm.c), 
                                     list(expSTOP = expSTOP, 
                                          COral = expCONC_Oral,
                                          CDermal = expCONC_Dermal, 
                                          Tinput = Tinput, 
                                          tinterval = tinterval)),
                   "Inhalation" = c(INHALATION_PARAMS(expAGE, expBW, sex = "M", base.parm.c), 
                                    list(expSTOP = expSTOP, 
                                         CLung = expCONC, 
                                         Tinput = Tinput, 
                                         tinterval = tinterval))
  )
  
  write.csv(parm.c, file = here(OUTPUT, "ModelParameters.csv"), row.names = FALSE)
  
  
  # Add state
  A_init <- switch (exposure_type,
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
  
  
  # Load R file with the PBK model
  source(here("Script", "PBK_model.R"), local = TRUE)  # Load model functions
  
  
  # Run the PBK model
  output_PFOA <- switch(exposure_type,
                        "Oral" = ORAL_PBK_RUN(A_init, parm.c, TIME = seq(Tstart,Tstop,by=Dt)),
                        "Dermal" = DERMAL_PBK_RUN(A_init, parm.c, TIME = seq(Tstart,Tstop,by=Dt)),
                        "Oral_Dermal" = ORAL_DERMAL_PBK_RUN(A_init, parm.c, TIME = seq(Tstart, Tstop, by=Dt)),
                        "Inhalation" = INHALATION_PBK_RUN(A_init, parm.c, TIME = seq(Tstart,Tstop,by=Dt))
  )
  
  
  # Creating data frames for data-analysis
  output.df <- output_PFOA %>% mutate(time = time/365)  # time in years
  write.csv(output.df, file = here(OUTPUT, "PFOA_PBKoutput.csv"), row.names = FALSE)
  
  C_organs.df <- switch (exposure_type,
                         "Oral" = output.df %>% 
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
                         "Dermal" = output.df %>% 
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
                         "Oral_Dermal" = output.df %>% 
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
                         "Inhalation" = output.df %>% 
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
  )
  
  # Check Mass Balance/Error 
  MB.df <- output.df %>% select(time, Ain, Atot, MB)  # change days to years if needed
  MB.df$MB <- round(MB.df$MB, 10) # rounding significance points
  MB.df$ERROR <- (MB.df$Ain - MB.df$Atot) / MB.df$Atot * 100
  MB.df$ERROR <- round(MB.df$ERROR, 10) # rounding significance points
  MB_plot <- ggplot(data = MB.df)+
    geom_line(aes(x = time, y = ERROR, color = "ERROR")) +
    geom_line(aes(x = time, y = MB, color = "MB")) +
    scale_color_manual(values = c("ERROR" = "blue", "MB" = "black"), 
                       name = NULL) +
    # ylim(0,1) +
    theme_minimal() +
    theme(
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 12)
    )+
    ylab("MB / ERROR")
  MB_plot
  ggsave(filename = here(OUTPUT,"MB_ERROR.png"), 
         dpi = 300,
         width = 17,      
         height = 8,      
         units = "cm")
  
  
  # Plot organ concentrations
  Plot_C_organs <- C_organs.df %>% 
    ggplot(aes(time, Concentration)) +
    geom_path(linewidth = 0.5) +
    facet_wrap(~Organ) +
    labs(title = "PFOA organ concentrations",
         x = "Time (years)", # check that time is indeed in days and not years
         y = "Concentration (ng/ml)") +
    theme_minimal()+
    theme(
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 12)
    )
  Plot_C_organs
  ggsave(filename = here(OUTPUT, "Plot_C_organ.png"), 
         dpi = 300,
         width = 17,      
         height = 8,      
         units = "cm")
  
  Plot_C_plasma <- C_organs.df %>% 
    filter(Organ == "Plasma") %>% 
    ggplot(aes(time, Concentration)) +
    geom_path(linewidth = 0.5) +
    labs(title = "PFOA plasma concentration",
         x = "Time (years)", # check that time is indeed in days and not years
         y = "Concentration (ng/ml)") +
    theme_minimal()+
    theme(
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 12)
    )
  Plot_C_plasma
  ggsave(filename = here(OUTPUT, "Plot_C_plasma.png"), 
         dpi = 300,
         width = 17,      
         height = 8,      
         units = "cm")
  
  
  # Calculate AUC
  AUC <- trapz(output_PFOA[ , "time"], output_PFOA[ , "CP"])  # ug*day/L
  print(AUC)
  
  
  # Calculate Half life
  time <- output.df[ , "time"] # years
  conc <- output.df[ , "CP"] # ug/L or ng/ml
  Cmax <- max(conc)
  Tmax <- time[which.max(conc)]
  tlast <- max(time[conc > 0])
  half_life <- pk.calc.half.life(
    conc,
    time,
    Tmax,
    tlast
  )
  HalfLife <- half_life$half.life  # half-life in years
  print(HalfLife)
  
  # Plot observed vs predicted half-life
  ObsHalfLifes <- read_csv(here("Input", "HalfLifes.csv"))
  
  Observed.df <- ObsHalfLifes %>%
    filter(species == "human",
           chemical == "pfoa",
           parameter== "HalfLife") %>%
    select(c(value_average,n)) %>%
    rename(HalfLife = value_average) %>%
    mutate(value = 1,
           Origin = "Observed")
  Observed.df$HalfLife <- as.numeric(Observed.df$HalfLife) # years
  Observed.df$n <- as.numeric(Observed.df$n)
  
  Predicted.df <- data.frame(
    HalfLife = HalfLife,
    Origin = "Predicted",
    value = 1, n = 1)
  Observed.df <- data.frame(
    HalfLife = Observed.df$HalfLife,
    Origin = "Observed",
    value = 1,
    n = Observed.df$n)
  
  HalfLifes <- rbind(Predicted.df, Observed.df)
  
  range <- c(min(Observed.df$n), max(Observed.df$n))
  
  Plot_HalfLifes <- ggplot() +
    geom_violin(
      data = Observed.df,
      aes(value, HalfLife),
      color = "transparent",
      fill = "grey89") +
    geom_point(
      data = Observed.df,
      aes(value, HalfLife, size = n),  
      color = "black",
      alpha = 0.5,  
      shape = 20) +
    geom_point(
      data = Predicted.df,
      aes(value, HalfLife),
      color = "red",
      alpha = 0.7,
      size = 10,
      shape = 18) +
    labs(y = "Half life (years)") + 
    scale_size_continuous(range = c(1, 10), 
                          name = "Sample size") + 
    theme_minimal() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank(),
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 12),
      legend.position = "top"
    )
  Plot_HalfLifes
  ggsave(filename = here(OUTPUT, "ExpVsSimHalfLife.png"), 
         dpi = 300,
         width = 17,      
         height = 8,      
         units = "cm")
  
  return(list(
    Plot_C_organs,
    Plot_HalfLifes,
    data = output.df,
    AUC = AUC,
    HalfLife = HalfLife
  ))
  
  
}
