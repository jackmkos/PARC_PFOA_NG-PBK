# --------------------------------------------------------------------------- #
# SCRIPT FOR COMPARING SIMULATED AGAINST OBSERVED DATA
# By: Chrysanthi Pachoulide
# Note: This file contains the codes used to evaluate the model against observed data and to produce the relevant plots
# Date: 06-05-2025
# --------------------------------------------------------------------------- #

rm(list=ls()) # to clear out the global environment

# Packages
library(here)
library(tidyverse)
library(deSolve)
library(PKNCA)
library(pracma)
library(showtext)
font_add(family = "Garamond", regular = "GARA.TTF")
showtext_auto()

# Set output storage directory
OUTPUT <- here("Output", format(Sys.Date(), "%Y-%m-%d"), format(Sys.time(), "%H-%M-%S"))
dir.create(OUTPUT, recursive = TRUE)

CP_theme <- theme_minimal() +
  theme(
    axis.text = element_text(size = 62),
    axis.title = element_text(size = 70),
    plot.title = element_text(size = 65, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 62, hjust = 0.5),
    # legend.position = "none",
    plot.margin = margin(0.2, 0.2, 0.2, 0.2, "cm") 
  )

# 1. Evaluate against Abraham data and EFSA-2020 model####

# Input data
Abraham_data <- read.csv(here("Input", "AbrahamData.csv"))
load(here("Output", "Predicted.AbrahamSettings", "RESULTS_PFOA_PBK.RData"))
PredEFSA <- read.csv(here("Output", "Validation.Abraham.EFSA_2020", "RESULTS_EFSA 2020.csv"))

## Plot serum ####
Abraham_serum <- Abraham_data %>% filter(Matrix == "Serum")
ObsSerum <- Abraham_serum %>%
  rename(time = Time.days) %>%
  rename(CP = MPFOA) %>%    # ng/ml
  mutate(CP = CP - 0.130) %>% # substracting the pre-existing level of 0.130ng/ml from their previous study, as also done in the ref. article: https://doi.org/10.1016/j.envint.2024.109047 (table 3)
  mutate(Origin = "Observed") %>% 
  select(c(time, CP))

Predicted_serum <- RESULTS$OUT_RAW_data %>% select(time, CP)
PredSerum <- Predicted_serum %>%
  select(time, CP) %>%
  mutate(time = time*365) %>%
  filter(time <=450) %>% 
  mutate(Origin = "Predicted")
PredictedEFSA_serum <- PredEFSA %>% select(time, CA_PFOA) %>% 
  rename(CP = CA_PFOA) %>% 
  mutate(Origin = "PredictedEFSA")
PredEFSASerum <- PredictedEFSA_serum %>% 
  filter(time <=450) 

Plot_Serum <- ggplot() +
  geom_path(data = PredSerum, aes(x = time, y = CP, color = "Predicted"), linewidth = 0.3) +
  geom_point(data = ObsSerum, aes(x = time, y = CP, color = "Observed"), size = 0.1) +
  geom_path(data = PredEFSASerum, aes(x = time, y = CP, color = "Predicted EFSA 2020"), 
    linetype = "dashed", linewidth = 0.3) +
  scale_color_manual(
    name = "Legend",
    values = c("Predicted" = "blueviolet", 
               "Observed" = "black", 
               "Predicted EFSA 2020" = "darkblue")) +
  CP_theme +
  ylab("PFOA Serum Concentration (ng/ml)") +
  xlab("Time (days)") +
  theme(legend.position = c(0.9, 0.9),
        legend.text = element_text(size = 52),
        legend.title = element_blank(),
        legend.key.size = unit(0.5, "lines"),
        legend.spacing.y = unit(0, "cm"),
        legend.margin = margin(0, 0, 0, 0))
Plot_Serum
ggsave(
  filename = here(OUTPUT, "AbrahamVsPred_Serum.png"),
  plot = Plot_Serum,
  dpi = 1000,
  width = 9,
  height = 7,  
  units = "cm"
)

## Plot All Organs ####

Predicted_conc <- RESULTS$OUT_RAW_data %>% select(time, CP, CI, CPTT, CPTL, CRKT, CRKL, CL_ec, CL_ic, CA) %>%
  mutate(time = time*365,
         CK = CPTT + CPTL + CRKT + CRKL, 
         CL = CL_ic + CL_ec) %>%
  rename(Plasma = CP, Intestine = CI, Kidney = CK, Liver = CL, Adipose = CA) %>% 
  filter(time <=450) %>% 
  mutate(Origin = "Predicted") %>% 
  select(time, Plasma, Intestine, Kidney, Liver, Adipose, Origin)

PredictedEFSA_conc <- PredEFSA %>% select(time, CA_PFOA, CK_PFOA, CG_PFOA, CL_PFOA, CF_PFOA) %>% 
  rename(Plasma = CA_PFOA, Intestine = CG_PFOA, Kidney = CK_PFOA, Liver = CL_PFOA, Adipose = CF_PFOA) %>% 
  mutate(Origin = "PredictedEFSA")%>% 
  filter(time <=450) %>% 
  select(time, Plasma, Intestine, Kidney, Liver, Adipose, Origin)

ComparePredicted <- rbind(Predicted_conc, PredictedEFSA_conc) %>% 
  pivot_longer(values_to = "Concentration", names_to = "Organ", Plasma:Adipose)

ComparePredicted_plot <- ComparePredicted %>% 
  ggplot() +
  geom_path(aes(x = time, y = Concentration, colour = Origin, linetype = Origin), linewidth = 0.3) +
  facet_wrap(~Organ)+
  scale_color_manual(
    values = c("Predicted" = "blueviolet",
               "PredictedEFSA" = "darkblue")) +
  scale_linetype_manual(
    values = c("Predicted" = "solid",
               "PredictedEFSA" = "dashed")
  )+
  CP_theme +
  ylab("PFOA Organ Concentration (ng/ml)") +
  xlab("Time (days)") +
  theme(legend.position = c(0.8, 0.3),
        legend.text = element_text(size = 42),
        legend.title = element_blank(),
        legend.key.size = unit(0.5, "lines"),
        legend.spacing.y = unit(0, "cm"),
        legend.margin = margin(0, 0, 0, 0),
        strip.text = element_text(size = 54)  
  )
ComparePredicted_plot
ggsave(
  filename = here(OUTPUT, "ComparePredicted.png"),
  plot = ComparePredicted_plot,
  dpi = 1000,
  width = 11,
  height = 7,  
  units = "cm"
)



## Plot urine ####
Abraham_urine <- Abraham_data %>% filter(Matrix == "Urine")
ObsUrine <- Abraham_urine %>%
  rename(time = Time.days) %>%
  rename(AU = MPFOA) %>% 
  mutate(
         AU = AU*0.985, # ng (correcting as in Abraham they mention that 1.5% of PFOA was from a contamination)
         Origin = "Observed") %>% 
  select(c(time, AU)) %>% 
  filter(time <= 10)

Predicted_urine <- RESULTS$OUT_RAW_data %>% select(time, AUr)
# According to the ICRP89, an adult male looses about 1600ml of urine per day
PredUrine <- Predicted_urine %>%
  select(time, AUr) %>%
  mutate(time = time*365) %>% #,
         # AUr.cumulative = AUr,
         # AUr = AUr.cumulative - lag(AUr.cumulative, default = 0),
         # CU = AUr/1600) %>%
  filter(time <=10) %>% 
  mutate(Origin = "Predicted")

PredictedEFSA_urine <- PredEFSA %>% select(time, Aurine_PFOA) 
PredEFSAUrine <- PredictedEFSA_urine %>%
  # mutate(AUr.cumulative = Aurine_PFOA,
  #        AUr = AUr.cumulative - lag(AUr.cumulative, default = 0),
  #        CU = AUr/1600) %>%
  filter(time <=10) %>% 
  mutate(Origin = "Predicted")

Plot_Urine <- ggplot() +
  geom_path(data = PredUrine, aes(x = time, y = AUr, color = "Predicted"), linewidth = 0.3) +
  geom_point(data = ObsUrine, aes(x = time, y = AU, color = "Observed"), size = 0.1) +
  geom_path(data = PredEFSAUrine, aes(x = time, y = Aurine_PFOA, color = "Predicted EFSA 2020"), 
            linetype = "dashed", linewidth = 0.3) +
  scale_color_manual(
    name = "Legend",
    values = c("Predicted" = "blueviolet", 
               "Observed" = "black", 
               "Predicted EFSA 2020" = "darkblue")) +
  CP_theme +
  ylab("PFOA Amount in Urine (ng)") +
  xlab("Time (days)") +
  theme(legend.position = c(0.2, 0.9),
        legend.text = element_text(size = 52),
        legend.title = element_blank(),
        legend.key.size = unit(0.5, "lines"),
        legend.spacing.y = unit(0, "cm"),
        legend.margin = margin(0, 0, 0, 0))
Plot_Urine
ggsave(
  filename = here(OUTPUT, "AbrahamVsPred_Urine.png"),
  plot = Plot_Urine,
  dpi = 1000,
  width = 9,
  height = 7,  
  units = "cm"
)

## Plot faeces ####
Abraham_feces <- Abraham_data %>% filter(Matrix == "Feces")
ObsFeces <- Abraham_feces %>%
  rename(time = Time.days) %>%
  rename(AF = MPFOA) %>% 
  mutate(Origin = "Observed") %>% 
  select(c(time, AF)) %>% 
  filter(time <=450) 
  

Predicted_feces <- RESULTS$OUT_RAW_data %>% select(time, AFe)
# According to the ICRP89, an adult male looses about 150g of feces per day
PredUFeces <- Predicted_feces %>%
  select(time, AFe) %>%
  mutate(time = time*365) %>% #,
         # AFe.cumulative = AFe,
         # AFe = AFe.cumulative - lag(AFe.cumulative, default = 0),
         # CF = AFe/150) %>%
  filter(time <=450) %>% 
  mutate(Origin = "Predicted")

Plot_Feces <- ggplot() +
  geom_path(data = PredUFeces, aes(x = time, y = AFe, color = "Predicted"), linewidth = 0.3)+
  geom_point(data = ObsFeces, aes(x = time, y = AF, color = "Observed"), size = 0.1)+
  scale_color_manual(
    name = "Legend",
    values = c("Predicted" = "blueviolet", 
               "Observed" = "black")) +
  CP_theme +
  ylab("PFOA Amount in Faeces (ng)")+
  xlab("Time (days)")+
  theme(legend.position = c(0.2, 0.9),
        legend.text = element_text(size = 52),
        legend.title = element_blank(),
        legend.key.size = unit(0.5, "lines"),
        legend.spacing.y = unit(0, "cm"),
        legend.margin = margin(0, 0, 0, 0))

Plot_Feces
ggsave(
  filename = here(OUTPUT, "AbrahamVsPred_Feces.png"),
  plot = Plot_Feces,
  dpi = 1000,
  width = 9,
  height = 7,  
  units = "cm"
)





# Calculate AUCs
AUC.Pred.EFSA <- trapz(PredEFSA[ , "time"], PredEFSA[ , "CA_PFOA"])  # ug*day/L
print(AUC.Pred.EFSA)

# Calculate Half life
time <- PredEFSA[ , "time"]/365 # years
conc <- PredEFSA[ , "CA_PFOA"] # ug/L or ng/ml
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

AUCsandHLs <- data.frame(
  Origin = c("Abraham", "Predicted", "Predicted EFSA 2020"),
  HL = c(4.016, RESULTS$ANALYSED_data$HalfLife,HalfLife),
  AUC = c(NA, RESULTS$ANALYSED_data$AUC, AUC.Pred.EFSA)
)



# 
# Observed.df <- ObsHalfLifes %>%
#   filter(species == "human",
#          chemical == "pfoa",
#          parameter== "HalfLife") %>%
#   select(c(value_average,n)) %>%
#   rename(HalfLife = value_average) %>%
#   mutate(value = 1,
#          Origin = "Observed")
# Observed.df$HalfLife <- as.numeric(Observed.df$HalfLife) # years
# Observed.df$n <- as.numeric(Observed.df$n)
# 
# Predicted.df <- data.frame(
#   HalfLife = HalfLife, #RESULTS$HalfLife,
#   Origin = "Predicted",
#   value = 1, n = 1)
# Observed.df <- data.frame(
#   HalfLife = Observed.df$HalfLife,
#   Origin = "Observed",
#   value = 1,
#   n = Observed.df$n)
# 
# HalfLifes <- rbind(Predicted.df, Observed.df)
# 
# range <- c(min(Observed.df$n), max(Observed.df$n))
# 
# Plot_HalfLifes <- ggplot() +
#   geom_violin(
#     data = Observed.df,
#     aes(value, HalfLife),
#     color = "transparent",
#     fill = "grey89") +
#   geom_point(
#     data = Observed.df,
#     aes(value, HalfLife, size = n),
#     color = "black",
#     alpha = 0.5,
#     shape = 20) +
#   geom_point(
#     data = Predicted.df,
#     aes(value, HalfLife),
#     color = "red",
#     alpha = 0.7,
#     size = 10,
#     shape = 18) +
#   labs(y = "Half life (years)") +
#   scale_size_continuous(range = c(1, 10),
#                         name = "Sample size") +
#   theme_minimal() +
#   theme(
#     axis.text.x = element_blank(),
#     axis.ticks.x = element_blank(),
#     axis.title.x = element_blank(),
#     axis.text = element_text(size = 10),
#     axis.title = element_text(size = 10),
#     legend.position = "top"
#   )
# Plot_HalfLifes
# ggsave(filename = here(OUTPUT, "ExpVsSimHalfLife.png"),
#        dpi = 300,
#        width = 17,
#        height = 8,
#        units = "cm")

# 2. Evaluate simulated half life against the half lives observed in Human Biomonitoring (HBM) data ####
# Should be used together with the INPUT_dummy.csv file and simulation results after running it

ObsHalfLifes <- read_csv(here("Input", "HalfLifes.csv"))
# Oral.F <- RESULTS$ANALYSED_data[[4]]$HalfLife
# Oral.M <- RESULTS$ANALYSED_data[[5]]$HalfLife
# Dermal.F <- RESULTS$ANALYSED_data[[6]]$HalfLife
# Inhalation.F <- RESULTS$ANALYSED_data[[2]]$HalfLife

# Prepare observed data

Observed.df <- ObsHalfLifes %>%
  filter(species == "human",
         chemical == "pfoa",
         parameter == "HalfLife") %>%
  select(c(value_average, n)) %>%
  rename(HalfLife = value_average) %>%
  mutate(value = 0.75,
         Origin = "Observed") %>%
  mutate(HalfLife = as.numeric(HalfLife),  # years
         n = as.numeric(n))

Observed2.df <- Observed.df %>%
  mutate(value = 1.5)

Observed.df <- rbind(Observed.df, Observed2.df)

# Prepare predicted data for each Exposure and sex
Predicted.df <- data.frame(
  # Exposure = c(rep("Oral", length(c(Oral.F, Oral.M))),
  #              rep("Dermal", length(Dermal.F)),
  #              rep("Inhalation", length(Inhalation.F))),
  # Sex = c(rep("F", length(Oral.F)),
  #         rep("M", length(Oral.M)),
  #         rep("F", length(Dermal.F)),
  #         rep("F", length(Inhalation.F))), 
  HalfLife = c(Oral.F, Oral.M, Dermal.F, Inhalation.F),
  Origin = "Predicted") %>% 
  mutate(
    HalfLife = as.numeric(str_remove(HalfLife, "_years")),
    value = case_when(
      Sex == "F" ~ 0.75,  
      Sex == "M" ~ 1.5   
    )
  )


# Plot
Halflifeplot <- 
  ggplot() +
  geom_violin(data = Observed.df, aes(x = 0.75, y = HalfLife), 
              fill = "grey89", color = NA, width = 0.5, trim = FALSE) +
  geom_point(data = Observed.df, aes(x = 0.75, y = HalfLife, size = n),
             color = "grey70", alpha = 0.5, position = position_jitter(width = 0.05)) +
  # geom_point(data = filter(Predicted.df), 
  #            aes(x = 0.75, y = HalfLife, color = Exposure), 
  #            shape = 18, size = 5, alpha = 0.9, position = position_jitter(width = 0.25)) +
  # 
  # geom_violin(
  #   data = Observed.df, aes(x = 1.5, y = HalfLife),
  #   fill = "grey89", color = NA, width = 0.5, trim = FALSE) +
  # geom_point(data = Observed.df, aes(x = 1.5, y = HalfLife, size = n),
  #            color = "grey70", alpha = 0.5, position = position_jitter(width = 0.05)) +
  # geom_point(data = filter(Predicted.df, Sex == "M"),
  #            aes(x = 1.5, y = HalfLife, color = Exposure),
  #            shape = 18, size = 5,  alpha = 0.9, position = position_jitter(width = 0.01)) +
  # 
  # scale_color_manual(values = c("Oral" = "#8934AA",
  #                               "Dermal" = "#238EFF", 
  #                               "Inhalation" = "#F5D475")) +
  scale_x_continuous(breaks = c(0.75, 1.5),       
                     labels = c("Female", "Male")) + 
  scale_size_continuous(range = c(1, 5)) +
  
  labs(x = "", y = "Half life (years)") +
  guides(size = guide_legend(title = "HBM sample size")) +
  theme_minimal() +
  theme(axis.title.x = element_text(size = 12),
        axis.text.x = element_text(size = 11),
        legend.position = "right")
Halflifeplot
ggsave(filename = here(OUTPUT, "NoLifestageHalf.life.png"), 
       dpi = 300,
       width = 12,      
       height = 8,      
       units = "cm")


# 3. Evaluate against Olsen data ####
# This should be done after performing reverse dosimetry to define which exposure concentration is needed to reach the measured plasma concentration for each participant of the study of Olsen et al. 
# Concentration at Olsen experiment start
# Results of reverse dosimetry and exposure scenario is found in OlsenData.csv and can be used directly as input to the model
InputData <- read.csv(here("Input", "INPUT_dummy.csv"))

CPredictedObserved.df <- data.frame(
  Idcode = INPUT_dummy$Idcode,
  exp = INPUT_dummy$exp,
  expSTOP = INPUT_dummy$expSTOP, # this is the time of the stop of exposure in days
  CP_measured = INPUT_dummy$CP_measured, # measured plasma PFOA concentration
  samplingT = INPUT_dummy$samplingT # time at which the plasma concentration was measured
)

## Regression on the plasma concentration 

# Import RESULTS from the PBK simulation
PBK_OUT <- RESULTS$OUT_RAW_data

CPredictedObserved.df <- CPredictedObserved.df %>%
  # Fint the predicted concentration at the time sampling time (samplingT)
  mutate(
    CP_predicted = map2_dbl( 
      PBK_OUT,samplingT,
      ~ {
        idx <- which.min(abs(.x$time - .y))
        .x$CP[idx]
      }
    )
  ) 

CPRegression <- lm(CP_predicted~CP_measured, data=PredictedObserved.df) #lm(y~x) (y is the dependent variable)
summary(CPRegression)


## Regression on the half life 

ANALYSED_data <- RESULTS$ANALYSED_data

PredictedObserved.df <- PredictedObserved.df %>% mutate(
  Idcode = seq_along(ANALYSED_data),
  HalfLife = sapply(ANALYSED_data, function(x) x$HalfLife)
  ) %>%
  separate(col = HalfLife, into = c("HL_predicted", "unit"), sep = "_") %>%
  mutate(
    HL_predicted = as.numeric(HL_predicted),
    HL_predicted = round(HL_predicted,1)
  ) %>% 
  mutate(
    HL_observed = c(3.6,
                 2.3,
                 2.8,
                 3.6,
                 3.3,
                 2.3,
                 3.3,
                 6.9,
                 3.8,
                 3,
                 4.2,
                 1.5,
                 3.5,
                 2.8,
                 9.1,
                 4.8,
                 3.8,
                 1.6,
                 7,
                 2.9,
                 4,
                 3.4,
                 3.7,
                 4.6,
                 3.3,
                 2.9)
  ) 

HLRegression <- lm(HL_observed~HL_predicted, data=PredictedObserved.df)
summary(HLRegression)

## Plots 
PlotHLRegression <- HLRegression %>%
  ggplot(aes(log(HL_predicted), log(HL_observed))) +
  geom_smooth(method='lm', color = "black", se = TRUE) +
  geom_abline(intercept = 0, slope = 1, linetype = "solid", linewidth = 0.5, color = "grey50") +  
  geom_abline(intercept = log(1.1), slope = 1, linetype = "dashed", linewidth = 0.5, color = "grey50") +  # +10% line
  geom_abline(intercept = log(0.9), slope = 1, linetype = "dashed", linewidth = 0.5, color = "grey50") +  # -10% line
  geom_abline(intercept = log(2), slope = 1, linetype = "dotted", linewidth = 0.5, color = "grey50") +  # 2-fold upper
  geom_abline(intercept = log(0.5), slope = 1, linetype = "dotted", linewidth = 0.5, color = "grey50") +  # 2-fold lower
  geom_point(color = "darkred", size = 1) +
  theme_minimal() +
  labs(title = "Predicted Vs Observed Halflife (log years)",
       x = "Predicted", y = "Observed") +
  theme(plot.title = element_text(size = 10, margin = margin(b = 20)),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 8))
PlotHLRegression
ggsave(filename = here(OUTPUT, "PlotHLRegression.png"), 
       dpi = 300,
       width = 12,      
       height = 8,      
       units = "cm")


PlotCRegression <- CRegression %>%
  ggplot(aes(x = log(CP_final_predicted), y = log(CP_final))) +
  geom_smooth(method = 'lm', color = "black", se = TRUE) +
  geom_abline(intercept = 0, slope = 1, linetype = "solid", linewidth = 0.5, color = "grey50") +  
  geom_abline(intercept = log(1.1), slope = 1, linetype = "dashed", linewidth = 0.5, color = "grey50") +  # +10% line
  geom_abline(intercept = log(0.9), slope = 1, linetype = "dashed", linewidth = 0.5, color = "grey50") +  # -10% line
  geom_abline(intercept = log(2), slope = 1, linetype = "dotted", linewidth = 0.5, color = "grey50") +  # 2-fold upper
  geom_abline(intercept = log(0.5), slope = 1, linetype = "dotted", linewidth = 0.5, color = "grey50") +  # 2-fold lower
  labs(title = "Predicted Vs Observed Serum Concentration (log ng/ml)",
       x = "Predicted", y = "Observed") +
  geom_point(color = "darkred", size = 1) +
  theme_minimal()+
  theme(plot.title = element_text(size = 10, margin = margin(b = 20)),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 8))
PlotCRegression
ggsave(filename = here(OUTPUT, "PlotCRegression.png"), 
       dpi = 300,
       width = 12,      
       height = 8,      
       units = "cm")


# 
# ObsHalfLifes <- read_csv(here("Input", "HalfLifes.csv"))
# ObsPlasmaConc <- read_csv(here("Input", "ObservedPFOA_CPlasma.csv"))
# 
