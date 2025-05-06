# --------------------------------------------------------------------------- #
# SCRIPT FOR COMPARING SIMULATED AGAINST OBSERVED DATA
# By: Chrysanthi Pachoulide
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


# Validate against Olsen data ####
# This should be done after performing reverse dosimetry to define which exposure concentration is needed to reach the measured plasma concentration for each participant of the study of Olsen et al. 
# Concentration at Olsen experiment start
# Results of reverse dosimetry and exposure scenario is found in OlsenData.csv and can be used directly as input to the model
OlsenData <- read.csv(here("Input", "OlsenData.csv"))

PredictedObserved.df <- data.frame(
  Idcode = OlsenData$Idcode,
  expSTOP = OlsenData$expSTOP/365,
  CP_initial = OlsenData$CP_initial, # plasma PFOA concentration at the begining of the study
  CP_final = OlsenData$CP_final,
  expCONC = OlsenData$expCONC
)

## Regression on the plasma concentration ####

# Import RESULTS from the PBK simulation
PBK_OUT <- RESULTS$OUT_RAW_data

PredictedObserved.df <- PredictedObserved.df %>%
  mutate(
    CPatexpSTOP = map2_dbl( # Find the predicted concentration at the time of the stop of exposure to check if it's the same as the one measured
      PBK_OUT,
      expSTOP,
      ~ {
        idx <- which.min(abs(.x$time - .y))
        .x$CP[idx]
      }
    )
  ) %>%
  mutate(
    CP_final_predicted = map_dbl( # Find the predicted concentration at the end of the study
      PBK_OUT,
      ~ {
        last_row <- nrow(.x)
        .x$CP[last_row]
      }
    )
  )

CRegression <- lm(CP_final~CP_final_predicted, data=PredictedObserved.df)
summary(CRegression)


## Regression on the half life ####

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

## Plots ####
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
# # Plot Halflife ####
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
# 
# 
# # ## Plot Concentration over time ####
# # 
# # ObsPlasma <- ObsPlasmaConc %>%
# #   filter(Timedays <= 450.00) %>% 
# #   mutate(Timedays = Timedays/365) %>% #to years
# #   rename(time = Timedays) %>% 
# #   rename(CP = MPFOAugperL) %>%    # ug/L or ng/ml
# #   mutate(CP = CP - 0.130) %>%     # substracting the pre-existing level of 0.130ug/L from their previous study, as also done in the ref. article: https://doi.org/10.1016/j.envint.2024.109047 (table 3)
# #   mutate(Origin = "Observed")
# # 
# # SimData <- RESULTS$data
# # SimPlasma <- SimData %>% 
# #   select(time, CP) %>% 
# #   mutate(time = time) %>% 
# #   mutate(Origin = "Predicted")
# # 
# # Plot_Plasma <- ggplot() +
# #   geom_path(data = SimData, aes(x = time, y = CP), color = "red", linewidth = 1.5)+
# #   geom_point(data = ObsPlasma, aes(x = time, y = CP), color = "black")+
# #   theme_minimal()+
# #   ylab("Plasma (ng/ml)")+
# #   xlab("Time (years)")
# # Plot_Plasma
# # ggsave(here(OUTPUT, "ObsVsSimConcOverTime.png"), dpi = 300)
# # 
