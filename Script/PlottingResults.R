# --------------------------------------------------------------------------- #
# SCRIPT FOR COMPARING SIMULATED AGAINST OBSERVED DATA
# By: Chrysanthi Pachoulide
# Date: 17-06-2025
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

# Input files ####
Input <- read.csv(here("Input", "INPUT_dummy.csv")) %>% 
  filter(Study == "dummy") # Input file for observed data, select the correct study
Input <- Input %>% filter(!is.na(expAGE)) # needed if lifestage is used
Lit.HalfLifes <- read_csv(here("Input", "HalfLifes.csv")) # Half lives reported in literature
load(here("Output", "2025-06-17", "17-12-04", "RESULTS_PFOA_PBK.RData")) # Change for the correct file

PBK_OUT <- RESULTS$OUT_RAW_data # PBK model results per subject
ANALYSED_data <- RESULTS$ANALYSED_data # PBK model results per subject, analysed (contains AUC and half lives)

# Evaluate against HBM data ####
Observed <- data.frame(
  Idcode = Input$Idcode,
  exp = Input$exp,
  expSTOP = Input$expSTOP, # this is the time of the stop of exposure in days
  CP_measured = Input$CP_measured, # measured plasma PFOA concentration
  samplingT = Input$samplingT, # time at which the plasma concentration was measured
  HL_measured = Input$HL_observed, # half life reported in the HBM study
  sex = Input$sex
) 

## Regression on the plasma concentration ####
CP_PredictedObserved <- Observed %>%
  # Fint the predicted concentration at the time sampling time (samplingT)
  mutate(
    CP_predicted = map2_dbl( 
      PBK_OUT,samplingT, ~ {
        idx <- which.min(abs(.x$time - .y))
        .x$CP[idx]
      }
    )
  ) %>% select(CP_measured, CP_predicted) 

CP_Regression <- lm(CP_predicted~CP_measured, data=CP_PredictedObserved) #lm(y~x) (y is the dependent variable)
summary(CP_Regression)

CP_Regression_plot <- CP_Regression %>%
  ggplot(aes(log(CP_predicted), log(CP_measured))) +
  # geom_smooth(method='lm', color = "black", se = TRUE) +
  geom_abline(intercept = 0, slope = 1, linetype = "solid", linewidth = 0.3, color = "grey50") +  
  geom_abline(intercept = log(2), slope = 1, linetype = "dashed", linewidth = 0.2, color = "grey50") + #2 fold
  geom_abline(intercept = log(0.5), slope = 1, linetype = "dashed", linewidth = 0.2, color = "grey50") + #2 fold
  geom_abline(intercept = log(3), slope = 1, linetype = "dotted", linewidth = 0.1, color = "grey50") +  #3 fold
  geom_abline(intercept = log(1/3), slope = 1, linetype = "dotted", linewidth = 0.1, color = "grey50") + #3 fold
  geom_point(color = "black", size = 0.2) +
  CP_theme +
  labs(x = "log10(Predicted Serum Concentration) (ng/ml)", y = "log10(Observed Serum Concentration) (ng/ml)")
CP_Regression_plot

ggsave(filename = here(OUTPUT, "CP_Regression_plot.png"),
       plot = CP_Regression_plot,
       dpi = 1000,
       width = 9,height = 7, units = "cm")


## Regression on the half life ####
HL_PredictedObserved <- Observed %>% mutate(
  Idcode = seq_along(ANALYSED_data),
  HalfLife = sapply(ANALYSED_data, function(x) x$HalfLife)) %>%
  separate(col = HalfLife, into = c("HL_predicted", "unit"), sep = "_") %>%
  mutate(HL_predicted = as.numeric(HL_predicted), 
         HL_predicted = round(HL_predicted,1)) %>% 
  select(HL_predicted, HL_measured)

HL_Regression <- lm(HL_predicted~HL_measured, data=HL_PredictedObserved)
summary(HL_Regression)

HL_Regression_plot <- HL_Regression %>%
  ggplot(aes(log10(HL_predicted), log10(HL_measured))) +
  # geom_smooth(method='lm', color = "black", se = TRUE) +
  geom_abline(intercept = 0, slope = 1, linetype = "solid", linewidth = 0.3, color = "grey50") +  
  geom_abline(intercept = log10(2), slope = 1, linetype = "dashed", linewidth = 0.2, color = "grey50") +  #2 fold
  geom_abline(intercept = log10(0.5), slope = 1, linetype = "dashed", linewidth = 0.2, color = "grey50") + #2 fold
  geom_abline(intercept = log10(3), slope = 1, linetype = "dotted", linewidth = 0.1, color = "grey50") +  #3 fold
  geom_abline(intercept = log10(1/3), slope = 1, linetype = "dotted", linewidth = 0.1, color = "grey50") + #3 fold
  geom_point(color = "black", size = 0.2) +
  CP_theme +
  labs(x = "log10(Predicted Half life) (years)", y = "log10(Observed Half life) (years)")

HL_Regression_plot

ggsave(filename = here(OUTPUT, "HL_Regression_plot.png"),
       plot = HL_Regression_plot,
       dpi = 1000,
       width = 9,height = 7, units = "cm")

## Plot predicted half lives over those reported in literature ####


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
  HL_predicted = HL_PredictedObserved$HL_predicted, 
  sex = sapply(ANALYSED_data, function(x) x$sex)) %>% 
  mutate(sex = str_remove(sex, "_unitless"), 
         value = case_when(sex == "F" ~ 0.75,  
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
  
  labs(x = "", y = "Half life (years)") +
  guides(size = guide_legend(title = "HBM sample size")) +
  CP_theme +
  theme(legend.position = "botom")
HL_violin_plot

ggsave(filename = here(OUTPUT, "HL_violin_plot.png"),
       plot = HL_violin_plot,
       dpi = 1000,
       width = 9,height = 7, units = "cm")
