# --------------------------------------------------------------------------- #
# SCRIPT FOR COMPARING SIMULATED AGAINST OBSERVED DATA
# By: Chrysanthi Pachoulide
# Date: 14-04-2025
# --------------------------------------------------------------------------- #

rm(list=ls()) # to clear out the global environment

# Packages
library(here)
library(tidyverse)
library(ggplot2)


# Set output storage directory
OUTPUT <- here("Output", format(Sys.Date(), "%Y-%m-%d"), format(Sys.time(), "%H-%M-%S"))
dir.create(OUTPUT, recursive = TRUE)

# Load input files
load(here("Output", "2025-04-14", "23-52-28", "Oral", "RESULTS_PFOA_PBK.RData"))
ObsHalfLifes <- read_csv(here("Input", "HalfLifes.csv"))
ObsPlasmaConc <- read_csv(here("Input", "ObservedPFOA_CPlasma.csv"))

# Plot Halflife ####

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
  HalfLife = HalfLife, #RESULTS$HalfLife,
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
    axis.title = element_text(size = 10),
    legend.position = "top"
  )
Plot_HalfLifes
ggsave(filename = here(OUTPUT, "ExpVsSimHalfLife.png"), 
       dpi = 300,
       width = 17,      
       height = 8,      
       units = "cm")


# ## Plot Concentration over time ####
# 
# ObsPlasma <- ObsPlasmaConc %>%
#   filter(Timedays <= 450.00) %>% 
#   mutate(Timedays = Timedays/365) %>% #to years
#   rename(time = Timedays) %>% 
#   rename(CP = MPFOAugperL) %>%    # ug/L or ng/ml
#   mutate(CP = CP - 0.130) %>%     # substracting the pre-existing level of 0.130ug/L from their previous study, as also done in the ref. article: https://doi.org/10.1016/j.envint.2024.109047 (table 3)
#   mutate(Origin = "Observed")
# 
# SimData <- RESULTS$data
# SimPlasma <- SimData %>% 
#   select(time, CP) %>% 
#   mutate(time = time) %>% 
#   mutate(Origin = "Predicted")
# 
# Plot_Plasma <- ggplot() +
#   geom_path(data = SimData, aes(x = time, y = CP), color = "red", linewidth = 1.5)+
#   geom_point(data = ObsPlasma, aes(x = time, y = CP), color = "black")+
#   theme_minimal()+
#   ylab("Plasma (ng/ml)")+
#   xlab("Time (years)")
# Plot_Plasma
# ggsave(here(OUTPUT, "ObsVsSimConcOverTime.png"), dpi = 300)
# 
