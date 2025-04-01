# Validation file
# Experimental Vs Predicted

# Input

PredData.sixty.once <- read_csv("C:/Users/pacho003/OneDrive - Wageningen University & Research/CP_L_R/PARC_PFOA_mechanistic/Output/2025-03-28/2025-03-28 19-33-03.554412/output.csv") %>% 
  rename(time = Days)
PredData.twenty.chron <- read_csv("C:/Users/pacho003/OneDrive - Wageningen University & Research/CP_L_R/PARC_PFOA_mechanistic/Output/2025-03-28/2025-03-28 22-57-39.398024/output.csv") %>% 
  rename(time = Days)
PredData.fourty.chron <- read_csv("C:/Users/pacho003/OneDrive - Wageningen University & Research/CP_L_R/PARC_PFOA_mechanistic/Output/2025-03-28/2025-03-28 23-27-29.094606/output.csv") %>% 
  rename(time = Days)

## AUC and Half life ####

AUC <- trapz(PredData.fourty.chron$time, PredData.fourty.chron$CP)   # ug*day/L

# Calculate predicted half-life
time <- PredData.fourty.chron$time                                   # days
conc <- PredData.fourty.chron$CP                                     # ug/L or ng/ml
Cmax <- max(conc)
Tmax <- time[which.max(conc)]
tlast <- max(time[conc > 0])

half_life <- pk.calc.half.life(
  conc,
  time,
  Tmax,
  tlast
)

HalfLife.sixty.once <- half_life$half.life/365                          # half-life in years
HalfLife.twenty.chron <- half_life$half.life/365                          # half-life in years
PredData.fourty.chron <- half_life$half.life/365

# Experimental
ExpData <- read_csv("C:/Users/pacho003/OneDrive - Wageningen University & Research/CP_L_R/PARC_PFOA_mechanistic/Input/HalfLifes.csv")

Experimental.df <- ExpData %>%
  filter(species == "human",
         chemical == "pfoa",
         parameter== "HalfLife") %>%
  select(c(value_average,n)) %>%
  rename(HalfLife = value_average) %>% 
  mutate(value = 1, 
         Origin = "Experimental")
Experimental.df$HalfLife <- as.numeric(Experimental.df$HalfLife) # years
Experimental.df$n <- as.numeric(Experimental.df$n)

Predicted.df <- data.frame(
  HalfLife = PredData.fourty.chron,
  Origin = "Predicted",
  value = 1, n = 1)

Plot_HalfLifes <- ggplot() +
  geom_violin(
    data = Experimental.df, 
    aes(value, HalfLife),
    color = "transparent",
    fill = "grey89"
  ) +
  geom_point(
    data = Experimental.df,
    aes(value, HalfLife, size = n),  # Ensure 'n' is numeric!
    color = "grey70",
    shape = 20  
  ) +
  geom_point(
    data = Predicted.df,
    aes(value, HalfLife),
    color = "slateblue3",
    size = 5,
    shape = 18
  ) +
  ylab("Half life (years)") +
  scale_size_continuous(range = c(1, 10)) +  # Customize size range
  theme_CP() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank()
  )

Plot_HalfLifes
ggsave("ExpVsSimHalfLife.png", dpi = 300)


## Experimental Vs Simulated ####
ExpPlasma <- read_excel("C:/Users/pacho003/OneDrive - Wageningen University & Research/CP_L_R/PARC_PFOA_mechanistic/Input/Experimental.Plasma.PFOA.xlsx", 
                        col_types = c("numeric", "numeric"))


ExpPlasma <- ExpPlasma %>% 
  rename(Days = Time_days) %>% 
  rename(CP = MPFOA_µg_per_L) %>%    # ug/L or ng/ml
  mutate(CP = CP - 0.130) %>%        # substracting the pre-existing level of 0.130ug/L from their previous study, as also done in the ref. article: https://doi.org/10.1016/j.envint.2024.109047 (table 3)
  filter(Days <=TSTOP)


Plot_Plasma <- ggplot()+
  geom_path(data = output.PFOA.df, aes(x = Days, y = CP), color = "aquamarine", linewidth = 1.5)+
  geom_point(data = ExpPlasma, aes(x = Days, y = CP), color = "black")+
  theme_CP()+
  ylab("Plasma (ng/ml)")
Plot_Plasma
ggsave("PlasmaExpVsPredicted.png", dpi = 300)
# 
# library(tidyverse)
# library(pkr)
# library(readxl)
# library(caTools)
# 
# # 1. Define reusable functions ----
# read_simulation_data <- function(file_path) {
#   read_csv(file_path) %>% 
#     rename(time = Days)
# }
# 
# calculate_half_life <- function(data) {
#   conc <- data$CP
#   time <- data$time
#   Cmax <- max(conc)
#   Tmax <- time[which.max(conc)]
#   tlast <- max(time[conc > 0])
#   
#   half_life <- pk.calc.half.life(conc, time, Tmax, tlast)
#   half_life$half.life / 365  # Convert to years
# }
# 
# # 2. Read all simulation data at once ----
# simulation_files <- c(
#   "C:/Users/pacho003/OneDrive - Wageningen University & Research/CP_L_R/PARC_PFOA_mechanistic/Output/2025-03-28/2025-03-28 19-33-03.554412/output.csv",
#   "C:/Users/pacho003/OneDrive - Wageningen University & Research/CP_L_R/PARC_PFOA_mechanistic/Output/2025-03-28/2025-03-28 22-57-39.398024/output.csv",
#   "C:/Users/pacho003/OneDrive - Wageningen University & Research/CP_L_R/PARC_PFOA_mechanistic/Output/2025-03-28/2025-03-28 23-27-29.094606/output.csv"
# )
# 
# simulation_data <- map(simulation_files, read_simulation_data)
# names(simulation_data) <- c("sixty.once", "twenty.chron", "fourty.chron")
# 
# # 3. Calculate half-lives ----
# predicted_half_lives <- map_dbl(simulation_data, calculate_half_life)
# 
# # 4. Prepare experimental data ----
# exp_data_path <- "C:/Users/pacho003/OneDrive - Wageningen University & Research/CP_L_R/PARC_PFOA_mechanistic/Input/HalfLifes.csv"
# exp_plasma_path <- "C:/Users/pacho003/OneDrive - Wageningen University & Research/CP_L_R/PARC_PFOA_mechanistic/Input/Experimental.Plasma.PFOA.xlsx"
# 
# experimental_df <- read_csv(exp_data_path) %>%
#   filter(species == "human",
#          chemical == "pfoa",
#          parameter == "HalfLife") %>%
#   transmute(
#     HalfLife = as.numeric(value_average),
#     n = as.numeric(n),
#     value = 1, 
#     Origin = "Experimental"
#   )
# 
# predicted_df <- tibble(
#   HalfLife = predicted_half_lives,
#   Origin = "Predicted",
#   value = 1, 
#   n = 1
# )
# 
# # 5. Create plots ----
# # Half-life comparison plot
# plot_half_lives <- ggplot() +
#   geom_violin(
#     data = experimental_df, 
#     aes(value, HalfLife),
#     color = "transparent",
#     fill = "grey89"
#   ) +
#   geom_point(
#     data = experimental_df,
#     aes(value, HalfLife, size = n),
#     color = "grey70",
#     shape = 20  
#   ) +
#   geom_point(
#     data = predicted_df,
#     aes(value, HalfLife),
#     color = "slateblue3",
#     size = 5,
#     shape = 18
#   ) +
#   ylab("Half life (years)") +
#   scale_size_continuous(range = c(1, 10)) +
#   theme_CP() +
#   theme(
#     axis.text.x = element_blank(),
#     axis.ticks.x = element_blank(),
#     axis.title.x = element_blank()
#   )
# 
# ggsave("ExpVsSimHalfLife.png", plot_half_lives, dpi = 300)
# 
# # Plasma concentration plot
# exp_plasma <- read_excel(exp_plasma_path, col_types = c("numeric", "numeric")) %>% 
#   rename(Days = Time_days, CP = MPFOA_µg_per_L) %>% 
#   mutate(CP = CP - 0.130) %>% 
#   filter(Days <= TSTOP)
# 
# plot_plasma <- ggplot() +
#   geom_path(data = output.PFOA.df, aes(Days, CP), color = "aquamarine", linewidth = 1.5) +
#   geom_point(data = exp_plasma, aes(Days, CP), color = "black") +
#   theme_CP() +
#   ylab("Plasma (ng/ml)")
# 
# ggsave("PlasmaExpVsPredicted.png", plot_plasma, dpi = 300)
