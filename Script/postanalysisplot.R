library(here)
library(tidyverse)
library(deSolve)
library(PKNCA)
library(pracma)

ObsHalfLifes <- read_csv(here("Input", "HalfLifes.csv"))
Oral.F <- RESULTS$ANALYSED_data[[4]]$HalfLife
Oral.M <- RESULTS$ANALYSED_data[[5]]$HalfLife
Dermal.F <- RESULTS$ANALYSED_data[[6]]$HalfLife
Inhalation.F <- RESULTS$ANALYSED_data[[2]]$HalfLife

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
  Exposure = c(rep("Oral", length(c(Oral.F, Oral.M))),
               rep("Dermal", length(Dermal.F)),
               rep("Inhalation", length(Inhalation.F))),
  Sex = c(rep("F", length(Oral.F)),
        rep("M", length(Oral.M)),
        rep("F", length(Dermal.F)),
        rep("F", length(Inhalation.F))), 
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
NoLifestageHalf <- 
  ggplot() +
  geom_violin(data = Observed.df, aes(x = 0.75, y = HalfLife), 
              fill = "grey89", color = NA, width = 0.5, trim = FALSE) +
  geom_point(data = Observed.df, aes(x = 0.75, y = HalfLife, size = n),
    color = "grey70", alpha = 0.5, position = position_jitter(width = 0.05)) +
  geom_point(data = filter(Predicted.df, Sex == "F"), 
             aes(x = 0.75, y = HalfLife, color = Exposure), 
             shape = 18, size = 5, alpha = 0.9, position = position_jitter(width = 0.25)) +
  
  geom_violin(
    data = Observed.df, aes(x = 1.5, y = HalfLife),
    fill = "grey89", color = NA, width = 0.5, trim = FALSE) +
  geom_point(data = Observed.df, aes(x = 1.5, y = HalfLife, size = n),
    color = "grey70", alpha = 0.5, position = position_jitter(width = 0.05)) +
  geom_point(data = filter(Predicted.df, Sex == "M"),
    aes(x = 1.5, y = HalfLife, color = Exposure),
    shape = 18, size = 5,  alpha = 0.9, position = position_jitter(width = 0.01)) +
  
  scale_color_manual(values = c("Oral" = "#8934AA",
                                "Dermal" = "#238EFF", 
                                "Inhalation" = "#F5D475")) +
  scale_x_continuous(breaks = c(0.75, 1.5),       
                     labels = c("Female", "Male")) + 
  scale_size_continuous(range = c(1, 5)) +
  
  labs(x = "", y = "Half life (years)") +
  guides(size = guide_legend(title = "HBM sample size")) +
  theme_minimal() +
  theme(axis.title.x = element_text(size = 12),
        axis.text.x = element_text(size = 11),
        legend.position = "right")
NoLifestageHalf
ggsave(filename = here(OUTPUT, "NoLifestageHalf.life.png"), 
       dpi = 300,
       width = 12,      
       height = 8,      
       units = "cm")
