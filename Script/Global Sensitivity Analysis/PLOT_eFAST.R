# --------------------------------------------------------------------------- #
# SCRIPT FOR PLOTING eFAST RESULTS
# By: Chrysanthi Pachoulide
# Date: 15-04-2025
# --------------------------------------------------------------------------- #

rm(list=ls()) # to clear out the global environment


# Packages
library(here)
library(tidyverse)
library(scales)
library(showtext)
font_add(family = "Garamond", regular = "GARA.TTF")
showtext_auto()

# Set output storage directory
OUTPUT <- here("Output", format(Sys.Date(), "%Y-%m-%d"), format(Sys.time(), "%H-%M-%S"))
dir.create(OUTPUT, recursive = TRUE)


# load eFAST Rdata file
load(here("Output", "2025-04-8", "2025-04-08 14-57-10", "ExperienceResults.RData"))
Parameters <- read_csv(here("Input", "Parameters.GSA.csv"))

P <- Parameters %>% select(Abbreviation, Distribution, Initial_value, Binf, Bsup, Mean, sdlog, zscore) %>%
  filter(Abbreviation %in% c("VILc", "fup", "COral", "R_PTL", "QT", "SF_OAT", "Km_OAT4c", "BW", "GFRc", "QC")) # "SA_SI", "QKc", "SF_OATP1B1", "Vmax_OAT4c", "Km_OATP1B1c", "tco", "Papp_SI", "VPc", "VKc", "VPTc", "Hct", "PAc", "Vmax_OATP1B1c", "Km_OATP1B3c"))
GSA.parms <- P$Abbreviation

# Results
Variance <- eFAST$V # Total variance
names(Variance) <- GSA.parms

Done <- eFAST$D1
names(Done) <- GSA.parms

Dt <- eFAST$Dt
names(Dt) <- GSA.parms

first_order <- Done / Variance  # D1 is the estimated Variance of the Conditional Expectation (VCE) with respect to each factor, normalized by total variance
names(first_order) <- GSA.parms

total_order <- Dt / eFAST$V  # Dt is the estimated VCE with respect to each factor complementary set of factors ("all but Xi"), normalized by total variance
names(total_order) <- GSA.parms

lowry_data <- data.frame(
  Parameter = GSA.parms,
  Main.Effect = first_order,
  Interaction = total_order) %>% 
  mutate(
    Main.Effect = Main.Effect/sum(Main.Effect), # to normalise
    Interaction = Interaction/sum(Interaction)  # to normalise
  ) 

write.csv(lowry_data, file = here(OUTPUT, "eFASTresults.csv"), row.names = FALSE)

ordered_data <- lowry_data %>% 
  arrange(desc(Main.Effect)) %>%  
  mutate(
    Parameter = factor(Parameter, levels = Parameter), # parameter to factor, fixing factor levels to the current order (descending Main.Effect)
    Total.effect = Main.Effect + Interaction,
    Cumulative.main = cumsum(Main.Effect) # cumulative sum of Main.Effects
  )

long_data <- ordered_data %>% pivot_longer(c(Main.Effect, Interaction))

lowry_plot <- ggplot(long_data) +
  geom_col(aes(x = Parameter, y = value, fill = name),
           position = position_dodge2(padding = -0.4), # padding adds overlap
           width = 0.7,
           alpha = 0.6  # Increased transparency
  ) +
  geom_ribbon(
    data = ordered_data,
    aes(x = as.numeric(Parameter),
        ymin = lag(Cumulative.main, default = 0),
        ymax = Cumulative.main),
    fill = "grey50", alpha = 0.2
  ) +
  scale_fill_manual(
    values = c("Main.Effect" = "darkblue", "Interaction" = "blueviolet"), 
    labels = c("Main Effect", "Interaction")
  ) +
  # Proper 100% scale
  scale_y_continuous(
    limits = c(0, 1),
    expand = c(0, 0),
    labels = percent_format(),
    name = "Sensitivity Index"
  ) +
  labs(x = "Parameters") +
  theme_minimal(base_size = 12) +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top",
    legend.title = element_blank()
  )

lowry_plot
ggsave(filename = here(OUTPUT, "lowry_plot.png"), 
       dpi = 300,
       width = 17,      
       height = 8,      
       units = "cm"     
       )