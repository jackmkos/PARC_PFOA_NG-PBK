# --------------------------------------------------------------------------- #
# SCRIPT FOR AGE-DEPENDENT GFR CALCULATION
# By: Jack Koster
# Date: 16-04-2025 (Updated: 22-07-2025) (Updated 27-02-2026)
# --------------------------------------------------------------------------- #

# Flow GFR calculation (for reference and plotting) ####
GFRc <- 0.18             # L/d, 18% of total renal plasma flow [ICRP 89]
QTotc <- 0.988           # Placeholder, replace with actual value if needed

# Parameters
QKc <- Physio_params[[paste0("Q_kidneyFraction", suffix)]] 
QC <- Physio_params[[paste0("CardOut", suffix)]]              
QK <- QKc/QTotc * QC    # Kidney flow
  
GFR_flow <- GFRc * QK

Variables_df <- Variables_df %>%
  mutate(
    GFR_M_flow = GFRc * (Q_kidneyFraction_M / QTotc * CardOut_M),
    GFR_F_flow = GFRc * (Q_kidneyFraction_F / QTotc * CardOut_F)
  )

# Model Settings (for switch) ####
  
# Choose between age GFR ("Age") or renal plasma flow GFR ("Flow")
GFR_type = "Age"
  
GFR  # age GFR
  
# PBK model code
GFR <- if (GFR_type == "Age") {
  GFR  # age GFR
} else {
  GFRc * QK  # renal plasma flow GFR
}

# Age GFR original with citations (18-03-2025) ####

Variables_df <- Variables_df %>%
  
  # Skin surface area 
  mutate(SkbTarea_F = 9.1*((BW_F*1000)^0.666)) %>%  #Skin barrier (epidermis = stratum corneum and viable epidermis) area; assumed to be the same as the total area of the skin (cm^2) from Husoy
  
  # eGFR = 107.3/[SCr/Q]                      age ≤ 40
  # eGFR = 107.3/[SCr/Q] * 0.988^(age-40)     age > 40  (Pottel 2017, https://doi.org/10.1093/ckj/sfx026)
  
  # mGFR of n=633 potential kidney donors used to validate equation and reference intervals.
  # 107.3 mL/min * (1440/1000) = 154.5 L/d #5% different from current (163 L/d).
  
  # Scr normally distributed, but Q (avg or med Scr for a pop) changes uniquely depending on age and gender. LB and UB of [Scr/Q] = 0.67, 1.33 (1 for most healthy adults)
  # Baseline/neonatal GFR should be 20.0 mL/min: Qi = (20*0.9)/107.3 = 0.1678 mg/dL (Smeets 2022, https://doi.org/10.1681/ASN.2021101326)
  # females: Q increases linearly from 0.1678 mg/dL at birth to 0.70 mg/dL at age 18.
  # males: Q increases linearly from 0.1678 at birth to 0.90 mg/dL at age 18.
  mutate(age = if_else(age > 18, 18, age)) %>%
  mutate(
    SCr_F = if_else(age < 18, 0.1678 + ((0.70 - 0.1678) / 18) * age, 0.70),
    SCr_M = if_else(age < 18, 0.1678 + ((0.90 - 0.1678) / 18) * age, 0.90)
  ) %>%
  # Calculate adult baseline GFR in mL/min
  # GFR_i = 107.3 / (Scr/Q), with Scr fixed at 0.9.
  # Convert to L/d
  mutate( # Note: this was before adding age-adjusted population SCr, and this should be individual/population not population/individual
    GFRi_F = (107.3 * 1.44 * (SkbTarea_F*1e-4) / (0.7/SCr_F)), # (ml/min/1.73m^2 -> L/d ) # m^2 = 1.73 for constant, or SkbTarea*1e-4 for variable
    GFRi_M = (107.3 * 1.44 * (SkbTarea_M*1e-4) / (0.9/SCr_M)) 
  ) %>%
  # Apply adult adjustment # For age <= 40, use the baseline GFR
  # For age > 40, apply an exponential decline with factor 0.988^(age - 40).
  mutate(
    GFR_F = if_else(age <= 40, GFRi_F, GFRi_F * 0.988^(age - 40)),
    GFR_M = if_else(age <= 40, GFRi_M, GFRi_M * 0.988^(age - 40))
  )


# Final Update (22-07-2025) ####

# Estimated neonatal SCr values (mg/dl)
SCr_n_M <- 0.1678
SCr_n_F <- 0.1305

# Individual adult SCr values
SCr_i_M <- 0.90 # Calculate or input bio monitoring data
SCr_i_F <- 0.70

# Average population SCr values (Q)
SCr_p_M <- 0.90
SCr_p_F <- 0.70

# Function to calculate age-dependent GFR
calc_gfr <- function(df) {
  df %>%
    # Age-adjusted individual and populations SCr 
    mutate( 
      SCra_i_M = if_else(age < 18, SCr_n_M + ((SCr_i_M - SCr_n_M) / 18) * age, SCr_i_M),
      SCra_i_F = if_else(age < 18, SCr_n_F + ((SCr_i_F - SCr_n_F) / 18) * age, SCr_i_F),
      SCra_p_M = if_else(age < 18, SCr_n_M + ((SCr_p_M - SCr_n_M) / 18) * age, SCr_p_M),
      SCra_p_F = if_else(age < 18, SCr_n_F + ((SCr_p_F - SCr_n_F) / 18) * age, SCr_p_F)
    ) %>%
    # Baseline GFRs, conversion, scaling # (ml/min/1.73m^2 -> L/day)
    mutate(
      GFRbase_M = (107.3 * 1.44 * (BSA_M)/1.73) / (SCra_i_M/SCra_p_M),
      GFRbase_F = (107.3 * 1.44 * (BSA_F)/1.73) / (SCra_i_F/SCra_p_F)
    ) %>%
    # Exponential decline after age 40
    mutate(
      GFR_M = if_else(age <= 40, GFRbase_M, GFRbase_M * 0.988^(age - 40)),
      GFR_F = if_else(age <= 40, GFRbase_F, GFRbase_F * 0.988^(age - 40))
    )
}

write.csv(Variables_df, here("Input", "PhysioVariables.csv"), row.names = FALSE)

## Update for variable Scr (22-07-2025) ####
  
Variables_df <- Variables_df %>%
    
  mutate(
    SCrQ_M = if_else(age < 18,
                         0.1678 + ((0.90 - 0.1678)/18)*age,
                         0.90),
  # simulate individual SCr around the age/sex mean
    CV = 0.2,
    SCr_M = exp(rnorm(n(),
                        mean = log(Q_GFRi_M),
                        sd   = log(1 + CV))),
  # compute baseline GFR with true SCr
    GFR_M_base = (107.3 * 1.44 * (BSA_M)/1.73) / (SCr_M/SCrQ_M),
  # apply age-related decline
    GFR_M = if_else(age <= 40, GFR_M_base, GFR_M_base * 0.988^(age - 40))
  )


# Plots ####
p_GFR <- ggplot() + 
  geom_path(data = Variables_df, aes(age, GFR_M/1.44, colour = "Male Age")) +
  geom_path(data = Variables_df, aes(age, GFR_F/1.44, colour = "Female Age")) +
  geom_path(data = Variables_df, aes(age, GFR_M_flow/1.44, colour = "Male Flow")) +
  geom_path(data = Variables_df, aes(age, GFR_F_flow/1.44, colour = "Female Flow")) +
  scale_colour_manual(values = c("Male Age" = "cornflowerblue",
                                 "Female Age" = "darkorchid3",
                                 "Male Flow" = "darkgreen",
                                 "Female Flow" = "brown"),
                      name = "") +
  theme_minimal(base_size = 10)+
  ylab("GFR (ml/min)") +
  xlab("Age (years)")

p_GFR
