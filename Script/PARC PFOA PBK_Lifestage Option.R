# --------------------------------------------------------------------------- #
# PBK MODEL FOR PFOA, TO BE USED TOGETHER WITH THE LATEST HBM DATA
# Model File
# CP, 10-11-2024
# --------------------------------------------------------------------------- #

rm(list=ls()) # to clear out the global environment

# Set working directory

HOME <- "C:/Users/pacho003/OneDrive - Wageningen University & Research/C Channel/R/PARC_PFAS_PBPKmodel/Codes_PFAS_models/Chrysanthi_Pachoulide_2024_PFOA_oral_dermal_blood"
# HOME <- "/home/westerj"
setwd(HOME)


# Set output storage directory
workingtime <- gsub(":", "-", Sys.time())
OUTPUT <- file.path("Output/Data", Sys.Date(), workingtime)
dir.create(OUTPUT, recursive = TRUE)
setwd(OUTPUT)


# Load packages

library(lubridate)
library(ggplot2)
library(deSolve)
library(writexl)
library(ggpubr)
library(tidyverse)

# INPUT ####
# ------------------------------------------------------ #

# Input variables
Final_variables_M_df = read_csv("C:/Users/pacho003/OneDrive - Wageningen University & Research/C Channel/R/PARC_PFAS_PBPKmodel/Codes_PFAS_models/Chrysanthi_Pachoulide_2024_PFOA_oral_dermal_blood/Input/2024-11-14/Final_variables_M.csv") %>% as.data.frame()


# EXPOSURE SCENARIO ####
# ------------------------------------------------------ #

exposure_stop <- 50*365         # days
sim_stop <- 50*365 # 80*365     # days

TSTART <- 0
TSTOP <- sim_stop               # years in days
DT <- 1                         # days in hours
TIME <- seq(TSTART,TSTOP,by=DT)


Oralconc <- 0.000187 # ug/kg/day [EFSA 2020; page 143]
Dermconc <- 0.000542 # ug/kg/day; mean of #as.numeric(SumExpPFOA_LB_val[i,14])

# ## Dosing for lifestage simulation ####
# Dosing <- Final_variables_M_df %>% 
#   select(TIME, age, BDW_M, AbsPFOA) %>% 
#   
#   # Oral dose
#   mutate(
#     Oraldose = Oralconc*BDW_M,
#     Oraldose = if_else(age <= exposure_stop,Oraldose,0)
#   ) %>% 
#   
#   # Dermal does
#   mutate(Dermaldose = AbsPFOA*Dermconc*BDW_M) %>% 
#   mutate(Dermaldose = if_else(age <= exposure_stop,Dermaldose,0))
# 
## Dosing for NOT lifestage simulation ####
# Choose age accordingly
BDW <- Final_variables_M_df %>%
  filter(age == 20) %>% pull(BDW_M) 

AbsPFOA <- Final_variables_M_df %>%
  filter(age == 20) %>% pull(AbsPFOA) 

Oraldose <- Oralconc*BDW # ug/day
Dermaldose <- AbsPFOA*Dermconc*BDW # ug/day


# PBK MODEL PARAMETERS ####
# ---------------------------------------------------------------------------- #

# ## For lifestage simulations ####
# parms <- Final_variables_M_df %>% 
#   select(
#     TIME,
#     # Not time dependend variables
#     kAbl, # affinity constant basolateral this is about OAT1 and OAT3 which have affinity to uptake (movement from plasma to cells; this is fitted value for now; kAbl = 0.01 is driving the equilibrium towards uptake into the proximal tubule cells
#     kAap, # affinity constant apical this is about OAT4 which has affinity to re-abs (movement from filtrate to cells; this is fitted value for now; kAap = 0.01 is driving the equilibrium towards re-absorption into the proximal tubule cells
#     fup,  # Unbound fraction in plasmal Fischer et al. 2024 https://doi.org/10.1021/acs.est.3c07415;
#     PL,   # Plasma/liver partition coefficient; Rat tissue data (Kudo et al. 2007)
#     PA,   # Plasma/fat partition coefficient; Rat tissue data (Kudo et al. 2007)
#     PK,   # Plasma/kidney partition coefficient; Rat tissue data (Kudo et al. 2007)
#     PSk,  # Plasma/skin partition coefficient; Rat tissue data (Kudo et al. 2007)
#     PR,   # Plasma/rest of the body partition coefficient; Rat tissue data (Kudo et al. 2007)
#     PG,
#     # Time dependent variables
#     V_skin_M,          
#     V_skinPlasma_M,
#     V_skinTissue_M,
#     V_kidney_M,
#     V_kidneyPlasma_M,
#     V_kidneyTissue_M,
#     V_filtrate_M,
#     V_gut_M,
#     V_liver_M,
#     V_adipose_M,       
#     V_rest_M, 
#     V_plasma_M,
#     Q_skin_M,
#     Q_kidney_M,
#     Q_gut_M,
#     Q_liver_M,
#     Q_hepatic_M, 
#     Q_adipose_M,
#     Q_rest_M, 
#     QUr_M,
#     GFR_M,
#     CLdermalabs_M,      
#     CL_PltPT_M,         
#     CL_FiltPT_M,        
#     CLbiliary_M,
#     CLfecal_M 
#   )

## For NOT lifestage simulations ####

Indep_parms <- Final_variables_M_df %>% select(kAbl, kAap, fup, PL, PF, PK, PSk, PR, PG) %>% distinct()

Age_parms <- Final_variables_M_df %>%
  filter(age == 20) %>% # Select age accordingly
  select(
    V_skin_M,
    V_skinBarrier_M,
    # V_skinPlasma_M,
    # V_skinTissue_M,
    V_kidney_M,
    V_kidneyPlasma_M,
    V_kidneyTissue_M,
    V_filtrate_M,
    V_gut_M,
    V_liver_M,
    V_adipose_M,
    V_rest_M,
    V_plasma_M,
    Q_skin_M,
    Q_kidney_M,
    Q_gut_M,
    Q_liver_M,
    Q_hepatic_M,
    Q_adipose_M,
    Q_rest_M,
    # Q_lungs_M, # Needs to be incorporated when inhalation exposure is included
    QUr_M,
    GFR_M,
    CLdermalabs_M,
    CL_PltPT_M,
    CL_FiltPT_M,
    CLbiliary_M,
    CLfecal_M
  ) 

parms <- c(
  kAbl = Indep_parms$kAbl, # affinity constant basolateral; this is about OAT1 and OAT3 which have affinity for the uptake (plasma to cells)
  kAap = Indep_parms$kAap,  # affinity constant apical; this is about OAT4 which has affinity for the re-abs (movement from filtrate to cells; this is fitted value for now; kAap = 0.01 is driving the equilibrium towards re-absorption into the proximal tubule cells
  fup = Indep_parms$fup,   # Unbound fraction in plasmal Fischer et al. 2024 https://doi.org/10.1021/acs.est.3c07415;
  PL = Indep_parms$PL,     # Plasma/liver partition coefficient; Rat tissue data (Kudo et al. 2007)
  PA = Indep_parms$PF,     # Plasma/fat partition coefficient; Rat tissue data (Kudo et al. 2007)
  PK = Indep_parms$PK,     # Plasma/kidney partition coefficient; Rat tissue data (Kudo et al. 2007)
  PSk = Indep_parms$PSk,   # Plasma/skin partition coefficient; Rat tissue data (Kudo et al. 2007)
  PR = Indep_parms$PR,     # Plasma/rest of the body partition coefficient; Rat tissue data (Kudo et al. 2007)
  PG = Indep_parms$PG,     # Plasma/gut partition coefficient; Rat tissue data (Kudo et al. 2007)
  VSk = Age_parms$V_skin_M,          
  VSkb = Age_parms$V_skinBarrier_M,
  # VSkP = Age_parms$V_skinPlasma_M, For future use when updating skin again
  # VSkT = Age_parms$V_skinTissue_M, For future use when updating skin again
  VK = Age_parms$V_kidney_M,
  VKP = Age_parms$V_kidneyPlasma_M,
  VKT = Age_parms$V_kidneyPlasma_M,
  VFil = Age_parms$V_filtrate_M,
  VG = Age_parms$V_gut_M,
  VL = Age_parms$V_liver_M,
  VA = Age_parms$V_adipose_M,        
  VR = Age_parms$V_rest_M,  
  VP = Age_parms$V_plasma_M, 
  QSk = Age_parms$Q_skin_M, 
  QK = Age_parms$Q_kidney_M, 
  QG = Age_parms$Q_gut_M, 
  QL = Age_parms$Q_liver_M, # Liver artery
  QH = Age_parms$Q_hepatic_M, # Includes Portal (gut) and Liver arteries 
  QA = Age_parms$Q_adipose_M, 
  QR = Age_parms$Q_rest_M,
  # QP = Age_parms$Q_lungs_M, # Needs to be incorporated when inhalation exposure is included
  QUr = Age_parms$QUr_M, 
  GFR = Age_parms$GFR_M, # Used to be called CFil
  CLdermalabs = Age_parms$CLdermalabs_M,       
  CL_PltPT = Age_parms$CL_PltPT_M,          
  CL_FiltPT = Age_parms$CL_FiltPT_M,         
  CLbiliary = Age_parms$CLbiliary_M, 
  CLfecal = Age_parms$CLfecal_M, 
  kfil = Age_parms$GFR_M, # parameter from Trine; in the original code: kfil = 0.2*QK  # Clearance from the kidney to the filtrate compartment (L/h); 20% of bloodstream to QK is cleared for 
  Tm = 5842.308*24*BDW^0.75, # ug/d from ug/h/kg^0.75 parameter from Trine; transporter maximum
  Kt = 0.055, # ug/L; changed from Trine who had it as 55ug # Resorption affinity, changed from 0.055 in the original Loccisano 2011 model (ug)
  CLurine = 0.00000183*24*BDW^(-0.25)
)


# PBK MODEL ####
# ---------------------------------------------------------------------------- #

PBPKmodPFOA_M <- function(t, state, parameters){
  with(as.list(c(state, parameters)), {
    
    # ## Dose ----
    Oraldose <- if_else(t <= exposure_stop, Oralconc * BDW, 0)
    Dermaldose <- if_else(t <= exposure_stop, AbsPFOA * Dermconc * BDW, 0)
    
    
    # ## Parameters ----
    # # Ignore section if lifestage not used
    # 
    # kAbl <- parms$kAbl[t]
    # kAap <- parms$kAap[t]
    # 
    # fup <-  parms$fup[t]
    # 
    # PL <-   parms$PL[t]
    # PA <-   parms$PF[t]
    # PK <-   parms$PK[t]
    # PSk <-  parms$PSk[t]
    # PR <-  parms$PR[t]
    # PG <- parms$PG[t]
    # 
    # VSk <- parms$V_skin_M[t]          
    # VSkb <- parms$V_skinBarrier_M[t]
    # # VSkP <- parms$V_skinPlasma_M[t] For future use when updating skin again
    # # VSkT <- parms$V_skinTissue_M[t] For future use when updating skin again
    # VK <- parms$V_kidney_M[t]
    # VKP <- parms$V_kidneyPlasma_M[t]
    # VKT <- parms$V_kidneyPlasma_M[t]
    # VFil <- parms$V_filtrate_M[t]
    # VG <- parms$V_gut_M[t]
    # VL <- parms$V_liver_M[t]
    # VA <- parms$V_adipose_M[t]        
    # VR <- parms$V_rest_M[t]  
    # VP <- parms$V_plasma_M[t] 
    # QSk <- parms$Q_skin_M[t] 
    # QK <- parms$Q_kidney_M[t] 
    # QG <- parms$Q_gut_M[t] 
    # QL <- parms$Q_liver_M[t] 
    # QH <- parms$Q_hepatic_M[t]  
    # QA <- parms$Q_adipose_M[t] 
    # QR <- parms$Q_rest_M[t]
    # # QP <- parms$Q_lungs_M[t] # Needs to be incorporated when inhalation exposure is included
    # QUr <- parms$QUr_M[t] 
    # GFR <- parms$GFR_M[t] # Used to be called CFil
    # CLdermalabs <- parms$CLdermalabs_M[t]       
    # CL_PltPT <- parms$CL_PltPT_M[t]          
    # CL_FiltPT <- parms$CL_FiltPT_M[t]         
    # CLbiliary <- parms$CLbiliary_M[t] 
    # CLfecal <- parms$CLfecal_M[t] 
    # 
    # 
    ## Concentrations ----
    
    # Organ concentrations (ug/L); these are TOTAL concentrations
    CP <- AP/VP  # Concentration in plasma (ug/L)
    CG <- AG/VG  # Concentration in gut (ug/L)
    CL <- AL/VL  # Concentration in liver (ug/L)
    CA <- AA/VA  # Concentration in adipose (ug/L)
    CR <- AR/VR  # Concentration in rest the body (ug/L)
    
    CSk <- ASk/VSk  # Concentration in skin (ug/L)
    # CSkb <- ASkb/VSkb # Concentration in skin barrier (ug/L) 
    # CSkP <- ASkP/VSkP # Concentration in skin plasma (ug/L) # For future use when updating skin again
    # CSkT <- ASkT/VSkT # Concentration in skin tissue (ug/L) # For future use when updating skin again
    
    CK <- AK/VK  # Concentration in kidney (ug/L)
    # CKP <- AKP/VKP # Concentration in kidney blood
    # CKT <- AKT/VKT # Concentration in kidney tissue
    CFil <- AFil/VFil # Concentration in filtrate compartment
    # CUr <- AUr/VUr # Concentration in urine (ug/L) # Not used
    
    CVG <- CG/PG # Concentration leaving gut (ug/L)
    CVL <- CL/PL  # Concentration leaving liver (ug/L)
    CVA <- CA/PA # Concentration leaving adipose (ug/L)
    CVSk <- CSk/PSk  # Concentration leaving skin (ug/L) 
    # CVSkP <- CSkP/PSk # Concentration leaving skin plasma (ug/L) # For future use when updating skin again
    CVR <- CR/PR  # Concentration leaving rest of the body (ug/L)
    
    CVK <- CK/PK # Concentration leaving the kidney (ug/L)
    # CVKP <- CKP/PK # Concentration leaving kidney blood (ug/L)
    
    
    
    ## Differential equations ----
    # adipose compartment
    dAA <- QA*(CP-CVA) # (ug/h)
    
    # Rest compartment
    dAR <- QR*(CP-CVR) # (ug/h)
    
    # Gut compartment: plasma to gut then to liver. Biliary clearance to gut, based on free concentration in the liver; faecal clearance based on total gut concentration; not only the free fraction as partitioning is not needed
    dAG <- Oraldose + QG*CP - QG*CVG + CLbiliary*CL*fup - CLfecal*CG #(ug/h)
    
    # Excretion fecal: cumulative
    dAEx_feces <- CLfecal*CG #ug/h, used to be CVG but I think this is wrong as it's not the tissue plasma partition that defines fecal excretion
    
    # Liver compartment
    dAL <- QL*CP + QG*CVG - (QL+QG)*CVL - CLbiliary*CL*fup #input from the hepatic artery
    # Rate of PFOA amount change in the liver (ug/h)
    
    ### OLD Kidney compartment ----
    # Kidney compartment
    dAK <- QK*(CP -CVK) + Tm*CFil/(Kt+CFil) - kfil*CK #from Trine's model ;Qfil*CK*Free was introduced to reflect clearance to filtrate compartment
    
    # Filtrate compartment
    dAFil <- kfil*(CK-CFil) - Tm*CFil/(Kt+CFil)
    
    # Storage compartment for urine
    dAdelay <- kfil*CFil - CLurine*Adelay # ug/h
    
    # Urine
    dAurine <- CLurine*Adelay # ug/h
    
    ### Updated Kidney compartment STEP 1 ----
    # ### Kidney: this is basically the kidney blood compartment
    # ### I think we're missing active secretion
    # dAK <- QK*(CP-CVK) + (Vmax*CFil)/(Km+CFil) - GFR*CK*fup
    # #                    re-absorption     ultrafiltration
    # 
    # ### Filtrate compartment: rate of formation of the filtrate(=urine) in the lumen; this is the urinary tract from Aude's lifestage equations
    # dAFil <- GFR*CK*fup - (Vmax*CFil)/(Km+CFil) - QUr*CFil
    # #        ultrafiltration  REABSORPTION  urine flow rate to the bladder
    # 
    # ### Remove this compartment
    # ### Urine: urine in the bladder
    # ### ??? why do we need this compartment? I would try having a simple compartment first; i.e not having an excretion compartment
    # dAUr <- QUr*CFil - kUr*AUr #(ug/h)
    # #       urineflow  urine excretion
    # 
    # ### Excretion urinary: cumulative PFOA concentration in the urine
    # dAEx_urine <- kUr*AUr #(ug/h)
    
    # ### Updated Kidney compartment STRP 2 ----
    # # Kidney Blood
    # dAKP <- QK*(CP-CVKP) - GFR*fup*CKP - CL_PltPT*fup*(CKP - (CKT*kAbl)) #  
    # 
    # # Filtrate
    # dAFil <- GFR*fup*CKP - CL_FiltPT*fup*(CFil- (CKT*kAap)) - QUr*CFil #
    # 
    # # Urine
    # dAUr <- QUr*CFil
    # 
    # # Kidney Tissue
    # dAKT <- CL_PltPT*fup*(CKP - (CKT*kAbl)) + CL_FiltPT*fup*(CFil- (CKT*kAap))
    # 
    ### OLD Skin compartment ----
    dASk <- QSk*(CP-CVSk) + Dermaldose # Rate of PFOA amount change in skin
    
    # ### Updated Skin compartment ----
    # # Skin barrier
    # dASkb <- Dermaldose - CLdermalabs*CSkb #(ug/h)
    # # Skin tissue
    # dASk <- CLdermalabs*CSkb + QSk*(CP - CVSk)
    # 
    ## Plasma compartment ----
    dAP <- - (QSk + QG + QL + QA + QK + QR)*CP +
      QSk*CVSk + (QL+QG)*CVL + QA*CVA + QK*CVK + QR*CVR #(ug/h) #QK*CVKP
    
    
    ## Mass Balance ----
    Atot <- AP + ASk + AG + AL + AA + AK + AFil + Adelay + Aurine + AR + AEx_feces # + AEx_urine + AKP + AKT + AFil + AUr + ASkb
    dInput <- Oraldose + Dermaldose
    MB = Input - Atot
    
    
    
    ## End ----
    
    list(c(dAP, 
           # dASkb, 
           dASk, 
           dAG, 
           dAL, 
           dAA, 
           dAK,
           dAFil,
           dAdelay,
           dAurine,
           # dAKP, 
           # dAKT, 
           # dAFil, 
           # dAUr, 
           dAR, 
           dAEx_feces, 
           dInput), #dAEx_urine
         c(CP = CP, 
           CSk = CSk, #CSkb = CSkb
           CG = CG, CVG = CVG, 
           CL = CL, CVL = CVL, 
           CA = CA, CVA = CVA,
           CK = CK, #CKP = CKP, CKT = CKT, 
           CFil = CFil, CVK = CVK, #CVKP = CVKP, 
           CR = CR, CVR = CVR,
           Atot = Atot,MB = MB)
    )
  })
}

## Initials ####

A_init <- c(AP = 0, 
            # ASkb = 0, 
            ASk = 0, 
            AG = 0, 
            AL = 0, 
            AA = 0, 
            AK = 0,
            AFil = 0,
            Adelay = 0,
            Aurine = 0,
            # AKP = 0,
            # AKT = 0,
            # AFil = 0, 
            # AUr = 0, 
            AR = 0, 
            AEx_feces = 0, 
            #AEx_urine = 0,
            Input = 0)


## Solving the model ####
output_PFOA <- lsoda(y = A_init, 
                     times = TIME, 
                     func = PBPKmodPFOA_M, 
                     parms = parms)
output.PFOA.df <- as.data.frame(output_PFOA) 
output.df <- Final_variables_M_df %>%
  rename(time = TIME) %>%
  left_join(output.PFOA.df) %>% 
  mutate(time = time/365) %>% 
  rename(Years = time) %>% 
  filter(Years <= 50) #according to sim_stop

# RESULTS ####
# ---------------------------------------------------------------------------- #

## Plotting per simulation output.df ####
# Plot_PFOA_doses <- ggplot()+
#   geom_line(data = Variables_df, aes(x = output.df, y = Oraldose_M),col="red")+
#   geom_line(data = Variables_df, aes(x = output.df, y = Dermaldose_M),col="blue")+
#   theme(axis.text.x = element_text(size = 7),axis.text.y = element_text(size = 7), axis.title = element_text(size = 8))+
#   #scale_colour_hue()+
#   ylab("Dose (ug)")
# 
# Plot_PFOA_doses


# PFOA in plasma of one individual
Plot_PFOA_Plasma <- ggplot()+
  geom_path(data = output.df, aes(x = Years, y = CP))+
  theme(axis.text.x = element_text(size = 7),axis.text.y = element_text(size = 7), axis.title = element_text(size = 8))+
  ylab("Plasma (ng/ml)")

Plot_PFOA_Plasma
ggsave("PlasmaConcentration.1.png", dpi = 300)


# PFOA in Kidney of one individual
Plot_PFOA_Kidney <- ggplot()+
  geom_path(data = output.df, aes(x = Years, y = CK, color="Kidney"))+
  # geom_path(data = output.df, aes(x = Years, y = CKT, color="Kidney tissue"))+
  geom_path(data = output.df, aes(x = Years, y = CFil, color="Filtrate"))+
  theme(axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7),
        axis.title = element_text(size = 8)) +
  ylab("Kidney Compartments (ng/ml)") +
  scale_color_manual(values = c("Kidney" = "lightblue",
                                # "Kidney tissue" = "blue",
                                "Filtrate" = "darkblue"))

Plot_PFOA_Kidney
ggsave("KidneyConcentrations.png", dpi = 300)


# PFOA in all organs of one individual
Plot_PFOA_All <- output.df %>% 
  # mutate(CK = CKP + CKT) %>% 
  select(Years, CK, CSk, CL, CA, CG, CR, CP) %>% 
  rename(Kidney = CK, Skin = CSk, Liver = CL, Adipose = CA, Gut = CG, Rest = CR, Plasma = CP) %>% 
  pivot_longer(names_to = "Organ", values_to = "Concentration", Skin:Plasma) %>% 
  ggplot()+
  geom_path(aes(x = Years, y = Concentration, color = Organ)) +
  facet_wrap(~ Organ)+
  theme(axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7), 
        axis.title = element_text(size = 8)) +
  ylab("Organ Concentration (ng/ml)")

Plot_PFOA_All
ggsave("OrganConcentrations.png", dpi = 300)




# ## Plotting per age ####
# # Plot_PFOA_doses <- ggplot()+
# #   geom_line(data = Variables_df, aes(x = age, y = Oraldose_M),col="red")+
# #   geom_line(data = Variables_df, aes(x = age, y = Dermaldose_M),col="blue")+
# #   theme(axis.text.x = element_text(size = 7),axis.text.y = element_text(size = 7), axis.title = element_text(size = 8))+
# #   #scale_colour_hue()+
# #   ylab("Dose (ug)")
# # 
# # Plot_PFOA_doses
# 
# 
# # PFOA in plasma of one individual
# Plot_PFOA_Plasma <- ggplot()+
#   geom_path(data = output.df, aes(x = age, y = CP))+
#   theme(axis.text.x = element_text(size = 7),axis.text.y = element_text(size = 7), axis.title = element_text(size = 8))+
#   ylab("Plasma (ng/ml)")
# 
# Plot_PFOA_Plasma
# ggsave("PlasmaConcentration.1.png", dpi = 300)
# 
# 
# # PFOA in Kidney of one individual
# Plot_PFOA_Kidney <- ggplot()+
#   geom_path(data = output.df, aes(x = age, y = CKP, color="Kidney blood"))+
#   geom_path(data = output.df, aes(x = age, y = CKT, color="Kidney tissue"))+
#   geom_path(data = output.df, aes(x = age, y = CFil, color="Filtrate"))+
#   theme(axis.text.x = element_text(size = 7),
#         axis.text.y = element_text(size = 7), 
#         axis.title = element_text(size = 8)) +
#   ylab("Kidney Compartments (ng/ml)") +
#   scale_color_manual(values = c("Kidney blood" = "lightblue", 
#                                 "Kidney tissue" = "blue", 
#                                 "Filtrate" = "darkblue"))
# 
# Plot_PFOA_Kidney
# ggsave("KidneyConcentration.2.png", dpi = 300)















# PBK code as by TRINE ####

# Special params
parms <- c(
  kfil = 0.1*Age_parms$GFR_M,#0.2*Age_parms$Q_kidney_M, #0.035, #Age_parms$GFR_M, # parameter from Trine; in the original code: kfil = 0.2*QK  # Clearance from the kidney to the filtrate compartment (L/h); 20% of bloodstream to QK is cleared for 
  Tm = 5842.308*24*BDW^0.75, # ug/d from ug/h/kg^0.75 parameter from Trine; transporter maximum
  Kt = 55, # ug/L; changed from Trine who had it as 55ug # Resorption affinity, changed from 0.055 in the original Loccisano 2011 model (ug)
  CLurine = 0.00000183*24*BDW^(-0.25)
)

PBPKmodPFOA_M <- function(t, state, parameters){
  with(as.list(c(state, parameters)), {
    
    # ## Dose ----
    Oraldose <- if_else(t <= exposure_stop, Oralconc * BDW, 0)
    Dermaldose <- if_else(t <= exposure_stop, AbsPFOA * Dermconc * BDW, 0)
    
    
    ## Concentrations ----
    
    # Organ concentrations (ug/L); these are TOTAL concentrations
    CP <- AP/VP  # Concentration in plasma (ug/L)
    CG <- AG/VG  # Concentration in gut (ug/L)
    CL <- AL/VL  # Concentration in liver (ug/L)
    CA <- AA/VA  # Concentration in adipose (ug/L)
    CR <- AR/VR  # Concentration in rest the body (ug/L)
    
    CSk <- ASk/VSk  # Concentration in skin (ug/L)
    # CSkb <- ASkb/VSkb # Concentration in skin barrier (ug/L) 
    # CSkP <- ASkP/VSkP # Concentration in skin plasma (ug/L) # For future use when updating skin again
    # CSkT <- ASkT/VSkT # Concentration in skin tissue (ug/L) # For future use when updating skin again
    
    CK <- AK/VK  # Concentration in kidney (ug/L)
    # CKP <- AKP/VKP # Concentration in kidney blood
    # CKT <- AKT/VKT # Concentration in kidney tissue
    CFil <- AFil/VFil # Concentration in filtrate compartment
    # CUr <- AUr/VUr # Concentration in urine (ug/L) # Not used
    
    CVG <- CG/PG # Concentration leaving gut (ug/L)
    CVL <- CL/PL  # Concentration leaving liver (ug/L)
    CVA <- CA/PA # Concentration leaving adipose (ug/L)
    CVSk <- CSk/PSk  # Concentration leaving skin (ug/L) 
    # CVSkP <- CSkP/PSk # Concentration leaving skin plasma (ug/L) # For future use when updating skin again
    CVR <- CR/PR  # Concentration leaving rest of the body (ug/L)
    
    CVK <- CK/PK # Concentration leaving the kidney (ug/L)
    # CVKP <- CKP/PK # Concentration leaving kidney blood (ug/L)
    
    
    
    ## Differential equations ----
    # adipose compartment
    dAA <- QA*(CP-CVA) # (ug/h)
    
    # Rest compartment
    dAR <- QR*(CP-CVR) # (ug/h)
    
    # Gut compartment: plasma to gut then to liver. Biliary clearance to gut, based on free concentration in the liver; faecal clearance based on total gut concentration; not only the free fraction as partitioning is not needed
    dAG <- Oraldose + QG*CP - QG*CVG + CLbiliary*CL*fup - CLfecal*CG #(ug/h)
    
    # Excretion fecal: cumulative
    dAEx_feces <- CLfecal*CG #ug/h, used to be CVG but I think this is wrong as it's not the tissue plasma partition that defines fecal excretion
    
    # Liver compartment
    dAL <- QL*CP + QG*CVG - (QL+QG)*CVL - CLbiliary*CL*fup #input from the hepatic artery
    # Rate of PFOA amount change in the liver (ug/h)
    
    ### OLD Kidney compartment ----
    # Kidney compartment
    dAK <- QK*(CP -CVK) + Tm*CFil/(Kt+CFil) - kfil*CK #from Trine's model ;Qfil*CK*Free was introduced to reflect clearance to filtrate compartment
    
    # Filtrate compartment
    dAFil <- kfil*(CK-CFil) - Tm*CFil/(Kt+CFil)
    
    # Storage compartment for urine
    dAdelay <- kfil*CFil - CLurine*Adelay # ug/h
    
    # Urine
    dAurine <- CLurine*Adelay # ug/h
    
    ### OLD Skin compartment ----
    dASk <- QSk*(CP-CVSk) + Dermaldose # Rate of PFOA amount change in skin
    
    ## Plasma compartment ----
    dAP <- - (QSk + QG + QL + QA + QK + QR)*CP +
      QSk*CVSk + (QL+QG)*CVL + QA*CVA + QK*CVK + QR*CVR #(ug/h) #QK*CVKP
    
    
    ## Mass Balance ----
    Atot <- AP + ASk + AG + AL + AA + AK + AFil + Adelay + Aurine + AR + AEx_feces # + AEx_urine + AKP + AKT + AFil + AUr + ASkb
    dInput <- Oraldose + Dermaldose
    MB = Input - Atot
    
    
    
    ## End ----
    
    list(c(dAP, 
           # dASkb, 
           dASk, 
           dAG, 
           dAL, 
           dAA, 
           dAK,
           dAFil,
           dAdelay,
           dAurine,
           # dAKP, 
           # dAKT, 
           # dAFil, 
           # dAUr, 
           dAR, 
           dAEx_feces, 
           dInput), #dAEx_urine
         c(CP = CP, 
           CSk = CSk, #CSkb = CSkb
           CG = CG, CVG = CVG, 
           CL = CL, CVL = CVL, 
           CA = CA, CVA = CVA,
           CK = CK, #CKP = CKP, CKT = CKT, 
           CFil = CFil, CVK = CVK, #CVKP = CVKP, 
           CR = CR, CVR = CVR,
           Atot = Atot,MB = MB)
    )
  })
}
