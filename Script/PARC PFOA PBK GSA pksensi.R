# --------------------------------------------------------------------------- #
# PBK MODEL FOR PFOA, TO BE USED TOGETHER WITH THE LATEST HBM DATA
# Global Sensitivity Analysis
# CP, 30-12-2024
# --------------------------------------------------------------------------- #


rm(list=ls()) # to clear out the global environment

# Set working directory

HOME <- "C:/Users/pacho003/OneDrive - Wageningen University & Research/C Channel/CP_L_R/PARC_PFAS_PBPKmodel/Codes_PFAS_models/Chrysanthi_Pachoulide_2024_PFOA_oral_dermal_blood"
# HOME <- "/home/westerj"
setwd(HOME)


# Set output storage directory
workingtime <- gsub(":", "-", Sys.time())
OUTPUT <- file.path("Output", Sys.Date(), workingtime)
dir.create(OUTPUT, recursive = TRUE)
setwd(OUTPUT)


# Load packages

library(ggplot2)
library(deSolve)
library(tidyverse)
library(data.table)
library(openxlsx)
library(writexl)
library(sensitivity)
library(pksensi)
library(PKNCA)


# PBK MODEL PARAMETERS ####
# ---------------------------------------------------------------------------- #

### Constants ####
#### Physiological  ####
Physio_params <- read_csv("C:/Users/pacho003/OneDrive - Wageningen University & Research/C Channel/CP_L_R/PARC_PFAS_PBPKmodel/Codes_PFAS_models/Chrysanthi_Pachoulide_2024_PFOA_oral_dermal_blood/Input/PhysioVariables.csv")

# 50 years old adult, male
Physio_params <- Physio_params %>% 
        filter(age == 50) %>% 
        select(ends_with('_M')) %>% 
        mutate(BloodFlowSum = rowSums(select(., starts_with("Q_")))) %>% # 0.9935; total blood flow as the sum of the fractional blood flows of all organs on which we have data
        mutate(VolumesSum = rowSums(select(., starts_with("V_")))) # 0.96; total volume as the sum of the fractional organ volumes of all organs on which we have data

BW <- Physio_params$BDW_M 
QC <- Physio_params$CardOut_M # This is corrected for hematocrit already so it's plasma


# Fractional organ volumes
VAc <- Physio_params$V_adiposeFraction_M   
MassAc <- Physio_params$AdiposeMass_M      #adipose mass
VGc <- Physio_params$V_gutFraction_M
VKc <- Physio_params$V_kidneyFraction_M
VLc <- Physio_params$V_liverFraction_M
VPc <- Physio_params$V_plasmaFraction_M
Hct <- Physio_params$Hct_M
VSkc <- Physio_params$V_skinFraction_M
VTotc <- Physio_params$VolumesSum 

# Fractional organ blood flows
QAc <- Physio_params$Q_adiposeFraction_M/Physio_params$BloodFlowSum
QGc <- Physio_params$Q_gutFraction_M/Physio_params$BloodFlowSum 
QKc <- Physio_params$Q_kidneyFraction_M/Physio_params$BloodFlowSum 
QLc <- Physio_params$Q_liverFraction_M/Physio_params$BloodFlowSum 
QSkc <- Physio_params$Q_skinFraction_M/Physio_params$BloodFlowSum 


#### Chemical Specific ####
PFOA_params <- read_csv("C:/Users/pacho003/OneDrive - Wageningen University & Research/C Channel/CP_L_R/PARC_PFAS_PBPKmodel/Codes_PFAS_models/Chrysanthi_Pachoulide_2024_PFOA_oral_dermal_blood/Input/PFOAParams.csv")

MW <- PFOA_params$MW 

fup <- PFOA_params$fup 

# Partition coefficients
# Calculate Plasma/Rest of the body partition coefficient
# KpRe = partition coefficient of each of the lumped organs * fractional volume of the respective organ / sum of the fractional volume of all these organs
KpRe <- (PFOA_params$KpBr*Physio_params$V_brainFraction_M +
                PFOA_params$KpHe*Physio_params$V_heartFraction_M + 
                PFOA_params$KpLu*Physio_params$V_lungFraction_M +
                PFOA_params$KpMu*Physio_params$V_muscleFraction_M +
                PFOA_params$KpSp*Physio_params$V_spleenFraction_M +
                PFOA_params$KpGo*Physio_params$V_reproFraction_M +
                PFOA_params$KpBo*Physio_params$V_boneFraction_M) / (
                        Physio_params$V_brainFraction_M +
                                Physio_params$V_heartFraction_M +
                                Physio_params$V_lungFraction_M +
                                Physio_params$V_muscleFraction_M +
                                Physio_params$V_spleenFraction_M +
                                Physio_params$V_reproFraction_M +
                                Physio_params$V_boneFraction_M)  
PFOA_params$KpRe <- KpRe

# Choose between calculated or initial ones (rat)
# In PFOA_params, PF, PG, PK, PL, PSK and PR are the initial PCs from Kudo 2007 (rat). These were recalculated to total tissue (/fup) to enable differentiating the fup for the sensitivity analysis
PAc <- PFOA_params$KpAd #PF # Adipose
PGc <- PFOA_params$KpGu #PG # Gut
PKc <- PFOA_params$KpKi #PK # Kidney
PLc <- PFOA_params$KpLi #PL # Liver
PSkc <- PFOA_params$KpSk #PSk # Skin
PRc <- PFOA_params$KpRe #PR #Rest


# Renal Clearance
CL_OAT4 <- PFOA_params$CL_OAT4 # L/d/kg protein; in vitro PFOA clearance in HEK239 cells overexpressing OAT4, Louisse et al. 2024
REF_OAT4 <- PFOA_params$REF_OAT4 # 1; as we don't have data on the in vitro expression of OAT4

# Biliary and Fecal Clearances
CLbiliaryc <- PFOA_params$CLbiliaryc #L/d/kg Fujii et al 2015 
CLfecalc <- PFOA_params$CLfaecesc # L/d/kg Fujii et al 2015 

### Final Parameters ####

#### Physiological  ####

# Body volumes (L)
VA <- VAc * BW/0.9 + MassAc/0.9
VG <- VGc * BW
VK <- VKc * BW
VL <- VLc * BW
VP <- VPc * (1-Hct) * BW
VSk <- VSkc * BW
VR <- (VTotc*BW) - (VA + VG + VK + VL + VP + VSk)

VFil <- VK * 0.05 # Kidney filtrate compartment, corresponds to the volume of the collecting system in [ICRP 89 page 149] http://www.icrp.org/publication.asp?id=ICRP%20Publication%2089

# Body flows (L/d)
QA <- QAc * QC
QG <- QGc * QC
QK <- QKc * QC
QL <- QLc * QC
QSk <- QSkc * QC
QR <- QC - (QA + QG + QK + QL + QSk)

QUr = 0.022 * BW # L/d, Urine flow rate to the bladder 22 mL/kg BW/d [ICRP 89 page 161]
GFR = 0.18 * QK # L/d, Glomerular filtration rate 18% of total renal plasma flow [ICRP 89 page 159] http://www.icrp.org/publication.asp?id=ICRP%20Publication%2089

#### Compound specific ####

# Partition coefficients
PA <- PAc * fup # Adipose
PG <- PGc * fup # Gut
PK <- PKc * fup # Kidney
PL <- PLc * fup # Liver
PSk <- PSkc * fup # Skin
PR <- PRc * fup #Rest

# Renal clearance
CL_FiltPT <- CL_OAT4 * REF_OAT4 * 0.17 * 0.7 * VK # L/d scaling for protein content instead of cell content; 17% of kidney is protein and 70% of kidney is cortex [ICRP 89], assuming that all the kidney protein is found in the cortex; this is an overestimation though; double ref for 17% protein Ruark 2020: DOI: https://doi.org/10.1016/B978-0-12-818596-4.00006-0

CLbiliary <- CLbiliaryc * BW #L/d  
CLfecal <- CLfecalc * BW # L/d 


parms <- unlist(c(data.frame(BW, #1
                             QC, #2
                             VA, #3
                             VG, #4
                             VK, #5
                             VL, #6
                             VP, #7
                             VSk, #8
                             VR, #9
                             VFil, #10
                             QA, #11
                             QG, #12 
                             QK, #13
                             QL, #14
                             QSk, #15
                             QR, #16
                             QUr, #17
                             GFR, #18
                             PA, #19
                             PG, #20
                             PK, #21
                             PL, #22
                             PSk, #23
                             PR, #24
                             CL_FiltPT, #25
                             CLbiliary, #26
                             CLfecal, #27
                             fup, #28
                             VAc, #29
                             VGc, #30
                             VKc, #31
                             VLc, #32
                             VPc, #33
                             VSkc, #34
                             QAc, #35
                             QGc, #36
                             QKc, #37
                             QLc, #38
                             QSkc, #39
                             PAc, #40
                             PGc, #41
                             PKc, #42
                             PLc, #43
                             PSkc, #44
                             PRc, #45
                             CL_OAT4, #46 
                             REF_OAT4, #47
                             CLbiliaryc, #48
                             CLfecalc #49
                             )))

parms


# EXPOSURE SCENARIO ####
# ------------------------------------------------------ #

exposure_stop <- 1 #50*365         # days
sim_stop <- (80*365) + 0.01 #150.01  # days

## Commment Chrysa 05-12-2024: I'm changing these for the sensitivity analysis. Does this affect the SA result?
TSTART <- 0.01
TSTOP <- sim_stop               # 24.01 years in days
DT <- 1                         # days
TIME <- seq(TSTART,TSTOP,by=DT)

Oraldose <- 3.96 # ug https://doi.org/10.1016/j.envint.2024.109047 # Oralconc*BDW # ug/day
# Dermaldose <- 0 # AbsPFOA*Dermconc*BDW # ug/day #AbsPFOA


# PBK MODEL ####
# ---------------------------------------------------------------------------- #

PBPKmodPFOA_M <- function(t, state, parameters){
        with(as.list(c(state, parameters)), {
                ## To make sure that volume and flow balances are correct
          
          QR <- QC - (QA + QG + QK + QL + QSk)
          VR <- (VTotc*BW) - (VA + VG + VK + VL + VP + VSk)
          
                ## Dose ----
                # Oraldose <- if_else(t <= exposure_stop, Oralconc * BDW, 0)
                # Dermaldose <- if_else(t <= exposure_stop, AbsPFOA * Dermconc * BDW, 0)
                
                ## Concentrations ----
                
                # Organ concentrations (ug/L); these are TOTAL concentrations
                CP <- AP/VP  # Concentration in plasma (ug/L)
                CG <- AG/VG  # Concentration in gut (ug/L)
                CL <- AL/VL  # Concentration in liver (ug/L)
                CA <- AA/VA  # Concentration in adipose (ug/L)
                CR <- AR/VR  # Concentration in rest the body (ug/L)
                
                CSk <- ASk/VSk  # Concentration in skin (ug/L)
                
                CK <- AK/VK  # Concentration in kidney (ug/L)
                CFil <- AFil/VFil # Concentration in filtrate compartment
                
                CVG <- CG/PG # Concentration leaving gut (ug/L)
                CVL <- CL/PL  # Concentration leaving liver (ug/L)
                CVA <- CA/PA # Concentration leaving adipose (ug/L)
                CVSk <- CSk/PSk  # Concentration leaving skin (ug/L) 
                CVR <- CR/PR  # Concentration leaving rest of the body (ug/L)
                
                CVK <- CK/PK # Concentration leaving the kidney (ug/L)
                
                
                ## Differential equations ----
                
                # Adipose compartment
                dAA <- QA*(CP-CVA) # (ug/h)
                
                # Rest compartment
                dAR <- QR*(CP-CVR) # (ug/h)
                
                # Gut compartment: 
                dAG <- + QG*CP - QG*CVG + CLbiliary*CL*fup - CLfecal*CG #Oraldose
                
                # Excretion fecal: cumulative
                dAEx_feces <- CLfecal*CG 
                
                # Liver compartment
                dAL <- QL*CP + QG*CVG - (QL+QG)*CVL - CLbiliary*CL*fup
                
                # Kidney compartment
                dAK <- QK*(CP -CVK) - GFR*CK*fup + CL_FiltPT*CFil*fup
                
                # Filtrate compartment
                dAFil <- GFR*CK*fup - QUr*CFil - CL_FiltPT*CFil*fup 
                
                
                # Urine compartment
                dAUr <- QUr*CFil
                
                # Skin compartment
                dASk <- QSk*(CP-CVSk) #+ Dermaldose 
                
                # Plasma compartment
                dAP <- - (QSk + QG + QL + QA + QK + QR)*CP +
                        QSk*CVSk + (QL+QG)*CVL + QA*CVA + QK*CVK + QR*CVR 
                
                
                # Mass Balance
                Atot <- AP + ASk + AG + AL + AA + AK + AFil + AR + AEx_feces + AUr #AKP + AKT
                # # dInput <- Oraldose + Dermaldose
                # # MB = Input - Atot
                MB = Oraldose - Atot
                
                # End
                
                list(c(dAP, 
                       dASk, 
                       dAG, 
                       dAL, 
                       dAA,
                       dAK,
                       dAFil,
                       dAUr, 
                       dAR, 
                       dAEx_feces #, 
                       # dInput
                ), 
                c(CP = CP, 
                  CSk = CSk, 
                  CG = CG, CVG = CVG, 
                  CL = CL, CVL = CVL, 
                  CA = CA, CVA = CVA,
                  CK = CK,
                  CFil = CFil, 
                  CVK = CVK,
                  CR = CR, CVR = CVR,
                  Atot = Atot, MB = MB
                  )
                )
        })
}



## Initials ####


A_init <- unlist(c(data.frame(
        AP = 0, 
        ASk = 0, # Dermaldose, 
        AG = Oraldose,
        AL = 0,
        AA = 0,
        AK = 0,
        AFil = 0,
        AUr = 0,
        AR = 0,
        AEx_feces = 0 #,
        # Input = 0
)))


## Solving the model ####
output_PFOA <- ode(y = A_init,
                   times = TIME,
                   func = PBPKmodPFOA_M, 
                   parms = parms,
                   method="lsoda"
                   )
output.PFOA.df <- as.data.frame(output_PFOA) 
plot(output.PFOA.df)

# Sensitivity analysis using pksensi ####
# ---------------------------------------------------------------------------- #

outputs = c("CP") # Variable in test uncertainty/sensitivity


## Define the distribution of the parameters that you will analyse in the sensitivity test
q <- c( "qunif", #1
        "qunif", #2
        "qunif", #3
        "qunif", #4
        "qunif", #5
        "qunif", #6
        "qunif", #7
        "qunif", #8
        "qunif", #9
        "qunif", #10,
        "qunif", #11, 
        "qunif", #12
        "qunif", #13
        "qunif", #14
        "qunif", #15
        "qunif"#, #16
        # "qunif", #17
        # "qunif", #18
        # "qunif", #19
        # "qunif", #20
        # "qunif", #21
        # "qunif", #22
        # "qunif", #23
        # "qunif", #24
        # "qunif", #25
        # "qunif", #26
        # "qunif", #27
        # "qunif", #28
        # "qunif"#, #29
        # "qunif", #30
        # "qunif", #31
        # "qunif", #32
        # "qunif", #33
        # "qunif", #34
        # "qunif", #35
        # "qunif", #36
        # "qunif", #37
        # "qunif", #38
        # "qunif", #39
        # "qunif", #40
        # "qunif", #41
        # "qunif", #42
        # "qunif", #43
        # "qunif", #44
        # "qunif", #45
        # "qunif", #46 
        # "qunif", #47
        # "qunif", #48
        # "qunif" #49
)

## Set parameter distribution ##
# we use 10% change in all parameters

LL <- 0.9 # 10% lower limit
UL <- 1.1 # 10% upper limit


q.arg <- list(list(min = parms["BW"]*LL, max= parms["BW"]*UL), #1
              # list(min = parms["QC"]*LL, max= parms["QC"]*UL), #2
              # list(min = parms["VA"]*LL, max = parms["VA"]*UL), #3
              # list(min = parms["VG"]*LL, max = parms["VG"]*UL), #4
              list(min = parms["VK"]*LL, max = parms["VK"]*UL), #5
              # list(min = parms["VL"]*LL, max = parms["VL"]*UL), #6
              # list(min = parms["VP"]*LL, max = parms["VP"]*UL), #7
              # list(min = parms["VSk"]*LL, max = parms["VSk"]*UL), #8
              # list(min = parms["VR"]*LL, max = parms["VR"]*UL), #9
              list(min = parms["VFil"]*LL, max = parms["VFil"]*UL), #10
              # list(min = parms["QA"]*LL, max = parms["QA"]*UL), #11
              # list(min = parms["QG"]*LL, max = parms["QG"]*UL), #12
              list(min = parms["QK"]*LL, max = parms["QK"]*UL), #13
              # list(min = parms["QL"]*LL, max = parms["QL"]*UL), #14
              # list(min = parms["QSk"]*LL, max = parms["QSk"]*UL), #15
              # list(min = parms["QR"]*LL, max = parms["QR"]*UL), #16
              list(min = parms["QUr"]*LL, max = parms["QUr"]*UL), #17
              list(min = parms["GFR"]*LL, max = parms["GFR"]*UL), #18
              list(min = parms["PA"]*LL, max = parms["PA"]*UL), #19
              list(min = parms["PG"]*LL, max = parms["PG"]*UL), #20
              list(min = parms["PK"]*LL, max = parms["PK"]*UL), #21
              list(min = parms["PL"]*LL, max = parms["PL"]*UL), #22
              list(min = parms["PSk"]*LL, max = parms["PSk"]*UL), #23
              list(min = parms["PR"]*LL, max = parms["PR"]*UL), #24
              list(min = parms["CL_FiltPT"]*LL, max = parms["CL_FiltPT"]*UL), #25
              list(min = parms["CLbiliary"]*LL, max = parms["CLbiliary"]*UL), #26
              list(min = parms["CLfecal"]*LL, max = parms["CLfecal"]*UL), #27
              list(min = parms["fup"]*LL, max = parms["fup"]*UL)#28
              # list(min = parms["VAc"]*LL, max = parms["VAc"]*UL), #29
              # list(min = parms["VGc"]*LL, max = parms["VGc"]*UL), #30
              # list(min = parms["VKc"]*LL, max = parms["VKc"]*UL), #31
              # list(min = parms["VLc"]*LL, max = parms["VLc"]*UL), #32
              # list(min = parms["VPc"]*LL, max = parms["VPc"]*UL), #33
              # list(min = parms["VSkc"]*LL, max = parms["VSkc"]*UL), #34
              # list(min = parms["QAc"]*LL, max = parms["QAc"]*UL), #35
              # list(min = parms["QGc"]*LL, max = parms["QGc"]*UL), #36
              # list(min = parms["QKc"]*LL, max = parms["QKc"]*UL), #37
              # list(min = parms["QLc"]*LL, max = parms["QLc"]*UL), #38
              # list(min = parms["QSkc"]*LL, max = parms["QSkc"]*UL), #39
              # list(min = parms["PAc"]*LL, max = parms["PAc"]*UL), #40
              # list(min = parms["PGc"]*LL, max = parms["PGc"]*UL), #41
              # list(min = parms["PKc"]*LL, max = parms["PKc"]*UL), #42
              # list(min = parms["PLc"]*LL, max = parms["PLc"]*UL), #43
              # list(min = parms["PSkc"]*LL, max = parms["PSkc"]*UL), #44
              # list(min = parms["PRc"]*LL, max = parms["PRc"]*UL), #45
              # list(min = parms["CL_OAT4"]*LL, max = parms["CL_OAT4"]*UL), #46
              # list(min = parms["REF_OAT4"]*LL, max = parms["REF_OAT4"]*UL), #47
              # list(min = parms["CLbiliaryc"]*LL, max = parms["CLbiliaryc"]*UL), #48
              # list(min = parms["CLfecalc"]*LL, max = parms["CLfecalc"]*UL) #49
              # 
)


## Create parameter matrix ##
set.seed(1234)

params <- c(
        "BW", #1
        # "QC", #2
        # "VA", #3
        # "VG", #4
        "VK", #5
        # "VL", #6
        # "VP", #7
        # "VSk", #8
        # "VR", #9
        "VFil", #10
        # "QA", #11
        # "QG", #12 
        "QK", #13
        # "QL", #14
        # "QSk", #15
        # "QR", #16
        "QUr", #17
        "GFR", #18
        "PA", #19
        "PG", #20
        "PK", #21
        "PL", #22
        "PSk", #23
        "PR", #24
        "CL_FiltPT", #25
        "CLbiliary", #26
        "CLfecal", #27
        "fup"#, #28
        # "VAc", #29
        # "VGc", #30
        # "VKc", #31
        # "VLc", #32
        # "VPc", #33
        # "VSkc", #34
        # "QAc", #35
        # "QGc", #36
        # "QKc", #37
        # "QLc", #38
        # "QSkc"#, #39
        # "PAc", #40
        # "PGc", #41
        # "PKc", #42
        # "PLc", #43
        # "PSkc", #44
        # "PRc", #45
        # "CL_OAT4", #46 
        # "REF_OAT4", #47
        # "CLbiliaryc", #48
        # "CLfecalc" #49
        )

length(params) == length(q)

# Create the sequences for each parameter by eFAST
x <- rfast99(params = params, n = 200, q = q, q.arg = q.arg, replicate = )

dim(x$a) # the array of c(model evaluation, replication, parameters)

## Conduct simulation ##
out <- solve_fun(x,
                 time = TIME,
                 func = PBPKmodPFOA_M,
                 initState = A_init,
                 outnames = outputs)

saveRDS(object=out, file="out_scaled.rds")
out <- readRDS("out_scaled.rds")

# Output of the Uncertainty analysis ##
pdf("out_scaled.pdf")
pksim(out) # PK plot of the outputs based on the given parameter (Uncertainty analysis)
dev.off()

# Output from the sensitivity analysis ##
pdf("out_scaled.pdf")
plot(out) # plot Time-dependent sensitivity (with 95 % CI)


dev.off()

output_PFOASI <- as.data.frame(print(out["tSI"]))
output_PFOASI$Times <- rownames(output_PFOASI)
output_PFOASI$Times <- as.numeric(as.character(output_PFOASI$Times))
output_PFOASI$Year <- output_PFOASI$Times/(365*24)

write.csv(output_PFOASI, "output_PFOASI.csv", col.names = TRUE)

check(out) # Check sensitivity measurement for parameter fixing
pdf("heat_check_CI_scaled.pdf")
heat_check(out, index = "CI") # heat_check Create heatmap to overview the result of GSA
dev.off()

lkjwqpdf("heat_check_all.pdf")
heat_check(out, show.all = TRUE)
dev.off()


