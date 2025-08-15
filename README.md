A mechanistic PBK model for PFOA that can be used together with biomonitoring data with information regarding oral, dermal (or inhalation) exposure. 

Git-hub branches structure:
- **main** branch is the branch to be used together with HBM data
- All other branches will be archived uppon publication
  
# Git-hub file structure:
## Input folder contains:
- **INPUT_dummy.csv**: file which is an example of how the input file from biomonitoring data should look like. Needed for running the _DO_PFOA_PBK.R_ file. In the file, Idcode should be a numeric value (it's only used to link PBK results of that number to the relvant input information from this file)
- **HalfLives.csv**: file which contains the observed half-lifes from human biomonitoring data, details regarding sex, number of subjects per study, reference and other study information are found in the file. Needed for compairing predicted with observed halflife in _RUN_and_OUTPUT.R_ and _PLOT_ObservedVSPredicted.R_ files.
- **PhysioVariables.csv**: file which contains the calculated physiological constants, based on the lifestage equations from Ratier et al. 2024. File calculated in _CALC_Lifestage_Constants.R_ file. Needed for calculating the PBK model input parameters in _CALC_Parameters.R_ file.
- **TissueComposition.csv**: file which contains the updated tissue composition (fractional volume of membrane lipids, albumine, structural proteins, water, fatty acids binding protein) per tissue. Needed for calculating the tissue-plasma partition coefficients in _CALC_Parameters.R_ file.

## Script folder contains:
- **DO_PFOA_PBK.R**: **the only file that a user needs**
- **RUN_and_OUTPUT.R**: file which contains the code for calculating the input parameters, running the PBK model and analysing the results depending on exposure type
- **PBK_model.R**: file which contains the actual PFOA PBK model code for all three types of exposure
- **CALC_Parameters.R**: file which contains the code for calculating the input parameters for the PBK model
- **CALC_Lifestage_Constants.R**: file with all the lifestage calculations, to calculate the physiological constants (volumes and flows). For computation efficiency the file is actually not used every time but instead a csv file containing all the constants for all ages is used instead (_Input/PhysioVariables.csv_)
- **PBK_Results_Report.qmd**: file that generates the final report containing the results of the simulation

# Dependencies for running the DO_PFOA_PBK.R file:

## Needed packages:

1. here
2. tidyverse
3. deSolve
4. PKNCA
5. pracma
6. glue
7. patchwork
8. quarto
9. showtext
10. DT
11. plotly
12. tinytex
13. webshot2

## Inputs:

The user should choose between:
- **Lifestage > "Yes":** to re-calculate human physiological input data per simulated year. **The data is re-calculated based on age, therefore expAGE is a mandatory input to run this option**
  - **Lifestage > "No":** indicates that the model assumes the same physiology throughout the whole simulation time. (For example, for a simulated person that is 30 years old with a bodyweight of 60kg, if the simulation is run for 50 year it is assumed that the bodyweight stays at 60kg during the 50 years.)
    
- **Population > "Yes":** if multuple people are run (for example in the case of running exposure data from HBM studies)
  - **Population > "No":** indicates that the user will input data related to a one person exposure

- **Test_study**: manaually add the contents of the "Study" column (column 18) from the input_file (used to filter the relevant input information), or add a name to the test.

Required input to run the model, provided in the *INPUT_dummy.csv* file, or mannualy in the document if a one-person simulation is run:
- **sex**: (M or F)subject sex
- **expBW**: (kg) if available the bodyweight of the subject at the time of exposure
- **expAGE**: (number in years) if available the age of the subject at the begining of exposure
- **exposure_type**: (characters) Oral if exposure was only oral, Dermal if exposure was only dermal, Inhalation if exposure was only via inhalation, Oral_Dermal if exposure was both oral and dermal simultaneously
- **exp**: (number in ug/kg/day) PFOA concentration at the begining of exposure, the model assumes that the exposure concentration is constant from the start until the end of exposure duration. In case of simultaneous oral and dermal exposures it should be the sum of those two. 
- **exp_Oral**: (number in ug/kg/day), only relevant in case exposure is modelled via both Oral and Dermal routes, if not, it can be left as NA or can contain the same value as expCONC, if exposure is Oral. PFOA concentration at the begining of exposure, the model assumes that the exposure concentration is constant from the start until the end of exposure duration.
- **exp_Dermal**: (number in ug/kg/day), only relevant in case exposure is modelled via both Oral and Dermal routes, if not, it can be left as NA or can contain the same value as expCONC, if exposure is Dermal. PFOA concentration at the begining of exposure, the model assumes that the exposure concentration is constant from the start until the end of exposure duration.
- **expSTOP**: (number in days), time after which the exposure should stop
- **Tinput**: (number in days), time that it takes for the exposure to happen, or where each exposure lasts (for example the dermal application of a make-up product is 8hours) (for oral exposure via water or food it's minimal and can be assumed 1)
- **tinterval**: (number in days), time between each exposure (for example the interval between each dermal application of a make-up product is 24h) (for oral exposure it can be assumed 1)
- **Tstart**: (time in days), time when the exposure started, default is 0
- **Tstop**: (time in days), time when the **simulation** should stop, should be at least equal to "expSTOP", can be longer
- **Dt**: (time in days), this is the iteration steps, default is 1, decreasing them would improve the accuracy of the model but will compromise computing speed
- **samplingT**: (time in days), time at which the observed serum samples were taken, usually equal to either "expSTOP" or "Tstop", cannot be longer than Tstop, used for plotting observed vs predicted _(not needed for running the one-person simulation)_
- **CP_measured**: (in ug/ml), measured plasma/serum concentration, optional, used for plotting observed vs predicted _(not needed for running the one-person simulation)_
- **HL_measured**: (in years), half life determined in the reference study, optional, used for plotting observed vs predicted _(not needed for running the one-person simulation)_
- **Study**: (character), an indicative study name that will be used to filter out the test population

## Defaults, calculations and assumptions  
- To run with the lifestage option, exposure age (expAGE) is requred. Noting also that the bodyweight is predicted based on the age so even a bodyweight is provided it will be over-run by the predicted one. 
- Physiological constants (fractional organ volumes, blood flows, surface areas, GFR, kidney tubular flow, intestinal transit time) are selected based on the provided age, or bodyweight, in _CALC_Parameters.R_. The actual physiological parameters are calculated based on bodyweight in _PBK_model.R_.
- If expAGE is not provided, then the physiological constants are selected based on the provided bodyweight. Physiological parameters are also calculated based on the provided bodyweight.
- If expAGE is provided, then the physiological constants are selected based on the provided age. If expBW is also provided, then physiological parameters are calculated based on the actual bodyweight.
- If neither expAGE and expBW are provided then the default bodyweight of 70kg is used.
- If sex is not provided, then sex is assumed to be male (M). If sex is used then physiological constants are selected based on sex.
- The main model parameters are the same for all exposure scenarios. If dermal exposure is used, then the skin is added with the relevant parameters. If inhalation exposure is used, then the lungs are added and plasma is changed to arterial and venous plasma. The "Rest" organ and relevant parameters are updated accordingly.
- **Default initial PFOA amounts in each organ are 0**, if background PFOA organ concentrations are available, these can be directly added in the _RUN_and_OUTPUT.R_ file.
- The model assumes **plasma-limited perfusion**. As PFOA is highly bound to plasma proteins, the blood volume and flow rates were directly converted to plasma values, by correcting with the hematocrit.
- Passive PFOA kinetics are driven from the tissue-plasma partition coefficients and PFOA fraction unbound in plasma. Tissue-plasma partition coefficients are calculated for all organs based on PFOA distribution coefficients to organ components (Allendorf et al. 2021) and the fractional volume of each component in each organ (found in _TissueComposition.csv_). 
- PFOA excretion to the bile is assumed to be the same as that of biliary acids.

# Results
Results are saved in an Output folder which is not synchronised in github. To change that update the git.ignore file.
An html file is generated automatically to show a summary of the results and the plots.

# The PBK model

## Base structure
The basic PBK model specifies the organs that are necessary to describe the toxicokinetics of PFOA.

(Semi-permeability limited liver and intestine) Entero-hepatic circulation is modelled as:
- Passive/active uptake from the intestinal lumen to the intestine (vascular and cellular).
- Transfer to the liver extracellular space (vascular and interstitial space) via the portal plasma flow.
- Active, OATP transporter mediated and albumin facilitated, uptake to hepatocytes.
- Active, BSEP mediated excretion to the intestinal lumen from the hepatocytes. Physiologically, a fraction of billiary secretion is stored in the gallbladder and released after food consumption. For model simplification and as simulations are on the yearly scale, it was decided to assume direct excretion back to the intestinal lumen.
- Fecal excretion based on physiological colonic transit time of luminal contents.

(Sequential, semi-permeability limited kidney) Renal secretion and re-absorption is modelled as:
- Glomerular filtration of the fraction unbound of PFOA in plasma to the proximal tubuly lumen, forming the primary urine.
- Active, OAT4 transported mediated and albumin facilitated, re-uptake to the proximal tubule (vascular and cellular). 
- Passive (kidney plasma flow driven) flow of PFOA from the proximal tubule to the rest of the kidney structures (loop of Henle, Distal tubule, Collecting duct), and then back to the systemic circulation.
- Passive (tubular fluid flow driven) flow of PFOA from the proximal tubule lumen to the lumen of the rest of the kidney structures.
- Renal excretion, from the lumen of the rest of the kindey structures to the urine, based on the physiological urinary flow rate.

PFOA excretion via menstrual plasma clearance is also simulated, as well as serum albumin variation throughout lifestage.

Given the affinity of PFOA to different type of lipids, a simple adipose compartment is also added. All other organs are lumped together in a "rest" compartment.
<img width="3022" height="2488" alt="Schematic_representation_PBK_model" src="https://github.com/user-attachments/assets/55afc791-05ac-4b71-9bd6-ae763b5c77a5" />

## Oral exposure
The PBK model for oral exposure is exactly that of the base structure. The initial amount of PFOA from oral dose is added to the intestinal lumen.

## Dermal exposure
The PBK model for dermal exposure has two additional compartments: skin and skin barrier (epidermis).

- PFOA absorption from the skin is modelled as a clearance from the epidermis to the skin (vascular and cellular).
- Distribution to the rest of the body from the skin is plasma flow driven.

## Oral and dermal exposure
As described above

## Inhalation
A lung compartment is added. Plasma is devided in arterial and venous. 
Exposure is directly added to the lung compartment. In future versions a more mechanistic description of the inhalation exposure will be added, based on available mechanistic studies. 

