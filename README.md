A mechanistic PBK model for PFOA that can be used together with biomonitoring data with information regarding oral, dermal (or inhalation) exposure. 

Git-hub branches structure:
- **main** branch is the branch to be used together with HBM data
- **GSA** branch is the branch with all information about the global sensitivity analysis
- All other branches will be archived uppon publication
  
# Git-hub file structure:
## Input folder contains:
- **INPUT_dummy.csv**: file which is an example of how the input file from biomonitoring data should look like. Needed for running the _Do_PFOA_PBK_4Population.R_ file.
- **HalfLives.csv**: file which contains the observed half-lifes from human biomonitoring data, details regarding sex, number of subjects per study, reference and other study information are found in the file. Needed for compairing predicted with observed halflife in _RUN_and_OUTPUT.R_ and _PLOT_ObservedVSPredicted.R_ files.
- **PhysioVariables.csv**: file which contains the calculated physiological constants, based on the lifestage equations from Ratier et al. 2024. File calculated in _CALC_Lifestage_Constants.R_ file. Needed for calculating the PBK model input parameters in _CALC_Parameters.R_ file.
- **TissueComposition.csv**: file which contains the updated tissue composition (fractional volume of membrane lipids, albumine, structural proteins, water, fatty acids binding protein) per tissue. Needed for calculating the tissue-plasma partition coefficients in _CALC_Parameters.R_ file.

## Script folder contains:
- **DO_PFOA_PBK_4Population.R**: the only file that the user needs to use. The file should be used together with human biomonitoring data. It has a loop to run the PBK model and save the PBK predictions per person(subject).
- **DO_PFOA_PBK.R**: same file as above, but the user should manually add the relevant input data. The file only runs for one subject.
- **RUN_and_OUTPUT.R**: file which contains the code for calculating the input parameters, running the PBK model and analysing the results depending on exposure type
- **PBK_model.R**: file which contains the actual PFOA PBK model code for all three types of exposure
- **CALC_Parameters.R**: file which contains the code for calculating the input parameters for the PBK model
- **CALC_Lifestage_Constants.R**: file with all the lifestage calculations, to calculate the physiological constants (volumes and flows). For computation efficiency the file is actually not used every time but instead a csv file containing all the constants for all ages is used instead (_Input/PhysioVariables.csv_)
- **PLOT_ObservedVSPredicted.R**: file with the code for ploting observed vs predicted concentration over time curves and observed vs predicted halflife distributions.


# Requirements for use:

## Needed packages:

- here
- tidyverse
- deSolve
- PKNCA
- pracma

## Dependencies for running the DO_PFOA_PBK_4Population.R file:
The user should create and INPUT.csv file, in the **same structure** as the _INPUT_dummy.csv_ file.

- **Idcode**: (number) subject Id (used to link PBK results to the input information)
- **sex**: (M or F)subject sex
- **expBW**: (kg) if available the bodyweight of the subject at the time of exposure
- **expAGE**: (number in years) if available the age of the subject at the begining of exposure
- **exposure_type**: (characters) Oral if exposure was only oral, Dermal if exposure was only dermal, Inhalation if exposure was only via inhalation, Oral_Dermal if exposure was both oral and dermal simultaneously
- **expCONC**: (number in ug/kg/day) PFOA concentration at the begining of exposure, the model assumes that the exposure concentration is constant from the start until the end of exposure duration. In case of simultaneous oral and dermal exposures it should be the sum of those two. 
- **expCONC_Oral**: (number in ug/kg/day), only relevant in case exposure is modelled via both Oral and Dermal routes, if not, it can be left as NA or can contain the same value as expCONC, if exposure is Oral. PFOA concentration at the begining of exposure, the model assumes that the exposure concentration is constant from the start until the end of exposure duration.
- **expCONC_Dermal**: (number in ug/kg/day), only relevant in case exposure is modelled via both Oral and Dermal routes, if not, it can be left as NA or can contain the same value as expCONC, if exposure is Dermal. PFOA concentration at the begining of exposure, the model assumes that the exposure concentration is constant from the start until the end of exposure duration.
- **expSTOP**: (number in days), time after which the exposure should stop
- **Tinput**: (number in days), time that it takes for the exposure to happen, or where each exposure lasts (for example the dermal application of a make-up product is 8hours) (for oral exposure via water or food it's minimal and can be assumed 1)
- **tinterval**: (number in days), time between each exposure (for example the interval between each dermal application of a make-up product is 24h) (for oral exposure it can be assumed 1)

The user should also provide the simulation relevant information directly in the R file:

- **Tstart**: (time in days), time when the exposure started, default is 0
- **Tstop**: (time in days), time when the simulation should stop, can be longer or shorter than the exposure time
- **Dt**: (time in days), this is the iteration steps, default is 1, decreasing them would improve the accuracy of the model but will compromise computing speed

## Defaults, calculations and assumptions  
- Physiological constants (fractional organ volumes, blood flows, surface areas, GFR, kidney tubular flow, intestinal transit time) are selected based on the provided age, or bodyweight, in _CALC_Parameters.R_. The actual physiological parameters are calculated based on bodyweight in _PBK_model.R_.
- If expAGE is not provided, then the physiological constants are selected based on the provided bodyweight. Physiological parameters are also calculated based on the provided bodyweight.
- If expAGE is provided, then the physiological constants are selected based on the provided age. If expBW is also provided, then physiological parameters are calculated based on the actual bodyweight.
- If neither expAGE and expBW are provided then the default bodyweight of 70kg is used.
- If sex is not provided, then sex is assumed to be male (M). If sex is used then physiological constants are selected based on sex.
- The main model parameters are the same for all exposure scenarios. If dermal exposure is used, then the skin is added with the relevant parameters. If inhalation exposure is used, then the lungs are added and plasma is changed to arterial and venous plasma. The "Rest" organ and relevant parameters are updated accordingly.
- **Default initial PFOA amounts in each organ are 0**, if background PFOA organ concentrations are available, these can be directly added in the _RUN_and_OUTPUT.R_ file.
- The model assumes **plasma-limited perfusion**. As PFOA is highly bound to plasma proteins, the blood volume and flow rates were directly converted to plasma values, by correcting with the hematocrit.
- PFOA excretion to the bile is assumed to be the same as that of biliary acids.

# Results

Results are saved in an Output folder which is not synchronised in github. To change that update the git.ignore file.

# The PBK model

## Base structure
The basic PBK model specifies the organs that are necessary to describe the toxicokinetics of PFOA.

Entero-hepatic circulation is modelled as:

- Active/Passive uptake from the intestinal lumen to the intestine (vascular and cellular).
- Transfer to the liver extracellular space (vascular and interstitial space) via the portal plasma flow.
- Active, OATP transporter driven and albumin facilitated, uptake to hepatocytes.
- Active, BSEP mediated excretion to the intestinal lumen. [1]Physiologically, a fraction of billiary secretion is stored in the gallbladder and released after food consumption. For model simplification and as simulations are on the yearly scale, it was decided to assume direct excretion back to the intestinal lumen.
- 
![PFOA_PBK_model](https://github.com/user-attachments/assets/8ca3f246-8781-4aba-b35e-36149908a123)
