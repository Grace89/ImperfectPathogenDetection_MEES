The data used in the analyses and figures can be found on the dryad data repository searching for:
Imperfect pathogen detection from non-invasive skin swabs biases disease inference
DiRenzo, Graziella; Grant, Evan; Longo, Ana; Che-Castaldo, Christian; Zamudio, Kelly; Lips, Karen


Model> Model code >
Double_Swab_Imperfect_sampling.R
  # Formats the data
  # Write the imperfect sampling detection model
  # Bundles the data
  # Runs the model 
Double_Swab_IGNORING_imperfect_sampling.R
  # Formats the data
  # Write the imperfect sampling detection model
  # Bundles the data
  # Runs the model 


Model> Model code > Model output
adjusted.rda 
   # JAGS output for model with imperfect sampling detection
NOTadjusted.rda 
   # JAGS output for model IGNORING imperfect sampling detection
Miller.rda
   # JAGS output for model replicating model output for Miller et al. 2012


Model>  Simulate Miller et al. 2012 data
Simulate_Miller_2012.R
  #This file simulates the data for Miller et al. 2012, runs the model, and saves the model output


Figures
Bayesianpvalue.R
  # Creates a figure to visualize the Bayesian p-value 
  # Calculates the Bayesian p-value 
Fig1_detection_Figure.R
  # Creates Figure 1 in the main text
Fig2_Pstar_combo.R
  # Creates Figure 2 in the main text


Parameter comparisons
Comparisons.R
  # Compares the parameter estimates between the imperfect pathogen detection model and the unadjusted model



