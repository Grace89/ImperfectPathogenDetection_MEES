# Citation:
# DiRenzo, Graziella; Grant, Evan; Longo, Ana; Che-Castaldo, Christian; Zamudio, Kelly; Lips, Karen. 2017.  Imperfect pathogen detection from non-invasive skin swabs biases disease inference. Methods in Ecology & Evolution.

# This file:
  # Compares the parameter estimates between the imperfect pathogen detection model and the unadjusted model

# Load libraries
library(coda)

# Load data
load("adjusted.rda")

load("NOTadjusted.rda")

#--- Prevalence dry vs wet
mean(out24$sims.list$alpha_psi_dry < out24$sims.list$alpha_psi_wet)
mean(out22$sims.list$alpha_psi_dry < out22$sims.list$alpha_psi_wet)

mean(out24$sims.list$alpha_psi_stream > out24$sims.list$alpha_psi_trail)
mean(out22$sims.list$alpha_psi_stream > out22$sims.list$alpha_psi_trail)


mean(out24$sims.list$alpha_psi_wet < out22$sims.list$alpha_psi_wet)
mean(out24$sims.list$alpha_psi_dry < out22$sims.list$alpha_psi_dry)
mean(out24$sims.list$alpha_psi_stream < out22$sims.list$alpha_psi_stream)
mean(out24$sims.list$alpha_psi_trail < out22$sims.list$alpha_psi_trail)

#---- Infection intensity

#---- Wet vs Dry
mean((out24$sims.list$alpha_N_wet) > out24$sims.list$alpha_N_dry)
mean((out22$sims.list$alpha_N_wet) > (out22$sims.list$alpha_N_dry))

mean((out24$sims.list$alpha_N_dry) > (out22$sims.list$alpha_N_dry))
mean((out24$sims.list$alpha_N_wet) > (out22$sims.list$alpha_N_wet))


#--- Stream vs trail
mean((out24$sims.list$alpha_N_stream) > (out24$sims.list$alpha_N_trail))
mean((out22$sims.list$alpha_N_stream) > (out22$sims.list$alpha_N_trail))

mean((out24$sims.list$alpha_N_stream) > (out22$sims.list$alpha_N_stream))
mean((out24$sims.list$alpha_N_trail) > (out22$sims.list$alpha_N_trail))
