# Citation:
# DiRenzo, Graziella; Grant, Evan; Longo, Ana; Che-Castaldo, Christian; Zamudio, Kelly; Lips, Karen. 2017.  Imperfect pathogen detection from non-invasive skin swabs biases disease inference. Methods in Ecology & Evolution.

# This file:
  # Formats the data
  # Write the model ignoring imperfect sampling detection
  # Bundles the data
  # Runs the model 

###########################
Det <- read.table( "~/Dropbox/PhD_cope/Side_projects/Detection_probability_paper/Data/Cope_2012_2013.txt", sep = "\t", header = TRUE)

str(Det)

Swab <- Det[is.na(Det$LargestZooValue) == FALSE, ]

levels(Swab$Species)

#--------------------- Formating the data

stream <- numeric(nrow(Swab))
stream[grep("stream", Swab$type)] <- 1

trail <- numeric(nrow(Swab))
trail[grep("trail", Swab$type)] <- 1

wet <- numeric(nrow(Swab))
wet[grep("wet", Swab$season)] <- 1

dry <- numeric(nrow(Swab))
dry[grep("dry", Swab$season)] <- 1

raw <-cbind(Swab$LargestZooValue, Swab$SecondZoo)

w <- round(log(cbind(Swab$LargestZooValue, Swab$SecondZoo)+0.01), dig = 3)

yy <- round(log(cbind(Swab$LargestZooValue, Swab$SecondZoo)+0.01), dig = 3)

yy[raw > 0] <- 1
yy[raw == 0] <- 0


#### Model

sink("modelDSNA.txt")
cat("
    model { 
    
### Infection Intensity model

for (i in 1:n){

    mu[i] <- z[i] * mu.delta[i]

    mu.delta[i] <-  alpha_N_dry * dry[i] + 
                    alpha_N_wet * wet[i] +
                    alpha_N_trail * trail[i] +
                    alpha_N_stream * stream[i] 
                    
    
    x[i] ~ dnorm(mu[i], sd)

    x.new[i] ~ dnorm(mu[i], sd)
    
    }
    

    for (i in 1:n){
    
    z[i] ~ dbern(psi[i])
    
logit(psi[i]) <- alpha_psi_dry * dry[i] +
                 alpha_psi_wet * wet[i] +
                 alpha_psi_stream * stream[i] +
                 alpha_psi_trail * trail[i]

    }
    
### priors

alpha_N_dry ~ dnorm(0, 0.01)T(-10, 10)
alpha_N_wet ~ dnorm(0, 0.01)T(-10, 10)
alpha_N_stream ~ dnorm(0, 0.01)T(-10, 10)
alpha_N_trail ~ dnorm(0, 0.01)T(-10, 10)

alpha_psi_dry ~ dnorm(0, 0.368)T(-5, 5)
alpha_psi_wet ~ dnorm(0, 0.368)T(-5, 5)
alpha_psi_stream ~ dnorm(0, 0.368)T(-5, 5)
alpha_psi_trail ~ dnorm(0, 0.368)T(-5, 5)

tau <- 1/ (sd * sd)
sd ~ dgamma(0.1, 0.1)

#------- Bayesian Predictive check

    
for(i in 1:n){

      eval[i] <- abs(mu[i])

      E[i] <-  pow(pow(abs(x[i]), 0.5) - pow(eval[i], 0.5), 2)

      E.new[i] <- pow(pow(abs(x.new[i]), 0.5) - pow(eval[i], 0.5), 2)

} #is


zzzfit 		<- sum(E[]) 	
zzzfit.new  <- sum(E.new[])


#------- Derived quantities




}
",fill=TRUE)
sink()

# Collapse the data to ignore the double swabs
# Infection intensity
Nst <- apply(w, 1, max, na.rm = T)
Nst[Nst == -4.605] <- NA

# Presence/absence
zst <- apply(yy, 1, max, na.rm = T)
zst[zst > 0] <- 1

# Look at the data
round(cbind(Nst, zst, Swab$LargestZooValue), dig = 2)

Y <- rowSums(yy, na.rm = T)

missing <- numeric(length(zst))
missing[is.na(y[,2]) == TRUE] <- 1

# Calculate the total number of samples
samples <- rep(2, times = length(zst))
samples <- samples - missing

# Vector with idividual ID
n <- length(zst)
individual <- rep(1:n, each = 2)

#---------------- Bundling the data

inits = function() {
  list(
    alpha_N_wet = mean(Nst[Swab$season == "wet"], na.rm = T),
    alpha_N_dry = mean(Nst[Swab$season == "dry"], na.rm = T), 
    alpha_N_stream = mean(Nst[Swab$type == "stream"], na.rm = T),
    alpha_N_trail = mean(Nst[Swab$type == "trail"], na.rm = T),
    
    sd = runif(1, 0, 1),
    
    alpha_psi_dry= runif(1, 0, 1),
    alpha_psi_wet= runif(1, 0, 1),
    alpha_psi_stream= runif(1, 0, 1),
    alpha_psi_trail= runif(1, 0, 1)
  )}

params = c("alpha_N_dry",
           "alpha_N_wet",
           "alpha_N_trail",
           "alpha_N_stream",
           
           "alpha_psi_dry",
           "alpha_psi_wet",
           "alpha_psi_stream",
           "alpha_psi_trail",
           
           "sd", 
           
           "zzzfit",
           "zzzfit.new"
)

data.new = list(
  # Prevalence
  z = zst, 
  
  # Infection intensity
  x = Nst, 
  
  stream = stream,
  trail = trail,
  
  wet = wet,
  dry = dry,
  
  n = n		# total number of sites (R)
  
)

# MCMC settings
ni <- 50000
nb <- 10000
nt <- 50
nc <- 3			

library("jagsUI")

out24 <- jags(data.new, inits, params, "modelDSNA.txt", 
              n.chains = nc, n.thin = nt, n.iter = ni, 
              n.burnin = nb, parallel = TRUE)

# Set working directory
setwd("/Users/Cici/Desktop/")
# Save the model output
save(out24, file = "NOTAdjusted.rda")


