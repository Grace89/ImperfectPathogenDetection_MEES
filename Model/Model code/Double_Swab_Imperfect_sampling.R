# Citation:
# DiRenzo, Graziella; Grant, Evan; Longo, Ana; Che-Castaldo, Christian; Zamudio, Kelly; Lips, Karen. 2017.  Imperfect pathogen detection from non-invasive skin swabs biases disease inference. Methods in Ecology & Evolution.

# This file:
  # Formats the data
  # Write the imperfect sampling detection model
  # Bundles the data
  # Runs the model 

###########################
# Load the data
Det <- read.table( "~/Dropbox/PhD_cope/Side_projects/Detection_probability_paper/Data/Cope 2012-2013.txt", sep = "\t", header = TRUE)

Swab <- Det[is.na(Det$LargestZooValue) == FALSE, ]

#--------------------- Formating the data

wet_stream <- numeric(nrow(Swab))
wet_stream[which(Swab$season== "wet" & Swab$type == "stream") ] <- 1

wet_trail <- numeric(nrow(Swab))
wet_trail[which(Swab$season== "wet" & Swab$type == "trail") ] <- 1

dry_stream <- numeric(nrow(Swab))
dry_stream[which(Swab$season== "dry" & Swab$type == "stream") ] <- 1

dry_trail <- numeric(nrow(Swab))
dry_trail[which(Swab$season== "dry" & Swab$type == "trail") ] <- 1

#--- Observation of infection intensity
w <- cbind(Swab$LargestZooValue, Swab$SecondZoo)
wNEW <- c(t((w)))
wNEW[wNEW == 0] <- NA
wNEW <- log(wNEW)

#--- Observed infection
yy <- cbind(Swab$LargestZooValue, Swab$SecondZoo)

yy[yy > 0] <- 1
yy[yy == 0] <- 0


#### Model

sink("modelDS.txt")
cat("
    model { 
    
#------------- Infection intensity
  
for (i in 1:n){

    mu.delta[i] <-  alpha_N_dry_trail *  dry_trail[i] + 
              alpha_N_dry_stream * dry_stream[i] + 
              alpha_N_wet_trail *  wet_trail[i] + 
              alpha_N_wet_stream * wet_stream[i]
    
    x[i] ~ dnorm(mu[i], tau)

    mu[i] <- z[i] * mu.delta[i] + 0.001

}


#--------------- Prevalence

for (i in 1:n){

    z[i] ~ dbern(psi[i])
    
    logit(psi[i]) <- alpha_psi_dry_stream * dry_stream[i] + 
                     alpha_psi_wet_stream * wet_stream[i] + 
                     alpha_psi_dry_trail * dry_trail[i] + 
                     alpha_psi_wet_trail * wet_trail[i]

}


#------------------ Detection & Measurement Error

for (i in 1:n){

  for(j in 1:D[i]){

    logit(p[i,j]) <- alpha + x[i]* beta
    
    P[i, j] <- p[i, j] * z[i]
    
    y[i, j] ~ dbern(P[i, j])
  }
    
}

for (i in 1:R){

  w[i] ~ dnorm(x[individual[i]], tau.error)

  w.new[i] ~ dnorm(x[individual[i]], tau.error)

}


#------------------ Priors

alpha_N_dry_stream ~ dnorm(0, 0.01)T(-10, 10)
alpha_N_wet_stream ~ dnorm(0, 0.01)T(-10, 10)
alpha_N_dry_trail ~ dnorm(0, 0.01)T(-10, 10)
alpha_N_wet_trail ~ dnorm(0, 0.01)T(-10, 10)

tau <- 1/ (sd * sd)
tau.error <- 1/ (sd.error * sd.error)

sd.error ~ dgamma(0.1, 0.1)
sd ~ dgamma(0.1, 0.1)

alpha_psi_dry_stream ~ dnorm(0, 0.368)T(-5, 5)
alpha_psi_wet_stream ~ dnorm(0, 0.368)T(-5, 5)
alpha_psi_dry_trail ~ dnorm(0, 0.368)T(-5, 5)
alpha_psi_wet_trail ~ dnorm(0, 0.368)T(-5, 5)

alpha ~ dnorm(0, 0.368)T(-5, 5)
beta ~ dnorm(0, 0.368)T(-5, 5)


#-------------------- Bayesian Posterior Predictive check

    
for(i in 1:n){
      eval[i] <- abs(x[i])
}  

for(i in 1:R){

E[i] <- pow(pow(abs(w[individual[i]]), 0.5) - pow(eval[individual[i]], .5), 2)


E.new[i] <- pow(pow(abs(w.new[individual[i]]), 0.5) - pow(eval[individual[i]], 0.5), 2)

} #is

zzzfit 		<- sum(E[]) 	
zzzfit.new  <- sum(E.new[])

}
",fill=TRUE)
sink()


#------ Initial values
Nst <- apply(log(w+0.001), 1, max, na.rm = TRUE)+1

zst <- apply(yy, 1, max, na.rm = TRUE)

missing <- numeric(length(zst))
missing[is.na(yy[,2]) == TRUE] <- 1

samples <- rep(2, times = length(zst))
samples <- samples - missing

n <- length(zst)

individual <- rep(1:n, each = 2)


#---------------- Bundling the data

inits = function() {
  list(
    # Prevalence
    z = zst, 
    
    # Infection intensity
    x = Nst, 
    
    alpha_N_wet_stream = mean(Nst[Swab$season == "wet" & Swab$type == "stream"], na.rm = TRUE)+1,
    alpha_N_dry_stream = mean(Nst[Swab$season == "dry" & Swab$type == "stream"], na.rm = TRUE)+1, 
    alpha_N_dry_trail = mean(Nst[Swab$season == "dry" &Swab$type == "trail"], na.rm = TRUE)+1,
    alpha_N_wet_trail = mean(Nst[Swab$season == "wet" & Swab$type == "trail"], na.rm = TRUE)+1,
    
    # Detection
    alpha = runif(1, 1, 2), 
    beta = runif(1, 1, 2),
    
    sd = runif(1, 0.2, 0.5),
    sd.error = runif(1, 0.2, 0.5),
    
    alpha_psi_dry_trail= runif(1, 0, 1),
    alpha_psi_dry_stream= runif(1, 0, 1),
    alpha_psi_wet_stream= runif(1, 0, 1),
    alpha_psi_wet_trail= runif(1, 0, 1)
  )}

params = c("alpha_N_dry_stream",
           "alpha_N_wet_stream",
           "alpha_N_dry_trail",
           "alpha_N_wet_trail",
           
           "alpha_psi_dry_stream",
           "alpha_psi_wet_stream",
           "alpha_psi_dry_trail",
           "alpha_psi_wet_trail",
           
           "sd",
           "sd.error",
           
           "alpha","beta",
           "zzzfit",
           "zzzfit.new"
)



data.new = list(
  
  y = yy,
  
  wet_stream = wet_stream,
  wet_trail = wet_trail,
  
  dry_stream = dry_stream,
  dry_trail = dry_trail,
  
  n = n,
  
  R = length(individual),
  
  D = samples,
  
  w = wNEW, 		
  
  individual = individual
  
)

# MCMC settings
ni <- 50000
nb <- 10000
nt <- 50
nc <- 3			

library("jagsUI")

out22 <- jags(data.new, inits, params, "modelDS.txt", 
              n.chains = nc, n.thin = nt, n.iter = ni, 
              n.burnin = nb, parallel = TRUE)

# Set working directory
setwd("/Users/Cici/Desktop/")
# Save model output
save(out22, file = "adjusted.rda")


