
####################################################
##########     Table of contents  ##################
####################################################

# 1. Load libraries
# 2. Define parameter estimates
# 2. Model true population parameters (prevalence & infection intensity)
# 3. Model data collected from swabs
# 4. Model data multiple qPCR runs from each swab collected 
# 5. Write file with results

####################################################
##########  1. Load libraries  #########
####################################################

library(foreach)
library(doParallel)
c1 <- makeCluster(4)

registerDoParallel(c1)

####################################################
##########  2. Define parameter values  #########
####################################################
param <- read.csv(file = "/home/gracedirenzo/params.csv")[,-1]

args <- commandArgs(trailingOnly = TRUE)

job_ID <- as.numeric(args[1])

n.swab <- param$n.swab[job_ID]     # Number of swabs collected per individual
n.PCR  <- param$n.PCR[job_ID]     # Number of qPCR runs per swab

B1 <- param$B1[job_ID]      # Slope and intercept of swab model
B2 <- param$B2[job_ID]

G1 <- param$G1[job_ID]  # Slope and intercept of qPCR model
G2 <- param$G2[job_ID]

lambda <- param$lambda[job_ID]             # True average pathogen infection intensity
psi    <- param$psi[job_ID]           # True pathogen prevalence
sigma1  <- 1    # [1] Variance of population infection intensity
sigma2  <- 1    # [2] Variance in sampling error of swabbing
sigma3  <- 1    # [3] Variance in qPCR error

# Define parameters and survey conditions
n.ind  <- 500   # Number of individuals sampled

# Number of simulated data set per parameter combination
n.sims <- 50

#--- indicate which parameters were used for this run
true_params <- matrix(rep(unlist(param[job_ID,]), times = n.sims), 
                      nrow = n.sims, ncol = ncol(param), 
                      byrow = TRUE)

colnames(true_params) <- paste("true_", colnames(param), sep = "")

####################################################
########## 3.  Simulate the data  #########
####################################################

results <- foreach(a = 1:n.sims, .combine = rbind) %dopar%{  #simulations

library(coda, lib.loc = "/sw/csc/R-3.2.3/library/")
library(rjags, lib.loc = "/sw/csc/R-3.2.3/library/")
library(jagsUI, lib.loc = "/sw/csc/R-3.2.3/library/")
                     
# Each parameter combination

# True infection status
z <- rbinom(n.ind, 1, prob = psi)
# 1 = infected
# 0 = uninfected

# True infection intensity
N <- rlnorm(n.ind, meanlog = log(lambda), sdlog = sigma1)

# Assign 0 to individuals with no infection  
N <- N * z

####################################################
##### 3.  Model data collected from swabs  #########
####################################################

# Create empty matrix to hold data
m <- w <- array(NA, dim = c(n.ind, n.swab))

# Calculate pathogen detection probability on swabs according to infection intensity
p_swab <- plogis(B1 + B2 * log(N))

for(i in 1:n.ind){
  for(j in 1:n.swab){
    
    # Determine if the pathogen was detected on a swab or not
    w[i, j] <- rbinom(1, 1, prob = p_swab[i])
    
    # If the pathogen was detected, determine what the infection load was
    m[i, j] <- rlnorm(1, meanlog = log(N[i]), sdlog = sigma2)
    
  }
}

# Assign 0 to individuals that did not have pathogen detected
m <- m * w

####################################################
##### 4.  Model data from multiple qPCR runs  ######
####################################################

# Create empty array to hold data
x <- y <- array(NA, dim = c(n.ind, n.swab, n.PCR))

p_PCR <- plogis(G1 + G2 * log(m))

for(i in 1:n.ind){
  for(j in 1:n.swab){
    for(k in 1:n.PCR){
      # Determine if the pathogen was detected on a swab or not
      y[i, j, k] <- rbinom(1, 1, prob = p_PCR[i, j])
      
      # If the pathogen was detected, determine what the infection load was
      x[i, j, k] <- rlnorm(1, meanlog = log(m[i, j]), sdlog = sigma3)
    }
  }
}

# Assign 0 to individuals that did not have pathogen detected
x <- x * y

#---- Replace 0's with NA's
N[N == 0] <- NA
m[m == 0] <- NA 
x[x == 0] <- NA 

# Data
win.data <- list(
  y = y,
  x = (x+0.001), 
  n.ind = dim(y)[1],
  n.swab = dim(y)[2],
  n.PCR = dim(y)[3]
)

# inits
zinit <- z
ninit <- N
winit <- w
minit <- m

inits <- function() {list(
  B1 = B1,
  B2 = B2,
  
  G1 = G1,
  G2 = G2,
  
  sigma1 = sigma1,
  sigma2 = sigma2,
  sigma3 = sigma3,
  
  psi = psi,
  lambda = lambda,
  z = zinit,  
  w = winit,
  
  m = minit+0.001, 
  N = ninit+0.001
)}

# Parameters

params2 <- c("B1","B2",
            "G1", "G2", 
            "sigma1","sigma2","sigma3", 
            "psi", "lambda")

# MCMC settings
ni <- 20000
nb <- 2000
nt <- 20
nc <- 3

#- Run the model

out <- jags(win.data, inits, params2, 
            "model_sampling.txt", 
            n.chains = nc, 
            n.thin = nt, 
            n.iter = ni, 
            n.burnin = nb, 
            parallel = TRUE)

c(out$mean$B1, 
  out$mean$B2,
  out$mean$G1,
  out$mean$G2,
  out$mean$sigma1,
  out$mean$sigma2,
  out$mean$sigma3, 
  out$mean$psi,
  out$mean$lambda
  )

}

colnames(results) <- c("B1","B2",
                       "G1", "G2", 
                       "sigma1","sigma2","sigma3", 
                       "psi", "lambda")
  
dat <- cbind(true_params, results)

####################################################
##########  5. Write file with results #########
####################################################

write.csv(dat, file = paste("/home/gracedirenzo/", "dat_high_", job_ID, ".csv", sep =""))

