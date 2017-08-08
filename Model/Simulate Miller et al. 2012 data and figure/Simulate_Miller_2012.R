# This file simulates the data for Miller et al. 2012, runs the model, and saves the model output

expit <- function(x){ exp(x)/(1+exp(x)) }
	# Anti-logist function = to plogis()

### simulation parameters which can be adjusted
n <- 400 
### total number of sampled individuals
ni <- n/2

samples <- 4 
### number of surveys 

sigma.load <- 2.74
### sd of intensity for the population

sigma.error <- 1.17
#### sampling error

alpha <- 0.77
#### intercept
beta <- 0.32 
#### slope of intensity detection relationship

alpha_psi <- 1

alpha_N = 10


### build data set
individual <- rep(1:n, each = samples)

# -------- Prevalence/ Occupancy

psi <- plogis(alpha_psi)

z <- rbinom(n, size = 1, prob = psi)

# ---------- Infection intensity

N <- alpha_N

x <- rnorm(n, mean = N, sd = sigma.load)
	# this is a vector, and when z = 1, then you simulate the true infection intensity with mean mu.load and sd of sigma.load

x[ z ==0] <- NA


#---------- Detection 

p <- plogis(alpha + beta*x)
	# Model detection probability as a function of true infection intensity

p[z == 0]<- 0

y <- t(sapply(p, function(r){rbinom(n = samples, size = 1, prob = r)}))
	# Simulate the observation (detection/ non-detection process)

w <- y * t(sapply(x, function(r){rnorm(n = samples, mean = r, sd = sigma.error)}))
	# y in a vector of detection/ non-detection
	# x is a vector of true infection intensity
	# simulate the uncertainity in observation of infection intensity


which(z == 1 & x == 0)
ww <- w

w[w==0] <- NA
w[w=="NaN"] <- NA
w <- c(t(w))
Y <- rowSums(y)

#----------- Bundle the data

Infection <- apply(ww, 1, max, na.rm = T)
	Infection[Infection == "-Inf"] <- 0

zst <- Infection

zst[zst > 0] <- 1

inits = function() {
 list(
   # Prevalence
   z = zst, 
   
   alpha_psi = runif(1, -3, 3),

   # Infection intensity
   x = Infection, 
   
   alpha_N = 5 ,

   sigma.error = 0.2,
   
   # Detection
   alpha = 0, beta = 0
   )}

params = c("alpha_N",
   
         "alpha_psi",
       
          "sigma.delta",
          "sigma.error",
          
          "alpha","beta")

data = list(Y = Y,	# Total number of times the host was found (Y)
            
          n = n,		# total number of sites (R)

          samples = samples, 	# total number of surveys (T)

          w = w, 			# total number of individuals found at a site

          individual = individual	# vector with survey and ID numbers
)

#### Model

sink("model.txt")
cat("
model { 

### intensity model
for (i in 1:n){
  
  mu.delta[i] <- alpha_N

  mu[i] <- mu.delta[i] * z[i]

	x[i] ~ dnorm(mu[i], tau.delta)

}

for (i in 1:(n*2)){

    w[i] ~ dnorm(x[individual[i]], tau.error)

}

### detection model including relationship of intensity to detection

for (i in 1:n){
 
  	z[i] ~ dbern(psi[i])
	  
    logit(psi[i]) <- alpha_psi

    logit(p[i]) <- alpha + beta * x[i]
    
    P[i] <- p[i] * z[i]

	Y[i] ~ dbin(P[i], samples)

}

### priors

alpha_N ~ dnorm(0, 0.368)

alpha_psi ~ dnorm(0, 0.368)

tau.delta <- 1/(sigma.delta*sigma.delta)
sigma.delta ~ dunif(0, 10)

tau.error <- 1/(sigma.error*sigma.error)
sigma.error ~ dunif(0, 10)

alpha ~ dnorm(0.77, 13.22)
  # 2*SD =0.55
  # SD = 0.275
  # 1/(0.275)^2 = 13.22

beta ~ dnorm(0.32, 110.80)
  # 2*SD =0.19
  # SD = 0.095
  # 1/(0.095)^2 = 110.80


}
",fill=TRUE)
sink()


# MCMC settings
ni <- 10000
nt <- 10
nb <- 1000
nc <- 3			

library("jagsUI")

output <- jags(data, inits, params, "model.txt", 
              n.chains = nc, n.thin = nt, n.iter = ni, 
              n.burnin = nb, parallel = TRUE)


save(output, file = "Miller.rda")

##--------------


Model <- out$BUGSoutput$summary[c(1:10, 12, 13),c(1,2,3,7)]

TRUTH <- c(alpha,

alpha_N_dry,

alpha_N_stream,

alpha_N_trail,

alpha_N_wet,

alpha_psi_dry,

alpha_psi_stream,

alpha_psi_trail,

alpha_psi_wet,

beta,

sigma.load,

sigma.error

)

 
cbind(TRUTH, Model)



#---------------- Plot output
NAMES <- rownames(Model)

fores <- as.data.frame(cbind(NAMES, Model))

for (i in 2:ncol(fores)){
  fores[,i] <- as.numeric(as.character(fores[,i]))
}

colnames(fores) <- c("x", "y", "SD", "ylo", "yhi")



credplot.gg <- function(d){
  # d is a data frame with 4 columns
  # d$x gives variable names
  # d$y gives center point
  # d$ylo gives lower limits
  # d$yhi gives upper limits
  require(ggplot2)
  ggplot(d, aes(x=x, y=y, ymin=ylo, ymax=yhi))+geom_pointrange()+
    geom_hline(aes(x=0), lty=2) + geom_hline(aes(x=1), lty=2) + xlab('Parameter Name') +   
    coord_flip() + ylab('Probability') +
    scale_y_continuous(limits = c(-5,70)) + theme_bw() + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
}


credplot.gg(fores) + geom_point(aes(y = TRUTH, x = NAMES), size = 2, col = "red")


