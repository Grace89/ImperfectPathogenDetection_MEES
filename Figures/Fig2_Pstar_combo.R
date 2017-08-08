# Citation:
# DiRenzo, Graziella; Grant, Evan; Longo, Ana; Che-Castaldo, Christian; Zamudio, Kelly; Lips, Karen. 2017.  Imperfect pathogen detection from non-invasive skin swabs biases disease inference. Methods in Ecology & Evolution.

# This file:
  # Creates Figure 2 in the main text


# Load libraries
library(coda)
library(plyr)
# Load data
load("adjusted.rda")

load("Miller.rda")

# combine chains into 1 matrix
allchains2 <- rbind(as.matrix(out22$samples[[1]]), 
                    as.matrix(out22$samples[[2]]), 
                    as.matrix(out22$samples[[3]]))

# Make it into an mcmc object
allchains2 <- as.mcmc(allchains2)

# Calculate column means and 95% CI
AA <- cbind(colMeans(allchains2), HPDinterval(allchains2))

# Total number of MCMC draws
mcmc.sample <- nrow(allchains2)
y <- nrow(allchains2)

# Number of samples to simulate (collect 1 swab, 2 swabs, etc.)
samp <- 5

# Create 3 matrices
  # 1 zoospore
  # 3 zoospores
  # 5 zoospores
# Dimensions = # of MCMC draws by 5 samples collected (this is the x axis)
Pstar1 <- Pstar3 <- Pstar5 <- Pstar10 <- array(NA, dim = c(y, samp))

# Create an x variable used later
x <- rep(1:samp, each = y)

# Estimate P* for each situation
for(i in 1:y){
  for(j in 1:samp){
    Pstar1[i, j] <- 1 - (1 - plogis(out22$sims.list$alpha[i] + out22$sims.list$beta[i] * log(1)))^j
    Pstar3[i, j] <- 1 - (1 - plogis(out22$sims.list$alpha[i] + out22$sims.list$beta[i] * log(3)))^j
    Pstar5[i, j] <- 1 - (1 - plogis(out22$sims.list$alpha[i] + out22$sims.list$beta[i] * log(5)))^j
    Pstar10[i, j] <- 1 - (1 - plogis(out22$sims.list$alpha[i] + out22$sims.list$beta[i] * log(10)))^j
  }
}

spore1 <- rep(1, times = nrow(Pstar1)*ncol(Pstar1))
samples <- rep(1:5, each = nrow(Pstar1))
values1 <- 
c(Pstar1[,1], 
  Pstar1[,2], 
  Pstar1[,3], 
  Pstar1[,4], 
  Pstar1[,5])

new_dat1 <- cbind(spore1, samples, values1)

spore3 <- rep(3, times = nrow(Pstar3)*ncol(Pstar3))
samples <- rep(1:5, each = nrow(Pstar3))
values3 <- 
  c(Pstar3[,1], 
    Pstar3[,2], 
    Pstar3[,3], 
    Pstar3[,4], 
    Pstar3[,5])

new_dat3 <- cbind(spore3, samples, values3)

spore5 <- rep(5, times = nrow(Pstar5)*ncol(Pstar5))
samples <- rep(1:5, each = nrow(Pstar5))
values5 <- 
  c(Pstar5[,1], 
    Pstar5[,2], 
    Pstar5[,3], 
    Pstar5[,4], 
    Pstar5[,5])
new_dat5 <- cbind(spore5, samples, values5)


spore10 <- rep(10, times = nrow(Pstar10)*ncol(Pstar10))
samples <- rep(1:5, each = nrow(Pstar10))
values10 <- 
  c(Pstar10[,1], 
    Pstar10[,2], 
    Pstar10[,3], 
    Pstar10[,4], 
    Pstar10[,5])
new_dat10 <- cbind(spore10, samples, values10)


ND <- rbind(new_dat1, new_dat3, new_dat5, new_dat10)
ND <- as.data.frame(ND)
str(ND)
colnames(ND) <- c("Zoospore", "Swabs", "Probability")

dat2 <- ddply(.data = ND, .variables = c("Zoospore", "Swabs"),
              .fun = summarize,
              est_mean = mean(Probability),
              CI.lower = HPDinterval(as.mcmc(Probability))[1],
              CI.upper = HPDinterval(as.mcmc(Probability))[2])



library(ggplot2)
swabs <- ggplot(data = dat2, aes(x = as.factor(Swabs), y = est_mean, ymin = CI.lower, ymax = CI.upper, col = as.factor(Zoospore)))+ 
  geom_pointrange(position=position_dodge(width=0.7), fatten = 1.5, size = 0.7)+
  scale_colour_manual(values = c("black", "gray41", "gray67", "gray87"), name = expression(paste(italic(Bd)," zoospores", sep = ""))) +
  ylab(expression(paste(italic(Bd), " detection probability", sep = ""))) + xlab("Number of swabs collected") + 
  theme_bw()+ ylim(c(0, 1))+
  theme(axis.text.x = element_text(size = 13, color = "black"), 
        axis.text.y = element_text(size = 13, color = "black"), 
        axis.title.y = element_text(size = 13, color = "black"), 
        axis.title.x =element_text(size = 13, color = "black"),
        legend.title =element_blank(),
        legend.text =element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "n") +
  geom_hline(yintercept = 0.95, lty = 2)


################### qPCR

y <- output$mcmc.info$n.samples
# Dimensions = # of MCMC draws by 5 samples collected (this is the x axis)
Pstar1 <- Pstar3 <- Pstar5 <- Pstar10 <- array(NA, dim = c(y, samp))

# Create an x variable used later
x <- rep(1:samp, each = y)

for(i in 1:y){
  for(j in 1:samp){
    Pstar1[i, j] <- 1 - (1 - plogis(output$sims.list$alpha[i] + output$sims.list$beta[i] * log(1)))^j
    Pstar3[i, j] <- 1 - (1 - plogis(output$sims.list$alpha[i] + output$sims.list$beta[i] * log(3)))^j
    Pstar5[i, j] <- 1 - (1 - plogis(output$sims.list$alpha[i] + output$sims.list$beta[i] * log(5)))^j
    Pstar10[i, j] <- 1 - (1 - plogis(output$sims.list$alpha[i] + output$sims.list$beta[i] * log(10)))^j
  }
}

spore1 <- rep(1, times = nrow(Pstar1)*ncol(Pstar1))
samples <- rep(1:5, each = nrow(Pstar1))
values1 <- 
  c(Pstar1[,1], 
    Pstar1[,2], 
    Pstar1[,3], 
    Pstar1[,4], 
    Pstar1[,5])

new_dat1 <- cbind(spore1, samples, values1)

spore3 <- rep(3, times = nrow(Pstar3)*ncol(Pstar3))
samples <- rep(1:5, each = nrow(Pstar3))
values3 <- 
  c(Pstar3[,1], 
    Pstar3[,2], 
    Pstar3[,3], 
    Pstar3[,4], 
    Pstar3[,5])

new_dat3 <- cbind(spore3, samples, values3)

spore5 <- rep(5, times = nrow(Pstar5)*ncol(Pstar5))
samples <- rep(1:5, each = nrow(Pstar5))
values5 <- 
  c(Pstar5[,1], 
    Pstar5[,2], 
    Pstar5[,3], 
    Pstar5[,4], 
    Pstar5[,5])

new_dat5 <- cbind(spore5, samples, values5)

spore10 <- rep(10, times = nrow(Pstar10)*ncol(Pstar10))
samples <- rep(1:5, each = nrow(Pstar10))
values10 <- 
  c(Pstar10[,1], 
    Pstar10[,2], 
    Pstar10[,3], 
    Pstar10[,4], 
    Pstar10[,5])
new_dat10 <- cbind(spore10, samples, values10)


ND <- rbind(new_dat1, new_dat3, new_dat5, new_dat10)
ND <- as.data.frame(ND)
str(ND)
colnames(ND) <- c("Zoospore", "Swabs", "Probability")

dat2 <- ddply(.data = ND, .variables = c("Zoospore", "Swabs"),
              .fun = summarize,
              est_mean = mean(Probability),
              CI.lower = HPDinterval(as.mcmc(Probability))[1],
              CI.upper = HPDinterval(as.mcmc(Probability))[2])


qPCR <- ggplot(data = dat2, aes(x = as.factor(Swabs), y = est_mean, ymin = CI.lower, ymax = CI.upper, col = as.factor(Zoospore)))+ 
  geom_pointrange(position=position_dodge(width=0.7), fatten = 1.5, size = 0.7)+
  scale_colour_manual(values = c("black", "gray41", "gray67", "gray87"), name = expression(paste(italic(Bd)," zoospores", sep = ""))) +
  ylab(expression(paste(italic(Bd), " detection probability", sep = ""))) + xlab("Number of qPCR runs") + 
  theme_bw()+ ylim(c(0, 1))+
  theme(axis.text.x = element_text(size = 13, color = "black"), 
        axis.text.y = element_blank(), 
        axis.title.y = element_blank(), 
        axis.title.x =element_text(size = 13, color = "black"),
        legend.title =element_text(size = 13, color = "black"),
        legend.text =element_text(size = 13, color = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  geom_hline(yintercept = 0.95, lty = 2)

#----

require(cowplot)

plot_grid(swabs, qPCR, rel_widths = c(1, 1.3), labels = "AUTO", vjust = 2.10, hjust = -1)


#setwd("/Users/Cici/Desktop/June 2017 Double swab/")
#ggsave("Pstar_Fig_16_June.pdf", width = 8, height = 5)







