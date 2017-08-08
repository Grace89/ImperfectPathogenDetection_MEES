# Citation:
# DiRenzo, Graziella; Grant, Evan; Longo, Ana; Che-Castaldo, Christian; Zamudio, Kelly; Lips, Karen. 2017.  Imperfect pathogen detection from non-invasive skin swabs biases disease inference. Methods in Ecology & Evolution.

# This file:
  # Creates Figure 1 in the main text

load("adjusted.rda")
load("Miller.rda")

library(coda)

allchains2 <- rbind(as.matrix(out22$samples[[1]]), 
                    as.matrix(out22$samples[[2]]), 
                    as.matrix(out22$samples[[3]]))

allchains2 <- as.mcmc(allchains2)
#-----------

AA <- cbind(colMeans(allchains2), HPDinterval(allchains2))

#--------

mcmc.sample <- nrow(allchains2)

# We scaled the log(II +0.01)

parasite_density <- seq(from = -5, to = 5, by = 0.1)

# Stream
mean <- plogis(AA["alpha",1] + AA["beta",1]* parasite_density)

upper <- plogis(AA["alpha",2] + AA["beta",2]* parasite_density)
lower <- plogis(AA["alpha",3] + AA["beta",3]* parasite_density)

array.pred <- array(NA, dim = c(length(parasite_density), mcmc.sample))

for(i in 1:mcmc.sample){
  
  array.pred[,i] <- plogis(allchains2[i, grep("alpha", colnames(allchains2))[9]] +allchains2[i, grep("beta", colnames(allchains2))] * parasite_density)

}


sub.set <- sort(sample(1:mcmc.sample, size = 500))


par(mar = c(5,5,5,5))
plot(mean ~ parasite_density, 
     xlab = expression(paste(italic(Bd), " infection intensity", sep = "")), 
     las = 1, ylab = expression(paste(italic(Bd), " detection probability", sep = "")),
     lwd = 4, type = "l", col = "black", main = "",
     cex.axis = 1.5, cex.lab = 1.5, ylim = c(0, 1),
     xlim = c(min(parasite_density), max(parasite_density)),
     xaxt = "n"
     )

axis(1, at = c(log(0.01), log(0.1), log(1), log(10), log(100), log(1000)), lab = c(0.01, 0.1, 1, 10, 100, 1000), cex.axis = 1.5)

for(i in sub.set){
  
  lines(parasite_density, array.pred[,i], type = "l", lwd = 1, col = "gray60")
  
}

lines(parasite_density, mean, type = "l", lwd = 4, col = "black")

cbind(parasite_density, mean)


#######------------- Miller

allchains21 <- rbind(as.matrix(output$samples[[1]]), 
                    as.matrix(output$samples[[2]]), 
                    as.matrix(output$samples[[3]]))

allchains21 <- as.mcmc(allchains21)
#-----------

AA1 <- cbind(colMeans(allchains21), HPDinterval(allchains21))

#--------

mcmc.sample <- nrow(allchains21)

# We scaled the log(II +0.01)

parasite_density <- seq(from = -5, to = 5, by = 0.1)

# Stream
mean1 <- plogis(AA1["alpha",1] + AA1["beta",1]* parasite_density)

array.pred1 <- array(NA, dim = c(length(parasite_density), mcmc.sample))

for(i in 1:mcmc.sample){
  
  array.pred1[,i] <- plogis(allchains21[i, grep("alpha", colnames(allchains21))[3]] +allchains21[i, grep("beta", colnames(allchains21))] * parasite_density)
  
}


sub.set <- sort(sample(1:mcmc.sample, size = 500))


par(mar = c(5,5,5,5))
plot(mean1 ~ parasite_density, xlab = "Ln(Infection intensity + 0.01)", 
     las = 1, ylab = "Pathogen detection probability",
     lwd = 4, type = "l", col = "black", main = "",
     cex.axis = 1.5, cex.lab = 1.5, ylim = c(0, 1)
)


for(i in sub.set){
  
  lines(parasite_density, array.pred1[,i], type = "l", lwd = 1, col = "gray60")
  
}

lines(parasite_density, mean1, type = "l", lwd = 4, col = "forestgreen")


legend(c(4.5, 0.05),
       c("qPCR", "Swabbing"), col = c("forestgreen", "black"), lwd = c(4, 4), bty = "n", cex = 1.5)




############################# Combo figures
dev.off()
postscript(file = "RPlot.eps", onefile = TRUE,
           width = 9, height = 7, horizontal = FALSE)
par(mar = c(5,5,5,5))

plot(mean1 ~ parasite_density, 
     xlab = expression(paste(italic(Bd), " infection intensity", sep = "")), 
     las = 1, ylab = expression(paste(italic(Bd), " detection probability", sep = "")),
     lwd = 4, type = "l", col = "black", main = "",
     cex.axis = 1.5, cex.lab = 1.5, ylim = c(0, 1),
     xlim = c(min(parasite_density), max(parasite_density)),
     xaxt = "n"
)


for(i in sub.set){
  lines(parasite_density, array.pred1[,i], type = "l", lwd = 1, col = "gray77")
}

lines(parasite_density, mean1, type = "l", lwd = 4, col = "black")


axis(1, at = c(log(0.01), log(0.1), log(1), log(10), log(100), log(1000)), lab = c(0.01, 0.1, 1, 10, 100, 1000), cex.axis = 1.5)

for(i in sub.set){
  lines(parasite_density, array.pred[,i], type = "l", lwd = 1, col = "gray77")
}

lines(parasite_density, mean, type = "l", lwd = 4, col = "gray37")

legend(2.5, 0.1,
       c("qPCR", "Swabbing"), col = c("black", "gray37"), lwd = c(4, 4), bty = "n", cex = 1.1)

dev.off()
