# Citation:
# DiRenzo, Graziella; Grant, Evan; Longo, Ana; Che-Castaldo, Christian; Zamudio, Kelly; Lips, Karen. 2017.  Imperfect pathogen detection from non-invasive skin swabs biases disease inference. Methods in Ecology & Evolution.

# This file:
  # Creates a figure to visulize the Bayesian p-value 
  # Calculates the Bayesian p-value 


# Load the data
load("adjusted.rda")
load("NOTadjusted.rda")

#----- Visulaize Model fit

library(plyr)

allchains2 <- rbind(as.matrix(out22$samples[[1]]),
                    as.matrix(out22$samples[[2]]), 
                    as.matrix(out22$samples[[3]]))

z.act <- grep("zzzfit", colnames(allchains2))[1]
z.new <- grep("zzzfit", colnames(allchains2))[2]

p2 <- round(mean(allchains2[,z.act]>allchains2[,z.new]),2)
m3 <- round_any(min(allchains2[,z.new],allchains2[,z.act]), 10, f = floor)
m4 <- round_any(max(allchains2[,z.new],allchains2[,z.act]), 10, f = ceiling)
plot(allchains2[,z.act], allchains2[,z.new], xlab = expression(T^{obs}), ylab=expression(T^{rep}), cex.lab=1, cex.axis=1, xlim = c(m3,m4), ylim = c(m3,m4), las=1)
abline(0, 1, lwd=2)
mtext(paste("Bayesian p-value =", p2, sep = " "), side = 3, line = -2, at=(m4-m3)*.1 + m3)



library(plyr)

allchains2 <- rbind(as.matrix(out24$samples[[1]]),
                    as.matrix(out24$samples[[2]]), 
                    as.matrix(out24$samples[[3]]))

z.act <- grep("zzzfit", colnames(allchains2))[1]
z.new <- grep("zzzfit", colnames(allchains2))[2]

p2 <- round(mean(allchains2[,z.act]>allchains2[,z.new]),2)
m3 <- round_any(min(allchains2[,z.new],allchains2[,z.act]), 10, f = floor)
m4 <- round_any(max(allchains2[,z.new],allchains2[,z.act]), 10, f = ceiling)
plot(allchains2[,z.act], allchains2[,z.new], xlab = expression(T^{obs}), ylab=expression(T^{rep}), cex.lab=1, cex.axis=1, xlim = c(m3,m4), ylim = c(m3,m4), las=1)
abline(0, 1, lwd=2)
mtext(paste("Bayesian p-value =", p2, sep = " "), side = 3, line = -2, at=(m4-m3)*.1 + m3)
