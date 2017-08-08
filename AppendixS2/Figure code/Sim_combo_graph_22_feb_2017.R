dat <- numeric()

for(i in c(1:128)){
  dat <- rbind(dat, read.csv(paste("/Users/Cici/Desktop/Double Swab revisions/SuperComputer simulation/17 Feb run/dat_", i, ".csv", sep = "")))
}

low <- seq(from = 1, to = 128, by = 1)
missing <- c(1, 2, 3, 4, 17, 33, 34, 49, 50, 65)
  
low <- low[-missing]  

for(i in low){
  dat <- rbind(dat, read.csv(paste("/Users/Cici/Desktop/Double Swab revisions/SuperComputer simulation/22 Mar 2017 scripts/dat_low_", i, ".csv", sep = "")))
}

library(ggplot2)

cols <- c("1" = "grey0",
          "2" = "grey27",
          "3" = "grey51", 
          "4" = "grey71")

dat2 <- dat

dat2$true_B1_new <- as.factor(dat2$true_B1)
levels(dat2$true_B1_new) <- c(
  expression(paste("logit(", beta[0], ") = 0.2", sep = "")), 
  expression(paste("logit(", beta[0], ") = 0.8", sep = "")))

dat2$true_G1_new <- as.factor(dat2$true_G1)
levels(dat2$true_G1_new) <- c(
  expression(paste("logit(", gamma[0], ") = 0.2", sep = "")), 
  expression(paste("logit(", gamma[0], ") = 0.8", sep = "")))
  
dat2$true_psi_new <- as.factor(dat2$true_psi)
levels(dat2$true_psi_new) <- c(
  expression(paste(psi, " = 0.2", sep = "")), 
  expression(paste(psi, " = 0.8", sep = "")))

dat2$true_lambda_new <- as.factor(dat2$true_lambda)
levels(dat2$true_lambda_new) <- c(
  expression(paste(lambda, " = 2", sep = "")), 
  expression(paste(lambda, " = 4", sep = "")))


#--- Average occupancy probability
ggplot(data = dat2, aes(x = as.factor(true_n.swab), 
                        y = psi, 
                        fill = as.factor(true_n.PCR))) + 
  scale_fill_manual("Number of \ndiagnostic runs", values=cols)+
  scale_x_discrete(breaks = c(1, 2, 3, 4),
                   labels = c(1, 2, 3, 4))+
  geom_boxplot(notch = TRUE) + 
  facet_grid(as.factor(true_psi_new) * as.factor(true_lambda_new)~ as.factor(true_B1_new) * as.factor(true_G1_new), labeller = label_parsed)+
  geom_hline(data= dat2, aes(yintercept = true_psi), lty = 2) +
  xlab("Number of samples collected")+
  ylab(expression(paste("Estimated occupancy probability (", psi, ")", sep = "")))+
  theme_bw()+ 
  theme(axis.text.x = element_text(size = 17, color = "black"), 
        axis.text.y = element_text(size = 17, color = "black"), 
        axis.title.y = element_text(size = 17, color = "black"), 
        axis.title.x =element_text(size = 17, color = "black"),
        legend.title =element_text(size = 17, color = "black"),
        legend.text =element_text(size = 17, color = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        strip.text.x = element_text(size=14),
        strip.text.y = element_text(size=14),
        strip.background = element_rect(colour="black", fill="white", 
                                        size=1.5, linetype="solid")) 
  

#----- Average infection intensity
ggplot(data = dat2, aes(x = as.factor(true_n.swab), y = lambda, fill = as.factor(true_n.PCR))) + 
  scale_fill_manual("Number of \ndiagnostic runs", values=cols)+
  scale_x_discrete(breaks = c(1, 2, 3, 4),
                   labels = c(1, 2, 3, 4))+
  geom_boxplot(notch = TRUE) + 
  facet_grid(as.factor(true_G1_new) ~ as.factor(true_B1_new), labeller = label_parsed)+
  geom_hline(data= dat2, aes(yintercept = true_lambda), lty = 2) +
  xlab("Number of samples collected")+
  ylab(expression(paste("Estimated average infection intensity (", lambda, ")", sep = "")))+
  theme_bw()+ 
  theme(axis.text.x = element_text(size = 17, color = "black"), 
        axis.text.y = element_text(size = 17, color = "black"), 
        axis.title.y = element_text(size = 17, color = "black"), 
        axis.title.x =element_text(size = 17, color = "black"),
        legend.title =element_text(size = 17, color = "black"),
        legend.text =element_text(size = 17, color = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        strip.text.x = element_text(size=14),
        strip.text.y = element_text(size=14),
        strip.background = element_rect(colour="black", fill="white", 
                                        size=1.5, linetype="solid")) 


#----- Detection of samples
ggplot(data = dat2, aes(x = as.factor(true_n.swab), y = plogis(B1), fill = as.factor(true_n.PCR))) + 
  scale_fill_manual("Number of \ndiagnostic runs", values=cols)+
  scale_x_discrete(breaks = c(1, 2, 3, 4),
                   labels = c(1, 2, 3, 4))+
  geom_boxplot(notch = TRUE) + 
  facet_grid(as.factor(true_G1_new) ~ as.factor(true_B1_new), labeller = label_parsed)+
  geom_hline(data= dat2, aes(yintercept = plogis(true_B1)), lty = 2) +
  xlab("Number of samples collected")+
  ylab(expression(paste("Estimated ", beta[0], sep = "")))+
  theme_bw()+ 
  theme(axis.text.x = element_text(size = 17, color = "black"), 
        axis.text.y = element_text(size = 17, color = "black"), 
        axis.title.y = element_text(size = 17, color = "black"), 
        axis.title.x =element_text(size = 17, color = "black"),
        legend.title =element_text(size = 17, color = "black"),
        legend.text =element_text(size = 17, color = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        strip.text.x = element_text(size=14),
        strip.text.y = element_text(size=14),
        strip.background = element_rect(colour="black", fill="white", 
                                        size=1.5, linetype="solid")) 

#----- Detection of diagnostic test
ggplot(data = dat2, aes(x = as.factor(true_n.swab), y = plogis(G1), fill = as.factor(true_n.PCR))) + 
  scale_fill_manual("Number of \ndiagnostic runs", values=cols)+
  scale_x_discrete(breaks = c(1, 2, 3, 4),
                   labels = c(1, 2, 3, 4))+
  geom_boxplot(notch = TRUE) + 
  facet_grid(as.factor(true_G1_new) ~ as.factor(true_B1_new), labeller = label_parsed)+
  geom_hline(data= dat2, aes(yintercept = plogis(true_G1)), lty = 2) +
  xlab("Number of samples collected")+
  ylab(expression(paste("Estimated ", gamma[0], sep = "")))+
  theme_bw()+ 
  theme(axis.text.x = element_text(size = 17, color = "black"), 
        axis.text.y = element_text(size = 17, color = "black"), 
        axis.title.y = element_text(size = 17, color = "black"), 
        axis.title.x =element_text(size = 17, color = "black"),
        legend.title =element_text(size = 17, color = "black"),
        legend.text =element_text(size = 17, color = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        strip.text.x = element_text(size=14),
        strip.text.y = element_text(size=14),
        strip.background = element_rect(colour="black", fill="white", 
                                        size=1.5, linetype="solid")) 