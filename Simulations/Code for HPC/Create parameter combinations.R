n.swab <- c(1, 2, 3, 4)
n.PCR <- c(1, 2, 3, 4)
B1 <- c(-1.386294, 1.386294)
B2 <- 1
G1 <- c(-1.386294, 1.386294)
G2 <- 1
psi <- c(0.2, 0.8)
lambda <- 4

params <- expand.grid(n.swab, n.PCR, B1, B2, G1, G2, psi, lambda)

colnames(params) <- c("n.swab", "n.PCR", "B1", "B2", "G1", "G2", "psi", "lambda")

write.csv(params, file = "/Users/Cici/Desktop/Double Swab revisions/params.csv")
