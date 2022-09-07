# Example of three N(0,1) PIBD with X+Y+Z NOT a N(0,3)

# Set seed
set.seed(1)

#Sample size
n <- 1000000

#Generate a sample of (X,Y,Z)
X <- rnorm(n)
Y <- rnorm(n)
W <- rnorm(n)
Z <- abs(W)*sign(X*Y)

#Sum of the PIBD Normals
S <- X + Y + Z

# Plot (and save) density and CDF (with, as comparison, those of a N(0, sqrt(3)))
pdf(file = "density_CDF_S_romano.pdf", width = 14, height = 7)
par(mfrow=c(1,2))
par(mar=c(3,5,1,1))

# Density 
my.d <- density(S)
plot(my.d, lwd= 4, col = 'brown3', cex.axis = 1.5, cex.lab = 1.75, main = "", 
     xlab = "",  bty="n")
curve(dnorm(x, 0, sqrt(3)), col="black", lty = 1, lwd=2, add=TRUE)

# CDF (we take a subset of "only" 100000 points)
plot(ecdf(S[1:100000]), col='brown3', lwd=4, xlab ="", 
     ylab = "CDF", cex.axis = 1.5, cex.lab = 1.75, main = "", col.01line = NULL)
curve(pnorm(x, 0, sqrt(3) ), lwd=2,add=TRUE)
legend("bottomright", NULL, ncol=1, cex = 1.75, legend=c("S", "N(0, 3)"), 
       col=c("brown3", "black"), lty = c(1,1), lwd = c(4,2))
dev.off()

