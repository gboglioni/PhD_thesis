# 0-correlation, but VaR(X + Y) different to VaR(X.ind + Y.ind)
set.seed(1)
n       <- 5000
X       <- rexp(n)
Y       <- rnorm(n, 1, X)
X.indep <- rexp(n)

# Save scatterplot of X versus Y
pdf(file = "0_correl_scatter.pdf", width = 7, height = 7)
par(mar=c(5,5,1,1))
plot(rank(X)/(n+1), rank(Y)/(n+1), cex = 0.5, cex.axis = 1.25, cex.lab = 1.5)
dev.off()

# Theoretical values, 
# (numbers obtained from calculations in Mathematica)
VaR.X         <- qexp(0.995)
VaR.Y         <- 6.607 
VaR.sum.indep <- 8.351
VaR.sum.dep   <- 10.405
sum.VaR       <- VaR.X + VaR.Y

# Diversification under independence
1 - VaR.sum.indep/sum.VaR
# Diversification under dependence
1 - VaR.sum.dep/sum.VaR
