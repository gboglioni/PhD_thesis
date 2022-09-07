# Plot a shifted Log-norm with different values of sigma 
# (standardised to have mean 0, sd 1)
# This highlights what it means for 'r' to be close to 0 or close to 1
Erf    <- function(x){2 * pnorm(x * sqrt(2)) - 1}
# Find r parameter (as function of sigma)
r.lnorm      <- function(sig){Erf(sig/sqrt(2))/sqrt(exp(sig^2)-1)}
med.lnorm    <- function(sig){(exp(sig^2)-1)^(-1/2) * (exp(-sig^2/2)-1)}

# Function to graph density, as function of 'sig'
graph.logN <- function(sig = sqrt(0.1480239), legend = T){

# 'mu' parameter of the Log-normal (as function of 'sig') (such that Var[X]=1)
mu     <- (- sig^2 -log(exp(sig^2)-1))/2
# EX (also function of 'sig'), such that (X - EX) has mean 0 and unit variance
EX      <- sqrt(1/(exp(sig^2)-1)) 
# median
med          <- exp(mu) - EX
r            <- r.lnorm(sig)
# because we deal with a standard distribution, mu_V = r
mu_V  <- r
x      <- seq(-2, 2.5, by = 0.002)
df     <- dlnorm(x+EX, meanlog = mu, sdlog = sig)
cdf    <- plnorm(x+EX, meanlog = mu, sdlog = sig)

par(mar = c(2, 5, 2, 1))
plot(x, df, type = "l", col="brown3", lwd=3, 
     main = paste("r =", round(r,3)), ylim = c(-max(df)/15, max(df)),
     xlab="", ylab="Density", cex.lab = 2.25, cex.axis = 2)
# lines at the median and mean
abline(v = med, lwd = 2, col = 'darkgoldenrod1')
abline(v = 0, lty = 3, lwd = 2)
abline(h=0)

ps <- matrix(c(0, 0, mu_V, 0), ncol = 2, byrow = T)
points(ps, col=c("steelblue4"), pch=c("|"), cex=1.7)
text(ps, col=c("steelblue4"), labels= c("", "r"), adj = c(.25,-.6), cex = 2.5)
arrows(x0 = 0, y0 = 0, x1 = mu_V, y1 = 0, length = 0.2, code = 2, lwd = 3, 
       col = 'steelblue4')
ifelse(legend, legend("topright", NULL, ncol = 1, 
               cex = 1.75, legend=c("median", "mean"), 
               col = c("darkgoldenrod1", "black"), lty = c(1,3), lwd = 5), "")
}

par(mfrow=c(3,1))

graph.logN(sig = sqrt(0.1480239))
graph.logN(sig = sqrt(0.6733226), legend = F)
graph.logN(sig = sqrt(1.596642),  legend = F)


