# Density and CDF of S (for diffent values of r and \ell)
# They are found via the convolution: 
# S = Normal(0, 1-r^2) + Gamma(shape = (ell-1)/2, scale = r*sqrt(2/(ell-1)) )
#     - r*sqrt((ell-1)/2)
f.s <- function(s){
  s <- s + r*sqrt((ell-1)/2)
  k = (ell-1)/2
  theta = r*sqrt(2/(ell-1))
  b <- (2*(1-r^2))
  I.S <- function(x){x^(k-1)*exp(x*(2*s/b - 1/theta) - x^2/b)}
  return(exp(-s^2/b)/(sqrt(pi*b)*gamma(k)*theta^k) 
         * integrate(I.S, lower=0, upper=Inf, abs.tol=10^(-20))$value)
}

F.S <- function(s){
  s <-  s + r*sqrt((ell-1)/2)
  k = (ell-1)/2
  theta = r*sqrt(2/(ell-1))
  I.S <- function(x){x^(k-1)*exp(-x/theta)*pnorm((s-x)/sqrt(1-r^2))}
  1/(gamma(k)*theta^k) * integrate(I.S, lower=0, upper=Inf, abs.tol=10^(-20))$value
}  

# Plot df and CDF for many values of 'r' or many values of 'ell'
par(mfrow=c(1,2))
par(mar=c(3,6,1,1))
sd <- 5
x <- matrix(seq(-2*sd/3, sd, by = 0.005), ncol = 1)

# Creating figures: either fix 'r' and change 'ell', or the other way around
#r   <- .95
r    <- .9
ell <- 3
hx.r095    <- apply(x, 1, f.s)
hxCDF.r095 <- apply(x, 1, F.S)
#r   <- .8
ell <- 6
hx.r08    <- apply(x, 1, f.s)
hxCDF.r08 <- apply(x, 1, F.S)
#r   <- .6
ell <- 15
hx.r06    <- apply(x, 1, f.s)
hxCDF.r06 <- apply(x, 1, F.S)

#Density
plot(x, hx.r095, type="l", lty=1, col = "darkorange", lwd = 4, xlab="", 
     ylab="Density", cex.lab = 2.25, cex.axis = 2)
lines(x, hx.r08, type="l", lty=3, col = "blueviolet", lwd = 5)
lines(x, hx.r06, type="l", lty=6, col = "cornflowerblue", lwd = 4)
curve(dnorm(x, 0, 1), col="black", lty = 1, lwd=2, add=TRUE)

#CDF
plot(x, hxCDF.r095, type="l", lty=1, col = "darkorange", lwd = 4, xlab="", 
     ylab="CDF", cex.lab = 2.25, cex.axis = 2)
lines(x, hxCDF.r08, type="l", lty=3, col = "blueviolet", lwd = 5)
lines(x, hxCDF.r06, type="l", lty=6, col = "cornflowerblue", lwd = 4)
curve(pnorm(x, 0, 1), col="black", lty = 1, lwd=2, add=TRUE)
legend("bottomright", NULL, ncol = 1, cex = 2, 
       legend=c("l = 3", "l = 6", "l = 15", "N(0,1)"), 
       col=c("darkorange", "blueviolet", "cornflowerblue", "black"), 
       lty = c(1,3,6, 1), lwd = c(7,7,7,3))


