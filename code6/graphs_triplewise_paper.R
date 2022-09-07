# density function for the VG(n, 0, 1, 0) distribution
density_VG <- function (x, n) {
  return (1 / (sqrt(pi) *  gamma(n / 2)) * (abs(x) / 2) ^ ((n - 1) / 2) * 
            besselK(abs(x), (n - 1) / 2))
}

# density of the convolution in the main theorem
f.s <- function (s) {
  convolution_density <- function(y) {dnorm(y / sqrt(1 - r ^ 2)) / sqrt(1 - r ^ 2)*
      density_VG((s - y) / (r / sqrt(ell - 1)), ell - 1) / (r / sqrt(ell - 1))}
  return (integrate(convolution_density, lower = -Inf, upper = s, 
                    abs.tol = 10 ^ (-20))$value + 
            integrate(convolution_density, lower = s, upper = Inf, 
                      abs.tol = 10 ^ (-20))$value)
}

# cdf of the convolution in the main theorem
F.S <- function (s) {
  convolution_cdf <- function(y) {pnorm(y / sqrt(1 - r ^ 2)) * 
      density_VG((s - y) / (r / sqrt(ell - 1)), ell - 1) / (r / sqrt(ell - 1))}
  return (integrate(convolution_cdf, lower = -Inf, upper = s, 
                    abs.tol = 10 ^ (-20))$value + 
            integrate(convolution_cdf, lower = s, upper = Inf, 
                      abs.tol = 10 ^ (-20))$value)
}  


# Plot df and CDF for many values of 'r' or many values of 'ell'
par(mfrow = c(1, 2))
par(mar = c(3, 6, 1, 1))
x <- matrix(seq(-3, 3, by = 0.005), ncol = 1)

# To create Figures: either fix r and change 'ell', or the other way around
r <- 0.99
ell <- 2
hx.r099    <- apply(x, 1, f.s)
hxCDF.r099 <- apply(x, 1, F.S)

#r <- 0.8
ell <- 4
hx.r08    <- apply(x, 1, f.s)
hxCDF.r08 <- apply(x, 1, F.S)

#r <- 0.6
ell <- 6
hx.r06    <- apply(x, 1, f.s)
hxCDF.r06 <- apply(x, 1, F.S)

#Density
plot(x, hx.r099, type = "l", lty = 1, col = "darkorange", lwd = 4, 
     xlab="", ylab="Density", cex.lab = 2.25, cex.axis = 2)
lines(x, hx.r08, type = "l", lty = 3, col = "blueviolet", lwd = 5)
lines(x, hx.r06, type = "l", lty = 6, col = "cornflowerblue", lwd = 4)
curve(dnorm(x, 0, 1), col="black", lty = 1, lwd=2, add=TRUE)

#CDF
plot(x, hxCDF.r099, type = "l", lty = 1, col = "darkorange", lwd = 4, 
     xlab="", ylab="CDF", cex.lab = 2.25, cex.axis = 2)
lines(x, hxCDF.r08, type = "l", lty = 3, col = "blueviolet", lwd = 5)
lines(x, hxCDF.r06, type = "l", lty = 6, col = "cornflowerblue", lwd = 4)
curve(pnorm(x, 0, 1), col="black", lty = 1, lwd=2, add=TRUE)
legend("bottomright", NULL, ncol = 1, cex = 1.1, 
       legend=c("r = 0.99", "r = 0.8", "r = 0.6", "N(0,1)"), 
       col=c("darkorange", "blueviolet", "cornflowerblue", "black"), 
       lty = c(1,3,6,1), lwd = c(4,5,4,2))

