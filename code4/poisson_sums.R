# Theoretical PMF of the sum S = X_1 + ... X_n
pois.S <- function(s,m){
  p.s <- 0
  for (k in 0:m){
    p.k <- (2*k^2+m^2-2*m*k-m)/2
    j   <- 0:p.k
    partial.sum <- choose(m,k)*ifelse(p.k<=s,1,0)*sum((-1)^j*(p.k-j)^s*choose(p.k,j))
    p.s = p.s + partial.sum
  }
  return(log(2)^s/factorial(s)*(1/2)^m*p.s)
}

# Function to  plot the PMF of S, compared to that of a Pois(n*log(2))
plot.pois <-function(m=3, max=m^2, legend=F, col.2='brown3', col.1='darkgoldenrod1'){
  domain <- 0:max
  probs  <- sapply(domain, pois.S, m = m)
  probs.pois <- sapply(domain, dpois, lambda = (m^2-m)/2*log(2))
  plot(domain,probs,type="h", lwd = 4, col= col.1,xlab="s",ylab="p(s)", 
       cex.lab=1.5, cex.axis=1.5)
  points(domain,probs,col=col.1, cex = 1.5, pch = 19)
  points(domain,probs.pois, col=col.2,cex = 2.5, pch=18)
  abline(h=0,col='black')
  if(legend){legend("topright", NULL, ncol = 1, cex = 2, legend=c("S", "Poisson"), 
         col=c(col.1, col.2), lwd = c(4,0), lty = c(1,0), pch = c(19,18), 
         pt.cex = c(3,3.5))}
}

# Create (and save) three plots (for different sample sizes)
pdf(file = "poisson_pmf.pdf", width = 8, height = 11)
par(mfrow=c(3,1))
par(mar=c(3,6,1,1))
plot.pois(m=3, 7)
plot.pois(m=5, 18)
plot.pois(m=7, 34, legend=T)
dev.off()

# Density and CDF of the standardised sample mean (if sample size goes to infinity)
f.s <- function(s){
  s <- s + r*sqrt((ell-1)/2)
  k = (ell-1)/2
  theta = r*sqrt(2/(ell-1))
  b <- (2*(1-r^2))
  I.S <- function(x){x^(k-1)*exp(x*(2*s/b - 1/theta) - x^2/b)}
  return(exp(-s^2/b)/(sqrt(pi*b)*gamma(k)*theta^k)
         *integrate(I.S, lower = 0, upper = Inf, abs.tol = 10^(-20))$value)
}

F.S <- function(s){
  s <-  s + r*sqrt((ell-1)/2)
  k = (ell-1)/2
  theta = r*sqrt(2/(ell-1))
  I.S <- function(x){x^(k-1)*exp(-x/theta)*pnorm((s-x)/sqrt(1-r^2))}
  1/(gamma(k)*theta^k)*integrate(I.S, lower=0, upper = Inf, abs.tol=10^(-20))$value
}  

sd <- 5
x <- matrix(seq(-2*sd/3, sd, by = 0.005), ncol = 1)
r    <- sqrt(log(2))
ell <- 2
hx.r08    <- apply(x, 1, f.s)
hxCDF.r08 <- apply(x, 1, F.S)

# Plot density and CDF, also compared to a N(0,1)
pdf(file = "poisson_CLT.pdf", width = 16, height = 8)
par(mfrow=c(1,2))
par(mar=c(3,6,1,1))
#Density
plot(x, hx.r08, type="l", lty=1, col = "brown3", lwd = 4, xlab="", 
     ylab="Density", cex.lab = 1.5, cex.axis = 1.25)
curve(dnorm(x, 0, 1), col="black", lty = 1, lwd=2, add=TRUE)
#CDF
plot(x, hxCDF.r08, type="l", lty=1, col = "brown3", lwd = 4, xlab="", 
     ylab="CDF", cex.lab = 1.5, cex.axis = 1.25)
curve(pnorm(x, 0, 1), col="black", lty = 1, lwd=2, add=TRUE)
legend("bottomright", NULL, ncol = 1, cex = 1.5, legend=c("Z", "N(0,1)"), 
       col=c("brown3", "black"), lty = c(1, 1), lwd = c(7,3))
dev.off()
