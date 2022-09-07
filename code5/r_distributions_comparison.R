# Analysing the "r" parameter as measure of tail-heaviness for common distributions
# The analysis is done for r* (i.e. 1-r) so that larger r* imply heavier-tail

# Log Normal(mu, s) (s is like sigma^2)
# min r: 1-sqrt(2/pi) = 0.2022
# max r: 1
erf        <- function(x){2 * pnorm(x * sqrt(2)) - 1}
r.logN     <- function(s){1-erf(sqrt(s/2))/sqrt(exp(s)-1)}
r.logN.inv <- function(r){uniroot(function(x){r.logN(x) - r}, lower = 0.00001, 
                                  upper=100)$root}
k.logN     <- function(s){exp(4*s) + 2*exp(3*s) + 3*exp(2*s) - 6  +3}
k.logN.inv <- function(k){uniroot(function(x){k.logN(x) - k}, lower = 0.0001, 
                                  upper=100)$root}
RQW.logN   <- function(s,q=7/8){(qlnorm(0.5+q/2, sdlog=sqrt(s))
                                 +qlnorm(1-q/2, sdlog=sqrt(s))
                                 -2*qlnorm(3/4, sdlog=sqrt(s)))/(qlnorm(0.5+q/2, 
                                  sdlog=sqrt(s))-qlnorm(1-q/2, sdlog=sqrt(s)))}
RQW.logN.inv <- function(k){uniroot(function(x){RQW.logN(x) - k}, 
                                    lower = 0.0001, upper=1000)$root}

# Gamma with shape parameter a
# min r: 1-sqrt(2/pi) = 0.2030 (seems to have numerical problems to evaluate at lower r)
# max r: 1
r.Gam <- function(a){
  mx <- qgamma(1/2, shape = a, rate = 1)
  1-2*mx^a*exp(-mx)/(gamma(a)*sqrt(a))
}
r.Gam.inv <- function(r){uniroot(function(x){r.Gam(x) - r}, 
                                 lower = .0001, upper=130)$root}
k.Gam     <- function(a){6/a + 3}
k.Gam.inv <- function(k){uniroot(function(x){k.Gam(x) - k}, 
                                 lower = .0001, upper=130)$root}
RQW.Gam   <- function(a,q=7/8){(qgamma(0.5+q/2, shape=a)+qgamma(1-q/2, shape=a)
                                -2*qgamma(3/4, shape=a))/(qgamma(0.5+q/2, shape=a)
                                                          -qgamma(1-q/2, shape=a))}
RQW.Gam.inv <- function(r){uniroot(function(x){RQW.Gam(x) - r}, 
                                   lower = .001, upper=500)$root}

# Weibull with shape parameter k
# min r: 0.191 (roughly, numerically)
# max r: 1
r.Wei     <- function(k){1-(2*igamma(1+1/k,log(2)) - gamma(1+1/k))/sqrt(gamma(1+2/k) 
                                                                        - gamma(1+1/k)^2)}
r.Wei.inv <- function(r){uniroot(function(x){r.Wei(x) - r}, lower = 0.05, upper=2.9)$root}
G1        <- function(k){gamma(1+1/k)}
G2        <- function(k){gamma(1+2/k)}
G3        <- function(k){gamma(1+3/k)}
G4        <- function(k){gamma(1+4/k)}
k.Wei     <- function(k){(-6*G1(k)^4+12*G1(k)^2*G2(k)-3*G2(k)^2-4*G1(k)*G3(k)
                          +G4(k))/(G2(k)-G1(k)^2)^2 + 3}
k.Wei.inv <- function(k){uniroot(function(x){k.Wei(x) - k}, lower = 0.05, upper=2.9)$root}
RQW.Wei   <- function(k,q=7/8){(qweibull(0.5+q/2, shape=k)+qweibull(1-q/2, shape=k)
                                -2*qweibull(3/4, shape=k))/(qweibull(0.5+q/2, shape=k)
                                                            -qweibull(1-q/2, shape=k))}
RQW.Wei.inv <- function(r){uniroot(function(x){RQW.Wei(x) - r}, lower = 0.05,
                                   upper=2.9)$root}

# Frechet
# min r: 0.2453 (roughly, numerically)
# max r: 1
# (kurtosis does not exist for r > 0.4013)
igamma <- function(a, x){gamma(a) * (1 - pgamma(x,a,1))}
r.Fre  <- function(a){1-(gamma(1-1/a) -
                      2*igamma(1-1/a, log(2)))/sqrt(gamma(1-2/a)-gamma(1-1/a)^2)}
r.Fre.inv <- function(r){uniroot(function(x){r.Fre(x) - r}, lower = 2.0001, 
                                 upper=10000)$root}
k.Fre  <- function(a){(gamma(1-4/a)-4*gamma(1-3/a)*gamma(1-1/a)
                       +3*gamma(1-2/a)^2)/(gamma(1-2/a)-gamma(1-1/a)^2)^2-6  +3}
k.Fre.inv <- function(k){uniroot(function(x){k.Fre(x) - k}, lower = 4.0001, 
                                 upper=1000)$root}
qfrechet <- function(q,a){log(1/q)^{-1/a}}
RQW.Fre   <- function(a,q=7/8){(qfrechet(0.5+q/2, a)+qfrechet(1-q/2, a)
                                -2*qfrechet(3/4, a))/(qfrechet(0.5+q/2, a)
                                                      -qfrechet(1-q/2, a))}
RQW.Fre.inv <- function(rqw){uniroot(function(x){RQW.Fre(x) - rqw}, lower = 0.2,
                                     upper=1000)$root}

# Find a in Pareto(a, lambda) (a>2)
# min r: 1-log(2) = 0.3069
# max r: 1
# (Kurtosis does not exist past r > 0.4648)
r.Par     <- function(a){1-sqrt(a*(a-2))*(2^(1/a)-1) }
r.Par.inv <- function(r){uniroot(function(x){r.Par(x) - r}, lower = 2.00001, 
                                 upper=10000)$root}
k.Par     <- function(a){6*(a^3+a^2-6*a-2)/(a*(a-3)*(a-4))  +3}
k.Par.inv <- function(k){uniroot(function(x){k.Par(x) - k}, lower = 4.00001,
                                 upper=10000)$root}
qpareto   <- function(q,a){(1-q)^(-1/a)-1}
RQW.Par   <- function(a,q=7/8){(qpareto (0.5+q/2, a)+qpareto (1-q/2, a)
                                -2*qpareto (3/4, a))/(qpareto (0.5+q/2, a)
                                                      -qpareto(1-q/2, a))}
RQW.Par.inv <- function(r){uniroot(function(x){RQW.Par(x) - r}, lower = .1, 
                                   upper=1000)$root}

# Student (v) (v > 2)
# min r: 1-sqrt(2/pi) = 0.2028 (seems to have numerical problems to evaluate at lower r)
# max r: 1
# (Kurtosis does not exist past r > 0.2928)
r.t     <- function(v){1-2*sqrt((v-2)/pi)*gamma((v+1)/2)/((v-1)*gamma(v/2))}
r.t.inv <- function(r){uniroot(function(x){r.t(x) - r}, lower = 2.00001, 
                               upper=340)$root}
# Excess Kurtosis
k.t     <- function(v){6/(v-4)  +3}
k.t.inv <- function(k){uniroot(function(x){k.t(x) - k}, lower = 4.00001, 
                               upper=340)$root}
RQW.t   <- function(v,q=7/8){(qt(0.5+q/2, df=v)+qt(1-q/2, df=v)
                              -2*qt(3/4, df=v))/(qt(0.5+q/2, df=v)-qt(1-q/2, df=v))}
RQW.t.inv <- function(r){uniroot(function(x){RQW.t(x) - r}, lower = 0.05, 
                                 upper=340)$root}

# Sequence of values of Kurtosis (NOT excess kurtosis)
n          <- 1000
range.kurt <- seq(from = 10, to = 300, length.out=n)

# Sequence of values of r
r.end      <- 0.95
# maximal ranges for all distributions
range.Fre  <- seq(from=0.2453,  to = r.end, length.out=n)
range.LogN <- seq(from=0.2022,  to = r.end, length.out=n)
range.Gam  <- seq(from=0.2030,  to = r.end, length.out=n)
range.Wei  <- seq(from=0.191,  to = r.end, length.out=n)
range.t    <- seq(from=0.2028,  to = r.end, length.out=n)
range.Par    <- seq(from=0.3069,  to = r.end, length.out=n)

# additional ranges such that the kurtosis exists
range.2.Fre  <- seq(from=0.2453,  to = 0.4013, length.out=n)
range.2.t    <- seq(from=0.2028,  to = 0.2928, length.out=n)
range.2.Par    <- seq(from=0.3069,  to = 0.4648, length.out=n)

# range of RQW
range.RQW <- seq(from=0.3,  to = .995, length.out=n)

# Colors to be used for plots
my.col = c("brown3","cornflowerblue", "blueviolet", "darkgoldenrod1", "darkgreen",
           "pink")

# Initialise sequences of parameters and kurtosis
logN.para  <- rep(NA, n)
k.logN.vec <- rep(NA, n)
RQW.logN.vec <- rep(NA, n)

Wei.para  <- rep(NA, n)
k.Wei.vec <- rep(NA, n)
RQW.Wei.vec <- rep(NA, n)

Gam.para  <- rep(NA, n)
k.Gam.vec <- rep(NA, n)
RQW.Gam.vec <- rep(NA, n)

t.para  <- rep(NA, n)
k.t.vec <- rep(NA, n)
RQW.t.vec <- rep(NA, n)

Par.para  <- rep(NA, n)
k.Par.vec <- rep(NA, n)
RQW.Par.vec <- rep(NA, n)

Fre.para  <- rep(NA, n)
k.Fre.vec <- rep(NA, n)
RQW.Fre.vec <- rep(NA, n)

# Finding the "k:kurtosis" corresponding to specific values of "RQW"
for (i in 1:n){
  logN.para[i]  <- RQW.logN.inv(range.RQW[i])
  k.logN.vec[i] <- k.logN(logN.para[i])
  
  Wei.para[i]   <- RQW.Wei.inv(range.RQW[i])
  k.Wei.vec[i]  <- k.Wei(Wei.para[i])
  
  Gam.para[i]   <- RQW.Gam.inv(range.RQW[i])
  k.Gam.vec[i]  <- k.Gam(Gam.para[i])

}

# Plots: kurtosis versus RQW (log scale)
pdf(file = "k.vs.RQW0875.pdf", width = 8, height = 6) 
y.range <- c(min(log(k.logN.vec),log(k.Wei.vec), log(k.Gam.vec)), 
             max(log(k.logN.vec),log(k.Wei.vec), log(k.Gam.vec)))
par(mar = c(5, 5, 2, 1))
plot(range.RQW,   log(k.logN.vec), type = 'l', lwd=3, ylim = y.range,
     lty = 1, col=my.col[1], xlab="RQW", ylab="log(kurtosis)", cex.lab = 1.5,
     cex.axis = 1.25)
lines(range.RQW, log(k.Wei.vec), lwd = 3, lty=2, col=my.col[2])
lines(range.RQW, log(k.Gam.vec), lwd = 3, lty=3, col=my.col[3])
legend('topleft', ncol=1, legend=c("LogN","Weibull","Gamma"), col=my.col[1:3], 
       lty = 1:3, lwd = 5, cex = 1.5)
# Saving plot
dev.off()

# Finding the "k:kurtosis" corresponding to specific values of "r*"
for (i in 1:n){
  logN.para[i]  <- r.logN.inv(range.LogN[i])
  k.logN.vec[i] <- k.logN(logN.para[i])
  
  Wei.para[i]   <- r.Wei.inv(range.Wei[i])
  k.Wei.vec[i]  <- k.Wei(Wei.para[i])
  
  Gam.para[i]   <- r.Gam.inv(range.Gam[i])
  k.Gam.vec[i]  <- k.Gam(Gam.para[i])
  
  Par.para[i]   <- r.Par.inv(range.2.Par[i])
  k.Par.vec[i]  <- k.Par(Par.para[i])
  
  Fre.para[i]   <- r.Fre.inv(range.2.Fre[i])
  k.Fre.vec[i]  <- k.Fre(Fre.para[i]) 
  
  t.para[i]     <- r.t.inv(range.2.t[i])
  k.t.vec[i]    <- k.t(t.para[i])
}


# Plots: kurtosis versus r (log scale)
pdf(file = "k.vs.r.pdf", width = 8, height = 6) 

y.range <- c(min(log(k.logN.vec),log(k.Wei.vec), log(k.Gam.vec), log(k.Par.vec), 
                 log(k.Fre.vec)), 
             max(log(k.logN.vec),log(k.Wei.vec), log(k.Gam.vec), log(k.Par.vec), 
                 log(k.Fre.vec)))
par(mar = c(5, 5, 2, 1))
plot(range.LogN,   log(k.logN.vec), type = 'l', lwd=3, ylim = y.range,
     lty = 1, col=my.col[1], xlab="r*", ylab="log(kurtosis)", cex.lab = 1.5, 
     cex.axis = 1.25)
lines(range.Wei, log(k.Wei.vec), lwd = 3, lty=2, col=my.col[2])
lines(range.Gam, log(k.Gam.vec), lwd = 3, lty=3, col=my.col[3])
legend('topleft', ncol=1, legend=c("LogN","Weibull","Gamma"), col=my.col[1:5], 
       lty = 1:5, lwd = 5, cex = 1.5)
# Saving plot
dev.off()

# Finding the "RQW" corresponding to specific values of "r*"
for (i in 1:n){
  logN.para[i]    <- r.logN.inv(range.LogN[i])
  RQW.logN.vec[i] <- RQW.logN(logN.para[i])
  
  Wei.para[i]     <- r.Wei.inv(range.Wei[i])
  RQW.Wei.vec[i]  <- RQW.Wei(Wei.para[i])
  
  Gam.para[i]     <- r.Gam.inv(range.Gam[i])
  RQW.Gam.vec[i]  <- RQW.Gam(Gam.para[i])
  
  Par.para[i]     <- r.Par.inv(range.Par[i])
  RQW.Par.vec[i]  <- RQW.Par(Par.para[i])
  
  Fre.para[i]     <- r.Fre.inv(range.Fre[i])
  RQW.Fre.vec[i]  <- RQW.Fre(Fre.para[i]) 
  
  t.para[i]       <- r.t.inv(range.t[i])
  RQW.t.vec[i]    <- RQW.t(t.para[i])
}

# Plots: RQW versus r
pdf(file = "RQW.vs.r.pdf", width = 8, height = 6) 
y.range <- c(0.0,1)
par(mar = c(5, 5, 2, 1))
plot(range.LogN,   (RQW.logN.vec), type = 'l', lwd=3, ylim = y.range,
     lty = 1, col=my.col[1], xlab="r*", ylab="RQW", cex.lab = 1.5, cex.axis = 1.25)
lines(range.Wei, (RQW.Wei.vec), lwd = 3, lty=2, col=my.col[2])
lines(range.Gam, (RQW.Gam.vec), lwd = 3, lty=3, col=my.col[3])
lines(range.Par, (RQW.Par.vec), lwd = 3, lty=4, col=my.col[4])
lines(range.Fre, (RQW.Fre.vec), lwd = 3, lty=5, col=my.col[5]) 
lines(range.t, RQW.t.vec, lwd = 3, lty=6, col=my.col[6]) 
legend('topleft', ncol=1, 
       legend=c("LogN","Weibull","Gamma","Pareto","Frechet", "Student"), 
       col=my.col[1:6], lty = 1:6, lwd = 5, cex = 1.5)
# Saving plot
dev.off()

