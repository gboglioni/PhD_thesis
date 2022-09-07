# Generate a random sample as in the sequence in Janson1988 (from Remark 2)
janson.seq <- function(n=10){
e <- runif(1)
eta <- runif(1)
Z <- rep(NA,n)
for (j in 1:n){
    Z[j] <- (cos(2*pi*(eta+j*e))+1)/2
}
return(Z)
}

# Find distribution of the standardised maxima in the Janson1988 sequence
# We do this in a function, which returns the histogram and ECDF of -log(-M), 
# (compared to the 'prediction' from the Fisher-Tippett THM)
# We use the log scale, because the distribution itself is very heavy-tailed

# Arguments of function 'EVT.comparison'
#   n           Vector of sample sizes
#   B:          Number of simulations (i.e. number of samples generated)
#   low/up:     Lower and upper limits for the histogram
#   seed:       Seed to reproduce simulation results
#   br:         Number of breaks for the histogram
#   compare:    TRUE for a comparison of the CDF for two sample sizes 
#               (the last two values of vector n are used)
#   plot.means: TRUE to plot -log(-mean(H.n)) as function of log(n) 

EVT.comparison <- function(n=c(10,1000,10000), B=50000, low = -15, up=15, 
                           seed=1, br=200, compare=F, plot.means=T){
set.seed(seed)
# Number of different sample sizes
n.n <- length(n)
# Matrix to contain the samples
M   <- matrix(NA, nrow = B, ncol = n.n)

# Generation of the samples (of maxima)
# (different rows represent different samples)
# (different columns represent different sample sizes)
for (i in 1:B){
  U     <- janson.seq(tail(n,1))
  for (j in 1:n.n){
    M[i,j]  <- (max(U[1:n[j]])-1)*(2*n[j]/pi)^2
  }
}
# Log transformation
logH  <- -log(-M)

# Mean of M.n (for all sample sizes)
mean    <- apply(M, 2, mean)

# Plotting (and saving plot) log-transformed mean of M.n versus sample size
if(plot.means){
  pdf(file = paste("evt_mean_B", B, "n", n[n.n], ".pdf", sep=""), 
      width = 14, height = 8)
  par(mfrow=c(1,1))
  par(mar=c(5,5,1,1))
  plot(log(n), -log(-mean), type="b", lty=1, col = "brown3", lwd = 6, 
       xlab="log(n)", ylab=bquote(-log(-E(H[n]))), cex.lab = 1.5, cex.axis=1.25)
  dev.off()
}

# Plots of empirical distribution (histogram & ECDF)
# (vs. what is 'predicted' by Fisher-Tippet)
pdf(file = paste("evt_df_B", B, "n", n[n.n], ".pdf", sep=""), width=16, height=8)
par(mfrow=c(1,2))
par(mar=c(3,5,1,1))
hist(logH[,n.n], xlim = c(low, up), breaks = br, prob=T, 
     cex.lab = 1.5, cex.axis =1.25, main = "", xlab="")
curve(0.5*exp(-exp(-x/2))*exp(-x/2), lwd = 2, n = 500, 
      col = "brown3", add  = TRUE) # Density of -Log(-Weibull(1/2))

plot(ecdf(logH[,n.n]),  cex.lab = 1.5, cex.axis =1.25, do.points = F, 
     col.01line = NULL, xlim = c(low, up), main = "", ylab="CDF", xlab="")
if(compare){
  plot(ecdf(logH[,n.n-1]), do.points = F, add = T, col = 'darkblue', lty = 2, 
       xlim = c(low, up), main = "", ylab="CDF", xlab="")
}
else{
  curve(exp(-exp(-x/2)), lwd = 2, n = 500, 
        col = "brown3", add  = TRUE) # Density of -Log(-Weibull(1/2))
  legend("bottomright", NULL, ncol = 1, cex = 1.5, legend=c("F-T THM"), 
         col=c("brown3"), lty = 1, lwd = 3)
}
dev.off()
}

# (warning: running time is many hours)
EVT.comparison(n=10^(1:5), B=3000000, seed=1)
