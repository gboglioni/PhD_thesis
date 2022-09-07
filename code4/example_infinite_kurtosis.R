# Load packages
library('moments')

# Function to generate the PIBD sample
piid.seq <- function(m=10){
  e   <- runif(1)
  eta <- runif(1)
  Z <- rep(NA, m)
  for (j in 1:m){Z[j] <- (eta + j*e)%%1}
  return(Z)
}

# Function to generate samples of S, for various choices of sample sizes 'n'
samples.of.S <- function(n = c(10, 100, 1000, 10000), B = 10000){
  set.seed(1)
  # Matrix to contain all samples of S
  # (rows represent samples: from one row to the next, 'n' increases)
  S    <- matrix(NA, nrow = length(n), ncol = B)
  kurt <- rep(NA, length(n))
  for (i in 1:length(n)){ # for all rows
    for (j in 1:B){ # for all columns within that row
      U <- piid.seq(n[i])
      S[i,j] <- sum(U)
    }
    kurt[i] <- kurtosis(S[i,])
  }
  plot(log(n), log(kurt), type="b", lty=1, col = "brown3", lwd = 6, 
       xlab="log(n)", ylab="log(kurtosis)", cex.lab = 1.5, cex.axis = 1.25)
  return(summary(lm(kurt ~ n)))
}

# Launch function and save plot
# (warning: running time is many hours)
pdf(file = "increasing_kurtosis.pdf", width = 14, height = 8)
par(mar=c(5,5,1,1))
samples.of.S(n = c(10, 100, 1000, 10000, 100000), B = 3000000)
dev.off()
