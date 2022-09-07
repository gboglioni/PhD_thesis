# Source function 'piid.generator' (which allows to generate PIBD variables) 
source('piid-generator.R')

# Functions needed for the generation of PIBD Poisson(log(2))
randF <- function(m) rpois(m, log(2))
indA  <- function(x) ifelse(x >= 1, TRUE, FALSE)

# Histogram for the empirical distribution of q.hat for the PIBD sequence
m     <- 10 # Controls the sample size: sample size = m(m-1)/2
B     <- 10000 # Number of samples
q     <- rep(NA, B)
set.seed(1)
for (i in 1:B){
  sample <- piid.generator(randF = randF, indA = indA, ell = 2L, m = m)
  q[i]   <- sum(sample == 0)/length(sample)
}

# Plot histogram
pdf(file = "piid_poisson_q_dist.pdf", width = 14, height = 8)
par(mar=c(5,5,1,1))
hist(q, probability = T, breaks = 20, main = '', xlab =  bquote(hat(q)), 
     cex.lab = 1.75, cex.axis = 1.5, col = 'lightgreen')
abline(v = 1/2, col = 'black', lwd = 4, lty = 2)
legend("topleft", NULL, ncol = 1, cex = 2, legend=c("true q"), 
       col=c("black"), lwd = 4, lty = c(2))
dev.off()

# Function 'boot.piid' produces bootstrapped re-samples of a dataset ('data')
# It computes a statistic (% of 0's) from the bootstrap samples and creates a CI
# It also checks whether the value (q=0.5) is inside the bootstrap CI

# Arguments of function 'boot.piid':
#   m:          For simulated data, controls the sample size n, with n = m(m-1)/2
#   B:          Number of bootstrapped samples
#   data:       Original sample of data to use 
#               (if unspecified, a random mutually independent sample is generated)
#   histo:      TRUE to plot a histogram of the bootstrapped statistics
#   br:         Number of 'bins' for the histograms
#   alpha:      Level for the confidence interval

boot.piid <- function(m = 3, B = 500, data = '', 
                      histo = F, br=30, alpha = 0.10){
  
  # For unspecified 'data', generate a random sample of independent Pois(log(2))
  if(length(data)==1){data <- rpois(m*(m-1)/2, log(2))}
  
  # Sample size of bootstrapped samples
  boot.size=length(data)
  
  # Matrix of bootstrapped samples (col = observations, row = bootstrap replicates)
  S <- matrix(NA, ncol = boot.size, nrow = B)
  for (i in 1:B){
    S[i,] <- sample(data, replace = T, size = boot.size)
  }
  
  # Bootstrapped replicates of the statistic
  q <- rep(NA, B) # initiate vector of bootstrapped statistics
  for (i in 1:B){
      q[i] = sum(S[i,] == 0)/boot.size
  }
  
  # Indicator equal to 1 if the real value is outside the bootstrap CI (L, 1)
  ind <- ifelse(1/2 < quantile(q, alpha), 1,0)
  
  # Histogram of the bootstrapped statistic 'q'
  if (histo == T){
    par(mar=c(5,5,1,1))
    hist(q, probability = T, breaks = br, main = '', xlab = bquote(hat(q)), 
         cex.lab = 1.75, cex.axis = 1.5, col = 'lightgreen')
    abline(v = quantile(q, alpha), col = 'brown3', lwd = 4)
    abline(v = 1/2, col = 'black', lwd = 4, lty = 2)
    legend("topright", NULL, ncol = 1, cex = 2, legend=c("L", "true q"), 
           col=c("brown3", "black"), lwd = 4, lty = c(1,2))
  }
  return(ind)
  
} # End function

# Call function 'boot.piid' to obtain the histogram of one bootstrapped sample
# (and save the plot)
# IID sample
pdf(file = "bootstrapCI_iid.pdf", width = 14, height = 8)
set.seed(1)
boot.piid(m=10, B = 5000, histo=T)
dev.off()
# PIBD sample
pdf(file = "bootstrapCI_piid.pdf", width = 14, height = 8)
set.seed(1)
boot.piid(data = piid.generator(randF = randF, indA = indA, ell = 2L, m = 10), 
          B = 5000, histo = T)
dev.off()

# Function that calls 'boot.piid' a large number of times 
# (to estimate the empirical level of the BCI)
boot.level <- function(nb.trials=10000, B=5000, m=3, indep=TRUE, seed=1){
  set.seed(seed)
  count<-0
  for (i in 1:nb.trials){
    if(indep){sample <- rpois(m*(m-1)/2, lambda = log(2))} #iid
    else{
      sample  <- piid.generator(randF = randF, indA = indA, ell = 2L, m = m) #PIBD
    }
    count <- count + boot.piid(m = m, B = B, data = sample, alpha = 0.1)
  }
  count/nb.trials
}

# Launch simulations
# (warning: takes a few hours)

# IID case 
boot.level(nb=10000, B=5000, m=10)
boot.level(nb=10000, B=5000, m=45)

# PIBD case
boot.level(nb=10000, B=5000, m=10, indep=F)
boot.level(nb=10000, B=5000, m=45, indep=F)
