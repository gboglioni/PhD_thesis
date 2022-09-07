# Load external functions
source("3D_visualisation_functions.R")

# Set path to where results (MSE values and graphs) should be saved
graphs.direct <- "C:/Users/z3519303/Dropbox/Apps/Overleaf/PhD_Thesis/code3"
mse.direct    <- "C:/Users/z3519303/Dropbox/Apps/Overleaf/PhD_Thesis/code3"

# Function that runs simulations
# Arguments:
#     n.sim:    number of samples generated
#     my.n:     sample size
#     r1,r2:    parameters for the max and min values of the radius 'r'
#     spacing:  interval length between different values of 'r'
tuning.sim <- function(n.sim=100, my.n=200, r1=0.1, r2=0.2, spacing=0.002){
set.seed(1)
r <- seq(from = r1, to = r2, by = spacing)
l <- length(r)
# Matrix to contain MSE:
# rows represent different samples, columns represent different r's
MSE <- matrix(NA, nrow=n.sim, ncol=l) 

for (k in 1:n.sim){
  X <- runif(my.n, 0, 1)
  Y <- runif(my.n, 0, 1)
  Z <- runif(my.n, 0, 1)
  my.data <- as.data.frame(cbind(X,Y,Z))
  # For different value of 'r', estimate f(h) at every point in the sample
  # Then compute the MSE
  for (j in 1:l){
    fh        <- pts.concentration(my.data, rad=r[j], return.h = F)
    MSE[k,j]  <- sqrt(mean((fh - 1/2)^2))
  }  
}
# Mean MSE (for every r, averaged across all simulated samples)
mse     <- apply(MSE, 2, mean)
results <- cbind(r, mse, c(NA, 100*(mse[-1] - mse[-l])/mse[-l]))

# Save results
write.csv(results, 
          file = paste(mse.direct, "/", "nsim", n.sim, "n", my.n, ".csv", sep=""))

# Plots of MSE versus 'r'
pdf(paste(graphs.direct, "/", "nsim", n.sim, "n", my.n, ".pdf", sep = ""), 
    width = 12, height = 10)
par(mfrow=c(1,1))
par(mar = c(5, 5, 2, 2))
plot(r, results[,2], type = 'l', lty = 1, lwd = 3, col='brown3', xlab="r", 
     ylab="MSE", cex.lab = 2, cex.axis = 2)
dev.off()

return(list(results, r[which.min(results[,2])], which.min(results[,2])))
}

# Launch simulations (takes many hours)
tuning.sim(n.sim=2000, my.n = 200,   r1=0.15, r2=.35)
tuning.sim(n.sim=2000, my.n = 300,   r1=0.10, r2=.30)
tuning.sim(n.sim=2000, my.n = 400,   r1=0.10, r2=.30)
tuning.sim(n.sim=2000, my.n = 500,   r1=0.10, r2=.40)
tuning.sim(n.sim=2000, my.n = 1000,  r1=0.07, r2=.25)
tuning.sim(n.sim=1000, my.n = 2000,  r1=0.05, r2=.23)
tuning.sim(n.sim=1000, my.n = 3000,  r1=0.05, r2=.23)


