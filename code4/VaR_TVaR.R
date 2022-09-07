# Source needed function
source('generate_PI.R')

# Function to find the empirical TVaR of a sample ("data")
emp.TVaR <- function(data, alpha = c(0.75,0.95)){
  l    <- length(alpha)
  TVAR <- rep(NA, l)
  VAR  <- quantile(data, alpha)
  for (i in 1:l){TVAR[i] <- mean(data[data > VAR[i]])}
  return(TVAR)
}

# Function 'VaR.TVaR.S' estimates VaR, TVaR of X1+X2+X3 (for many PIBD examples)
# Arguments of this function are:
#   B          Number of simulated samples
#   alpha:     Vector of levels (for which VaR and TVaR are computed)
#   dep.type:  Name of the dependence example (as defined in function 'generate.PI')
#   param:     Value of the parameter within the example (if one is needed)
#   relative:  TRUE to divide all values by their corresponding values 
#              under mutual independence
#   seed:      seed of the random generation

VaR.TVaR.S <- function(B=200000, alpha = c(0.70, 0.90, 0.95, 0.99, 0.995), 
                       dep.type = 'indep', param=1, relative=T, seed=1){
  
  # Set seed
  set.seed(seed)
  
  # Data with uniform margins
  Udata <- as.matrix(generate.PI(n = B, type = dep.type, alpha = param))
  
  # Transform the data to different margins (always such that E[X]=1, Var[X]=1)
  U.data  <- qunif(Udata, min = 1-sqrt(3), max = 1+sqrt(3))
  N.data  <- qnorm(Udata, mean = 1, sd = 1)
  G.data  <- qgamma(Udata, shape = 1, rate = 1)
  LN.data <- qlnorm(Udata, -log(2)/2, sqrt(log(2)))
  
  # Get samples of "S", for all margins
  S.piid.U  <- apply(U.data, MARGIN = 1, sum)
  S.piid.N  <- apply(N.data, MARGIN = 1, sum)
  S.piid.G  <- apply(G.data, MARGIN = 1, sum)
  S.piid.LN <- apply(LN.data, MARGIN = 1, sum)

  # VaR
  VaR.1  <- quantile(S.piid.U, alpha)-3
  VaR.2  <- quantile(S.piid.N, alpha)-3
  VaR.3  <- quantile(S.piid.G, alpha)-3
  VaR.4  <- quantile(S.piid.LN, alpha)-3

  #TVaR
  TVaR.1 <- emp.TVaR(S.piid.U, alpha)-3
  TVaR.2 <- emp.TVaR(S.piid.N, alpha)-3
  TVaR.3 <- emp.TVaR(S.piid.G, alpha)-3
  TVaR.4 <- emp.TVaR(S.piid.LN, alpha)-3
  
  result <- c(VaR.1, VaR.2, VaR.3, VaR.4, TVaR.1, TVaR.2, TVaR.3, TVaR.4)
  ifelse(relative, return(round(result/VaR.TVaR.S(B, relative=F, seed=2),2)), 
         return(result))
}

# Run the function for all examples
# (warning: running time is many hours)
VaR.TVaR.S(B=10^7, dep.type = 'triangles')
VaR.TVaR.S(B=10^7, dep.type = 'bernstein')
VaR.TVaR.S(B=10^7, dep.type = 'tetra')
VaR.TVaR.S(B=10^7, dep.type = 'cosine')
VaR.TVaR.S(B=10^7, dep.type = 'copula', param = 1)
VaR.TVaR.S(B=10^7, dep.type = 'copula', param = -1)


