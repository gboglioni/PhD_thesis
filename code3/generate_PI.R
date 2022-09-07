# Function to simulate PIBD variables (U1,U2,U3) under various examples

# Arguments of function 'generate.PI '
#   n:         Sample size
#   type:      Key word for the name of the example
#   alpha:     Parameter (used in some examples only)

generate.PI <- function(n, type = 'indep', alpha=1){
  
  # Two independent samples of U(0,1)  
  U1 <- runif(n, 0, 1)
  U2 <- runif(n, 0, 1)
  
  # when "mixing" a PIBD structure with mutual independence,
  # "n.dep" is the number of observations stemming from the dependent structure
  if (type == 'bernstein' || type == 'triangles'){
    n.dep <- rbinom(1, n, alpha)
  }
  
  # mutual independence
  if (type == 'indep'){
    U3 <- runif(n, 0, 1)
  }
  
  # Examples 3.1 and 3.6
  if (type == 'triangles'){
    U3 <- c((U1[1:n.dep]+U2[1:n.dep])%%1, runif(n-n.dep))
  }
  
  # Example 3.2
  if (type == 'bernstein'){
    X <- rbinom(n.dep, 1, 1/2)
    Y <- rbinom(n.dep, 1, 1/2)
    Z <- -abs(X-Y)+1       
    U1 <- c((X+runif(n.dep))/2, runif(n-n.dep))
    U2 <- c((Y+runif(n.dep))/2, runif(n-n.dep))
    U3 <- c((Z+runif(n.dep))/2, runif(n-n.dep))
  }
  
  # Example 3.3
  if (type == 'tetra'){
    n1 <- rbinom(1, n, 1/2)
    U3 <- rep(NA, n)
    for (i in 1:n1){
      ifelse(U1[i]+U2[i] >= 1, U3[i]<- U1[i]+U2[i]-1, U3[i] <- 1-(U1[i]+U2[i]))
    }
    for (i in (n1+1):n){
      ifelse(U1[i]-U2[i] >= 0, U3[i]<- 1 - (U1[i]-U2[i]), U3[i] <- 1+(U1[i]-U2[i]))
    }
  }
  
  # Example 3.4
  if (type == 'cosine'){
    U3 <- pbeta((cos(2*pi*(U1+U2))+1)/2, 1/2, 1/2)
  }
  
  # Example 3.5
  if (type == 'copula'){
    Tval  <- runif(n)
    a  <- alpha*(1-2*U1)*(1-2*U2)
    U3 <- (1+a - sqrt((1+a)^2-4*a*Tval))/(2*a)
  }
  
  # Example 3.7
  if (type == 'bernstein4D'){
    M1 <- rmultinom(n, 1, rep(1/alpha, alpha))
    M2 <- rmultinom(n, 1, rep(1/alpha, alpha))
    M3 <- rmultinom(n, 1, rep(1/alpha, alpha))
    M4 <- rmultinom(n, 1, rep(1/alpha, alpha))
    X1 <- apply(M1 == M3, 2, prod)
    X2 <- apply(M1 == M4, 2, prod)
    X3 <- apply(M2 == M3, 2, prod)
    X4 <- apply(M2 == M4, 2, prod)
    U1 <- X1 * runif(n, (alpha-1)/alpha,1) + ifelse(X1==0,1,0)*runif(n,0,(alpha-1)/alpha)
    U2 <- X2 * runif(n, (alpha-1)/alpha,1) + ifelse(X2==0,1,0)*runif(n,0,(alpha-1)/alpha)
    U3 <- X3 * runif(n, (alpha-1)/alpha,1) + ifelse(X3==0,1,0)*runif(n,0,(alpha-1)/alpha)
    U4 <- X4 * runif(n, (alpha-1)/alpha,1) + ifelse(X4==0,1,0)*runif(n,0,(alpha-1)/alpha)
    data <- as.data.frame(cbind(U1,U2,U3,U4))
  }
  
  # Return data.frame
  if(type != 'bernstein4D'){
    data <- as.data.frame(cbind(U1,U2,U3))
  }
  return(data)
}
