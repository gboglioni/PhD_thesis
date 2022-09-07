# Generator of pairwise independent observations
piid.generator <- function(m = 3, randF = rnorm, 
                           indA = function(x) ifelse(x <= 0, FALSE, TRUE),
                           ell = 2) {
  
  # Check that the value of 'ell' provided is coherent 
  #  with the function 'indA' provided
  # This check is valid only for moderate values of 'ell'
  if (round(1 / mean(indA(randF(10 ^ 5)))) != ell) 
    warning("Is 'ell' consistent with your A?")
  
  # Sample size
  n <- choose(m, 2)
  
  # Generate the 'initial' multinomial sample of size m
  M <- rmultinom(m, 1, rep(1 / ell, ell))
  
  # Find all possible pairs out of the m multinomials
  # Use those pairs to create the 'n' D's in Eq. (2.4)
  combin <- combn(1:m, 2)
  D <- apply(M[, combin[1,]] == M[, combin[2,]], 2, all)
  D <- as.integer(D)
  # Compute the number of 1's among the D's
  pN <- sum(D)
  
  # Generate pN r.v.s with distribution F restricted to A 
  # and (n - pN) r.v.s with distribution F restricted to A^c
  nU <- 0
  nV <- 0
  U <- rep(NA, n - pN)
  V <- rep(NA, pN)
  while((nU < n - pN) | (nV < pN)) {
    W <- randF(1)
    indAW <- indA(W)
    if (indAW & (nV < pN)) {
      nV <- nV + 1
      V[nV] <- W
    } else if ((indAW == 0) & (nU < n - pN)) {
      nU <- nU + 1
      U[nU] <- W
    }
  }
  # Return the resulting random generated variables  
  X <- rep(NA, n)
  X[which(D == 0)] <- U
  X[which(D == 1)] <- V
  return(X)
}

