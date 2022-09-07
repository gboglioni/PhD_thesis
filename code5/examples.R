# Specific values of randF, indA and ell for all examples

# Example 5.1
randF <- function(m, beta) rlnorm(m, 0, beta)
indA  <- function(x) ifelse(x >= 1, TRUE, FALSE)
piid.generator(randF = function(m) randF(m, beta = 2), indA = indA, ell = 2L)

# Example 5.3
r.ex6 <- function(rand, ell) {
  if (rand < 1 / (2 * ell)) {
    res <- -1L
  } else if (rand < 1 / ell) {
    res <- 1L
  } else if (rand < 1 / (2 * ell) + 1 / 2) {
    res <- -2L
  } else {
    res <- 2L
  }
  return(res)
}
randF <- function(m, ell) sapply(runif(m), FUN = r.ex6, ell = ell)
indA  <- function(x) ifelse((x == 1L) | (x == -1L), TRUE, FALSE)
piid.generator(randF = function(m) randF(m, ell = 2L), indA = indA, ell = 2L)

# Example 5.4
randF <- function(m) runif(m, -1, 1)
indA  <- function(x, ell) ifelse((x >= -1 / ell) & (x <= 1 / ell), TRUE, FALSE)
piid.generator(randF = randF, indA = function(x) indA(x, ell = 2L), ell = 2L)

# Example 5.5
randF <- function(m, ell) as.integer(2 * rbinom(m, 1, 1 / ell) - 1)
indA  <- function(x) ifelse(x == 1L, TRUE, FALSE)
piid.generator(randF = function(m) randF(m, ell = 10L), indA = indA, ell = 10L)

# Example 5.6
F  <-  function(x, ell, sigma) sum(c(1 - 1 / ell, 1 / ell) *
                         pnorm(x, mean = c(-1 / ell, 1 - 1 / ell),
                               sd = c(sigma, sigma)))
Finv  <- function(p, ell, sigma = 1/(4*ell)){ 
  G = function(x) F(x, ell, sigma) - p
  return(uniroot(G, c(-100,100))$root)
}
r.ex9 <- function(rand, ell, sigma) {
  if (rand < 1 / ell) {
    res <- rnorm(1, 1 - 1 / ell, sigma)
  } else {
    res <- rnorm(1, -1 / ell, sigma)
  }
  return(res)
}
randF <- function(m, ell, sigma = 1/(4*ell)) 
          sapply(runif(m), FUN = r.ex9, ell = ell, sigma = sigma)
indA  <- function(x, ell, sigma = 1/(4*ell)) 
          ifelse(x >= Finv(1 - 1 / ell, ell = ell, sigma = sigma), TRUE, FALSE)
piid.generator(randF = function(m) randF(m, ell = 3L), 
                indA = function(x) indA(x, ell = 3L), ell = 3L)

# Example 5.7
randF <- function(m, mu, sigma) rnorm(m, mu, sigma)
indA  <- function(x, mu) ifelse(x >= mu, TRUE, FALSE)
piid.generator(randF = function(m) randF(m, mu = 2, sigma = 1),
               indA = function(x) indA(x, mu = 2), ell = 2L)
