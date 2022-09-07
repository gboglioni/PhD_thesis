###########################################
## Incidence matrix of projective planes ##
###########################################

require(nortest) # for the Anderson-Darling and Pearson tests at the end

## the graph is bipartite, (q+1)-regular, has girth 6, diameter 3 
## and 2*(q^2 + q + 1) vertices
## (q must be a prime power)

## Sigma matrix

sigma <- function (u, q) {
  res <- (u * (matrix(rep(0:(q - 1), q), nrow = q) 
            + t(matrix(rep(0:(q - 1), q), nrow = q)))) %% q
  return(replace(res, res == 0, q))
}

## Position matrix of A

pm <- function (A, q) {
  res <- matrix(as.integer(A == 1), nrow = q)
  for (i in 2:q) {
    res <- cbind(res, matrix(as.integer(A == i), nrow = q))
  }
  return(res)
}

## Position matrix of a family F

pmf <- function (q) {
  res <- pm(sigma(1, q), q)
  for (i in 2:(q - 1)) {
    res <- rbind(res, pm(sigma(i, q), q))
  }
  return(res)
}

## Incidence matrix of the bipartite graph

inc_mat <- function (q) {
  res1 <- rbind(pmf(q), do.call(cbind, replicate(q, diag(q), simplify=FALSE)), 
                pm(matrix(rep(1:q, q), nrow = q), q), rep(0, q ^ 2))
  res2 <- rbind(t(pm(matrix(rep(1:(q + 1), q), nrow = q), q + 1)), rep(1, q + 1))
  return(cbind(res1, res2))
}

## Adjacency matrix

adj_mat <- function (q) {
  n <- q ^ 2 + q + 1
  zz <- matrix(0, nrow = n, ncol = n)
  im <- inc_mat(q)
  return(rbind(cbind(zz, im), cbind(t(im), zz)))
}

## Computation of the standardized rv Z

stand_rand <- function (q) {
  n_vert <- 2 * (q ^ 2 + q + 1) ## number of vertices in the graph
  n_edges <- n_vert * (q + 1) / 2 ## number of edges in the graph
  vec <- rbinom(n_vert, 1, 1 / 2)
  mat <- matrix(rep(vec, n_vert), n_vert)
  res <- 1 - ((mat + t(mat)) %% 2)
  x <- res * adj_mat(q) ## generates rvs on the edges
  return((sum(x) / 2 - n_edges / 2) / sqrt(n_edges / 4))
}

## Monte-Carlo histogram

q <- 2 ^ 6
sim <- 5000
z <- rep(0, sim)
  
for (i in 1:sim) {
  z[i] <- stand_rand(q)
}
hist(z)
qqnorm(z, pch = 1, frame = FALSE)
qqline(z, col = "steelblue", lwd = 2)

shapiro.test(z)
ad.test(z)
pearson.test(z)

