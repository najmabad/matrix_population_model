# Author: Najma Bader <najmabader@hotmail.it>
# Code implementation to reproduce results of 
# Caswell, H. (2009). Stage, age and individual stochasticity in demography. 
# Oikos, 118(12), 1763–1782. https://doi.org/10.1111/j.1600-0706.2009.17620.x

library(zeallot)

# Model parameters
s1 <- 0.9026497
s2 <- 0.96769866
g2 <- 0.12517818
s3 <- 1
g3 <- 0.28986564
s4 <- 0.84956954
s5 <- 1
f <- 0.13


# Population growth
populationGrowthProjection <- function(times, C, I, M, MM, P, t, TOT, lambda, plot = FALSE, log = FALSE) {
  for (i in 1:times) {
    C[i+1] <- f * M[i]
    I[i+1] <- s1 * C[i] + s2 * (1 - g2) * I[i] 
    M[i+1] <- s2 * g2 * I[i] + s3 * (1 - g3) * M[i] + s5 * P[i] 
    MM[i+1] <-  s3 * g3 * M[i]
    P[i+1] <- s4 * MM[i]
    t[i+1] <- i
    TOT[i+1] <-  C[i+1] + I[i+1] + M[i+1] + MM[i+1] + P[i+1] 
    lambda[i+1] <-  TOT[i+1]/TOT[i]
  }
  if (plot) {
    if (log) {
      c(C, I, M, MM, P) %<-% list(log(C), log(I), log(M), log(MM), log(P))
      suffix <- "(log)"
    } else {
      suffix <- ""
    }
    plot(t, C, type="o", ylab=paste("Number of individuals", suffix), xlab="Years", pch=20, lwd=2, main= "Simulated NAW Population Growth")
    lines(t, I, type="o", col="red")
    lines(t, M,  type="o",col="blue")
    lines(t, MM,  type="o",col="green")
    lines(t, P,  type="o",col="purple")
    legend("top", legend = c(paste("Calf", suffix), paste("Immature", suffix), paste("Mature", suffix), paste("Mother", suffix), paste("Post-breading mother", suffix)),
           col = c("black", "red", "blue", "green", "purple"), lty = 1:2, cex = 0.8)
    plot(lambda, ylab = "Years", xlab = "Growth Rate")
  }
  return(list(C, I, M, MM, P, t, TOT, lambda))
}


# Scenario 1: we start with 500 calves and project over 50 years
populationGrowthProjection(
  times  = 50,
  C      = 500,
  I      = 0,
  M      = 0,
  MM     = 0,
  P      = 0,
  TOT    = 500,
  lambda = 0,
  t      = 0,
  plot   = TRUE,
  log    = TRUE
)


# Scenario 2: we start with an heterogeneous distribution of the initial population
# and project over 50 years
populationGrowthProjection(
  times  = 50,
  C      = 5,
  I      = 45,
  M      = 150,
  MM     = 50,
  P      = 250,
  TOT    = 500,
  lambda = 0,
  t      = 0,
  plot   = TRUE,
  log    = TRUE
)


# Scenario 3: we start with 500 calves and project over 500 years
populationGrowthProjection(
  times  = 500,
  C      = 5,
  I      = 45,
  M      = 150,
  MM     = 50,
  P      = 250,
  TOT    = 500,
  lambda = 0,
  t      = 0,
  plot   = TRUE,
  log    = TRUE
)


# Scenario 4: we start with an heterogeneous distribution of the initial population
# and project over 500 years

populationGrowthProjection(
  times  = 500,
  C      = 5,
  I      = 45,
  M      = 150,
  MM     = 50,
  P      = 250,
  TOT    = 500,
  lambda = 0,
  t      = 0,
  plot   = TRUE,
  log    = TRUE
)

################## Perturbation analysis ##################


theta <- c(s1, s2, s3, s4, g2, g3)

invertNonNull <- function(vector) {
  vectorInv <- 0
  for (i in 1:length(vector)){
    if (vector[i] != 0) {
      vectorInv[i] <- 1 / vector[i]
    } else {
      vectorInv[i] <- 0
    }
  }
  return(vectorInv)
}

## Fundamental matrix N 
# Sensitivity 
generateN <- function(theta, stacked = TRUE){
  U <- matrix(
    c(0, theta[1], 0, 0, 0,
      0, theta[2] * (1-theta[5]), theta[2] * theta[5], 0, 0,
      0, 0, theta[3] * (1-theta[6]), theta[3] * theta[6], 0,
      0, 0, 0, 0, theta[4],
      0, 0, theta[3], 0, 0 ), 
    nrow = 5,
    ncol = 5)
  N <- solve(diag(5) - U)
  if (stacked) {
    return(c(vec(N)))
  } else {
    return(N)
  }
}

sensitivityN <- jacobian(generateN, theta, method="Richardson")


# Elasticity
elasticityN <- diag(invertNonNull(generateN(theta)))  %*% sensitivityN  %*% diag(theta)


# Plot
plot(
  x    = c(1, 2, 3, 4, 5, 6),
  y    = c(elasticityN[4,]), 
  pch  = 21,  
  bg   = "black", 
  xaxt = "n", 
  xlab = "theta", 
  ylab = "Elasticity", 
  main = "Elasticity of lifetime reproductive events"
)
axis(1, at = c(1, 2, 3, 4, 5, 6), labels = expression(s1, s2, s3, s4, g2, g3))


## Variance of the fundamental matrix N
# Sensitivity 
generateV <- function(theta, stacked = TRUE){
  N <- generateN(theta, stacked = FALSE)
  V <- (2*diag(diag(N)) - diag(5)) %*% N - hadamard.prod(N, N)
  if (stacked) {
    return(c(vec(V)))
  } else {
    return(V)
  }
}

sensitivityV <- jacobian(generateV, theta, method="Richardson")


# Elasticity
elasticityV <- diag(invertNonNull(generateV(theta)))  %*% sensitivityV  %*% diag(theta)


# Plot 
plot(
  x    = c(1, 2, 3, 4, 5, 6),
  y    = c(elasticityV[4,]),
  pch  = 21,  
  bg   = "black", 
  xaxt = "n", 
  xlab = "theta", 
  ylab = "Elasticity", 
  main = "Elasticity of variance in lifetime reproductive events"
)
axis(1, at = c(1, 2, 3, 4, 5, 6), labels = expression(s1, s2, s3, s4, g2, g3))


## Standard deviation of the fundamental matrix N
stddevV <- sqrt(generateV(theta))


## Expected longevity
# Sensitivity
generateExpLongevity <- function(theta, stacked = TRUE){
  N <- generateN(theta, stacked = FALSE)
  expLongevity <- t(c(1,1,1,1,1)) %*% N
  if (stacked) {
    return(c(vec(expLongevity)))
  } else {
    return(expLongevity)
  }
}

sensitivityExpLongevity <- jacobian(generateExpLongevity, theta, method="Richardson")


# Elasticity
elasticityExpLongevity <- diag(invertNonNull(generateExpLongevity(theta)))  %*% sensitivityExpLongevity  %*% diag(theta)


# Plot
plot(
  x    = c(1, 2, 3, 4, 5, 6),
  y    = c(elasticityExpLongevity[1,]), 
  pch  = 21,  
  bg   = "black", 
  xaxt = "n", 
  xlab = "theta", 
  ylab = "Elasticity", 
  main = "Elasticity of life expectancy (calf)"
)
axis(1, at = c(1,2,3,4,5,6), labels = expression(s1, s2, s3, s4, g2, g3))


## Variance of longevity
N <- generateN(theta, stacked = FALSE)
varLongevity <-  t(c(1,1,1,1,1)) %*% N %*% (2*N - diag(5)) - hadamard.prod(expLongevity,expLongevity)

# Sensitivity
generateVarLongevity <- function(theta, stacked = TRUE){
  N <- generateN(theta, stacked = FALSE)
  expLongevity <- t(c(1,1,1,1,1)) %*% N
  varLongevity <-  t(c(1,1,1,1,1)) %*% N %*% (2*N - diag(5)) - hadamard.prod(expLongevity, expLongevity)
  if (stacked) {
    return(c(vec(varLongevity)))
  } else {
    return(varLongevity)
  }
}

sensitivityVarLongevity <- jacobian(generateVarLongevity, theta, method="Richardson")

# Elasticity
elasticityVarLongevity  <- diag(invertNonNull(generateVarLongevity(theta)))  %*% sensitivityVarLongevity %*% diag(theta)

# Plot
plot(
  x    = c(1, 2, 3, 4, 5, 6),
  y    = c(elasticityVarLongevity[1,]), 
  pch  = 21,
  bg   = "black", 
  xaxt = "n", 
  xlab = "theta", 
  ylab =" Elasticity", 
  main = "Elasticity of variance in life expectancy (calf)"
)
axis(1, at = c(1, 2, 3, 4, 5, 6), labels <- expression(s1, s2, s3, s4, g2, g3))


## Standard deviation of longevity
stddevLongevity <- sqrt(varLongevity)


## Net Reproductive Rate
# Sensitivity
generateNetReproductiveRate <- function(theta){
  F <- matrix(
    c(0, 0, 0, 0, 0,
      0, 0, 0, 0, 0,
      0.13, 0, 0, 0, 0, 
      0, 0, 0, 0, 0,
      0, 0, 0, 0, 0), 
    nrow = 5,
    ncol = 5
  ) 
  N <- generateN(theta, stacked = FALSE)
  FN <- F %*% N
  return(FN[1,1])
}

sensitivityNetReproductiveRate <- c(jacobian(generateNetReproductiveRate, theta, method="Richardson"))


# Plot
plot(
  x    = c(1, 2, 3, 4, 5, 6),
  y    = sensitivityNetReproductiveRate, 
  pch  = 21,  
  bg   = "black", 
  xaxt = "n", 
  xlab = "theta", 
  ylab = "Sensitivity", 
  main = "Sensitivity of the Net Reproductive Rate"
)
axis(1, at = c(1, 2, 3, 4, 5, 6), labels = expression(s1, s2, s3, s4, g2, g3))
