# Structural stability (in norme L1)

Omega.L1 <- function(alpha){
  if(nrow(alpha)==1){return(NA)}
  Omega <- NA
  alpha <- alpha %*% diag(1/colSums(alpha))
  try({
    Omega <- determinant(alpha,logarithm = TRUE)$modulus[1]
  },silent = TRUE)
  Omega <- Omega / log(10)
  #Omega <- sum(log(eigen(alpha)$values))
  return(Omega)
}


r.centroid.L1 <- function(alpha){
  if(nrow(alpha)==1){return(NA)}
  alpha.L1 <- alpha %*% diag(1/colSums(alpha))
  r.c <- rowSums(alpha.L1)
  r.c <- r.c / sum(r.c)
  return(r.c)  
}


theta.L1 <- function(alpha,r){
  if(nrow(alpha)==1){return(NA)}
  r.c <- r.centroid.L1(alpha)
  cos.theta <- sum(r * r.c) / (sqrt(sum(r^2))*sqrt(sum(r.c^2)))
  theta <- acos(cos.theta) * 180 / pi
  return(theta)  
}



# Structural stability (in norme L2)
library(mvtnorm)

Omega.L2 <- function(alpha){
  n <- nrow(alpha)
  Omega <- NA
  alpha <- alpha %*% diag(1/colSums(alpha)) * n
  try({
    Sigma <-solve(t(alpha) %*% alpha)
    d <- pmvnorm(lower = rep(0,n), upper = rep(Inf,n), mean = rep(0,n), sigma = Sigma)
    if ((attr(d,"msg") == "Normal Completion") && (attr(d,"error")/d[1] < 1)){
      Omega <- log10(d[1]) + n * log10(2)}
  },silent = TRUE)
  return(Omega)
}


r.centroid.L2 <- function(alpha){
  alpha.L2 <- alpha %*% diag(1/sqrt(colSums(alpha^2)))
  r.c <- rowSums(alpha.L2)
  r.c <- r.c / sqrt(sum(r.c^2))
  return(r.c)  
}


theta.L2 <- function(alpha,r){
  # browser()
  r.c <- r.centroid.L2(alpha)
  cos.theta <- sum(r * r.c) / (sqrt(sum(r^2))*sqrt(sum(r.c^2)))
  theta <- acos(cos.theta) * 180 / pi
  return(theta)  
}