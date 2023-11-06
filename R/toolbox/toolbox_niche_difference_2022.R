require(mvtnorm)


Omega_L2 <- function(alpha){
  n <- nrow(alpha)
  D <- diag(1/sqrt(diag(t(alpha)%*%alpha)))
  alpha <- alpha %*% D
  Sigma <-solve(t(alpha) %*% alpha)
  d <- pmvnorm(lower = rep(0,n), upper = rep(Inf,n), mean = rep(0,n), sigma = Sigma)
  out <- log10(d[1]) + n * log10(2)
  return(out) 
}


Omega_L1 <- function(alpha){
  D <- diag(1/colSums(alpha))
  alpha <- alpha %*% D
  out <- (determinant(alpha, logarithm = TRUE)$modulus[1]) / log(10)
  return(out)
}


rc_L2 <- function(alpha){
  D <- diag(1/sqrt(diag(t(alpha)%*%alpha)))
  alpha <- alpha %*% D
  out <- rowSums(alpha)
  out <- out / sqrt(sum(out^2))
  return(out)
}

rc_L1 <- function(alpha){
  D <- diag(1/colSums(alpha))
  alpha <- alpha %*% D
  out <- rowSums(alpha)
  out <- out / sum(out)
  return(out)
}


theta_L2 <- function(r,alpha){
  rc <- rc_L2(alpha)
  out <- acos(sum(rc*r)/(sqrt(sum(r^2))*sqrt(sum(rc^2))))*180/pi
  return(out)
}


theta_L1 <- function(r,alpha){
  rc <- rc_L1(alpha)
  out <- acos(sum(rc*r)/(sqrt(sum(r^2))*sqrt(sum(rc^2))))*180/pi
  return(out)
}


angle_to_border <- function(r,alpha){
  n <- length(r)
  if(n == 1){return(list(eta = NA, p = NA))}
  D <- diag(1/colSums(alpha)) * n
  alpha <- alpha %*% D
  r <- r / sum(r)
  p <- matrix(NA,nrow = n,ncol = n)
  eta <- rep(NA,n)
  for (i in 1:n){
    A <- alpha[,-i]
    p_t <- A%*% solve(t(A) %*% A, t(A)) %*% r
    cos_eta <- sum(r * p_t) / sqrt(sum(r^2)*sum(p_t^2))
    eta[i] <- acos(cos_eta) * 180 / pi
    p[,i] <- p_t
  }
  out <- list(eta = eta, p = p)
  return(out)
}


angle_to_corner <- function(r,alpha){
  n <- length(r)
  if(n == 1){return(NA)}
  psi <- rep(NA,n)
  for (i in 1:n){
    rc <- alpha[,i]
    cos_psi <- sum(r * rc) / sqrt(sum(r^2)*sum(rc^2))
    psi[i] <- acos(cos_psi) * 180 / pi
  }
  return(psi)
}
  

