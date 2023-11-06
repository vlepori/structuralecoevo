require(mvtnorm)
require(kernlab)
#require(cubature)

#input parameters:
#alpha = compeition strenght matrix (assumed to be positive and globally stable)
#r = vector of intrinsic growth rates


#solid angle in topology L2 of the frasibility domain ( = niche difference)
Omega_L2 <- function(alpha){
  n <- nrow(alpha)
  Sigma <-solve(t(alpha) %*% alpha)
  d <- pmvnorm(lower = rep(0,n), upper = rep(Inf,n), mean = rep(0,n), sigma = Sigma)
  out <- log10(d[1]) + n * log10(2)
  return(out) 
}


#solid angle in topology L1 of the frasibility domain ( = niche difference)
Omega_L1 <- function(alpha){
  n <- nrow(alpha)
  out <- determinant(alpha,logarithm = TRUE)$modulus[1] - sum(log(colSums(alpha)))
  out <- out / log(10)
  return(out)
}


#normalize the columns of the alpha matrix to unit norm L2 (L2 norm of columns = 1)
norm_L2 <- function(alpha){
  D <- diag(1/sqrt(diag(t(alpha)%*%alpha)))
  alpha_n <- alpha %*% D
  return(alpha_n)
}


#normalize the columns of the alpha matrix to unit norm L1 (L1 norm of columns = 1)
norm_L1 <- function(alpha){
  D <- diag(1/colSums(alpha))
  alpha_n <- alpha %*% D 
  return(alpha_n)
}


#vector defining the centroid of the feasibility domain in topology L2 (norm of vector = 1)
r_centroid_L2 <- function(alpha){
  # browser()
  n <- nrow(alpha)
  D <- diag(1/sqrt(diag(t(alpha)%*%alpha)))
  alpha_n <- alpha %*% D # rescaled matrix: l2 norm of cols = 1
  r_c <- rowSums(alpha_n)
  r_c <- r_c / sqrt((sum(r_c^2)))
  r_c <- t(t(r_c))
  return(r_c)
}

#vector defining the centroid of the feasibility domain in topology L1
r_centroid_L1 <- function(alpha){
  n <- nrow(alpha)
  D <- diag(1/colSums(alpha))
  alpha_n <- alpha %*% D # rescaled matrix: l1 norm of cols = 1
  r_c <- rowSums(alpha_n) /n 
  r_c <- t(t(r_c))
  return(r_c)
}


#deviation from centroid in topology L2 (in degree)
theta_L2 <- function(alpha,r){
  r_c <- r_centroid_L2(alpha)
  out <- acos(sum(r_c*r)/(sqrt(sum(r^2))*sqrt(sum(r_c^2))))*180/pi
  return(out)
}


#deviation from centroid in topology L1 (in degree)
theta_L1 <- function(alpha,r){
  r_c <- r_centroid_L1(alpha)
  out <- acos(sum(r_c*r)/(sqrt(sum(r^2))*sqrt(sum(r_c^2))))*180/pi
  return(out)
}

##################################################
### evenness
##################################################

evenness <- function(N){
  N <- N[N>0]
  n <- length(N)
  p <- N/sum(N)
  out <- - sum(p*log10(p))/log10(n)
  return(out)
}

##################################################
### alpha, r, and K
##################################################

alpha_matrix <- function(mu,sigma,a){
  n <- dim(mu)[1]
  d <- dim(mu)[2]
  Mdist <- as.matrix(dist(mu))
  Ma <- t(t(a)) %*% t(a)
  Msigma <- t(t(rep(1,n))) %*% t(sigma)
  Msigma2 <- Msigma^2 + t(Msigma^2)
  out <- ( Ma * Msigma^d * t(Msigma)^d * (2*pi)^d  / (sqrt(Msigma2))^d ) * exp(-Mdist^2 / (2*Msigma2))
  return(out)
}



alpha_matrix_bounded <- function(mu,sigma,a){
  n <- dim(mu)[1]
  d <- dim(mu)[2]
  Mdist <- as.matrix(dist(mu))
  Ma <- t(t(a)) %*% t(a)
  Msigma <- t(t(rep(1,n))) %*% t(sigma)
  Msigma2 <- Msigma^2 + t(Msigma^2)
  alpha <- Ma * exp(-Mdist^2 / (2*Msigma2))
  for (i in 1:n){
    for (j in i:n){
      sigmaij <- sqrt(sigma[i]^2*sigma[j]^2/Msigma2[i,j])
      muij <- (mu[i,] * sigma[j]^2 + mu[j,] * sigma[i]^2)/Msigma2[i,j]
      alpha[i,j] <- alpha[i,j] * prod(pnorm(rep(1,d),mean=muij,sd=sigmaij) - pnorm(rep(0,d),mean=muij,sd=sigmaij)) * (2*pi)^(d/2) * sigmaij^d
      alpha[j,i] <- alpha[i,j]
    }
  }
  return(alpha)
}



alpha_matrix_bounded_m_fast <- function(mu,sigma,a,alpha,sp){
  n <- dim(mu)[1]
  d <- dim(mu)[2]
  Mdist <- as.matrix(dist(mu))
  Ma <- t(t(a)) %*% t(a)
  Msigma <- t(t(rep(1,n))) %*% t(sigma)
  Msigma2 <- Msigma^2 + t(Msigma^2)
  alpha_f <- Ma * exp(-Mdist^2 / (2*Msigma2))
  for (i in 1:n){
    sigmaij <- sqrt(sigma[i]^2*sigma[sp]^2/Msigma2[i,sp])
    muij <- (mu[i,] * sigma[sp]^2 + mu[sp,] * sigma[i]^2)/Msigma2[i,sp]
    alpha[i,sp] <- alpha_f[i,sp] * prod(pnorm(rep(1,d),mean=muij,sd=sigmaij) - pnorm(rep(0,d),mean=muij,sd=sigmaij)) * (2*pi)^(d/2) * sigmaij^d
    alpha[sp,i] <- alpha[i,sp] 
  }
  return(alpha)
}


alpha_matrix_bounded_m <- function(mu,sigma,a,alpha,sp){
  n <- dim(mu)[1]
  d <- dim(mu)[2]
  ff <- function(x,mu1,mu2,sigma1,sigma2){
    exp(-sum((x-mu1)^2)/(2*sigma1^2)-sum((x-mu2)^2)/(2*sigma2^2))
  }
  for (i in 1:n){
    alpha[i,sp] <- a[i]*a[sp]*adaptIntegrate(ff, rep(0,d), rep(1,d), tol=1e-4, mu1 = mu[i,] , mu2 = mu[sp,], sigma1 = sigma[i], sigma2 = sigma[sp])$integral
    alpha[sp,i] <- alpha[i,sp] 
  }
  return(alpha)
}




r_eff_vector_bounded <- function(m,mu,sigma,a,r){
  n <- dim(mu)[1]
  d <- dim(mu)[2]
  out <- rep(NA,n)
  for (i in 1:n){
    out[i] <- -m[i] + r * a[i] * prod(pnorm(rep(1,d),mean=mu[i,],sd=sigma[i]) - pnorm(rep(0,d),mean=mu[i,],sd=sigma[i])) * (2*pi)^(d/2) * sigma[i]^d
  }
  return(out)
}


r_eff_vector_bounded_m <- function(m,mu,sigma,a,r,r_eff,sp){
  d <- dim(mu)[2]
  r_eff[sp] <- -m[sp] + r * a[sp] * prod(pnorm(rep(1,d),mean=mu[sp,],sd=sigma[sp]) - pnorm(rep(0,d),mean=mu[sp,],sd=sigma[sp])) * (2*pi)^(d/2) * sigma[sp]^d
  return(r_eff)
}


m_vector_bounded <- function(r_eff,mu,sigma,a,r){
  n <- dim(mu)[1]
  d <- dim(mu)[2]
  out <- rep(NA,n)
  for (i in 1:n){
    out[i] <- r * a[i] * prod(pnorm(rep(1,d),mean=mu[i,],sd=sigma[i]) - pnorm(rep(0,d),mean=mu[i,],sd=sigma[i])) * (2*pi)^(d/2) * sigma[i]^d - r_eff[i]
  }
  return(out)
}


K_vector <- function(alpha,r){
  out <- r / diag(alpha)
  return(out)
}

niche_overlap_alpha <- function(alpha){
  out <- diag(1/diag(alpha)) %*% alpha
  return(out)
}

#K_vector <- function(m,sigma,a,r){
#  d <- dim(mu)[2]
#  out <- (-m + r*a*(2*pi)^(d/2)*sigma^d) / (a^2 * pi^(d/2) * sigma^d)
  #out <- (-m + r*a*(2*pi)^(d/2)) / (a^2 * pi^(d/2) / sigma^d)
#  return(out)
#}


#m_vector <- function(K,sigma,a,r){
#  d <- dim(mu)[2]
#  out <- r*a*(2*pi)^(d/2)*sigma^d - K * (a^2 * pi^(d/2) * sigma^d)
  #out <- r*a*(2*pi)^(d/2) - K * (a^2 * pi^(d/2) / sigma^d)
#  return(out)
#}


##################################################
### LCP
##################################################

LCP <- function(alpha,r){
  
  ff <- sqrt(sum(r^2))
  
  r <- r / ff
  
  out <- solve(alpha,r)
  
  if ( prod(out > 0) == 0 ){
  
    n <- length(r)
  
    H <- 2*t(alpha)
    c <- -matrix(r)
  
    A <- alpha
    b <- matrix(r)
    l <- matrix(rep(0,n))
    u <- matrix(rep(10000000,n))
    r <- matrix(rep(10000000,n))
  
    sv <- ipop(c,H,A,b,l,u,r,sigf = 12)
    out <- t(t(sv@primal))
  }
    
  return(out*ff)

}

##################################################
### fitness equivalence
##################################################


r_invasion <- function(alpha,r){
  n <- length(r)
  out <- matrix(NA,n,1)
  for (i in 1:n){
    out[i] <- r[i] - alpha[i,-i] %*%  solve(alpha[-i,-i],r[-i])  
  }  
  return(out)
}



A_matrix <- function(alpha){
  n <- dim(alpha)[1]
  A <- matrix(1,nrow = n, ncol = n)
  for (i in 1:n){
    A[i,-i] <- -alpha[i,-i] %*% solve(alpha[-i,-i])
  }  
  return(A)
}

##################################################
### max-min function (angle)
##################################################

#angle in degree
min_deviation <- function(alpha,r){
  n <- length(r)
  if(n<3){return(NA)}
  theta_min <- rep(NA,n)
  D <- diag(1/sqrt(diag(t(alpha)%*%alpha)))
  alpha_n <- alpha %*% D
  alpha_n <- alpha
  for (i in 1:n){
    alpha2 <- alpha_n[,-i]
    alpha3 <- alpha2[,-1]
    v <- alpha2[,1]
    alpha4 <- alpha3 - v %*% matrix(1,1,n-2)
    
    angle <- function(x) { w <- v + alpha4 %*% x
                          out <- sum(w*r)/(sqrt(sum(w^2))*sqrt(sum(r^2)))
                          return(acos(out))}
    out <- optim(rep(1/n,n-2),angle, method = ifelse(n==3,"BFGS","Nelder-Mead"), control = list(abstol = 1e-8, reltol=1e-8) )
    theta_min[i] <- out$value*180/pi  
  }
  return(theta_min)
}

#angle in degree
max_deviation <- function(alpha,r){
  D <- diag(1/sqrt(diag(t(alpha)%*%alpha)))
  alpha_n <- alpha %*% D
  r_n <- r / sqrt(sum(r^2))
  theta_max <- acos(t(alpha_n) %*% r_n)*180/pi
  return(theta_max)
}


#in degree
min_center <- function(alpha){
  n <- dim(alpha)[1]
  D <- diag(1/sqrt(diag(t(alpha)%*%alpha)))
  alpha_n <- alpha %*% D
  sd_theta <- function(x){ v <- alpha_n %*% x
                      out <- min_deviation(alpha_n,v)
                      return(sd(out))}
  out <- optim(rep(1/n,n),sd_theta, control = list(abstol = 1e-8, reltol=1e-8))  
  x <- out$par
  vc <- alpha_n %*% x
  return(vc)
}
  

##################################################
### min-max function (distance on simplex)
##################################################

min_distance <- function(alpha,K){
  n <- dim(alpha)[1]
  if(n<3){return(NA)}
  distances <- rep(NA,n)
  alpha_n <- norm_L1(alpha)
  K <- K / sum(K)
  nv <- n_vectors(alpha_n)
  for (i in 1:n){
    vi <- alpha_n[,(i %% n)+1]
    dim(vi) <- c(n,1)
    ni <- nv[,i]
    dim(ni) <- c(n,1)
    distances[i] <- abs(sum((K-vi) * ni))
  }
  return(distances)
}


  delta_alpha <- function(alpha){
  n1 <- dim(alpha)[1]
  n2 <- dim(alpha)[2]
  v1 <- alpha[,1]
  dim(v1) <- c(n1,1)
  d_alpha <- alpha[,-1]
  dim(d_alpha) <- c(n1,n2-1)
  A <- rep(v1,n2-1)
  dim(A) <- c(n1,n2-1)
  d_alpha <- d_alpha - A
  return(d_alpha)
}


n_vectors <- function(alpha){
  n <- dim(alpha)[1]
  vo <- o_vector(alpha)
  vn <- matrix(NA,n,n)
  for (i in 1:n){
    vi <- alpha[,i]
    alpha_c <- alpha[,-i]
    d_alpha <- delta_alpha(alpha_c)
    A <- t(cbind(vo, d_alpha))
    lambda <- c(1,solve(A[,-1],-A[,1]))
    lambda <- lambda / sqrt(sum(lambda^2))
    lambda <- lambda * sign(sum((vi- alpha_c[,1]) * lambda))
    vn[,i] <- lambda
  }
  return(vn)
}


o_vector <- function(alpha){
  d_alpha <- delta_alpha(alpha)
  A <- t(d_alpha)
  lambda <- solve(A[,-1],-A[,1])
  vo <- c(1,lambda)
  vo <- vo / sqrt(sum(vo^2))
  return(vo)
}


in_center <- function(alpha){
  alpha_n <- norm_L1(alpha)
  n <- dim(alpha_n)[1]
  v1 <- alpha_n[,1]
  dim(v1) <- c(1,n)
  vo <- o_vector(alpha_n)
  vn <- n_vectors(alpha_n)
  d_alpha <- delta_alpha(alpha_n)
  A <- cbind(t(vn) %*% d_alpha, rep(-1,n))
  B <- rep(v1,n)
  dim(B) <- c(n,n)
  C <- alpha_n[,c(2:n,1)]
  D <- diag(t((C - B)) %*% vn)
  out <- solve(A,D)
  lambda <- out[1:(n-1)]
  r <- out[n]
  vin <- t(v1) + d_alpha %*% lambda
  return(as.vector(vin))
}



















##################################################



set_vector_on_triangle <- function(beta,v){
  
  S <- dim(beta)[1]
  d_beta <- delta_beta(beta)  
  
  v1 <- beta[,1]
  dim(v1) <- c(S,1)
  
  A <- cbind(v, -d_beta) 
  
  lambda <- solve(A,v1)[1]  
  v <- v * lambda
  
  return(v)
}


distance_farthest_border <- function(beta,v){
  
  S <- length(v)
  vs <- set_vector_on_triangle(beta,v)
  
  A <- rep(vs,S)
  dim(A) <- c(S,S)
  
  distances <- sqrt(colSums ((A - beta)^2))
  
  return(distances)
}











##################################################
### sample alpha matricies functions
##################################################


sample_niche_space_rho <- function(n,d,rho){
  p <- runif(n*d)
  dim(p) <- c(n,d)
  A <- as.matrix(dist(p))
  f <- function(lambda){
    b <- exp(-lambda*A)
    S <- dim(A)[1]
    out <- ((sum(b) - S)/((S-1)*S) - rho)^2
    return(out)
  }
  lambda <- optimize(f,interval = c(0,100),maximum = F,tol = 10^-14)$minimum
  alpha <- exp(-lambda*A)
  return(alpha)
}



sample_niche_space_OmegaL2 <- function(n,d,Omega_t){
  test <- TRUE
  while (test == TRUE){
    p <- runif(n*d)
    dim(p) <- c(n,d)
    A <- as.matrix(dist(p))
    f <- function(lambda){
      b <- exp(-lambda*A)
      Omega <- OmegaL2(b)
      out <- (Omega-Omega_t)^2
      return(out)
    }
    out_o <- optimize(f,interval = c(0.001,2000),maximum = F,tol = 10^-14)
    val <- out_o$objectiv
    test <- (val>1e-10)
  }
  lambda <- out_o$minimum
  alpha <- exp(-lambda*A)
  return(alpha)
}


sample_alpha_rho <- function(n,d,rho){
  det_alpha <- 1e-8
  while (det_alpha < 2e-8){
    p <- runif(n*d)
    dim(p) <- c(n,d)
    A <- as.matrix(dist(p))
    alpha_in <- runif(n) * 2 + 0.75
    f <- function(lambda){
      b <- diag(alpha_in) %*% exp(-lambda*A)
      S <- dim(A)[1]
      out <- ((sum(b) - S)/((S-1)*S) - rho)^2
      return(out)
    }
    lambda <- optimize(f,interval = c(0,100),maximum = F,tol = 10^-14)$minimum
    alpha <- diag(alpha_in) %*% exp(-lambda*A)
    det_alpha <- det(alpha)
  }
  return(alpha)
}


