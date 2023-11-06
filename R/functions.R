add_if_even <-function(x, e = 1e-3){
  
  y=rep(NA,length(x))
  
  for(i in 1:length(x)){
    y[i]=if(i%%2==0){x[i]+e}else{x[i]}
  }
  return(y)
  
}

make_mc = function(mc1,mc2){

    mc1 = mc1 %>%
      rename(step = "x1", time = "x2", x = "x3") %>%
      long2nested() %>%
      as_tibble %>%
      select(x)
    mc2 = mc2 %>%
      rename(step = "x1", time = "x2", x = "x3") %>%
      long2nested() %>%
      as_tibble %>%
      select(x)
    mc = bind_rows(mc2,mc1,.id = "id") #!
    return(mc)
}

duplicate_rows <- function(df){
  
  out = df[1,]
  
  for (i in 2:nrow(df)){
    out=  rbind(out,df[i-1,])
    out=  rbind(out,df[i,])
  }
  
  out =  rbind(out,df[nrow(df),])
  
  return(out)
}


niche_diff <- function(alpha){
  if(!all(dim(alpha)==c(2,2))){return(NA)}
  1 - sqrt((alpha[1,2]*alpha[2,1])/(prod(diag(alpha))))
}

fitness_diff <- function(alpha,r){
  if(!all(dim(alpha)==c(2,2))){return(NA)}
  if(!length(r)==2L){return(NA)}
  out = (r[1]/r[2])*(sqrt((alpha[2,2]*alpha[2,1])/(alpha[1,1]*alpha[1,2])))
 return(round(out,digits = 9))
}


save_plots_rda <- function(path){
  
  all = mget(ls(envir = globalenv()),
             envir = globalenv())
  
  isplot = lapply(all, is.ggplot) %>%  unlist
  
  out=all[isplot]
  
  saveRDS(out,file = path)
  
  message(paste("saving",length(out),"plot(s):"))
  for(i in 1:length(out)){
    message(names(out)[i])
  }
  
  return(NULL)
}


Nstar <- function(alpha,r){
  
  if(length(r)==1){
    return(r / alpha)
  }
  N = tryCatch(solve(alpha,r),error = function(e) NA)
  
  if(any(is.na(N))){return(NA)}
  if(all(N>1e-9)){
    return(N)
  }
  return(NA)
}

long2nested = function(data, key = "step"){
  # browser()
  idx = unique(unlist(data[key]))
  
  olist = list()
  times = rep(NA, length(idx))
  
  # browser()
  
  for(i in 1:length(idx)){
    tmp = filter(data, step == idx[i])
    # print(i)
    times[i] = unique(tmp$time)
    olist[[i]] =  tmp %>%  pull(x)
  }
  
  out = tibble::tibble(step = unique(data$step), 
                       time = times, 
                       x = olist)
  
  return(out)
}

plotmat <- function(d){
  d <- data.frame(d)
  d <- cbind(row = rownames(d), d)
  d <-
    d %>%
    tidyr::gather(key = "column", value = "value", 2:ncol(.)) %>%
    mutate(row = as.integer(as.character(row)),
           column = as.integer(substring(column, 2)))
  pl<-
    ggplot(d)+
    geom_tile(aes(x=column,y=row,fill = value))+
    coord_fixed()+
    scale_fill_gradient2()+
    scale_y_continuous(trans = "reverse"
    )
  
  pl  
}

df2list = function(df, return.vec = FALSE){
  
  stopifnot(nrow(df)==1L)
  
  names = colnames(df)
  
  values = as.vector(as.matrix(df[1,]))
  
  names(values) = names
  
  if(return.vec){return(values)}
  
  return(as.list(values))
  
}

mode_dist <- function(x) { # Mode function
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

lst2str<-function(l){paste(names(l),unlist(l), collapse = ", ", sep = "=")}
obj2str<-function(o){paste(deparse(substitute(o)),o, sep = "=")}

alpha_mat <- function(mu_v, alphafun, aa){
  
  l <- length(mu_v)
  # if(length(sigma_v) == 1L){sigma_v <- rep(sigma_v, l)} #recycle
  # stopifnot(length(mu_v) == length(sigma_v))
  
  out <- matrix(rep(NA, l^2), nrow = l)
  
  for(f in 1:l){
    for(c in 1:l){
      
      arg.list <- c(x1 = mu_v[f],x2 = mu_v[c], aa)
      afc <- do.call(alphafun,arg.list)
      out[f,c] <- afc
    }
  }
  return(out)
}

r_vec <- function(mu_v, sigma_v, rfun, aa = NULL){
  
  l <- length(mu_v)
  if(length(sigma_v) == 1L){sigma_v <- rep(sigma_v, l)}
  stopifnot(length(mu_v) == length(sigma_v))
  
  out <- rep(NA, l)
  
  for(x in 1:l){
    
    arg.list <- list(mu_1 = mu_v[x], sigma_1 = sigma_v[x], aa)
    out[x] <- do.call(rfun, arg.list)
    
  }
  return(out)
}


adiag <- function(m){
  
  stopifnot(dim(m)[1]==dim(m)[2])
  
  l = nrow(m)
  out = rep(NA,l)
  for(i in 1:l){
    out[i]= m[i,l-i+1]
  }
  
  return(out)
  
}


alpha_tilde <- function(alphamat){ # efficient
  alphamat/diag(alphamat)
}

alpha_tilde2 <- function(alphamat){ #sanity check
  d<- dim(alphamat)
  out <- matrix(NA, nrow = d[1], ncol = d[2])
  for (i in 1:d[1]) {
    for (j in 1:d[2]) {
      out[i,j] = alphamat[i,j]/alphamat[i,i]
    }
  }
  return(out)
}

alpha_bar <- function(alphamat, K) { # diagonal should always be 1
  d<- dim(alphamat)
  out <- matrix(NA, nrow = d[1], ncol = d[2])
  for (i in 1:d[1]) {
    for (j in 1:d[2]) {
      out[i,j] = (K[j]/K[i])*(alphamat[i,j]/alphamat[i,i])
    }
  }
  return(out)
}

#convenience functions

submat <- function(indices, matrx, ind2=NULL){
  
  if(is.null(ind2)){return(matrx[indices,indices])}
  
  return(matrx[indices,ind2])
  
}

subvect <- function(indices, vctr){
  vctr[indices]
}


avg_nn <- function(x){
  xs <- sort(x)
  l <- length(xs)
  
  if(l==1){return(NA)}
  
  nn <- rep(NA, l)
  
  nn[1] <- xs[2]-xs[1]
  for(i in 2:(l-1)){
    nn[i] <- min(xs[i]-xs[i-1],xs[i+1]-xs[i])
  }
  nn[l] <- xs[l]-xs[l-1]
  
  return(mean(nn))
}


PIP <- function(rs, alphas, mu_v, return.df=F, log = F){
  mn = length(mu_v)
  expand.grid(FI=1:mn, CI=1:mn) %>% 
    mutate(Fo = mu_v[FI], Co = mu_v[CI]) %>% 
    mutate(rj = rep(rs, mn),
           ri = rep(rs, each = mn)) %>% 
    mutate(alphaii = sapply(X = .$CI,FUN = submat, matrx= alphas),
           alphaij = mapply(FUN = submat, indices = .$FI, ind2 = .$CI, MoreArgs = list(matrx = alphas))) %>% 
    mutate(w = rj - (ri/alphaii)*alphaij) -> d.f
  
  if(return.df)return(d.f)
  
  d.f %>% 
    ggplot(aes(x = Co, y= Fo))+
    # geom_tile(aes(fill=w))+
    {if(log)geom_tile(aes(fill=alogv(w)))}+
    {if(!log)geom_tile(aes(fill=w))}+
    stat_contour(aes(z = w),breaks = 0)+
    labs(x="Resident",y="Mutant")+
    scale_fill_gradient2()+
    coord_equal()
  
  
}

PIP2 <- function(range,rfun=NULL,alphafun=NULL, res = 100){
  
  if(is.null(rfun)){rfun=rfun}
  if(is.null(alphafun)){rfun=alphafun}
  
  a=range[1]
  b=range[2]
  
  df<-
  expand.grid(FI=seq(a,b,length.out = res), 
              CI=seq(a,b,length.out = res))
  
  # browser()
  
  df$rf = sapply(df$FI, rfun)
  df$rc = sapply(df$CI, rfun)
  
  df$aFF = mapply(alphafun,df$FI,df$FI)
  df$aCC = mapply(alphafun,df$CI,df$CI)
  df$aFC = mapply(alphafun,df$FI,df$CI)
  df$aCF = mapply(alphafun,df$CI,df$FI)
  
  df$wm = df$rc - (df$rf/df$aFF)*df$aFC
  
    
  df %>% 
    ggplot(aes(x = FI, y= CI))+
    # geom_tile(aes(fill=w))+
    # {if(log)geom_tile(aes(fill=alogv(w)))}+
    geom_tile(aes(fill=wm))+
    stat_contour(aes(z = wm),breaks = 0, color = "black")+
    labs(x="Resident",y="Mutant",fill = "Mutant fitness")+
    scale_fill_gradient2()+
    coord_equal()+
    theme_classic()
  
}






Jeigs <- function(alphas, rs, Ns=NULL){
  
  ns = nrow(alphas)
  if(ns==1){return(NA)}
  if(is.null(Ns)) {Ns <- solve(alphas)%*%rs}
  if(any(is.na(Ns))){ return(rep(NA,ns))}
  D = matrix(data = 0, nrow = ns, ncol = ns)
  diag(D) = Ns
  J = D%*%-alphas
  lambdas <-eigen(J, only.values = T)$values # real part of eigenvalues < 0   <->   stability
  return(lambdas)
}



lt3 <- function(data, columns=1:ncol(data), color=NULL, keep = 'ltdiag') { # options lt, ltdiag, all.
  # data <- upgrade_scatmat_data(data)
  data.choose <- data[columns]
  dn <- data.choose[sapply(data.choose, is.numeric)]
  factor <- data[color]
  p <- ncol(dn)
  q <- nrow(dn)
  newdata <- as.data.frame(matrix(NA, nrow = q*p*p, ncol = 6+ncol(factor)))
  newdata[5:6]<-data.frame(character(q),character(q),stringsAsFactors = FALSE)
  r <-1
  for (i in 1:p) {
    for (j in 1:p) {
      # browser()
      newdata[r:(r+q-1),] <- cbind(dn[[i]], dn[[j]], i, j, colnames(dn)[i], colnames(dn)[j], factor,stringsAsFactors = FALSE)
      r <- r+q
    }
  }
  colnames(newdata) <- c("xvalue", "yvalue", "xslot", "yslot", "xlab", "ylab", colnames(factor))
  
  rp <- data.frame(newdata)
  
  rp$xvalue <- (as.numeric(as.character(rp$xvalue)))
  rp$yvalue <- (as.numeric(as.character(rp$yvalue)))
  rp$xlab <- as.factor(rp$xlab)
  rp$ylab <- factor(rp$ylab, levels = levels(rp$xlab))
  
  rp<-
    switch (keep,
            'ltdiag' = {dplyr::filter(rp, as.integer(xlab)<=as.integer(ylab))},
            'lt' = {dplyr::filter(rp, as.integer(xlab)<as.integer(ylab))},
            'all' = rp
    )
  
  return(rp)
}

leave_one_out_omega <- function(alpha, fun){
  S=nrow(alpha)
  if(S<3){return(NA)}
  out = rep(NA, S)
  for(i in 1:S){
    out[i]=fun(alpha[-i,-i])
  }
  return(out)
}



# Calculate L-p norm of a vector
vector_norm = function(x,p=2){
  (sum((abs(x))^p))^(1/p)
}

# Normalize vector to make its L-p norm = 1
normalize = function(x, p=2){
  x / vector_norm(x,p=p)
}

# euclidian_d = function(v1,v2){
#   vector_norm(v1-v2,p=2)
# }

fitness_landscape <- function(range, resident_x, log = F){
  
  ls_x = seq(range[1],range[2],length.out = 1000)
  ls_y = rep(NA,length(ls_x))
  
  x_end = resident_x
  N_end = Nstar(alpha = alpha_mat(mu_v = x_end,alphafun = alphafun, aa = list()),
                r = sapply(x_end, rfun))
  
  for(i in 1:length(ls_x)){
    alpha_mr = mapply(alphafun, rep(ls_x[i], length(x_end)), x_end)
    ls_y[i] = rfun(ls_x[i]) - sum(alpha_mr * N_end)
  }
  
  if(log){
    plot(ls_x,-log(abs(ls_y)))
  }else{
    plot(ls_x,ls_y)
  }
}

log_seq= function(from, to,length.out){
  f_ = log(from)
  t_ = log(to)
  out_ = seq(f_,t_,length.out = length.out)
  out = exp(out_)
  return(out)
}


# log_abs<-function(x){
#   
#   sign(x)*log(abs(x))
# }
# log_abs_inv<-function(x){
#   sqrt(x)
# }

# Calculates metrics for a given dataframe 


get_metrics = function(df, short = TRUE){
  
  # alphafun = get(settings$alpha)
  # rfun = get(settings$r)
  
  df$S = sapply(df$x, length)

  df$Nstars <- mapply(FUN = Nstar, alpha = df$alpha, r = df$rs)
  df$sum_n_stars <- sapply(df$Nstars, sum)
  df$omega <- sapply(FUN = tb$Omega.L1, X = df$alpha) # log10
  df$omega_rel <- 10^(df$omega/ sapply(df$rs, FUN = length)) # linear
  df$theta <-  mapply(tb$theta.L1, alpha=df$alpha, r=df$rs)
  
  df$alphaii <- lapply(df$alpha, diag)
  df$Ks <- mapply(FUN="/", df$rs, df$alphaii)
  df$sum_K <- sapply(df$Ks,sum)
  df$Y <- mapply(FUN="/",df$Nstars,df$Ks)
  df$sum_Y <- sapply(df$Y, sum)
  
  df$bd_effect <- purrr::map2_dbl(.x = df$Nstars, .y = df$Ks, .f = ~sum(.x) - sum(.y*(.x/sum(.x))))
  
  df$rel_bio <- purrr::map2_dbl(.x = df$Nstars, .y = df$Ks, .f = ~sum(.x / mean(.y)))
  
  df$n_evenness <- sapply(df$Nstars, evenness)
  
  df$k_evenness <- sapply(df$Ks, evenness)
  
  df$dist <- mapply(tb16$min_distance, alpha=df$alpha, K=df$rs)
  df$min_dist <- sapply(df$dist, min)
  
  df$eta = mapply(function(x,y){min(tb22$angle_to_border(x,y)$eta)}, x=df$rs, y=df$alpha)

  
  if(short){return(df)}
  
  # df$S = sapply(df$x, length)
  
  df$mean_nn <- sapply(df$x, avg_nn)
  df$alphaii <- lapply(df$alphas, diag)
  df$mean_alphaii <- sapply(df$alphaii, FUN = mean)
  df$alpha_tilde <- lapply(df$alphas, FUN = alpha_tilde) # is this right?
  df$mean_alpha_tilde <- sapply(df$alpha_tilde, FUN = mean)
  
  # sim$rhos <- lapply(sim$x, sapply, rhofun, pr) # watch out for positional matching. else use r_vec() or an anonymous function
  
  df$mean_r <- sapply(df$rs,mean)
  df$alpha_bar <- mapply(FUN = alpha_bar, alphamat = df$alphas, K = df$Ks)
  df$mean_alpha_bar <- sapply(df$alpha_bar, FUN = mean)
  df$mean_K <- sapply(df$Ks,mean)
  
  return(df)
}


get_metrics_r_alpha = function(df,settings,short = TRUE){
  
  # alphafun = get(settings$alpha)
  # rfun = get(settings$r)
  
  df$S = sapply(df$x, length)
  df$alphas <- lapply(df$x, FUN = alpha_mat, alphafun = alphafun, aa = list())
  df$rs <-  lapply(df$x, sapply, rfun)
  df$Nstars <- mapply(FUN = Nstar, alpha = df$alphas, r = df$rs)
  df$sum_n_stars <- sapply(df$Nstars, sum)
  df$omega <- sapply(FUN = tb$Omega.L1, X = df$alphas) # log10
  df$omega_rel <- 10^(df$omega/ sapply(df$rs, FUN = length)) # linear
  df$theta <-  mapply(tb$theta.L1, alpha=df$alphas, r=df$rs)
  
  if(short){return(df)}
  
  # df$S = sapply(df$x, length)
  
  df$mean_nn <- sapply(df$x, avg_nn)
  df$alphaii <- lapply(df$alphas, diag)
  df$mean_alphaii <- sapply(df$alphaii, FUN = mean)
  df$alpha_tilde <- lapply(df$alphas, FUN = alpha_tilde) # is this right?
  df$mean_alpha_tilde <- sapply(df$alpha_tilde, FUN = mean)
  
  # sim$rhos <- lapply(sim$x, sapply, rhofun, pr) # watch out for positional matching. else use r_vec() or an anonymous function
  
  df$mean_r <- sapply(df$rs,mean)
  df$Ks <- mapply(FUN="/", df$rs, df$alphaii)
  df$alpha_bar <- mapply(FUN = alpha_bar, alphamat = df$alphas, K = df$Ks)
  df$mean_alpha_bar <- sapply(df$alpha_bar, FUN = mean)
  df$mean_K <- sapply(df$Ks,mean)
  
  return(df)
}

evenness = function(y){
  
  if(length(y)<2)return(NA)
  p = y/sum(y)
  out = -sum(p * log(p))/log(length(y))
  return(out)
}

