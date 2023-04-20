
"Functions multivariate change point detection"

revcumsum <- function(x){
  return(rev(cumsum(rev(x))))
}



nullmodel <- function(y,invsigma){
  "Compute recursive components A, corresponding to the Gaussian kernel of the null model"
  
  n <- dim(y)[1]
  save <- c()
  for (i in 1:n){
    save <- c(save,-t(y[i,] %*% invsigma %*% y[i,]))
  }
  return(A=cumsum(save))
}

singlecp_poisson_1d <- function(y,exp_lam){
  
  
  "Univariate Poisson single change in mean"
  "NB: 
  -it assumes pi=1/T
  - the posterior expectation is conditioanlly on gammat=1"
  par_sum <- revcumsum(y)

  
}


singlecp_poisson_d <- function(y,exp_lamCA,b0,log_a_post,a_post){
  
  "Multivariate  DoublePoisson single change in mean"
  "NB: 
  -it assumes pi=1/T
  - the posterior expectation is conditioanlly on gammat=1"
  
  T <- dim(y)[1]
  p <- dim(y)[2]
  #A_all <- nullmodel(y,invsigma)
  savelogPt <- c()
  cumsum_lam <- apply(exp_lamCA,2,cumsum)
  revcumsum_lam <- apply(exp_lamCA,2,revcumsum)
  b_post <- b0 + revcumsum_lam #Posterior b for the gamma
  logPt <- -cumsum_lam + log_a_post - a_post * log(b_post)
  logPt <- rowSums(logPt)
  nc <- matrixStats::logSumExp(logPt)
  alpha <- exp(logPt -nc)
  exp_lam_l = apply(a_post/b_post*alpha,2,cumsum)+(1-cumsum(alpha))
  return(list(alpha=alpha,b_post=b_post,exp_lam_l=exp_lam_l))#Wrong exp
}

poisson_prod <- function(y,L,a0=0.001,b0=0.001,ite=1000){
  
  "Multiple changes in Poisson"
  
    alpha.mat <- matrix(rep(0,T*L),ncol=T)
    par_sum <- apply(y,2,revcumsum)
    a_post <- a0+par_sum
    log_a_post <- lgamma(a_post)
    exp_lam<- list()
    for (l in 1:L){exp_lam[[l]] <- matrix(rep(1,T*p),nrow=T)}
    #Start the coordinate ascent
    for (i in 1:ite){
      for (l in 1:L){
        exp_lamCA <- Reduce("*",exp_lam[-l])
        single <- singlecp_poisson_d(y,exp_lamCA,b0,log_a_post,a_post)
        alpha.mat[l,] <- single$alpha
        exp_lam[[l]] <- single$exp_lam_l #Wrong
      }
    }
  return(list(alpha.mat=alpha.mat))
}


singlecp_doublepoisson_d <- function(y,th_exp_lamCA,b0,log_a_post,a_post,theta){
  
  "Multivariate Double Poisson single change in mean"
  "NB: 
  -it assumes pi=1/T
  - the posterior expectation is conditioanlly on gammat=1"
  
  T <- dim(y)[1]
  p <- dim(y)[2]
  savelogPt <- c()
  cumsum_th_lam <- apply(th_exp_lamCA,2,cumsum)
  revcumsum_th_lam <- apply(th_exp_lamCA,2,revcumsum)
  b_post <- b0 + revcumsum_th_lam #Posterior b for the gamma
  logPt <- -cumsum_th_lam + log_a_post - a_post * log(b_post) #- revcumsum_th_ylogof
  #maxlogPt <- max(logPt)
  logPt <- rowSums(logPt)
  nc <- matrixStats::logSumExp(logPt)#-maxlogPt)
  alpha <- exp(logPt -nc)#-maxlogPt)
  exp_lam_l = apply(a_post/b_post*alpha,2,cumsum)+(1-cumsum(alpha))
  return(list(alpha=alpha,b_post=b_post,exp_lam_l=exp_lam_l))#Wrong exp
}




doublepoisson_prod <- function(y,theta=1,L,a0=0.001,b0=0.001,ite=1000,theta_infer=FALSE){
  
  "Multiple changes in Double Poisson"
  
  if (theta_infer==TRUE) theta <- 1
  alpha.mat <- matrix(rep(0,T*L),ncol=T)
  par_sum <- apply(y,2,revcumsum)
  a_post <- a0+theta*par_sum
  log_a_post <- lgamma(a_post)
  exp_lam<- list()
  for (l in 1:L){exp_lam[[l]] <- matrix(rep(1,T*p),nrow=T)}
  #Start the coordinate ascent
  for (i in 1:ite){
    for (l in 1:L){
      th_exp_lamCA <- theta*Reduce("*",exp_lam[-l])
      single <- singlecp_doublepoisson_d(y,th_exp_lamCA,b0,log_a_post,a_post,theta)
      alpha.mat[l,] <- single$alpha
      exp_lam[[l]] <- single$exp_lam_l 
    }
    if (theta_infer==TRUE){
      logofy <- ifelse(y == 0, 1, log(y))
      lam <- Reduce("*",exp_lam)
      theta <- 1/(2*(mean(y*logofy)-mean(y*log(lam))+mean(lam)-mean(y)))
      a_post <- a0+theta*par_sum
      log_a_post <- lgamma(a_post)
    }
  }
  return(list(alpha.mat=alpha.mat,theta=theta))
}








