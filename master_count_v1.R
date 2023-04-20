
setwd("~/My Drive/Statistics/change_points/count")

library(MASS)
library(prisca)

source("functions_count.R")


#This doe not work bB <- c(2,2,2) ,bC <- rep(-1,p)
T <- 150
p <- 1
nA <- 40
nB <- 70
theta <- 1/4
#yA <- rpois(nA, 1)
#yB <- rpois(nB, 10)
#yC <- rpois(T-nA-nB, 5)
yA <- rDPO(nA, mu = 1, sigma = 1/theta, max.value = 10000)
yB <- rDPO(nB, mu = 3, sigma = 1/theta, max.value = 10000)
yC <- rDPO(T-nA-nB, mu = 5, sigma = 1/theta, max.value = 10000)
y <- c(yA,yB,yC)
y <- matrix(y,ncol=p)
plot(y)
#y <- cbind(y,y)

var(yA)

a0 <- 0.0000001
b0 <- a0
#Initialize
L <- 2



########### POISSON ####################
out <- poisson_prod(y,L,a0,b0,ite=1000)
cs <- credible_set(out$alpha.mat,p=0.9)
cp.poi <- sort(unlist(sapply(cs,function(y) y$point)))
cp.poi

#matplot(t(alpha.mat),type="l")
#abline(v=41,col="red")
#abline(v=111,col="red")


########### DOUBLE POISSON ####################
#theta <- 1/2
#theta <- 1/theta
out <- doublepoisson_prod(y,theta,L,a0,b0,ite=1000,theta_infer = TRUE)
cs <- credible_set(out$alpha.mat,p=0.9)
cp.dopoi <- sort(unlist(sapply(cs,function(y) y$point)))
cp.dopoi



#matplot(t(out$alpha.mat),type="l")
#abline(v=41,col="red")
#abline(v=111,col="red")




#z <- sqrt(y+3/8) #anscombe transform
z <- 2*sqrt(y+1/4)

#WBS
library(wbs)
w <- wbs(z)
w.cpt <-  changepoints(w,penalty="bic.penalty")
cp.wbs <- sort(w.cpt$cpt.ic$bic.penalty)+1


#susie
library(susieR)
fitted <- susie(X, z,L = 2)
#fitted$sets
cp.wbs
cp.poi
cp.dopoi
out$theta

#PELT
library(changepoint)
cpt.mean(z,penalty="MBIC",method="PELT",test.stat="Normal",
         class=TRUE,param.estimates=TRUE)@cpts

 
##Inspect
library("InspectChangepoint")
library(RSpectra)
threshold <- compute.threshold(T,p)
alt <- inspect(t(y),threshold = threshold)
alt
plot(alt)

##Example inspect
# n <- 500; p <- 100; ks <- 30; zs <- c(125,250,375)
# varthetas <- c(0.2,0.4,0.6); overlap <- 0.5
# obj <- multi.change(n, p, ks, zs, varthetas, overlap)
# x <- obj$x
# threshold <- compute.threshold(n,p)
# ret <- inspect(x, threshold = threshold)
# ret
# summary(ret)
# plot(ret)



###What susie can get (1d attempt)
library(susieR)
covmatrix_simple <- function(T){
  
  X <- matrix(0,nrow=T,ncol=T)
  #bracketing
  for (j in 1:T){
    X[j,1:j] <- 1
  }
  
  return(X=X)
}
X <- covmatrix_simple(T)
fitted <- susie(X, z,L = 2)
susie_plot(fitted, y="PIP")
abline(v=41,col="red")
abline(v=111,col="red")
