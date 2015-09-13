## ---- monte-carlo-experiments ----

## 1. Set parameters
mciter <- 2 #500
niter <- 10 #400000
nodes <- 4

## 2. Setup parallel backend to use 4 processors
library(foreach); library(doSNOW)
cl <- makeCluster(nodes); registerDoSNOW(cl)

## 3. Define foreach loop function
mce.add <- function(mciter, niter, N, n, type, method){
  h <- foreach(i=1:mciter) %dopar% {
    library(matchingMarkets)
    mce(seed=i,niter, N, n, type, method)
  }
  do.call(rbind, h)
}

## 4. Run siumlations:

## 4-a. Benchmark study
exp.5.5.ols <- mce.add(mciter=mciter, niter=niter, N=5, n=5,
                       type="group.members", method="outcome")
exp.5.5.ntu <- mce.add(mciter=mciter, niter=niter, N=5, n=5,
                       type="group.members", method="NTU")

## 4-b. Experiment 1: randomly sampled group members
exp.6.5.ols <- mce.add(mciter=mciter, niter=niter, N=6, n=5,
                       type="group.members", method="outcome")
exp.6.5.ntu <- mce.add(mciter=mciter, niter=niter, N=6, n=5,
                       type="group.members", method="NTU")

## 4-c. Experiment 2: randomly sampled counterfactual groups
exp.6.6.ols <- mce.add(mciter=mciter, niter=niter, N=6, n=6,
                       type="counterfactual.groups", method="outcome")
exp.6.6.ntu <- mce.add(mciter=mciter, niter=niter, N=6, n=6,
                       type="counterfactual.groups", method="NTU")

## 5. Stop parallel backend
stopCluster(cl)




## ---- mc-text ----
## load data
library("matchingMarkets")
data(klein15b)

sd.b.ntu  <- round( sd(data.frame(klein15b$exp.5.5.ntu)$beta.wst.ieq), 2)
sd.e2.ntu <- round( sd(data.frame(klein15b$exp.6.6.ntu)$beta.wst.ieq), 2)




## ---- mc-table ----
## load data
library("matchingMarkets")
data(klein15b)

## define function to obtain the mode
mode <- function(x){
d <- density(x,bw="SJ")
formatC( round( d$x[which.max(d$y)], 3), format='f', digits=3)
}

## Benchmark study
b.ntu <- apply(klein15b$exp.5.5.ntu, 2, mode)
b.ols <- apply(klein15b$exp.5.5.ols, 2, mode)

## Experiment 1
e1.ntu <- apply(klein15b$exp.6.5.ntu, 2, mode)
e1.ols <- apply(klein15b$exp.6.5.ols, 2, mode)

## Experiment 2
e2.ntu <- apply(klein15b$exp.6.6.ntu, 2, mode)
e2.ols <- apply(klein15b$exp.6.6.ols, 2, mode)




## ---- mc-plots ----
library("matchingMarkets")
data(klein15b)

par(mfrow=c(3,3))
tpe <- c(rep("Benchmark",2), rep("Experiment 1",2), rep("Experiment 2",2))
par(mar=c(5.1,4.6,0.8,2.1))

for(i in seq(1,length(klein15b)-1,2)){
  ntu <- klein15b[[i]]
  ols <- klein15b[[i+1]]

  ntu <- ntu[,colnames(ntu) %in% c("alpha","beta.wst.ieq","delta")]
  ols <- ols[,colnames(ols) == "beta.wst.ieq"]

  plot(density(ntu[,1]), xlab=expression(hat(alpha)), ylab="density", main="", axes=FALSE, xlim=c(-1,2))
  axis(2,lwd=2,cex.axis=0.8); axis(1,lwd=2,cex.axis=0.8)
  legend("topleft","Struct.",lty=NULL,bty="n")
  #legend("topleft","Struct.",lty=1,bty="n")
  abline(v=1, lty=3)

  plot(density(ntu[,2]), xlab=expression(hat(beta)), ylab="density", main=tpe[i], axes=FALSE)
  axis(2,lwd=2,cex.axis=0.8); axis(1,lwd=2,cex.axis=0.8)
  points(density(ols), type="l", lty=2)
  legend("topleft","Struct.",lty=NULL,bty="n")
  legend("topright","OLS",lty=NULL,bty="n")
  #legend("topright",c("Struct.","OLS"),lty=c(1,2),bty="n")
  abline(v=-1, lty=3)

  plot(density(ntu[,3]), xlab=expression(hat(delta)), ylab="density", main="", axes=FALSE)
  axis(2,lwd=2,cex.axis=0.8); axis(1,lwd=2,cex.axis=0.8)
  legend("topleft","Struct.",lty=NULL,bty="n")
  #legend("topleft","Struct.",lty=1,bty="n")
  abline(v=0.5, lty=3)
}




## ---- figure-measurement-error-panel-1 ----

A <- seq(0.5,1,0.01)
B <- 1-A

random <- A^2 + B^2
matching <- (1 + (2*(A-.5))^2 + (2*B)^2)/2

par(mar=c(5.1,4.6,0.8,2.1))
par(lwd=2, cex.axis=1, cex=1.8)

plot(random ~ A,type="l",axes=F,xlab=expression(theta[A]),ylab=expression(list(X[t],tilde(X)[t])))
axis(2); axis(1)
polygon(c(A, rev(A)), c(matching, rev(random)),
     col = "grey90", border = NA)
grid()
points(matching ~ A,type="l",lty=2)
points(random ~ A,type="l")
legend("bottomright", c("assortative","random"), lty=c(2,1), bty="n")
#text(x=.65,y=.7,"measurement error")
text(x=.65,y=.7,"measurement")
text(x=.65,y=.65,"error")


## ---- figure-measurement-error-panel-2 ----

par(mar=c(5.1,4.6,0.8,2.1))
par(lwd=2, cex.axis=1, cex=1.8)

plot(matching-random ~ random,type="l",axes=F,lwd=2,xlab=expression(tilde(X)[t]),ylab=expression(X[t]-tilde(X)[t]))
axis(2); axis(1)
grid()






