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

#par(mfrow=c(3,3))
tpe <- c(rep("Benchmark",2), rep("Experiment 1",2), rep("Experiment 2",2))
#bottom, left, top, right
#par(mar=c(5.1,4.6,2.8,2.1))

for(i in seq(1,length(klein15b)-1,2)){
  ntu <- klein15b[[i]]
  ols <- klein15b[[i+1]]

  ntu <- ntu[,colnames(ntu) == "beta.wst.ieq"]
  ols <- ols[,colnames(ols) == "beta.wst.ieq"]
  
  if(i == 1){
    draws <- data.frame(Structural=ntu, OLS=ols, type=tpe[i]) #, stringsAsFactors=FALSE
  } else{
    draws <- rbind(draws, data.frame(Structural=ntu, OLS=ols, type=tpe[i]))
  }
}

library(lattice)
lattice.options(default.theme = standard.theme(color = FALSE))
keys <- list(text=c("Structural model","OLS"), space="top", columns=2, lines=TRUE)
densityplot( ~ Structural + OLS | type, plot.points=FALSE, auto.key=keys,
       data = draws, xlab = "coefficient draws", ylab = "density", type = "l",
       panel = function(x,...) {
         panel.densityplot(x,...)
         panel.abline(v=-1, lty=3)
       })




## ---- mf-gibbs-draws ----

load("~/Documents/Research/Microfinance/2_cleanData/klein15a.RData")

dd <- cbind(#beta.pi   = klein15a$draws$betadraws["pi.inv",seq(1,800000,1000)],
            beta.wst  = klein15a$draws$betadraws["wst.ieq",seq(1,800000,1000)],
            delta     = klein15a$draws$deltadraws[1,seq(1,800000,1000)],
            #alpha.pi  = klein15a$draws$alphadraws["pi.inv",seq(1,800000,1000)],
            alpha.wst = klein15a$draws$alphadraws["wst.ieq",seq(1,800000,1000)],
            iteration = 1:800)

library(tidyr)
dd.long <- gather(as.data.frame(dd), condition, measurement, beta.wst:alpha.wst)
dd.long$condition <- as.factor(dd.long$condition)

library(lattice)
lattice.options(default.theme = standard.theme(color = FALSE))
xyplot(measurement ~ iteration | condition, 
       data = dd.long, scales=list(relation="free"),
       xlab = "iterations in thousands",
       ylab = "coefficient draws", type = "l") #index.cond = list(c(3,4,5,1,2))




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






