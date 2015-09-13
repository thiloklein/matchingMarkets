\name{klein15b}
\docType{data}
\alias{klein15b}

\title{Results of Monte Carlo Simulations in Klein (2015b).}

\description{
   Results of Monte Carlo Simulations in Klein (2015b) for 40 two-group markets.
}

\usage{data(klein15b)}

\format{
   This list contains the following elements:
   \describe{
      \item{exp.5.5.ols}{Benchmark study, OLS: coefficient estimates for 40 markets with groups of 5. Data for all 5 group members is observed.}
      \item{exp.5.5.ntu}{Benchmark study, structural model.}
      \item{exp.6.5.ols}{Experiment 1, OLS: coefficient estimates for 40 markets with groups of 6. Only Data for 5 group members is observed.}
      \item{exp.6.5.ntu}{Experiment 1, structural model.}
      \item{exp.6.6.ols}{Experiment 2, OLS: coefficient estimates for 40 markets with groups of 6. Data for all 6 group members is observed but only a random sample of 250 of the 922 counterfactual groups is used in the analysis.}
      \item{exp.6.6.ntu}{Experiment 2, structural model.}
   }
}

\examples{
\dontrun{
#######################################
## Modes of posterior distributions  ##
#######################################

## load data
data(klein15b)

## define function to obtain the mode
mode <- function(x){
  d <- density(x,bw="SJ")
  formatC(round(d$x[which.max(d$y)], 3), format='f', digits=3)
}

## Benchmark study
apply(klein15b$exp.5.5.ntu, 2, mode)
apply(klein15b$exp.5.5.ols, 2, mode)

## Experiment 1
apply(klein15b$exp.6.5.ntu, 2, mode)
apply(klein15b$exp.6.5.ols, 2, mode)

## Experiment 2
apply(klein15b$exp.6.6.ntu, 2, mode)
apply(klein15b$exp.6.6.ols, 2, mode)

#####################################
## Plot of posterior distributions ##
#####################################

## load data
data(klein15b)

par(mfrow=c(3,3))
tpe <- c(rep("Benchmark",2), rep("Experiment 1",2), rep("Experiment 2",2))
par(mar=c(5.1,4.6,0.8,2.1))

for(i in seq(1,length(klein15b)-1,2)){
  ntu <- klein15b[[i]]
  ols <- klein15b[[i+1]]

  ntu <- ntu[,colnames(ntu) %in% c("alpha","beta.wst.ieq","delta")]
  ols <- ols[,colnames(ols) == "beta.wst.ieq"]

  plot(density(ntu[,1]), xlab=expression(hat(alpha)), 
    ylab="density", main="", axes=FALSE, xlim=c(-1,2))
  axis(2,lwd=2,cex.axis=0.8); axis(1,lwd=2,cex.axis=0.8)
  legend("topleft","Struct.",lty=1,bty="n")
  abline(v=1, lty=3)

  plot(density(ntu[,2]), xlab=expression(hat(beta)), 
    ylab="density", main=tpe[i], axes=FALSE)
  axis(2,lwd=2,cex.axis=0.8); axis(1,lwd=2,cex.axis=0.8)
  points(density(ols), type="l", lty=2)
  legend("topright",c("Struct.","OLS"),lty=c(1,2),bty="n")
  abline(v=-1, lty=3)

  plot(density(ntu[,3]), xlab=expression(hat(delta)), 
    ylab="density", main="", axes=FALSE)
  axis(2,lwd=2,cex.axis=0.8); axis(1,lwd=2,cex.axis=0.8)
  legend("topleft","Struct.",lty=1,bty="n")
  abline(v=0.5, lty=3)
}
}}

\references{
Klein, T. (2015a). \href{https://ideas.repec.org/p/cam/camdae/1521.html}{Does Anti-Diversification Pay? A One-Sided Matching Model of Microcredit}. \emph{Cambridge Working Papers in Economics}, #1521.

Klein, T. (2015b). \href{http://cran.r-project.org/package=matchingMarkets}{Analysis of stable matchings in R: Package matchingMarkets}. \emph{Vignette to R package matchingMarkets}, The Comprehensive R Archive Network.
}

\keyword{datasets}
