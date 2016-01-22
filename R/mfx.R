# ----------------------------------------------------------------------------
# R-code (www.r-project.org/) to obtain marginal effects for probit and matching models
#
# Copyright (c) 2013 Thilo Klein
#
# This library is distributed under the terms of the GNU Public License (GPL)
# for full details see the file LICENSE
#
# ----------------------------------------------------------------------------

#' @title Marginal effects for probit and matching models
#'
#' @description Marginal effects from regression coefficients for probit 
#' and matching models. 
#'
#' @param m an object returned from functions \code{stabit} or \code{stabit2}.
#' @param toLatex logical: if \code{TRUE} the result tables are printed in Latex format. The default setting is \code{FALSE}.
#' 
#' @export
#' 
#' @import stats
#' 
#' @author Thilo Klein 
#' 
#' @keywords summary
#' 
#' @references Klein, T. (2015a). \href{https://ideas.repec.org/p/cam/camdae/1521.html}{Does Anti-Diversification Pay? A One-Sided Matching Model of Microcredit}.
#' \emph{Cambridge Working Papers in Economics}, #1521.
#' 
#' @examples
#' ## 1. load results from Klein (2015a)
#'  data(klein15a)
#' 
#' ## 2. apply mfx function and print results
#'  mfx(m=klein15a)
mfx <- function(m,toLatex=FALSE){
  
  if(!is.null(m$coefs$alpha)){ ## Selectiom and Outcome Eqns

    ## model matrix
    X <- do.call(rbind.data.frame, m$model.list$X)
    eta <- c(m$coefs$eta, rep(0, length(m$model.list$X)-length(m$model.list$W)))
    X <- as.matrix(data.frame(X=X,eta=eta))
    
    ## valuation equation
    nrowX <- dim(do.call(rbind.data.frame, m$model.list$W))[1]
    sel <- mfxVal(postmean=m$coefs$alpha[,1], poststd=m$coefs$alpha[,2],
           nrowX=nrowX, toLatex=toLatex)
    
    ## structral model outcome
    out <- mfxOut(sims=10000, postmean=unlist(c(m$coefs$beta[,1], data.frame(delta=m$coefs$delta[1]))),
           poststd=c(m$coefs$beta[,2], m$coefs$delta[2]), X=X, toLatex=toLatex)

    return(list(mfx.selection=sel, mfx.outcome=out))

  } else{ ## Outcome Eqn only
    
    ## model matrix
    X <- do.call(rbind.data.frame, m$model.list$X)
    X <- as.matrix(X)
    
    ## model outcome    
    out <- mfxOut(sims=10000, postmean=m$coefs$beta[,1],
           poststd=m$coefs$beta[,2], X=X, toLatex=toLatex)
    
    return(list(mfx.outcome=out))
  }
}


mfxOut <- function(sims=10000,x.mean=TRUE,postmean,poststd,X,toLatex){
  ## source: http://researchrepository.ucd.ie/handle/10197/3404
  ## method: average of individual marginal effects at each observation
  ## interpretation: http://www.indiana.edu/~statmath/stat/all/cdvm/cdvm.pdf page 8
  set.seed(1984)
  if(x.mean==TRUE){
    ## marginal effects are calculated at the means of independent variables
    pdf <- dnorm(mean(X%*%postmean))
    pdfsd <- dnorm(sd(X%*%postmean))
  } else{
    ## marginal effects are calculated for each observation and then averaged
    pdf <- mean(dnorm(X%*%postmean))
    pdfsd <- sd(dnorm(X%*%postmean))
  }  
  mx <- pdf*postmean

  sim <- matrix(rep(NA,sims*length(postmean)), nrow=sims)
  for(i in 1:length(postmean)){
    sim[,i] <- rnorm(sims,postmean[i],poststd[i])
  }
  pdfsim <- rnorm(sims,pdf,pdfsd)
  sim.se <- pdfsim*sim
  s.e. <- apply(sim.se,2,sd)

  t.stat <- mx/s.e.
  p.val <- pt(-abs(t.stat),df=dim(X)[1]-length(postmean)+1)
  stars <- ifelse(p.val<0.001,"***",ifelse(p.val<0.01,"**",ifelse(p.val<0.05,"*",ifelse(p.val<0.10,".",""))))
  if(toLatex==FALSE){
    res <- data.frame(round(cbind(mx, s.e., t.stat, p.val),3), stars)
  } else{
      sign <- ifelse(mx>0,"~","")
      res <- data.frame("&", sign, round(mx,3), se=paste(paste("(",round(s.e.,3),sep=""),")",sep=""), stars, "\\")
  }
  return(res)
}


mfxVal <- function(postmean,poststd,nrowX,toLatex){

  ## Reference: Sorensen (2007, p. 2748)

  mx <- dnorm(0)*postmean/sqrt(2)
  s.e. <- dnorm(0)*poststd/sqrt(2)
  t.stat <- mx/s.e.
  p.val <- pt(-abs(t.stat),df=nrowX-length(postmean))
  stars <- ifelse(p.val<0.001,"***",ifelse(p.val<0.01,"**",ifelse(p.val<0.05,"*",ifelse(p.val<0.10,".",""))))
  if(toLatex==FALSE){
    res <- data.frame( round(cbind(mx, s.e., t.stat, p.val),3), stars)
  } else{
    sign <- ifelse(mx>0,"~","")
    res <- data.frame( "&", sign, round(mx,3), se=paste(paste("(",round(s.e.,3),sep=""),")",sep=""), stars, "\\")
  }
  return(res)
}
