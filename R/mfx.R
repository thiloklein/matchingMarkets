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
#' @param m EITHER: an object returned from function \code{smm}; OR: the file path to a \code{.mat} file after running \code{smm.m} in Matlab.
#' @param toLatex logical: if \code{TRUE} the result tables are printed in Latex format. The default setting is \code{FALSE}.
#' @export
#' @author Thilo Klein 
#' @keywords summary
#' @examples
#' ## Read model results from Klein (2014), Table 5
#' filepath <- system.file("scripts/TK_gibbsiter_ntu.mat", package="matchingMarkets")
#' 
#' ## Apply mfx function and print results
#' M <- mfx(m=filepath)
#' M$mfx.selection; M$mfx.outcome
mfx <- function(m,toLatex=FALSE){

  ## read results from Matlab
  if(is.character(m)){
    m <- readMat(m)
  }
  
  if(!is.null(m$eta) && !is.nan(m$eta[1])){
    ## variable names
    names(m$postmean) <- names(m$poststd) <- c(unlist(m$an), unlist(m$bn), "delta")
    
    ## model matrix
    X <- do.call(rbind.data.frame, m$X)
    names(X) <- unlist(m$bn)
    eta <- c(m$eta, rep(0, length(m$X)-length(m$W)))
    X <- as.matrix(cbind(X,eta))
    
    ## valuation equation
    nrowX <- sum(sapply(1:length(m$W), function(i) nrow(m$W[[i]])))
    sel <- mfxVal(postmean=m$postmean[1:length(unlist(m$an))]
      , poststd=m$poststd[1:length(unlist(m$an))], nrowX=nrowX, toLatex=toLatex)    
    ## structral model outcome
    out <- mfxOut(sims=10000,postmean=m$postmean[-(1:length(unlist(m$an)))]
      , poststd=m$poststd[-(1:length(unlist(m$an)))], X=X, toLatex=toLatex)

    return(list(X=data.frame(X),R=unlist(m$R),alpha=m$postmean[1:length(unlist(m$an))]
      ,alpha.se=m$poststd[1:length(unlist(m$an))],beta=m$postmean[-(1:length(unlist(m$an)))]
      ,beta.se=m$poststd[-(1:length(unlist(m$an)))],mfx.selection=sel,mfx.outcome=out))

  } else{
    ## variable names
    m$postmean <- na.omit(m$postmean)
    m$poststd <- na.omit(m$poststd)
    names(m$postmean) <- names(m$poststd) <- unlist(m$bn)
    
    ## model matrix
    X <- do.call(rbind.data.frame, m$X)
    names(X) <- unlist(m$bn)
    X <- as.matrix(X)
    
    ## model outcome    
    out <- mfxOut(sims=10000,postmean=m$postmean[1:length(m$postmean)]
      , poststd=m$poststd[1:length(m$poststd)], X=X, toLatex=toLatex)
    
    return(list(X=data.frame(X),R=unlist(m$R),beta=c(m$postmean),beta.se=c(m$poststd),mfx.outcome=out))
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
