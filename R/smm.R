# ----------------------------------------------------------------------------
# R-code (www.r-project.org/) for the Structural Matching Model
#
# Copyright (c) 2013 Thilo Klein
#
# This library is distributed under the terms of the GNU Public License (GPL)
# for full details see the file LICENSE
#
# ----------------------------------------------------------------------------

#' @title Structural Matching Model
#'
#' @description Corrects for sample selection when the selection process is a one-sided matching game.
#'
#' @param selection logical: if TRUE estimate structural model with selection and outcome equation; if FALSE use outcome equation only.
#' @param NTU logical: if TRUE use non-transferable utility (NTU) matching game; if FALSE use transferable utility (TU) matching game.
#' @param binary logical: if TRUE outcome variable is taken to be binary; if FALSE outcome variable is taken to be continuous.
#' @param offsetOut vector of integers indicating the indices of columns in X for which coefficients should be forced to 1. Use 0 for none.
#' @param offsetSel vector of integers indicating the indices of columns in W for which coefficients should be forced to 1. Use 0 for none.
#' @param marketFE logical: if TRUE use market-level fixed effects in outcome equation; if FALSE don't.
#' @param censored delta is 0:not censored, 1:from below, 2:from above
#' @param dropOnes one-group-markets exluded for estimation
#' @param repeatRun repeated run
#' @param niter number of iterations
#' @param data data. see \code{Details} below
#' @param software either \code{"matlab"} or \code{"octave"}
#' @export
#' @return
#' 'smm' returns a list with the following items.
#' \item{alphadraws}{matrix of dimension}
#' \item{betadraws}{matrix of dimension }
#' \item{deltadraws}{vector of length }
#' \item{postmean}{vector of lenght }
#' \item{poststd}{vector of lenght }
#' \item{tstat}{vector of lenght }
#' @author Thilo Klein \email{thilo@@klein.co.uk}
#' @references Klein, T. (2014). matchingMarkets: An R Package for the Analysis of Stable Matchings.
#' @examples
#' ############
#' ## OCTAVE ##
#' ############
#' 
#' ## 1. Load RcppOctave
#' library(RcppOctave)
#' 
#' ## 2. Run Gibbs sampler
#' mdata <- system.file("scripts/MatchingData.mat", package="matchingMarkets")
#' res <- smm(selection=1,NTU=1,binary=1,offsetOut=0,offsetSel=0,marketFE=0
#'   ,censored=0,dropOnes=0,repeatRun=0,niter=10,data=mdata,software="octave")
#'
#' ## 3. Get results
#' res$postmean
#' plot(res$alphadraws[1,], type="l")
#' 
#' ############
#' ## MATLAB ##
#' ############
#' 
#' ## 1. Load R.matlab
#' library(R.matlab)
#' 
#' ## 2. Start MATLAB
#' Matlab$startServer()
#' 
#' ## 3. Create a MATLAB client object used to communicate with MATLAB
#' matlab <- Matlab()
#' isOpen <- open(matlab) ## Connect to the MATLAB server
#' 
#' ## 4. Run Gibbs sampler
#' mdata <- system.file("scripts/MatchingData.mat", package="matchingMarkets")
#' smm(selection=1,NTU=1,binary=1,offsetOut=0,offsetSel=0,marketFE=0,censored=0,
#'   dropOnes=0,repeatRun=0,niter=10,data=mdata,software="matlab")
#' 
#' ## 5. Get results
#' res <- getVariable(matlab, c("alphadraws","betadraws","deltadraws",
#'   "postmean","poststd","tstat"))
#' res$postmean
#' plot(res$alphadraws[1,], type="l")
#' 
#' ## 6. Close MATLAB server
#' close(matlab)
smm <- function(selection,NTU,binary,offsetOut,offsetSel,marketFE,censored,dropOnes,repeatRun,niter,data,software){

  if(software == "matlab"){

    # 4.6 Create a function (M-file) on the MATLAB server
    mfile <- system.file("scripts/smm.m",package="matchingMarkets")
    code <- readLines(mfile)
    setFunction(matlab, code)

    # 4.7 Use the MATLAB function just created
    param <- paste(selection,NTU,binary,offsetOut,offsetSel,marketFE,censored,dropOnes,repeatRun,niter,sep=",")
    evaluate(matlab, paste("[alphadraws,betadraws,deltadraws,postmean,poststd,tstat] = smm(",param,",'",data,"');", sep=""))
    
  } else{

    ## clear all Octave session
    o_clear(all=TRUE)

    ## source Octave file
    mfile <- system.file("scripts/smm.m",package="matchingMarkets")
    o_source(mfile)

    ## call Octave
    res <- .CallOctave("smm",selection,NTU,binary,offsetOut,offsetSel,marketFE,censored,dropOnes,repeatRun,niter,data,argout=7)
    #names(res) <- c("alphadraws","betadraws","deltadraws","postmean","poststd","tstat")
    
  }
}
