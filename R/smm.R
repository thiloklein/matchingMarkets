# ----------------------------------------------------------------------------
# R-code (www.r-project.org/) for the Structural Matching Model
#
# Copyright (c) 2013 Thilo Klein
#
# This library is distributed under the terms of the GNU Public License (GPL)
# for full details see the file LICENSE
#
# ----------------------------------------------------------------------------

#' @title Structural Matching Model to correct for sample selection bias
#'
#' @description Please use \code{\link[matchingMarkets:matchingMarkets-package]{help(matchingMarkets)}}
#' to check the \emph{System requirements} in the package description before
#' using this function.
#' 
#' The function provides a Gibbs sampler for a structural matching model that corrects for sample 
#' selection bias when the selection process is a one-sided matching game; that is, 
#' group/coalition formation.
#'
#' The \emph{Selection Equation} determines which matches are observed (\eqn{D=1}) and which 
#' are not (\eqn{D=0}).
#' \deqn{ \begin{array}{lcl}
#'        D &= & 1[V \in \Gamma] \\
#'        V &= & W\alpha + \eta
#'        \end{array}
#'      }{ D = 1[V in \Gamma] with V = W\alpha + \eta
#'      }
#' Here, \eqn{V} is a vector of latent valuations of \emph{all feasible} matches, ie observed and unobserved, and \eqn{1[.]} is the Iverson 
#' bracket. A match is observed if its match valuation is in the set of valuations \eqn{\Gamma}
#' that satisfy the equilibrium condition (see Klein, 2014). This condition differs for matching
#' games with transferable and non-transferable utility and can be specified using the \code{NTU} 
#' argument. The match valuation \eqn{V} is
#' a linear function of \eqn{W}, a matrix of characteristics for \emph{all feasible} groups,
#' and \eqn{\eta}, a vector of random errors. \eqn{\alpha} is a paramter vector to be estimated.
#' 
#' The \emph{Outcome Equation} determines the outcome for \emph{observed} matches. The dependent
#' variable can either be continuous or binary, dependent on the value of the \code{binary}
#' argument. In the binary case, the dependent variable \eqn{R} is determined by a threshold 
#' rule for the latent variable \eqn{Y}.
#' \deqn{ \begin{array}{lcl}
#'        R &= & 1[Y > c] \\
#'        Y &= & X\beta + \epsilon
#'        \end{array}
#'      }{ R = 1[Y > c] with Y = X\beta + \epsilon
#'      }
#' Here, \eqn{Y} is a linear function of \eqn{X}, a matrix of characteristics for \emph{observed} 
#' matches, and \eqn{\epsilon}, a vector of random errors. \eqn{\beta} is a paramter vector to 
#' be estimated.
#' 
#' The structural model imposes a linear relationship between the error terms of both equations 
#' as \eqn{\epsilon = \delta\eta + \xi}, where \eqn{\xi} is a vector of random errors and \eqn{\delta}
#' is the covariance paramter to be estimated. If \eqn{\delta} were zero, the marginal distributions
#' of \eqn{\epsilon} and \eqn{\eta} would be independent and the selection problem would vanish.
#' That is, the observed outcomes would be a random sample from the population of interest.
#' 
#' @param selection logical: if \code{TRUE} estimate structural model with selection and outcome equation; if \code{FALSE} use outcome equation only.
#' @param NTU logical: if \code{TRUE} use non-transferable utility (NTU) group formation game; if \code{FALSE} use transferable utility (TU) group formation game.
#' @param binary logical: if \code{TRUE} outcome variable is taken to be binary; if \code{FALSE} outcome variable is taken to be continuous.
#' @param offsetOut vector of integers indicating the indices of columns in \code{X} for which coefficients should be forced to 1. Use 0 for none.
#' @param offsetSel vector of integers indicating the indices of columns in \code{W} for which coefficients should be forced to 1. Use 0 for none.
#' @param marketFE logical: if \code{TRUE} market-level fixed effects are used in outcome equation; if \code{FALSE} no market fixed effects are used.
#' @param censored draws of the \code{delta} parameter that estimates the covariation between the error terms in selection and outcome equation are 0:not censored, 1:censored from below, 2:censored from above.
#' @param gPrior logical: if \code{TRUE} the g-prior (Zellner, 1986) is used for the variance-covariance matrix.
#' @param dropOnes logical: if \code{TRUE} one-group-markets are exluded from estimation.
#' @param interOut two-colum matrix indicating the indices of columns in \code{X} that should be interacted in estimation. Use 0 for none.
#' @param interSel two-colum matrix indicating the indices of columns in \code{W} that should be interacted in estimation. Use 0 for none.
#' @param repeatRun logical: if \code{TRUE} repeated run. This is useful when .
#' @param niter number of iterations to use for the Gibbs sampler.
#' @param data an object from \code{design.matrix}.
#' @param software which software to use for the Gibbs sampler. EITHER: \code{"matlab"}; OR: \code{"octave"}.
#' @export
#' @return
#' \code{smm} returns a list with the following items.
#' \item{alphadraws}{matrix of dimension \code{ncol(W)} \code{x} \code{niter} comprising all paramter draws for the selection equation.}
#' \item{betadraws}{matrix of dimension \code{ncol(X)} \code{x} \code{niter} comprising all paramter draws for the outcome equation.}
#' \item{deltadraws}{vector of length \code{niter} comprising all draws for the \code{delta} parameter.}
#' \item{X}{list with \code{X[[t]][G,]} containing characteristics of group \code{G} in market \code{t} (equilibrium groups only).}
#' \item{W}{list with \code{W[[t]][G,]} containing characteristics of group \code{G} in market \code{t} (all feasible groups).}
#' \item{R}{list of group-level outcomes.}
#' \item{eta}{vector containing the mean of all \code{eta} draws for each observed group.}
#' \item{postmean}{vector comprising the coefficient estimates.}
#' \item{poststd}{vector comprising the coefficient standard errors.}
#' \item{an}{vector comprising the coefficient names for the selection equation.}
#' \item{bn}{vector comprising the coefficient names for the outcome equation.}
#' \item{sigmasquareximean}{variance estimate of the error term \code{xi} in the outcome equation.}
#' @author Thilo Klein 
#' @keywords regression
#' @import R.matlab 
#' @references Klein, T. (2014). Stable matching in microcredit: Implications for market design & econometric analysis, PhD thesis, \emph{University of Cambridge}.
#' @references Zellner, A. (1986). \emph{On assessing prior distributions and Bayesian regression analysis with g-prior distributions}, volume 6, pages 233--243. North-Holland, Amsterdam.
#' @examples
#' ############
#' ## OCTAVE ##
#' ############
#' 
#' \dontrun{
#' ## 1. Load data and run Gibbs sampler
#' mdata <- system.file("scripts/matchingDataNTU.mat", package="matchingMarkets")
#' fit <- smm(selection=TRUE, NTU=TRUE, binary=TRUE, gPrior=TRUE, niter=200, 
#'        data=mdata, software="octave")
#'
#' ## 2. Get results
#' fit$postmean
#' plot(fit$alphadraws[1,], type="l")}
#' 
#' ############
#' ## MATLAB ##
#' ############
#' 
#' \dontrun{
#' ## 1. Load data and run Gibbs sampler
#' mdata <- system.file("scripts/matchingDataNTU.mat", package="matchingMarkets")
#' fit <- smm(selection=TRUE, NTU=TRUE, binary=TRUE, gPrior=TRUE, niter=200, 
#'        data=mdata, software="matlab")
#' 
#' ## 2. Get results
#' fit$postmean
#' plot(fit$alphadraws[1,], type="l")}
smm <- function(selection, NTU=TRUE, binary, offsetOut=FALSE, offsetSel=FALSE, 
  marketFE=FALSE, censored=FALSE, gPrior=FALSE, dropOnes=FALSE, interOut=FALSE, 
  interSel=FALSE, repeatRun=FALSE, niter, data, software){

  ## If 'data' is not a file path, then safe data.mat to temporary directory and obtain path
  if(!is.character(data)){
    fName <- paste(tempdir(),"/",sep="")
    toMatLab(x=data, folder=fName)
    data <- paste(fName,"data.mat",sep="")
  }

  ## Convert logical statements TRUE/FALSE to 1/0 for compatibility with MATLAB
  arglist <- mget(names(formals()),sys.frame(sys.nframe()))
  arglist$data <- data
  arglist <- lapply(arglist, function(i) ifelse(i==TRUE,1, ifelse(i==FALSE,0,i)))
  for(i in 1:length(arglist)){
    assign(names(arglist)[i], arglist[[i]])
  } 
  
  if(software == "matlab"){
    
    ## 1. Load R.matlab and start MATLAB
    Matlab$startServer()

    ## 2. Create a MATLAB client object used to communicate with MATLAB
    matlab <- Matlab()
    isOpen <- open(matlab) ## Connect to the MATLAB server

    ## 3. Create a function (M-file) on the MATLAB server
    mfile <- system.file("scripts/smm.m", package="matchingMarkets")
    code <- readLines(mfile)
    matlab <- get("matlab")
    setFunction(matlab, code)
    #setOption(matlab, "readResult/interval", 10)
    #setOption(matlab, "readResult/maxTries", 30*(60/10))

    ## 4. Use the MATLAB function just created
    inputs <- paste(selection, NTU, binary, offsetOut, offsetSel, marketFE, censored,
      gPrior, dropOnes, interOut, interSel, repeatRun, niter, sep=",")
    outputs <- "alphadraws,betadraws,deltadraws,X,W,R,eta,postmean,poststd,an,bn,sigmasquareximean"
    evaluate(matlab, paste("[",outputs,"]=smm(",inputs,",'",data,"');", sep=""))
    
    ## 5. Get results
    fit <- getVariable(matlab, strsplit(outputs, ",")[[1]])
    
    ## 6. Close MATLAB server and return results
    close(matlab)
    fit
    
  } else{

    if(require(RcppOctave)){
      ## clear all Octave session
      o_clear(all=TRUE)

      ## source Octave file
      mfile <- system.file("scripts/smm.m", package="matchingMarkets")
      o_source(mfile)

      ## call Octave
      fit <- .CallOctave("smm", selection, NTU, binary, offsetOut, offsetSel, marketFE, censored, 
        gPrior, dropOnes, interOut, interSel, repeatRun, niter, data, argout=12)      
    }    
    
  }
}
