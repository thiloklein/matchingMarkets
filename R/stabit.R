# ----------------------------------------------------------------------------
# R-code (www.r-project.org/) for the Structural Matching Model
#
# Copyright (c) 2013 Thilo Klein
#
# This library is distributed under the terms of the GNU Public License (GPL)
# for full details see the file LICENSE
#
# ----------------------------------------------------------------------------

#' @title Matching model and selection correction for group formation
#'
#' @description The function provides a Gibbs sampler for a structural matching model that 
#' estimates preferences and corrects for sample selection bias when the selection process 
#' is a one-sided matching game; that is, group/coalition formation.
#'
#' The input is individual-level data of all group members from one-sided matching marktes; that is, 
#' from group/coalition formation games. 
#' 
#' In a first step, the function generates a model matrix with characteristics of \emph{all feasible} 
#' groups of the same size as the observed groups in the market. 
#' 
#' For example, in the stable roommates problem with \eqn{n=4} students \eqn{\{1,2,3,4\}}{{1,2,3,4}} 
#' sorting into groups of 2, we have \eqn{ {4 \choose 2}=6 }{choose(4,2) = 6} feasible groups: 
#' (1,2)(3,4) (1,3)(2,4) (1,4)(2,3).
#' 
#' In the group formation problem with \eqn{n=6} students \eqn{\{1,2,3,4,5,6\}}{{1,2,3,4,5,6}} 
#' sorting into groups of 3, we have \eqn{ {6 \choose 3} =20}{choose(6,3) = 20} feasible groups. 
#' For the same students sorting into groups of sizes 2 and 4, we have \eqn{ {6 \choose 2} + 
#' {6 \choose 4}=30}{choose(6,2) + choose(6,4) = 30} feasible groups.
#'
#' The structural model consists of a selection and an outcome equation. The \emph{Selection Equation} 
#' determines which matches are observed (\eqn{D=1}) and which are not (\eqn{D=0}).
#' \deqn{ \begin{array}{lcl}
#'        D &= & 1[V \in \Gamma] \\
#'        V &= & W\alpha + \eta
#'        \end{array}
#'      }{ D = 1[V in \Gamma] with V = W\alpha + \eta
#'      }
#' Here, \eqn{V} is a vector of latent valuations of \emph{all feasible} matches, ie observed and 
#' unobserved, and \eqn{1[.]} is the Iverson bracket. 
#' A match is observed if its match valuation is in the set of valuations \eqn{\Gamma}
#' that satisfy the equilibrium condition (see Klein, 2015a). This condition differs for matching
#' games with transferable and non-transferable utility and can be specified using the \code{method} 
#' argument. 
#' The match valuation \eqn{V} is a linear function of \eqn{W}, a matrix of characteristics for 
#' \emph{all feasible} groups, and \eqn{\eta}, a vector of random errors. \eqn{\alpha} is a paramter 
#' vector to be estimated.
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
#' @param x data frame with individual-level characteristics of all group members including
#' market- and group-identifiers.
#' @param m.id character string giving the name of the market identifier variable. Defaults to \code{"m.id"}.
#' @param g.id character string giving the name of the group identifier variable. Defaults to \code{"g.id"}.
#' @param R dependent variable in outcome equation. Defaults to \code{"R"}.
#' @param selection list containing variables and pertaining operators in the selection equation. The format is 
#' \code{operation = "variable"}. See the Details and Examples sections.
#' @param outcome list containing variables and pertaining operators in the outcome equation. The format is 
#' \code{operation = "variable"}. See the Details and Examples sections.
#' @param simulation should the values of dependent variables in selection and outcome equations be simulated? Options are \code{"none"} for no simulation, \code{"NTU"} for non-transferable utility matching, \code{"TU"} for transferable utility or \code{"random"} for random matching of individuals to groups. Simulation settings are (i) all model coefficients set to \code{alpha=beta=1}; (ii) covariance between error terms \code{delta=0.5}; (iii) error terms \code{eta} and \code{xi} are draws from a standard normal distribution.
#' @param seed integer setting the state for random number generation if \code{simulation=TRUE}.
#' @param max.combs integer (divisible by two) giving the maximum number of feasible groups to be used for generating group-level characteristics.
#' @param method estimation method to be used. Either \code{"NTU"} or \code{"TU"} for selection correction using non-transferable or transferable utility matching as selection rule; \code{"outcome"} for estimation of the outcome equation only; or \code{"model.frame"} for no estimation.
#' @param binary logical: if \code{TRUE} outcome variable is taken to be binary; if \code{FALSE} outcome variable is taken to be continuous.
#' @param offsetOut vector of integers indicating the indices of columns in \code{X} for which coefficients should be forced to 1. Use 0 for none.
#' @param offsetSel vector of integers indicating the indices of columns in \code{W} for which coefficients should be forced to 1. Use 0 for none.
#' @param marketFE logical: if \code{TRUE} market-level fixed effects are used in outcome equation; if \code{FALSE} no market fixed effects are used.
#' @param censored draws of the \code{delta} parameter that estimates the covariation between the error terms in selection and outcome equation are 0:not censored, 1:censored from below, 2:censored from above.
#' @param gPrior logical: if \code{TRUE} the g-prior (Zellner, 1986) is used for the variance-covariance matrix.
#' @param dropOnes logical: if \code{TRUE} one-group-markets are exluded from estimation.
#' @param interOut two-colum matrix indicating the indices of columns in \code{X} that should be interacted in estimation. Use 0 for none.
#' @param interSel two-colum matrix indicating the indices of columns in \code{W} that should be interacted in estimation. Use 0 for none.
#' @param standardize numeric: if \code{standardize>0} the independent variables will be standardized by dividing by \code{standardize} times their standard deviation. Defaults to no standardization \code{standardize=0}. 
#' @param niter number of iterations to use for the Gibbs sampler.
#' @param verbose .
#' 
#' @export
#' 
#' @useDynLib matchingMarkets, .registration = TRUE 
#' 
#' @import partitions stats
#' @importFrom Rcpp evalCpp
#' @importFrom utils setTxtProgressBar txtProgressBar
#' 
#' @details 
#' Operators for variable transformations in \code{selection} and \code{outcome} arguments.
#' \describe{
#' \item{\code{add}}{sum over all group members and divide by group size.}
#' \item{\code{int}}{sum over all possible two-way interactions \eqn{x*y} of group members and divide by the number of those, given by \code{choose(n,2)}.}
#' \item{\code{ieq}}{sum over all possible two-way equality assertions \eqn{1[x=y]} and divide by the number of those.}
#' \item{\code{ive}}{sum over all possible two-way interactions of vectors of variables of group members and divide by number of those.}
#' \item{\code{inv}}{...}
#' \item{\code{dst}}{sum over all possible two-way distances between players and divide by number of those, where distance is defined as \eqn{e^{-|x-y|}}{exp(-|x-y|)}.}
#' }
#' 
#' @author Thilo Klein 
#' 
#' @keywords regression
#' 
#' @references Klein, T. (2015a). \href{https://ideas.repec.org/p/cam/camdae/1521.html}{Does Anti-Diversification Pay? A One-Sided Matching Model of Microcredit}.
#' \emph{Cambridge Working Papers in Economics}, #1521.
#' @references Zellner, A. (1986). \emph{On assessing prior distributions and Bayesian regression analysis with g-prior distributions}, 
#' volume 6, pages 233--243. North-Holland, Amsterdam.
#' 
#' @examples
#' \dontrun{
#' ## --- SIMULATED EXAMPLE ---
#' 
#' ## 1. Simulate one-sided matching data for 200 markets (m=200) with 2 groups
#' ##    per market (gpm=2) and 5 individuals per group (ind=5). True parameters 
#' ##    in selection equation is wst=1, in outcome equation wst=0. 
#' 
#' ## 1-a. Simulate individual-level, independent variables
#'  idata <- stabsim(m=200, ind=5, seed=123, gpm=2)
#'  head(idata)
#'  
#' ## 1-b. Simulate group-level variables 
#'  mdata <- stabit(x=idata, simulation="NTU", method="model.frame",
#'  selection = list(add="wst"), outcome = list(add="wst"), verbose=FALSE)
#'  head(mdata$OUT)
#'  head(mdata$SEL)
#' 
#' 
#' ## 2. Bias from sorting
#' 
#' ## 2-a. Naive OLS estimation
#'  lm(R ~ wst.add, data=mdata$OUT)$coefficients
#' 
#' ## 2-b. epsilon is correlated with independent variables
#'  with(mdata$OUT, cor(epsilon, wst.add))
#'  
#' ## 2-c. but xi is uncorrelated with independent variables
#'  with(mdata$OUT, cor(xi, wst.add))
#' 
#' ## 3. Correction of sorting bias when valuations V are observed
#' 
#' ## 3-a. 1st stage: obtain fitted value for eta
#' lm.sel <- lm(V ~ -1 + wst.add, data=mdata$SEL)
#' lm.sel$coefficients
#' 
#' eta <- lm.sel$resid[mdata$SEL$D==1]
#' 
#' ## 3-b. 2nd stage: control for eta
#'  lm(R ~ wst.add + eta, data=mdata$OUT)$coefficients
#' 
#' 
#' ## 4. Run Gibbs sampler
#'  fit1 <- stabit(x=idata, method="NTU", simulation="NTU", censored=1, 
#'                 selection = list(add="wst"), outcome = list(add="wst"), 
#'                 niter=2000, verbose=FALSE)
#' 
#' 
#' ## 5. Coefficient table
#'  summary(fit1)
#' 
#' 
#' ## 6. Plot MCMC draws for coefficients
#'  plot(fit1)
#' 
#' 
#' ## --- REPLICATION, Klein (2015a) ---
#' 
#' ## 1. Load data 
#'  data(baac00); head(baac00)
#'  
#' ## 2. Run Gibbs sampler
#'  klein15a <- stabit(x=baac00, selection = list(inv="pi",ieq="wst"), 
#'         outcome = list(add="pi",inv="pi",ieq="wst",
#'         add=c("loan_size","loan_size2","lngroup_agei")), offsetOut=1,
#'         method="NTU", binary=TRUE, gPrior=TRUE, marketFE=TRUE, niter=800000)
#' 
#' ## 3. Marginal effects
#'  summary(klein15a, mfx=TRUE)
#'  
#' ## 4. Plot MCMC draws for coefficients
#'  plot(klein15a)
#' }
stabit <- function(x, m.id="m.id", g.id="g.id", R="R", selection=NULL, outcome=NULL, 
                   simulation="none", seed=123, max.combs=Inf,
                   method="NTU", binary=FALSE, offsetOut=0, offsetSel=0, 
                   marketFE=FALSE, censored=0, gPrior=FALSE, dropOnes=FALSE, interOut=0, interSel=0, 
                   standardize=0, niter=10, verbose=FALSE){
  
  # -----------------------------------------------------------------------------
  # Generate design matrix.
  # -----------------------------------------------------------------------------
  assignment <- simulation
  simulation <- ifelse(simulation!="none", TRUE, FALSE)
  
  data <- design.matrix(x, m.id=m.id, g.id=g.id, R=R, selection=selection, outcome=outcome, 
                        simulation, assignment=assignment, seed=seed, max.combs=max.combs,
                        standardize=standardize, verbose=verbose)
  
  # -----------------------------------------------------------------------------
  # Obtain parameter estimates.
  # -----------------------------------------------------------------------------
  
  
  sel <- ifelse(method=="outcome", FALSE, TRUE)
  NTU <- ifelse(method=="NTU", TRUE, FALSE)
  
  # -----------------------------------------------------------------------------
  # Attach variables in 'data' object.
  # -----------------------------------------------------------------------------
  D = data$D
  R = data$R
  W = data$W
  X = data$X
  P = data$P
  combs = data$combs
  
  # -----------------------------------------------------------------------------
  # Preliminaries (1).
  # -----------------------------------------------------------------------------  
  ## replace NA's with 0 in combinations matrices
  combs <- lapply(combs, function(j) replace(j,is.na(j),0))
  
  ## add intercept to design matrix of outcome equation
  X <- lapply(X, function(i) data.frame(intercept=1,i))
  
  ## variable names
  an = colnames(W[[1]])
  bn = colnames(X[[1]])
  
  # -----------------------------------------------------------------------------
  # Market identifiers.
  # -----------------------------------------------------------------------------
  T = 1:length(R) # market identifiers.
  One = vector()
  Two = vector()
  for(t in T){
    if( length(D[[t]]) == 1 ){
      One = c(One, t) # one-group markets.
    } else{
      Two = c(Two, t) # two-group markets.
    }
  }
  
  # -----------------------------------------------------------------------------
  # Interactions.
  # -----------------------------------------------------------------------------
  if( !is.null(dim(interSel)) ){
    for( i in 1:dim(interSel)[1] ){
      h = dim(W[[1]])[2]
      for( j in 1:length(Two) ){
        W[[j]][,h+1] = W[[j]][,interSel[i,1]] * W[[j]][,interSel[i,2]]
      }
      an = c( an, paste(an[interSel[i,1]],an[interSel[i,2]],sep=":") )
    }
  }
  
  if( !is.null(dim(interOut)) ){
    for( i in 1:dim(interOut)[1] ){
      h = dim(X[[1]])[2]
      for( j in 1:length(T) ){
        X[[j]][,h+1] = X[[j]][,interOut[i,1]] * X[[j]][,interOut[i,2]]
      }
      bn = c( bn, paste(bn[interOut[i,1]],bn[interOut[i,2]],sep=":") )
    }
  }
  
  # -----------------------------------------------------------------------------
  # One-group-markets.
  # -----------------------------------------------------------------------------
  if(dropOnes == TRUE  | length(One) == 0){
    One = NULL # drop one-group markets ???
    T = Two
  } else{
    for(t in Two){
      X[[t]] = cbind(X[[t]], 0)
    }
    for(t in One){
      X[[t]] = cbind(X[[t]], 1)
      names(X[[t]]) = names(X[[1]])
    }
    bn = c(bn, "one")
  }
  
  # -----------------------------------------------------------------------------
  # Fixed effects.
  # -----------------------------------------------------------------------------
  if(marketFE == TRUE){
    for(t in T){
      X[[t]] = X[[t]][,-1] # drop intercept to avoid perfect multicollinearity
    }
    bn = bn[-1]
    if(binary == TRUE){
      count1 = 0
      for(t in Two){
        if( R[[t]][1] != R[[t]][2] ){ # both groups in market have different outcome.
          count1 = count1+1
        }   
      }
      count2 = 0
      for(t in Two){
        X[[t]] = cbind( X[[t]], matrix(0,nrow=2,ncol=count1) )
        if( R[[t]][1] != R[[t]][2] ){ # both groups in market have different outcome.
          end = dim(X[[t]])[2]
          X[[t]][,end-count2] = rbind(1,1)
          count2 = count2+1
          bn = c(bn, paste("d",count2,sep=""))
        } 
      }
      for(t in One){
        X[[t]] = cbind( X[[t]], matrix(0,nrow=1,ncol=count1) )
      }        
    } else{ # continuous outcome variable.
      count2 = 0
      for(t in Two){
        X[[t]] = cbind( X[[t]], matrix(0,nrow=2,ncol=length(Two)) )
        end = dim(X[[t]])[2]
        X[[t]][,end-count2] = rbind(1,1)
        count2 = count2+1
        bn = c(bn, paste("d",count2,sep=""))
      }
      for(t in One){
        X[[t]] = cbind( X[[t]], matrix(0,nrow=1,ncol=length(Two)) )
      }
    }
  } 
  
  # -----------------------------------------------------------------------------
  # Offset outcome equation.
  # -----------------------------------------------------------------------------
  offOut = list()
  if( sum(offsetOut) != 0 ){
    for(t in T){ # set offset and drop column from X
      offOut[[t]] = t(t(apply(as.matrix(X[[t]][,offsetOut]),1,sum)))
      X[[t]] = X[[t]][,-offsetOut]
    } 
    bn = bn[-offsetOut]
  } else{
    for(t in T){ # set offset to zero
      offOut[[t]] = matrix(0,nrow=dim(X[[t]])[1],ncol=1)
    }
  }
  
  # -----------------------------------------------------------------------------
  # Offset selection equation.
  # -----------------------------------------------------------------------------
  offSel = list()
  if( sum(offsetSel) != 0 ){
    for(t in Two){ # set offset and drop column from W
      offSel[[t]] = t(t(apply(as.matrix(W[[t]][,offsetSel]),1,sum)))
      W[[t]] = W[[t]][,-offsetSel]
    } 
    an = an[-offsetSel]
  } else{
    for(t in Two){ # set offset to zero
      offSel[[t]] = matrix(0,nrow=dim(W[[t]])[1],ncol=1)
    } 
  }
  
  # -----------------------------------------------------------------------------
  # Size of feasible groups.
  # -----------------------------------------------------------------------------
  Uneq = vector()
  l = list()
  p = list()
  for(t in Two){
    l[[t]] = length(D[[t]])
    p[[t]] = (l[[t]]/2)+1 # cut point between "ego" and "other" groups
    if( min(combs[[t]]) == 0 ){ # unequal group sizes in market.
      Uneq = c(Uneq, t)
    }
  }
  
  # -----------------------------------------------------------------------------
  # Indicator for group size difference.
  # -----------------------------------------------------------------------------
  if(NTU == TRUE){
    s = dim(W[[1]])[2]
    for(t in Two){
      W[[t]] = cbind( W[[t]], matrix(0,nrow=l[[t]],ncol=length(Uneq)) )
    }
    count = 0
    for(t in Uneq){
      count = count + 1
      W[[t]][,s+count] = matrix(c(1,0,rep(1,p[[t]]-2),rep(0,p[[t]]-2)),ncol=1) 
      an = c(an, paste("d",count,sep=""))
    }
  }
  
  # -----------------------------------------------------------------------------
  # Preliminaries (2).
  # -----------------------------------------------------------------------------
  data$X <- X
  data$W <- W
  
  W <- lapply(W,as.matrix)
  X <- lapply(X,as.matrix)
  
  if(!is.null(One) & !is.null(Two)){
    T = length(One) + length(Two) 
    One = c( One[1]-1, One[1]-1+length(One) )
    Two = c( 0, length(Two) )
  } else{
    if(is.null(One)){
      T = length(Two)
      One = c(0,0)
      Two = c( 0, length(Two) )
    }
  }
  
  n <- dim(do.call(rbind,X))[1]
  N <- dim(do.call(rbind,W))[1]
  
  sigmabarbetainverse <- ( t(do.call(rbind,X)) %*% do.call(rbind,X) ) / n
  sigmabaralphainverse <- ( t(do.call(rbind,W)) %*% do.call(rbind,W) ) / N
  
  if(method != "model.frame"){
    
    # -----------------------------------------------------------------------------
    # Source C++ script
    # -----------------------------------------------------------------------------    
    
    est <- stabitSel2(Xr=X, Rr=R, Wr=W, One=One, Two=Two, T=T, 
                      offOutr=offOut, offSelr=offSel,
                      sigmabarbetainverse=sigmabarbetainverse, sigmabaralphainverse=sigmabaralphainverse,
                      niter=niter, n=n, l=matrix(unlist(l),length(l),1), Pr=P, p=matrix(unlist(p),length(p),1),
                      binary=binary, selection=sel, censored=censored, gPrior=gPrior, ntu=NTU)
    
    # -----------------------------------------------------------------------------
    # Add names to coefficients.
    # ----------------------------------------------------------------------------- 
    
    ## variable names
    
    an <- colnames(X[[1]])
    bn <- colnames(W[[1]])
    
    # parameter draws
    
    rownames(est$alphadraws) = an
    if(sel==TRUE){
      est$alphadraws <- with(est, rbind(alphadraws, deltadraws))
      rownames(est$alphadraws) <- c(an,"delta")
      est$deltadraws <- NULL
      rownames(est$betadraws) = bn
    }
    
    # posterior means
    
    ## The last half of all draws are used in approximating the posterior means and the posterior standard deviations.
    niter <- ncol(est$alphadraws)
    startiter <- floor(niter/2)
    
    posMeans <- function(x){
      sapply(1:nrow(x), function(z) mean(x[z,startiter:niter]))
    }
    
    est$alpha <- posMeans(est$alphadraws) 
    names(est$alpha) <- rownames(est$alphadraws)
    
    if(sel==TRUE){
      est$beta <- posMeans(est$betadraws)
      names(est$beta) <- bn  
    } 
    if(binary==FALSE){
      est$sigma <- sqrt(posMeans(est$sigmasquarexidraws))
      est$sigmasquarexidraws <- NULL
    }
    
    ## error terms
    
    if(sel==TRUE){
      est$eta <- posMeans(est$etadraws)
      est$etadraws <- NULL
    }
    
    ## vcov
    
    est$vcov$alpha <- var(t(est$alphadraws))
    if(sel==TRUE){
      est$vcov$beta <- var(t(est$betadraws))
    }
    
    ## variables
    
    est$X <- do.call(rbind, X)
    est$Y <- do.call(c, R)
    if(sel==TRUE){
      est$X <- cbind(est$X, eta=c(est$eta, rep(0,nrow(est$X)-length(est$eta))) )
      est$W <- do.call(rbind, W)
      est$D <- do.call(c, D)
    }
    
    ## coefficients
    
    if(sel==TRUE){
      est$coefficients <- unlist(with(est, list(o=alpha, s=beta)))  
    } else{
      est$coefficients <- est$alpha
    }
    
    est$fitted.values <- with(est, as.vector(X %*% alpha))
    est$residuals <- est$Y - est$fitted.values
    est$df <- n-ncol(est$X)
    
    est$call <- match.call()
    est$method <- ifelse(sel==FALSE, "Outcome-only", "Selection")
    est$binary <- binary
    #est$formula <- ???
    
    ## consolidate
    
    h <- est[grep(pattern="draws",x=names(est))]
    est[grep(pattern="draws",x=names(est))] <- NULL
    est$draws <- h
    
    h <- est[which(names(est) %in% c("alpha","beta","delta"))]
    est[which(names(est) %in% c("alpha","beta","delta"))] <- NULL
    est$coefs <- h
    
    h <- est[which(names(est) %in% c("X","W","Y","D"))]
    est[which(names(est) %in% c("X","W","Y","D"))] <- NULL
    est$variables <- h
    
    h <- NULL
    
    # -----------------------------------------------------------------------------
    # Returns .
    # ----------------------------------------------------------------------------- 
    
    class(est) <- "stabit2"
    return(est)
    
  } else{ ## if method == "model.frame"
    
    # -----------------------------------------------------------------------------
    # Returns .
    # ----------------------------------------------------------------------------- 
    
    model.frame <- unlistData(x=data)
    list(OUT=model.frame$OUT, SEL=model.frame$SEL, combs=combs)
    #return(list(model.list=data, model.frame=model.frame))
    
  }
}




unlistData <- function(x){
  
  # --------------------------------------------------------------------
  # R-code (www.r-project.org) to unlist the matching data for OLS regression
  
  # The arguments of the function are:
  # x        : matching data after applying the oneSidedMatching() function
  
  # ## Examples:
  #
  # unlistData(x)
  # --------------------------------------------------------------------
  
  ## Outcome equation
  if(is.null(dim(x$X[[1]]))){
    countX <- unlist(lapply(x$X,length)) # count of groups per market
  } else{
    countX <- unlist(lapply(x$X,nrow)) # count of groups per market
  }
  g.id   <- unlist(c(sapply(countX, function(i) 1:i )))
  m.id   <- unlist(c(sapply(1:length(countX), function(i) rep(i,countX[i]))))
  if(length(unlist(x$xi)) == length(g.id)){ # if simulation = TRUE in model.matrix
    OUT    <- with(x, data.frame( m.id, g.id, do.call(rbind,X), R=do.call(c,R),
                                  xi=do.call(c,xi), epsilon=do.call(c,epsilon) ))
  } else{
    OUT    <- with(x, data.frame( m.id, g.id, do.call(rbind,X), R=do.call(c,R) ))    
  }
  
  ## Selection equation
  if(is.null(dim(x$W[[1]]))){ # only one variable in W
    countW <- unlist(lapply(x$W,length)) # count of groups per market
  } else{
    countW <- unlist(lapply(x$W,nrow)) # count of groups per market
  }
  g.id   <- unlist(c(sapply(countW, function(i) 1:i )))
  m.id   <- unlist(c(sapply(1:length(countW), function(i) rep(i,countW[i]))))
  h <- unlist(lapply(x$D, length))
  x$D <- x$D[which(h>1)] # for 2-group markets only
  if(length(unlist(x$V)) == length(g.id)){ # if simulation = TRUE in model.matrix
    if(is.null(dim(x$W[[1]]))){ # only one variable in W
      SEL <- with(x, data.frame( m.id, g.id, reg.id=na.omit(do.call(c,P)), do.call(c,W), 
                                 D=do.call(c,D), V=do.call(c,V), eta=do.call(c,eta) ))
    } else{
      SEL <- with(x, data.frame( m.id, g.id, reg.id=na.omit(do.call(c,P)), do.call(rbind,W), 
                                 D=do.call(c,D), V=do.call(c,V), eta=do.call(c,eta) ))      
    }
  } else{
    if(is.null(dim(x$W[[1]]))){ # only one variable in W
      SEL    <- with(x, data.frame( m.id, g.id, reg.id=na.omit(do.call(c,P)), do.call(c,W), D=do.call(c,D) ))    
    } else{
      SEL    <- with(x, data.frame( m.id, g.id, reg.id=na.omit(do.call(c,P)), do.call(rbind,W), D=do.call(c,D) ))    
    }
  }
  
  return(list(OUT=OUT, SEL=SEL))
}




design.matrix <- function(x, m.id="m.id", g.id="g.id", R="R", selection=NULL, outcome=NULL, 
                          simulation=FALSE, assignment="NTU", seed=123, max.combs=Inf,
                          standardize=0, verbose=verbose){  
  
  # --------------------------------------------------------------------
  # R-code (www.r-project.org) to obtain a design matrix for the analysis
  # of data from the coalition formation game.
  
  # Calls the following functions in the given order:
  # 1. listByMarket
  # 2. speci
  # 3. combmats.coalitions / combmats.roomates
  # 4. indexmat
  # 5. designmatrix (-> wrap -> combmats.interactions)
  # --------------------------------------------------------------------
  
  ## rename dependent variable, group and market identifiers
  names(x)[names(x) == R]    <- "R"
  names(x)[names(x) == g.id] <- "g.id"
  names(x)[names(x) == m.id] <- "m.id"
  
  ## list data by market id
  x <- listByMarket(x=x, m.id=m.id, g.id=g.id, R=R)
  
  ## create combinatorial matrices
  spec <- speci(x=x)
  
  CMATS <- combmats.coalitions(x=spec, dodrop=TRUE, max.combs=max.combs) 
  
  INDEXMAT <- indexmat(spec)
  
  ## create design matrix
  designmatrix(selection=selection, outcome=outcome, x=x,
               simulation=simulation, assignment=assignment, seed=seed, spec=spec, CMATS=CMATS,
               INDEXMAT=INDEXMAT, standardize=standardize, verbose=verbose)
}




combmats.coalitions <- function(x=NULL, n1, n2, dodrop=FALSE, max.combs=Inf){
  
  # --------------------------------------------------------------------
  # R-code (www.r-project.org) for obtaining all possible groups when 
  # partitioning the market into TWO groups of size 'n1' and 'n2'.
  # The resulting matrix has number of rows equal to all possible groups
  # with each each row containing the group members' indices.
  
  # The arguments of the function are either:
  # n1 : integer indicating the number of individuals in group 1
  # n2 : integer indicating the number of individuals in group 2
  #
  # OR:
  # x : vector indicating the number of individuals in group 1 and 2,
  #     or a matrix with two columns indicating the number of individuals 
  #     in group 1 and 2, and number of rows equal to the number of markets
  # 
  # AND:
  # dodrop    : logical, should the two observed groups be dropped?
  # max.combs : integer (divisible by two) giving the maximum number of feasible 
  #             groups to be used for generating group-level characteristics.
  
  
  # ## Examples:
  #
  # combmats.coalitions(n1=2, n2=2)
  # combmats.coalitions(n1=2, n2=2, dodrop=TRUE)
  # combmats.coalitions(x=c(2,3))
  # combmats.coalitions(x=matrix(1:4, ncol=2, nrow=2, byrow=TRUE))
  # --------------------------------------------------------------------
  
  #library(partitions)
  set.seed(123)
  
  ## If 'x' not given, obtain it from 'n1' and 'n2'
  if(is.null(x)){
    x <- c(n1,n2)
  }
  
  ## Consistency checks:
  if(is.matrix(x)==FALSE){
    if(length(x)>2){stop("combmats() does only handle two groups per market!")}
    x <- matrix(x, ncol=2)
  } else{
    if(dim(x)[2]>2){stop("combmats() does only handle two groups per market!")}
  }
  
  ## sufficient to obtain partitions for unique elements of 'x'  
  x <- t(sapply(1:dim(x)[1], function(j) sort(x[j,]) ))  # order doesn't matter
  x <- unique(x)  # drop duplicates
  
  lapply(1:dim(x)[1], function(i){
    
    if( choose(sum(x[i,]), x[i,1]) < max.combs ){ ## if number of feasible groups below limit
      L  <- listParts(x[i,])
      L1 <- t(sapply(1:length(L), function(z) c(L[[z]][2]$`2`,rep(NA,abs(diff(x[i,])))) ))
      L2 <- t(sapply(1:length(L), function(z) L[[z]][1]$`1` ))
      LC <- rbind(L1,L2)
      
      if(dodrop){
        ## prepare variables needed to drop indices of the two observed groups in the market
        numA <- x[i,1]; rangeA <- 1:numA
        numB <- x[i,2]; rangeB <- (numA+1):(numB+numA)      
        ncombs <- dim(LC)[1]
        num    <- c(rep(numA,ncombs/2), rep(numB,ncombs/2))
        
        todrop <- NA; s <- 1  # two write the rows of indices of the 2 observed groups to be dropped
        for(j in 1:ncombs){
          indices <- LC[j,1:num[j]]
          if( (all.equal(indices, rangeA)==TRUE) | (all.equal(indices, rangeB)==TRUE) ){
            todrop[s] <- j
            s <- s+1
          }
        }
        LC <- LC[-todrop,]
      } else{
        LC <- LC
      }
    } else{
      A <- t(sapply(1:max.combs, function(z) sort(sample(1:sum(x[i,]),x[i,1],replace=F))))
      A <- A[apply(A,1,function(z) all.equal(z, 1:x[i,1])!=TRUE),] ## drop equilibrium groups
      dim(unique(A))
      A <- unique(A)[1:(max.combs/2),]
      B <- t(apply(A,1,function(z) sort(c(1:sum(x[i,]))[-z])))
      LC <- rbind(cbind(A,matrix(NA,nrow=max.combs/2,ncol=dim(B)[2]-dim(A)[2])),B)
    }
  })
}




indexmat <- function(x=NULL){
  
  # --------------------------------------------------------------------
  # R-code (www.r-project.org) for setting up an index matrix for quick 
  # indexing to the list of combination matrices (combmats.coalitions)
  # based on the number and size of groups in a market.
  
  # The arguments of the function are:
  # x : Either: a matrix with columns indicating the number of individuals in
  #             group 1 and 2, and number of rows equal to the number of markets.
  #     Or: a vector of length equal to the number of markets, giving the 
  #         number of players in each market.
  
  # ## Examples:
  #
  # ## index matrices for coalition formation game
  # indexmat(x=matrix(c(1,2,2,1,3,3), ncol=2, nrow=3))
  # indexmat(x=matrix(1:6, ncol=2, nrow=3))
  # --------------------------------------------------------------------
  
  ## Consistency checks:
  if(is.matrix(x)){
    if(dim(x)[2]!=2){stop("'x' must have 2 columns!")}
    
    ## sufficient to obtain partitions for unique elements of 'x'  
    x <- t(sapply(1:dim(x)[1], function(j) sort(x[j,]) ))  # order doesn't matter
    x <- unique(x)  # drop duplicates  
    
    mat <- matrix(NA, nrow=max(x), ncol=max(x))
    for(i in 1:dim(x)[1]){ 
      mat[x[i,1],x[i,2]] <- i
      mat[x[i,2],x[i,1]] <- i
    }
    return(mat)
  } else{
    x <- unique(x)  # drop duplicates
    mat <- rep(NA,max(x))
    for(i in 1:length(x)){
      mat[x[i]] <- i
    }
    return(mat)
  }
}




listByMarket <- function(x, m.id, g.id, R){
  
  # --------------------------------------------------------------------
  # R-code (www.r-project.org) for creating a list of markets from 
  # one-sided matching data.
  
  # The arguments of the function are:
  # x    : data frame 
  # m.id : variable name of market identifier
  # g.id : variable name of group identifier
  # R    : dependent variable
  
  # ## Examples:
  #
  # listByMarket(x=data, g.id="g.id", m.id="m.id")
  # --------------------------------------------------------------------
  
  ## Consistency checks:
  if(is.null(g.id) | is.null(m.id)){stop("variables 'g.id' and/or 'm.id' missing!")}
  
  ## data.frame to list
  x <- split(x, x$m.id)
  
  ## add group size (g.size), number of groups (g.numb), and individual identifier (i.id)
  x <- lapply(x, function(i) data.frame(i, g.size=rep(table(i$g.id),table(i$g.id)),
                                        g.numb=length(unique(i$g.id))))
  
  ## normalise g.id to counts staring from 1,2,...
  x <- lapply(x, function(i) data.frame(g.id=as.numeric(factor(i$g.id)), i.id=1:dim(i)[1], m.id=i$m.id, 
                                        i[-which(names(i)%in%c("m.id","i.id","g.id"))]) )
  
  ## sort list such that two-group markets are listed first
  l <- unlist(lapply(x, function(i) length(unique(i$g.id))))
  x <- x[order(l,decreasing=TRUE)]
  
  ## sort groups in each market by group size
  x <- lapply(x, function(i) i[order(i$g.size,decreasing=FALSE),])
  
  return(x)
}




speci <- function(x){
  
  # --------------------------------------------------------------------
  # R-code (www.r-project.org) for creating either: a group size matrix with
  # number of rows equal to the number of markets and 2 columns indicating the
  # number of players in group 1 and 2 respectively; or: a vector of length 
  # equal to the number of markets with each element giving the market size.
  
  # The arguments of the function are:
  # x         : data list from listByMarket()
  
  # ## Examples:
  #
  # speci(x=data)
  # --------------------------------------------------------------------
  
  spec <- unlist(lapply(1:length(x), function(i){
    if(2 %in% x[[i]]$g.numb){  # only for 2-group markets
      table(x[[i]]$g.id )
    }
  }))
  spec <- matrix(spec, byrow=TRUE, ncol=2)
  spec <- t(sapply(1:dim(spec)[1], function(j) sort(spec[j,]) ))  # order doesn't matter
  return(spec)
}




combmats.interactions <- function(x=NULL){ 
  
  # --------------------------------------------------------------------
  # R-code (www.r-project.org) for obtaining all possible groups when 
  # partitioning a market into groups of size TWO.
  
  # The arguments of the function are:
  # x : integer indicating the indices of selected individuals in the market
  
  # ## Examples:
  #
  # combmats.interactions(x=c(1,3,5))
  # --------------------------------------------------------------------
  
  ## Consistency checks:
  if(is.null(dim(x))==FALSE){stop("'x' must be of dimension NULL!")}
  
  l <- length(x)
  s <- NULL
  for(j in 1:(l-1)){
    for(k in (1+j):l){
      s <- rbind(s, x[c(j,k)])
    }
  }
  return(s)
}




wrap <- function(thisdata, indices, num, denom, j, names.xw){
  
  # --------------------------------------------------------------------
  # R-code (www.r-project.org) to wrap the below group variable 
  # transformations.
  
  # The arguments of the function are:
  # thisdata : matching data for a particular market
  # indices  : indices of the group members
  # num      : group size (= number of indices)
  # denom    : number of possible ways to draw pairs from the group with size 'num'
  # j        : index running from 1 to number of possible groups in the market
  # names.xw : variable names in outcome (x) and selection (w) equation
  
  # ## Examples:
  # 
  # wrap()
  # --------------------------------------------------------------------
  
  comb <- combmats.interactions(x=indices)  # for pipj terms
  varq <- sapply(1:length(names.xw), function(i) strsplit(names.xw[i],split="@")[[1]])[1,]
  funq <- sapply(1:length(names.xw), function(i) strsplit(names.xw[i],split="@")[[1]])[2,]
  
  l <- 1:length(varq)
  sapply(l, function(a){
    if(funq[a]=="add"){
      sum(thisdata[indices,varq[a]], na.rm=TRUE)/num[j]
    } else if(funq[a]=="int"){
      sum( thisdata[comb[,1],varq[a]] * thisdata[comb[,2],varq[a]] ) / denom[j]
    } else if(funq[a]=="inv"){
      sum( thisdata[comb[,1],varq[a]] * (1-thisdata[comb[,2],varq[a]]) +
             thisdata[comb[,2],varq[a]] * (1-thisdata[comb[,1],varq[a]]) ) / (2*denom[j])
    } else if(funq[a]=="dst"){
      sum( exp( -1 * abs( thisdata[comb[,1],varq[a]] - thisdata[comb[,2],varq[a]] ))) / denom[j]
    } else if(funq[a]=="iln"){
      sum( log(thisdata[comb[,1],varq[a]] * thisdata[comb[,2],varq[a]]) ) / denom[j]
    } else if(funq[a]=="ieq"){
      sum( thisdata[comb[,1],varq[a]] == thisdata[comb[,2],varq[a]] ) / denom[j]
    } else if(funq[a]=="ive"){
      n <- nchar(names(thisdata))
      posis <- which( substr(names(thisdata), 1, n-1) == varq[a])
      sum( diag(as.matrix(thisdata[comb[,1],posis]) %*% t(as.matrix(thisdata[comb[,2],posis]))) ) / denom[j]
    } else if(funq[a]=="iem"){
      n <- nchar(names(thisdata))
      posis <- which( substr(names(thisdata), 1, n-1) == varq[a])
      fu <- function(a,b) apply(cbind(a,b), 1, function(z) z[1:(length(z)/2)] %in% z[(length(z)/2):length(z)])
      sum( fu(a=thisdata[comb[,1],posis], b=thisdata[comb[,2],posis]) ) / denom[j] / 2
    } else if(funq[a]=="sel"){
      thisdata[indices[1],varq[a]]
    } else if(funq[a]=="oth"){
      thisdata[indices[2],varq[a]]
    } else if(funq[a]=="avg"){
      mean(thisdata["<me>",varq[a]])
    } else if(funq[a]=="val"){
      n <- nchar(names(thisdata))
      ##posis = which( substr(names(thisdata), 1, n-1) == varq[a] )
      posis <- which( names(thisdata) %in% paste(varq[a],1:1000,sep="") )
      thisdata[indices[1], posis[indices[2]]]
    } else if(funq[a]=="val2"){
      n <- nchar(names(thisdata))
      posis <- which( names(thisdata) %in% paste(varq[a],1:1000,sep="") )
      thisdata[indices[1], posis[indices[2]]] + thisdata[indices[2], posis[indices[1]]]
    } else if(funq[a]=="ieg"){
      n <- nrow(thisdata)
      ( sum( ( table(c(1,1,2)) / n )^2 ) - 1/n ) / (1 - 1/n)
      ( sum( ( table(thisdata[,varq[a]]) / n )^2 ) - 1/n ) / (1 - 1/n)
    } else(stop("function must be either of 'add', 'int', 'ieq', 'ive' or 'iem'!"))
  })
}




designmatrix <- function(selection, outcome, x, simulation=FALSE, assignment, 
                         seed, spec, CMATS, INDEXMAT, standardize, verbose=verbose){ 
  
  # --------------------------------------------------------------------
  # R-code (www.r-project.org) to set up the design matrix for the analysis
  # of matching data.
  
  # The arguments of the function are:
  # selection : list ...
  # outcome   : list ...
  # x         : list from listByMarket
  # simulation: logical: if TRUE the dependent variables of selection and outcome equations are simulated
  # assignment: character string giving the assignment procedure to be used if simulation=TRUE.
  # seed      : integer setting the state for random number generation if simulation=TRUE.
  # spec      : object from function speci(). Either: a group size matrix with number of rows
  #             equal to the number of markets and 2 columns indicating the number of players
  #             in grou 1 and 2, respectively. Or: a vector of length equal to the number of
  #             markets with each element giving the market size.
  # CMATS     : combination matrix from function combmats.coalitions()
  # INDEXMAT  : index matrix from function indexmat() for quick indexing to the list of
  #             combination matrices based on the number and size of groups in a market
  # 
  # Options for variable transformations:
  # add : sum over group observations and divide by group size
  # int : sum all possible two-way interactions of group members
  #       and divide by the number of those (=choose(n,2))
  # ieq : sum over all possible two-way equality assertions and
  #       divide by the number of those (=choose(n,2))
  # ive : sum over all possible two-way interactions of vectors
  #       of variables of group members and divide by number of those
  # sel : variable for individual (peer effects only!)
  # oth : variable for other in the group (peer effects only!)
  # avg : variable average for others in the group (peer effects only, not yet implemented!)
  # dst : sum over all possible two-way distances between players and divide by
  #       number of those (= choose(n,2)) where distance is defined as exp(-abs(x1-x2))
  # iln : sum all possible two-way log-interactions of group members
  #       and divide by the number of those (=choose(n,2))
  # val : valuation of player 1 (row) for player 2 (column)
  # val2: sum of valuation of players 1 and 2
  
  # ## Examples:
  #
  # ## Thai group lending paper
  # designmatrix( selection = list(add="pi", int="pi"), 
  #   outcome = list(add="pi", int="pi", ieq="wst", ive="occ"), data=data)
  # ## Simulation for coalition formation game
  # designmatrix( selection = list(add="pi", int="pi"), 
  #   outcome = list(add="pi", int="pi"), data=data)
  # --------------------------------------------------------------------
  
  X.ind <- outcome
  W.ind <- selection
  
  ## Selected variables
  X.names <- unlist(sapply(1:length(X.ind), function(i) paste(X.ind[[i]], names(X.ind[i]), sep="@")))
  W.names <- unlist(sapply(1:length(W.ind), function(i) paste(W.ind[[i]], names(W.ind[i]), sep="@")))
  vars    <- unique(c(X.names, W.names))
  
  ## Prepare data frames
  numvills <- length(x)
  data.combs <- D <- R <- V <- P <- E <- combs <- xi <- eta <- epsilon <- 
    lapply(1:numvills, function(i) NA)
  
  cat("Generating group-level data for", numvills, "markets...","\n")
  if(verbose==TRUE){
    pb <- txtProgressBar(style = 3)
  }
  
  for(i in 1:numvills){
    
    if(i <= dim(spec)[1]){  ## 2-GROUP MARKETS
      
      ## DEFINE AUXILIARY VARIABLES
      
      # Note: numA <= numB
      numA <- spec[i,1];  rangeA <- 1:numA
      numB <- spec[i,2];  rangeB <- (numA+1):(numB+numA)
      
      thisdata <- x[[i]]
      # unobs + obs groups:
      thiscmat <- rbind( c(rangeA, rep(NA,numB-numA)), rangeB, CMATS[[INDEXMAT[numA, numB]]] )  
      ncombs   <- dim(thiscmat)[1]  
      data.combs[[i]]        <- data.frame(matrix(NA,ncol=length(vars),nrow=ncombs))  
      names(data.combs[[i]]) <- vars
      # if numA != numB, the first 1/2 contain NA:
      num   <- c(numA, numB, rep(numA,(ncombs-2)/2), rep(numB,(ncombs-2)/2))  
      denom <- c(choose(numA,2), choose(numB,2), rep(choose(numA,2),(ncombs-2)/2), rep(choose(numB,2),(ncombs-2)/2))
      
      ## OBSERVABILITY INDICATOR, 'D'
      D[[i]] <- c(rep(1,2), rep(0,ncombs-2))
      
      ## OUTCOME VARIABLE, 'R'
      #if( (abs(diff(range(thisdata[rangeA,"R"]))) > 0.1) | (abs(diff(range(thisdata[rangeB,"R"]))) > 0.1) ){
      #  stop("dependent variable must be group-level!")
      #}
      R[[i]] <- c(thisdata[rangeA,"R"][1], thisdata[rangeB,"R"][1])#, rep(NA,ncombs-2))
      
      ## OBSERVED GROUPS
      for(j in 1:2){
        indices <- thiscmat[j,1:num[j]] 
        #data.combs[[i]][j,X.names] <- wrap(thisdata, indices, num, denom, j, X.names)    
        data.combs[[i]][j,vars] <- wrap(thisdata, indices, num, denom, j, vars)    
      }
      
      ## UNOBSERVED GROUPS  
      for(j in 3:ncombs){
        indices <- thiscmat[j,1:num[j]] 
        data.combs[[i]][j,W.names] <- wrap(thisdata, indices, num, denom, j, W.names)    
      }
      
      ## Replace "@" by "." in variable names
      names(data.combs[[i]]) <- gsub("@",".",names(data.combs[[i]]))
      
      ## PARTNER GROUP INDEX IN SAME MARKET
      l             <- dim(thiscmat)[1]
      P[[i]]        <- c(2, 1, ((l/2)+2):l, 3:((l/2)+1))
      names(P[[i]]) <- NULL
      
      ## COMBINATION MATRICES
      combs[[i]] <- thiscmat
      
      
    } else{  ## 1-GROUP MARKETS
      
      ## DEFINE AUXILIARY VARIABLES
      thisdata <- x[[i]]
      indices  <- x[[i]]$i.id
      num      <- length(indices)
      denom    <- choose(num,2)
      
      data.combs[[i]]        <- data.frame(matrix(NA,ncol=length(vars),nrow=1))  
      names(data.combs[[i]]) <- vars
      
      ## SINGLE GROUP
      j <- 1
      data.combs[[i]][j,X.names] <- wrap(thisdata, indices, num, denom, j, X.names)
      
      ## Replace "@" by "." in variable names
      names(data.combs[[i]]) <- gsub("@",".",names(data.combs[[i]]))
      
      ## OBSERVABILITY INDICATOR, 'D', AND OUTCOME VARIABLE, 'R'
      D[[i]] <- 1
      R[[i]] <- thisdata[1,"R"]
      
    }
    if(verbose==TRUE){
      setTxtProgressBar(pb, i/numvills)
    }
  }
  
  
  #####################
  ## STANDARDIZATION ##
  #####################
  
  if(standardize > 0){
    
    ## standardize variance of exogneous variables to 1
    std <- apply(do.call(rbind, data.combs), 2, sd)
    for(i in 1:numvills){
      data.combs[[i]] <- data.combs[[i]] / (standardize*std)
    }
    
  }
  
  #################
  ## SIMULATIONS ##
  #################
  
  if(simulation == TRUE){
    
    set.seed(seed)
    
    for(i in 1:numvills){
      
      ncombs <- dim(combs[[i]])[1]
      thiscmat <- combs[[i]]
      l      <- dim(thiscmat)[1]
      
      if(i <= dim(spec)[1]){  ## 2-GROUP MARKETS
        
        # Note: numA <= numB
        numA <- spec[i,1];  rangeA <- 1:numA
        numB <- spec[i,2];  rangeB <- (numA+1):(numB+numA) 
        
        # if numA != numB, the first 1/2 contain NA:
        num   <- c(numA, numB, rep(numA,(ncombs-2)/2), rep(numB,(ncombs-2)/2))  
        denom <- c(choose(numA,2), choose(numB,2), rep(choose(numA,2),(ncombs-2)/2), rep(choose(numB,2),(ncombs-2)/2))
        
        ## EQUILIBRIUM GROUP SELECTION
        xi[[i]]      <- rnorm(ncombs)
        eta[[i]]     <- rnorm(ncombs)
        delta        <- 1
        epsilon[[i]] <- xi[[i]] + delta*eta[[i]]
        
        V[[i]] <- apply(data.combs[[i]], 1, sum) + eta[[i]]
        
        if(assignment=="NTU"){
          
          ## A: NON-TRANSFERABLE UTILITY
          equ1 <- which(V[[i]]==max(V[[i]]))[1]
          equ2 <- P[[i]][equ1]
          
        } else if(assignment=="TU"){
          
          ## B: TRANSFERABLE UTILITY
          market.value <- V[[i]] + V[[i]][P[[i]]]
          equ1 <- which(market.value == max(market.value))[1]
          equ2 <- P[[i]][equ1]
          
        } else if(assignment=="random"){
          
          ## C: RANDOM GROUP ASSIGNMENT
          equ1 <- sample(1:ncombs, 1)
          equ2 <- P[[i]][equ1]
          
        } else(stop("assignment must be either of 'NTU', 'TU' or 'random'!"))
        
        ## Swap groups at position 1 and 2 with equilibrium groups equ1 and equ2
        if(sum(equ1,equ2)!=3){ ## otherwise (equ1, equ2) are already in position (1,2)
          #data.combs[[i]] <- rbind( data.combs[[i]][c(equ1,equ2),], 
          #                          data.combs[[i]][1,], data.combs[[i]][(3:(l/2+1))[-(min(equ1,equ2)-2)],], 
          #                          data.combs[[i]][2,], data.combs[[i]][((l/2+2):l)[-(min(equ1,equ2)-2)],] )
          nvars <- length(vars)
          data.combs[[i]][1:dim(data.combs[[i]])[1],] <- rbind( as.matrix( data.combs[[i]][c(equ1,equ2),], ncol=nvars), 
                                                                data.combs[[i]][1,], 
                                                                as.matrix( data.combs[[i]][(3:(l/2+1))[-(min(equ1,equ2)-2)],], ncol=nvars), 
                                                                data.combs[[i]][2,], 
                                                                as.matrix( data.combs[[i]][((l/2+2):l)[-(min(equ1,equ2)-2)],], ncol=nvars) )
          #colnames(data.combs[[i]]) <- vars
          #colnames(data.combs[[i]]) <- gsub("@",".",names(data.combs[[i]]))
          
          V[[i]]          <- c( V[[i]][c(equ1,equ2)], 
                                V[[i]][1], V[[i]][(3:(l/2+1))[-(min(equ1,equ2)-2)]], 
                                V[[i]][2], V[[i]][((l/2+2):l)[-(min(equ1,equ2)-2)]] )
          
          eta[[i]]        <- c( eta[[i]][c(equ1,equ2)], 
                                eta[[i]][1], eta[[i]][(3:(l/2+1))[-(min(equ1,equ2)-2)]], 
                                eta[[i]][1], eta[[i]][((l/2+2):l)[-(min(equ1,equ2)-2)]] )
        }
        xi[[i]]      <- xi[[i]][c(equ1,equ2)]
        epsilon[[i]] <- epsilon[[i]][c(equ1,equ2)]
        R[[i]] <- 0*apply(as.matrix( data.combs[[i]][1:2,], nrow=2), 1, sum) + epsilon[[i]]
        #R[[i]] <- ifelse(R[[i]] > 1, 1, 0) # uncomment me!
        
        E[[i]] <- thiscmat[c(equ1,equ2),]
        
        
      } else{  ## 1-GROUP MARKETS
        
        xi[[i]]  <- rnorm(1)
        eta[[i]] <- rnorm(1)          
        delta        <- 0.5
        epsilon[[i]] <- delta*eta[[i]] + xi[[i]]
        
        R[[i]] <- sum(data.combs[[i]]) + epsilon[[i]]
        
      }
    }
  }
  
  ###################
  ## WRITE RESULTS ##
  ###################
  
  ## Replace "@" by "." in W.names and X.names
  W.names <- gsub("@",".",W.names)
  X.names <- gsub("@",".",X.names)
  
  W <- lapply(1:dim(spec)[1], function(i){
    h <- as.data.frame(data.combs[[i]][,W.names])
    names(h) <- W.names
    h
  })
  
  X <- lapply(1:length(data.combs), function(i){ 
    if(i <= dim(spec)[1]){
      h <- as.data.frame(data.combs[[i]][1:2,X.names])
      names(h) <- X.names
      h
    } else{
      h <- as.data.frame(data.combs[[i]][,X.names])
      names(h) <- X.names
      h
    }
  })
  
  return(list(D=D, R=R, W=W, X=X, V=V, P=P, epsilon=epsilon, eta=eta, xi=xi, combs=combs, E=E))
  
}
