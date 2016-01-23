# ----------------------------------------------------------------------------
# R-code (www.r-project.org/) for the Structural Matching Model
#
# Copyright (c) 2013 Thilo Klein
#
# This library is distributed under the terms of the GNU Public License (GPL)
# for full details see the file LICENSE
#
# ----------------------------------------------------------------------------

#' @title Structural Matching Model to correct for sample selection bias in two-sided matching markets
#'
#' @description The function provides a Gibbs sampler for a structural matching model that corrects 
#' for sample selection bias when the selection process is a two-sided matching game; i.e., 
#' a matching of students to colleges.
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
#' that satisfy the equilibrium condition (see Sorensen, 2007). 
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
#' @param OUT data frame with characteristics of all observed matches, including
#' market identifier \code{m.id}, college identifier \code{c.id} and student identifier \code{s.id}.
#' @param SEL optional: data frame with characteristics of all observed and unobserved matches, including 
#' market identifier \code{m.id}, college identifier \code{c.id} and student identifier \code{s.id}.
#' @param colleges character vector of variable names for college characteristics. These variables carry the same value for any college.
#' @param students character vector of variable names for student characteristics. These variables carry the same value for any student.
#' @param m.id character name of the market identifier variable. Defaults to \code{"m.id"}.
#' @param c.id character name of the college identifier variable. Defaults to \code{"c.id"}.
#' @param s.id character name of the student identifier variable. Defaults to \code{"s.id"}.
#' @param outcome formula for match outcomes.
#' @param selection formula for match valuations.
#' @param selection.college formula for match valuations of colleges. This argument is ignored when \code{selection} is provided.
#' @param selection.student formula for match valuations of students. This argument is ignored when \code{selection} is provided.
#' @param binary logical: if \code{TRUE} outcome variable is taken to be binary; if \code{FALSE} outcome variable is taken to be continuous.
#' @param niter number of iterations to use for the Gibbs sampler.
#' @param gPrior logical: if \code{TRUE} the g-prior (Zellner, 1986) is used for the variance-covariance matrix. (Not yet implemented)
#' @param censored draws of the \code{delta} parameter that estimates the covariation between the error terms in selection and outcome equation are 0:not censored, 1:censored from below, 2:censored from above.
#' @param seed integer setting the state for MCMC draws.
#' 
#' @export
#' 
#' @useDynLib matchingMarkets
#' 
#' @import partitions lpSolve stats
#' @importFrom Rcpp evalCpp
#' 
#' @aliases stabitCpp2
#' 
#' @return
#' \code{stabit} returns a list with the following items.
#' \item{model.list}{}
#' \item{draws}{}
#' \item{coefs}{}
#' 
#' @author Thilo Klein 
#' 
#' @keywords regression
#' 
#' @references Sorensen, M. (2007). How Smart is Smart Money? A Two-Sided Matching Model of Venture Capital.
#' \emph{Journal of Finance}, 62 (6): 2725-2762.
#' 
#' @examples
#' ## --- SIMULATED EXAMPLE ---
#' \dontrun{
#' ## 1. Simulate two-sided matching data for 20 markets (m=20) with 100 students
#' ##    (nStudents=100) per market and 20 colleges with quotas of 5 students, each
#' ##    (nSlots=rep(5,20)).
#' 
#' xdata <- stabsim2(m=20, nStudents=100, nSlots=rep(5,20), 
#'   colleges = "c1",
#'   students = "s1",
#'   outcome = ~ c1:s1 + eta +nu,
#'   selection = ~ -1 + c1:s1 + eta
#' )
#' head(xdata$OUT)
#' 
#' 
#' ## 2-a. Bias from sorting
#'  lm1 <- lm(y ~ c1:s1, data=xdata$OUT)
#'  summary(lm1)
#' 
#' ## 2-b. Cause of the bias
#'  with(xdata$OUT, cor(c1*s1, eta))
#' 
#' ## 2-c. Correction for sorting bias
#'  lm2a <- lm(V ~ -1 + c1:s1, data=xdata$SEL); summary(lm2a)
#'  etahat <- lm2a$residuals[xdata$SEL$D==1]
#'  
#'  lm2b <- lm(y ~ c1:s1 + etahat, data=xdata$OUT)
#'  summary(lm2b)
#' 
#' 
#' ## 3. Correction for sorting bias when match valuation V is unobserved
#' 
#' ## 3-a. Run Gibbs sampler (when SEL is given)
#'  fit2 <- stabit2(OUT = xdata$OUT, 
#'            SEL = xdata$SEL,
#'            outcome = y ~ c1:s1, 
#'            selection = ~ -1 + c1:s1,
#'            niter=1000
#'  )
#'
#' ## 3-b. Run Gibbs sampler (when SEL is not given)
#'  fit2 <- stabit2(OUT = xdata$OUT, 
#'            colleges = "c1",
#'            students = "s1",
#'            outcome = y ~ c1:s1, 
#'            selection = ~ -1 + c1:s1,
#'            niter=1000
#'  )
#'
#' ## 4-a. Get marginal effects (for linear model)
#'  fit2$coefs
#'  
#' ## 4-b. Get marginal effects (for probit)
#'  mfx(fit2)
#'  
#'  
#' ## 5. Plot MCMC draws for coefficients
#'  plot(fit2$draws$alphadraws[1,], type="l")
#'  plot(fit2$draws$betadraws[1,], type="l")
#'  plot(fit2$draws$deltadraws[1,], type="l")
#'  
#'  
#' ## 6. Obtain the model list used in estimation
#'  head(fit2$model.list)
#' }
stabit2 <- function(OUT, SEL=NULL,
                    colleges=NULL, students=NULL, m.id="m.id", c.id="c.id", s.id="s.id", outcome, selection=NULL,
                    selection.student=NULL, selection.college=NULL, binary=FALSE, niter, gPrior=FALSE, 
                    censored=1, seed=123){
  
  ## 1) creates an edgelist for all feasible matches taken as given an edgelist of either 
  ##    (i) equilibrium matches [as in daa()] or
  ##    (ii) equilibrium plus other matches [as in school choice data]
  ## 2) in edgelist of all feasible matches, bring equilibrim matches to top and split into two datasets
  ##    one for Xall, the other for Xmatch
  
  if(is.null(selection)){
    method <- "Klein"
  } else{
    method <- "Sorensen"
  }
  
  if(!is.null(SEL)){
    SEL <- split(SEL, SEL$m.id) 
  }
    
  #OUT = xdata$OUT 
  ##SEL = xdata$SEL
  #colleges = "c1"
  #students = "s1"
  #outcome = y ~ c1:s1 
  #selection = ~ -1 + c1:s1
  #niter=100

  
  Y=list(); Xmatch=list(); C=list(); Cmatch=list(); S=list(); Smatch=list(); D=list(); d=list(); M=list(); H=list()
  
  for(i in 1:length(unique(OUT[,m.id]))){
    
    X <- stabit2_inner(iter=i, data=OUT, SEL=SEL,
                       colleges=colleges, students=students, outcome=outcome, selection=selection,
                       selection.student=selection.student, selection.college=selection.college, method=method)

    Y[[i]]      <- as.matrix(X$Y)
    Xmatch[[i]] <- as.matrix(X$Xmatch)
    C[[i]]      <- as.matrix(X$C)
    Cmatch[[i]] <- as.matrix(X$Cmatch)
    D[[i]]      <- as.matrix(X$D)
    d[[i]]      <- as.matrix(X$d) - 1
    M[[i]]      <- as.matrix(X$M) - 1
    H[[i]]      <- as.matrix(X$H)
    if(method=="Klein"){
      S[[i]]      <- as.matrix(X$S)
      Smatch[[i]] <- as.matrix(X$Smatch)
    }
    
  }

  ## Preliminaries.
  T <- length(Y); #// Number of markets.
  nColleges <- nStudents <- XXmatch <- CC <- SS <- CCmatch <- SSmatch <- list()
  L <- studentIds <- collegeId <- rep(list(vector()),T)
  
  s <- 0
  for(i in 1:T){
    nColleges[[i]] <- nrow(H[[i]])
    nStudents[[i]] <- ncol(H[[i]])
    
    XXmatch[[i]]   <- t(Xmatch[[i]]) %*% Xmatch[[i]]
    CC[[i]]        <- t(C[[i]]) %*% C[[i]]
    CCmatch[[i]]   <- t(Cmatch[[i]]) %*% Cmatch[[i]]
    if(method=="Klein"){
      SSmatch[[i]]   <- t(Smatch[[i]]) %*% Smatch[[i]]
      SS[[i]]        <- t(S[[i]]) %*% S[[i]]
    }
    
    ## Record the id's of students matched to each college, and the id of the college matched to each student.
    for(j in 1:nColleges[[i]]){
      s <- s+1
      L[[i]][j] <- s - 1
      studentIds[[s]] <- which(H[[i]][j,] == 1) - 1
    }
    
    for(j in 1:nStudents[[i]]){
      collegeId[[i]][j] <- which(H[[i]][,j] == 1) - 1
    }
  }
  n <- sum(unlist(nStudents)) ## Total number of students/matches.
  N <- sum(unlist(nColleges)) ## Total number of students/matches.
  
  if(method=="Klein"){
    return( list(Y=Y, Xmatch=Xmatch, C=C, Cmatch=Cmatch, S=S, Smatch=Smatch, D=D, d=d, M=M, H=H, 
                 nColleges=unlist(nColleges), nStudents=unlist(nStudents), XXmatch=XXmatch, CC=CC, SS=SS, 
                 CCmatch=CCmatch, SSmatch=SSmatch, L=L, studentIds=studentIds, collegeId=collegeId, n=n, N=N,
                 binary=binary, niter=niter, T=T, gPrior=gPrior) )    
  } else if(method=="Sorensen"){
    fit2 <- list(Y=Y, Xmatch=Xmatch, C=C, Cmatch=Cmatch, D=D, d=d, M=M, H=H, 
                 nColleges=unlist(nColleges), nStudents=unlist(nStudents), XXmatch=XXmatch, CC=CC, 
                 CCmatch=CCmatch, L=L, studentIds=studentIds, collegeId=collegeId, n=n, N=N,
                 binary=binary, niter=niter, T=T, gPrior=gPrior)
    
    # -----------------------------------------------------------------------------
    # Source C++ script
    # -----------------------------------------------------------------------------    
    #sourceCpp("/home/thilo/Documents/Research/Matching/matchingMarkets2/src/stabitCpp2sorensen.cpp")
    
    res <- with(fit2, stabitCpp2(Yr=Y, Xmatchr=Xmatch, Cr=C, Cmatchr=Cmatch, 
                                 #Sr=S, Smatchr=Smatch, 
                                 Dr=D, dr=d,
                                 Mr=M, Hr=H, nCollegesr=nColleges, nStudentsr=nStudents, XXmatchr=XXmatch, CCr=CC, 
                                 #SSr=SS, SSmatchr=SSmatch,
                                 CCmatchr=CCmatch, Lr=L, studentIdsr=studentIds, collegeIdr=collegeId, n=n, 
                                 #Vcr=Vc, Vsr=Vs,
                                 N=N, binary=binary, niter=niter, T=T, censored=censored))
    
    # -----------------------------------------------------------------------------
    # Add names to coefficients.
    # ----------------------------------------------------------------------------- 
    
    an <- names(fit2$Cmatch)
    bn <- names(fit2$Xmatch)
    
    # parameter draws
    rownames(res$alphadraws) = an
    rownames(res$betadraws) = bn
    rownames(res$deltadraws) = "delta"
    # posterior means
    rownames(res$alpha) = an
    rownames(res$beta) = bn
    rownames(res$delta) = "delta"
    colnames(res$eta) = "eta"
    # vcov
    rownames(res$alphavcov) = colnames(res$alphavcov) = an
    rownames(res$betavcov) = colnames(res$betavcov) = bn
    #
    colnames(res$alpha) = colnames(res$beta) = colnames(res$delta)  = colnames(res$sigmasquarexi) = c("coef","s.e.")
    #
    if(binary==TRUE){
      rownames(res$sigmasquarexi) = "sigma"
      out <- list(draws=with(res,list(alphadraws=alphadraws,betadraws=betadraws,deltadraws=deltadraws)), 
                  coefs=with(res,list(eta=eta,alphavcov=alphavcov,betavcov=betavcov,alpha=alpha,beta=beta,
                                      delta=delta,sigmasquarexi=sigmasquarexi)))    
    } else{
      out <- list(draws=with(res,list(alphadraws=alphadraws,betadraws=betadraws,deltadraws=deltadraws)), 
                  coefs=with(res,list(eta=eta,alphavcov=alphavcov,betavcov=betavcov,alpha=alpha,beta=beta,
                                      delta=delta)))
    }  
    
    # -----------------------------------------------------------------------------
    # Returns .
    # ----------------------------------------------------------------------------- 
    #model.frame <- do.call(rbind, data)
    #return(list(model.list=data, model.frame=model.frame, draws=out$draws, coefs=out$coefs))
    
    if(method=="Sorensen"){
      
      model.list <- list(D=D, Y=Y, W=C, X=Xmatch)
      return(list(model.list=model.list, draws=out$draws, coefs=out$coefs))
      
    } else if(method=="Klein"){
      
      return(list(draws=out$draws, coefs=out$coefs))
      
    }
    
  }  
}


rFormula <- function(formula, data=list(), ...){
  mf <- model.frame(formula=formula, data=data)
  x <- model.matrix(attr(mf, "terms"), data=mf)
  y <- model.response(mf)
  as.data.frame(cbind(y,x))
}


stabit2_inner <- function(iter, data, SEL, colleges, students, m.id="m.id", c.id="c.id", s.id="s.id", outcome, selection, 
                          selection.student, selection.college, method){
  
  data <- data[data[,m.id]==iter,]
  
  if(is.null(SEL)){
    
    ## 
    c.id <- data[,c.id]
    s.id <- data[,s.id]
    
    ## unique student and college ids
    uColleges <- sort(unique(c.id))
    uStudents <- sort(unique(s.id))
    
    ## all feasible combinations
    combs <- function(uColleges, uStudents){
      nColleges <- length(uColleges)
      nStudents <- length(uStudents)
      data.frame( c.id = c(sapply(uColleges, function(i){ rep(i, nStudents) })), 
                  s.id = rep(uStudents, nColleges), 
                  stringsAsFactors=FALSE )
    }
    indices <- as.data.frame(combs(uColleges, uStudents))
    
    ## index equilibrium groups from x in indices
    indices$id <- paste(indices$c.id, indices$s.id, sep="_")
    data$id    <- paste(data$c.id, data$s.id, sep="_")
    
    ## sort indices such that observed matches come first
    indices$D <- ifelse(indices$id %in% data$id, 1, 0)
    indices <- indices[order(indices$D, decreasing=TRUE),]
    
    ## sort data such that it concordes with indices
    data    <- data[match(indices$id[indices$D==1], data$id), ] 
    
  } else{
    
    SEL <- SEL[[iter]]
    indices <- SEL
    
    ## unique student and college ids
    uColleges <- sort(unique(indices$c.id))
    uStudents <- sort(unique(indices$s.id))
    
  }
  
  ## --- look-up matrix of dim nColleges x nStudents ---
  ind.new       <- indices
  ind.new$index <- 1:nrow(indices)
  ind.org <- ind.new[order(ind.new$c.id,ind.new$s.id),]
  M <- matrix(ind.org$index, nrow=length(uColleges), ncol=length(uStudents), byrow=TRUE)
  H <- matrix(ind.org$D, nrow=length(uColleges), ncol=length(uStudents), byrow=TRUE)  
  
  ## Main effects...
  #all.vars(outcome)
  
  if(is.null(SEL)){
    
    c.data <- data[!duplicated(data$c.id),]
    C <- data.frame(apply(data.frame(c.data[,colleges]), 2, function(i) i[match(indices$c.id,c.data$c.id)]))
    names(C) <- colleges
    
    s.data <- data[!duplicated(data$s.id),]
    S <- data.frame(apply(data.frame(s.data[,students]), 2, function(i) i[match(indices$s.id,s.data$s.id)]))
    names(S) <- students
    
    ## ... and interaction effects
    Xmain <- data.frame(C, S)
    Xmatch <- rFormula(formula = outcome, data=data)
    
    if(method=="Klein"){
      
      C <- rFormula(formula = selection.college, data=Xmain)
      S <- rFormula(formula = selection.student, data=Xmain)    
      
    } else if(method == "Sorensen"){
      
      C <- rFormula(formula = selection, data=Xmain)
      
    }
    
  } else{
    
    Xmatch <- rFormula(formula = outcome, data=data)
    
    if(method=="Klein"){
      
      C <- rFormula(formula = selection.college, data=SEL) ## need to define SELc !!!
      S <- rFormula(formula = selection.student, data=SEL) ## need to define SELs !!!   
      
    } else if(method == "Sorensen"){
      
      C <- rFormula(formula = selection, data=SEL)
      
    }
    
    
  }
  


  D <- indices$D
  d <- which(D==1)
  
  ## --- OUTPUT for stabit2() ---
  
  if(method=="Klein"){
    
    return( list(Y=Xmatch[,1], Xmatch=Xmatch[,-1], C=C, Cmatch=C[D==1,], S=S, Smatch=S[D==1,], D=D, d=d, M=M, H=H) )
    
  } else if(method=="Sorensen"){
    
    return( list(Y=Xmatch[,1], Xmatch=Xmatch[,-1], C=C, Cmatch=C[D==1,], D=D, d=d, M=M, H=H) )
  }  
}  



