# ----------------------------------------------------------------------------
# R-code (www.r-project.org/) for the Structural Matching Model
#
# Copyright (c) 2016 Thilo Klein
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
#'        V &= & W\beta + \eta
#'        \end{array}
#'      }{ D = 1[V in \Gamma] with V = W\beta + \eta
#'      }
#' Here, \eqn{V} is a vector of latent valuations of \emph{all feasible} matches, ie observed and 
#' unobserved, and \eqn{1[.]} is the Iverson bracket. 
#' A match is observed if its match valuation is in the set of valuations \eqn{\Gamma}
#' that satisfy the equilibrium condition (see Sorensen, 2007). 
#' The match valuation \eqn{V} is a linear function of \eqn{W}, a matrix of characteristics for 
#' \emph{all feasible} matches, and \eqn{\eta}, a vector of random errors. \eqn{\beta} is a paramter 
#' vector to be estimated.
#' 
#' The \emph{Outcome Equation} determines the outcome for \emph{observed} matches. The dependent
#' variable can either be continuous or binary, dependent on the value of the \code{binary}
#' argument. In the binary case, the dependent variable \eqn{R} is determined by a threshold 
#' rule for the latent variable \eqn{Y}.
#' \deqn{ \begin{array}{lcl}
#'        R &= & 1[Y > c] \\
#'        Y &= & X\alpha + \epsilon
#'        \end{array}
#'      }{ R = 1[Y > c] with Y = X\alpha + \epsilon
#'      }
#' Here, \eqn{Y} is a linear function of \eqn{X}, a matrix of characteristics for \emph{observed} 
#' matches, and \eqn{\epsilon}, a vector of random errors. \eqn{\alpha} is a paramter vector to 
#' be estimated.
#' 
#' The structural model imposes a linear relationship between the error terms of both equations 
#' as \eqn{\epsilon = \kappa\eta + \nu}, where \eqn{\nu} is a vector of random errors and \eqn{\kappa}
#' is the covariance paramter to be estimated. If \eqn{\kappa} were zero, the marginal distributions
#' of \eqn{\epsilon} and \eqn{\eta} would be independent and the selection problem would vanish.
#' That is, the observed outcomes would be a random sample from the population of interest.
#' 
#' @param OUT data frame with characteristics of all observed matches, including
#' market identifier \code{m.id}, college identifier \code{c.id} and student identifier \code{s.id}.
#' @param SEL optional: data frame with characteristics of all observed and unobserved matches, including 
#' market identifier \code{m.id}, college identifier \code{c.id} and student identifier \code{s.id}.
#' @param colleges character vector of variable names for college characteristics. These variables carry the same value for any college.
#' @param students character vector of variable names for student characteristics. These variables carry the same value for any student.
#' @param outcome formula for match outcomes.
#' @param selection formula for match valuations. 
#' @param binary logical: if \code{TRUE} outcome variable is taken to be binary; if \code{FALSE} outcome variable is taken to be continuous.
#' @param niter number of iterations to use for the Gibbs sampler.
#' @param gPrior logical: if \code{TRUE} the g-prior (Zellner, 1986) is used for the variance-covariance matrix. (Not yet implemented)
#' @param censored draws of the \code{kappa} parameter that estimates the covariation between the error terms in selection and outcome equation are 0:not censored, 1:censored from below, 2:censored from above.
#' @param thin integer indicating the level of thinning in the MCMC draws. The default \code{thin=1} saves every draw, \code{thin=2} every second, etc.
#' @param ... .
#' 
# @param selection.college formula for match valuations of colleges. Is ignored when \code{selection} is provided.
# @param selection.student formula for match valuations of students. Is ignored when \code{selection} is provided.
# @param s.prefs list of matrices (one element for each market) of dimension \code{nColleges} \code{x} \code{nStudents} with the \code{j}th 
# column containing student \code{j}'s ranking over colleges in decreasing order of 
# preference (i.e. most preferred first).
# @param c.prefs list of matrices (one element for each market) of dimension \code{nStudents} \code{x} \code{nColleges} with the \code{i}th 
# column containing college \code{i}'s ranking over students in decreasing order of 
# preference (i.e. most preferred first).
# @param SELs optional: same as \code{SEL} but for student valuation when estimating separate selection equations for college and student preferences. Is ignored when \code{SEL} is provided.
# @param SELc optional: same as \code{SEL} but for college valuation when estimating separate selection equations for college and student preferences. Is ignored when \code{SEL} is provided.
#' 
#' @export
#' 
#' @useDynLib matchingMarkets
#' 
#' @import partitions stats
#' @importFrom Rcpp evalCpp
#' 
#' @aliases stabitCpp2 stabitCpp3
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
#' ##    (nSlots=rep(5,20)). True parameters in selection and outcome equations are 
#' ##    all equal to 1.
#' 
#' xdata <- stabsim2(m=20, nStudents=100, nSlots=rep(5,20), 
#'   colleges = "c1",
#'   students = "s1",
#'   outcome = ~ c1:s1 + eta + nu,
#'   selection = ~ -1 + c1:s1 + eta
#' )
#' head(xdata$OUT)
#' 
#' 
#' ## 2. Correction for sorting bias when match valuations V are observed
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
#' ## 3. Correction for sorting bias when match valuations V are unobserved
#' 
#' ## 3-a. Run Gibbs sampler (when SEL is given)
#'  fit2 <- stabit2(OUT = xdata$OUT, 
#'            SEL = xdata$SEL,
#'            outcome = y ~ c1:s1, 
#'            selection = ~ -1 + c1:s1,
#'            niter=1000
#'  )
#'
#' ## 3-b. Alternatively: Run Gibbs sampler (when SEL is not given)
#'  fit2 <- stabit2(OUT = xdata$OUT, 
#'            colleges = "c1",
#'            students = "s1",
#'            outcome = y ~ c1:s1, 
#'            selection = ~ -1 + c1:s1,
#'            niter=1000
#'  )
#'
#'
#' ## 4. Implemented methods
#'
#' ## 4-a. Get coefficients
#'  fit2
#'  
#' ## 4-b. Coefficient table
#'  summary(fit2)
#'  
#' ## 4-c. Get marginal effects
#'  summary(fit2, mfx=TRUE)
#'  
#' ## 4-d. Also try the following functions
#'  coef(fit2)
#'  fitted(fit2)
#'  residuals(fit2)
#'  predict(fit, newdata=NULL)
#'
#'    
#' ## 5. Plot MCMC draws for coefficients in outcome equation
#'  res <- as.data.frame(t(fit2$draws$alphadraws))
#'  res$iteration <- 1:nrow(res)
#'  library(tidyr)
#'  res.long <- gather(res, condition, measurement, 1:(ncol(res)-1))
#'  
#'  library(lattice)
#'  lattice.options(default.args=list(as.table=TRUE), 
#'                  default.theme=standard.theme(color=FALSE))
#'  xyplot(measurement ~ iteration | factor(condition), 
#'         data = res.long, scales=list(relation="free"),
#'         xlab = "iterations",
#'         ylab = "paramter draws", type = "l") 
#' }
stabit2 <- function(OUT, SEL=NULL, colleges=NULL, students=NULL, outcome, selection,
                    binary=FALSE, niter, gPrior=FALSE, 
                    censored=1, thin=1, ...) UseMethod("stabit2")

#' @export
print.stabit2 <- function(x, ...){
  
  cat("Call:\n")
  print(x$call)
  cat("\nCoefficients:\n")
  print(x$coefficients)
}

#' @export
stabit2.default <- function(OUT, SEL=NULL, colleges=NULL, students=NULL, outcome, selection,
                            binary=FALSE, niter, gPrior=FALSE, 
                            censored=1, thin=1, ...){
  
  ## ------------------------
  ## --- 1. Preliminaries ---
  
  ## select method based on arguments provided
  if(is.list(selection)){
    method <- "Klein" # separate selection equations for students and college preferences
    selection.college <- selection$college
    selection.student <- selection$student
    selection <- NULL
  } else{
    method <- "Sorensen" # single selection equation with equal sharing rule for student and college utility
  }
  
  if(!is.data.frame(SEL)){
    SELs <- SEL$SELs
    SELc <- SEL$SELc
    s.prefs <- SEL$s.prefs
    c.prefs <- SEL$c.prefs
    SEL <- NULL
  }
  
  ## split input data by market
  OUT <- split(OUT, OUT$m.id)
  if(!is.null(SEL)){
    SEL <- split(SEL, SEL$m.id) 
    SELs <- NULL
    SELc <- NULL
  } else if(!is.null(SELs) & !is.null(SELc)){
    SELs <- split(SELs, SELs$m.id) 
    SELc <- split(SELc, SELc$m.id) 
  }
  
  ## market and agent identifiers
  m.id <- "m.id"
  c.id <- "c.id"
  s.id <- "s.id"
  
  ## empty lists for results
  Y=list(); Xmatch=list(); C=list(); Cmatch=list(); S=list(); Smatch=list(); D=list(); d=list()
  M=list(); H=list(); H2=list(); copt.id=vector(); indices=list()
  
  for(i in 1:length(OUT)){
    
    ## ------------------------------------------------------------------------------------------
    ## --- 2. Bring data in consistent order and produce matrices for all equilibrium matches ---
    
    X <- stabit2_inner(iter=i, OUT=OUT[[i]], SEL=SEL[[i]], SELs=SELs[[i]], SELc=SELc[[i]],
                       colleges=colleges, students=students, outcome=outcome, selection=selection,
                       selection.student=selection.student, selection.college=selection.college, 
                       s.prefs=s.prefs, c.prefs=c.prefs,
                       method=method)
    
    ## continue to write the results per market
    Y[[i]]      <- as.matrix(X$Y)
    Xmatch[[i]] <- as.matrix(X$Xmatch)
    C[[i]]      <- as.matrix(X$C)
    Cmatch[[i]] <- as.matrix(X$Cmatch)
    D[[i]]      <- as.matrix(X$D)
    d[[i]]      <- as.matrix(X$d) - 1
    M[[i]]      <- as.matrix(X$M) - 1
    indices[[i]] <- X$indices
    if(method=="Klein"){
      S[[i]]      <- as.matrix(X$S)
      Smatch[[i]] <- as.matrix(X$Smatch)
      H[[i]]      <- X$H
      H2[[i]]     <- X$H2
      copt.id[i]  <- X$copt.id
    } else{
      H[[i]]      <- as.matrix(X$H)
    }
  }
  
  ## some preliminaries
  T <- length(Y); #// Number of markets.
  nColleges <- nStudents <- XXmatch <- CC <- SS <- CCmatch <- SSmatch <- list()
  L <- rep(list(vector()),T)
  if(method=="Klein"){
    studentIds <- collegeId <- rep(list(matrix()),T)
  } else{
    studentIds <- collegeId <- rep(list(vector()),T)
  }
  
  ## further results
  s <- 0
  for(i in 1:T){
    if(method=="Klein"){
      nColleges[[i]] <- nrow(H[[i]][,,1])
      nStudents[[i]] <- ncol(H[[i]][,,1])
    } else{
      nColleges[[i]] <- nrow(H[[i]])
      nStudents[[i]] <- ncol(H[[i]])
    }
    
    XXmatch[[i]]   <- t(Xmatch[[i]]) %*% Xmatch[[i]]
    CC[[i]]        <- t(C[[i]]) %*% C[[i]]
    CCmatch[[i]]   <- t(Cmatch[[i]]) %*% Cmatch[[i]]
    if(method=="Klein"){
      SSmatch[[i]]   <- t(Smatch[[i]]) %*% Smatch[[i]]
      SS[[i]]        <- t(S[[i]]) %*% S[[i]]
      collegeId[[i]] <- matrix(NA, nrow=nStudents[[i]], ncol=dim(H[[i]])[3]) ##
    }
    
    ## Record the id's of students matched to each college, and the id of the college matched to each student.
    for(j in 1:nColleges[[i]]){
      s <- s+1
      L[[i]][j] <- s - 1
      if(method=="Klein"){
        studentIds[[s]] <- matrix(NA, nrow=length(which(H[[i]][j,,1] == 1)), ncol=dim(H[[i]])[3])
        for(k in 1:dim(H[[i]])[3]){
          studentIds[[s]][,k] <- which(H[[i]][j,,k] == 1) - 1
        }
      } else{
        studentIds[[s]] <- which(H[[i]][j,] == 1) - 1
      }      
    }
    
    for(j in 1:nStudents[[i]]){
      if(method=="Klein"){
        for(k in 1:dim(H[[i]])[3]){
          collegeId[[i]][j,k] <- which(H[[i]][,j,k] == 1) - 1
        }
      } else{
        collegeId[[i]][j] <- which(H[[i]][,j] == 1) - 1
      }      
    }
  }
  n <- sum(unlist(nStudents)) # total number of students/matches
  N <- sum(unlist(nColleges)) # total number of colleges
  nEquilibs <- unlist(lapply(H, function(z) dim(z)[3])) # number of equilibria
  
  ## ---------------------------------------------------------------------------------------------
  ## --- 3. Mapping of hospital/residents problem (HR) to related stable marriage problem (SM) ---
  
  if(method=="Klein"){
    
    sopt.id <- 1 ## student-optimal matching is first in lists
    
    sopt2equ <- list()
    for(i in 1:length(H2)){ ## for each market i
      h <- apply(H2[[i]][,,sopt.id], 2, function(z) which(z==1))
      sopt2equ[[i]] <- -1 + sapply(1:dim(H2[[i]])[3], function(z){
        sapply(1:length(h), function(j){
          which(H2[[i]][h[j],,z]==1)
        })
      }) 
    }
    
    copt2equ <- list() ##!!!
    for(i in 1:length(H2)){
      h <- apply(H2[[i]][,,copt.id[i]], 2, function(z) which(z==1))
      copt2equ[[i]] <- -1 + sapply(1:dim(H2[[i]])[3], function(z){
        sapply(1:length(h), function(j){
          which(H2[[i]][h[j],,z]==1)
        })
      }) 
    }
    
    equ2sopt <- list()
    for(i in 1:length(H2)){
      h <- sapply(1:dim(H2[[i]])[3], function(z){
        apply(H2[[i]][,,z], 2, function(j){
          which(j==1)
        })
      })
      equ2sopt[[i]] <- -1 + sapply(1:dim(H2[[i]])[3], function(z){
        sapply(1:nrow(h), function(j){
          which(H2[[i]][h[j,z],,sopt.id]==1)
        })
      }) 
    }
    
    equ2copt <- list() ##!!!
    for(i in 1:length(H2)){
      h <- sapply(1:dim(H2[[i]])[3], function(z){
        apply(H2[[i]][,,z], 2, function(j){
          which(j==1)
        })
      })
      equ2copt[[i]] <- -1 + sapply(1:dim(H2[[i]])[3], function(z){
        sapply(1:nrow(h), function(j){
          which(H2[[i]][h[j,z],,copt.id[i]]==1)
        })
      }) 
    }
  }
  
  ## adjust index for college-optimal matching for C++
  copt.id <- copt.id -1
  
  ## ----------------------------
  ## --- 4. Run Gibbs sampler ---
  
  if(method=="Klein"){
    
    est <- stabitCpp3(Yr=Y, Xmatchr=Xmatch, Cr=C, Cmatchr=Cmatch, Sr=S, Smatchr=Smatch, Dr=D, dr=d,
                      Mr=M, Hr=H, nCollegesr=unlist(nColleges), nStudentsr=unlist(nStudents), XXmatchr=XXmatch, 
                      CCr=CC, SSr=SS, SSmatchr=SSmatch, CCmatchr=CCmatch, Lr=L,
                      studentIdsr=studentIds, collegeIdr=collegeId, nEquilibsr=nEquilibs,
                      equ2soptr=equ2sopt, sopt2equr=sopt2equ, equ2coptr=equ2copt, copt2equr=copt2equ,
                      coptidr=copt.id, n=n, N=N, binary=binary, niter=niter, thin=thin, T=T, censored=censored)
    
  } else if(method=="Sorensen"){
    
    est <- stabitCpp2(Yr=Y, Xmatchr=Xmatch, Cr=C, Cmatchr=Cmatch, Dr=D, dr=d,
                      Mr=M, Hr=H, nCollegesr=unlist(nColleges), nStudentsr=unlist(nStudents), XXmatchr=XXmatch, 
                      CCr=CC, CCmatchr=CCmatch, Lr=L, 
                      studentIdsr=studentIds, collegeIdr=collegeId, 
                      n=n, N=N, binary=binary, niter=niter, thin=thin, T=T, censored=censored)
  }
  
  # ------------------------------------
  # --- 5. Add names to coefficients ---
  
  ## variable names
  
  an <- colnames(Xmatch[[1]])
  bn <- colnames(C[[1]])
  if(method=="Klein"){
    cn <- colnames(S[[1]])
  }
  
  ## parameter draws
  
  rownames(est$betadraws) <- bn
  if(method=="Sorensen"){
    est$alphadraws <- with(est, rbind(alphadraws, kappadraws))
    rownames(est$alphadraws) <- c(an,"kappa")
    est$kappadraws <- NULL
  } else if(method=="Klein"){
    est$alphadraws <- with(est, rbind(alphadraws, kappadraws, lambdadraws))
    rownames(est$alphadraws) <- c(an,"kappa","lambda")
    est$kappadraws  <- NULL
    est$lambdadraws <- NULL
    rownames(est$gammadraws) = cn
  }
  
  ## posterior means
  
  ## The last half of all draws are used in approximating the posterior means and the posterior standard deviations.
  niter <- dim(est$alphadraws)[2]
  startiter <- floor(niter/2)
  
  posMeans <- function(x){
    sapply(1:nrow(x), function(z) mean(x[z,startiter:niter]))
  }
  
  est$alpha <- posMeans(est$alphadraws) 
  names(est$alpha) <- rownames(est$alphadraws)
  
  est$beta <- posMeans(est$betadraws)
  names(est$beta) <- bn
  
  if(method=="Klein"){
    est$gamma <- posMeans(est$gammadraws)
    names(est$gamma) <- cn
  }
  
  if(binary==FALSE){
    est$sigma <- sqrt(posMeans(est$sigmasquarenudraws))
    est$sigmasquarenudraws <- NULL
  }
  
  ## error terms
  
  est$eta <- posMeans(est$etadraws)
  est$etadraws <- NULL
  if(method=="Klein"){
    est$delta <- posMeans(est$deltadraws)
    est$deltadraws <- NULL
  }
  
  ## vcov
  
  est$vcov$alpha <- var(t(est$alphadraws))
  est$vcov$beta <- var(t(est$betadraws))
  if(method=="Klein"){
    est$vcov$gamma <- var(t(est$gammadraws))
  }
  
  ## variables
  
  est$X <- do.call(rbind, Xmatch)
  est$X <- cbind(est$X, eta=est$eta)
  if(method=="Klein"){
    est$X <- cbind(est$X, delta=est$delta)
    est$S <- do.call(rbind, S)  
  }
  est$C <- do.call(rbind, C)
  est$Y <- do.call(c, Y)
  est$D <- do.call(c, D)
  
  if(method=="Sorensen"){
    
    est$coefficients <- unlist(with(est, list(o=alpha, s=beta)))
    
  } else if(method=="Klein"){
    
    est$coefficients <- unlist(with(est, list(o=alpha, c=beta, s=gamma)))
  }
  est$fitted.values <- with(est, as.vector(X %*% alpha))
  est$residuals <- est$Y - est$fitted.values
  est$df <- n-ncol(est$X)
  
  est$call <- match.call()
  est$method <- method
  est$binary <- binary
  est$formula <- outcome
  
  ## consolidate
  
  h <- est[grep(pattern="draws",x=names(est))]
  est[grep(pattern="draws",x=names(est))] <- NULL
  est$draws <- h
  
  h <- est[which(names(est) %in% c("alpha","beta","gamma","kappa","lambda"))]
  est[which(names(est) %in% c("alpha","beta","gamma","kappa","lambda"))] <- NULL
  est$coefs <- h
  
  h <- est[which(names(est) %in% c("X","S","C","Y","D"))]
  est[which(names(est) %in% c("X","S","C","Y","D"))] <- NULL
  est$variables <- h
  
  h <- NULL
  
  # ------------------
  # --- 6. Returns ---
  
  class(est) <- "stabit2"
  return(est)
}




rFormula <- function(formula, data=list(), ...){
  mf <- model.frame(formula=formula, data=data)
  x <- model.matrix(attr(mf, "terms"), data=mf)
  y <- model.response(mf)
  as.data.frame(cbind(y,x))
}




stabit2_inner <- function(iter, OUT, SEL, SELs, SELc, colleges, students, 
                          m.id="m.id", c.id="c.id", s.id="s.id", 
                          outcome, selection, selection.student, selection.college, 
                          s.prefs=s.prefs, c.prefs=c.prefs, method){
  
  ## ----------------------------------------------------------------------------------------------
  ## --- 1. Sort datasets by colleges (c.id) and students (s.id), put equilibrium matches first ---
  
  if(!is.null(SEL)){
    
    ## sort SEL by c.id and s.id (based on OUT)
    idOUT <- with(OUT, paste(c.id, s.id, sep="_"))
    SEL$idSEL <- with(SEL, paste(c.id, s.id, sep="_"))
    
    ## which matches in SEL are in OUT?
    SEL_top <- SEL[SEL$idSEL %in% idOUT,]
    SEL_top <- SEL_top[match(idOUT, SEL_top$idSEL),]
    SEL_bot <- SEL[! SEL$idSEL %in% idOUT,]
    SEL <- rbind(SEL_top, SEL_bot)
    rm(SEL_top, SEL_bot)
    
    ## unique student and college ids
    indices <- SEL
    uColleges <- sort(unique(indices$c.id))
    uStudents <- sort(unique(indices$s.id))
    
  } else if(!is.null(SELs) & !is.null(SELc)){
    
    ## sort SELs/SELc by c.id and s.id (based on OUT)
    idOUT <- with(OUT, paste(c.id, s.id, sep="_"))
    SELs$idSEL <- with(SELs, paste(c.id, s.id, sep="_"))
    SELc$idSEL <- with(SELc, paste(c.id, s.id, sep="_"))
    
    ## which matches in SELs are in OUT?
    SEL_top <- SELs[SELs$idSEL %in% idOUT,]
    SEL_top <- SEL_top[match(idOUT, SEL_top$idSEL),]
    SEL_bot <- SELs[! SELs$idSEL %in% idOUT,]
    SELs <- rbind(SEL_top, SEL_bot)
    rm(SEL_top, SEL_bot)
    
    ## which matches in SELc are in OUT?
    SEL_top <- SELc[SELc$idSEL %in% idOUT,]
    SEL_top <- SEL_top[match(idOUT, SEL_top$idSEL),]
    SEL_bot <- SELc[! SELc$idSEL %in% idOUT,]
    SELc <- rbind(SEL_top, SEL_bot)
    rm(SEL_top, SEL_bot)
    
    ## unique student and college ids
    indices <- SELs
    uColleges <- sort(unique(indices$c.id))
    uStudents <- sort(unique(indices$s.id))
    
  } else{
    
    ## get college and student ids from OUT
    c.id <- OUT[,c.id]
    s.id <- OUT[,s.id]
    
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
    
    ## index equilibrium matches from x in indices
    indices$id <- paste(indices$c.id, indices$s.id, sep="_")
    OUT$id    <- paste(OUT$c.id, OUT$s.id, sep="_")
    
    ## sort indices such that observed matches come first
    indices$D <- ifelse(indices$id %in% OUT$id, 1, 0)
    indices <- indices[order(indices$D, decreasing=TRUE),]
    
    ## sort OUT such that it concordes with indices
    OUT    <- OUT[match(indices$id[indices$D==1], OUT$id), ] 
  }
  
  ## -----------------------------------------------------------------------------------------
  ## --- 2. Produce matrices indicating the matches (H) and their position in the data (M) ---
  
  ## matrices are of dimension nColleges x nStudents
  ind.new       <- indices
  ind.new$index <- 1:nrow(ind.new)
  ind.org <- ind.new[order(ind.new$c.id,ind.new$s.id),]
  M <- matrix(ind.org$index, nrow=length(uColleges), ncol=length(uStudents), byrow=TRUE)
  H <- matrix(ind.org$D, nrow=length(uColleges), ncol=length(uStudents), byrow=TRUE)  
  
  if(method == "Klein"){
    
    ## preference inputs
    s.prefs <- s.prefs[[iter]]
    c.prefs <- c.prefs[[iter]]
    nSlots <- rowSums(H)
    
    ## find all stable matchings
    res <- hri(s.prefs=s.prefs, c.prefs=c.prefs, nSlots=nSlots)$matchings
    copt.id <- unique(res$matching[res$cOptimal==1]) # college-optimal matching
    res <- split(res, as.factor(res$matching))
    
    ## create adjecency matrices for all equilibrium matchings
    myfun <- function(x, type){
      H <- array(0, dim=c(length(unique(x[[1]][,type])), nrow(x[[1]]), length(x)))
      for(j in 1:length(x)){
        for(z in 1:nrow(x[[1]])){
          H[x[[j]][z,type], x[[j]][z,"student"], j] <- 1
        }
      }
      return(H)
    }
    H1 <- myfun(x=res, type="college"); names(H1) <- NULL # college admissions problem
    H2 <- myfun(x=res, type="slots"); names(H2) <- NULL # related stable marriage problem
    
    ## consistency check
    if( length( table(c(H) == c(H1[,,1]))) > 1 ){
      stop(paste("Data provided is not the student-optimal matching obtained from the 
                 preference lists in market ", iter, ".", sep=""))
    } else{
      H <- H1 # replace matrix for student-optimal matching (H) with all matchings (H1)
      rm(H1)
    }
  }
  
  ## ---------------------------------------------------------------------
  ## --- 3. Produce the datasets to be used based on formulas provided ---
  
  if(is.null(SEL) & is.null(SELs) & is.null(SELc)){
    
    c.OUT <- OUT[!duplicated(OUT$c.id),]
    C <- data.frame(apply(data.frame(c.OUT[,colleges]), 2, function(i) i[match(indices$c.id,c.OUT$c.id)]))
    names(C) <- colleges
    
    s.OUT <- OUT[!duplicated(OUT$s.id),]
    S <- data.frame(apply(data.frame(s.OUT[,students]), 2, function(i) i[match(indices$s.id,s.OUT$s.id)]))
    names(S) <- students
    
    ## ... and interaction effects
    Xmain <- data.frame(C, S)
    Xmatch <- rFormula(formula = outcome, data=OUT)
    
    if(method=="Klein"){
      
      C <- rFormula(formula = selection.college, data=Xmain)
      S <- rFormula(formula = selection.student, data=Xmain)    
      
    } else if(method == "Sorensen"){
      
      C <- rFormula(formula = selection, data=Xmain)
    }
    
  } else{
    
    Xmatch <- rFormula(formula = outcome, data=OUT)
    
    if(method=="Klein"){
      
      C <- rFormula(formula = selection.college, data=SELc)
      S <- rFormula(formula = selection.student, data=SELs)
      
    } else if(method == "Sorensen"){
      
      C <- rFormula(formula = selection, data=SEL)
    }
  }
  
  D <- indices$D
  d <- which(D==1)
  
  ## -----------------------------
  ## --- 4. Return the results ---
  
  if(method=="Klein"){
    
    return( list(Y=Xmatch[,1], Xmatch=Xmatch[,-1], C=C, Cmatch=C[D==1,], S=S, Smatch=S[D==1,], 
                 D=D, d=d, M=M, H=H, H2=H2, copt.id=copt.id, indices=indices$id) )
    
  } else if(method=="Sorensen"){
    
    return( list(Y=Xmatch[,1], Xmatch=Xmatch[,-1], C=C, Cmatch=C[D==1,], 
                 D=D, d=d, M=M, H=H, indices=indices$id) )
  }  
}  




#' @export
summary.stabit2 <- function(object, mfx=FALSE, ...){
  
  ## function to produce regression tables
  tab <- function(coefs, vcov){
    se <- sqrt(diag(vcov))
    tval <- coefs / se
    TAB <- cbind(Estimate = coefs,
                 StdErr = se,
                 t.value = tval,
                 p.value = 2*pt(-abs(tval), df=object$df))
    TAB
  }
  
  res <- list()
  
  ## Selection equation(s)
  
  if(object$method!="Outcome-only"){
    
    if(mfx==FALSE){
      
      if(object$method=="Klein"){
        
        res$college <- tab(coefs=object$coefs$beta, vcov=object$vcov$beta)
        
        res$student <- tab(coefs=object$coefs$gamma, vcov=object$vcov$gamma)
        
      } else{
        
        res$selection <- tab(coefs=object$coefs$beta, vcov=object$vcov$beta)
      }
      
    } else if(mfx==TRUE){
      
      if(object$method=="Klein"){
        
        res$college <- mfxVal(coefs = object$coefs$beta, 
                              vcov = object$vcov$beta,
                              df = nrow(object$variables$C) - ncol(object$variables$C) )
        
        res$student <- mfxVal(coefs = object$coefs$gamma, 
                              vcov = object$vcov$gamma,
                              df = nrow(object$variables$S) - ncol(object$variables$S) )
        
      } else if(object$method=="Sorensen"){
        
        res$selection <- mfxVal(coefs = object$coefs$beta, 
                                vcov = object$vcov$beta,
                                df = nrow(object$variables$C) - ncol(object$variables$C) )
      }  else{
        
        res$selection <- mfxVal(coefs = object$coefs$beta, 
                                vcov = object$vcov$beta,
                                df = nrow(object$variables$W) - ncol(object$variables$W) )
      }
    }
  }
  
  
  ## Outcome equation
  
  if(mfx==FALSE){
    
    res$outcome <- tab(coefs=object$coefs$alpha, vcov=object$vcov$alpha)
    
  } else{
    
    if(object$binary==TRUE){
      
      res$outcome <- mfxOut(sims = 10000, X = object$variables$X,
                            coefs = object$coefs$alpha,
                            vcov = object$vcov$alpha,
                            df = nrow(object$variables$X) - ncol(object$variables$X) )
    } else{
      
      res$outcome <- tab(coefs = object$coefs$alpha,
                         vcov = object$vcov$alpha )
    }
  }
  
  res$call <- object$call
  res$method <- object$method
  res$mfx <- mfx
  
  class(res) <- "summary.stabit2"
  res
}


mfxOut <- function(sims=10000, x.mean=TRUE, coefs, vcov, X, df){
  
  ## source: http://researchrepository.ucd.ie/handle/10197/3404
  ## method: average of individual marginal effects at each observation
  ## interpretation: http://www.indiana.edu/~statmath/stat/all/cdvm/cdvm.pdf page 8
  
  set.seed(1984)
  se <- sqrt(diag(vcov))
  
  if(x.mean==TRUE){
    
    ## marginal effects are calculated at the means of independent variables
    pdf <- dnorm(mean(X%*%coefs))
    pdfsd <- dnorm(sd(X%*%coefs))
    
  } else{
    
    ## marginal effects are calculated for each observation and then averaged
    pdf <- mean(dnorm(X%*%coefs))
    pdfsd <- sd(dnorm(X%*%coefs))
  }  
  mx <- pdf*coefs
  
  sim <- matrix(rep(NA,sims*length(coefs)), nrow=sims)
  
  for(i in 1:length(coefs)){
    sim[,i] <- rnorm(sims,coefs[i],se[i])
  }
  
  pdfsim <- rnorm(sims,pdf,pdfsd)
  sim.se <- pdfsim*sim
  s.e. <- apply(sim.se,2,sd)
  
  tval <- mx / s.e.
  TAB <- cbind(Estimate = mx,
               StdErr = s.e.,
               t.value = tval,
               p.value = pt(-abs(tval), df=df))
  TAB
}


mfxVal <- function(coefs, vcov, df){
  
  ## Reference: Sorensen (2007, p. 2748)
  se <- sqrt(diag(vcov))
  mx <- dnorm(0)*coefs/sqrt(2)
  s.e. <- dnorm(0)*se/sqrt(2)
  tval <- mx / s.e.
  TAB <- cbind(Estimate = mx,
               StdErr = s.e.,
               t.value = tval,
               p.value = pt(-abs(tval), df=df))
  TAB
}




#' @export
print.summary.stabit2 <- function(x, ...){
  
  if(x$mfx==TRUE){
    cat("\nMarginal effects for multi-index sample selection model.")
  } else{
    cat("\nCoefficients for multi-index sample selection model.")
  }
  
  if(x$method=="Klein"){
    cat("\nMethod: Klein (2016), two-sided matching market\n")
  } else if(x$method=="Sorensen"){
    cat("\nMethod: Sorensen (2007), two-sided matching market\n")
  } else{
    cat("\nMethod: Klein (2015), one-sided matching market\n")
  }
  
  cat("\nCall:\n")
  print(x$call)
  
  if(x$method!="Outcome-only"){
    
    if(x$method=="Klein"){
      
      cat("\nSelection equation (Valuation over colleges):")
      cat("\n")
      printCoefmat(x$college, P.values=TRUE, has.Pvalue=TRUE, signif.legend=FALSE)
      
      cat("\nSelection equation (Valuation over students):")
      cat("\n")
      printCoefmat(x$student, P.values=TRUE, has.Pvalue=TRUE, signif.legend=FALSE)
      
    } else{
      
      cat("\nSelection equation:")
      cat("\n")
      printCoefmat(x$selection, P.values=TRUE, has.Pvalue=TRUE, signif.legend=FALSE) 
    }
  }
  
  cat("\nOutcome equation:")
  cat("\n")
  printCoefmat(x$outcome, P.values=TRUE, has.Pvalue=TRUE, signif.legend=TRUE)
  
}




#' @export
predict.stabit2 <- function(object, newdata=NULL, ...){
  if(is.null(newdata))
    y <- fitted(object)
  else{
    if(!is.null(object$formula)){
      ## model has been fitted using formula interface
      x <- model.matrix(object$formula, newdata)
    }
    else{
      x <- newdata
    }
    y <- as.vector(x %*% object$coefs$alpha)
  }
  y
}


