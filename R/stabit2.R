# ----------------------------------------------------------------------------
# R-code (www.r-project.org/) for the Structural Matching Model
#
# Copyright (c) 2016 Thilo Klein
#
# This library is distributed under the terms of the GNU Public License (GPL)
# for full details see the file LICENSE
#
# ----------------------------------------------------------------------------
#' @title Matching model and selection correction for college admissions
#'
#' @description The function provides a Gibbs sampler for a structural matching model that 
#' estimates preferences and corrects for sample selection bias when the selection process 
#' is a two-sided matching game; i.e., a matching of students to colleges.
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
#' @param nCores number of cores to be used in parallel Gibbs sampling.
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
#' @useDynLib matchingMarkets, .registration = TRUE 
#' 
#' @import stats lattice parallel RcppProgress
#' @importFrom Rcpp evalCpp
#' @importFrom graphics par plot
#' 
#' @author Thilo Klein 
#' 
#' @keywords regression
#' 
#' @references Sorensen, M. (2007). How Smart is Smart Money? A Two-Sided Matching Model of Venture Capital.
#' \emph{Journal of Finance}, 62 (6): 2725-2762.
#' 
#' @examples
#' \dontrun{
#' ## --- SIMULATED EXAMPLE ---
#' 
#' ## 1. Simulate two-sided matching data for 20 markets (m=20) with 100 students
#' ##    (nStudents=100) per market and 20 colleges with quotas of 5 students, each
#' ##    (nSlots=rep(5,20)). True parameters in selection and outcome equations are 
#' ##    all equal to 1.
#' 
#' xdata <- stabsim2(m=20, nStudents=100, nSlots=rep(5,20), verbose=FALSE,
#'   colleges = "c1", students = "s1",
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
#'  #coef(fit2)
#'  #fitted(fit2)
#'  #residuals(fit2)
#'  #predict(fit2, newdata=NULL)
#'
#'    
#' ## 5. Plot MCMC draws for coefficients
#'  plot(fit2)
#' }
stabit2 <- function(OUT=NULL, SEL=NULL, colleges=NULL, students=NULL, outcome=NULL, selection,
                    binary=FALSE, niter, gPrior=FALSE, 
                    censored=1, thin=1, nCores=max(1,detectCores()-1), ...) UseMethod("stabit2")

#' @export
stabit2.default <- function(OUT=NULL, SEL=NULL, colleges=NULL, students=NULL, outcome=NULL, selection,
                            binary=FALSE, niter, gPrior=FALSE, 
                            censored=1, thin=1, nCores=max(1,detectCores()-1), ...){
  
  ## ------------------------
  ## --- 1. Preliminaries ---
  
  ## select method based on arguments provided
  if(is.list(selection)){
    method <- "Klein" # separate selection equations for students and college preferences
    selection.college <- selection$college
    selection.student <- selection$student
    selection <- NULL
    if(is.null(OUT)){
      method <- "Klein-selection" # no outcome equation
      OUT <- SEL$SELs[SEL$SELs$D==1,]
      OUT$Y <- rep(0,nrow(OUT))
    }
  } else{
    method <- "Sorensen" # single selection equation with equal sharing rule for student and college utility
  }
  
  s.prefs <- NULL
  c.prefs <- NULL
  
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
  
  ## number of cores need not exceed number of markets
  nCores <- min(length(SELc), nCores)
  
  ## market and agent identifiers
  m.id <- "m.id"
  c.id <- "c.id"
  s.id <- "s.id"
  
  ## empty lists for results
  Y=list(); Xmatch=list(); C=list(); Cmatch=list(); S=list(); Smatch=list(); D=list(); d=list()
  M=list(); H=list(); indices=list()
  c.better=list(); c.worse=list(); s.better=list(); s.worse=list()
  c.betterNA=list(); c.worseNA=list(); s.betterNA=list(); s.worseNA=list()
  
  for(i in 1:length(OUT)){
    
    ## ------------------------------------------------------------------------------------------
    ## --- 2. Bring data in consistent order and produce matrices for all equilibrium matches ---
    
    X <- stabit2_inner(iter=i, OUT=OUT[[i]], SEL=SEL[[i]], SELs=SELs[[i]], SELc=SELc[[i]],
                       colleges=colleges, students=students, outcome=outcome, selection=selection,
                       selection.student=selection.student, selection.college=selection.college, 
                       s.prefs=s.prefs[[i]], c.prefs=c.prefs[[i]],
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
    if(method=="Klein" | method=="Klein-selection"){
      S[[i]]      <- as.matrix(X$S)
      Smatch[[i]] <- as.matrix(X$Smatch)
      H[[i]]      <- X$H
    } else{
      H[[i]]      <- as.matrix(X$H)
    }
    c.better[[i]]  <- X$c.better
    c.worse[[i]]   <- X$c.worse
    s.better[[i]]  <- X$s.better
    s.worse[[i]]   <- X$s.worse
    c.betterNA[[i]]  <- X$c.betterNA
    c.worseNA[[i]]   <- X$c.worseNA
    s.betterNA[[i]]  <- X$s.betterNA
    s.worseNA[[i]]   <- X$s.worseNA
  }
  
  ## some preliminaries
  T <- length(Y); #// Number of markets.
  nColleges <- nStudents <- XXmatch <- CC <- SS <- CCmatch <- SSmatch <- list()
  L <- rep(list(vector()),T)
  studentIds <- collegeId <- rep(list(vector()),T)
  
  ## further results
  s <- 0
  for(i in 1:T){
    nColleges[[i]] <- nrow(H[[i]])
    nStudents[[i]] <- ncol(H[[i]])
    XXmatch[[i]]   <- t(Xmatch[[i]]) %*% Xmatch[[i]]
    CC[[i]]        <- t(C[[i]]) %*% C[[i]]
    CCmatch[[i]]   <- t(Cmatch[[i]]) %*% Cmatch[[i]]
    if(method=="Klein" | method=="Klein-selection"){
      SSmatch[[i]]   <- t(Smatch[[i]]) %*% Smatch[[i]]
      SS[[i]]        <- t(S[[i]]) %*% S[[i]]
      #collegeId[[i]] <- matrix(NA, nrow=nStudents[[i]], ncol=dim(H[[i]])[3]) ##
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
  n <- sum(unlist(nStudents)) # total number of students/matches
  N <- sum(unlist(nColleges)) # total number of colleges
  
  ## ----------------------------
  ## --- 3. Run Gibbs sampler ---
  
  #save.image("~/Desktop/play.RData")
  #load("~/Desktop/play.RData")
  
  if(nCores > 1){
    
    ## split market ids (Ts) and Ns evenly to the number of cores (nCores)
    Ts <- split(1:T, cut(1:T, breaks=nCores, labels=FALSE))
    Ns <- list()
    for(i in 1:length(Ts)){
      if(i==1){
        Ns[[i]] <- 1:sum(unlist(nColleges[ Ts[[i]] ]))
      } else{
        Ns[[i]] <- max(Ns[[i-1]]) + ( 1:sum(unlist(nColleges[ Ts[[i]] ])) )
        L[ Ts[[i]] ] <- lapply(L[ Ts[[i]] ], function(z){
          z - min(unlist(L[ Ts[[i]] ])) 
        })
      }
    }  
    
    ## setup cluster
    cl <- makeCluster(nCores) 
    clusterEvalQ(cl, library(matchingMarkets))
    
    if(method == "Klein" & (!is.null(s.prefs) & !is.null(c.prefs)) ){
      
      ## split 
      parObject <- lapply(1:nCores, function(i){
        list(Y=Y[Ts[[i]]], Xmatch=Xmatch[Ts[[i]]], C=C[Ts[[i]]], Cmatch=Cmatch[Ts[[i]]], 
        S=S[Ts[[i]]], Smatch=Smatch[Ts[[i]]], D=D[Ts[[i]]], d=d[Ts[[i]]], M=M[Ts[[i]]], 
        H=H[Ts[[i]]], nColleges=unlist(nColleges[Ts[[i]]]), nStudents=unlist(nStudents[Ts[[i]]]), 
        XXmatch=XXmatch[Ts[[i]]], CC=CC[Ts[[i]]], SS=SS[Ts[[i]]], SSmatch=SSmatch[Ts[[i]]], 
        CCmatch=CCmatch[Ts[[i]]], L=L[Ts[[i]]], studentIds=studentIds[Ns[[i]]], 
        collegeId=collegeId[Ts[[i]]], n=sum(unlist(nStudents[Ts[[i]]])), N=sum(unlist(nColleges[Ts[[i]]])), 
        binary=binary, niter=niter, thin=thin, T=length(Ts[[i]]), censored=censored,
        c.better=c.better[Ts[[i]]], c.worse=c.worse[Ts[[i]]], s.better=s.better[Ts[[i]]], s.worse=s.worse[Ts[[i]]],
        c.betterNA=c.betterNA[Ts[[i]]], c.worseNA=c.worseNA[Ts[[i]]], s.betterNA=s.betterNA[Ts[[i]]], s.worseNA=s.worseNA[Ts[[i]]])
      })
      
      cat(paste("Running parallel Gibbs sampler on", nCores, "cores...", sep=" "))
      res <- parLapply(cl, parObject, function(d){
        
        with(d, stabit2Sel0(Yr=Y, Xmatchr=Xmatch, Cr=C, Cmatchr=Cmatch, 
                            Sr=S, Smatchr=Smatch, Dr=D, dr=d, Mr=M, 
                            Hr=H, nCollegesr=nColleges, nStudentsr=nStudents, 
                            XXmatchr=XXmatch, CCr=CC, SSr=SS, SSmatchr=SSmatch, 
                            CCmatchr=CCmatch, Lr=L, studentIdsr=studentIds, 
                            collegeIdr=collegeId, n=n, N=N, 
                            binary=binary, niter=niter, thin=thin, T=T, censored=censored,
                            cbetterr=c.better, cworser=c.worse, sbetterr=s.better, sworser=s.worse,
                            cbetterNAr=c.betterNA, cworseNAr=c.worseNA, sbetterNAr=s.betterNA, sworseNAr=s.worseNA))
      })
      
    } else if(method == "Klein" & (is.null(s.prefs) | is.null(c.prefs)) ){
      
      ## split 
      parObject <- lapply(1:nCores, function(i){
        list(Y=Y[Ts[[i]]], Xmatch=Xmatch[Ts[[i]]], C=C[Ts[[i]]], Cmatch=Cmatch[Ts[[i]]], 
             S=S[Ts[[i]]], Smatch=Smatch[Ts[[i]]], D=D[Ts[[i]]], d=d[Ts[[i]]], M=M[Ts[[i]]], 
             H=H[Ts[[i]]], nColleges=unlist(nColleges[Ts[[i]]]), nStudents=unlist(nStudents[Ts[[i]]]), 
             XXmatch=XXmatch[Ts[[i]]], CC=CC[Ts[[i]]], SS=SS[Ts[[i]]], SSmatch=SSmatch[Ts[[i]]], 
             CCmatch=CCmatch[Ts[[i]]], L=L[Ts[[i]]], studentIds=studentIds[Ns[[i]]], 
             collegeId=collegeId[Ts[[i]]], n=sum(unlist(nStudents[Ts[[i]]])), N=sum(unlist(nColleges[Ts[[i]]])), 
             binary=binary, niter=niter, thin=thin, T=length(Ts[[i]]), censored=censored)
      })
      
      cat(paste("Running parallel Gibbs sampler on", nCores, "cores...", sep=" "))
      res <- parLapply(cl, parObject, function(d){
        
        with(d, stabit2Sel1(Yr=Y, Xmatchr=Xmatch, Cr=C, Cmatchr=Cmatch, 
                            Sr=S, Smatchr=Smatch, Dr=D, dr=d, Mr=M, 
                            Hr=H, nCollegesr=nColleges, nStudentsr=nStudents, 
                            XXmatchr=XXmatch, CCr=CC, SSr=SS, SSmatchr=SSmatch, 
                            CCmatchr=CCmatch, Lr=L, studentIdsr=studentIds, 
                            collegeIdr=collegeId, n=n, N=N, 
                            binary=binary, niter=niter, thin=thin, T=T, censored=censored))
      })
      
    } else if(method == "Klein-selection" & (!is.null(s.prefs) & !is.null(c.prefs)) ){
      
      ## split
      parObject <- lapply(1:nCores, function(i){
        list(C=C[Ts[[i]]], Cmatch=Cmatch[Ts[[i]]], S=S[Ts[[i]]], Smatch=Smatch[Ts[[i]]], D=D[Ts[[i]]], 
             d=d[Ts[[i]]], M=M[Ts[[i]]], H=H[Ts[[i]]], nColleges=unlist(nColleges[Ts[[i]]]), 
             nStudents=unlist(nStudents[Ts[[i]]]), CC=CC[Ts[[i]]], SS=SS[Ts[[i]]], 
             SSmatch=SSmatch[Ts[[i]]], CCmatch=CCmatch[Ts[[i]]], L=L[Ts[[i]]], 
             studentIds=studentIds[Ns[[i]]], collegeId=collegeId[Ts[[i]]], n=sum(unlist(nStudents[Ts[[i]]])), 
             N=sum(unlist(nColleges[Ts[[i]]])), niter=niter, thin=thin, T=length(Ts[[i]]),
             c.better=c.better[Ts[[i]]], c.worse=c.worse[Ts[[i]]], s.better=s.better[Ts[[i]]], s.worse=s.worse[Ts[[i]]],
             c.betterNA=c.betterNA[Ts[[i]]], c.worseNA=c.worseNA[Ts[[i]]], s.betterNA=s.betterNA[Ts[[i]]], s.worseNA=s.worseNA[Ts[[i]]])
      })
      
      cat(paste("Running parallel Gibbs sampler on", nCores, "cores...", sep=" "))
      res <- parLapply(cl, parObject, function(d){  
        
        with(d, stabit2Mat0(Cr=C, Cmatchr=Cmatch, Sr=S, Smatchr=Smatch, Dr=D, 
                            dr=d, Mr=M, Hr=H, nCollegesr=nColleges, 
                            nStudentsr=nStudents, CCr=CC, SSr=SS, 
                            SSmatchr=SSmatch, CCmatchr=CCmatch, Lr=L, 
                            studentIdsr=studentIds, collegeIdr=collegeId, n=n, 
                            N=N, niter=niter, thin=thin, T=T,
                            cbetterr=c.better, cworser=c.worse, sbetterr=s.better, sworser=s.worse,
                            cbetterNAr=c.betterNA, cworseNAr=c.worseNA, sbetterNAr=s.betterNA, sworseNAr=s.worseNA))
      })

    } else if(method == "Klein-selection" & (is.null(s.prefs) | is.null(c.prefs)) ){
      
      ## split
      parObject <- lapply(1:nCores, function(i){
        list(C=C[Ts[[i]]], Cmatch=Cmatch[Ts[[i]]], S=S[Ts[[i]]], Smatch=Smatch[Ts[[i]]], D=D[Ts[[i]]], 
             d=d[Ts[[i]]], M=M[Ts[[i]]], H=H[Ts[[i]]], nColleges=unlist(nColleges[Ts[[i]]]), 
             nStudents=unlist(nStudents[Ts[[i]]]), CC=CC[Ts[[i]]], SS=SS[Ts[[i]]], 
             SSmatch=SSmatch[Ts[[i]]], CCmatch=CCmatch[Ts[[i]]], L=L[Ts[[i]]], 
             studentIds=studentIds[Ns[[i]]], collegeId=collegeId[Ts[[i]]], n=sum(unlist(nStudents[Ts[[i]]])), 
             N=sum(unlist(nColleges[Ts[[i]]])), niter=niter, thin=thin, T=length(Ts[[i]]))
      })
      
      cat(paste("Running parallel Gibbs sampler on", nCores, "cores...", sep=" "))
      res <- parLapply(cl, parObject, function(d){  
        
        with(d, stabit2Mat1(Cr=C, Cmatchr=Cmatch, Sr=S, Smatchr=Smatch, Dr=D, 
                            dr=d, Mr=M, Hr=H, nCollegesr=nColleges, 
                            nStudentsr=nStudents, CCr=CC, SSr=SS, 
                            SSmatchr=SSmatch, CCmatchr=CCmatch, Lr=L, 
                            studentIdsr=studentIds, collegeIdr=collegeId, n=n, 
                            N=N, niter=niter, thin=thin, T=T))
      })
            
    } else if(method == "Sorensen"){

      ## split
      parObject <- lapply(1:nCores, function(i){
        list(Y=Y[Ts[[i]]], Xmatch=Xmatch[Ts[[i]]], C=C[Ts[[i]]], Cmatch=Cmatch[Ts[[i]]], 
             D=D[Ts[[i]]], d=d[Ts[[i]]], M=M[Ts[[i]]], H=H[Ts[[i]]], 
             nColleges=unlist(nColleges[Ts[[i]]]), nStudents=unlist(nStudents[Ts[[i]]]), 
             XXmatch=XXmatch[Ts[[i]]], CC=CC[Ts[[i]]], CCmatch=CCmatch[Ts[[i]]], L=L[Ts[[i]]], 
             studentIds=studentIds[Ns[[i]]], collegeId=collegeId[Ts[[i]]], 
             n=sum(unlist(nStudents[Ts[[i]]])), N=sum(unlist(nColleges[Ts[[i]]])), 
             binary=binary, niter=niter, thin=thin, T=length(Ts[[i]]), censored=censored)
      })
      
      cat(paste("Running parallel Gibbs sampler on", nCores, "cores...", sep=" "))
      res <- parLapply(cl, parObject, function(d){  
        
        with(d, stabit2Sel2(Yr=Y, Xmatchr=Xmatch, Cr=C, Cmatchr=Cmatch, 
                            Dr=D, dr=d, Mr=M, Hr=H, 
                            nCollegesr=nColleges, nStudentsr=nStudents, 
                            XXmatchr=XXmatch, CCr=CC, CCmatchr=CCmatch, Lr=L, 
                            studentIdsr=studentIds, collegeIdr=collegeId, 
                            n=n, N=N, 
                            binary=binary, niter=niter, thin=thin, T=T, censored=censored))
      })
    }
    
    rm(parObject, Ts, Ns)
    stopCluster(cl)
    
    ## split results into draws for paramters (res1) and error terms (res2)
    res1 <- lapply(res, function(z) z[!names(z) %in% c("etadraws", "deltadraws")])
    res2 <- lapply(res, function(z) z[names(z) %in% c("etadraws", "deltadraws")])
    rm(res)
    
    ## combine results from different cores for error terms
    est2 <- lapply(1:length(res2[[1]]), function(z){
      do.call(rbind, lapply(res2, function(d) d[[z]] ))
    })
    names(est2) <- names(res2[[1]]); rm(res2)
    
    ## combine results from different cores for parameter draws
    toarray <- function(x, y){
      array(unlist(lapply(x, function(z) z[[y]] )), 
            dim = c(nrow(res1[[1]][[y]]), ncol(res1[[1]][[y]]), length(res1)) )
    }
    est1 <- lapply(1:length(res1[[1]]), function(z){
      toarray( res1, names(res1[[1]])[z] )
    })
    names(est1) <- names(res1[[1]]); rm(res1)
    est1 <- lapply(est1, consensusMC)
    
    ## combine results for paramters (est1) and error terms (est2)
    est <- c(est1, est2); rm(est1, est2)
    
  } else{
    #Ts <- 1:T
    #Ns <- 1:N  
    
    if(method=="Klein" & !is.null(s.prefs) & !is.null(c.prefs) ){
      
      est <- stabit2Sel0(Yr=Y, Xmatchr=Xmatch, Cr=C, Cmatchr=Cmatch, Sr=S, Smatchr=Smatch, Dr=D, dr=d,
                         Mr=M, Hr=H, nCollegesr=unlist(nColleges), nStudentsr=unlist(nStudents), XXmatchr=XXmatch, 
                         CCr=CC, SSr=SS, SSmatchr=SSmatch, CCmatchr=CCmatch, Lr=L,
                         studentIdsr=studentIds, collegeIdr=collegeId, n=n, N=N, binary=binary, niter=niter, 
                         thin=thin, T=T, censored=censored,
                         cbetterr=c.better, cworser=c.worse, sbetterr=s.better, sworser=s.worse,
                         cbetterNAr=c.betterNA, cworseNAr=c.worseNA, sbetterNAr=s.betterNA, sworseNAr=s.worseNA)
      
    } else if(method=="Klein" & (is.null(s.prefs) | is.null(c.prefs)) ){
      
      est <- stabit2Sel1(Yr=Y, Xmatchr=Xmatch, Cr=C, Cmatchr=Cmatch, Sr=S, Smatchr=Smatch, Dr=D, dr=d,
                         Mr=M, Hr=H, nCollegesr=unlist(nColleges), nStudentsr=unlist(nStudents), XXmatchr=XXmatch, 
                         CCr=CC, SSr=SS, SSmatchr=SSmatch, CCmatchr=CCmatch, Lr=L,
                         studentIdsr=studentIds, collegeIdr=collegeId, n=n, N=N, binary=binary, niter=niter, 
                         thin=thin, T=T, censored=censored)
      
    } else if(method=="Klein-selection" & !is.null(s.prefs) & !is.null(c.prefs) ){
      
      est <- stabit2Mat0(Cr=C, Cmatchr=Cmatch, Sr=S, Smatchr=Smatch, Dr=D, dr=d,
                         Mr=M, Hr=H, nCollegesr=unlist(nColleges), nStudentsr=unlist(nStudents),
                         CCr=CC, SSr=SS, SSmatchr=SSmatch, CCmatchr=CCmatch, Lr=L,
                         studentIdsr=studentIds, collegeIdr=collegeId, n=n, N=N, niter=niter, thin=thin, T=T,
                         cbetterr=c.better, cworser=c.worse, sbetterr=s.better, sworser=s.worse,
                         cbetterNAr=c.betterNA, cworseNAr=c.worseNA, sbetterNAr=s.betterNA, sworseNAr=s.worseNA)

    } else if(method=="Klein-selection" & (is.null(s.prefs) | is.null(c.prefs)) ){
      
      est <- stabit2Mat1(Cr=C, Cmatchr=Cmatch, Sr=S, Smatchr=Smatch, Dr=D, dr=d,
                         Mr=M, Hr=H, nCollegesr=unlist(nColleges), nStudentsr=unlist(nStudents),
                         CCr=CC, SSr=SS, SSmatchr=SSmatch, CCmatchr=CCmatch, Lr=L,
                         studentIdsr=studentIds, collegeIdr=collegeId, n=n, N=N, niter=niter, thin=thin, T=T)
      
    } else if(method=="Sorensen"){
      
      est <- stabit2Sel2(Yr=Y, Xmatchr=Xmatch, Cr=C, Cmatchr=Cmatch, Dr=D, dr=d,
                         Mr=M, Hr=H, nCollegesr=unlist(nColleges), nStudentsr=unlist(nStudents), XXmatchr=XXmatch, 
                         CCr=CC, CCmatchr=CCmatch, Lr=L, 
                         studentIdsr=studentIds, collegeIdr=collegeId, 
                         n=n, N=N, binary=binary, niter=niter, thin=thin, T=T, censored=censored)
    }
  }
  
  # ------------------------------------
  # --- 5. Add names to coefficients ---
  
  ## variable names
  
  if(method!="Klein-selection"){
    an <- colnames(Xmatch[[1]])
  }
  bn <- colnames(C[[1]])
  if(method=="Klein" | method=="Klein-selection"){
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
  } else if(method=="Klein-selection"){
    rownames(est$gammadraws) = cn
  }
  
  ## posterior means
  
  ## The last half of all draws are used in approximating the posterior means and the posterior standard deviations.
  niter <- dim(est$betadraws)[2]
  startiter <- floor(niter/2)
  
  posMeans <- function(x){
    sapply(1:nrow(x), function(z) mean(x[z,startiter:niter]))
  }
  
  if(method!="Klein-selection"){
    est$alpha <- posMeans(est$alphadraws) 
    names(est$alpha) <- rownames(est$alphadraws)
  }
  
  est$beta <- posMeans(est$betadraws)
  names(est$beta) <- bn
  
  if(method=="Klein" | method=="Klein-selection"){
    est$gamma <- posMeans(est$gammadraws)
    names(est$gamma) <- cn
  }
  
  if(method!="Klein-selection" & binary==FALSE){
    est$sigma <- sqrt(posMeans(est$sigmasquarenudraws))
    est$sigmasquarenudraws <- NULL
  }
  
  ## error terms
  
  est$eta <- posMeans(est$etadraws)
  est$etadraws <- NULL
  if(method=="Klein" | method=="Klein-selection"){
    est$delta <- posMeans(est$deltadraws)
    est$deltadraws <- NULL
  }
  
  ## vcov
  
  if(method!="Klein-selection"){
    est$vcov$alpha <- var(t(est$alphadraws))
  }
  est$vcov$beta <- var(t(est$betadraws))
  if(method=="Klein" | method=="Klein-selection"){
    est$vcov$gamma <- var(t(est$gammadraws))
  }
  
  ## variables
  
  est$C <- do.call(rbind, C)
  est$D <- do.call(c, D)
  
  if(method!="Klein-selection"){
    est$Y <- do.call(c, Y)
    est$X <- do.call(rbind, Xmatch)
    
    if(method=="Klein"){
      est$X <- cbind(est$X, eta=est$eta, delta=est$delta)
      est$S <- do.call(rbind, S)  
    } else if(method=="Sorensen"){
      est$X <- cbind(est$X, eta=est$eta)
    }
  } else if(method=="Klein-selection"){
    est$S <- do.call(rbind, S)
  } 
  
  ## other output
  
  if(method=="Sorensen"){
    
    est$coefficients <- unlist(with(est, list(o=alpha, s=beta)))
    
  } else if(method=="Klein"){
    
    est$coefficients <- unlist(with(est, list(o=alpha, c=beta, s=gamma)))
    
  } else if(method=="Klein-selection"){
    
    est$coefficients <- unlist(with(est, list(c=beta, s=gamma)))
    
  }
  
  if(method!="Klein-selection"){
    est$fitted.values <- with(est, as.vector(X %*% alpha))
    est$residuals <- est$Y - est$fitted.values
    est$df <- n - ncol(est$X)
    est$binary <- binary
    est$formula <- outcome
  } else{
    #est$fitted.values <- with(est, as.vector(X %*% alpha))
    #est$residuals <- est$Y - est$fitted.values
    est$df <- n*N - ncol(est$C) - ncol(est$S)
    #est$binary <- binary
    #est$formula <- outcome
  }
  
  est$call <- match.call()
  est$method <- method
  
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




consensusMC <- function(subchain){
  
  ddata = length(dim(subchain))
  
  d          <- dim(subchain)[1]  
  sampletotT <- dim(subchain)[2]    
  M          <- dim(subchain)[3]
  
  if ( M==1 ){ 
    theta <- array(subchain[,,1],c(d,sampletotT))
    return (theta)
  }
  
  ## compute sigmahatm & sigmahatm.inverse (= W_s)  
  sigmahatm     <- array(NA,c(d,d,M))    
  sigmahatM     <- matrix(NA,d,d)
  sigmahatM.pre <- matrix(NA,d,d)        
  if (d==1){
    for (k in 1:M) { sigmahatm[1,,k] <- var(subchain[1,,k]) }    
  } else{
    for (k in 1:M) { sigmahatm[,,k] <- cov(t(subchain[,,k])) }
  }
  
  ## compute inverses of covariance matrices (with try()):
  sigmahatm.inverse <- array(NA,dim=c(d,d,M))
  for (k in 1:M){    
    res <- try( sigmahatm.inverse[,,k] <- solve(sigmahatm[,,k]), silent=TRUE)  
  }            
  sigmahatM <- solve(rowSums(sigmahatm.inverse, dims=2)) # (sum W_s)^(-1)
  
  ## compute unified posterior samples
  theta <- matrix(NA,nrow=d,ncol=sampletotT) # the resulting posterior samples
  wvec  <- array(NA, c(d,1))
  for (i in 1:sampletotT){
    wvec <- rep(0,d)
    for (s in 1:M){
      wvec <- wvec + sigmahatm.inverse[,,s] %*% subchain[,i,s] 
    }
    theta[,i] <- sigmahatM %*% wvec        
  }
  return(theta)
}




stabit2_inner <- function(iter, OUT, SEL, SELs, SELc, colleges, students, 
                          m.id="m.id", c.id="c.id", s.id="s.id", 
                          outcome, selection, selection.student, selection.college, 
                          s.prefs, c.prefs, method){
  
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

  ## -----------------------------------------------------------------------------------------
  ## --- 3. Create upper/lower bounds on latent match valuations based on rank order lists ---
  
  if(!is.null(s.prefs) & !is.null(c.prefs)){
    
    cprefs_better <- rbind(NA, c.prefs[-dim(c.prefs)[1],])
    cprefs_worse  <- rbind(c.prefs[-1,], NA)
    
    sprefs_better <- rbind(NA, s.prefs[-dim(s.prefs)[1],])
    sprefs_worse  <- rbind(s.prefs[-1,], NA)
    
    mapper <- function(clmn, direction, perspective){
      
      ## clmn       : column in rank order list (cprefs/sprefs)
      ## direction  : worse or better
      ## perspective: student or college perspective
      
      if(perspective == "college"){
        
        ## creates lookup matrix c.better (c.worse) of dimension CxS such that cell (c,s) gives
        ## student that college c ranks just above (below) student s.
        
        if(direction == "worse"){

          h <- merge(x=data.frame(x=1:max(uStudents)), 
                     y=data.frame(cbind(x=c.prefs[,clmn], y=cprefs_worse[,clmn])), 
                     by="x", all.x=TRUE); h[order(h$x),"y"]
        } else{ ## "better"

          h <- merge(x=data.frame(x=1:max(uStudents)), 
                     y=data.frame(cbind(x=c.prefs[,clmn], y=cprefs_better[,clmn])), 
                     by="x", all.x=TRUE); h[order(h$x),"y"]
          
        } ## end direction
        
      } else{ ## "student"
        
        ## creates lookup matrix s.better (s.worse) of dimension CxS such that cell (c,s) gives
        ## college that student s ranks just above (below) college c.
        
        if(direction == "worse"){

          h <- merge(x=data.frame(x=1:max(uColleges)), 
                     y=data.frame(cbind(x=s.prefs[,clmn], y=sprefs_worse[,clmn])), 
                     by="x", all.x=TRUE); h[order(h$x),"y"]
        } else{ ## "better"
          
          h <- merge(x=data.frame(x=1:max(uColleges)), 
                     y=data.frame(cbind(x=s.prefs[,clmn], y=sprefs_better[,clmn])), 
                     by="x", all.x=TRUE); h[order(h$x),"y"]
        } ## end direction
      } ## end perspective
    }
    
    c.worse    <- t(sapply(1:ncol(c.prefs), function(x) mapper(clmn=x, direction="worse", perspective="college"))) -1
    c.better   <- t(sapply(1:ncol(c.prefs), function(x) mapper(clmn=x, direction="better", perspective="college"))) -1
    c.worseNA  <- c.worse;  c.worseNA[!is.na(c.worseNA)]   <- 0; c.worseNA[is.na(c.worseNA)]   <- 1
    c.betterNA <- c.better; c.betterNA[!is.na(c.betterNA)] <- 0; c.betterNA[is.na(c.betterNA)] <- 1
    c.worse[is.na(c.worse)] <- 0
    c.better[is.na(c.better)] <- 0
    
    s.worse    <- sapply(1:ncol(s.prefs), function(x) mapper(clmn=x, direction="worse", perspective="student")) -1
    s.better   <- sapply(1:ncol(s.prefs), function(x) mapper(clmn=x, direction="better", perspective="student")) -1
    s.worseNA  <- s.worse;  s.worseNA[!is.na(s.worseNA)]   <- 0; s.worseNA[is.na(s.worseNA)]   <- 1
    s.betterNA <- s.better; s.betterNA[!is.na(s.betterNA)] <- 0; s.betterNA[is.na(s.betterNA)] <- 1
    s.worse[is.na(s.worse)] <- 0
    s.better[is.na(s.better)] <- 0
    
  } else{
    
    s.better <- NULL; c.better <- NULL; s.worse <- NULL; c.worse <- NULL
    s.betterNA <- NULL; c.betterNA <- NULL; s.worseNA <- NULL; c.worseNA <- NULL
  }
  
  ## ---------------------------------------------------------------------
  ## --- 4. Produce the datasets to be used based on formulas provided ---
  
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
    
    if(method!="Klein-selection"){
      Xmatch <- rFormula(formula = outcome, data=OUT)
    } else{
      Xmatch <- rFormula(formula = selection.student, data=OUT)
    }
    
    if(method=="Klein" | method=="Klein-selection"){
      
      C <- rFormula(formula = selection.college, data=SELc)
      S <- rFormula(formula = selection.student, data=SELs)
      
    } else if(method == "Sorensen"){
      
      C <- rFormula(formula = selection, data=SEL)
    }
  }
  
  D <- indices$D
  d <- which(D==1)
  
  ## -----------------------------
  ## --- 5. Return the results ---
  
  if(method=="Klein" | method=="Klein-selection"){
    
    return( list(Y=Xmatch[,1], Xmatch=Xmatch[,-1], C=C, Cmatch=C[D==1,], S=S, Smatch=S[D==1,], 
                 D=D, d=d, M=M, H=H, indices=indices$id, 
                 c.better=c.better, c.worse=c.worse, s.better=s.better, s.worse=s.worse,
                 s.betterNA=s.betterNA, c.betterNA=c.betterNA, s.worseNA=s.worseNA, c.worseNA=c.worseNA) )
    
  } else if(method=="Sorensen"){
    
    return( list(Y=Xmatch[,1], Xmatch=Xmatch[,-1], C=C, Cmatch=C[D==1,], 
                 D=D, d=d, M=M, H=H, indices=indices$id) )
  }  
}  



