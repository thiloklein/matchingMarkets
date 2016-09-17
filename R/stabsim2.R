# ----------------------------------------------------------------------------
# R-code (www.r-project.org/) for simulating data for
# all players in the market.
#
# Copyright (c) 2015 Thilo Klein
#
# This library is distributed under the terms of the GNU Public License (GPL)
# for full details see the file LICENSE
#
# ----------------------------------------------------------------------------

#' @title Simulated data for college admissions problem
#'
#' @description Simulate data for two-sided matching markets. In the simulation for the 
#' Sorensen (2007) model with one selection equation, an equal sharing rule of 
#' \eqn{\lambda = 0.5} is used.
#'
#' @param m integer indicating the number of markets to be simulated.
#' @param nStudents integer indicating the number of students per market.
#' @param nColleges integer indicating the number of colleges per market.
#' @param nSlots vector of length \code{nColleges} indicating the number of places at each college, i.e. 
#' the college's quota.
#' @param outcome formula for match outcomes.
#' @param selection formula for match valuations.
# @param selection.college formula for match valuations of colleges. Is ignored when \code{selection} is provided.
# @param selection.student formula for match valuations of students. Is ignored when \code{selection} is provided.
#' @param colleges character vector of variable names for college characteristics. These variables carry the same value for any college.
#' @param students character vector of variable names for student characteristics. These variables carry the same value for any student.
#' @param binary logical: if \code{TRUE} outcome variable is binary; if \code{FALSE} outcome variable is continuous.
#' @param seed integer setting the state for random number generation. Defaults to \code{set.seed(123)}.
#' @param verbose .
#' 
#' @export
#' 
#' @import stats
#' @importFrom utils setTxtProgressBar txtProgressBar
#' 
#' @return
#' 
#' @return
#' \code{stabsim2} returns a list with the following items.
#' \item{OUT}{}
#' \item{SEL}{}
#' \item{SELc}{}
#' \item{SELs}{}
#' 
#' @author Thilo Klein 
#' 
#' @keywords generate
#' 
#' @examples
#' 
#' ## Simulate two-sided matching data for 2 markets (m=2) with 10 students
#' ## (nStudents=10) per market and 3 colleges (nColleges=3) with quotas of
#' ## 2, 3, and 5 students, respectively.
#' 
#' xdata <- stabsim2(m=2, nStudents=10, nSlots=c(2,3,5), verbose=FALSE,
#'   colleges = "c1", students = "s1",
#'   outcome = ~ c1:s1 + eta + nu,
#'   selection = ~ -1 + c1:s1 + eta
#' )
#' head(xdata$OUT)
#' head(xdata$SEL)
#' 
stabsim2 <- function(m, nStudents, nColleges=length(nSlots), nSlots, 
                     colleges, students, outcome, selection, 
                     binary=FALSE, seed=123, verbose=TRUE){

  set.seed(seed)
  
  ## select method based on arguments provided
  if(is.list(selection)){
    method <- "Klein"
    selection.college <- selection$college
    selection.student <- selection$student
    selection <- NULL
  } else{
    method <- "Sorensen"
  }
  
  # --------------------------------------------------------------------
  # R-code (www.r-project.org) for simulating purely random (!) data for
  # all players in the market.
  
  stabsim2_inner <- function(mi, nStudents, nColleges=length(nSlots), nSlots, colleges, students, 
                             outcome, selection, selection.student, selection.college){

    ## unique student and college ids
    uColleges <- 1:nColleges
    uStudents <- 1:nStudents
    
    ## all feasible combinations
    combs <- function(uColleges, uStudents){
      data.frame( c.id = c(sapply(uColleges, function(i){ rep(i, nStudents) })), 
                  s.id = rep(uStudents, nColleges), 
                  stringsAsFactors=FALSE )
    }
    indices <- as.data.frame(combs(uColleges, uStudents))
    
    ## create random, exogenous variables
    C <- data.frame(matrix(rnorm(nColleges*length(colleges), sd=sqrt(0.5)), nrow=nColleges, ncol=length(colleges)))
    names(C) <- colleges
    
    S <- data.frame(matrix(rnorm(nStudents*length(students), sd=sqrt(0.5)), nrow=nStudents, ncol=length(students)))
    names(S) <- students
    
    ## Main and interaction effects
    Cexp <- as.data.frame(apply(C, 2, function(i) i[indices$c.id]))
    Sexp <- as.data.frame(apply(S, 2, function(i) i[indices$s.id]))
    Xexp <- cbind(Cexp, Sexp)
        
    
    ## ---- Match valuations ----    
    
    if(method == "Klein"){
      
      student.terms <- attr( attr(terms.formula(selection.student), "factors"), "dimnames")[[1]]
      student.terms <- student.terms[! student.terms %in% c(colleges, students)]
      college.terms <- attr( attr(terms.formula(selection.college), "factors"), "dimnames")[[1]]
      college.terms <- college.terms[! college.terms %in% c(colleges, students)]
      outcome.terms <- attr( attr(terms.formula(outcome), "factors"), "dimnames")[[1]]
      outcome.terms <- outcome.terms[! outcome.terms %in% c(colleges, students, college.terms, student.terms)]
      terms <- c(college.terms, student.terms, outcome.terms)
      Mexp <- data.frame(matrix(rnorm(length(terms)*nrow(Xexp), sd=1), ncol=length(terms))) ## !!! sd.
      names(Mexp) <- terms
      
      if(dim(Mexp)[1] >0 ){
        Xexp <- cbind(Xexp, Mexp)
      }
      
      Smtch <- rFormula(formula = selection.student, data=Xexp)
      Cmtch <- rFormula(formula = selection.college, data=Xexp)
      Xmtch <- rFormula(formula = outcome, data=Xexp)
      
      Vc <- apply(Cmtch, 1, sum)
      Vs <- apply(Smtch, 1, sum)
      
    } else if(method == "Sorensen"){

      selection.terms <- attr( attr(terms.formula(selection), "factors"), "dimnames")[[1]]
      selection.terms <- selection.terms[! selection.terms %in% c(colleges, students)]
      outcome.terms <- attr( attr(terms.formula(outcome), "factors"), "dimnames")[[1]]
      outcome.terms <- outcome.terms[! outcome.terms %in% c(colleges, students, selection.terms)]
      terms <- c(selection.terms, outcome.terms)
      Mexp <- data.frame(matrix(rnorm(length(terms)*nrow(Xexp), sd=1), ncol=length(terms))) ## !!! sd.
      names(Mexp) <- terms
      
      if(dim(Mexp)[1] >0 ){
        Xexp <- cbind(Xexp, Mexp)
      }
      
      Wmtch <- rFormula(formula = selection, data=Xexp)
      Xmtch <- rFormula(formula = outcome, data=Xexp)
      
      Vc <- Vs <- apply(Wmtch, 1, sum)*0.5 ## equal sharing rule
    }
    
    ## convert preferences to ranks in matrix format
    c.prefs <- matrix(Vs, ncol=nColleges, nrow=nStudents, byrow=FALSE)
    s.prefs <- matrix(Vc, ncol=nStudents, nrow=nColleges, byrow=TRUE)
    
    c.prefs <- apply(-1*c.prefs, 2, order)
    s.prefs <- apply(-1*s.prefs, 2, order)
    
    ## obtain student-optimal matching
    matching <- hri(s.prefs=s.prefs, c.prefs=c.prefs, nSlots=nSlots)$matchings
    matching <- matching[matching$sOptimal==1,]
    
    ## obtain equilibrium identifier 'd'
    matching$id <- paste(matching$college, matching$student, sep="_")
    indices$id  <- paste(indices$c.id, indices$s.id, sep="_")
    indices$D <- ifelse(indices$id %in% matching$id, 1, 0) # identifier for s-optimal matches
    
    d <- which(indices$id %in% matching$id)
    indx <- c(d, (1:nrow(indices))[-d]) # s-optimal matchings come first in output
    indices$id <- NULL
    
    ## --- OUTPUT for stabsim2() ---
    
    Xmatch <- Xmtch[d,]
    y <- apply(Xmatch, 1, sum)
    OUT <- cbind(m.id=mi, y, Xmatch, Xexp[d,setdiff(names(Xexp),names(Xmtch))], c.id=indices$c.id[d], s.id=indices$s.id[d]) ## add instrumental variables that only enter selection eqn
    
    if(method == "Klein"){
      
      ## for colleges: add main effects for all interation terms
      college.terms <- attr( attr(terms.formula(selection.college), "factors"), "dimnames")[[1]]
      college.terms <- college.terms[! college.terms %in% names(Cmtch)]
      addvars <- as.data.frame(Xexp[indx,college.terms]); names(addvars) <- college.terms
      
      SELc <- cbind(m.id=mi, Vc=Vc[indx], Cmtch[indx,], addvars, indices[indx,])
      
      ## for students: add main effects for all interation terms
      student.terms <- attr( attr(terms.formula(selection.student), "factors"), "dimnames")[[1]]
      student.terms <- student.terms[! student.terms %in% names(Smtch)]
      addvars <- as.data.frame(Xexp[indx,student.terms]); names(addvars) <- student.terms
      
      SELs <- cbind(m.id=mi, Vs=Vs[indx], Smtch[indx,], addvars, indices[indx,])
      
      list(OUT=OUT, SELc=SELc, SELs=SELs)
            
    } else if(method == "Sorensen"){
      
      selection.terms <- attr( attr(terms.formula(selection), "factors"), "dimnames")[[1]]
      selection.terms <- selection.terms[! selection.terms %in% names(Wmtch)]
      addvars <- as.data.frame(Xexp[indx,selection.terms]); names(addvars) <- selection.terms
      
      SEL <- cbind(m.id=mi, V=Vc[indx]+Vs[indx], Wmtch[indx,], addvars, indices[indx,]) #Xexp[indx,setdiff(names(Xexp),names(Xmtch))]
      
      list(OUT=OUT, SEL=SEL)
    }
  }
  
  if(method == "Klein"){
    
    RETURN <- list(OUT=list(), SELs=list(), SELc=list())
    
    cat("Generating data for", m, "matching markets...","\n")
    if(verbose==TRUE){
      pb <- txtProgressBar(style = 3)
    }
    
    for(i in 1:m){
      X <- stabsim2_inner(mi=i, nStudents=nStudents, nSlots=nSlots, 
                          colleges=colleges, students=students, outcome=outcome, selection=selection,
                          selection.student=selection.student, selection.college=selection.college)  
      RETURN$OUT[[i]]  <- X$OUT
      RETURN$SELs[[i]] <- X$SELs
      RETURN$SELc[[i]] <- X$SELc
      if(verbose==TRUE){
        setTxtProgressBar(pb, i/m)
      }
    }
       
  } else if(method == "Sorensen"){
    
    RETURN <- list(OUT=list(), SEL=list())  
    
    cat("Generating data for", m, "matching markets...","\n")
    if(verbose==TRUE){
      pb <- txtProgressBar(style = 3)
    }
    
    for(i in 1:m){
      X <- stabsim2_inner(mi=i, nStudents=nStudents, nSlots=nSlots, 
                          colleges=colleges, students=students, outcome=outcome, selection=selection,
                          selection.student=selection.student, selection.college=selection.college)  
      RETURN$OUT[[i]] <- X$OUT
      RETURN$SEL[[i]] <- X$SEL
      if(verbose==TRUE){
        setTxtProgressBar(pb, i/m)
      }
    }
  }
  
  RETURN <- lapply(RETURN, function(i){
    h <- do.call(rbind,i)
    rownames(h) <- 1:dim(h)[1]
    h
  })
  
  if(binary == TRUE){
    RETURN$OUT$y <- ifelse(RETURN$OUT$y > median(RETURN$OUT$y), 1, 0)
  }
  
  return(RETURN)
}




rFormula <- function(formula, data=list(), ...){
  mf <- model.frame(formula=formula, data=data)
  x <- model.matrix(attr(mf, "terms"), data=mf)
  y <- model.response(mf)
  as.data.frame(cbind(y,x))
}

