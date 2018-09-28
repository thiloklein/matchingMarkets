# ----------------------------------------------------------------------------
# Copyright (c) 2018 Sven Giegerich, Thilo Klein
#
# This library is distributed under the terms of the GNU Public License (GPL)
# for full details see the file LICENSE
#
# ----------------------------------------------------------------------------
#
#' @title Resident-optimal matching in the hospital/residents problem with couples
#' @description Implements the Roth Peranson matching algorithm for the \href{https://en.wikipedia.org/wiki/National_Resident_Matching_Program}{hospital/residents problem with couples} as described in Roth and Peranson (1999). The function is based on an adoption of Bacchus (2018). 
#' @param nStudents integer indicating the number of students (in the college admissions problem) 
#' or men (in the stable marriage problem) in the market. Defaults to \code{ncol(s.prefs)}.
#' @param nColleges integer indicating the number of colleges (in the college admissions problem) 
#' or women (in the stable marriage problem) in the market. Defaults to \code{ncol(c.prefs)}.
#' @param nCouples  integer indicating the number of couples (in the college admissions problem) 
#' or men (in the stable marriage problem) in the market. Defaults to \code{ncol(co.prefs)}
#' @param nSlots vector of length \code{nColleges} indicating the number of places (i.e. 
#' quota) of each college. Defaults to \code{rep(1,nColleges)} for the marriage problem.
#' @param s.prefs matrix of dimension \code{nColleges} \code{x} \code{nStudents} with the \code{j}th 
#' column containing student \code{j}'s ranking over colleges in decreasing order of 
#' preference (i.e. most preferred first).
#' @param c.prefs matrix of dimension \code{nStudents} \code{x} \code{nColleges} with the \code{i}th 
#' column containing college \code{i}'s ranking over students in decreasing order of 
#' preference (i.e. most preferred first).
#' @param co.prefs matrix of dimension \code{4} \code{x} \code{nCouplesPrefs} in long format with the \code{1}th and \code{2}th
#' columns containing student couple id's; \code{3}th and \code{4}th is a 2-tuple ranking over college preference for the couple (coupleStudent1.pref, coupleStudent2.pref) in decreasing order of 
#' preference by rows (i.e. most preferred first).
#' @param randomization determines at which level and in which order random lottery numbers for student priorities are drawn. The default is \code{randomization = "multiple"}, where a student's priority is determined by a separate lottery at each college (i.e. local tie-breaking). For the second variant, \code{randomization = "single"}, a single lottery number determines a student's priority at all colleges (i.e. global tie breaking). A third variant is common in the context of course allocation, where a "couple" represents a student who submits a preference ranking over single courses (first course) and combinations of courses (first and second course). Here, the option \code{randomization = "single-course-first"} gives applications for a student's single courses strictly higher priority than for course combinations. This ensures the fairness criterion that a student is only assigned a second course after single course applications of all students have been considered.
#' @param seed integer setting the state for random number generation. 
#' @param ... .
#' 
#' @export
#' 
#' @useDynLib matchingMarkets, .registration = TRUE 
#' 
#' @importFrom Rcpp evalCpp
#' 
#' @section Minimum required arguments:
#' \code{hri2} requires the following combination of arguments, subject to the matching problem.
#' \describe{
#' \item{\code{nStudents, nColleges}}{Residence hospital problem without couples and random preferences}
#' \item{\code{nStudents, nColleges, nCouples, nSlots}}{Residence hospital problem with couples and random preferences.}
#' \item{\code{s.prefs, c.prefs, co.prefs, nSlots}}{Residence hospital problem with couples and given preferences.}
#' }
#' 
#' @return
#' \code{hri2} returns a list of the following elements:
#' \item{matchings}{List of matched students and colleges.}
#' \item{summary}{Detailed report of the matching result, including futher information on ranks.}
#' 
#' @author Sven Giegerich, Thilo Klein 
#' 
#' @keywords algorithms, matching
#' 
#' @references Bacchus, F. (2018). Stable matching suite. GitHub repository.
#' 
#' Gale, D. and L.S. Shapley (1962). College admissions and the stability 
#' of marriage. \emph{The American Mathematical Monthly}, 69(1):9--15.
#' 
#' Roth, A. E., & Peranson, E. (1999). The redesign of the matching market for American physicians: Some engineering aspects of economic design. \emph{American economic review}, 89(4), 748-780.
#' 
#' Kojima, F., Pathak, P. A., & Roth, A. E. (2013). Matching with couples: Stability and incentives in large markets. \emph{The Quarterly Journal of Economics}, 128(4), 1585-1632.
#' 
#' @examples
#' ## Example with given preferences
#' (s.prefs <- matrix(c(4,2,3,5, 2,1,3,NA, 1,2,3,4), 4,3))
#' (c.prefs <- matrix(rep(1:5,5), 5,5))
#' (co.prefs <- matrix(c(rep(4,3), rep(5,3), 3,3,NA, 3,NA,3), 3,4))
#' res <- hri2(s.prefs=s.prefs, c.prefs=c.prefs, co.prefs=co.prefs, nSlots=rep(1,5))
#' res$matchings
#' # summary(res)
#' 
#' ## Example with random preferences
#' nStudents <- 50
#' nColleges <- 30
#' nCouples <- 4
#' nSlots <- sample(1:nStudents, nColleges)
#' res <- hri2(nStudents=nStudents, nColleges=nColleges, nCouples=nCouples, nSlots=nSlots)
#' res$matchings
#' # summary(res)

hri2 <- function(nStudents=ncol(s.prefs), nColleges=ncol(c.prefs), nSlots=rep(1,nColleges), nCouples=ncol(co.prefs), 
                  s.prefs=NULL, c.prefs=NULL, co.prefs=NULL, randomization="multiple", seed=NULL, ...) UseMethod("hri2")

#' @export
hri2.default <- function(nStudents=ncol(s.prefs), nColleges=ncol(c.prefs), nSlots=rep(1,nColleges), nCouples=ncol(co.prefs), 
                        s.prefs=NULL, c.prefs=NULL, co.prefs=NULL, randomization="multiple", seed=NULL, ...){
  
  ## -------------------------------------------------------
  ## --- 1. consistency checks:  ---------------------------
  if (missing(nStudents) && missing(s.prefs)) {
    stop("Required arguments are missing")
  }
  
  if (missing(nColleges) && missing(c.prefs) && missing(nSlots)) {
    stop("Required arguments are missing")
  }
  
  ## ------------------------
  ## --- 2. Preliminaries ---
  
  ## set seed for random preference draws
  if(!is.null(seed)){
    set.seed(seed)
  }
  
  ## if 'nColleges' not given, obtain it from nSlots
  if(is.null(nColleges)){
    nColleges <- length(nSlots)
  }
  
  if(is.null(nCouples) && !is.null(co.prefs)){
    nCouples <- length(co.prefs)
  } else if (is.null(nCouples) && is.null(co.prefs)) {
    message("normal hri problem without couples")
    nCouples <- 0
  }
  
  ## if no resident prefs given, generate random ranking:
  if(is.null(s.prefs)){  
    s.prefs <- replicate(n=nStudents,sample(seq(from=1,to=nColleges,by=1)))
  }
  if(is.null(c.prefs)){  
    if(randomization == "single"){ 
      
      c.prefs <- matrix(sample(seq(from=1, to=nStudents+2*nCouples, by=1)), nrow=nStudents+2*nCouples, ncol=nColleges) 

    } else if(randomization == "single-course-first"){
      
      c.prefs <- matrix( c(sample(c(1:nStudents, seq(nStudents+1, nStudents+2*nCouples, by=2))), 
                 sample(seq(nStudents+2, nStudents+2*nCouples, by=2))),
                 nrow=nStudents+2*nCouples, ncol=nColleges)
      
    } else{ # if(randomization == "multiple")
      
      c.prefs <- replicate(n=nColleges ,sample(seq(from=1, to=nStudents+2*nCouples, by=1))) 
    }
  }
  if(is.null(co.prefs) && nCouples > 0){
    co.prefs <- matrix(ncol = nCouples, nrow = 2+2*nColleges)
    co <- 1
    for (i in 1:nCouples) {
      co.prefs[1,i] <- nStudents+co
      co.prefs[2, i] <- nStudents+co+1
      co.prefs[3:(2+2*nColleges), i] <- c(sample(seq(from=1,to=nColleges,by=1)), sample(seq(from=1,to=nColleges,by=1)))
      co <- co+2
    }
    co.prefs <- t(co.prefs)
  }
  
  ##co.prefs in wide or long format? -> to wide
  if (!missing(co.prefs) && is.matrix(co.prefs) && dim(co.prefs)[2] == 4) {
    co.prefs <- copL2copW(co.prefs)
  }
  
  ##missing prefernces need to be tranlated to -1
  if (!missing(co.prefs) && is.matrix(co.prefs)) {
    co.prefs[is.na(co.prefs)] <- 0
  }
  
  ## -------------------------------------------------------
  ## --- 3. Prepare preference matrices and apply solver ---
  
  ## prepare and write preference matrices
  c.matrix <- sapply(1:ncol(c.prefs), function(z) paste("p", z - 1, if(is.na(nSlots[z])){"0"} else{nSlots[z]}, paste( t(c.prefs)[z,][!is.na(t(c.prefs)[z,])] -1, collapse = " ")))
  s.matrix <- sapply(1:ncol(s.prefs), function(z) paste("r", z - 1, paste( t(s.prefs)[z,][!is.na(t(s.prefs)[z,])] -1, collapse = " ")))
  if (nCouples > 0) {
    co.matrix <- sapply(1:nrow(co.prefs), function(z) paste("c", z - 1, paste(co.prefs[z,][!is.na(co.prefs[z,])] -1, collapse = " ")))
  } else {
    co.matrix <- c("")
  }
  
  matchResult <- runMatch(s.matrix, c.matrix, co.matrix)
  
  matchResult$matchings <- cbind(matchResult$matchings$ResidentID, matchResult$matchings$matchResultResident)
  colnames(matchResult$matchings) <- c("student", "college")
  
  ## drop unmatched students and colleges
  matchResult$matchings <- matchResult$matchings[(matchResult$matchings[,1] != 0) & (matchResult$matchings[,2] != 0),]
  
  class(matchResult) <- "hri2"
  return(matchResult)
}

##Couples pref from long format to wide
copL2copW <- function(couplesPref) {
  #Example
  ##couplesPref <- matrix(c(5,5,5, 6,6,6, 1,1,0, 1,0,1), 3,4)
  ##couplesPref <- data.frame(c(5,5,5,7), c( 6,6,6,8), c( 1,1,0,1), c(1,0,1,1))
  ##5 6 1 1
  ##5 6 1 0
  ##5 6 0 1
  ##7 8 1 1
  ##7 8 1 0
  ##7 8 0 1
  ## -> to
  # 5 6 1 1 0 1 1 0
  # 7 8 1 1 1 0 0 1
  couplesPref <- as.data.frame(couplesPref)
  couplesPref.split <- split(couplesPref, list(couplesPref[,1], couplesPref[,2]), drop = TRUE)
  
  spreadCouplesPrefs <- function(co.prefs, x, cID) {
    x <- data.frame(x)
    co.prefs[cID, 1] <- x[1,1]
    co.prefs[cID, 2] <- x[1,2]
    c <- 3
    for ( i in 1:(dim(x)[1]) ) {
      co.prefs[cID, c] <- x[i, 3]
      co.prefs[cID, c+1] <- x[i, 4]
      c <- c+2
    }
    return(co.prefs)
  }
  
  co.prefs <- data.frame()
  coupleID <- 1
  for (j in 1:length(couplesPref.split)) {
    co.prefs <- spreadCouplesPrefs(co.prefs, couplesPref.split[j], coupleID)
    coupleID <- coupleID+1
  }
  return(as.matrix(co.prefs))
}