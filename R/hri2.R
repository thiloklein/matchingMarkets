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
#' @param check_consistency Performs additional consicentcy checks if the preference matrices are given by characters. Defaults to \code{FALSE}. Set to \code{FALSE} to reduce run-time.
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
#' @author Sven Giegerich, Thilo Klein, Alexander Sauer
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
#' 
#' ## Example with characters in the preferences matrices
#' s.prefs <- matrix(c("Micro1", NA, NA,
#'                     "Micro2", "Micro1", "Macro",
#'                     "Macro",NA ,NA), 
#'                     ncol = 3)
#' colnames(s.prefs) <- c('Lea', 'Mia', 'Kai')
#' c.prefs <- matrix(c("Niklas", "Kai", "Mia", "Anna",
#'                     "Lea", "Kai", "Anna",NA,
#'                     "Kai", "Mia", "Lea",NA), 
#'                     ncol = 3)
#' colnames(c.prefs) <- c('Micro1', 'Micro2', 'Macro')
#' col1 <- c(rep("Niklas",4),rep("Anna",5))
#' col2 <- c(rep("Jan",4),rep("Lisa",5))
#' col3 <- c("Micro1","Macro","Micro1",NA,"Macro",
#'           NA,"Micro2","Micro2","Macro")
#' col4 <- c("Micro2","Micro1",NA,"Macro","Macro",
#'           "Micro1","Micro2","Macro",NA)
#' co.prefs <- matrix(c(col1,col2,col3,col4), ncol = 4)
#' res <- hri2(s.prefs=s.prefs, c.prefs=c.prefs, co.prefs=co.prefs, 
#'             nSlots=c(2,1,1))                     
#' res$matching
#' 
#' ## Example if students are allowed to apply and be accepted by two courses   
#' col12 <- c(rep(c(rep("Niklas",4),rep("Anna",2)),2))
#' col3 <- c("Micro1","Macro","Micro1","Macro","Macro","Macro")
#' col4 <- c("Micro2","Micro1",NA,NA,"Micro1","Micro2")
#' co.prefs <- matrix(c(col12,col3,col4), ncol = 4)
#' res <- hri2(s.prefs=s.prefs, c.prefs=c.prefs, co.prefs=co.prefs, 
#'             nSlots=c(2,1,1))                     
#' res$matching                                                                                     


hri2 <- function(nStudents=ncol(s.prefs), nColleges=ncol(c.prefs), nSlots=rep(1,nColleges), nCouples=ncol(co.prefs), 
                 s.prefs=NULL, c.prefs=NULL, co.prefs=NULL, randomization="multiple", seed=NULL, check_consistency=TRUE, ...) UseMethod("hri2")

#' @export
hri2.default <- function(nStudents=ncol(s.prefs), nColleges=ncol(c.prefs), nSlots=rep(1,nColleges), nCouples=ncol(co.prefs), 
                         s.prefs=NULL, c.prefs=NULL, co.prefs=NULL, randomization="multiple", seed=NULL, check_consistency=TRUE,...){
  
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
  
  ## Check if preferences are given with 'true'/character names
  prefs_as_char <- FALSE
  if(is.character(s.prefs) || is.character(c.prefs)  || is.character(co.prefs)){
    if(!(is.character(s.prefs) && is.character(c.prefs) && is.character(co.prefs))) {
      stop('Prefs must all be as characters or all numeric')
    }
    prefs_as_char <- TRUE
    
    # Check for consistency
    if(check_consistency){consistency_check(s.prefs, c.prefs, co.prefs)}
    
    #Check if the co.prefs matrix represents real couples or preferences over two subjects
    kurs_pref <- all(co.prefs[,1] == co.prefs[,2])
    
    # All student names (single + couples)
    # If the couples preferences are no real couples, the name of the person is duplicated
    if(kurs_pref){
      s.names <- 1:(ncol(s.prefs) + 2*length(unique(as.character(co.prefs[,c(1,2)]))))
      names(s.names) <- c(colnames(s.prefs), as.character(vapply(unique(as.character(co.prefs[,c(1,2)])), rep, FUN.VALUE = character(2), times = 2)))
    } else{
      s.names <- 1:(ncol(s.prefs) + length(unique(as.character(co.prefs[,c(1,2)]))))
      names(s.names) <- c(colnames(s.prefs), unique(as.character(co.prefs[,c(1,2)])))
    }
    
    c.names <- 1:ncol(c.prefs)
    names(c.names) <- colnames(c.prefs)
    
    ## Student prefs
    # Replace names with identifiers in preferences matrices
    s.prefs <- apply(s.prefs, 2, function(pref){
      return(c.names[pref])
    })
    colnames(s.prefs) <- NULL
    rownames(s.prefs) <- NULL
    
    ## College prefs
    # If the couples prefs are not real couples, the preferences of the college over these people have to include both IDs
    if(kurs_pref){
      matched_pref <- apply(c.prefs,2, function(pref){
        return(unlist(sapply(pref, function(e){which(names(s.names) == e)}, USE.NAMES = FALSE)))  # Match names to IDs
      })
      # Transform list to matrix:
      max_n_row <- max(sapply(matched_pref, length))
      c.prefs <- sapply(matched_pref, function(x) {
        length(x) <- max_n_row
        return(x)
        })
    }else{
      c.prefs <-  apply(c.prefs, 2, function(pref){
        return(s.names[pref])
      })
    }
    colnames(c.prefs) <- NULL
    rownames(c.prefs) <- NULL
    
    ## Couples prefs
    co.prefs <- sapply(1:4,  function(col){
      pref <- co.prefs[,col]
      if(col <= 2){
        if(kurs_pref && col == 2){ # If the co.prefs correspond to one individual, then use that the second column is the same individual with the ID + 1
          return(s.names[pref] + 1)
        }
        return(s.names[pref])
      }
      else{
        return(c.names[pref])
      }
    })
    colnames(co.prefs) <- NULL
    rownames(co.prefs) <- NULL
  }
  ######## Transformation of characters prefs finished.
  
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
  
  # Match back identifiers to string-names
  if(prefs_as_char){
    # matchResult$matchings$student <- names(s.names)[matchResult$matchings$student]
    # matchResult$matchings$college <- names(c.names)[matchResult$matchings$college]
    # First column represents students/ second colleges
    matchResult$matchings[,1] <- names(s.names)[as.numeric(matchResult$matchings[,1])]
    matchResult$matchings[,2] <- names(c.names)[as.numeric(matchResult$matchings[,2])]
  }
  
  
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


consistency_check <- function(s.prefs, c.prefs, co.prefs){
  # Check if student names are unique:
  if(length(unique(colnames(s.prefs))) != ncol(s.prefs)) {stop('Student names not unique')}
  if(length(intersect(colnames(s.prefs), unique(as.character(co.prefs[,c(1,2)])))) != 0) { stop('Student appears in s.prefs and co.prefs')}
  
  # Check if colleges are unique:
  if(length(unique(colnames(c.prefs))) != ncol(c.prefs)) {stop('College/Course names not unique')}
  
  # Check if a student applied for a course/college that is not in c.prefs
  applied_colleges <- unique(c(as.character(s.prefs), as.character(co.prefs[,c(3,4)])))
  applied_colleges <- applied_colleges[!is.na(applied_colleges)]
  
  if(length(setdiff(applied_colleges,colnames(c.prefs))) != 0){
    missing_college <- setdiff(applied_colleges,colnames(c.prefs))
    missing_college <- paste('Someone applied to a college (', missing_college, ') that has no ranking', sep = '')
    stop(missing_college)
  }
  
  
  # Check if a college ranked someone who is not in s.prefs
  ranked_stud <- unique(c(as.character(c.prefs)))
  ranked_stud <- ranked_stud[!is.na(ranked_stud)]
  
  if(length(setdiff(ranked_stud,c(colnames(s.prefs), unique(as.character(co.prefs[,c(1,2)]))))) != 0){
    missing_stud <- setdiff(ranked_stud,c(colnames(s.prefs), unique(as.character(co.prefs[,c(1,2)]))))
    missing_stud <- paste('A course/college ranked to a student (', missing_stud, ') that has no ranking', sep = '')
    stop(missing_stud)
  }
  
  # Check if co.prefs are real couples or individuals how rank two subjects, but not a mix!
  if(!all(co.prefs[,1] == co.prefs[,2])){
    # Check if co1 and co2 are different in all rows
    if( !all(apply(co.prefs, 1, function(row){return(row[1] != row[2]) })
    )) { stop('Some rows in co.prefs have the same student as participant one and two.')}
  }
  print('Input passed consistency test!')
}


