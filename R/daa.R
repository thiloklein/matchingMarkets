# ----------------------------------------------------------------------------
# R-code (www.r-project.org/) for the Deferred Acceptance Algorithm
#
# Copyright (c) 2013 Thilo Klein
#
# This library is distributed under the terms of the GNU Public License (GPL)
# for full details see the file LICENSE
#
# ----------------------------------------------------------------------------

#' @title Deferred Acceptance Algorithm
#'
#' @description Finds the student (men) optimal matching in the
#' \href{http://en.wikipedia.org/wiki/Hospital_resident}{college admissions} 
#' (\href{http://en.wikipedia.org/wiki/Stable_matching}{stable marriage}) problem. 
#' Uses the Gale-Shapley (1962) Deferred Acceptance Algorithm with student (male) 
#' offer based on given or randomly generated preferences.
#'
#' @param nStudents integer indicating the number of students (in the college admissions problem) or men (in the stable marriage problem) in the market
#' @param nColleges integer indicating the number of colleges (in the college admissions problem) or women (in the stable marriage problem) in the market
#' @param nSlots vector of length \code{nColleges} indicating the number of places (i.e. 
#' quota) of each college. Defaults to \code{rep(1,nColleges)} for the marriage problem.
#' @param s.prefs matrix of dimension \code{nColleges x nStudents} with the \code{i}th 
#' column containing student \code{i}'s ranking over colleges in decreasing order of 
#' preference (i.e. most preferred first)
#' @param c.prefs matrix of dimension \code{nStudents x nColleges} with the \code{j}th 
#' column containing college \code{j}'s rankimatg over students in decreasing order of 
#' preference (i.e. most preferred first) 
#' @export
#' @section Minimum required arguments:
#' Subject to the matching problem, 'daa' requires the following arguments.
#' Marriage problem with random preferences: \code{nStudents, nColleges}.
#' Marriage problem with given preferences: \code{s.prefs, c.prefs}.
#' College admissions problem with random preferences: \code{nStudents, nSlots}.
#' College admissions problem with given preferences: \code{s.prefs, c.prefs, nSlots}.
#' @section Value: 
#' 'daa' returns a list with the following items.
#' \describe{
#' \item{\code{s.prefs}}{students' preference matrix.}
#' \item{\code{c.prefs}}{colleges' preference matrix.}
#' \item{\code{iterations}}{number of interations required to find the stable matching.}
#' \item{\code{matches}}{identifier of students (men) assigned to colleges (women).}
#' \item{\code{match.mat}}{matching matrix of dimension \code{nStudents x nColleges}.}
#' \item{\code{singles}}{identifier of single/unmatched students (men).}
#' }
#' @author Thilo Klein <\email{thilo@@klein.co.uk}>
#' @references Gale, D. and Shapley, L.S. (1962). College admissions and the stability 
#' of marriage. The American Mathematical Monthly, 69(1):9--15.
#' @examples
#' ## Marriage problem (3 men, 2 women) with random preferences:
#' daa(nStudents=3, nColleges=2)
#' 
#' ## Marriage problem (3 men, 2 women) with given preferences:
#' s.prefs <- matrix(c(1,2, 1,2, 1,2), 2,3)
#' c.prefs <- matrix(c(1,2,3, 1,2,3), 3,2)
#' daa(s.prefs=s.prefs, c.prefs=c.prefs)
#' 
#' ## College admission problem (7 students, 2 colleges 
#' ## with 3 slots each) with random preferences:
#' daa(nStudents=7, nSlots=c(3,3))
#' 
#' ## College admission problem (7 students, 2 colleges 
#' ## with 3 slots each) with given preferences:
#' s.prefs <- matrix(c(1,2, 1,2, 1,2, 1,2, 1,2, 1,2, 1,2), 2,7)
#' c.prefs <- matrix(c(1,2,3,4,5,6,7, 1,2,3,4,5,6,7), 7,2)
#' daa(s.prefs=s.prefs, c.prefs=c.prefs, nSlots=c(3,3))
daa <- function(nStudents=dim(s.prefs)[2], nColleges=dim(c.prefs)[2], nSlots=rep(1,nColleges), s.prefs=NULL, c.prefs=NULL){

  ## If 'nColleges' not given, obtain it from nSlots
  if (is.null(nColleges)){
    nColleges <- length(nSlots)
  }
  ## If no prefs given, make them randomly:
  if (is.null(s.prefs)){  
    s.prefs <- replicate(n=nStudents,sample(seq(from=1,to=nColleges,by=1)))
  }
  if (is.null(c.prefs)){    
    c.prefs <- replicate(n=nColleges,sample(seq(from=1,to=nStudents,by=1)))
  }
  
  ## Consistency checks:
  if ( dim(s.prefs)[1] != dim(c.prefs)[2] | dim(s.prefs)[2] != dim(c.prefs)[1] | 
       dim(s.prefs)[2] != nStudents | dim(c.prefs)[2] != nColleges | 
       dim(c.prefs)[1] != nStudents | dim(s.prefs)[1] != nColleges ){
    stop("'s.prefs' and 'c.prefs' must be of dimensions 'nColleges x nStudents' and 'nStudents x nColleges'!")
  }
  if ( length(nSlots) != nColleges | length(nSlots) != dim(c.prefs)[2] ){
    stop("length of 'nSlots' must equal 'nColleges' and the number of columns of 'c.prefs'!")
  }

  iter = 0

  s.hist    <- rep(0,length=nStudents)  # number of proposals made
  c.hist    <- lapply(nSlots, function(x) rep(0,length=x))  # current students
  s.singles <- 1:nStudents

  s.mat <- matrix(data=1:nStudents,nrow=nStudents,ncol=nColleges,byrow=F)

  while (max(s.hist) < nColleges){  # there are as many rounds as maximal preference orders
    # look at market: all unassigned students
    # if history not full (been rejected by all colleges in his prefs)
    # look at unassigned students' history
    # propose to next college on list
    iter = iter + 1
    print(paste("Iteration: ",iter))
    offers <- NULL

    for (i in 1:length(s.singles)){
      s.hist[s.singles[i]] <- s.hist[s.singles[i]] + 1  # set history of student i one up.
    offers[i] <- s.prefs[s.hist[s.singles[i]],s.singles[i]]  # offer if unassigned i is index of current round college
	}

	approached   <- unique(offers)	# index of colleges who received offers
	temp.singles <- s.singles
	s.singles    <- NULL	        # reset unassigned students
	
	for (j in approached){
	  proposers   <- temp.singles[offers==j]
	  stay.single <- temp.singles[offers==0]		# students who prefer remaining unassigned at current history

	  for (k in 1:length(proposers)){
	    if (0 %in% (c.hist[[j]] & any(c.prefs[ ,j]==proposers[k]))){  # if no history and proposer is on preference list
		  c.hist[[j]][c.hist[[j]]==0][1] <- proposers[k]			  # then accept
		} else if (TRUE %in% (match(c.prefs[c.prefs[ ,j]==proposers[k],j],c.prefs[ ,j]) < match(c.prefs[c.prefs[ ,j] %in% c.hist[[j]], j],c.prefs[ ,j]))){   # if proposer better than any current student
		  worst = max(match(c.prefs[c.prefs[ ,j] %in% c.hist[[j]], j], c.prefs[ ,j])) # determine worst current student
		  s.singles <- c(s.singles,c.prefs[worst,j])   # reject worst current student
		  c.hist[[j]][c.hist[[j]] == c.prefs[worst,j]] <- proposers[k]	# and take proposer on
		} else {
		  s.singles <- c(s.singles,proposers[k])	# otherwise k stays unassigned
		}
	  }	
    }

    s.singles <- sort(c(s.singles,stay.single))
	if (length(s.singles)==0){	# if no unassigned students left: stop
	  current.match   <- sapply(1:nColleges, function(x) s.mat[,x] %in% c.hist[[x]])
	  return(list(s.prefs=s.prefs,c.prefs=c.prefs,iterations=iter,matches=c.hist,match.mat=current.match,singles=s.singles))
  	  break
	}
	current.match   <- sapply(1:nColleges, function(x) s.mat[,x] %in% c.hist[[x]])
  }

  return(list(s.prefs=s.prefs,c.prefs=c.prefs,iterations=iter,matches=c.hist,match.mat=current.match,singles=s.singles))
}