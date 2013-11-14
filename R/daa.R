#' Deferred Acceptance Algorithm
#'
#' This function implements the student-proposing version of the Gale-Shapley algorithm. It takes student and college preference matrices and returns the unique stable assignment.
#'
#' @param nStudents an integer indicating the number of students in the market
#' @param nSlots a vector of length nColleges that indicates the quota of each college
#' @param s.prefs a nColleges x nStudents matrix with each column containing student i's ranking over colleges in decreasing order of preference (most preferred first)
#' @param c.prefs a nStudents x nColleges matrix with each column containing college j's ranking over students in decreasing order of preference (most preferred first) 
#' @export
daa <- function(nStudents, nSlots, s.prefs=NULL, c.prefs=NULL){

  nColleges = length(nSlots)
  iter = 0

  if (is.null(s.prefs)){  # if no prefs given, make them randomly
    s.prefs <- replicate(n=nStudents,sample(seq(from=1,to=nColleges,by=1)))
    c.prefs <- replicate(n=nColleges,sample(seq(from=1,to=nStudents,by=1)))
  }

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