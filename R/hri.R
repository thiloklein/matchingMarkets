# ----------------------------------------------------------------------------
# R-code (www.r-project.org/) for a constraint model to produce all stable matchings in the
# hospital/residents problem with incomplete lists
#
# Copyright (c) 2016 Thilo Klein
#
# This library is distributed under the terms of the GNU Public License (GPL)
# for full details see the file LICENSE
#
# ----------------------------------------------------------------------------

#' @title Constraint model to produce all stable matchings in the hospital/residents problem with 
#' incomplete lists
#'
#' @description Finds \emph{all} stable matchings in either the 
#' \href{http://en.wikipedia.org/wiki/Hospital_resident}{hospital/residents} problem (a.k.a. college 
#' admissions problem) or the related 
#' \href{http://en.wikipedia.org/wiki/Stable_matching}{stable marriage} problem. 
#' Dependent on the problem, the results comprise the student and college-optimal or 
#' the men and women-optimal matchings. The implementation allows for \emph{incomplete preference 
#' lists} (some agents find certain agents unacceptable) and \emph{unbalanced instances} (unequal 
#' number of agents on both sides). The function uses the Prosser (2014) constraint encoding based on 
#' either given or randomly generated preferences.
#'
#' @param nStudents integer indicating the number of students (in the college admissions problem) 
#' or men (in the stable marriage problem) in the market. Defaults to \code{ncol(s.prefs)}.
#' @param nColleges integer indicating the number of colleges (in the college admissions problem) 
#' or women (in the stable marriage problem) in the market. Defaults to \code{ncol(c.prefs)}.
#' @param nSlots vector of length \code{nColleges} indicating the number of places (i.e. 
#' quota) of each college. Defaults to \code{rep(1,nColleges)} for the marriage problem.
#' @param s.prefs matrix of dimension \code{nColleges} \code{x} \code{nStudents} with the \code{j}th 
#' column containing student \code{j}'s ranking over colleges in decreasing order of 
#' preference (i.e. most preferred first).
#' @param c.prefs matrix of dimension \code{nStudents} \code{x} \code{nColleges} with the \code{i}th 
#' column containing college \code{i}'s ranking over students in decreasing order of 
#' preference (i.e. most preferred first).
#' @param seed integer setting the state for random number generation. 
#' @param s.min integer, when specified produces incomplete preference lists with the length of each 
#' student's list randomly sampled from the range [s.min, nrow(s.prefs)].
#' @param c.min integer, when specified produces incomplete preference lists with the length of each 
#' college's list randomly sampled from the range [c.min, nrow(c.prefs)].
#' @export
#' @import stats rJava
#' @section Minimum required arguments:
#' \code{hri} requires the following combination of arguments, subject to the matching problem.
#' \describe{
#' \item{\code{nStudents, nColleges}}{Marriage problem with random preferences.}
#' \item{\code{s.prefs, c.prefs}}{Marriage problem with given preferences.}
#' \item{\code{nStudents, nSlots}}{College admissions problem with random preferences.}
#' \item{\code{s.prefs, c.prefs, nSlots}}{College admissions problem with given preferences.}
#' }
#' @return
#' \code{hri} returns a list of the following elements.
#' \item{s.prefs.smi}{student-side preference matrix for the stable marriage problem with incomplete lists (SMI).}
#' \item{c.prefs.smi}{college-side preference matrix for the stable marriage problem with incomplete lists (SMI).}
#' \item{s.prefs.hri}{student-side preference matrix for the college admissions problem (a.k.a. hospital/residents problem) with incomplete lists (HRI).}
#' \item{c.prefs.hri}{college-side preference matrix for the college admissions problem (a.k.a. hospital/residents problem) with incomplete lists (HRI).}
#' \item{matchings}{edgelist of matched students and colleges, inculding the number of the match
#' (\code{matching}) and two variables that indicate the student-optimal match (\code{sOptimal}) and 
#' college-optimal match (\code{cOptimal})}.
#' @author Thilo Klein 
#' @keywords algorithms
#' @references Gale, D. and L.S. Shapley (1962). College admissions and the stability 
#' of marriage. \emph{The American Mathematical Monthly}, 69(1):9--15.
#' 
#' Prosser, P. (2014). Stable Roommates and Constraint Programming. \emph{Lecture Notes in Computer Science, CPAIOR 2014 Edition}. 
#' Springer International Publishing, 8451: 15--28.
#' 
#' @examples
#' ######################
#' ## Marriage problem ## 
#' ######################
#' 
#' ## 3 men, 2 women, random preferences:
#' hri(nStudents=7, nColleges=6, seed=4)
#' 
#' ## 3 men, 2 women, given preferences:
#' s.prefs <- matrix(c(1,2, 1,2, 1,2), 2,3)
#' c.prefs <- matrix(c(1,2,3, 1,2,3), 3,2)
#' hri(s.prefs=s.prefs, c.prefs=c.prefs)
#' 
#' ###############################
#' ## College admission problem ##
#' ###############################
#' 
#' ## 7 students, 2 colleges with 3 slots each, random preferences:
#' hri(nStudents=7, nSlots=c(3,3), seed=21)
#' 
#' ## 7 students, 2 colleges with 3 slots each, given preferences:
#' s.prefs <- matrix(c(1,2, 1,2, 1,NA, 1,2, 1,2, 1,2, 1,2), 2,7)
#' c.prefs <- matrix(c(1,2,3,4,5,6,7, 1,2,3,4,5,NA,NA), 7,2)
#' hri(s.prefs=s.prefs, c.prefs=c.prefs, nSlots=c(3,3))
#' 
hri <- function(nStudents=ncol(s.prefs), nColleges=ncol(c.prefs), nSlots=rep(1,nColleges), 
                s.prefs=NULL, c.prefs=NULL, seed=NULL, s.min=NULL, c.min=NULL){
  
  ## ------------------------
  ## --- 1. Preliminaries ---
  
  ## set seed for random preference draws
  if(!is.null(seed)){
    set.seed(seed)
  }
  
  ## if 'nColleges' not given, obtain it from nSlots
  if(is.null(nColleges)){
    nColleges <- length(nSlots)
  }
  
  ## if no prefs given, generate random ranking:
  if(is.null(s.prefs)){  
    s.prefs <- replicate(n=nStudents,sample(seq(from=1,to=nColleges,by=1)))
  }
  if(is.null(c.prefs)){    
    c.prefs <- replicate(n=nColleges,sample(seq(from=1,to=nStudents,by=1)))
  }
  
  ## consistency checks:
  if( dim(s.prefs)[1] != dim(c.prefs)[2] | dim(s.prefs)[2] != dim(c.prefs)[1] | 
       dim(s.prefs)[2] != nStudents | dim(c.prefs)[2] != nColleges | 
       dim(c.prefs)[1] != nStudents | dim(s.prefs)[1] != nColleges ){
    stop("'s.prefs' and 'c.prefs' must be of dimensions 'nColleges x nStudents' and 'nStudents x nColleges'!")
  }
  if( length(nSlots) != nColleges | length(nSlots) != dim(c.prefs)[2] ){
    stop("Length of 'nSlots' must equal 'nColleges' and the number of columns of 'c.prefs'!")
  }
  
  ## ---------------------------------------------------------
  ## --- 2. Incomplete preferences (and consistency check) ---
  
  ## make complete preference lists incomplete
  if(!is.null(s.min)){
    s.range <- c(s.min, nrow(s.prefs))
    thre.s <- sample(x=s.range, size=ncol(s.prefs), replace=TRUE)
    s.prefs <- sapply(1:ncol(s.prefs), function(z) c(s.prefs[1:thre.s[z],z], rep(NA,s.range[2]-thre.s[z])) )
  }
  if(!is.null(c.min)){
    c.range <- c(c.min, nrow(c.prefs))
    thre.c <- sample(x=c.range, size=ncol(c.prefs), replace=TRUE)
    c.prefs <- sapply(1:ncol(c.prefs), function(z) c(c.prefs[1:thre.c[z],z], rep(NA,c.range[2]-thre.c[z])) )
  }
  
  ## make preference lists consistent (include only mutual preferences)
  c.prefs <- sapply(1:ncol(c.prefs), function(z){
    x <- na.omit(c.prefs[,z])
    y <- sapply(x, function(i) z %in% s.prefs[,i])
    return( c(x[y], rep(NA,nrow(c.prefs)-length(x[y]))) )
  })
  s.prefs <- sapply(1:ncol(s.prefs), function(z){
    x <- na.omit(s.prefs[,z])
    y <- sapply(x, function(i) z %in% c.prefs[,i])
    return( c(x[y], rep(NA,nrow(s.prefs)-length(x[y]))) )
  })
  
  ## consistency checks
  drop <- which( apply(s.prefs, 2, function(z) sum(!is.na(z))) == 0)
  if( length(drop)>0 ){
    stop(paste("Need to drop s.prefs column(s):",drop))
  }
  drop <- which( apply(c.prefs, 2, function(z) sum(!is.na(z))) == 0)
  if( length(drop)>0 ){
    stop(paste("Need to drop c.prefs column(s):",drop))
  }
  
  ## -------------------------------------------------------
  ## --- 3. Prepare preference matrices and apply solver ---
  
  ## create stable marriage instance from given hospital residents instance
  if(sum(nSlots) != length(nSlots)){ # hospital-residents problem
    
    ## set aside preferences for hospital residents instance (to be returned below)
    s.prefs.hri <- s.prefs
    c.prefs.hri <- c.prefs
    
    ## apply the transformation
    sm <- c2m(s.prefs, c.prefs, nSlots)
    s.prefs <- sm$s.prefs
    c.prefs <- sm$c.prefs
    collegeSlots <- sm$collegeSlots 
    rm(sm)
  }
  
  ## prepare and write preference matrices
  c.matrix <- sapply(1:nrow(t(c.prefs)), function(z) paste(t(c.prefs)[z,][!is.na(t(c.prefs)[z,])],collapse=" "))
  s.prefs2 <- s.prefs + ncol(s.prefs)
  s.matrix <- sapply(1:nrow(t(s.prefs2)), function(z) paste(t(s.prefs2)[z,][!is.na(t(s.prefs2)[z,])],collapse=" "))
  instance <- paste( c(length(s.matrix)+length(c.matrix),s.matrix,c.matrix), collapse="n" )
  
  ## call java jar file with choco solver
  hjw <- .jnew("smi") # create instance of smi class
  out <- .jcall(hjw, "S", "sayHello", .jarray(instance)) # invoke sayHello method
  
  ## transform the results in matrix format
  out <- as.list(strsplit(out, split="\n\n")[[1]])[[1]]
  out <- gsub(",", " ", out)
  out <- data.frame(matrix(scan(text=out, what=integer(), quiet=TRUE), ncol=3, byrow=TRUE))
  names(out) <- c("matching","student","college")
  out$college <- out$college - ncol(s.prefs)
  
  ## -----------------------------------------------------------------
  ## --- 4. Identify student-optimal and college-optimal matchings ---
  
  ## obtain sum of colleges' (and students') rankings over assigned matches in each equailibrium
  A <- split(out, out$matching)
  for(j in 1:length(A)){
    A[[j]]$srank <- sapply(1:nrow(A[[j]]), function(z) which(s.prefs[,A[[j]]$student[z]] == A[[j]]$college[z]) )
    A[[j]]$crank <- sapply(1:nrow(A[[j]]), function(z) which(c.prefs[,A[[j]]$college[z]] == A[[j]]$student[z]) )
  }
  sums <- unlist(lapply(A, function(z) sum(z$srank))) # sum students' preference rankings
  sopt.id <- which(sums == min(sums)) # s-optimal matching minimises students' ranks (i.e. maximises utility)
  sums <- unlist(lapply(A, function(z) sum(z$crank))) # # sum colleges' preference rankings
  copt.id <- which(sums == min(sums)) # c-optimal matching minimises colleges' ranks (i.e. maximises utility)
  
  ## add sOptimal and cOptimal to results
  out <- with(out, data.frame(out, sOptimal=0, cOptimal=0))
  out$sOptimal[out$matching==sopt.id] <- 1
  out$cOptimal[out$matching==copt.id] <- 1
  
  ## ----------------------------------------------------------
  ## --- 5. For hospital-residents problem: add college ids ---
  
  if(sum(nSlots) != length(nSlots)){ # hospital-residents problem
    out <- merge(x=out, y=collegeSlots, by.x="college", by.y="slots")
    out <- with(out, data.frame(matching=matching, college=colleges, slots=college, 
                                student=student, sOptimal=sOptimal, cOptimal=cOptimal))
  } else{
    out <- with(out, data.frame(matching=matching, college=college,  
                                student=student, sOptimal=sOptimal, cOptimal=cOptimal))
  }
  
  ## ----------------------------------
  ## --- 6. Sort and return results ---
  
  ## sort by matching id and college id
  out <- with(out, out[order(matching,college,student),])
  rownames(out) <- NULL
  
  ## return results
  if(sum(nSlots) != length(nSlots)){ # hospital-residents problem
    return(list(s.prefs.smi=s.prefs, c.prefs.smi=c.prefs, s.prefs.hri=s.prefs.hri,
                c.prefs.hri=c.prefs.hri, matchings=out))
  } else{ # stable-marriage problem
    return(list(s.prefs.smi=s.prefs, c.prefs.smi=c.prefs, matchings=out))
  }
}




c2m <- function(s.prefs, c.prefs, nSlots){
  
  ## expand college preferences
  c.prefs <- lapply(1:length(nSlots), function(z){
    matrix(rep(c.prefs[,z], nSlots[z]), ncol=nSlots[z])
  })
  c.prefs <- do.call(cbind,c.prefs)
  c.prefs <- rbind(c.prefs, matrix(NA, nrow=max(ncol(c.prefs)-nrow(c.prefs),0), ncol=ncol(c.prefs)))
  
  ## expand student preferences 
  fun1 <- function(i){
    x <- sapply(1:nrow(s.prefs), function(z){
      if(!is.na(s.prefs[z,i])){
        paste( rep( s.prefs[z,i], nSlots[s.prefs[z,i]] ), 1:nSlots[s.prefs[z,i]], sep=".")  
      }
    })
    if(is.list(x)){
      x <- do.call(c, x)
      #c(x, rep(NA, ncol(s.prefs)-length(x)))
      c(x, rep(NA, sum(nSlots)-length(x)))
    } else{
      c(x)
    }    
  }
  s.prefs <- sapply(1:ncol(s.prefs), function(i) fun1(i) )
  
  ## map college slots to integers
  collegeSlots_chr <- paste(rep(1:length(nSlots), nSlots), unlist(sapply(nSlots, function(z) 1:z)), sep=".")
  s.prefs <- matrix(match(c(s.prefs), collegeSlots_chr), nrow=nrow(s.prefs))
  #s.prefs <- s.prefs[-which(rowSums(is.na(s.prefs))==ncol(s.prefs)),]
  collegeSlots <- data.frame(colleges=unlist(lapply(strsplit(collegeSlots_chr,"[.]"), function(i) i[1])), 
                             slots=1:nrow(s.prefs))
  collegeSlots$colleges <- as.numeric(as.character(collegeSlots$colleges))
  
  return(list(s.prefs=s.prefs, c.prefs=c.prefs, collegeSlots=collegeSlots))
  
}

