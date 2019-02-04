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

#' @title All stable matchings in the hospital/residents problem with incomplete lists
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
#' @param s.range range of two intergers \code{s.range = c(s.min, s.max)}, where \code{s.min < s.max}. 
#' Produces incomplete preference lists with the length of each student's list randomly sampled from 
#' the range \code{[s.min, s.max]}. Note: interval is only correct if either c.range or s.range is used.
#' @param c.range range of two intergers \code{c.range = c(c.min, c.max)}, where \code{c.min < c.max}. 
#' Produces incomplete preference lists with the length of each college's list randomly sampled from 
#' the range \code{[c.min, c.max]}. Note: interval is only correct if either c.range or s.range is used.
#' @param randomization determines at which level random lottery numbers for student priorities are drawn. The default is \code{randomization = "multiple"}, where a student's priority is determined by a separate lottery at each college (i.e. local tie-breaking). For the second variant, \code{randomization = "single"}, a single lottery number determines a student's priority at all colleges (i.e. global tie breaking). 
#' @param check_consistency Performs consicentcy checks (Checks if there are columns in the preference matrices that only contains zeros and drops them and checks the matrixes for consistencies if they are given by characters). Defaults to \code{TRUE} but changing it to \code{FALSE} might reduce the running-time for large problems.
#' @param ... .
#' 
#' @export
#' 
#' @import stats rJava
#' @importFrom grDevices rgb
#' @importFrom graphics legend
#' 
#' @section Minimum required arguments:
#' \code{hri} requires the following combination of arguments, subject to the matching problem.
#' \describe{
#' \item{\code{nStudents, nColleges}}{Marriage problem with random preferences.}
#' \item{\code{s.prefs, c.prefs}}{Marriage problem with given preferences.}
#' \item{\code{nStudents, nSlots}}{College admissions problem with random preferences.}
#' \item{\code{s.prefs, c.prefs, nSlots}}{College admissions problem with given preferences.}
#' }
#' 
#' @return
#' \code{hri} returns a list of the following elements.
#' \item{s.prefs.smi}{student-side preference matrix for the stable marriage problem with incomplete lists (SMI).}
#' \item{c.prefs.smi}{college-side preference matrix for the stable marriage problem with incomplete lists (SMI).}
#' \item{s.prefs.hri}{student-side preference matrix for the college admissions problem (a.k.a. hospital/residents problem) with incomplete lists (HRI).}
#' \item{c.prefs.hri}{college-side preference matrix for the college admissions problem (a.k.a. hospital/residents problem) with incomplete lists (HRI).}
#' \item{matchings}{edgelist of matched students and colleges, inculding the number of the match
#' (\code{matching}) and two variables that indicate the student-optimal match (\code{sOptimal}) and 
#' college-optimal match (\code{cOptimal})}.
#' 
#' @author Thilo Klein
#' 
#' @keywords algorithms
#' 
#' @references Gale, D. and L.S. Shapley (1962). College admissions and the stability 
#' of marriage. \emph{The American Mathematical Monthly}, 69(1):9--15.
#' 
#' Morizumi, Y., T. Hayashi and Y. Ishida (2011). A network visualization of stable matching in the stable 
#' marriage problem. \emph{Artificial Life Robotics}, 16:40--43.
#' 
#' Prosser, P. (2014). Stable Roommates and Constraint Programming. \emph{Lecture Notes in Computer Science, CPAIOR 2014 Edition}. 
#' Springer International Publishing, 8451: 15--28.
#' 
#' @examples
#' \dontrun{
#' ## -----------------------
#' ## --- Marriage problem 
#' 
#' ## 7 men, 6 women, random preferences:
#'  hri(nStudents=7, nColleges=6, seed=4)
#' 
#' ## 3 men, 2 women, given preferences:
#'  s.prefs <- matrix(c(1,2, 1,2, 1,2), 2,3)
#'  c.prefs <- matrix(c(1,2,3, 1,2,3), 3,2)
#'  hri(s.prefs=s.prefs, c.prefs=c.prefs)
#' 
#' ## 3 men, 2 women, given preferences:
#'  s.prefs <- matrix(c("x","y", "x","y", "x","y"), 2,3)
#'  colnames(s.prefs) <- c("A","B","C")
#'  c.prefs <- matrix(c("A","B","C", "A","B","C"), 3,2)
#'  colnames(c.prefs) <- c("x","y")
#'  hri(s.prefs=s.prefs, c.prefs=c.prefs)
#' 
#' ## --------------------------------
#' ## --- College admission problem 
#' 
#' ## 7 students, 2 colleges with 3 slots each, random preferences:
#'  hri(nStudents=7, nSlots=c(3,3), seed=21)
#' 
#' ## 7 students, 2 colleges with 3 slots each, given preferences:
#'  s.prefs <- matrix(c(1,2, 1,2, 1,NA, 1,2, 1,2, 1,2, 1,2), 2,7)
#'  c.prefs <- matrix(c(1,2,3,4,5,6,7, 1,2,3,4,5,NA,NA), 7,2)
#'  hri(s.prefs=s.prefs, c.prefs=c.prefs, nSlots=c(3,3))
#'  
#' ## 7 students, 2 colleges with 3 slots each, given preferences:
#'  s.prefs <- matrix(c("x","y", "x","y", "x",NA, "x","y", 
#'                      "x","y", "x","y", "x","y"), 2,7)
#'  colnames(s.prefs) <- c("A","B","C","D","E","F","G")
#'  c.prefs <- matrix(c("A","B","C","D","E","F","G", 
#'                      "A","B","C","D","E",NA,NA), 7,2)
#'  colnames(c.prefs) <- c("x","y")
#'  hri(s.prefs=s.prefs, c.prefs=c.prefs, nSlots=c(3,3))
#'  
#' ## 7 students, 3 colleges with 3 slots each, incomplete preferences:
#'  hri(nStudents=7, nSlots=c(3,3,3), seed=21, s.range=c(1,3))
#'  
#'  s.prefs <- matrix(c('S1', 'S2', NA,
#'                      'S3', 'S1', NA,
#'                      'S1', NA, NA,
#'                       NA, NA,NA,
#'                      'S2', 'S1', 'S5'),
#'                    nrow = 3, ncol = 5)
#'  # Note that we explicitly allow for the existence of entries refering to colleges
#'  # that do not exist. A warning is generated and the entry is ignored.
#'  colnames(s.prefs) <- c('A', 'B', 'C', 'D', 'E')
#'  c.prefs <- matrix(c('B', 'C','D', 'A',
#'                      'C', 'D', NA, NA, 
#'                      'D', 'B', 'A', 'E'),
#'                    nrow = 4, ncol = 3)
#'  colnames(c.prefs) <- c('S1', 'S2', 'S3')
#'  hri(s.prefs=s.prefs, c.prefs=c.prefs, nSlots=c(3,3,3), check_consistency = TRUE)
#' ## --------------------
#' ## --- Summary plots
#' 
#' ## 200 students, 200 colleges with 1 slot each
#'  res <- hri(nStudents=200, nColleges=200, seed=12)
#'  plot(res)
#'  plot(res, energy=TRUE)
#' }
hri <- function(nStudents=ncol(s.prefs), nColleges=ncol(c.prefs), nSlots=rep(1,nColleges), 
                s.prefs=NULL, c.prefs=NULL, s.range=NULL, c.range=NULL, randomization=NULL, seed=NULL, check_consistency = TRUE, ...) UseMethod("hri")

#' @export
hri.default <- function(nStudents=ncol(s.prefs), nColleges=ncol(c.prefs), nSlots=rep(1,nColleges), 
                        s.prefs=NULL, c.prefs=NULL, s.range=NULL, c.range=NULL, randomization='multiple', seed=NULL, check_consistency = TRUE, ...){
  
  
  ###############################################################  
  #print('Section 1a')
  #print(Sys.time())
  ## ------------------------
  ## --- 1-a. Preliminaries ---
  
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
    if(randomization == "single"){ 
      
      c.prefs <- matrix(sample(seq(from=1, to=nStudents, by=1)), nrow=nStudents, ncol=nColleges) 
      
    } else{ # if(randomization == "multiple")
      
      c.prefs <- replicate(n=nColleges, sample(seq(from=1, to=nStudents, by=1)))
    }
  }
  
  
  if( length(nSlots) != nColleges | length(nSlots) != dim(c.prefs)[2] ){
    stop("Length of 'nSlots' must equal 'nColleges' and the number of columns of 'c.prefs'!")
  }
  # Print a warning if colnames are given but the pref matrices are numeric
  if( !is.null(colnames(s.prefs)) || !is.null(colnames(c.prefs))){
    if(is.numeric(c.prefs) || is.numeric(s.prefs)){
      warning('Colnames are given but preferences matrizes are numeric. Make sure that they are consistent!')
    }
  }
  
  # Are preferences given as characters?
  prefs_char <- FALSE
  if(is.character(s.prefs) || is.character(c.prefs)){
    if(!(is.character(s.prefs) && is.character(c.prefs))){
      stop('Both prefs must be in characters')
    }
    prefs_char <- TRUE
  } 
  
  ## add colnames index to s.prefs and c.prefs if colnames are NULL 
  if(is.null(colnames(s.prefs))){
    colnames(s.prefs) <- 1:ncol(s.prefs)
  }
  if(is.null(colnames(c.prefs))){
    colnames(c.prefs) <- 1:ncol(c.prefs)
  }
  
  
  # Generated s.names and c.names since they are necessary for the consistency check
  c.names <- colnames(c.prefs)
  s.names <- colnames(s.prefs)
  # c.prefs_named will be created later but initialize it here. In this way it will not 
  # be generated in the global environment.
  c.prefs_named <- NULL
  s.prefs_named <- NULL 
  
  # Perform consistency checks
  if(check_consistency){
    
    # Check if student names are unique:
    if(length(unique(colnames(s.prefs))) != ncol(s.prefs)) {stop('Student names not unique')}
    
    # Check if colleges are unique:
    if(length(unique(colnames(c.prefs))) != ncol(c.prefs)) {stop('College/Course names not unique')}
    
    # Check if a student applied for a course/college that is not in c.prefs
    applied_colleges <- unique(as.character(s.prefs))
    applied_colleges <- applied_colleges[!is.na(applied_colleges)]
    if(length(setdiff(applied_colleges,colnames(c.prefs))) != 0){
      missing_college <- setdiff(applied_colleges,colnames(c.prefs))
      missing_college <- paste('Someone applied to a college (named ', missing_college, ') that has no ranking. This preference entries will be deleted! \n', sep = '')
      warning(missing_college)
      
      # Delete the preference entries that refer to non existing columns in the c.pref matrix:
      
      s.prefs <- sapply(s.names, function(z){
        x <- c(na.omit(s.prefs[,z]))
        if(length(x) == 0) return(rep(NA, length(s.prefs[,z])))
        y <- sapply(x, function(i) i %in% colnames(c.prefs))
        return( c(x[y], rep(NA,nrow(s.prefs)-length(x[y]))) )
      })
      
    }
    
    # Check if a college ranked someone who is not in s.prefs
    ranked_stud <- unique(c(as.character(c.prefs)))
    ranked_stud <- ranked_stud[!is.na(ranked_stud)]
    if(length(setdiff(ranked_stud,colnames(s.prefs))) != 0){
      missing_stud <- setdiff(ranked_stud,colnames(s.prefs))
      missing_stud <- paste('A course/college ranked a student (named ', missing_stud, ') that has no ranking.  This preference entries will be deleted! \n', sep = '')
      warning(missing_stud)
    }
    # Delete the preference entries that refer to non existing columns in the s.pref matrix:
    c.prefs <- sapply(c.names, function(z){
      x <- c(na.omit(c.prefs[,z]))
      if(length(x) == 0) return(rep(NA, length(c.prefs[,z])))
      y <- sapply(x, function(i) i %in% colnames(s.prefs))
      return( c(x[y], rep(NA,nrow(c.prefs)-length(x[y]))) )
    })
    #print('Input passed consistency test!')
  }
  
  
  
  ###############################################################  
  #print('Section 2')
  #print(Sys.time())
  ## ---------------------------------------------------------
  ## --- 2-a. Incomplete preferences (and consistency check) ---
  
  
  ## make complete preference lists incomplete
  if(!is.null(s.range)){
    #s.range <- c(s.min, nrow(s.prefs))
    thre.s <- sample(x=s.range[1]:s.range[2], size=ncol(s.prefs), replace=TRUE)
    s.prefs <- sapply(1:ncol(s.prefs), function(z) c(s.prefs[1:thre.s[z],z], rep(NA,nrow(s.prefs)-thre.s[z])) )
    colnames(s.prefs) <- s.names
    rownames(s.prefs) <- NULL
  }
  if(!is.null(c.range)){
    #c.range <- c(c.min, nrow(c.prefs))
    thre.c <- sample(x=c.range[1]:c.range[2], size=ncol(c.prefs), replace=TRUE)
    c.prefs <- sapply(1:ncol(c.prefs), function(z) c(c.prefs[1:thre.c[z],z], rep(NA,nrow(c.prefs)-thre.c[z])) )
    colnames(c.prefs) <- c.names
    rownames(c.prefs) <- NULL
  }
  
  
  
  ## --- 2-b. Make prefs consistent and map names to IDs---
  
  # Drop make preferences consistent and drop NA columns

  # Check prefs for mutual consistency and delete all unneccessary entries
  c.prefs <- sapply(c.names, function(z){
    x <- c(na.omit(c.prefs[,z]))
    if(length(x) == 0) return(rep(NA, length(c.prefs[,z])))
    y <- sapply(x, function(i) z %in% s.prefs[,i])
    return( c(x[y], rep(NA,nrow(c.prefs)-length(x[y]))) )
  })
  s.prefs <- sapply(s.names, function(z){
    x <- c(na.omit(s.prefs[,z]))
    if(length(x) == 0) return(rep(NA, length(s.prefs[,z])))
    y <- sapply(x, function(i) z %in% c.prefs[,i])
    return( c(x[y], rep(NA,nrow(s.prefs)-length(x[y]))) )
  })
  

  # drop columns that contain only missings and update nSlots!
  drop <- which( apply(s.prefs, 2, function(z) all(is.na(z))))
  if( length(drop)>0 ){
    s.prefs <- matrix(s.prefs[,-drop], nrow=nrow(s.prefs))
    colnames(s.prefs) <- s.names[-drop]
    print(paste("Dropped s.prefs column(s):", paste(s.names[drop], collapse=", ")))
  }
  drop <- which( apply(c.prefs, 2, function(z) all(is.na(z))))
  if( length(drop)>0 ){
    c.prefs <- matrix(c.prefs[,-drop], nrow=nrow(c.prefs))
    colnames(c.prefs) <- c.names[-drop]
    print(paste("Dropped c.prefs column(s):", paste(c.names[drop], collapse=", ")))
    
    # Update nSlots:
    nSlots <- nSlots[-drop]
  }
  
  # Replace character names with numeric IDs
  # Probably this should be done regardless of the type of prefs (numeric or character)
  # since even in the numeric case the column numbers might not correspond to the ids since some
  # columns might have been deleted (only nas after mutually consistent)

  # Replace college/student names with ids
  
  # Save prefs
  s.prefs_named <-  s.prefs;
  c.prefs_named <-  c.prefs;
  
  # Replace names with IDs
  s.names <- 1:ncol(s.prefs)
  names(s.names) <- colnames(s.prefs)
  c.names <- 1:ncol(c.prefs)
  names(c.names) <- colnames(c.prefs)
  s.names <-  s.names
  c.names <-  c.names
  
  s.prefs <- apply(s.prefs, 2, function(pref){
    c.names[as.character(pref)]
  })
  colnames(s.prefs) <- s.names[colnames(s.prefs)]
  s.prefs <-  s.prefs
  
  c.prefs <- apply(c.prefs, 2, function(pref){
    s.names[as.character(pref)]
  })
  colnames(c.prefs) <- c.names[colnames(c.prefs)]
  c.prefs <-  c.prefs
  
  ###############################################################  
  #print('Section 3')
  #print(Sys.time())
  ## -------------------------------------------------------
  ## --- 3. Prepare preference matrices and apply solver ---
  
  ## create stable marriage instance from given hospital residents instance
  if(sum(nSlots) != length(nSlots)){ # hospital-residents problem
    
    ## set aside preferences for hospital residents instance (to be returned below)
    if(!prefs_char){
      s.prefs.hri <- s.prefs
      c.prefs.hri <- c.prefs
    }
    
    ## apply the transformation
    sm <- c2m(s.prefs, c.prefs, nSlots)
    s.prefs <- sm$s.prefs
    c.prefs <- sm$c.prefs
    collegeSlots <- sm$collegeSlots 
    rm(sm)
  }
  
  ## prepare and write preference matrices
  c.matrix <- sapply(1:nrow(t(c.prefs)), function(z) paste(t(c.prefs)[z,][!is.na(t(c.prefs)[z,])],collapse=" "))
  s.prefs1 <- s.prefs + ncol(s.prefs)
  s.matrix <- sapply(1:nrow(t(s.prefs1)), function(z) paste(t(s.prefs1)[z,][!is.na(t(s.prefs1)[z,])],collapse=" "))
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
  
  
  ###############################################################
  #print('Section 4')
  #print(Sys.time())
  ## -----------------------------------------------------------------
  ## --- 4. Identify student-optimal and college-optimal matchings ---
  
  ## obtain sum of colleges' (and students') rankings over assigned matches in each equailibrium
  A <- split(out, out$matching)
  for(j in 1:length(A)){
    A[[j]]$sRank <- sapply(1:nrow(A[[j]]), function(z) which(s.prefs[,A[[j]]$student[z]] == A[[j]]$college[z]) )
    A[[j]]$cRank <- sapply(1:nrow(A[[j]]), function(z) which(c.prefs[,A[[j]]$college[z]] == A[[j]]$student[z]) )
  }
  sums <- unlist(lapply(A, function(z) sum(z$sRank))) # sum students' preference rankings
  sopt.id <- which(sums == min(sums)) # s-optimal matching minimises students' ranks (i.e. maximises utility)
  sums <- unlist(lapply(A, function(z) sum(z$cRank))) # # sum colleges' preference rankings
  copt.id <- which(sums == min(sums)) # c-optimal matching minimises colleges' ranks (i.e. maximises utility)
  
  ## add sOptimal and cOptimal to results
  out <- with(out, data.frame(out, sOptimal=0, cOptimal=0))
  out$sOptimal[out$matching==sopt.id] <- 1
  out$cOptimal[out$matching==copt.id] <- 1
  
  ## add student/college ranking (sRank, cRank)
  out <- cbind(out, do.call("rbind",A)[,c("sRank","cRank")])
  
  ###############################################################
  #print('Section 5')
  #print(Sys.time())
  ## ----------------------------------------------------------
  ## --- 5. For hospital-residents problem: add college ids ---
  
  if(sum(nSlots) != length(nSlots)){ # hospital-residents problem
    
    out <- merge(x=out, y=collegeSlots, by.x="college", by.y="slots")
    
    # Reverse first ids
    out$student <- colnames(s.prefs_named)[out$student]
    out$colleges <- colnames(c.prefs_named)[out$colleges]
    
    out <- with(out, data.frame(matching=matching, college=colleges, slots=college,
                                student=student, sOptimal=sOptimal, cOptimal=cOptimal,
                                sRank=sRank, cRank=cRank, stringsAsFactors=FALSE))
    
    
  } else {   # stable-marriage problem
    
    # Reverse first ids
    out$student <- colnames(s.prefs_named)[out$student]
    out$college <- colnames(c.prefs_named)[out$college]
    
    out <- with(out, data.frame(matching=matching, college=college, slots=college,
                                student=student, sOptimal=sOptimal, cOptimal=cOptimal,
                                sRank=sRank, cRank=cRank, stringsAsFactors=FALSE))
    
  }
  
  
  ###############################################################
  #print('Section 6')
  #print(Sys.time())
  ## ----------------------------------
  ## --- 6. Sort and return results ---
  
  ## sort by matching id and college id
  out <- with(out, out[order(matching,college,student),])
  rownames(out) <- NULL
  
  ## return results
  if(sum(nSlots) != length(nSlots)){ # hospital-residents problem
    
    # Prefs as character
    
    res <- list(
      
      s.prefs.smi = s.prefs[rowSums(is.na(s.prefs))<ncol(s.prefs),], 
      c.prefs.smi = c.prefs[rowSums(is.na(c.prefs))<ncol(c.prefs),], 
      s.prefs.hri = s.prefs_named[rowSums(is.na(s.prefs_named))<ncol(s.prefs_named),],
      c.prefs.hri = c.prefs_named[rowSums(is.na(c.prefs_named))<ncol(c.prefs_named),],
      
      matchings = out
    )
    
    
  } else { # stable-marriage problem
    
    # Prefs as characters
    res <- list( 
      s.prefs.smi = s.prefs_named[rowSums(is.na(s.prefs))<ncol(s.prefs),],
      c.prefs.smi =c.prefs_named[rowSums(is.na(c.prefs))<ncol(c.prefs),],
      matchings = out
    )
  }
  #print('Section return')
  #print(Sys.time())
  #res$call <- match.call()
  class(res) <- "hri"
  return(res)
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
      #c(x)
      c(x, rep(NA, sum(nSlots)-length(x)))
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



#' @export
plot.hri <- function(x, energy=FALSE, ...){
  
  ## obtain edgelist of stable matchings
  x <- x$matchings
  
  ## add satisfaction and energy measure
  n <- nrow(x[x$matching==1,])
  x <- with(x, data.frame(x, sSatisf = n + 1 - sRank, cSatisf = n +1 - cRank))
  x$energy <- with(x, sSatisf*cSatisf)
  
  ## aggregate by matching
  x <- x[,-which(names(x) %in% c("college","slots","student"))]  
  x <- aggregate(x, by=list(x$matching), sum)
  x$csSatisf <- with(x, sSatisf - cSatisf)
  
  ## sort by student net satisfaction
  x <- with(x, x[order(csSatisf),])
  
  ## set graphic parameters
  par(mar = c(5.1, 4.1, 1.8, 0.5))
  lightgray <- rgb(84,84,84, 50, maxColorValue=255) 
  
  add_legend <- function(...) {
    opar <- par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0), 
                mar=c(0, 0, 0, 0), new=TRUE)
    on.exit(par(opar))
    plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
    legend(...)
  }
  
  ## plots
  if(energy==TRUE){
    
    with(x, plot(energy ~ csSatisf, col="gray80", type="l", 
                 xlab=expression(paste("net satisfaction: ",P[s]-P[c])), ylab=expression(paste("energy: ",E))))
    with(x, points(energy ~ csSatisf, col=lightgray, pch=16))
    with(x, points(energy ~ csSatisf, col=lightgray))
    with(x, points(energy ~ csSatisf, data=x[sOptimal>0,], pch=16, col="black"))
    with(x, points(energy ~ csSatisf, data=x[sOptimal>0,], col="black"))
    with(x, points(energy ~ csSatisf, data=x[cOptimal>0,], pch=16, col="white"))
    with(x, points(energy ~ csSatisf, data=x[cOptimal>0,], col="black"))
    
  } else{
    
    with(x, plot(cSatisf ~ sSatisf, col="gray80", type="l", 
                 xlab=expression(paste("student satisfaction: ",P[s])), 
                 ylab=expression(paste("college satisfaction: ",P[c]))))
    with(x, points(cSatisf ~ sSatisf, col=lightgray, pch=16))
    with(x, points(cSatisf ~ sSatisf, col=lightgray))
    with(x, points(cSatisf ~ sSatisf, data=x[sOptimal>0,], pch=16, col="black"))
    with(x, points(cSatisf ~ sSatisf, data=x[sOptimal>0,], col="black"))
    with(x, points(cSatisf ~ sSatisf, data=x[cOptimal>0,], pch=16, col="white"))
    with(x, points(cSatisf ~ sSatisf, data=x[cOptimal>0,], col="black"))
  }
  
  ## legend
  add_legend("topright", legend=c("c-optimal", "s-optimal", "stable"), 
             pch=c(16,16,16), col=c("white","black",lightgray), 
             text.col="white", horiz=TRUE, bty='n')
  add_legend("topright", legend=c("c-optimal", "s-optimal", "stable"), 
             pch=c(1,1,1), col=c("black","black",lightgray), 
             horiz=TRUE, bty='n')
  
  ## reset R default
  par(mar=c(5.1, 4.1, 4.1, 2.1)) 
  
}


