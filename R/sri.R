# ----------------------------------------------------------------------------
# R-code (www.r-project.org/) for a constraint model to produce all stable matchings in the
# stable roommates problem with incomplete lists
#
# Copyright (c) 2016 Thilo Klein
#
# This library is distributed under the terms of the GNU Public License (GPL)
# for full details see the file LICENSE
#
# ----------------------------------------------------------------------------

#' @title All stable matchings in the stable roommates problem with incomplete lists 
#'
#' @description Finds all stable matchings (if one exists) in the 
#' \href{http://en.wikipedia.org/wiki/Stable_roommates_problem}{stable roommates problem} with incomplete 
#' lists using the Prosser (2014) constraint encoding based on either given or randomly generated preferences.
#'
#' @param nAgents integer that gives the number of players in the market.
#' @param prefs valuation matrix of dimension \code{nAgents x nAgents} that gives column-players' 
#' ranking over other players in decreasing order of preference (i.e. most preferred first).
#' @param seed integer setting the state for random number generation. 
#' @param p.range range of two intergers \code{p.range = c(p.min, p.max)}, where \code{p.min < p.max}. 
#' Produces incomplete preference lists with the length of each player's list randomly sampled from 
#' the range \code{[p.min, p.max]}. Note: interval is always too small because non-permissible matches 
#' are automatically deleted.
#' 
#' @export
#' 
#' @return
#' \code{sri} returns a list with the following items.
#' \item{prefs}{agents' preference list.}
#' \item{matching}{edgelist of matched pairs, inculding the number of the match (\code{matching}).}
#' 
#' @author Thilo Klein 
#' 
#' @keywords algorithms
#' 
#' @import stats rJava
#' 
#' @references Gusfield, D.M. and R.W. Irving (1989). The Stable Marriage Problem: 
#' Structure and Algorithms, MIT Press.
#' 
#' Prosser, P. (2014). Stable Roommates and Constraint Programming. \emph{Lecture Notes in Computer Science, CPAIOR 2014 Edition}. 
#' Springer International Publishing, 8451: 15--28.
#' 
#' Irving, R.W. and S. Scott (2007). The stable fixtures problem: A many-to-many extension of stable roommates. 
#' \emph{Discrete Applied Mathematics}, 155: 2118--2129.
#' 
#' @examples
#' ## Roommate problem with 10 players, given preferences:
#'  prefs <- matrix(rep(1:10, 10), 10, 10)
#'  sri(prefs=prefs)
#' 
#' ## Roommate problem with 10 players, random preferences:
#'  sri(nAgents=10, seed=1)
#' 
#' ## Roommate problem with no equilibrium matching:
#'  sri(nAgents=10, seed=2)
#' 
#' ## Roommate problem with 3 equilibria:
#'  sri(nAgents=10, seed=3)
#' 
sri <- function(prefs=NULL, nAgents=NULL, seed=NULL, p.range=NULL){

  ## ------------------------
  ## --- 1. Preliminaries ---
  
  ## Simulate preferences (or use given preferences)
  if(is.null(prefs)==TRUE){
    if(!is.null(seed)){
      set.seed(seed)
    }
    prefs <- replicate(n=nAgents, sample(seq(from=1, to=nAgents, by=1)))
  } else{
    nAgents <- ncol(prefs)
  }

  ## Consistency checks
  if(dim(prefs)[1] != dim(prefs)[2]){stop("preference matrix must be symmetric!")}
  
  ## ---------------------------------------------------------
  ## --- 2. Incomplete preferences (and consistency check) ---
  
  ## make complete preference lists incomplete
  if(!is.null(p.range)){
    #p.range <- c(p.min, nrow(prefs))
    thres <- sample(x=p.range[1]:p.range[2], size=ncol(prefs), replace=TRUE)
    prefs <- sapply(1:ncol(prefs), function(z) c(prefs[1:thres[z],z], rep(NA,nrow(prefs)-thres[z])) )
  }
  
  ## make preference lists consistent (include only mutual preferences)
  prefs <- sapply(1:ncol(prefs), function(z){
    x <- na.omit(prefs[,z])
    y <- sapply(x, function(i) z %in% prefs[,i])
    return( c(x[y], rep(NA,nrow(prefs)-length(x[y]))) )
  })
  
  ## -------------------------------------------------------
  ## --- 3. Prepare preference matrices and apply solver ---
  
  ## prepare and write preference matrices
  prefs <- sapply(1:nAgents, function(z){ # drop "self-preferences"
    x <- prefs[,z][prefs[,z] != z]
    c(x, rep(NA, nAgents-length(x)))  
  })
  prefs <- prefs[rowSums(is.na(prefs))<ncol(prefs),]
  v.matrix <- sapply(1:nrow(t(prefs)), function(z) paste(t(prefs)[z,][!is.na(t(prefs)[z,])],collapse=" "))
  instance <- paste( c(length(v.matrix),v.matrix), collapse="n" )
  
  ## call java jar file with choco solver
  hjw <- .jnew("smi") # create instance of smi class
  out <- .jcall(hjw, "S", "sayHello", .jarray(instance)) # invoke sayHello method
  
  ## -------------------------
  ## --- 4. Return results ---
  
  if(substr(out,1,12) != "\nsolutions:0" ){ ## solutions exist
  
    ## transform the results in matrix format
    out <- as.list(strsplit(out, split="\n\n")[[1]])[[1]]
    out <- gsub(",", " ", out)
    out <- data.frame(matrix(scan(text=out, what=integer(), quiet=TRUE), ncol=3, byrow=TRUE))
    names(out) <- c("matching","playerA","playerB")
    
    list(prefs=prefs, matchings=out)
    
  } else{ ## solutions don't exist
    
    list(prefs=prefs, matchings="No stable matching exists for this instance.")
  }
}

