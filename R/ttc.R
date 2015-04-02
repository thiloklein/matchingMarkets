# ----------------------------------------------------------------------------
# R-code (www.r-project.org/) for the Top-Trading-Cycles Algorithm
#
# Copyright (c) 2013 Thilo Klein
#
# This library is distributed under the terms of the GNU Public License (GPL)
# for full details see the file LICENSE
#
# ----------------------------------------------------------------------------

#' @title Top-Trading-Cycles Algorithm for the house allocation problem
#'
#' @description Finds the stable matching in the \href{http://en.wikipedia.org/wiki/Herbert_Scarf#8._The_Housing_Market}{house allocation problem} with existing tenants.
#' Uses the Top-Trading-Cycles Algorithm proposed in Abdulkadiroglu and Sonmez (1999).
#'
#' @param P list of individuals' preference rankings over objects.
#' @param X 2-column-matrix of objects (\code{obj}) and their owners (\code{ind}).
#' @return \code{ttc} returns a 2-column matrix of the stable matching solution for the housing market problem based on the Top-Trading-Cycles algorithm.
#' @author Thilo Klein 
#' @keywords algorithms
#' @references Abdulkadiroglu, A. and Sonmez, T. (1999). House Allocation with Existing Tenants. \emph{Journal of Economic Theory}, 88(2):233--260.
#' @examples
#' ## generate list of individuals' preference rankings over objects
#' P <- list()
#' P[[1]] <- c(2,5,1,4,3)    # individual 1
#' P[[2]] <- c(1,5,4,3,2)    # individual 2
#' P[[3]] <- c(2,1,4,3,5)    # individual 3
#' P[[4]] <- c(2,4,3,1,5)    # individual 4
#' P[[5]] <- c(4,3,1,2,5); P # individual 5
#' 
#' ## generate 2-column-matrix of objects ('obj') and their owners ('ind')
#' X <- data.frame(ind=1:5, obj=1:5); X
#' 
#' ## find assignment based on TTC algorithm
#' ttc(P=P,X=X)
#' @export
ttc <- function(P=NULL,X=NULL){
  
  ## 2-column-matrix of home objects ('obj') and their owners ('ind')
  Y <- data.frame(ind=NULL, obj=NULL)
  
  for(z in 1:length(unique(X$obj))){

    ## 1. Find cycle
    Cycle <- findCycle(P=P,X=X)
    
    ## 2. Add objects in this cycle to 'home territory'
    Y <- rbind(Y,Cycle)

    ## 3. Remove objects in this cycle from tradable objects
    X <- X[-which(X$obj %in% Cycle$obj),]
    for(i in 1:length(P)){
      P[[i]] <- P[[i]][!P[[i]] %in% Y$ind]
    }

    ## 4. Process ends if no tradable objects remain
    if(nrow(X)==0){
      Y <- rbind(Y,X)
      return(Y)
      break
    }
  }
}
findCycle <- function(P=NULL,X=NULL){
  Cycle   <- data.frame(ind=NA, obj=NA)
  thisind <- X$ind[1] # start with first individual in line
  for(j in 1:length(unique(X$ind))){
    Cycle[j,] <- c(thisind,P[[thisind]][1]) # id and top-ranked object of the individual in line
    thisind   <- X[X$obj == P[[thisind]][1],"ind"] # individual whose object is requested
    if(Cycle[j,1] == Cycle[j,2]){ # if individual points to own object
      return(Cycle[j,])
      break
    }
    if(thisind %in% Cycle$ind){ # if this individual completes a cycle
      return(Cycle)
      break
    } 
  }
}