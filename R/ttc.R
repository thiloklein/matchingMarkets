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
#' @description Finds the stable matching in the 
#' \href{https://en.wikipedia.org/wiki/Top_trading_cycle}{house allocation problem} 
#' with existing tenants. Uses the Top-Trading-Cycles Algorithm proposed in Abdulkadiroglu and Sonmez (1999).
#'
#' @param P matrix of individuals' preference rankings (Rank Order Lists) over objects.
#' @param X 2-column-matrix of objects (\code{obj}) and their owners (\code{ind}).
#' @return \code{ttc} returns a 2-column matrix of the matching solution for the housing market problem based on the Top-Trading-Cycles algorithm.
#' @author Thilo Klein 
#' @keywords algorithms
#' @references Abdulkadiroglu, A. and Sonmez, T. (1999). House Allocation with Existing Tenants. \emph{Journal of Economic Theory}, 88(2):233--260.
#' 
#' Shapley, L. and H. Scarf (1974). On Cores and Indivisibility. \emph{Journal of Mathematical Economics}, 1(1):23--37.
#' @examples
#' ## generate matrix of individuals' preference rankings over objects, a.k.a. Rank Order Lists (ROL)
#' P <- matrix(c(2, 5, 1, 4, 3,  # ROL of individual 1
#'               1, 5, 4, 3, 2,  # ind 2
#'               2, 1, 4, 3, 5,  # ind 3
#'               2, 4, 3,NA,NA,  # ind 4
#'               4, 3, 1, 2,NA), # ind 5
#'             byrow=FALSE, nrow=5); P 
#' 
#' ## generate 2-column-matrix of objects ('obj') and their owners ('ind')
#' X <- data.frame(ind=1:5, obj=5:1); X
#' 
#' ## find assignment based on TTC algorithm
#' ttc(P=P,X=X)
#' @export
ttc <- function(P=NULL,X=NULL){
  
  ## convert matrix to list and drop NAs
  P <- as.list(data.frame(P))
  P <- lapply(P, function(i) i[!is.na(i)])
  
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
      P[[i]] <- P[[i]][!P[[i]] %in% Y$obj]
    }

    ## 4. Process ends if no tradable objects remain
    if(nrow(X)==0){
      Y <- rbind(Y,X)
      return(Y[order(Y$ind),])
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
    
    if(thisind %in% Cycle$ind){ # if this individual (thisind) completes a cycle
      h <- which(Cycle$ind == thisind)
      return(Cycle[h:length(Cycle$ind),])
      break
    } 
  }
}