# ----------------------------------------------------------------------------
# R-code (www.r-project.org/) to unlist the matching data for OLS regression
#
# Copyright (c) 2013 Thilo Klein
#
# This library is distributed under the terms of the GNU Public License (GPL)
# for full details see the file LICENSE
#
# ----------------------------------------------------------------------------

#' @title Unlist one-sided matching data for simple regression analysis
#'
#' @description Unlist an object from \code{design.matrix} for further analysis.
#'
#' @param x matching data after applying the \code{design.matrix} function.
#' @param roommates logical: if \code{TRUE} data is assumed to come from a roomate game. This means that groups are of size two and the model matrix is prepared for individual-level analysis (peer-effects estimation). If \code{FALSE} (which is the default) data is assumed to come from a group/coalition formation game and the model matrix is prepared for group-level analysis.
#' @export
#' @return
#' \code{unlistData} returns a list with the following elements.
#' \item{SEL}{data frame comprising variables in selection equation and number of observations equal to the number of feasible groups.}
#' \item{OUT}{data frame comprising variables in outcome equation and number of observations equal to the number of equilibrium groups.}
#' @author Thilo Klein 
#' @keywords generate
#' @examples
#' ## Load one-sided matching data
#' data(mdata, package="matchingMarkets")
#' 
#' ## Unlist the data
#' mdata <- unlistData(x=mdata)
#' str(mdata$SEL); str(mdata$OUT)
unlistData <- function(x, roommates=FALSE){
  
  # --------------------------------------------------------------------
  # R-code (www.r-project.org) to unlist the matching data for OLS regression
  
  # The arguments of the function are:
  # x        : matching data after applying the oneSidedMatching() function
  # roommates: 
  
  # ## Examples:
  #
  # unlistData(x)
  # --------------------------------------------------------------------
  
  ## Outcome equation
  countX <- unlist(lapply(x$X,nrow)) # count of groups per market
  if(roommates==TRUE){ # if roommates = TRUE in model.matrix
    g.id   <- unlist(c(sapply(countX/2, function(i) rep(1:i,2) )))
  } else{
    g.id   <- unlist(c(sapply(countX, function(i) 1:i )))
  }
  m.id   <- unlist(c(sapply(1:length(countX), function(i) rep(i,countX[i]))))
  if(length(unlist(x$xi)) == length(g.id)){ # if simulation = TRUE in model.matrix
    OUT    <- with(x, data.frame( m.id, g.id, do.call(rbind,X), R=do.call(c,R),
                xi=do.call(c,xi), epsilon=do.call(c,epsilon) ))
  } else{
    OUT    <- with(x, data.frame( m.id, g.id, do.call(rbind,X), R=do.call(c,R) ))    
  }

  ## Selection equation
  countW <- unlist(lapply(x$W,nrow)) # count of groups per market
  if(roommates==TRUE){ # if roommates = TRUE in model.matrix
    g.id   <- unlist(c(sapply(1:length(countW), function(i){
      c( rep(1:(countX[i]/2), 2), rep((countX[i]/2+1):(countW[i]/2), 2) )
    } )))
  } else{
    g.id   <- unlist(c(sapply(countW, function(i) 1:i )))
  }
  m.id   <- unlist(c(sapply(1:length(countW), function(i) rep(i,countW[i]))))
  h <- unlist(lapply(x$D, length))
  x$D <- x$D[which(h>1)] # for 2-group markets only
  if(length(unlist(x$V)) == length(g.id)){ # if simulation = TRUE in model.matrix
    SEL    <- with(x, data.frame( m.id, g.id, do.call(rbind,W), D=do.call(c,D),
                V=do.call(c,V), eta=do.call(c,eta) ))
  } else{
    SEL    <- with(x, data.frame( m.id, g.id, do.call(rbind,W), D=do.call(c,D) ))    
  }
  
  return(list(OUT=OUT, SEL=SEL))
}
