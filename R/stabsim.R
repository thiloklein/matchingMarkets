# ----------------------------------------------------------------------------
# R-code (www.r-project.org/) for simulating individual-level data for
# all players in the market.
#
# Copyright (c) 2013 Thilo Klein
#
# This library is distributed under the terms of the GNU Public License (GPL)
# for full details see the file LICENSE
#
# ----------------------------------------------------------------------------

#' @title Simulate individual-level data for one-sided matching markets
#'
#' @description Simulate individual-level data for one-sided matching markets.
#'
#' @param m integer indicating the number of markets to be simulated.
#' @param ind integer (or vector) indicating the number of individuals per group.
#' @param seed integer setting the state for random number generation. Defaults to \code{set.seed(123)}.
#' @param singles integer giving the number of one-group markets.
#' @param gpm integer giving the number of groups per market.
#' 
#' @export
#' 
#' @import stats
#' 
#' @return
#' \code{stabsim} returns a data frame with the randomly generated variables 
#' mimicking those in dataset \code{\link{baac00}}.
#' \item{m.id}{categorical: market identifier.}
#' \item{g.id}{categorical: group identifier.}
#' \item{wst}{binary: indicator taking the value 1 if last year was worse than the year before; 0 otherwise.}
#' \item{R}{NA: group outcome is not simulated. It can be obtained using the  \code{simulation} argument 
#' in function \code{stabit}.}
#' 
#' @author Thilo Klein 
#' 
#' @keywords generate
#' 
#' @examples
#' ## Coalitions [gpm := 2 !]
#' ## Simulate one-sided matching data for 4 markets (m=4) with 2 groups
#' ## per market (gpm=2) and 2 to 4 individuals per group (ind=2:4)
#'  idata <- stabsim(m=4, ind=2:4, seed=124, singles=2, gpm=2)  
#' 
#' ## Rommmates [ind := 2 !]
#' ## Simulate one-sided matching data for 3 markets (m=3) with 3 groups
#' ## per market (gpm=3) and 2 individuals per group (ind=2)
#'  idata <- stabsim(m=3, ind=2, seed=124, gpm=3)
stabsim <- function(m, ind, seed=123, singles=NULL, gpm=2){
  
  # --------------------------------------------------------------------
  # R-code (www.r-project.org) for simulating purely random (!) data for
  # all players in the market.
  
  # The arguments of the function are:
  # m       : integer indicating the number of markets to be simulated
  # ind     : integer indicating the number of individuals per group
  #           or vector
  # seed    : seed, defaults to set.seed(123)
  # singles : number of 1-group markets
  # gpm     : number of groups per market
  
  # ## Examples:
  #
  # stabsim(m=30, ind=2:4, seed=124, singles=5, gpm=2)  # coalitions
  # stabsim(m=3, ind=2, seed=124, gpm=3)  # rommmates
  # --------------------------------------------------------------------
  
  set.seed(seed)
  if(length(ind)==1){
    g.s <- rep(ind, gpm*m)
  } else{
    g.s <- sample(ind, gpm*m, replace=TRUE)
  }
  g.id  <- c(unlist( sapply(1:length(g.s), function(x) rep(x,g.s[x])) ))
  m.id  <- c(unlist( sapply(1:(length(g.s)/gpm), function(x) rep(x, sum(g.s[(1:gpm)+(x*gpm-gpm)]) )) )) 
  Nrows <- length(m.id)
  Ncols <- max(table(m.id))
  i.id  <- 1:Nrows

  #pi  <- runif(n=Nrows,min=0.5,max=1)
  
  #wst <- unlist(c(by(m.id, m.id, function(i){
  #  l <- length(i) # market size
  #  s <- 0.3 # share
  #  r <- round(s*l,0) # size of smaller share
  #  sample(c(rep(0,r),rep(1,l-r)),l,replace=FALSE)
  #})))
  wst <- sample(0:1,Nrows,replace=TRUE)
  
  #occ1 <- runif(Nrows); occ2 <- runif(Nrows,max=1-occ1); occ3 <- runif(Nrows,max=1-occ1-occ2)
  #sat <- rnorm(Nrows); mot <- rnorm(Nrows)
  R <- rep(NA,Nrows)
  
  #xi.i  <- rnorm(n=Nrows, sd=sqrt(max(ind)))
  #eta.i <- rnorm(n=Nrows, sd=sqrt(max(ind)))
  
  x <- data.frame(m.id, g.id, wst, R)
  
  if(is.null(singles)){
    return(x)
  } else {  ## make some 1-group markets by dropping the second group for the first 'singles' markets
    x <- x[-which((x$m.id %in% (m-singles):m) & (x$g.id %in% seq(2*(m-singles)+2,2*m,2))),]
    return(x)
    ## Two notes for roommates game (to be considered in future version):
    ## 1. there may be more than two groups per market.
    ## 2. roommates' outcomes are not group-level and can differ
  }
}
