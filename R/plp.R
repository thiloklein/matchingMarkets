# ----------------------------------------------------------------------------
# R-code (www.r-project.org/) for the Partitioning Linear Programme
#
# Copyright (c) 2013 Thilo Klein
#
# This library is distributed under the terms of the GNU Public License (GPL)
# for full details see the file LICENSE
#
# ----------------------------------------------------------------------------

#' @title Partitioning Linear Programme for the stable roommates problem 
#'
#' @description Finds the stable matching in the \href{http://en.wikipedia.org/wiki/Stable_roommates_problem}{stable roommates problem} 
#' with transferable utility. 
#' Uses the Partitioning Linear Programme formulated in Quint (1991).
#'
#' @param N integer (divisible by 2) that gives the number of players in the market.
#' @param V valuation matrix of dimension \code{NxN} that gives row-players valuation 
#' over column players (or vice versa).
#' @export
#' @return
#' \code{plp} returns a list with the following items.
#' \item{Valuation.matrix}{input values of V.}
#' \item{Assignment.matrix}{upper triangular matrix of dimension \code{NxN} with entries of 1 for equilibrium pairs and 0 otherwise.}
#' \item{Equilibrium.groups}{matrix that gives the \code{N/2} equilibrium pairs and equilibrium partners' mutual valuations.}
#' @author Thilo Klein 
#' @keywords algorithms
#' @import lpSolve stats
#' @references Quint, T. (1991). Necessary and sufficient conditions for balancedness 
#' in partitioning games. \emph{Mathematical Social Sciences}, 22(1):87--91.
#' @examples
#' ## Roommate problem with 10 players, transferable utility and random preferences:
#' plp(N=10)
#' 
#' ## Roommate problem with 10 players, transferable utility and given preferences:
#' V <- matrix(rep(1:10, 10), 10, 10)
#' plp(V=V)
plp <- function(V=NULL,N=NULL){

  ## Simulate preferences or use given preferences?
  if(is.null(V)==TRUE){
    V = matrix(rnorm(N*N), N, N)
    rownames(V) = colnames(V) = 1:N
    VV = V + t(V)  # mutual match valuations
  } else{
    VV = V + t(V)  # mutual match valuations
    N = dim(VV)[1]
    rownames(VV) = colnames(VV) = 1:N
  }

  ## Consistency checks
  if(N %% 2 == 1){stop("market size needs to be an even number!")}
  if(dim(V)[1] != dim(V)[2]){stop("valuation matrix must be symmetric!")}

  ## Prepare linear system
  f = VV[upper.tri(VV, diag=FALSE)]  # valuation of all (n^2 - n)/2 possible groups.
  b = rep(1,N)  # vector of ones. every of the N players can only be matched once! 
  dir = rep("=",N)
  A = matrix(NA,N,length(f))  # matrix of dim nxf indicating all possible groups player i can belong to
  for(count in 1:N){
    M = matrix(0,N,N)
    M[count,] = 1
    M[,count] = 1
    A[count,] = M[upper.tri(M, diag=FALSE)]
  }

  ## Integer LP solver
  #library(lpSolve)
  so = lp(direction="max", objective.in=f, const.mat=A, const.dir=dir, const.rhs=b, all.bin=TRUE)

  ## Equilibrium groups
  S = matrix(NA, N, N)
  S[upper.tri(S,diag=FALSE)] = so$solution
  G = which(S==1,arr.ind=TRUE)
  G = data.frame(player.A=rownames(VV)[G[,1]], player.B=rownames(VV)[G[,2]], mutual.valuation=VV[G])
  rownames(S) = colnames(S) = colnames(V)

  list(Valuation.matrix=V, Assignment.matrix=S, Equilibrium.groups=G)

}