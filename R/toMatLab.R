# ----------------------------------------------------------------------------
# R-code (www.r-project.org/) to write the data to MAT v5 format using
# the R.matlab package
#
# Copyright (c) 2013 Thilo Klein
#
# This library is distributed under the terms of the GNU Public License (GPL)
# for full details see the file LICENSE
#
# ----------------------------------------------------------------------------

#' @title Write one-sided matching data to Matlab .mat file
#'
#' @description Write an object from \code{design.matrix} to MAT v5 format using the \code{\link[R.matlab:R.matlab-package]{R.matlab}} package
#'
#' @param x one-sided matching data after applying the \code{design.matrix()} function.
#' @param folder path that .mat file will be saved to. Defaults to current working directory \code{getwd()}.
#' @param id number or character to be added to file name. The default is none.
#' @export
#' @author Thilo Klein 
#' @keywords generate
#' @import R.matlab
#' @examples
#' \dontrun{
#' ## Load one-sided matching data
#' data(mdata, package="matchingMarkets")
#' 
#' ## Write data to file mdata123.mat on Desktop
#' toMatLab(x=mdata, folder="~/Desktop/", id=123)}
toMatLab <- function(x, folder=getwd(), id=""){
  
  # --------------------------------------------------------------------
  # R-code (www.r-project.org) to write the data to MAT v5 format using
  # the R.matlab package
  
  # The arguments of the function are:
  # x : matching data after applying the OneSidedMatching function
  # folder : path that .mat file will be saved to
  # id : number or character to be added to file name
  
  # ## Examples:
  #
  # toMatLab(x, folder="~/Desktop/", id=001)
  # --------------------------------------------------------------------
  
  #library(R.matlab)
  filepath <- paste(folder,deparse(substitute(x)),id,".mat",sep="")
  
  ## replace NA's with 0 in combinations matrices
  x$combs <- lapply(x$combs, function(j) replace(j,is.na(j),0))

  ## add intercept to design matrix of outcome equation
  x$X <- lapply(x$X, function(i) data.frame(intercept=1,i))

  ## variable names
  an <- colnames(x$W[[1]])
  bn <- colnames(x$X[[1]])

  ## rename the nth market with 'n' for all vars
  for(i in names(x)){
    names(x[[i]]) <- as.character(1:length(x[[i]]))
  }

  ## write to MAT v5 format
  writeMat(filepath, an=an, bn=bn, D=x$D, R=x$R, W=x$W, X=x$X, P=x$P, combs=x$combs)

}
