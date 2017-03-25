#' @title Predict method for fitted matching models
#' 
#' @description Calculate predicted values for matching models fitted with functions \code{stabit} and \code{stabit2}. 
#' 
#' @param object a fitted object of class \code{stabit}
#' @param newdata optionally, a data frame in which to look for variables with which to 
#' predict. If omitted, the fitted linear predictors or the fitted response values are returned.
#' @param ... .
#' 
#' @export
#' 
#' @author Thilo Klein 
#' 
#' @keywords summary
#' 
#' @references Klein, T. (2015a). \href{https://ideas.repec.org/p/cam/camdae/1521.html}{Does Anti-Diversification Pay? A One-Sided Matching Model of Microcredit}.
#' \emph{Cambridge Working Papers in Economics}, #1521.
#' 
#' @examples
#' 
#' ## load the results from Klein (2015) paper
#'  data(klein15a)
#'  
#' ## predict the latent outcome variable
#'  predict(klein15a)
#' 
predict.stabit2 <- function(object, newdata=NULL, ...){
  
  if(object$method=="Klein-selection"){
    stop("Prediction method not yet implemented for matching equations!")
  }
  
  if(is.null(newdata))
    y <- fitted(object)
  else{
    if(!is.null(object$formula)){
      ## model has been fitted using formula interface
      x <- model.matrix(object$formula, newdata)
    }
    else{
      x <- newdata
    }
    y <- as.vector(x %*% object$coefs$alpha)
  }
  y
}



