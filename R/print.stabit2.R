#' @export
print.stabit2 <- function(x, ...){
  
  cat("Call:\n")
  print(x$call)
  cat("\nCoefficients:\n")
  print(x$coefficients)
}



