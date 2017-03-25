#' @export
print.summary.stabit2 <- function(x, ...){
  
  if(x$method=="Klein-selection"){
    
    if(x$mfx==TRUE){
      cat("\nMarginal effects for two-sided matching model.")
    } else{
      cat("\nCoefficients for two-sided matching model.")
    }  
    
  } else{
    
    if(x$mfx==TRUE){
      cat("\nMarginal effects for multi-index sample selection model.")
    } else{
      cat("\nCoefficients for multi-index sample selection model.")
    }
  }
  
  
  if(x$method=="Klein-selection"){
    cat("\nMethod: Klein (2016)\n")
  } else if(x$method=="Klein"){
    cat("\nMethod: Klein (2016), two-sided matching market\n")
  } else if(x$method=="Sorensen"){
    cat("\nMethod: Sorensen (2007), two-sided matching market\n")
  } else{
    cat("\nMethod: Klein (2015), one-sided matching market\n")
  }
  
  cat("\nCall:\n")
  print(x$call)
  
  if(x$method!="Outcome-only"){
    
    if(x$method=="Klein" | x$method=="Klein-selection"){
      
      cat("\nSelection equation (Valuation over colleges):")
      cat("\n")
      printCoefmat(x$college, P.values=TRUE, has.Pvalue=TRUE, signif.legend=FALSE)
      
      cat("\nSelection equation (Valuation over students):")
      cat("\n")
      printCoefmat(x$student, P.values=TRUE, has.Pvalue=TRUE, signif.legend=FALSE)
      
    } else{
      
      cat("\nSelection equation:")
      cat("\n")
      printCoefmat(x$selection, P.values=TRUE, has.Pvalue=TRUE, signif.legend=FALSE) 
    }
  }
  
  if(x$method!="Klein-selection"){
    
    cat("\nOutcome equation:")
    cat("\n")
    printCoefmat(x$outcome, P.values=TRUE, has.Pvalue=TRUE, signif.legend=FALSE)  
  }
  
  cat("---\nSignif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1")
}



