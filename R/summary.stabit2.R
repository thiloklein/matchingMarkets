#' @export
summary.stabit2 <- function(object, mfx=FALSE, ...){
  
  ## function to produce regression tables
  tab <- function(coefs, vcov){
    se <- sqrt(diag(vcov))
    tval <- coefs / se
    TAB <- cbind(Estimate = coefs,
                 StdErr = se,
                 t.value = tval,
                 p.value = 2*pt(-abs(tval), df=object$df))
    TAB
  }
  
  res <- list()
  
  ## Selection equation(s)
  
  if(object$method!="Outcome-only"){
    
    if(mfx==FALSE){
      
      if(object$method=="Klein" | object$method=="Klein-selection"){
        
        res$college <- tab(coefs=object$coefs$beta, vcov=object$vcov$beta)
        
        res$student <- tab(coefs=object$coefs$gamma, vcov=object$vcov$gamma)
        
      } else{
        
        res$selection <- tab(coefs=object$coefs$beta, vcov=object$vcov$beta)
      }
      
    } else if(mfx==TRUE){
      
      if(object$method=="Klein" | object$method=="Klein-selection"){
        
        res$college <- mfxVal(coefs = object$coefs$beta, 
                              vcov = object$vcov$beta,
                              df = nrow(object$variables$C) - ncol(object$variables$C) )
        
        res$student <- mfxVal(coefs = object$coefs$gamma, 
                              vcov = object$vcov$gamma,
                              df = nrow(object$variables$S) - ncol(object$variables$S) )
        
      } else if(object$method=="Sorensen"){
        
        res$selection <- mfxVal(coefs = object$coefs$beta, 
                                vcov = object$vcov$beta,
                                df = nrow(object$variables$C) - ncol(object$variables$C) )
      } else{
        
        res$selection <- mfxVal(coefs = object$coefs$beta, 
                                vcov = object$vcov$beta,
                                df = nrow(object$variables$W) - ncol(object$variables$W) )
      }
    }
  }
  
  
  ## Outcome equation
  
  if(object$method!="Klein-selection"){
    
    if(mfx==FALSE){
      
      res$outcome <- tab(coefs=object$coefs$alpha, vcov=object$vcov$alpha)
      
    } else{
      
      if(object$binary==TRUE){
        
        res$outcome <- mfxOut(sims = 10000, X = object$variables$X,
                              coefs = object$coefs$alpha,
                              vcov = object$vcov$alpha,
                              df = nrow(object$variables$X) - ncol(object$variables$X) )
      } else{
        
        res$outcome <- tab(coefs = object$coefs$alpha,
                           vcov = object$vcov$alpha )
      }
    } 
  }
  
  res$call <- object$call
  res$method <- object$method
  res$mfx <- mfx
  
  class(res) <- "summary.stabit2"
  res
}


mfxOut <- function(sims=10000, x.mean=TRUE, coefs, vcov, X, df){
  
  ## source: http://researchrepository.ucd.ie/handle/10197/3404
  ## method: average of individual marginal effects at each observation
  ## interpretation: http://www.indiana.edu/~statmath/stat/all/cdvm/cdvm.pdf page 8
  
  set.seed(1984)
  se <- sqrt(diag(vcov))
  
  if(x.mean==TRUE){
    
    ## marginal effects are calculated at the means of independent variables
    pdf <- dnorm(mean(X%*%coefs))
    pdfsd <- dnorm(sd(X%*%coefs))
    
  } else{
    
    ## marginal effects are calculated for each observation and then averaged
    pdf <- mean(dnorm(X%*%coefs))
    pdfsd <- sd(dnorm(X%*%coefs))
  }  
  mx <- pdf*coefs
  
  sim <- matrix(rep(NA,sims*length(coefs)), nrow=sims)
  
  for(i in 1:length(coefs)){
    sim[,i] <- rnorm(sims,coefs[i],se[i])
  }
  
  pdfsim <- rnorm(sims,pdf,pdfsd)
  sim.se <- pdfsim*sim
  s.e. <- apply(sim.se,2,sd)
  
  tval <- mx / s.e.
  TAB <- cbind(Estimate = mx,
               StdErr = s.e.,
               t.value = tval,
               p.value = pt(-abs(tval), df=df))
  TAB
}


mfxVal <- function(coefs, vcov, df){
  
  ## Reference: Sorensen (2007, p. 2748)
  se <- sqrt(diag(vcov))
  mx <- dnorm(0)*coefs/sqrt(2)
  s.e. <- dnorm(0)*se/sqrt(2)
  tval <- mx / s.e.
  TAB <- cbind(Estimate = mx,
               StdErr = s.e.,
               t.value = tval,
               p.value = pt(-abs(tval), df=df))
  TAB
}



