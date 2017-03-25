#' @export
plot.stabit2 <- function(x, ...){
  
  mcmcPlots <- function(dat, main){
    
    res <- as.data.frame(t(dat))
    res$iteration <- 1:nrow(res)
    res.long <- data.frame( iteration=rep(1:nrow(res),ncol(res)-1), 
                            condition=c(sapply(names(res)[-ncol(res)], function(z) rep(z,nrow(res)))), 
                            measurement=c(unlist(res[,1:(ncol(res)-1)])) )
    
    lattice.options(default.args=list(as.table=TRUE), 
                    default.theme=standard.theme(color=FALSE))
    keys <- list(text=main, space="top", columns=1, lines=TRUE, points=FALSE)
    plotOut <- xyplot(measurement ~ iteration | factor(condition, levels=unique(condition)), 
                      data = res.long, scales=list(relation="free"), auto.key=keys,
                      xlab = "iterations", ylab = "paramter draws", type = "l") 
    plotOut
  }
  
  
  if(x$method=="Klein-selection"){
    
    plot( mcmcPlots(dat = x$draws$betadraws, main = "MCMC draws in college selection equation") )
    par(ask=TRUE)
    plot( mcmcPlots(dat = x$draws$gammadraws, main = "MCMC draws in student selection equation") )
    par(ask=FALSE)
    
  } else if(x$method=="Klein"){
    
    plot( mcmcPlots(dat = x$draws$alphadraws, main = "MCMC draws in outcome equation") )
    par(ask=TRUE)
    plot( mcmcPlots(dat = x$draws$betadraws, main = "MCMC draws in college selection equation") )
    par(ask=TRUE)
    plot( mcmcPlots(dat = x$draws$gammadraws, main = "MCMC draws in student selection equation") )
    par(ask=FALSE)
    
  } else if(x$method %in% c("Sorensen","Selection")){
    
    plot( mcmcPlots(dat = x$draws$alphadraws, main = "MCMC draws in outcome equation") )
    par(ask=TRUE)
    plot( mcmcPlots(dat = x$draws$betadraws, main = "MCMC draws in selection equation") )
    par(ask=FALSE)
    
  } else if(x$method=="Outcome-only"){
    
    plot( mcmcPlots(dat = x$draws$alphadraws, main = "MCMC draws in outcome equation") )
    
  }
}



