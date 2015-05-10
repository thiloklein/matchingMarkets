# ----------------------------------------------------------------------------
# R-code (www.r-project.org/) for Monte Carlo Experiments
#
# Copyright (c) 2015 Thilo Klein
#
# This library is distributed under the terms of the GNU Public License (GPL)
# for full details see the file LICENSE
#
# ----------------------------------------------------------------------------

#' @title MC Experiments
#'
#' @description MC Experiments
#'
#' @param seed seed value for Monte Carlo Experiment
#' @param niter number of draws in estimation
#' @param N group size (population)
#' @param n group size (sample)
#' @param m number of markets
#' @param type type of the MC Experiment. Either \code{group.members} for randomly sampled group members or \code{counterfactual.groups} for randomly sampled number of counterfactual (or feasible) groups in selection equation (capped at limit max.combs=250)
#' @param method either \code{group.members} or \code{counterfactual.groups}
#' 
#' @export
#' 
#' @author Thilo Klein 
#' 
#' @examples
#' \dontrun{
#' ## 1. Set parameters
#' mciter <- 2 #500
#' niter <- 10 #400000
#' nodes <- 4
#' 
#' ## 2. Setup parallel backend to use 4 processors
#' library(foreach); library(doSNOW)
#' cl <- makeCluster(4); registerDoSNOW(cl)
#' 
#' ## 3. Define foreach loop function
#' mce.add <- function(mciter, niter, N, n, m, type, method){
#'   h <- foreach(i=1:mciter) %dopar% {
#'     library(matchingMarkets)
#'     mce(seed=i,niter, N, n, m, type, method)
#'   }
#'   do.call(rbind, h)
#' }
#' 
#' ## 4. Run siumlations:
#' 
#' ## 4-a. Benchmark study
#' exp.5.5.ols <- mce.add(mciter=mciter, niter=niter, N=5, n=5, m=40,
#'                        type="group.members", method="outcome")
#' exp.5.5.ntu <- mce.add(mciter=mciter, niter=niter, N=5, n=5, m=40, 
#'                        type="group.members", method="NTU")
#' 
#' ## 4-b. Experiment 1: randomly sampled group members
#' exp.6.5.ols <- mce.add(mciter=mciter, niter=niter, N=6, n=5, m=40,
#'                        type="group.members", method="outcome")
#' exp.6.5.ntu <- mce.add(mciter=mciter, niter=niter, N=6, n=5, m=40,
#'                        type="group.members", method="NTU")
#' 
#' ## 4-c. Experiment 2: randomly sampled counterfactual groups
#' exp.6.6.ols <- mce.add(mciter=mciter, niter=niter, N=6, n=6, m=40, 
#'                        type="counterfactual.groups", method="outcome")
#' exp.6.6.ntu <- mce.add(mciter=mciter, niter=niter, N=6, n=6, m=40, 
#'                        type="counterfactual.groups", method="NTU")
#' 
#' ## 5. Stop parallel backend
#' stopCluster(cl)
#' }
mce <- function(seed, niter, N, n, m, type, method){
    
    ######################
    ## A. Generate data ##
    ######################
    
    #library(matchingMarkets)
    l <- N-n # number of individuals excluded by random sampling of 'n' of 'N' members
    
    ## Individual-level independent variables are fixed (seed=123) in repeated samples
    idata <- stabsim(m=m, ind=N, seed=123, gpm=2)
    
    if(l == 0 & type != "counterfactual.groups"){ 
      ## If sample group size 'n' equals population group size 'N'
      ## and experiment is not to drop counterfactual groups
      ## -> run benchmark experiment
      
      ## Group outcomes and equilibrium group constellation are dependent on fixed 
      ## individual variables of group members and draws -- set.seed(seed) -- of unobservables    
      #! set.seed(seed) # set seed for draw of error terms
      #! idata[,c("xi.i","eta.i")] <- matrix(rnorm(2*dim(idata)[1]), ncol=2) # draw new xi and eta
      mdata <- stabit(x=idata, selection = list(ieq="wst"), 
                      outcome = list(ieq="wst"), simulation="NTU", 
                      method=method, niter=niter, seed=seed, censored=1)
      
      h <- mdata$coefs[!names(mdata$coefs) %in% "eta"]
      h <- unlist(lapply(h, function(x) x[,1]))
      
    } else{
      ## If sample group size 'n' is lower than population group size 'N'
      ## -> 'l' individuals dropped
      
      ## Group outcomes and equilibrium group constellation are dependent on fixed 
      ## individual variables of group members and draws -- set.seed(seed) -- of unobservables    
      #! set.seed(seed) # set seed for draw of error terms
      #! idata[,c("xi.i","eta.i")] <- matrix(rnorm(2*dim(idata)[1]), ncol=2) # draw new xi and eta
      mdata <- stabit(x=idata, selection = list(ieq="wst"), 
                      outcome = list(ieq="wst"), simulation="NTU", 
                      method="model.frame", seed=seed)
      
      ####################################
      ## B-1. Define equilibrium groups ##
      ####################################
      
      for(j in 1:length(mdata$model.list$E)){ ## for all markets
        
        ## add individual indices (E), group indicators (g.id) and outcomes (R)
        
        # rearrange individuals by equilibrium group affiliation
        idata[idata$m.id==j,]       <- idata[idata$m.id==j,][c(t(mdata$model.list$E[[j]])),]  
        # add simulated group outcome
        idata[idata$m.id==j,"R"]    <- c(rep(mdata$model.list$R[[j]][1],N) ,rep(mdata$model.list$R[[j]][2],N))  
        # add new group id
        idata[idata$m.id==j,"g.id"] <- sort(idata[idata$m.id==j,"g.id"])  
        
        #############################
        ## B-2. Drop group members ##
        #############################
        
        if(type == "group.members"){
          
          ## drop randomly chosen l=1,2,3 members per group
          for(k in unique(idata$g.id[idata$m.id==j])){
            dodrop  <- which(idata$g.id==k)[sample(N)>n]
            idata <- idata[-dodrop,]
          } 
        }
      }        
      
      ###################
      ## C. Estimation ##
      ###################
      
      if(type == "group.members"){
        
        ## Generates independent variables for sampled group members, 
        ## keeping outcomes fixed (simulation="none")
        mdata <- stabit(x=idata, selection = list(ieq="wst"), 
                        outcome = list(ieq="wst"), simulation="none", 
                        method=method, niter=niter, censored=1)
        
      } else if(type == "counterfactual.groups"){
        
        ## Generates independent variables for sampled group members, 
        ## keeping outcomes fixed (simulation="none"). 
        ## Samples from counterfactual groups (max.combs)
        mdata <- stabit(x=idata, selection = list(ieq="wst"), 
                        outcome = list(ieq="wst"), simulation="none", 
                        method=method, niter=niter, max.combs=250, seed=123, censored=1)
      }
      
      h <- mdata$coefs[!names(mdata$coefs) %in% "eta"]
      h <- unlist(lapply(h, function(x) x[,1]))
    }
}

