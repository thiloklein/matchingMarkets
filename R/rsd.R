#' @title Random serial dictatorship mechanism
#'
#' @description Implements the \href{https://en.wikipedia.org/wiki/Random_serial_dictatorship}{random serial dictatorship} algorithm algorithm for a fair division of indivisible objects among individuals. The mechanism takes individuals' prioirty order as an input or alternatively draws a random permutation of the agents form the uniform distribution. Individuals are then successively assigned an object in that order (so the first agent in the ordering gets the first pick and so on).
#
#' @param nIndividuals integer indicating the number of individuals in the matching problem. Defaults to ncol(prefs).
#' @param nObjects integer indication the number of objects in the matching problem. Defaults to nrow(prefs).
#' @param prefs matrix of dimension \code{nObjects} x \code{nIndividuals} with the jth column containing students j's ranking over the objects in decreasing order of preference (i.e. most preferred first).
#' @param priority (Optional) vector of length \code{nIndividuals} indicating the priority of the individuals in decreasing order (i.e. highest individuals first). If none is given, a random order is chosen.
#' @param seed (Optional) integer setting the state for random number generation. Defaults to seed = 123.
#' @param nSlots (Optional) vector of length \code{nObjects} indicating the owners/slots possible/available for each object. (If nSlots only consists of ones this means that every object can only have one owner, if the numbers are higher interpreting objects for example as schools might be more intuitive). Defaults to a vector of length \code{nObjects} filled with ones.
#'
#' @export
#'
#' @return \code{rsd} returns a data frame containing the final matching of individuals (ind) to objects (obj).
#' @author Thilo Klein, Alexander Sauer
#' @keywords algorithms
#' @examples
#' ## Generate preference-matrix
#' prefs <- matrix(c(1,2,3,
#'                   3,1,2,
#'                   1,3,2),
#'                   byrow = FALSE, ncol = 3)
#'
#' priority <- c(1,2,3)
#' nSlots <- c(1,1,1)
#'
#' rsd(prefs = prefs, priority = priority, nSlots = nSlots)




rsd <- function(nIndividuals = ncol(prefs), nObjects = nrow(prefs) ,prefs, priority, seed = NULL, nSlots = rep(1, nObjects)){
  
  if(missing(priority)){
    if(!is.null(seed)){ set.seed(seed) }
    priority <- sample(1:nIndividuals)   # Assign random priority ordering if none is given
  }
  
  #Check dimensions
  if(!(nIndividuals == ncol(prefs) && length(priority)== nIndividuals && length(nSlots)== nObjects)){
    stop('Dimensions do not match')
  }
  
  Res <- data.frame('ind' = NULL, 'obj' = NULL)
  for(i in 1:nIndividuals){
    ind <- priority[i]
    index <- match(TRUE, nSlots[prefs[,ind]] >= 1)
    if(is.na(index)){
      #print(paste('Warning: No Object for individual', as.character(ind)))
      next
    }
    obj <- prefs[,ind][index]
    
    Res <- rbind(Res, data.frame('ind' = ind, 'obj' = obj))
    nSlots[obj] <- nSlots[obj] -1
  }
  match_return <- Res[order(Res$ind),]
  rownames(match_return) <- NULL
  return(match_return)
}
