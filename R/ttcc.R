#' @title Top-Trading-Cycles and Chains Algorithm
#'
#' @description Implements the Top Trading Cycle and Chains algorithm proposed by Roth et al. (2004) for the kidney exchange problem.
#' The algorithm requires a rule to determine which chain will be used if there is more than one possibility. The chosen rule is
#' to search for the longest chain and remove it from the problem (even the first kidney which was unassigned).
#'
#' @param nPatients integer indicating the number of patient/donor-pairs in the matching problem. Defaults to  \code{ncol(prefs)}.
#' @param prefs matrix of dimension (\code{nPatients} + 1) x \code{nPatients} with column j containg patients jth ranking over kidneys in decreasing order of preferences (i.e. most preferred first). An entry with value (\code{nPatients} +1) indicates that the patient prefers the waiting list to all kidney below in his ranking (therefore they do not matter and can be neglected/NA).
#' @param priority (Optional) Vector of length \code{nPatients}. Gives the prioirty ordering of the patient, if nothing is specified a random ordering is chosen.
#' @param seed (Optional) integer setting the state for random number generation. Defaults to seed = 123.
#'
#' @export
#'
#' @return \code{ttcc} returns a list with the matching and a vector containing the patients who are assigned to the waiting list. The matching comprises a data frame of the operations to be performed between patient-donor pairs (ind-obj).
#' @author Thilo Klein, Alexander Sauer
#' @keywords algorithms
#' @references Roth, A.; T. Sonmez; U. Unver (2004). Kidney Exchange. \emph{Quarterly Journal of Economics}, 119 (2): 457-488.
#' @examples
#' \dontrun{
#' ## Compare Example 1 from Roth et al. (2004) on page 469 - 475
#' ## generate matrix of patients' preference rankings over kidneys, a.k.a. Rank Order Lists (ROL)
#'
#' prefs <- matrix(c( 9,10, 1,NA,NA,NA,NA,
#'                   11, 3, 5, 6, 2,NA,NA,
#'                    2, 4, 5, 6, 7, 8,13,
#'                    5, 9, 1, 8,10, 3,13,
#'                    3, 7,11, 4, 5,NA,NA,
#'                    3, 5, 8, 6,NA,NA,NA,
#'                    6, 1, 3, 9,10, 1,13,
#'                    6, 4,11, 2, 3, 8,NA,
#'                    3,11,13,NA,NA,NA,NA,
#'                   11, 1, 4, 5, 6, 7,13,
#'                    3, 6, 5,11,NA,NA,NA,
#'                   11, 3, 9, 8,10,12,NA),
#'               byrow = FALSE, ncol = 12); prefs
#' priority <- 1:12
#' ttcc(prefs = prefs, priority = priority)
#' ## The final matching differs slightly because in Round 3 another chain is chosen due to a different
#' ## decision rule (compare Figure 3, p472. Here W1 instead of W2 is chosen)
#' }

ttcc <- function(nPatients = ncol(prefs), prefs, priority = NULL, seed = 123){
  # Chain Rule:   Searches for the longest chain and removes it from the problem (even the first kidney which was unassigned)

  if(missing(priority)){
    set.seed(seed)
    priority = sample(1:nPatients)
  }
  prefs <- prefs[, priority]

  Res <- data.frame('ind' = NULL, 'obj' = NULL)
  matched_kid <- rep(FALSE, nPatients+1)
  waiting_list <- rep(FALSE, nPatients+1)

  while(!all(matched_kid | append(waiting_list[1:(length(waiting_list)-1)], TRUE))){   # The last entry of the waiting list has to be changed to true, to avoid an infinite loop

    cycle_chain <- find_cycle_chain(prefs, matched_kid, waiting_list, priority)

    if(cycle_chain$type == 'cycle'){
      #Update Result
      Res <- rbind(Res, cycle_chain$obj )
      #Update matched kidneys
      matched_kid[cycle_chain$obj[,1]] <- TRUE
    }
    if(cycle_chain$type == 'chain'){   #Only a chain was found
      chain <- cycle_chain$obj

      # Last individual on the wating list
      last <- tail(chain[,1],1)
      waiting_list[last] <- TRUE

      if(dim(chain)[1] > 1){
        # Update matchings of the chain:
        other <- chain[1:(dim(chain)[1]-1),]

        Res <- rbind(Res, other)
        matched_kid[other[,1]] <- TRUE
      }
    }
  }
  match_return <- Res[order(Res$ind),]
  rownames(match_return) <- NULL
  Result = list('matching' = match_return, 'waiting_list' = which(waiting_list[1:(length(waiting_list)-1)] == TRUE))
  return(Result)
}



find_cycle_chain <- function(prefs, matched, waiting_list, priority){
  # Takes a preference matrix, a vector indicating which patients have already been matched and one indicating which have already been assigned to the waiting list
  # Go through all individuals and try to find a cycle or a chain
  # If a cycle is found, return it and stop
  # If a chain is found, check if it is longer than all previous found chains, and start with the next individual and try to find a longer chain
  # Return the cycle or the longest chain as the second element of a list, first element indicates whether it is a cycle or chain

  nPatients <- dim(prefs)[2]

  Chain <- data.frame('ind' = NA, 'obj' = NA)
  for(k in priority){

    this_ind <- k
    # Go through priority list until we have a unmatched individuum that is not already on the waiting list
    if(matched[this_ind] == TRUE | waiting_list[this_ind] == TRUE){ next }

    Cycle <- data.frame('ind' = NA, 'obj' = NA)

    for(i in 1:nPatients){
      wanted_kidney_index <- match(FALSE, matched[prefs[,this_ind]] | waiting_list[prefs[,this_ind]])
      if(is.na(wanted_kidney_index)){
        print(k)
        print(matched)
        print(waiting_list)
        stop('NO MATCH FOUND')
      }
      wanted_kidney <- prefs[,this_ind][wanted_kidney_index]

      Cycle[i,] <- c(this_ind, wanted_kidney)

      # Patient points to the waiting list -> Chain
      if(wanted_kidney == (nPatients+1)){
        # Is it a longer chain?
        if(dim(Cycle)[1] >= dim(Chain)[1]){
          Chain <- Cycle
        }
        break
      }

      this_ind <- wanted_kidney

      if(this_ind %in% Cycle$ind){
        h <- which(this_ind == Cycle$ind)
        Result <- list('type' = 'cycle', 'obj' = Cycle[h:length(Cycle$ind),])
        return(Result)
        break
      }
    }
  }
  Result <- list('type' = 'chain', 'obj' = Chain)
  return(Result)
}

