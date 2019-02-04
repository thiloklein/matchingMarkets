#' @title Top-Trading-Cycles Algorithm with existing tenants
#'
#' @description Implements an algorithm for the
#' \href{https://en.wikipedia.org/wiki/Top_trading_cycle}{house allocation problem} proposed by Abdulkadiroglu and Sonmez (1999) for a matching problem in which there are both vacant houses and existing tenants.
#'
#' @param nStudents integer indicating the number of students. Defaults to \code{ncol(s.prefs)}.
#' @param nHouses integer indicating the number of houses. Defaults to \code{length(houses)}.
#' @param s.prefs matrix of dimension \code{nHouses} x \code{nStudents} with column j containing student jth ranking over houses in decreasing order of preferences (i.e. most preferred first).
#' @param houses vector of length \code{nHouses} which represents the occupation of the houses. Entry in \code{k} contains \code{j} if student \code{j} is living in house \code{k} and \code{NA} if house \code{k} is vacant.
#' @param priority (Optional) vector of length \code{nStudents}. Gives the prioirity ordering of the students in the search for cycles (Do not confuse it with the preferences!), if nothing is specified a random ordering is chosen.
#' @param seed (Optional) integer setting the state for random number generation. Defaults to seed = NULL
#'
#' @export
#'
#' @return \code{ttc} returns a data frame of the matching of students (int) to houses (obj)  for the house allocation problem based on the Top-Trading-Cycles algorithm.
#' @author Thilo Klein, Alexander Sauer
#' @keywords algorithms
#' @references Abdulkadiroglu, A. and T. Sonmez (1999). House Allocation with Existing Tenants. \emph{Journal of Economic Theory},  88 (2): 233-260.
#' @references Shapley, L. and H. Scarf (1974). On Cores and Indivisibility. \emph{Journal of Mathematical Economics}, 1(1): 23-37.
#' @examples
#' ##\dontrun{
#' ## 1-a. Generate matrix of individuals' preference rankings over objects,
#' ## a.k.a. Rank Order Lists (ROL).
#' s.prefs <- matrix(c(3,2,4,1,        # ROL of student 1
#'                    3,5,6, NA,
#'                    3,1, NA,NA,
#'                    2,5,6,4,
#'                    1,3,2,NA,
#'                    2,4,5,6), nrow = 4, ncol = 6, byrow = FALSE); s.prefs
#'
#' ## 1-b. Generate vector of house occupation objects ('obj') and their owners ('ind')
#' (houses <- 1:6)
#'
#' ## 1-c. Find assignment based on TTC algorithm
#' ttc(s.prefs = s.prefs, houses = houses, nHouses = 6, priority = 1:6)
#'
#' ## 2-a.Compare the example in the paper Abdulkadiroglu et al. (1999)
#' ## on page 246-248 (section 5.1 An Example):
#' ## generate matrix of students' preference rankings over houses, a.k.a. Rank Order Lists (ROL)
#' s.prefs <- matrix(c(2,6,5,1,4,3,7,NA,
#'                  7,1,6,5,4,3,2,NA,
#'                  2,1,4,7,3,6,5,NA,
#'                  2,4,3,6,1,7,5,NA,
#'                  4,3,7,1,2,5,6,NA), byrow = FALSE, ncol= 5); s.prefs
#'
#' ## 2-b. Generate house occupation, so student 1 lives in house 1, ..., student 4 lives in house 4
#' ## and the other houses are vacant.
#' houses <- c(1,2,3,4,NA,NA,NA,NA); houses
#'
#' ## 2-c. Generate priority ordering
#' priority <- 1:5
#'
#' ## 2-d. Find assigment
#' ttc(s.prefs = s.prefs, houses = houses, priority = priority)
#' ##}

ttc <- function(nStudents = ncol(s.prefs), nHouses = length(houses), s.prefs, houses, priority = NULL, seed = NULL){
  
  # Check if prioirty is given:
  if(missing(priority)){
    set.seed(seed)
    priority <- sample(1:nStudents)   # Assign random priority ordering if none is given
  }
  
  #Check dimensions
  if(!(nStudents == ncol(s.prefs)  )){
    stop('Number of students does not match the dimension of s.prefs')
  }
  if(!(nStudents == length(priority) && nHouses == length(houses))){
    #print(nStudents)
    #print(length(priority))
    #print(nHouses)
    #print(length(houses))
    stop('Number of students/houses does not match the length of priority ordering or houses')
  }
  
  matched_ind = rep(FALSE, nStudents)
  matched_houses = rep(FALSE, nHouses)
  Res <- data.frame('ind' = NULL, 'house' = NULL)
  
  # Run TTC:
  while(!(all(matched_ind) || all(matched_houses))){  #Are all individuals matched or no houses left?
    
    # 1. find cycle:
    Cycle <- find_cycle_extenants(P = s.prefs, H = houses, ordering = priority, matched_ind = matched_ind, matched_houses = matched_houses)
    
    # 2. Update Result
    Res <- rbind(Res,Cycle)
    
    # 3. Update matched individuals and houses
    matched_houses[Cycle[,2]] <- TRUE
    matched_ind[Cycle[,1]] <- TRUE
  }
  return(Res[order(Res$ind),])
}


find_cycle_extenants <- function(P, H, ordering, matched_ind, matched_houses){
  Cycle   <- data.frame(ind=NA, obj=NA)
  
  this_ind_orderindex <- match(FALSE, matched_ind[ordering]) # Individual that is highest in the priority list and not already matched
  this_ind <- ordering[this_ind_orderindex]
  
  for(j in 1:dim(P)[1]){
    # Highest ranked house that is not already matched and most wanted by this_ind
    index_wanted_house <- match(FALSE, matched_houses[P[,this_ind]])
    wanted_house <- P[,this_ind][index_wanted_house]
    
    Cycle[j,] <- c(this_ind, wanted_house)
    
    # Does the house have a tenant? Or is the tenant already matched?
    owner <- H[wanted_house]
    if(is.na(owner) == FALSE && matched_ind[owner] == FALSE){
      this_ind <- owner   # If the house owner has not been matched, the house points to its tentant
    }else{
      # Choose the indvidual with the highest priority that is not already matched
      index_ind <- match(FALSE, matched_ind[ordering])
      this_ind <- ordering[index_ind]
    }
    
    if(this_ind %in% Cycle$ind){
      h <- which(this_ind == Cycle$ind)
      return(Cycle[h:length(Cycle$ind),])
      break
    }
  }
}