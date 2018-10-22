#' @title Top-Trading-Cycles Algorithm for a two sided matching problem
#'
#' @description Implements the school matching algorithm proposed in Abdulkadiroglu and Sonmez (2003) for a matching problem
#' in which both sides have preferences. Missing preferences are handled in the following ways: Suppose that a student only ranked colleges that are already matched
#' to other students. This student is removed from the matching process and a list with all unmatchable students is printed.
#' If \code{full_return} is set to \code{TRUE}, a vector with this students is returned as well.
#' Now suppose during the matching process a student points to a college that still has capacities but does not rank any more students.
#' We assume now that the college is indifferent over all other students (so we do not allow for free capacieties) and we match the student who wants to go there to the college.
#'
#' @param nStudents integer indicating the number of students in the matching problem. Defaults to \code{ncol(s.prefs)}.
#' @param nColleges integer indicating the number of colleges in the matching problem. Defaults to \code{ncol(c.prefs)}.
#' @param s.prefs matrix of dimension \code{nColleges} x \code{nStudents} with the jth column containing student j's ranking over colleges in decreasing order of preference (i.e. most preferred first).
#' @param c.prefs matrix of dimension \code{nStudents} x \code{nColleges} with the ith column containing college i's ranking over students in decreasing order of preference (i.e. most preferred first).
#' @param nSlots vector of length \code{nColleges} indicating the number of places (i.e. quota) of each college.
#' @param priority (Optional) vector of length \code{nStudents}. Gives the prioirity ordering of the students, if nothing is specified a random ordering is chosen.
#' @param seed (Optional) integer setting the state for random number generation. Defaults to seed =123.
#' @param full_return (Optinal) If \code{TRUE} the return value is a list with the matching, the remaining seats and the unmatchable students is returned. Defaults to \code{FALSE} and only the matching is returned.
#'
#' @export
#'
#' @return \code{ttc2} returns a data frame of the matching of students (ind) to colleges (obj) for the school market problem based on the Top-Trading-Cycles algorithm.
#' @author Thilo Klein, Alexander Sauer
#' @keywords algorithms
#' @references Abdulkadiroglu, A. and T. Sonmez (2003). School Choice: A Mechanism Design Approach. \emph{American Economic Review}, 93 (3): 729-747.
#' @examples
#' \dontrun{
#' ## 1-a. Compare example from the Abdulkadiroglu et al. (2003) (in the Appendix, page 742-744)
#' ## 1-b. Generate matrix of students' preference rankings over schools, a.k.a. Rank Order Lists (ROL)
#' s.prefs <- matrix(c(
#'                   2,1,3,4,
#'                   1,2,3,4,
#'                   3,2,1,4,
#'                   3,4,1,2,
#'                   1,3,4,2,
#'                   4,1,2,3,
#'                   1,2,3,4,
#'                   1,2,4,3),
#'                   byrow = FALSE, ncol = 8); s.prefs
#'
#' ## 1-c. Generate matrix of schools' preference rankings over students, a.k.a. Rank Order Lists (ROL)
#' c.prefs <- matrix(c(
#'                   1,2,3,4,5,6,7,8,
#'                   3,5,4,8,7,2,1,6,
#'                   5,3,1,7,2,8,6,4,
#'                   6,8,7,4,2,3,5,1),
#'                   byrow = FALSE, ncol = 4); c.prefs
#'
#' ## 1-d. Generate capacities
#' nSlots <- c(2,2,3,3)
#'
#' ## 1-e. Find assignment based on TTC algorithm
#' ttc2(s.prefs = s.prefs, c.prefs = c.prefs, nSlots = nSlots)
#'
#' ## 2-a. Generate college preferences with college 1 only ranking student 1
#' c.prefs <- matrix(c(
#'                    1,rep(NA,7),
#'                    3,5,4,8,7,2,1,6,
#'                    5,3,1,7,2,8,6,4,
#'                    6,8,7,4,2,3,5,1),
#'                    byrow = FALSE, ncol = 4); c.prefs
#'
#' ## 2-b. Find assignment based on TTC algorithm
#' ttc2(s.prefs = s.prefs, c.prefs = c.prefs, nSlots = nSlots, priority = 1:8)
#'
#' ## If all schools have the same preferences the two sided ttc and the serial dictator yield
#' ## the same outcome if the preferences are taken to be the prioirty order for the serial dictator
#'
#' # Preferences are the same for all schools:
#' c.prefs <- matrix(c(
#'                   5,3,1,7,2,8,6,4,
#'                   5,3,1,7,2,8,6,4,
#'                   5,3,1,7,2,8,6,4,
#'                   5,3,1,7,2,8,6,4),
#'                   byrow = FALSE, ncol = 4)
#' priority <- c.prefs[,1]
#'
#' match_ttc <- ttc2(s.prefs = s.prefs, c.prefs = c.prefs, nSlots = nSlots); match_ttc
#' match_sd <- rsd(prefs = s.prefs, priority = priority, nSlots = nSlots); match_sd
#' all(match_ttc == match_sd)
#' }

ttc2 <- function( nStudents = ncol(s.prefs), nColleges = ncol(c.prefs), s.prefs = NULL,  c.prefs = NULL,nSlots = NULL, priority = NULL , seed = 123, full_return = FALSE){

  ## To Do
  ## Slots in find_cycle aktualisieren?

  # Check priority
  if(missing(priority)){
    set.seed(seed)
    priority <- sample(1:nStudents)   # Assign random priority ordering if none is given
  }
  # Check dimensions
  if(!(nStudents == ncol(s.prefs) && nColleges == ncol(c.prefs))){
    stop('Number of students/colleges do not match dimensions')
  }
  if(!(length(nSlots) == nColleges && nrow(c.prefs == ncol(s.prefs) ))){
    stop('Dimensions of nSlots and/or preference matrices do not match!')
  }


  # Start matching:
  matched_stud = rep(FALSE, nStudents)
  unmatchable_stud = c()
  Res <- data.frame('ind' = NULL, 'obj' = NULL)

  # Loop as long as there is an unmatched student
  while(!all(matched_stud)){

    if(all(nSlots <= 0)){
      print('Not enough capacity!')
      break
    }
    # Find Cycle
    Cycle <- find_cycle_two_side(s.prefs = s.prefs, c.prefs = c.prefs, matched_stud = matched_stud, nSlots = nSlots, priority = priority)

    # Check if find_cyle returned that there was an individual which could not be matched: In this case find_cycle_two_side returns the index of the unmatchable stud
    if(is.numeric(Cycle)){
      # Cycle is the index of the student that was not matchable
      unmatchable_stud <- append(unmatchable_stud, Cycle)
      # Remove student from the following matching process
      matched_stud[Cycle] <- TRUE
      next
      # Start again and try to find a new cycle
    }
    # Update Matching
    Res <- rbind(Res, Cycle)

    # Update Capacity and matched students
    nSlots[Cycle[,2]] <- nSlots[Cycle[,2]] -1
    matched_stud[Cycle[,1]] <- TRUE
  }

  if(length(unmatchable_stud) > 0){
    print('It was not possible to match the following student(s) according to their preference rankings:')
    print(unmatchable_stud)
  }

  match_return <- Res[order(Res$ind),]
  rownames(match_return) <- NULL

  if(full_return){
    output <- list('matching' = match_return, 'nSlots' = nSlots, 'unmatchable_students' = unmatchable_stud )
    return(output)
  }

  return(match_return)
}



find_cycle_two_side <-function(s.prefs = NULL, c.prefs = NULL, matched_stud = NULL, nSlots, priority){
  # Takes Preferences of students and schools, vector with remaining capacities as well as a vector indicating which students have already been matched.
  # Returns data frame with Cycle of Students and schools OR a numeric indicating the student which could not be matched
  #####
  current_slots <- nSlots

  Cycle <- data.frame('ind' = NA, 'obj' = NA)
  this_ind_prior_index <- match(FALSE, matched_stud[priority]) # First Student that still has to be matched in order of the given priority
  this_ind <- priority[this_ind_prior_index]

  if(is.na(this_ind)){
    stop('All students are matched!')
  }

  for(j in 1:dim(s.prefs)[1]){

    # Find college with the highest preference ranking that still has capacitie
    index_wanted_school <- match(TRUE, nSlots[s.prefs[,this_ind]]>= 1) # Order the nSlots-vector by the preferences of this ind and take the first instance with free capacity
    #index_wanted_school <- match(TRUE, current_slots[s.prefs[,this_ind]]>= 1) # Order the current_Slots-vector by the preferences of this ind and take the first instance with free capacity
    wanted_school = s.prefs[,this_ind][index_wanted_school]

    # If index_wanted_school is NA, then there was no school ranked by the current individual that still has capacities
    if(is.na(index_wanted_school)){
      return(this_ind)
    }

    Cycle[j,] <- c(this_ind, wanted_school)

    # Adjust slots
    current_slots[wanted_school] <- current_slots[wanted_school] - 1

    # Find student with highest priority that is not already matched
    index_wanted_stud <- match(FALSE, matched_stud[c.prefs[,wanted_school]])
    this_ind <- c.prefs[,wanted_school][index_wanted_stud]

    # Suppose a student points to a college that still has capacities but does not rank any more students.
    # Then we assume that the college is indifferent over all other students. Therefore we can assume that it points to the first
    # individual in the cycle and the cycle is therefore closed:
    if(is.na(index_wanted_stud)){
      return(Cycle)
    }

    if(this_ind %in% Cycle$ind){
      h <- which(this_ind == Cycle$ind)
      ##########
      #print(Cycle)
      #print(this_ind)
      return(Cycle[h:length(Cycle$ind),])
      break
    }
  }
}

