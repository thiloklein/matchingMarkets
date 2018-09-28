#' @title Stability-Check
#'
#' @description Checks a given two sided matching for blocking pairs.
#'
#' @param nColleges integer indicating the number of colleges
#' @param nStudents integer indicating the number of students
#' @param matching data frame or matrix of dimension (\code{min[nColleges, nStudents]}) x 2 containing in column 1 the colleges and in column 2 the students with each row forming a couple.
#' @param c.prefs matrix of dimension \code{nStudents} x \code{nColleges} with column j containing college j'th ranking over students in decreasing order of preferences.
#' @param s.prefs matrix of dimension \code{nColleges} x \code{nStudents} with column j containing student j'th ranking over colleges in decreasing order of preferences.
#'
#' @export
#'
#' @return \code{stabchk} returns a data frame with as many rows as blocking pairs were found. Column 1 indicates the college and column 2 indicate the student of the blocking pairs. Returns \code{NULL} if no blocking pair is found.
#' @author Thilo Klein, Alexander Sauer
#' @keywords algorithms
#' @examples
#'
#'
#' ## 1-a. Generate preferences for colleges
#' c.prefs = matrix(c(1,2,3,
#'                    3,2,1,
#'                    3,2,1),
#'                     byrow = FALSE, ncol = 3); c.prefs
#'
#' ## 1-b. Generate preferences for students
#' s.prefs = matrix(c(1,2,3,
#'                    3,2,1,
#'                    2,1,3),
#'                  byrow = FALSE, ncol = 3);s.prefs
#'
#' ## 1-c. Generate matching
#' matching = matrix(c(1,2,
#'                     2,1,
#'                     3,3),
#'                   byrow = TRUE, ncol = 2); matching
#'
#' ## 1-d. Check stability
#' stabchk(matching = matching, c.prefs = c.prefs, s.prefs = s.prefs)
#'
#' ## 2-a. Generate new matching without blocking pairs as a data frame
#' matching = data.frame('colleges' = c(1,2,3), 'student' = c(1,3,2))
#' stabchk(matching = matching, c.prefs = c.prefs, s.prefs = s.prefs)
#'
#' ## 3-a. Example with missing values:
#' matching  <- matrix(c(1,1,2,2,3,3), byrow = FALSE, ncol = 2)
#' c.prefs <- matrix(c(1,1,3,rep(NA, 6)), byrow = TRUE, ncol = 3)
#' s.prefs <- matrix(c(2,2,3,rep(NA, 6)), byrow = TRUE, ncol = 3)
#' stabchk(matching = matching, c.prefs = c.prefs, s.prefs = s.prefs)


stabchk <- function(matching, c.prefs, s.prefs, nColleges = ncol(c.prefs), nStudents = ncol(s.prefs)){

  #####
  # TO do:
  # Namen der Spalten der Input dataframes kontrollieren
  # Bsp Formatieren und umschreiben

  blocking_pairs = matrix(nrow = 0, ncol = 2)

  # For every college determine the lowest rank that he was assigned to
  college_ranks <- sapply(1:nColleges, function(college){
    matched_stud <- matching[which(matching[,1]==college),2]
    ranks <- sapply(matched_stud, function(stud){
      match(stud, c.prefs[, college])
    })
    if(all(is.na(ranks))){
      return(nStudents)
    }
    worst_rank <- max(ranks, na.rm = TRUE)
    return(worst_rank)
  })

  for(college in 1:nColleges){
    for(stud in c.prefs[1:college_ranks[college],college]){
      matched_school <- matching[match(stud, matching[,2]),1]
      current_rank_stud <- match(matched_school, s.prefs[, stud])
      possible_rank_stud <- match(college, s.prefs[,stud])

      if(is.na(current_rank_stud) && is.na(possible_rank_stud)){
        next
      }
      if(is.na(current_rank_stud)){
        current_rank_stud <- nColleges
      }
      if(is.na(possible_rank_stud)){
        possible_rank_stud <- nColleges
      }

      if(possible_rank_stud < current_rank_stud){
        blocking_pairs <- rbind(blocking_pairs, c(college, stud))
      }

    }
  }
  if(nrow(blocking_pairs) == 0){
    print('No blocking pairs!')
    return(NULL)
  }
  blocking_pairs_df <- data.frame(blocking_pairs)
  colnames(blocking_pairs_df) <- c('college', 'student')


  blocking_pairs_df <- blocking_pairs_df [!duplicated(blocking_pairs_df),]
  return(blocking_pairs_df)
}

