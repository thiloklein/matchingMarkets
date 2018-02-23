#' @export
summary.hri2 <- function(object, ...) {
  message("Matching Summary Stats:")
  message("# Unmatched Singles: ", object$stats$unmatchedSingles)
  message("# Unmatched Couples: ", object$stats$unmatchedCouples)
  message("# Unmatched Program slots: ", object$stats$unmatchedProgrammSlots)
  message("# Ave Resident Rank of their matching: ", object$stats$aveRRank)
  message("# Num Residents getting their top rank: ", object$stats$numRTop)
  message("# Num Couples getting their top rank: ", object$stats$numCTop)
  message("# Ave Program Rank of their matched residents: ", object$stats$avePRank)
  message("# Num Programs getting their top rank: ", object$stats$numPTop)
}