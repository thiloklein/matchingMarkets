\name{mfx}
\alias{mfx}
\title{Marginal effects for probit and matching models}
\usage{
mfx(m, toLatex = FALSE)
}
\arguments{
  \item{m}{EITHER: an object returned from function
  \code{smm}; OR: the file path to a \code{.mat} file after
  running \code{smm.m} in Matlab.}

  \item{toLatex}{logical: if \code{TRUE} the result tables
  are printed in Latex format. The default setting is
  \code{FALSE}.}
}
\description{
Marginal effects from regression coefficients for probit
and matching models.
}
\examples{
## Read model results from Klein (2014), Table 5
filepath <- system.file("scripts/TK_gibbsiter_ntu.mat", package="matchingMarkets")

## Apply mfx function and print results
M <- mfx(m=filepath)
M$mfx.selection; M$mfx.outcome
}
\author{
Thilo Klein
}
\keyword{summary}

