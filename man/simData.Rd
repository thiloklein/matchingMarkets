\name{simData}
\alias{simData}
\title{Simulate individual-level data for one-sided matching markets}
\usage{
simData(m, ind, seed = 123, singles = NULL, gpm = 2)
}
\arguments{
  \item{m}{integer indicating the number of markets to be
  simulated.}

  \item{ind}{integer (or vector) indicating the number of
  individuals per group.}

  \item{seed}{integer setting the state for random number
  generation. Defaults to \code{set.seed(123)}.}

  \item{singles}{integer giving the number of one-group
  markets.}

  \item{gpm}{integer giving the number of groups per
  market.}
}
\value{
\code{simData} returns a data frame with the randomly
generated variables mimicking those in dataset
\code{\link{baac00}}. \item{m.id}{categorical: market
identifier.} \item{g.id}{categorical: group identifier.}
\item{pi}{continuous: uniformly distributed project success
probability in [0,1].} \item{wst}{binary: indicator taking
the value 1 if last year was worse than the year before; 0
otherwise.} \item{occ1}{continuous: percentage of revenue
from income group 1.} \item{occ2}{continuous: percentage of
revenue from income group 2.} \item{occ3}{continuous:
percentage of revenue from income group 3.} \item{R}{NA:
group outcome is not simulated. It can be obtained using
the \code{design.matrix} function.}
}
\description{
Simulate individual-level data for one-sided matching
markets.
}
\examples{
## Coalitions
idata <- simData(m=4, ind=2:4, seed=124, singles=2, gpm=2)

## Rommmates
idata <- simData(m=3, ind=2, seed=124, gpm=3)
}
\author{
Thilo Klein
}
\keyword{generate}

