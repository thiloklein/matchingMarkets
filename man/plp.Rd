\name{plp}
\alias{plp}
\title{Partitioning Linear Programme}
\usage{
  plp(V = NULL, N = NULL)
}
\arguments{
  \item{N}{number of players in the market}

  \item{V}{valuation matrix of dimension \code{NxN} that
  gives row-players valuation over column players (or vice
  versa)}
}
\value{
  'plp' returns a list with the following items.
  \item{Assignment.matrix}{upper triangular matrix of
  dimension \code{NxN} with entries of 1 for equilibrium
  pairs and 0 otherwise.} \item{Equilibrium.groups}{matrix
  that gives the \code{N/2} equilibrium pairs and
  equilibrium partners' mutual valuations.}
}
\description{
  Finds the stable matching in the
  \href{http://matchingmarkets.org/sites/en.wikipedia.org/wiki/Stable_roommates_problem}{stable
  roommates problem} with transferable utility. Uses the
  Partitioning Linear Programme formulated in Quint (1991).
}
\examples{
## Roommate problem with 10 players, transferable utility and random preferences:
plp(N=10)

## Roommate problem with 10 players, transferable utility and given preferences:
V <- matrix(rep(1:10, 10), 10, 10)
plp(V=V)
}
\author{
  Thilo Klein \email{thilo@klein.co.uk}
}
\references{
  Quint, T. (1991). Necessary and sufficient conditions for
  balancedness in partitioning games. Mathematical Social
  Sciences, 22(1):87--91.
}

