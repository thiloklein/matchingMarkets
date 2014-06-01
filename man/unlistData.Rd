\name{unlistData}
\alias{unlistData}
\title{Unlist one-sided matching data for simple regression analysis}
\usage{
unlistData(x, roommates = FALSE)
}
\arguments{
  \item{x}{matching data after applying the
  \code{design.matrix} function.}

  \item{roommates}{logical: if \code{TRUE} data is assumed
  to come from a roomate game. This means that groups are
  of size two and the model matrix is prepared for
  individual-level analysis (peer-effects estimation). If
  \code{FALSE} (which is the default) data is assumed to
  come from a group/coalition formation game and the model
  matrix is prepared for group-level analysis.}
}
\value{
\code{unlistData} returns a list with the following
elements. \item{SEL}{data frame comprising variables in
selection equation and number of observations equal to the
number of feasible groups.} \item{OUT}{data frame
comprising variables in outcome equation and number of
observations equal to the number of equilibrium groups.}
}
\description{
Unlist an object from \code{design.matrix} for further
analysis.
}
\examples{
## Load one-sided matching data
data(mdata, package="matchingMarkets")

## Unlist the data
mdata <- unlistData(x=mdata)
str(mdata$SEL); str(mdata$OUT)
}
\author{
Thilo Klein
}
\keyword{generate}

