\name{toMatLab}
\alias{toMatLab}
\title{Write one-sided matching data to Matlab .mat file}
\usage{
toMatLab(x, folder = getwd(), id = "")
}
\arguments{
  \item{x}{one-sided matching data after applying the
  \code{design.matrix()} function.}

  \item{folder}{path that .mat file will be saved to.
  Defaults to current working directory \code{getwd()}.}

  \item{id}{number or character to be added to file name.
  The default is none.}
}
\description{
Write an object from \code{design.matrix} to MAT v5 format
using the \code{\link[R.matlab:R.matlab-package]{R.matlab}}
package
}
\examples{
\dontrun{
## Load one-sided matching data
data(mdata, package="matchingMarkets")

## Write data to file mdata123.mat on Desktop
toMatLab(x=mdata, folder="~/Desktop/", id=123)}
}
\author{
Thilo Klein
}
\keyword{generate}

