\name{ttc}
\alias{ttc}
\title{Top-Trading-Cycles Algorithm}
\usage{
  ttc(P = NULL, X = NULL)
}
\arguments{
  \item{P}{list of individuals' preference rankings over
  objects}

  \item{X}{2-column-matrix of objects ('obj') and their
  owners ('ind')}
}
\description{
  Finds the stable matching in the house allocation problem
  with existing tenants. Uses the Top-Trading-Cycles
  Algorithm proposed in Abdulkadiroglu and Sonmez (1999).
}
\section{Value}{
  'ttc' returns a stable matching solution for the housing
  market problem based on the Top-Trading-Cycles algorithm.
}
\examples{
## generate list of individuals' preference rankings over objects
P <- list()
P[[1]] <- c(2,5,1,4,3)    # individual 1
P[[2]] <- c(1,5,4,3,2)    # individual 2
P[[3]] <- c(2,1,4,3,5)    # individual 3
P[[4]] <- c(2,4,3,1,5)    # individual 4
P[[5]] <- c(4,3,1,2,5); P # individual 5

## generate 2-column-matrix of objects ('obj') and their owners ('ind')
X <- data.frame(ind=1:5, obj=1:5); X

## find assignment based on TTC algorithm
ttc(P=P,X=X)
}
\author{
  Thilo Klein <\email{thilo@klein.co.uk}>
}
\references{
  Abdulkadiroglu, A. and Sonmez, T. (1999). House
  Allocation with Existing Tenants. Journal of Economic
  Theory, 88(2):233--260.
}

