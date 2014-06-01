\name{mdata}
\docType{data}
\alias{mdata}

\title{Simulated one-sided matching data}

\description{
   One-sided matching data for 200 markets with 2 groups of size 5 each. The \code{mdata} data frame was generated using \code{simData()} and \code{design.matrix(...,simulation=TRUE)} which returns a list of lists, whose elements each represent a matching market.
}

\usage{data(mdata)}

\format{
   This list contains the following elements:
   \describe{
\item{D}{vector that indicates -- for all feasible groups in the market -- whether a group is observed in the data \code{D=1} or not \code{D=0}.}
\item{V}{vector of group valuations for all feasible groups in the market.}
\item{R}{vector of group outcomes for observed groups in the market only.}
\item{W}{data frame of group characteristics for all feasible groups in the market.}
\item{X}{data frame of group characteristics for observed groups in the market only.}
\item{P}{vector that gives for each group the index of the group comprised of residual individuals in the market (for 2-group markets).}
\item{combs}{partitions matrix that gives all feasible partitions of the market into groups of the observed sizes.}
\item{eta}{standard normally distributed errors in the selection equation.}
\item{xi}{tandard normally distributed component of the errors in the selection equation that is independent of \code{eta}.}
\item{epsilon}{errors in the outcome equation, given by \code{delta*eta + xi}.}
\item{E}{matrix that gives the indices of equilibrium group members for each group in the market.}
   }
}

\keyword{datasets}

