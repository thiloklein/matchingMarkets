\name{design.matrix}
\alias{design.matrix}
\title{Design matrix for one-sided matching markets}
\usage{
design.matrix(x, m.id = "m.id", g.id = "g.id", R = "R",
  selection = NULL, outcome = NULL, roommates = FALSE,
  simulation = FALSE, assignment = "NTU", seed = 123, max.combs = Inf)
}
\arguments{
  \item{x}{data frame with individual-level characteristics
  of all group members including market- and
  group-identifiers.}

  \item{m.id}{character string giving the name of the
  market identifier variable. Defaults to \code{"m.id"}.}

  \item{g.id}{character string giving the name of the group
  identifier variable. Defaults to \code{"g.id"}.}

  \item{R}{dependent variable in outcome equation. Defaults
  to \code{"R"}.}

  \item{selection}{list containing variables and pertaining
  operators in the selection equation. The format is
  \code{operation = "variable"}. See the Details and
  Examples sections.}

  \item{outcome}{list containing variables and pertaining
  operators in the outcome equation. The format is
  \code{operation = "variable"}. See the Details and
  Examples sections.}

  \item{roommates}{logical: if \code{TRUE} data is assumed
  to come from a roomate game. This means that groups are
  of size two and the model matrix is prepared for
  individual-level analysis (peer-effects estimation). If
  \code{FALSE} (which is the default) data is assumed to
  come from a group/coalition formation game and the model
  matrix is prepared for group-level analysis.}

  \item{simulation}{logical: if \code{TRUE} the dependent
  variables of selection and outcome equations are
  simulated as follows: (i) model coefficients in both
  equations \code{alpha=beta=1}; (ii) covariance between
  error terms \code{delta=0.5}; (iii) error terms
  \code{eta} and \code{xi} are draws from a standard normal
  distribution.}

  \item{assignment}{character string giving the assignment
  procedure to be used if \code{simulation=TRUE}. Either
  \code{"NTU"} for non-transferable utility (which is the
  default value), \code{"TU"} for transferable utility or
  \code{"random"} for random assignment of individuals to
  groups.}

  \item{seed}{integer setting the state for random number
  generation if \code{simulation=TRUE}.}

  \item{max.combs}{integer (divisible by two) giving the
  maximum number of feasible groups to be used for
  generating group-level characteristics.}
}
\value{
\code{design.matrix} returns a list with the following
lists, whose elements each represent a matching market.
\item{D}{vector that indicates -- for all feasible groups
in the market -- whether a group is observed in the data
\code{D=1} or not \code{D=0}.} \item{V}{vector of group
valuations for all feasible groups in the market.}
\item{R}{vector of group outcomes for observed groups in
the market only.} \item{W}{data frame of group
characteristics for all feasible groups in the market.}
\item{X}{data frame of group characteristics for observed
groups in the market only.} \item{P}{vector that gives for
each group the index of the group comprised of residual
individuals in the market (for 2-group markets).}
\item{combs}{partitions matrix that gives all feasible
partitions of the market into groups of the observed
sizes.} \item{eta}{if \code{simulation=TRUE}, the standard
normally distributed errors in the selection equation.}
\item{xi}{if \code{simulation=TRUE}, the standard normally
distributed component of the errors in the selection
equation that is independent of \code{eta}.}
\item{epsilon}{if \code{simulation=TRUE}, the errors in the
outcome equation, given by \code{delta*eta + xi}.}
\item{E}{matrix that gives the indices of equilibrium group
members for each group in the market. Only differs from the
first two rows in \code{combs} if \code{simulation=TRUE}.}
}
\description{
Returns a model matrix with characteristics of \emph{all
feasible} groups of the same size as the observed groups in
the market.

For example, in the stable roommates problem with \eqn{n=4}
students \eqn{\{1,2,3,4\}}{{1,2,3,4}} sorting into groups
of 2, we have \eqn{{4 \choose 2}=6}{"4 choose 2" = 6}
feasible groups: (1,2)(3,4) (1,3)(2,4) (1,4)(2,3).

In the group formation problem with \eqn{n=6} students
\eqn{\{1,2,3,4,5,6\}}{{1,2,3,4,5,6}} sorting into groups of
3, we have \eqn{{6 \choose 3}=20}{"6 choose 3" = 20}
feasible groups. For the same students sorting into groups
of sizes 2 and 4, we have \eqn{{6 \choose 2} + {6 \choose
4}=30}{"6 choose 2" + "6 choose 4" = 30} feasible groups.

The input is individual-level data of all group members
from one-sided matching marktes; that is, from
group/coalition formation games.
}
\details{
Operators for variable transformations in \code{selection}
and \code{outcome} arguments. \describe{
\item{\code{add}}{sum over all group members and divide by
group size.} \item{\code{int}}{sum over all possible
two-way interactions \eqn{x*y} of group members and divide
by the number of those, given by \code{choose(n,2)}.}
\item{\code{ieq}}{sum over all possible two-way equality
assertions \eqn{1[x=y]} and divide by the number of those.}
\item{\code{ive}}{sum over all possible two-way
interactions of vectors of variables of group members and
divide by number of those.} \item{\code{dst}}{sum over all
possible two-way distances between players and divide by
number of those, where distance is defined as
\eqn{e^{-|x-y|}}{exp(-|x-y|)}.} \item{\code{sel}}{for
\code{roommates=TRUE} only: variable for individual (for
peer effects estimation).} \item{\code{oth}}{for
\code{roommates=TRUE} only: variable for other in the group
(for peer effects estimation).} }
}
\examples{
###############
## ROOMMATES ##
###############

## Simulate one-sided matching data for 3 markets (m=3) with 3 groups
## per market (gpm=3) and 2 individuals per group (ind=2)
idata <- simData(m=3, ind=2, gpm=3)

## Generate the design matrix
mdata <- design.matrix(x=idata,
  selection = list(add="pi", ieq="wst"),
  outcome = list(add="pi", ieq="wst"),
  assignment="TU", roommates=TRUE, simulation=TRUE)

################
## COALITIONS ##
################

## EITHER: Simulate one-sided matching data for 3 markets (m=3) with 2 groups
## per market (gpm=2) and 2 to 4 individuals per group (ind=2:4)
idata <- simData(m=3, ind=2:4, gpm=2)

## Generate the design matrix
mdata <- design.matrix(x=idata,
  selection = list(add="pi", ieq="wst"),
  outcome = list(add="pi", ieq="wst"),
  assignment="NTU", roommates=FALSE, simulation=TRUE)

\dontrun{
## OR: Load baac00 data from the Townsend Thai project
data(baac00)

## Generate the design matrix
mdata <- design.matrix(x=baac00,
  selection = list(add="pi", int="pi", ieq="wst", ive="occ"),
  outcome = list(add="pi", int="pi", ieq="wst", ive="occ",
    add=c("loan_size","loan_size2","lngroup_agei")),
  simulation=FALSE)}
}
\author{
Thilo Klein
}
\keyword{generate}

