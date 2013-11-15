\name{daa}
\alias{daa}
\title{Deferred Acceptance Algorithm}
\usage{
  daa(nStudents = dim(s.prefs)[2],
    nColleges = dim(c.prefs)[2],
    nSlots = rep(1, nColleges), s.prefs = NULL,
    c.prefs = NULL)
}
\arguments{
  \item{nStudents}{integer indicating the number of
  students (in the college admissions problem) or men (in
  the stable marriage problem) in the market}

  \item{nColleges}{integer indicating the number of
  colleges (in the college admissions problem) or women (in
  the stable marriage problem) in the market}

  \item{nSlots}{vector of length \code{nColleges}
  indicating the number of places (i.e. quota) of each
  college. Defaults to \code{rep(1,nColleges)} for the
  marriage problem.}

  \item{s.prefs}{matrix of dimension \code{nColleges x
  nStudents} with the \code{i}th column containing student
  \code{i}'s ranking over colleges in decreasing order of
  preference (i.e. most preferred first)}

  \item{c.prefs}{matrix of dimension \code{nStudents x
  nColleges} with the \code{j}th column containing college
  \code{j}'s rankimatg over students in decreasing order of
  preference (i.e. most preferred first)}
}
\description{
  Finds the student (men) optimal matching in the
  \href{http://en.wikipedia.org/wiki/Hospital_resident}{college
  admissions}
  (\href{http://en.wikipedia.org/wiki/Stable_matching}{stable
  marriage}) problem. Uses the Gale-Shapley (1962) Deferred
  Acceptance Algorithm with student (male) offer based on
  given or randomly generated preferences.
}
\section{Minimum required arguments}{
  Subject to the matching problem, 'daa' requires the
  following arguments. Marriage problem with random
  preferences: \code{nStudents, nColleges}. Marriage
  problem with given preferences: \code{s.prefs, c.prefs}.
  College admissions problem with random preferences:
  \code{nStudents, nSlots}. College admissions problem with
  given preferences: \code{s.prefs, c.prefs, nSlots}.
}

\section{Value}{
  'daa' returns a list with the following items. \describe{
  \item{\code{s.prefs}}{students' preference matrix.}
  \item{\code{c.prefs}}{colleges' preference matrix.}
  \item{\code{iterations}}{number of interations required
  to find the stable matching.}
  \item{\code{matches}}{identifier of students (men)
  assigned to colleges (women).}
  \item{\code{match.mat}}{matching matrix of dimension
  \code{nStudents x nColleges}.}
  \item{\code{singles}}{identifier of single/unmatched
  students (men).} }
}
\examples{
## Marriage problem (3 men, 2 women) with random preferences:
daa(nStudents=3, nColleges=2)

## Marriage problem (3 men, 2 women) with given preferences:
s.prefs <- matrix(c(1,2, 1,2, 1,2), 2,3)
c.prefs <- matrix(c(1,2,3, 1,2,3), 3,2)
daa(s.prefs=s.prefs, c.prefs=c.prefs)

## College admission problem (7 students, 2 colleges
## with 3 slots each) with random preferences:
daa(nStudents=7, nSlots=c(3,3))

## College admission problem (7 students, 2 colleges
## with 3 slots each) with given preferences:
s.prefs <- matrix(c(1,2, 1,2, 1,2, 1,2, 1,2, 1,2, 1,2), 2,7)
c.prefs <- matrix(c(1,2,3,4,5,6,7, 1,2,3,4,5,6,7), 7,2)
daa(s.prefs=s.prefs, c.prefs=c.prefs, nSlots=c(3,3))
}
\author{
  Thilo Klein <\email{thilo@klein.co.uk}>
}
\references{
  Gale, D. and Shapley, L.S. (1962). College admissions and
  the stability of marriage. The American Mathematical
  Monthly, 69(1):9--15.
}

