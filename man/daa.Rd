\name{daa}
\alias{daa}
\title{Deferred Acceptance Algorithm}
\usage{
  daa(nStudents = ncol(s.prefs), nColleges = ncol(c.prefs),
    nSlots = rep(1, nColleges), s.prefs = NULL,
    c.prefs = NULL)
}
\arguments{
  \item{nStudents}{integer indicating the number of
  students (in the college admissions problem) or men (in
  the stable marriage problem) in the market. Defaults to
  \code{ncol(s.prefs)}.}

  \item{nColleges}{integer indicating the number of
  colleges (in the college admissions problem) or women (in
  the stable marriage problem) in the market. Defaults to
  \code{ncol(c.prefs)}.}

  \item{nSlots}{vector of length \code{nColleges}
  indicating the number of places (i.e. quota) of each
  college. Defaults to \code{rep(1,nColleges)} for the
  marriage problem.}

  \item{s.prefs}{matrix of dimension \code{nColleges} x
  \code{nStudents} with the \code{i}th column containing
  student \code{i}'s ranking over colleges in decreasing
  order of preference (i.e. most preferred first).}

  \item{c.prefs}{matrix of dimension \code{nStudents} x
  \code{nColleges} with the \code{j}th column containing
  college \code{j}'s ranking over students in decreasing
  order of preference (i.e. most preferred first).}
}
\value{
  'daa' returns a list with the following items.
  \item{s.prefs}{students' preference matrix.}
  \item{c.prefs}{colleges' preference matrix.}
  \item{iterations}{number of interations required to find
  the stable matching.} \item{matches}{identifier of
  students (men) assigned to colleges (women).}
  \item{match.mat}{matching matrix of dimension
  \code{nStudents x nColleges}.} \item{singles}{identifier
  of single/unmatched students (men).}
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
  'daa' requires the following combination of arguments,
  subject to the matching problem. \describe{
  \item{\code{nStudents, nColleges}}{Marriage problem with
  random preferences.} \item{\code{s.prefs,
  c.prefs}}{Marriage problem with given preferences.}
  \item{\code{nStudents, nSlots}}{College admissions
  problem with random preferences.} \item{\code{s.prefs,
  c.prefs, nSlots}}{College admissions problem with given
  preferences.} }
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
  Thilo Klein \email{thilo@klein.co.uk}
}
\references{
  Gale, D. and Shapley, L.S. (1962). College admissions and
  the stability of marriage. The American Mathematical
  Monthly, 69(1):9--15.
}

