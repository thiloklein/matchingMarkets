\name{daa}
\alias{daa}
\title{Deferred Acceptance Algorithm}
\usage{
  daa(nStudents, nSlots, s.prefs = NULL, c.prefs = NULL)
}
\arguments{
  \item{nStudents}{an integer indicating the number of
  students in the market}

  \item{nSlots}{a vector of length nColleges that indicates
  the quota of each college}

  \item{s.prefs}{a nColleges x nStudents matrix with each
  column containing student i's ranking over colleges in
  decreasing order of preference (most preferred first)}

  \item{c.prefs}{a nStudents x nColleges matrix with each
  column containing college j's ranking over students in
  decreasing order of preference (most preferred first)}
}
\description{
  This function implements the student-proposing version of
  the Gale-Shapley algorithm. It takes student and college
  preference matrices and returns the unique stable
  assignment.
}

