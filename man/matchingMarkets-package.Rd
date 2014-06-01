\name{matchingMarkets-package}
\alias{matchingMarkets}
\title{An R package for the analysis of stable matchings.}
\docType{package}




\description{
  The \code{matchingMarkets} package contains \R and Matlab code for the estimation 
  of structural models that correct for the sample selection bias of observed outcomes in 
  matching markets.

  Matching is concerned with who transacts with whom, and how. For example, who 
  works at which job, which students go to which school, who forms a workgroup 
  with whom, and so on.

  The empirical analysis of matching markets is naturally subject to sample 
  selection problems. If agents match assortatively on characteristics 
  unobserved to the analyst but correlated with both the exogenous variable and 
  the outcome of interest, regression estimates will generally be biased.

  The \code{matchingMarkets} package comes with
  \enumerate{
  \item \emph{Estimators}. The estimators implemented in function \code{\link{smm}} 
  correct for the selection bias from endogenous matching. In addition, \code{\link{mfx}} 
  produces marginal effects from coefficients in binary outcome and selection equations 
  and \code{\link{khb}} implements a statistical test of confounding due to sample selection.

  This version of the package provides solutions for two commonly observed matching 
  processes: 
  (i) the \emph{group formation problem} with fixed group sizes and 
  (ii) the \emph{roommates problem} with transferable utility. 
  Work on the college admissions problem is still in progress.
  These processes determine which matches are observed -- 
  and which are not -- and this is a sample selection problem.

  \item \emph{Design matrices}. The estimators require independent variables for all feasible, 
  i.e., observed and counterfactual, matches in the market. Generating the characteristics
  of all counterfactual matches is a combinatorial problem. The function \code{\link{design.matrix}} 
  generates a design matrix based on data for observed matches and pre-specified 
  transformations to generate counterfactual matches. 

  \item \emph{Algorithms}. Furthermore, the package contains several matching 
  algorithms such as the deferred acceptance algorithm (\code{\link{daa}}) for \href{http://en.wikipedia.org/wiki/Stable_matching}{stable marriages} and \href{http://en.wikipedia.org/wiki/Hospital_resident}{college admissions}, the top-trading-cycles algorithm (\code{\link{ttc}}) for \href{http://en.wikipedia.org/wiki/Herbert_Scarf\#8._The_Housing_Market}{house allocation} and a partioning linear program (\code{\link{plp}}) for the \href{http://en.wikipedia.org/wiki/Stable_roommates_problem}{stable roommates problem}. These can be used to obtain stable matchings from simulated or 
  real preference data. 

  \item \emph{Data}. In addition to the \code{\link{baac00}} dataset from borrowing groups in Thailands largest agricultural lending program, the package provides functions to simulate one's own data from matching markets. \code{\link{simData}} generates individual-level data and \code{\link{design.matrix}} has an argument \code{simulation} which, if set to \code{TRUE}, generates group-level data and determines which groups are observed in equilibrium based on an underlying linear stochastic model.
  }
}




\author{
  Thilo Klein
}




\references{
  Bresnahan, T. and Reiss, P. (1991). Empirical models of discrete games. \emph{Journal of Econometrics}, 48(1-2):57--81.

  Gale, D. and Shapley, L.S. (1962). College admissions and the stability of marriage. \emph{The American Mathematical Monthly}, 69(1):9--15.

  Heckman, J. (1979). Sample selection bias as a specification error. \emph{Econometrica}, 47(1):153--161.

  Pycia, M. (2012). Stability and preference alignment in matching and coalition formation. \emph{Econometrica}, 80(1):323--362.

  Sorensen, M. (2007). How smart is smart money? A two-sided matching model of venture capital. \emph{The Journal of Finance}, 62(6):2725--2762.
}

\section{System requirements}{
  This package is implemented in plain \R. To speed up the iterations of the Gibbs sampling algorithm, the estimator in function \code{\link{smm}} uses an interface to either \href{http://www.mathworks.co.uk/products/matlab/}{Matlab} or the open source Matlab clone \href{http://www.gnu.org/software/octave/}{Octave} which requires the \R packages \code{\link[R.matlab:R.matlab-package]{R.matlab}} or \code{\link[RcppOctave:RcppOctave]{RcppOctave}}, respectively.

  To connect to Matlab, Matlab v6 or higher is required. It does not work with Matlab v5 or before (because those versions do not support Java).

  For Octave on Linux machines, the only requirement is to have Octave >= 3.2.4 and its development files installed. To install Octave (>= 3.6) under Linux, run the following at the command line.
  \preformatted{
  # install octave and development files
  sudo apt-get install octave liboctave-dev
  # open octave 
  sudo octave
  # install packages 'io' and 'statistics' from the octave prompt
  pkg install -forge io
  pkg install -forge statistics}

  For Octave under Windows and Mac OS X, the \code{\link[RcppOctave:RcppOctave]{RcppOctave}} package supports Octave versions \href{http://wiki.octave.org/Octave_for_Windows}{\emph{minGW}} and \href{http://wiki.octave.org/Octave_for_MacOS_X\#Homebrew}{\emph{Homebrew}}, respectively. See the \href{http://cran.r-project.org/web/packages/RcppOctave/vignettes/RcppOctave.pdf}{vignette} for more details.
}




\section{Frequently Asked Questions}{
  \itemize{
  \item \emph{Why can I not use the classic Heckman correction?}

        Estimators such as the Heckman (1979) correction (in package \href{http://cran.r-project.org/web/packages/sampleSelection/index.html}{\code{sampleSelection}}) or double selection models are 
        inappropriate for this class of selection problems. To see this, note that 
        a simple first stage discrete choice model assumes that an observed match 
        reveals match partners' preferences over each other. In a matching market, 
        however, agents can only choose from the set of partners who would be 
        willing to form a match with them and we do not observe the players' 
        relevant choice sets. 

  \item \emph{Do I need an instrumental variable to estimate the model?}

        Short answer: No. Long answer: The characteristics of other agents in the 
        market serve as the source of exogenous variation necessary to identify 
        the model. The identifying exclusion restriction is that characteristics 
        of all agents in the market affect the matching, i.e., who matches with 
        whom, but it is only the characteristics of the match partners that affect 
        the outcome of a particular match once it is formed. No additional instruments 
        are required for identification (Sorensen, 2007). 

  \item \emph{What are the main assumptions underlying the estimator?}

        The approach has certain limitations rooted in its restrictive economic assumptions. 
        \enumerate{
        \item The matching models are \emph{complete information} models. That is, agents are assumed to have a complete knowledge of the qualities of other market participants.

        \item The models are \emph{static equilibrium} models. This implies that (i) the observed matching must be an equilibrium, i.e., no two agents would prefer to leave their current partners in order to form a new match (definition of pairwise stability), and (ii) the equilibrium must be unique for the likelihood function of the model to be well defined (Bresnahan and Reiss, 1991).

        \item Uniqueness results can be obtained in two ways. First, as is common in the industrial organization literature, by imposing suitable \emph{preference restrictions}. A suitable restriction on agents' preferences that guarantees a unique equilibrium is alignment (Pycia, 2012). In a group formation model, (pairwise) preference alignment states that any two agents who belong to the same groups must prefer the same group over the other. Second, by choosing a \emph{market assigment} based on matching algorithms that produce a unique stable matching, such as the well-studied Gale and Shapley (1962) deferred acceptance algorithm.

       \item Finally, the models assume \emph{bivariate normality} of the errors in selection and outcome equation. If that assumption fails, the estimator is generally inconsistent and can provide misleading inference in small samples.
       }
  }
}




\section{How to cite this package}{
  Whenever using this package, please cite as

  Klein, T. (2014). matchingMarkets: An R package for the analysis of stable matchings. R package version 0.1.
}




\seealso{
  \href{http://cran.r-project.org/web/packages/sampleSelection/index.html}{sampleSelection}
}



\section{Work in progress}{
  \enumerate{
  \item Extension of the estimator in \code{\link{smm}} and the \code{\link{design.matrix}} function to allow for \eqn{n \ge 2} groups per market for the group formation problem. Currently, \code{\link{design.matrix}} supports \eqn{n \ge 2} groups per market for the roommates problem only.
  \item Extension of the estimator in the \code{\link{smm}} function to allow for two-sided matching data from the college admissions problem. Currently, only one-sided matching data from the group formation and stable roommates problem is supported.
  \item Rewrite of the code from Matlab to C/C++ using the GSL libraries.
%  \item S3 methods for functions \code{smm} and \code{mfx}. At least \code{print} and \code{summary}.
%  \item In \code{smm} function the arguments "simulation" and "NTU" are redundant. Replace with \code{simulation} options \code{"TU"; "NTU"; "random"; "none"}.
%  \item Add more consistency checks.
  }
}

\keyword{package}


