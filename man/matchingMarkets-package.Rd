\name{matchingMarkets-package}
\alias{matchingMarkets}
\title{An R package for the analysis of stable matchings.}
\docType{package}




\description{
  The \code{matchingMarkets} package contains \R and \code{C++} code for the estimation 
  of structural models that correct for the sample selection bias of observed outcomes in 
  matching markets.

  Matching is concerned with who transacts with whom, and how. For example, who 
  works at which job, which students go to which school, who forms a workgroup 
  with whom, and so on.

  The empirical analysis of matching markets is naturally subject to sample 
  selection problems. If agents match assortatively on characteristics 
  unobserved to the analyst but correlated with both the exogenous variable and 
  the outcome of interest, regression estimates will generally be biased.

  The \code{matchingMarkets} package comprises
  \enumerate{
  \item \emph{Bayes estimators}. The estimators implemented in function \code{\link{stabit}} 
  and \code{\link{stabit2}} correct for the selection bias from endogenous matching. 

  The current package version provides solutions for two commonly observed matching 
  processes: 
  (i) the \emph{group formation problem} with fixed group sizes and
  (ii) the \emph{college admissions problem}.
  These processes determine which matches are observed -- 
  and which are not -- and this is a sample selection problem.
  
  \item \emph{Post-estimation tools}. Function \code{\link{mfx}} 
  computes marginal effects from coefficients in binary outcome and selection equations 
  and \code{\link{khb}} implements the Karlson-Holm-Breen test for confounding due to sample selection.

  \item \emph{Design matrix generation}. The estimators are based on independent variables for all feasible, 
  i.e., observed and counterfactual, matches in the market. Generating the characteristics
  of all feasible matches from individual-level data is a combinatorial problem. 
  The \code{\link{stabit}} function has an argument \code{method="model.frame"} 
  that returns a design matrix based on pre-specified transformations to generate counterfactual matches. 

  \item \emph{Algorithms}. The package also contains several matching 
  algorithms such as the deferred acceptance algorithm (\code{\link{daa}}) for \href{http://en.wikipedia.org/wiki/Stable_matching}{stable marriages} and \href{http://en.wikipedia.org/wiki/Hospital_resident}{college admissions}, the top-trading-cycles algorithm (\code{\link{ttc}}) for \href{http://en.wikipedia.org/wiki/Herbert_Scarf#8._The_Housing_Market}{house allocation} and a partioning linear program (\code{\link{plp}}) for the \href{http://en.wikipedia.org/wiki/Stable_roommates_problem}{stable roommates problem}. These can be used to obtain stable matchings from simulated or 
  real preference data. 

  \item \emph{Data}. In addition to the \code{\link{baac00}} dataset from borrowing groups in Thailands largest agricultural lending program, the package provides functions to simulate one's own data from matching markets. \code{\link{stabsim}} generates individual-level data and the \code{\link{stabit}} function has an argument \code{simulation} which generates group-level data and determines which groups are observed in equilibrium based on an underlying linear stochastic model.
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




\section{Frequently Asked Questions}{
  \itemize{
  \item \emph{Why can I not use the classic Heckman correction?}

        Estimators such as the Heckman (1979) correction (in package \href{http://cran.r-project.org/package=sampleSelection}{\code{sampleSelection}}) or double selection models are 
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

  Klein, T. (2015). matchingMarkets: Structural Estimator and Algorithms for the Analysis of Stable Matchings. R package version 0.1-8.
}




\seealso{
  \href{http://cran.r-project.org/package=sampleSelection}{sampleSelection}
}


\keyword{package}


