# matchingMarkets
> Analysis of Stable Matchings in R

[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/matchingMarkets)](http://cran.r-project.org/package=matchingMarkets)
[![CRAN_Downloads](http://cranlogs.r-pkg.org/badges/grand-total/matchingMarkets?color=brightgreen)](http://cran.r-project.org/package=matchingMarkets)


## Functions

`matchingMarkets` comes with two estimators:

* `stabit`: Implements a Bayes estimator that corrects for sample selection in matching markets when the selection process is a one-sided matching game (i.e. group formation).

* `stabit2`: Implements the Bayes estimator for a two-sided matching game (i.e. the [college admissions](http://en.wikipedia.org/wiki/Stable_marriage_problem#Similar_problems) and [stable marriage](http://en.wikipedia.org/wiki/Stable_marriage_problem) problems).

and three algorithms that can be used to simulate matching data:

* `daa`: Gale-Shapley Deferred Acceptance Algorithm. Finds stable matchings in two-sided matching markets. Implemented for both the [stable marriage problem](http://en.wikipedia.org/wiki/Stable_marriage_problem) (one-to-one matching) and the [college admissions problem](http://en.wikipedia.org/wiki/Stable_marriage_problem#Similar_problems) (many-to-one matching).

* `plp`: Partitioning Linear Programme. Finds stable matchings in the [roommates problem](https://en.wikipedia.org/wiki/Stable_roommates_problem) (one-sided matching market) with transferable utility.

* `ttc`: Top-Trading-Cycles Algorithm. Finds stable matchings in the [housing market problem](http://en.wikipedia.org/wiki/Herbert_Scarf#8._The_Housing_Market).


## Installation

To get the latest *stable version* from [CRAN](http://cran.r-project.org/package=matchingMarkets):

```R
install.packages("matchingMarkets")
library(matchingMarkets)
```

Under Linux, the dependency package `gmp` requires that you have GNU MP (> 4.1.4) installed, see http://gmplib.org.

To get the most recent *development version* from [GitHub](https://github.com/thiloklein/matchingMarkets):

```R
install.packages("devtools")
devtools::install_github("thiloklein/matchingMarkets")
library(matchingMarkets)
```
or from [R-Forge](https://r-forge.r-project.org/R/?group_id=1906):

```R
install.packages("matchingMarkets", repos="http://R-Forge.R-project.org")
library(matchingMarkets)
```


## Documentation

Package [documentation](http://cran.r-project.org/web/packages/matchingMarkets/matchingMarkets.pdf) and [vignette](https://cran.r-project.org/web/packages/matchingMarkets/vignettes/matching.pdf) are available from the [CRAN page](http://cran.r-project.org/package=matchingMarkets). An application of the estimator in function `stabit` is in [Klein (2015)](https://ideas.repec.org/p/cam/camdae/1521.html).

