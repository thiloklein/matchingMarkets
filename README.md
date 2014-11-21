# matchingMarkets

> An R package for the analysis of stable matchings.

# Functions

`matchingMarkets` currently comes with one estimator:

* `stabit`: Implements a Bayes estimator that corrects for sample selection in matching markets when the selection process is a one-sided matching game.

and three algorithms that can be used to simulate matching data:

* `daa`: Gale-Shapley Deferred Acceptance Algorithm. Finds stable matchings in two-sided matching markets. Implemented for both the Stable Marriage Problem (one-to-one matching) and the College Admissions Problem (many-to-one matching).

* `plp`: Partitioning Linear Programme. Finds stable matchings in the Roommates Problem (one-sided matching market) with transferable utility.

* `ttc`: Top-Trading-Cycles Algorithm. Finds stable matchings in the Housing Market Probem.

# Installation

To get the current development version from GitHub:

```R
install.packages("devtools")
devtools::install_github("thiloklein/matchingMarkets")
library(matchingMarkets)
```
or from R-Forge:

```R
install.packages("matchingMarkets", repos="http://R-Forge.R-project.org")
library(matchingMarkets)
```

To browse the package documentation:

```R
?matchingMarkets
```
