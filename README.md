# matchingMarkets

> An R package for the analysis of data from stable matchings.

# Functions

`matchingMarkets` currently comes with three algorithms:

* `daa`: Gale-Shapley Deferred Acceptance Algorithm. Finds stable matchings in two-sided matching markets. Implemented for both the Stable Marriage Problem (one-to-one matching) and the College Admissions Problem (many-to-one matching).

* `plp`: Partitioning Linear Programme. Finds stable matchings in the Roommates Problem (one-sided matching market) with transferable utility.

* `ttc`: Top-Trading-Cycles Algorithm. Finds stable matchings in the Housing Market Probem.

# Installation

To get the current development version from github:

```R
install.packages("devtools")
devtools::install_github(repo="matchingMarkets", username="thiloklein")
library(matchingMarkets)
```

