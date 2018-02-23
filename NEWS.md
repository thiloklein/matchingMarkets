# Changelog of package `matchingMarkets`

Please note that only the most significant changes are reported here.
A full ChangeLog is available in the log messages of the SVN repository
on [R-Forge](https://r-forge.r-project.org/scm/viewvc.php/?root=matchingmarkets)
and on [GitHub](https://github.com/thiloklein/matchingMarkets).


## matchingMarkets 0.3-5 (2018-02-23)

This is a minor update

* Added R wrapper for Roth-Peranson Algorithm in function `hri2`.

## matchingMarkets 0.3-3 (2017-03-27)

This is a minor update

* Implemented multi-core parrallel processing for estimators in function `stabit2`, which can be specified using the `nCores` argument.
* Updated immediate acceptance algorithm `iaa` and top-trading-cycles `ttc` functions. Thanks to [Sándor Sóvágó](http://www.tinbergen.nl/phd-student/sandor-s-sovago-2/) at Tinbergen Institute and [Kevin Breuer](http://economicdesign.uni-koeln.de/en/home/researchers/kevin-breuer/) at University of Cologne for the reports.

## matchingMarkets 0.3-1 (2016-08-21)

This is a major update

* Replaced stable matching algorithms with constraint programming model implemented for hospital/residents problem `hri` and stable roommates `sri` with incomplete lists.
* Added `plot` and `summary` methods for estimators.
* Allowed for thinning in `stabit2` function.

## matchingMarkets 0.2-2 (2016-02-07)

This is a minor update

* Allowed for two selection equations in `stabit2` function for two-sided matching markets. 


## matchingMarkets 0.2-1 (2016-01-30)

This is a major update

* Added `stabit2` function for two-sided matching markets. 


## matchingMarkets 0.1-6 (2015-10-05)

This is a major update

* Added [package vignette](https://CRAN.R-project.org/package=matchingMarkets/vignettes/matching.pdf). 


## matchingMarkets 0.1-5 (2015-07-08)

This is a minor update.

* Fixed `daa` function for college admissions problems when number of students exceeds number of colleges. Thanks to [Jan Tilly](http://jtilly.io/) at University of Pennsylvania for the report.

## matchingMarkets 0.1-1 (2015-07-08)

Initial commit.

