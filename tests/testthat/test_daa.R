library(matchingMarkets)

test_that("check if daa matching is stable", {
  s.prefs <- matrix(c(1,2, 1,2, 1,2), 2,3)
  c.prefs <- matrix(c(1,2,3, 1,2,3), 3,2)
  expect_equal( daa(s.prefs=s.prefs, c.prefs=c.prefs)$singles, 3)
})

