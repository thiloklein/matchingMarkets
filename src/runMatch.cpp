#include <Rcpp.h>
#include <problem.h>
#include <rpmatcher.h>
#include <params.h>

using namespace Rcpp;

static DAmatcher* dam{};

// [[Rcpp::export]]
List runMatch(Rcpp::List students, Rcpp::List programs, Rcpp::List couples) {
  //matchrp using Roth Peranson 1999 algorithm with re-randomization of couple ordering
  params.readOptions();

  Problem prob {};
  if(!prob.readProblem(students, couples, programs)) {
    Rcpp::Rcout  << prob.getError();
  }
  
  //cpu limit -> 10 sec
  rlimit rl;
  rl.rlim_cur = 10;
  
  dam = new RPmatcher {};
  auto match = dam->match(prob);
  
  return List::create(Named("matchings") = prob.returnMatch(match),
                      //instead of long description print short stats
                      //Named("summary") = prob.summaryMatch(match),
                      Named("stats") = prob.getStats(match));
}
