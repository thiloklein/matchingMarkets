/***********[rpmatcher.cc]

-------

Implementation of "rpmatcher" class. 
This class implements the rp matching algorithm. It alters updates the matching stored in
a problem object. 

------

Copyright (c) 2014, Fahiem Bacchus

Permission is hereby granted, free of charge, to any person obtaining a
copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be included
in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

***********/

#include <assert.h>
#include <iostream>
#include <vector>
#include <utility>
#include <algorithm>
//SVEN
#include <random>
#include <algorithm>
#include <iterator>

#include <params.h>
#include <io.h>
#include <rpmatcher.h>
#include <System.h>

using std::ostream;
using std::vector;
using std::pair;
//using std::cout;
//using std::cerr;
//using std::random_shuffle;
//Sven
using std::shuffle;

bool RPmatcher::match_(Problem& prob) {
  vector<Rid> rMatched {};
  rpCplsMatched = 0;
/*  for(auto& r : prob.Res()) {
    if(!r.inCouple()) {
      rMatched.push_back(r);
      if(!extendMatch(r, rMatched)) {
	return false;
      }
    }
    }*/

  //So simpler initial matching of singles.
  //Use rmatched as a stack of single residents to match.
  for(auto& r: prob.Res()) 
    if(!r.inCouple()) 
      rMatched.push_back(r);

  while(!rMatched.empty()) {
    auto r = rMatched.back();
    rMatched.pop_back();
    processResident(r, rMatched);
  }
  //Now reset rMatched to be set of matched singles
  for(auto& r: prob.Res()) 
    if(!r.inCouple()) 
      rMatched.push_back(r);


  if(params.verbosity > 0) {
    Rcpp::Rcout  << "#Matched singles\n";
    if(params.verbosity > 2) {
      Rcpp::Rcout  << "#Current match:\n";
      prob.printMatch(true, false);
    }
  }

  vector<int> randOrder {};
  for(size_t i = 0; i < prob.Cpl().size(); ++i)
    randOrder.push_back(i);

  while(1) {
    size_t i = 0;
    for(  ; i < prob.Cpl().size(); ++i) {
      auto c = prob.Cpl()[randOrder[i]];
      rMatched.push_back(c.r1());
      if(!extendMatch(c.r1(), rMatched)) 
	break;
      ++rpCplsMatched;
      if(params.verbosity > 1) {
	Rcpp::Rcout  << "#Matched couple: " << c.id() << "\n";
      }
    }
    if(i == prob.Cpl().size())
      return true;
    if (params.rnd) {
      if(params.verbosity > 1) {
	Rcpp::Rcout  << "#Reordering Couples\n";
      }
      ++rndReorderings;
      for(auto cpl : prob.Cpl()) {
	unmatch(cpl.r1());
	unmatch(cpl.r2());
      }
      //std::random_shuffle(randOrder.begin(), randOrder.end());
      //SVEN, http://en.cppreference.com/w/cpp/algorithm/random_shuffle
      std::random_device rd;
      std::mt19937 g(rd());
      std::shuffle(randOrder.begin(), randOrder.end(), g);
    }
    else {
      Rcpp::Rcout  << "Failed at Round " << i << "\n";
      return false;
    }
  }
}

bool RPmatcher::extendMatch(Rid newRes, vector<Rid>& matched) {
  //Invariant:
  // r on rToProcess stack:
  //  if r is a member of a couple then r == r1() the first
  //     resident of the couple. 
  //  There are no duplicates on the rToProcess stack
  //     Nor on the pToProcess stack.

  for(auto rid : matched)
    rid->clearMatchCounts();
  Resident::maxMatchCount = 0;

  if(params.verbosity > 2) {
    Rcpp::Rcout  << "#LOG: extending match to " << newRes;
    if(newRes->inCouple()) 
      Rcpp::Rcout  << " couple = " << newRes->couple();
    Rcpp::Rcout  << "\n";
  }

  vector<Rid> rToProcess {newRes}; //residents to process
  while(!rToProcess.empty() || !chkMatch(rToProcess, matched)) {
    if(params.cclim >= 0 && Resident::maxMatchCount > params.cclim)
      return false;
    if(Resident::maxMatchCount > maxRepeatedMatches)
      maxRepeatedMatches = Resident::maxMatchCount;

    if(!rToProcess.empty()) {
	auto r = rToProcess.back();
	rToProcess.pop_back();
	processResident(r, rToProcess);
    }
    else {
      Rcpp::Rcerr << "ERROR: extendMatch looping with empty rToProcess stack\n";
    }
  }
  return true;
}

template<typename T> 
//Add to r to vector only if r is not already in vector.
void addUnique(T r, vector<T>& v) {
  for(auto x : v)
    if(x == r) 
      return;
  v.push_back(r);
}

void RPmatcher::processResident(Rid r, vector<Rid>& rs) {
  if(params.verbosity > 2) {
    Rcpp::Rcout  << "#LOG: processResident " << r;
    if(r->inCouple()) {
      Rcpp::Rcout  << " couple = " << r->couple();
      if(r->couple()->isMatched()) 
	Rcpp::Rcout  << " current match = " << r->couple()->matchedTo();
    }
    else if(r->isMatched()) {
	Rcpp::Rcout  << " current match = " << r->matchedTo();
    }
    Rcpp::Rcout  << "\n";
  }

  if(r->inCouple())
    processCouple(r->couple(), rs);
  else {
    for(auto p : r->rol()) {
      if(p == r->matchedTo())
	break;
      if(!p->willAccept(r))
	continue;
      //accepted
      if(params.verbosity > 2)
	Rcpp::Rcout  << "#LOG: matching resident " << r << " into "
	     << "program " << p << "\n";
      if(r->isMatched())
	++totalRUnMatches;
      unmatch(r);
      for(auto x : p->match(r))
	processBumped(x, rs);
      r->match(p);
      ++totalRMatches;
      break;
    }
  }
}

void RPmatcher::processCouple(Cid c, vector<Rid>& rs) {
  for(auto pp : c->rol()) {
    if(pp == c->matchedTo())
       break;
    if(!willAccept(c, pp.first, pp.second))
      continue;

    //accepted
    if(c->isMatched())
      ++totalCUnMatches;
    unmatch(c->r1());
    unmatch(c->r2());
    if(pp.first == pp.second) {
      assert(pp.first != nilPid);
      for(auto x: pp.first->match(c))
	processBumped(x, rs);
    }
    else {
      if(pp.first != nilPid)
	for(auto x : pp.first->match(c->r1()))
	  processBumped(x, rs);
      if(pp.second != nilPid)
	for(auto x : pp.second->match(c->r2()))
	  processBumped(x, rs);
    }
    c->match(pp);
    ++totalCMatches;
    break;
  }
}

void RPmatcher::unmatch(Rid r) {
  if(r->isMatched()) {
    auto p = r->matchedTo();
    p->unmatch(r);
    r->unmatch();
  }
}

void RPmatcher::processBumped(Rid r, vector<Rid>& res) {
  r->unmatch();
  if(r->inCouple()) {
    unmatch(r->partner());
    addUnique(r->couple()->r1(), res);
    ++totalCUnMatches;
  }
  else {
    ++totalRUnMatches;
    addUnique(r, res);
  }
}

bool RPmatcher::chkMatch(vector<Rid>& rToProcess, vector<Rid>& matched) {
  //Roth Peranson maintain a program stack of programs that could
  //generate unstable pairings. These are the programs from which a
  //member of a couple was bumped during the matching process.  This
  //is a senseless "optimization" as it requires tracking things at
  //every stage of the matching. Instead we simply check all residents
  //and couples after matching all of them to see if any instability
  //exists. All unstable residents/couples back on the rToProcess
  //stack to be rematched.
  
  size_t initSize = rToProcess.size();
  for(auto rid : matched) {
    if(params.verbosity > 2)
      Rcpp::Rcout  << "#LOG: chkMatch processing resident " << rid << " current match = "
	   << rid->matchedTo() << "\n";

    if(!rid->inCouple()) {
      for(auto pid : rid->rol()) {
	if(pid == rid->matchedTo()) {
	  if(params.verbosity > 2) 
	    Rcpp::Rcout  << "#LOG: chkMatch resident " << rid << " in stable match ("
		 << rid->matchedTo() << ")\n";
	  break;
	}
	if(pid->willAccept(rid)) {
	  if(params.verbosity > 2)
	    Rcpp::Rcout  << "#LOG: chkMatch resident " << rid << " in unstable match "
		 << " prefers program " << pid << " to current match "
		 << rid->matchedTo() << "\n";
	  rToProcess.push_back(rid);
	  break;
	}
      }
    }
    else {
      auto cid = rid->couple();
      if(rid != cid->r1()) {
        Rcpp::Rcerr << "ERROR: residents to process stack contains non-r1 member of couple "
	     << rid << " couple = " << cid << "\n";
	continue;
      }
      for(auto ppid : cid->rol()) {
	if(ppid == cid->matchedTo()) {
	  if(params.verbosity > 2) 
	    Rcpp::Rcout  << "#LOG: chkMatch couple " << cid << " in stable match ("
		 << cid->matchedTo() << ")\n";
	  break;
	}
	if(willAccept(cid, ppid.first, ppid.second)) {
	  if(params.verbosity > 2)
	    Rcpp::Rcout  << "#LOG: chkMatch couple " << cid << " in unstable match "
		 << " prefers program pair " << ppid << " to current match "
		 << cid->matchedTo() << "\n";
	  rToProcess.push_back(rid);
	  break;
	}
      }
    }
  }
  return (rToProcess.size() == initSize);
}

