/***********[kprmatcher.cc]

-------

Implementation of "kprmatcher" class. 
This class implements the kpr matching algorithm. It alters updates the matching stored in
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

#include <params.h>
#include <io.h>
#include <kprmatcher.h>
#include <System.h>

using std::ostream;
using std::vector;
using std::pair;
//using std::Rcpp::Rcout ;
using std::cerr;
using std::random_shuffle;

template<typename T> 
//Add to r to deque only if r is not already in deque
void addUnique(T r, deque<T>& q) {
  for(auto x : q)
    if(x == r) 
      return;
  q.push_back(r);
}

void KPRmatcher::initData(Problem& prob) {
  for(auto &r: prob.Res()) {
    resNxtP.push_back(0);
    if(!r.inCouple()) {
      vector<int> counts(r.rol().size(), 0);
      resNApps.push_back(counts);
      resProcessQ.push_back(r); //auto convert to Rid
    }
    else if(r.couple()->isR1(r)) {
      vector<int> counts(r.couple()->rol().size(), 0);
      resNApps.push_back(counts);
      cplProcessQ.push_back(r);  //auto convert to Rid
    }
    else
      resNApps.push_back(vector<int>{});
  }
  maxNapps = 0;
}

void KPRmatcher::processBumped(Rid r) {
  //Resident or couple bumped from program 
  //Mark as unmatched add to the proper queue and
  //unmatch couple's partner
  r->unmatch();
  if(r->inCouple()) {
    unmatch(r->partner());
    addUnique(r->couple()->r1(), cplProcessQ);
    ++totalCUnMatches;
  }
  else {
    ++totalRUnMatches;
    addUnique(r, resProcessQ);
  }
}


void KPRmatcher::matchSingles(bool frmLast) {
  //Step 1, 3 and 6 of the algorithm.  Match Singles using deferred
  //match. If frmLast is true singles propose from their last proposal
  //(step 3), else they propose from the top of their ROL (step 6).
  //Note that step 1 (a standard GS matching of the singles) is also
  //achieved with this routine under the setting frmLast ==
  //true. (With only singles, resident cannot get back into a higher
  //ranked program---under the program's ranking matches into it are
  //monotonically non-dereasing. Hence we can always proceed downwards
  //in each resident's ROL.)
  
  if(params.verbosity > 1) 
    Rcpp::Rcout  << "#LOG: matchSingles(" << frmLast << ") "
	 << resProcessQ.size() << " singles\n";
  while(!resProcessQ.empty()) {
    Rid r = resProcessQ.front();
    resProcessQ.pop_front();
    
    if(params.verbosity > 3) 
      Rcpp::Rcout  << "#LOG: processing single " << r << " Next apply to #"
	   << resNxtP[r] << " = program " << r->rol()[resNxtP[r]] << "\n";
    for(size_t i = frmLast ? resNxtP[r] : 0;
	i < r->rol().size(); ++i) {
      Pid p = r->rol()[i];
      resNxtP[r] = i+1;     //r moves down its ROL
      
      if(p == r->matchedTo()) 
	break;
      
      resNApps[r][i] += 1; //r applies to prog #i once more
      if(resNApps[r][i] > maxNapps) {
	maxNapps = resNApps[r][i];
	if(params.verbosity > 2)
	  Rcpp::Rcout  << "#LOG: maxNapps increased to " << maxNapps << " (resident " 
	       << r << ", program " << p << ", rol index " << i << ")\n";

	if(params.kpr_cclim > 0 && maxNapps >= params.kpr_cclim) {
	  return;
	  resProcessQ.push_front(r); //not strictly needed but avoids "dangling" r
	}
      }

      if(!p->willAccept(r)) 
	continue;
      
      //p accepts r's application
      if(params.verbosity > 3)
	Rcpp::Rcout  << "#LOG: matching resident " << r << " into "
	     << "program " << p << " rol index = " << i
	     << " nxtApp " << resNxtP[r] << "\n";
      
      if(r->isMatched())
	++totalRUnMatches;
      unmatch(r); //in first phase r is not currently matched.
      for(Rid x: p->match(r))
	processBumped(x);
      r->match(p);
      ++totalRMatches;
      break;
    }
  }
}

void KPRmatcher::matchCouples(bool frmLast) {
  //Step 2 and 6 of the algorithm. Match couples. frmLast as in matchSingles
  if(params.verbosity > 1) 
    Rcpp::Rcout  << "#LOG: matchCouples(" << frmLast << ") "
	 << cplProcessQ.size() << " couples\n";

  while(!cplProcessQ.empty()) {
    Rid r1 = cplProcessQ.front();
    Cid c = r1->couple();
    cplProcessQ.pop_front();
    
    if(params.verbosity > 3) {
      Rcpp::Rcout  << "#LOG: processing couple " << c  << " Next apply to #"
	   << resNxtP[r1] << " = programs " << c->rol()[resNxtP[r1]] << "\n";
    }
    
    for(size_t i = frmLast ? resNxtP[r1] : 0;
	i < c->rol().size(); ++i) {
      auto pp = c->rol()[i];
      resNxtP[r1] = i+1;    //c moves down its ROL

      if(pp == c->matchedTo())
	break;

      resNApps[r1][i] += 1; //c applies to prog #i once more
      if(resNApps[r1][i] > maxNapps) {
	maxNapps = resNApps[r1][i];

	if(params.verbosity > 2)
	  Rcpp::Rcout  << "#LOG: maxNapps increased to " << maxNapps << " (couple " 
	       << c << ", program " << pp << ", rol index " << i << ")\n";

	if(params.kpr_cclim > 0 && maxNapps >= params.kpr_cclim) {
	  cplProcessQ.push_front(r1); //not strictly needed but avoids "dangling" c
	  return;
	}
      }
      
      if(!willAccept(c, pp.first, pp.second))
	continue;
      
      //pp accepts c's application
      if(params.verbosity > 3)
	Rcpp::Rcout  << "#LOG: matching couple " << c << " into "
	     << "program pair " << pp << " rol index = " << i
	     << " nxtApp " << resNxtP[r1] << "\n";

      if(c->isMatched())
	++totalCUnMatches;;
      unmatch(c->r1());
      unmatch(c->r2());
      if(pp.first == pp.second) {
	assert(pp.first != nilPid);
	for(auto x: pp.first->match(c))
	  processBumped(x);
      }
      else {
	if(pp.first != nilPid)
	  for(auto x : pp.first->match(c->r1()))
	    processBumped(x);
	if(pp.second != nilPid)
	  for(auto x : pp.second->match(c->r2()))
	    processBumped(x);
      }
      c->match(pp);
      ++totalCMatches;
      break;
    }
  }
}


bool KPRmatcher::match_(Problem& prob) {
  initData(prob);
  //Step 1
  matchSingles(true);
  if(params.verbosity > 0) {
    Rcpp::Rcout  << "#Initial DA match of singles completed\n";
    if(params.verbosity > 4) {
      Rcpp::Rcout  << "#Current match:\n";
      prob.printMatch(true, false);
    }
  }
  //Step 4. Loop over steps 2 and 3. Start everyone off with proposals to their top rank.
  for(size_t i = 0; i < resNxtP.size(); i++)
    resNxtP[i] = 0;

  if(params.verbosity > 0)
    Rcpp::Rcout  << "#Step 4 (iternate 2 and 3)\n";

  while(cplProcessQ.size() + resProcessQ.size() > 0) {
    int cplN = cplProcessQ.size();
    matchCouples(true);
    if(cplN > 0 && params.verbosity > 2) {
      Rcpp::Rcout  << "#Matched Couples\n";
      Rcpp::Rcout  << "#Number of bumped singles = " << resProcessQ.size() << "\n";
      if(params.verbosity > 4) {
	Rcpp::Rcout  << "#Current match:\n";
	prob.printMatch(true, false);
      }
    }
    int resN = resProcessQ.size();
    matchSingles(true);
    if(resN > 0 && params.verbosity > 2) {
      Rcpp::Rcout  << "#Matched of Singles\n";
      Rcpp::Rcout  << "#Number of bumped couples = " << cplProcessQ.size() << "\n";
      if(params.verbosity > 4) {
	Rcpp::Rcout  << "#Current match:\n";
	prob.printMatch(true, false);
      }
    }
  }

  if(params.verbosity > 0) {
    Rcpp::Rcout  << "#Initial round robin match (step 2&3) completed (maxNapps = "
	 << maxNapps << ")\n";
    if(params.verbosity > 4) {
      Rcpp::Rcout  << "#Current match:\n";
      prob.printMatch(true, false);
    }
  }
  
  int iters {0};
  //step 6.
  while(!chkMatch(prob)) {
    ++iters;
    if(params.verbosity > 1) 
      Rcpp::Rcout  << "#LOG: Unstable match. Iteration " << iters
	   << " maxNapps = " << maxNapps << "(limit = " << params.kpr_cclim << ")\n";

    for(auto& r: prob.Res()) {
      if(!r.inCouple())
	resProcessQ.push_back(r);
      else if(r.couple()->isR1(r))
	cplProcessQ.push_back(r);
    }

    if(params.verbosity > 2)
      Rcpp::Rcout  << "#LOG: Processing " << resProcessQ.size() << " residents and "
	   << cplProcessQ.size() << " couples\n";

    while(cplProcessQ.size() +resProcessQ.size() > 0) {
      int cplN = cplProcessQ.size();
      matchCouples(false); //start from top for each match
      if(params.kpr_cclim > 0 && maxNapps >= params.kpr_cclim) {
	if(params.verbosity > 0)
	  Rcpp::Rcout  << "#Max applications exceeded. No match found\n";
	return false;
      }

      if(params.verbosity > 2 && cplN > 0) {
	Rcpp::Rcout  << "#Matched Couples\n";
	Rcpp::Rcout  << "#Number of bumped singles = " << resProcessQ.size() << "\n";
	if(params.verbosity > 4) {
	  Rcpp::Rcout  << "#Current match:\n";
	  prob.printMatch(true, false);
	}
      }
      int resN = resProcessQ.size();
      matchSingles(false); //start from top for each match
      if(params.kpr_cclim > 0 && maxNapps >= params.kpr_cclim) {
	if(params.verbosity > 0)
	  Rcpp::Rcout  << "#Max applications exceeded. No match found\n";
	return false;
      }

      if(params.verbosity > 2 && resN > 0) {
	Rcpp::Rcout  << "#Match of Singles\n";
	Rcpp::Rcout  << "#Number of bumped couples = " << cplProcessQ.size() << "\n";
	if(params.verbosity > 4) {
	  Rcpp::Rcout  << "#Current match:\n";
	  prob.printMatch(true, false);
	}
      }
    }
  }
   
  return true;
}

void KPRmatcher::unmatch(Rid r) {
  if(r->isMatched()) {
    auto p = r->matchedTo();
    p->unmatch(r);
    r->unmatch();
  }
}

bool KPRmatcher::chkMatch(Problem& prob) {
  //check if current match is stable
  if(params.verbosity > 2) 
    Rcpp::Rcout  << "#LOG: checking match\n";

  for(auto& r : prob.Res()) {
    if(params.verbosity > 3)
      Rcpp::Rcout  << "#LOG: chkMatch processing resident " << r.id() << " current match = "
	   << r.matchedTo() << "\n";

    if(!r.inCouple()) {
      for(auto pid : r.rol()) {
	if(pid == r.matchedTo()) {
	  if(params.verbosity > 3) 
	    Rcpp::Rcout  << "#LOG: chkMatch resident " << r.id() << " in stable match ("
		 << r.matchedTo() << ")\n";
	  break;
	}
	if(pid->willAccept(r)) {
	  if(params.verbosity > 2)
	    Rcpp::Rcout  << "#LOG: chkMatch resident " << r.id() << " in unstable match "
		 << " prefers program " << pid << " to current match "
		 << r.matchedTo() << "\n";
	  return false;
	}
      }
    }
    else if (r.couple()->isR1(r)) {
      auto cid = r.couple();
      for(auto ppid : cid->rol()) {
	if(ppid == cid->matchedTo()) {
	  if(params.verbosity > 3) 
	    Rcpp::Rcout  << "#LOG: chkMatch couple " << cid << " in stable match ("
		 << cid->matchedTo() << ")\n";
	  break;
	}
	if(willAccept(cid, ppid.first, ppid.second)) {
	  if(params.verbosity > 2)
	    Rcpp::Rcout  << "#LOG: chkMatch couple " << cid << " in unstable match "
		 << " prefers program pair " << ppid << " to current match "
		 << cid->matchedTo() << "\n";
	  return false;
	}
      }
    }
    else 
      continue;
  }
  if(params.verbosity > 1) 
    Rcpp::Rcout  << "#LOG: match stable\n\n";
  return true;
}

