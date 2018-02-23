/***********[rcp.cc]

-------

Implementation of Resident, Couple, and Program classes.
These classes maintain information about the preferences and the current matching.

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
#include <algorithm>
#include <rcp.h>
#include <io.h>
//using std::cout;
//using std::cerr;
using std::ostream;

//Residents
Resident::Resident() :
  _id{nilRid}, _rol{}, matchCount{}, pid2rank{},
  _couple{nilCid}, _match{nilPid}
{}

Resident::Resident(Rid ident, vector<Pid> rankedPrograms, int couple) :
    _id {ident},
    _rol {rankedPrograms},
    matchCount(_rol.size()),
    pid2rank {},
    _couple {couple},
    _match {nilPid}
  { 
    assert(ident.id() >= 0);
    for(size_t i=0; i < rankedPrograms.size(); ++i) 
      pid2rank[rankedPrograms[i]] = i;
  }

int Resident::rankOf(Pid p) const {
  if (p == nilPid) return _rol.size();
  auto it = pid2rank.find(p);
  if(it == pid2rank.end())
    return std::numeric_limits<int>::max();      
  else
    return it->second;
}

void Resident::match(Pid pid) { //does not check ranking
  assert(inCouple() || ranks(pid)); //couple residents don't have their own rols

  _match = pid;
  if(!inCouple() && params.verbosity > 4)
    Rcpp::Rcout  << "#LOG: matching resident " << _id << " to program " << pid << "\n";
  
  //Cycle Checking
  if(pid != nilPid) {
    //don't count matches to null in cycle checking
    size_t i = rankOf(pid);
    if(i < _rol.size()) {
      matchCount[i] += 1;
      if(matchCount[i] > maxMatchCount) {
	if(params.verbosity > 4) 
	  Rcpp::Rcout  << "#LOG: new maxMatchcount = " << maxMatchCount << "\n";
	maxMatchCount = matchCount[i];
      }
    }
  }
}

void Resident::unmatch() {
  match(nilPid); 
  if(!inCouple() && params.verbosity > 4)
    Rcpp::Rcout  << "#LOG: unmatching resident " << _id << "\n";
}

Rid Resident::partner() const {
    if(!inCouple()) return nilRid;
    return (couple()->isR1(_id)) ? couple()->r2() : couple()->r1();
}

void Resident::clearMatchCounts() {
  for(size_t i = 0; i < matchCount.size(); i++)
    matchCount[i] = 0;
}

int Resident::getMatchCount(Pid pid) {
  size_t i = rankOf(pid);
  if(i < matchCount.size())
    return matchCount[i];
  else
    return 0;
}

//Define static class members
int Resident::maxMatchCount;

//Couples
Couple::Couple() :
  _id{nilCid}, _rol{}, pid2rank{}, _r1{nilRid}, _r2{nilRid} 
  {}

Couple::Couple(Cid ident, Rid res1, Rid res2, vector<PidPair> rankedPrograms) :
  _id {ident},
  _rol {rankedPrograms},
  pid2rank {},
  _r1 {res1},
  _r2 {res2}
  {
    assert(ident.id() >= 0);
    for(size_t i=0; i < rankedPrograms.size(); ++i) 
      pid2rank[rankedPrograms[i]] = i;
  }

int Couple::rankOf(PidPair p) const {
  if(p == nilPPid)
    return _rol.size();
  auto it = pid2rank.find(p);
  if(it == pid2rank.end())
    return std::numeric_limits<int>::max();      
  else
    return it->second;
}

bool Couple::ranksRi(Pid p, Rid r) const { //ranked by r[1/2]
  if(r == _r1) {
    auto ranksP = [p](PidPair x) {return x.first == p;};
    auto it = std::find_if(_rol.begin(), _rol.end(), ranksP);
    return it != _rol.end();
  }
  else {
    auto ranksP = [p](PidPair x) {return x.second == p;};
    auto it = std::find_if(_rol.begin(), _rol.end(), ranksP);
    return it != _rol.end();
  }
}

void Couple::match(PidPair p) { 
  if(params.verbosity > 4)
    Rcpp::Rcout  << "#LOG: matching couple " << _id << " to programs " << p << "\n";
  assert(ranks(p));
  _r1->match(p.first);
  _r2->match(p.second);
}

void Couple::unmatch() { 
 if(params.verbosity > 4)
    Rcpp::Rcout  << "#LOG: unmatching couple " << _id << "\n";
 match({nilPid, nilPid});
}

//Programs
Program::Program() :
  _id {nilPid}, _quota {std::numeric_limits<int>::max()},
  _rol{}, rid2rank{}, 
  _accepted{RidCmp{_id}} 
 {}

Program::Program(Pid ident, int quota, vector<Rid> rankedResidents) :
  _id {ident},
  _quota {quota},
  _rol {rankedResidents},
  rid2rank {},
  _accepted{RidCmp{_id}}
  {
    assert(_id >= 0);
    assert(_quota > 0);
    for(size_t i=0; i < _rol.size(); ++i) 
      rid2rank[_rol[i]] = i;
  }

int Program::rankOf(Rid r) const {
  if(r == nilRid) 
    return _rol.size();

  auto it = rid2rank.find(r);
  if(it == rid2rank.end())
    return std::numeric_limits<int>::max();      
  else
    return it->second;
}

Rid Program::minRes() const {
  if(_quota <= 0)
    return nilRid;
  if(nmatched() == _quota)
    return *(--_accepted.end());
  return nilRid;
}

Rid Program::min2ndRes() const {
  if(_quota <= 1)
    return nilRid;
  if(nmatched() == _quota)
    return *(--(--_accepted.end()));
  if(nmatched() == _quota-1)
    return *(--_accepted.end());
  return nilRid;
}

bool Program::inProgram(const Rid r) const {
  return (_accepted.find(r) != _accepted.end());
}

bool Program::willAccept(Rid r) const {
  //Note: r is acceptable if they are preferred some currently 
  //accepted resident, *OR* if they are already in the program.
  return _quota > 0 && rankOf(r) <= rankOf(minRes());
}

bool Program::willAccept(Cid c) const {
  //Note: c is acceptable if the program would prefer to accept both
  //members of the couple to the current acceptances *OR* c is
  //already in the program.
  if(inProgram(c->r1()) && inProgram(c->r2()))
    return true;
  auto rm2 = rankOf(min2ndRes());
  return _quota > 1 && rankOf(c->r1()) <= rm2 && rankOf(c->r2()) <= rm2;
}

vector<Rid> Program::match(Rid r) {
  //Raw match. r should already be unmatched and not nil
  vector<Rid> bumped {};
  assert(r != nilRid);
  assert(!r->isMatched());

  if(params.verbosity > 4)
    Rcpp::Rcout  << "#LOG: placing resident " << r << " in program " << _id << "\n";
  if(_quota == 0) {
    Rcpp::Rcerr << "ERROR: resident placement failed quota is zero\n";
    return bumped;
  }

  if(nmatched() == _quota) {
    auto it = --_accepted.end();
    bumped.push_back(*it);
    _accepted.erase(it);
  }
  _accepted.insert(r);
  if(params.verbosity > 4 && !bumped.empty())
    Rcpp::Rcout  << "#LOG: placement bumped residents " << bumped << "\n";
  return bumped;
}

vector<Rid> Program::match(Cid c) {
  //Raw Match. Both members of c should already be unmatched.
  assert(c->r1() != nilRid && c->r2() != nilRid);
  assert(!c->isMatched());

  vector<Rid> bumped {};
  if(params.verbosity > 4)
    Rcpp::Rcout  << "#LOG: placing couple " << c << "in program " << _id << "\n"; 
  if(_quota < 2) {
    if(params.verbosity > 4)
      Rcpp::Rcout  << "#LOG: couple placement failed quota is < 2\n";
    return bumped;
  }
  auto it = _accepted.end();
  while(nmatched() >= _quota-1){
    bumped.push_back(*(--it));
    _accepted.erase(it);
  }
  _accepted.insert(c->r1());
  _accepted.insert(c->r2());

  if(params.verbosity > 4 && !bumped.empty())
    Rcpp::Rcout  << "#LOG: placement bumped residents " << bumped << "\n";
  return bumped;
}

void Program::unmatch(Rid r) {
  if(params.verbosity > 4)
    Rcpp::Rcout  << "#LOG: removing resident " << r << " from program "
	 << _id << "\n";
  _accepted.erase(r);
}

vector<Rid> Program::accepted() const 
{
  return vector<Rid>(_accepted.begin(), _accepted.end());
}

//process a pair of programs
bool willAccept(Cid c, Pid p1, Pid p2)
{
  if(p1 == nilPid && p2 == nilPid)
    return true;
  if(p1 == p2) 
    return p1->willAccept(c);
  return 
    (p1 == nilPid || p1->willAccept(c->r1()))
    && (p2 == nilPid || p2->willAccept(c->r2()));
}

//Indexer statics.

Problem* Rid::prob;
Problem* Cid::prob;
Problem* Pid::prob;

// IO
ostream& operator<<(ostream& os, const Resident& r) {
  os << "#Resident " << r.id() << ". ";
  os << " match = " << r.matchedTo() << " ";
  if(r.inCouple())
    os << "in couple " << r.couple();
  else {
    os << "Not in couple ";
    os << "ROL = " << r.rol() << " ";
    os << "pid2rank = ";
    for(auto item : r.pid2rank)
      os << "[" << item.first << "," << item.second << "] ";
  }
  os << "\n";
  return os;
}

ostream& operator<<(ostream& os, const Couple& c) {
  os << "#Couple " << c.id() << ". ";
  os << "r1 = " << c.r1() << " r2 = " << c.r2();
  os << " match = " << c.matchedTo() << " ";
  os << "ROL = " << c.rol() << " ";
  os << "pid2rank = ";
  for(auto item : c.pid2rank)
    os << "[" << item.first << "," << item.second << "] ";
  os << "\n";
  return os;
}


ostream& operator<<(ostream& os, const Program& p) {
  os << "#Program " << p.id() << ". ";
  os << "quota = " << p.quota() << " ";
  os << "accepted  = " << p.accepted() << " ";
  os << "ROL = " << p.rol() << " ";
  os << "rid2rank = ";
  for(auto item : p.rid2rank)
    os << "[" << item.first << "," << item.second << "] ";
  os << "\n";
  return os;
}
