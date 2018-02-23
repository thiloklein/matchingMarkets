/***********[rcp.h]

-------

Class header the Resident, Couple, and Program classes.
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
#ifndef RCP_h
#define RCP_h

#include <iostream>
#include <vector>
#include <map>
#include <set>
#include <assert.h>
#include "params.h"
#include "problem.h"

using std::set;
using std::map;
using std::vector;
using std::cout;
using std::ostream;

class Problem;
class Resident;
class Couple;



//Integer oncrete classes that index into the underlying objects, via
//a pointer to the a problem object.
class Rid;
class Pid;
class Cid;
typedef std::pair<Pid, Pid> PidPair;

class Rid {
public:
  explicit Rid(int i) : _id {i} {}
  Rid() : _id {} {}
  operator size_t() const { return static_cast<size_t>(_id); }
  bool operator<(Rid o) const { return _id < o._id; }
  int  id() const { return _id; }
  Resident* operator->() const { return prob->ithRes(_id); }
  friend class Problem;
  friend ostream& operator<<(ostream& os, const Rid& r);

private:
  int _id;
  static Problem* prob;
};

class Cid {
public:
  explicit Cid(int i) : _id {i} {}
  Cid() : _id {} {}
  operator size_t() const { return static_cast<size_t>(_id); }
  bool operator<(Cid o) const { return _id < o._id; }
  int id() const { return _id; };
  Couple* operator->() const { return prob->ithCpl(_id); }

  friend class Problem;
  friend ostream& operator<<(ostream& os, const Cid& c);

private:
  int _id;
  static Problem* prob;
};

class Pid {
public:
  explicit Pid(int i) : _id {i} {}
  Pid() : _id {} {}

  operator size_t() const { return static_cast<size_t>(_id); }
  bool operator<(Pid o) const { return _id < o._id; }
  int id() const { return _id; };
  Program* operator->() const { return prob->ithProg(_id); }

  friend class Problem;
  friend ostream& operator<<(ostream& os, const Pid& p);

private:
  int _id;
  static Problem* prob;
};

//Constants
const Pid nilPid {-1};
const PidPair nilPPid {nilPid, nilPid};
const Rid nilRid {-1};
const Cid nilCid {-1};

//PidPair processing
bool willAccept(Cid c, Pid p1, Pid p2);

inline bool operator<(const PidPair& x, const PidPair& y) 
{
  if(x.first == y.first)
    return x.second < y.second;
  return x.first < y.first;
}

//Classes holding actual data

//Classes 
class Resident {
 public:
  Resident();
  Resident(Rid ident, vector<Pid> rankedPrograms, int couple = nilCid);

  operator Rid() const { return _id; } //auto and explicit converters to Rid.
  Rid id() const { return _id; }

  int rankOf(Pid p) const;
  bool prefers(Pid p1, Pid p2) const { return rankOf(p1) < rankOf(p2); }
  bool ranks(Pid p) const { return rankOf(p) <= static_cast<int>(_rol.size()); }
  Pid matchedTo() const { return _match; }
  bool isMatched() const { return _match != nilPid; }
  void match(Pid pid); //forces a match
  void unmatch();
  bool inCouple() const { return _couple != nilCid; }
  Cid couple() const { return _couple; }
  Rid partner() const;
  const vector<Pid>& rol() const { return _rol; }
  //Cycle Checking
  static int maxMatchCount;
  void clearMatchCounts();
  int getMatchCount(Pid pid);

  friend class Problem;
  friend ostream& operator<<(ostream& os, const Resident& r);


 private:
  Rid _id;
  vector<Pid> _rol;
  vector<int> matchCount;
  map<Pid, int> pid2rank;
  Cid _couple;
  Pid _match;
};

class Couple {
 public:
  Couple();
  Couple(Cid ident, Rid res1, Rid res2, vector<PidPair> rankedPrograms);

  operator Cid() const { return _id; } 
  Cid id() const { return _id; }
  int rankOf(PidPair p) const;
  bool prefers(PidPair p1, PidPair p2) const  { return rankOf(p1) < rankOf(p2); }
  bool ranks(PidPair p) const { return rankOf(p) <= static_cast<int>(_rol.size()); }
  bool isR1(Rid r) const { return r == _r1; }
  bool isR2(Rid r) const { return r == _r2; }
  bool ranksRi(Pid p, Rid r) const; //p is ranked by r[1/2]
  PidPair matchedTo() const { return std::make_pair(_r1->matchedTo(), _r2->matchedTo()); }
  bool isMatched() const { return _r1->isMatched() || _r2->isMatched(); }
  void match(PidPair p); //force a match
  void unmatch();

  Rid r1() const { return _r1; }
  Rid r2() const { return _r2; }
  const vector<PidPair>& rol() const { return _rol; }

  friend class Problem;
  friend ostream& operator<<(ostream& os, const Couple& c);
 
private:
  Cid _id;
  vector<PidPair> _rol;
  map<PidPair, int> pid2rank;
  Rid _r1;
  Rid _r2;
};

class Program {
public:
  Program();
  Program(Pid ident, int quota, vector<Rid> rankedResidents);

  operator Pid() const { return _id; } 
  Pid id() const { return _id; }
  int rankOf(Rid r) const;
  bool prefers(Rid r1, Rid r2) const { return rankOf(r1) < rankOf(r2); }
  bool ranks(Rid r) const { return rankOf(r) <= static_cast<int>(_rol.size()); }
  Rid minRes() const; 
  Rid min2ndRes() const;
  bool inProgram(const Rid r) const;
  bool willAccept(Rid r) const;
  bool willAccept(Cid c) const;
  int nmatched() const { return _accepted.size(); }
  vector<Rid> match(Rid r);
  vector<Rid> match(Cid c);
  void unmatch(Rid r);

  int quota() const { return _quota; }
  vector<Rid> accepted() const;
  const vector<Rid>& rol() const { return _rol; }

  struct RidCmp {
    Pid _p;
     bool operator()(Rid r1, Rid r2) const {
     return _p->prefers(r1, r2); 
     }
     RidCmp(Pid p) : _p{p} {}
   };

  friend class Problem;
  friend ostream& operator<<(ostream& os, const Program& p);

private:
  Pid _id;
  int _quota;
  vector<Rid> _rol;
  map<Rid, int> rid2rank;
  set<Rid, RidCmp> _accepted;
};

//IO
ostream& operator<<(ostream& os, const Resident& r);
ostream& operator<<(ostream& os, const Couple& c);
ostream& operator<<(ostream& os, const Program& p);

//IO
inline ostream& operator<<(ostream& os, const Rid& r) {
  os << r._id ;
  return os;
}

inline ostream& operator<<(ostream& os, const Pid& p) {
  os << p._id ;
  return os;
}

inline ostream& operator<<(ostream& os, const Cid& c) {
  os << c._id << " [" << c->r1() << "," << c->r2() << "]";
  return os;
}


#endif
