/***********[problem.cc]

-------

"Problem" class implementation

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
#include <fstream>
#include <sstream>
#include <problem.h>
#include <params.h>
#include <io.h>
#include <rcp.h>
#include <list>

#include <Rcpp.h>
using namespace Rcpp;

using std::istringstream;
//using std::cout;

Problem::Problem() : 
    errMsg {}, probOK {true},
    resIDs {}, progIDs {}, cplIDs {}, progsRanked {}, resRanked {},
    residents {}, programs {}, couples {}
{
  Rid::prob = Cid::prob = Pid::prob = this;
}

/*bool Problem::readProblem(string filename) {
  std::ifstream in {filename};
  vector<string> ilines;
  string line;
  while(getline(in, line))
    ilines.push_back(line);
  for(auto l : ilines) {
    if(l.size() == 0)
      continue;
    switch( l[0] ) {
    case ' ':
      break;
    case '#':
      break;
    case 'r':
      readResident(l);
      break;
    case 'c':
      readCouple(l);
      break;
    case 'p':
      readProgram(l);
      break;
    default:
      postError("Input ERROR: line \"" + l + "\" from input is invalid\n");
    }
  }
  furtherInputChecks();
  clearErrVecs();
  if(ok())
     postProcess();
  return ok();
}*/

//S
bool Problem::readProblem(Rcpp::List residents, Rcpp::List couples, Rcpp::List programs) {
 
  //residents
  for(auto i = 0; i < residents.size(); i++) {
    //Rcpp::Rcout << "resident + " << Rcpp::as<std::string>(residents[i]) << std::endl;
    readResident(Rcpp::as<std::string>(residents[i]));
  }
  
  //programms
  for(auto i = 0; i < programs.size(); i++) {
    //Rcpp::Rcout << "program + " << Rcpp::as<std::string>(programs[i]) << std::endl;
    readProgram(Rcpp::as<std::string>(programs[i]));
  }

  //couples
  for(auto i = 0; i < couples.size(); i++) {
    //Rcpp::Rcout << "couple + " << Rcpp::as<std::string>(couples[i]) << std::endl;
    //test if couples are part of the problem
    string p = as<std::string>(couples[0]);
    if (!p.empty()) {
      readCouple(Rcpp::as<std::string>(couples[i]));
    }
  }
 
  furtherInputChecks();
  clearErrVecs();
  if(ok())
    postProcess();
  return ok(); 
}

bool Problem::chkID(int rid, unordered_set<int>& Ids, string errmsg) {
  auto fnd = Ids.find(rid);
  if(fnd != Ids.end()) {
    postError(errmsg);
    return false;
  }
  else {
    Ids.insert(rid);
    return true;
  }
}

void Problem::readResident(string l) {
  //Format:
  //"r <resident id==int> <rol>"
  //Where rol is a sequence of program ids most prefered first.
  istringstream iss {l};
  char c;
  int rid, pid;
  vector<int> pids;

  iss >> c >> rid;
  while(iss >> pid) {
    pids.push_back(pid);
    progsRanked.push_back(pid);
  }

  if(rid < 0) {
    postError("Input ERROR: negative Resident ID in resident spec.\n");
    return;
  }
  if(!chkID(rid, resIDs, "Input ERROR: Duplicate resident ID in resident specs.\n"))
    return;

  if(static_cast<int>(residents.size()) <= rid)
    residents.resize(rid+1);

  vector<Pid> _pids;
  for(auto p : pids) 
    _pids.push_back(Pid{p});
  residents[rid] = Resident(Rid{rid}, _pids);
}
    
void Problem::readCouple(string l) {
  //Format:
  //"c <couple id> <r1id> <r2id> <rol>"
  //The couple id, followed by the two resident ids.
  //Neither resident can have appeared before.
  //rol is an even number of program ids. They are regarded 
  //as pairs with the most preferred pair comming as the first
  //two program ids. 
  //Note -1 is a legit progam id, denoting the null program. 
  //(e.g., the pair 3 -1 indicates a preference for res1 matching
  //to program 3 and res2 having a null match. 
  istringstream iss {l};
  char c;
  int r1id, r2id, cid, pid;
  vector<int> pids;

  iss >> c >> cid >> r1id >> r2id;
  while(iss >> pid) {
    pids.push_back(pid);
    progsRanked.push_back(pid);
  }

  if(pids.size() % 2 != 0) {
    postError("Input ERROR: Couple input had odd number or programs specified (not pairs)\n");
    return;
  }

  if(!chkID(r1id, resIDs, "Input ERROR: Duplicate resident ID in couple spec.\n"))
    return;
  if(r1id != r2id &&
     !chkID(r2id, resIDs, "Input ERROR: Duplicate resident ID in couple spec.\n"))
    return;
  if(!chkID(cid, cplIDs, "Input ERROR: Duplicate couple ID in couple specs.\n"))
    return;
  if(r1id < 0 || r2id < 0) {
    postError("Input ERROR: negative resident ID in couple spec\n");
    return;
  }
  
  if(static_cast<int>(couples.size()) <= cid)
    couples.resize(cid+1);
  vector<PidPair> ppairs {};
  for(size_t i = 0; i < pids.size(); i += 2) 
    ppairs.push_back({Pid{pids[i]}, Pid{pids[i+1]}});
  couples[cid] = Couple(Cid{cid}, Rid{r1id}, Rid{r2id}, ppairs);

  if(static_cast<int>(residents.size()) <= r1id)
    residents.resize(r1id+1);
  residents[r1id] = Resident(Rid{r1id}, {}, Cid{cid});
  if(static_cast<int>(residents.size()) <= r2id)
    residents.resize(r2id+1);
  residents[r2id] = Resident(Rid{r2id}, {}, Cid{cid});
}

void Problem::readProgram(string l) {
  //Format
  //"p <program id> <quota> <rol>"
  istringstream iss {l};
  char c;
  int pid, rid, quota;
  vector<int> rids;

  iss >> c >> pid >> quota;
  while(iss >> rid) {
    rids.push_back(rid);
    resRanked.push_back(rid);
  }

  if(!chkID(pid, progIDs, "Input ERROR: Duplicate program ID in program specs.\n"))
    return;

  if(static_cast<int>(programs.size()) <= pid)
    programs.resize(pid+1);
  
  vector<Rid> _rids;
  for(auto r : rids) 
    _rids.push_back(Rid{r});

  programs[pid] = Program(Pid{pid}, quota, _rids);
}

void Problem::furtherInputChecks() {
  for(auto pid : progsRanked) {
    if(pid != -1) {
      auto fnd = progIDs.find(pid);
      if(fnd == progIDs.end() ) 
	postError("Input ERROR: Resident or Couple ranked unspecified program.\n");
    }
  }
  for(auto rid : resRanked) {
    auto fnd = resIDs.find(rid);
    if(fnd == resIDs.end()) 
      postError("Input ERROR: Program unspecified resident.\n");
  }
}

void Problem::clearErrVecs() {
  resIDs = unordered_set<int> {};
  progIDs = unordered_set<int> {};
  cplIDs = unordered_set<int> {};
  progsRanked = vector<int> {};
  resRanked = vector<int> {};
}

void Problem::postProcess() {
  //remove rankings that are unreciprocated
  //friend of all classes so can access rol vectors.
  for(auto& r : residents) {
    size_t j {0};
    for(auto p : r._rol) {
      if(p->ranks(r)) 
	r._rol[j++] = p; //p ranks r
    }
    r._rol.resize(j);
    r.pid2rank.clear();
    for(size_t i=0; i < r._rol.size(); ++i) 
      r.pid2rank[r._rol[i]] = i;
  }
  for(auto& c : couples) {
    size_t j {0};
    for(auto p : c._rol)
      if((p.first == nilPid || p.first->ranks(c.r1())) 
	 && (p.second == nilPid || p.second->ranks(c.r2())))
	c._rol[j++] = p;
    c._rol.resize(j);
    c.pid2rank.clear();
    for(size_t i=0; i < c._rol.size(); ++i) 
      c.pid2rank[c._rol[i]] = i;


  }
  for(auto& p : programs) {
    size_t j {0};
    for(auto r : p._rol) {
      if(r->inCouple()) {
	if(r->couple()->ranksRi(p, r))
	  p._rol[j++] = r;
      }
      else {
	if(r->ranks(p))
	  p._rol[j++] = r;
      }
    }
    p._rol.resize(j);
    p.rid2rank.clear();
    for(size_t i=0; i < p._rol.size(); ++i) 
      p.rid2rank[p._rol[i]] = i;
  }
}

//Access to objects via IDs

Program* Problem::ithProg(const int id) {
  assert(id >= 0 && static_cast<size_t>(id) < programs.size());
  return &programs[id];
}

Resident* Problem::ithRes(const int id) {
  assert(id >= 0 && static_cast<size_t>(id) < residents.size());
  return &residents[id];
}

Couple* Problem::ithCpl(const int id) {
  assert(id >= 0 && static_cast<size_t>(id) < couples.size());
  return &couples[id];
}

//SVEN
Rcpp::List Problem::returnMatch(bool fndMatch) {
  List matchList;
  
  //Create vector for matching
  Rcpp::NumericVector matchResult(Res().size());
  Rcpp::NumericVector matchResidents(Res().size());
  int i; i = 0;
  for (auto& r : Res()) {
    matchResult[i] = r.matchedTo()+1;
    matchResidents[i] = r.id()+1;
    i++;
  }
  
  return Rcpp::List::create(
    Rcpp::Named("matchResultResident", matchResult),
    Rcpp::Named("ResidentID", matchResidents)
  );
}

//SVEN
string Problem::summaryMatch(bool fndMatch) {
  std::stringstream buffer;
  
  int resNotMatched {0};
  int cplNotMatched {0};
  int progSpareCap {0};
  int nSingRes {0};
  int resGotTopRank {0};
  int cplGotTopRank {0};
  int prgGotTopRank {0};
  
  double resAveRank {0};
  double cplAveRank {0};
  double prgAveRank {0};
  
  buffer << "#Single Residents:" << std::endl;
  
  for(const auto&r : Res()) 
    if(!r.inCouple()) {
      buffer << "#Resident " << r.id() << ": ";
      buffer << "match = " << r.matchedTo() << " Res Ranking = "
           << r.rankOf(r.matchedTo()) << "/" << r.rol().size() << std::endl;
      nSingRes++;
      if(!r.isMatched()) {
        ++resNotMatched;
      }
      else
        resAveRank += r.rankOf(r.matchedTo());
      
      if(r.rankOf(r.matchedTo()) == 0)
        ++resGotTopRank;
      
    }
    
    buffer << "\n#Couples:" << std::endl;
    for(const auto&c : Cpl()) {
      buffer << "#Couple " << c.id() << ": ";
      buffer << "match = " << c.matchedTo() << " Cpl Ranking = "
           << c.rankOf(c.matchedTo()) << "/" << c.rol().size() << std::endl;
      if(!c.isMatched())
        ++cplNotMatched;
      else
        cplAveRank += c.rankOf(c.matchedTo());      
      
      if(c.rankOf(c.matchedTo()) == 0)
        ++cplGotTopRank;
    }
    
    int matchedProgs {0};
    buffer << "\n#Programs:" << std::endl;
    for(const auto&p : Prog()) {
      buffer << "#Program " << p.id() << ": ";
      buffer << "spares = " << p.quota() - p.accepted().size();
      progSpareCap += p.quota() - p.accepted().size();
      buffer << " accepted = [";
      for(const auto &res : p.accepted()) 
        buffer << res << " ";
      buffer << "] " << " Prog rankings = [";
      double aveRank {0};
      for(const auto &res : p.accepted()) {
        buffer << p.rankOf(res) << " ";
        aveRank += p.rankOf(res);
        if(p.rankOf(res) == 0)
          ++prgGotTopRank;
      }
      buffer << "]" << "/" << p.rol().size();
      if(p.accepted().size() > 0) {
        prgAveRank += aveRank/p.accepted().size();
        ++matchedProgs;
        buffer << " ave Prog rank (accepted) = " << aveRank/p.accepted().size();
      }
      buffer << std::endl;
    }
    
    buffer << "\n#Matching Summary Stats:" << std::endl;
    buffer << "#Unmatched Singles: " << resNotMatched << std::endl;
    buffer << "#Unmatched Couples: " << cplNotMatched << std::endl;
    buffer << "#Unmatched Program slots: " << progSpareCap << std::endl;
    
    if(nSingRes-resNotMatched > 0)
      buffer << "#Ave Resident Rank of their matching = "
           << resAveRank/(nSingRes-resNotMatched) << std::endl;
      buffer << "#Num Residents getting their top rank = "
           << resGotTopRank << std::endl;
      
      if(Cpl().size() - cplNotMatched > 0)
        buffer << "#Ave Couple Rank of their matching = "
             << cplAveRank/(Cpl().size()-cplNotMatched) << std::endl;
        buffer << "#Num Couples getting their top rank = "
             << cplGotTopRank << std::endl;
        
        if(matchedProgs > 0) 
          buffer << "#Ave Program Rank of their matched residents "
               << prgAveRank/matchedProgs << std::endl;
          buffer << "#Num Programs getting their top rank = "
               << prgGotTopRank;
          buffer << std::endl;
  
  return buffer.str();
}

//SVEN
Rcpp::List Problem::getStats(bool fndMatch) {
  int resNotMatched {0};
  int cplNotMatched {0};
  int progSpareCap {0};
  int nSingRes {0};
  int resGotTopRank {0};
  int cplGotTopRank {0};
  int prgGotTopRank {0};
  
  double resAveRank {0};
  double cplAveRank {0};
  double prgAveRank {0};
  
  for(const auto&r : Res()) 
    if(!r.inCouple()) {
      nSingRes++;
      if(!r.isMatched()) {
        ++resNotMatched;
      }
      else
        resAveRank += r.rankOf(r.matchedTo());
      
      if(r.rankOf(r.matchedTo()) == 0)
        ++resGotTopRank;
      
    }
    
    for(const auto&c : Cpl()) {
      if(!c.isMatched())
        ++cplNotMatched;
      else
        cplAveRank += c.rankOf(c.matchedTo());      
      
      if(c.rankOf(c.matchedTo()) == 0)
        ++cplGotTopRank;
    }
    
    int matchedProgs {0};
    for(const auto&p : Prog()) {
      progSpareCap += p.quota() - p.accepted().size();
      double aveRank {0};
      for(const auto &res : p.accepted()) {
        aveRank += p.rankOf(res);
        if(p.rankOf(res) == 0)
          ++prgGotTopRank;
      }
      if(p.accepted().size() > 0) {
        prgAveRank += aveRank/p.accepted().size();
        ++matchedProgs;
      }
    }
  
  return Rcpp::List::create(Rcpp::Named("unmatchedSingles", resNotMatched),
                            Rcpp::Named("unmatchedCouples", cplNotMatched),
                            Rcpp::Named("unmatchedProgrammSlots", progSpareCap),
                            Rcpp::Named("aveRRank", resAveRank/(nSingRes-resNotMatched)),
                            Rcpp::Named("aveCRank", cplAveRank/(Cpl().size()-cplNotMatched)),
                            Rcpp::Named("avePRank", prgAveRank/matchedProgs),
                            Rcpp::Named("numRTop", resGotTopRank),
                            Rcpp::Named("numCTop", cplGotTopRank),
                            Rcpp::Named("numPTop", prgGotTopRank));
}

//IO
void Problem::printMatch(bool fndMatch, bool outputmatch) {
  if(!fndMatch) {
   Rcpp::Rcout << "#Match not found (resource bounds exceeded)\n";
  }
  int resNotMatched {0};
  int cplNotMatched {0};
  int progSpareCap {0};
  int nSingRes {0};
  int resGotTopRank {0};
  int cplGotTopRank {0};
  int prgGotTopRank {0};

  double resAveRank {0};
  double cplAveRank {0};
  double prgAveRank {0};

 Rcpp::Rcout << "#Single Residents:\n";
  
  for(const auto&r : Res()) 
    if(!r.inCouple()) {
     Rcpp::Rcout << "#Resident " << r.id() << ": ";
     Rcpp::Rcout << "match = " << r.matchedTo() << " Res Ranking = "
	   << r.rankOf(r.matchedTo()) << "/" << r.rol().size() << "\n";
      nSingRes++;
      if(!r.isMatched()) {
	++resNotMatched;
      }
      else
	resAveRank += r.rankOf(r.matchedTo());

      if(r.rankOf(r.matchedTo()) == 0)
	++resGotTopRank;

    }

 Rcpp::Rcout << "\n#Couples:\n";
  for(const auto&c : Cpl()) {
   Rcpp::Rcout << "#Couple " << c.id() << ": ";
   Rcpp::Rcout << "match = " << c.matchedTo() << " Cpl Ranking = "
	 << c.rankOf(c.matchedTo()) << "/" << c.rol().size() << "\n";
    if(!c.isMatched())
      ++cplNotMatched;
    else
      cplAveRank += c.rankOf(c.matchedTo());      

    if(c.rankOf(c.matchedTo()) == 0)
      ++cplGotTopRank;
  }

  int matchedProgs {0};
 Rcpp::Rcout << "\n#Programs:\n";
  for(const auto&p : Prog()) {
   Rcpp::Rcout << "#Program " << p.id() << ": ";
   Rcpp::Rcout << "spares = " << p.quota() - p.accepted().size();
    progSpareCap += p.quota() - p.accepted().size();
   Rcpp::Rcout << " accepted = [";
    for(const auto &res : p.accepted()) 
     Rcpp::Rcout << res << " ";
   Rcpp::Rcout << "] " << " Prog rankings = [";
    double aveRank {0};
    for(const auto &res : p.accepted()) {
     Rcpp::Rcout << p.rankOf(res) << " ";
      aveRank += p.rankOf(res);
      if(p.rankOf(res) == 0)
	++prgGotTopRank;
    }
   Rcpp::Rcout << "]" << "/" << p.rol().size();
    if(p.accepted().size() > 0) {
      prgAveRank += aveRank/p.accepted().size();
      ++matchedProgs;
     Rcpp::Rcout << " ave Prog rank (accepted) = " << aveRank/p.accepted().size();
    }
   Rcpp::Rcout << "\n";
  }

 Rcpp::Rcout << "\n#Matching Summary Stats:\n";
 Rcpp::Rcout << "#Unmatched Singles: " << resNotMatched << "\n";
 Rcpp::Rcout << "#Unmatched Couples: " << cplNotMatched << "\n";
 Rcpp::Rcout << "#Unmatched Program slots: " << progSpareCap << "\n";

  if(nSingRes-resNotMatched > 0)
   Rcpp::Rcout << "#Ave Resident Rank of their matching = "
	 << resAveRank/(nSingRes-resNotMatched) << "\n";
 Rcpp::Rcout << "#Num Residents getting their top rank = "
       << resGotTopRank << "\n";

  if(Cpl().size() - cplNotMatched > 0)
   Rcpp::Rcout << "#Ave Couple Rank of their matching = "
	 << cplAveRank/(Cpl().size()-cplNotMatched) << "\n";
 Rcpp::Rcout << "#Num Couples getting their top rank = "
       << cplGotTopRank << "\n";
    
  if(matchedProgs > 0) 
   Rcpp::Rcout << "#Ave Program Rank of their matched residents "
	 << prgAveRank/matchedProgs << "\n";
 Rcpp::Rcout << "#Num Programs getting their top rank = "
       << prgGotTopRank;
 Rcpp::Rcout << "\n";

  //Output match spec for verification.
  if(outputmatch) {
    if(!fndMatch) 
     Rcpp::Rcout << "m 0\n";
    else{
     Rcpp::Rcout << "m 1\n";
      for(auto& r : Res()) 
     Rcpp::Rcout << "r " << r.id () << " " << r.matchedTo() << "\n";
    }
  }
}

ostream& operator<<(ostream& os, const Problem& prob) {
  os << "#Problem Spec\n#Residents:\n";
  for(auto& res : prob.residents) 
    os << res;
  os << "\n#Couples:\n";
  for(auto& cpl : prob.couples)
    os << cpl;
  os << "\n#Programs:\n";
  for(auto& prog: prob.programs)
    os << prog;
  os << "\n";
  return os;
}