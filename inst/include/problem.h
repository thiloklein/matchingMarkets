/***********[problem.h]

-------

Class header for "Problem" class. 
This class can read and store a matching problem with residents, couples, and programs.

To support an API additional functions to storing a initializing
a matching problem can be added. 

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
#ifndef PROBLEM_h
#define PROBLEM_h

#include <vector>
#include <string>
#include <list>
#include <unordered_set>

#include <Rcpp.h>
using namespace Rcpp;

using std::string;
using std::unordered_set;
using std::vector;
using std::ostream;

class Resident;
class Couple;
class Program;

using namespace std;

class Problem {
public:
  Problem();
  //bool readProblem(string filename);  
  bool readProblem(Rcpp::List residents, Rcpp::List couples, Rcpp::List programs);
  string& getError() { return errMsg; }  
  void printMatch(bool haveMatch, bool outputmatch = true);
  //Sven
  Rcpp::List returnMatch(bool haveMatch);
  string summaryMatch(bool haveMatch);
  Rcpp::List getStats(bool haveMatch);

  Program* ithProg(const int id);
  Resident* ithRes(const int id);
  Couple* ithCpl(const int id);

  const vector<Resident>& Res() const { return residents; }
  const vector<Program>& Prog() const { return programs; }
  const vector<Couple>& Cpl() const { return couples; }

private:
  string errMsg;
  bool probOK;

  //Problem IO and error processing
  unordered_set<int> resIDs;
  unordered_set<int> progIDs;
  unordered_set<int> cplIDs;
  vector<int> progsRanked;
  vector<int> resRanked;

  void readResident(string l);
  void readCouple(string l);
  void readProgram(string l);
  //  Post and Error processing
  void postProcess();
  bool chkID(int rid, unordered_set<int>& Ids, string errmsg);
  void furtherInputChecks();
  void clearErrVecs();
  void postError(string msg) {
    errMsg += msg;
    probOK = false;
  }
  bool ok() { return probOK; }

  friend ostream& operator<<(ostream& os, const Problem& prob);
  
  vector<Resident> residents;
  vector<Program> programs;
  vector<Couple> couples;
};

//IO

ostream& operator<<(ostream& os, const Problem& prob);

#endif
