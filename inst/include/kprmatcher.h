/***********[kprmatcher.h]

-------

Class header for "kprmatcher" class.  This class implements the kpr
(Kojima, Pathak & Roth appendix B.2) matching algorithm. It alters
updates the matching stored in a problem object.

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
#ifndef KPRMATCHER_h
#define KPRMATCHER_h

#include "damatcher.h"
#include "problem.h"
#include "rcp.h"
#include <deque>
using std::deque;

class KPRmatcher : public DAmatcher {
public:
 KPRmatcher(): maxNapps {0} {}
  virtual ~KPRmatcher() {}
  
protected:
  bool match_(Problem&);
  void initData(Problem &);
  void processBumped(Rid);
  void matchSingles(bool frmLast);
  void matchCouples(bool frmLast);
  bool chkMatch(Problem &);
  void unmatch(Rid);


  //couples are dealt with indexing the data using the first resident of the couple

  deque<Rid> resProcessQ;
  deque<Rid> cplProcessQ; 
  vector<int> resNxtP;         //resident has applied to all programs in their ROL up to but not including resNxtP[r]
  vector<vector<int>> resNApps; //Number of proposals to each program 
  int maxNapps;
};

#endif
