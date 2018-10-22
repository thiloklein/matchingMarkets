/***********[damatcher.h]

-------

Class header for "matcher" class. 
This class is an abstract base class defining the interface for various implementations of 
Deferred Acceptance Matching Algorithms.

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
#ifndef DAMATCHER_h
#define DAMATCHER_h

#include "problem.h"
#include "System.h"
#include "io.h"

using std::ostream;
using std::cout;
using std::cerr;
using Minisat::cpuTime;
using Minisat::memUsedPeak;

class DAmatcher {
public:
  DAmatcher() : curStatsPrinted {false}, maxRepeatedMatches {0},
  totalRMatches {0}, totalRUnMatches {0}, 
  totalCMatches {0}, totalCUnMatches {0}, rndReorderings {0},
  rpCplsMatched {0}
    {}
  virtual ~DAmatcher() {}

  bool match(Problem& p) {
    cpu_time_start = cpuTime();
    curStatsPrinted = false;
    rndReorderings = 0;
    return match_(p);
  }
    
  void printStatsAndExit(int signum, int exitType) {
    cout << "#Signal Number = " << signum << "\n";
    printStats();
    cout << "#Match not found (interrupted)\n";
    cout << "m 0\n";
    cout.flush();
    cerr.flush();
    _exit(exitType);
  }

  void printStats() {
    if (curStatsPrinted)
      return;
    curStatsPrinted = true;
    double cpu_time_used = cpuTime() - cpu_time_start;
    //mem_used does not track memory used for only the current problem.
    //it tracks peak memory total run of the program---which could
    //involve multiple invocations of the matcher. Perhaps fix later.
    double mem_used = memUsedPeak();
    
    cout << "#CPU = " << cpu_time_used << "\n";
    cout << "#MEM = " << mem_used << "MB\n"; 
    cout << "#Total random reorders of couples " << rndReorderings << "\n";
    cout << "#Max Repeated Matches = " << maxRepeatedMatches << "\n"; 
    cout << "#Total Matches Made = " << totalRMatches + totalCMatches << "\n"; 
    cout << "#Total UnMatches Made = " << totalRUnMatches + totalCUnMatches << "\n"; 
    cout << "#Total Resident Matches Made = " << totalRMatches << "\n"; 
    cout << "#Total Resident UnMatches Made = " << totalRUnMatches << "\n"; 
    cout << "#Total Couple Matches Made = " << totalCMatches << "\n"; 
    cout << "#Total Couple UnMatches Made = " << totalCUnMatches << "\n"; 
    cout << "#RP Couples Matched = " << rpCplsMatched << "\n";
  }

protected:
  virtual bool match_(Problem& p) = 0;
  bool curStatsPrinted;
  double cpu_time_start;
  int maxRepeatedMatches;
  int totalRMatches;
  int totalRUnMatches;
  int totalCMatches;
  int totalCUnMatches;
  int rndReorderings;
  int rpCplsMatched;
};

#endif
