/***********[rpmatcher.h]

-------

Class header for "rpmatcher" class. 
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
#ifndef RPMATCHER_h
#define RPMATCHER_h

#include "damatcher.h"
#include "problem.h"
#include "rcp.h"

class RPmatcher : public DAmatcher {
public:
  RPmatcher() {}
  virtual ~RPmatcher() {}

protected:
  bool match_(Problem& p);
  bool extendMatch(Rid, vector<Rid>&);
  void processResident(Rid, vector<Rid>&);
  void processCouple(Cid, vector<Rid>&);
  void processBumped(Rid, vector<Rid>&);
  bool chkMatch(vector<Rid>& rStack, vector<Rid>& matched);
  void unmatch(Rid);
};

#endif
