/***********[io.h]

--- 
Generic (Templated) IO functions
---


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
#ifndef IO_h
#define IO_h
#include <utility>
#include <vector>
#include <iostream>
using std::pair;
using std::ostream;
using std::vector;

//Generic Output specializations
template<typename T>
ostream& operator<<(ostream& os, const vector<T>& v) {
  os << "[ ";
  for(const auto& i : v) 
    os << i << " ";
  os << "] (" << v.size() << ")";
  return os;
}

template<typename T>
ostream& operator<<(ostream& os, const pair<T,T>& p) {
  os << "(" << p.first << ", " << p.second << ")";
  return os;
}

#endif
