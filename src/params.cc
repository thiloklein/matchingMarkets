/***********[params.cc]

-------

"params" class implementation
This class holds various parameters for the program.

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
#include <Options.h>
#include <limits>
#include <params.h>

#include <Rcpp.h>

using namespace Minisat;

/*
static IntOption opt_verb("MAIN", "verb", "Verbosity level (0=silent, 1=some,"
			  " 2=more, 3=debuging info, 4=more debugging, 5=prepare to be swamped).", 0,
			  IntRange(0,5));

static IntOption opt_cclim("MAIN", "CClim", 
		"Abort of resident matched to same program more than this number"
		" of times, -1 == no limit",
		-1);

static IntOption opt_algo("MAIN", "alg", "Pick matching algorithms: 0 = Roth Peranson 1999, "
                          "1 = Roth Peranson plus random permute of couple ordering "
			  "(same as \"-rand\" option), 2 = Kojima Pathak & Roth 2013 (Appendix B.2)", 
			  0, IntRange(0,2));

static BoolOption opt_rnd("MAIN", "rnd",
			   "When a match is repeated more than  \"rndCClim\" times"
			   " randomly permute the order we attempt to match couples",
			  false);

static IntOption opt_rnd_cclim("MAIN", "rndCClim",
			 "When \"rnd\" is used we randomize the order of the couples"
			 " after this many matches of a resident into the same program."
			 ". (-1 == limit but this effectively turns off \"rnd\")",
			 200, IntRange(0, std::numeric_limits<int>::max()));

static IntOption opt_kpr_cclim("MAIN", "kprCClim",
			       "Terminate KPR algorithm (alg=2) when a resident applies this many times"
			       " to the same program (-1 = no limit)\n", 100);
*/

Params::Params() {}

void Params::readOptions(int opt_algo) {
  verbosity = 0;
  
  algo = 2;
}

Params params;
