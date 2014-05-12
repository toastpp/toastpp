// ==========================================================================
// Module libfe
// File timespec.cc
// Definition of class TimeSpec
// ==========================================================================

#define FELIB_IMPLEMENTATION

#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "mathlib.h"
#include "felib.h"
#include "timespec.h"

using namespace std;

istream& operator>> (istream& i, TimeSpec& tspec)
{
    char cbuf[200];

    while (i >> cbuf && strncmp (cbuf, "TimeSpec", 8));
    if (!i) return i;
    i >> tspec.nsteps >> tspec.dtim >> tspec.theta >> tspec.skip;
    return i;
}

ostream& operator<< (ostream& o, TimeSpec& tspec)
{
    o << "TimeSpec " << tspec.nsteps << " " << tspec.dtim << " "
	<< tspec.theta << " " << tspec.skip << endl;
    return o;
}
