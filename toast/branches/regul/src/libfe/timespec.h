// ==========================================================================
// Module libfe
// File timespec.h
// Declaration of class TimeSpec
// ==========================================================================

#ifndef __TIMESPEC_H
#define __TIMESPEC_H

// ==========================================================================
// class TimeSpec

class TimeSpec {
public:
    int nsteps;			// number of time steps
    int skip;			// interleave factor
    double dtim;		// time step interval length [ps]
    double theta;		// time step coupling factor
    friend std::istream& operator>> (std::istream& i, TimeSpec& tspec);
    friend std::ostream& operator<< (std::ostream& o, TimeSpec& tspec);
};

#endif // !__TIMESPEC_H
