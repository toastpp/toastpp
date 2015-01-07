// ==========================================================================
// Module libfe
// File ptsource.cc
// Definition of class PointSource
// ==========================================================================

#define FELIB_IMPLEMENTATION

#include <iostream>
#include <sstream>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "mathlib.h"
#include "felib.h"

using namespace std;

PointSource::PointSource ()
: Point ()
{
    strength = 0;
}

PointSource::PointSource (const PointSource &ps)
: Point (ps)
{
    strength = ps.strength;
}

ostream &operator<< (ostream &os, const PointSource &ps)
{
    os << "PointSource " << ps.strength << endl;
    for (int i = 0; i < ps.Dim(); i++) os << ps[i] << ' ';
    os << endl;
    return os;
}

istream &operator>> (istream &is, PointSource &ps)
{
    char cbuf[200];
    do {
	is.getline (cbuf, 200);
    } while (strncmp (cbuf, "PointSource", 11));
    std::istringstream iss (cbuf+11);
    iss >> ps.strength;
    for (int i = 0; i < ps.Dim(); i++) is >> ps[i];
    return is;
}
