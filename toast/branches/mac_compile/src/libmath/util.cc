// ============================================================================
// TOAST                                              (c) Martin Schweiger 2000
// Library: libmath
// File:    util
//
// General-purpose utility routines
// ============================================================================

#define MATHLIB_IMPLEMENTATION

#include <ctype.h>
#include <stdlib.h>
#include "mathlib.h"

using namespace std;

// ============================================================================
// Converts 0-terminated string 'str' to uppercase
// Returns str

char *ToupperString (char *str)
{
    for (char *s = str; *s; s++) *s = toupper (*s);
    return str;
}

// ============================================================================
// This takes 0-terminated input string 'str' of the form
//
// <string-a> = <string-b>
//
// and returns pointers to the beginnings of <string-a> and <string-b> in
// 'cat' and 'val', respectively. Both <string-a> and <string-b> are
// 0-terminated on exit, and have leading and trailing whitespace removed.
// The original string is modified.

void SplitEqString (char *str, char **cat, char **val)
{
    char *c;

    // find beginning of cat
    for (*cat = str; **cat && (**cat == ' ' || **cat == '\t'); (*cat)++);

    // find '='
    for (*val = *cat; **val && **val != '='; (*val)++);

    // cut whitespace before '='
    for (c = *val-1; c >= *cat && (*c == ' ' || *c == '\t'); c--)
        *c = '\0';

    // find beginning of val
    if (**val) {
        *((*val)++) = '\0';
	for (; **val && (**val == ' ' || **val == '\t'); (*val)++);

	// cut trailing whitespace
	for (c = *val; *c; c++);
	while (--c >= *val && (*c == ' ' || *c == '\t')) *c = '\0';
    }
}

// ============================================================================
// Write an image (stored in 1-D data array 'img') of dimension xdim,ydim
// to a PPM pixmap file. Image is scaled to range scalemin-scalemax.
// If scalemin and/or scalemax is NULL, then autscaling is used.

void WritePPM (const RVector &img, int xdim, int ydim,
    double *scalemin, double *scalemax, char *fname)
{
    typedef struct {
        unsigned char r,g,b;
    } RGB;

    int i, dim = xdim*ydim;
    unsigned char *pixmap = new unsigned char[dim];
    double imgmin, imgmax, scale;

    static RGB colmap[256];
    static bool have_colmap = false;

    if (!have_colmap) {
        int r, g, b;
        char colormap[256];
	strcpy (colormap, getenv ("TOASTDIR"));
	strcat (colormap, "/scales/fire2.pal");
        ifstream ifs (colormap);
	for (i = 0; i < 256; i++) {
	    ifs >> r >> g >> b;
	    colmap[i].r = (unsigned char)r;
	    colmap[i].g = (unsigned char)g;
	    colmap[i].b = (unsigned char)b;
	}
	have_colmap = true;
    }

    // rescale image
    if (scalemin) {
        imgmin = *scalemin;
    } else {
	for (i = 0, imgmin = 1e100; i < dim; i++)
	    if (img[i] < imgmin) imgmin = img[i];
    }
    if (scalemax) {
        imgmax = *scalemax;
    } else {
      for (i = 0, imgmax = -1e100; i < dim; i++)
	  if (img[i] > imgmax) imgmax = img[i];
    }
    scale = 256.0/(imgmax-imgmin);

    for (i = 0; i < dim; i++) {
        int v = (int)((img[i]-imgmin)*scale);
	if      (v < 0  ) v = 0;
	else if (v > 255) v = 255;
	pixmap[i] = (unsigned char)v;
    }
    ofstream ofs(fname);
    ofs << "P6" << endl;
    ofs << "# CREATOR: TOAST (test_grid)" << endl;
    ofs << xdim << ' ' << ydim << endl;
    ofs << 255 << endl;
    for (i = 0; i < dim; i++)
        ofs << colmap[pixmap[i]].r
	    << colmap[pixmap[i]].g
	    << colmap[pixmap[i]].b;
    ofs << endl;

    delete []pixmap;
}

