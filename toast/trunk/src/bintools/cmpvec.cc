// compare two toast vectors in files

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include "mathlib.h"

using namespace std;

bool is_silent = false;
const char *vecname[2] = {0,0};

void print_usage_and_exit (int errcode = 0);
void print_error_and_exit (const char *msg, int errcode = 0);

int main (int argc, char *argv[])
{
    int i, n, nvec = 0;
    double maxerr = 0;
    bool is_complex = false;
    bool is_verbose = false;
    bool vec_diff = false;
    bool is_reldiff;

    for (i = 1; i < argc; i++) {
	const char *arg = argv[i];
	if (arg[0] == '-') {
	    switch (arg[1]) {
	    case 'c':
		is_complex = true;
		break;
	    case 'e':
		if (arg[2] != '=' || sscanf(arg+3, "%lf", &maxerr) != 1
		    || maxerr < 0)
		    print_usage_and_exit (4);
		is_reldiff = false;
		break;
	    case 'r':
		if (arg[2] != '=' || sscanf(arg+3, "%lf", &maxerr) != 1
		    || maxerr < 0)
		    print_usage_and_exit (4);
		is_reldiff = true;
		break;
	    case 'h':
		print_usage_and_exit (0);
		break;
	    case 's':
		is_silent = true;
		break;
	    case 'v':
		is_verbose = true;
		break;
	    }
	} else {
	    if (nvec < 2) vecname[nvec++] = arg;
	    else print_error_and_exit ("Syntax error", 4);
	}
    }

    if (nvec != 2) print_error_and_exit ("Syntax error: two files containing toast vectors are required", 4);

    ifstream ifs1(vecname[0]);
    ifstream ifs2(vecname[1]);

    if (!ifs1 || !ifs2) {
        char cbuf[256];
	sprintf (cbuf, "Could not open file: %s", vecname[!ifs1 ? 0 : 1]);
	print_error_and_exit(cbuf, 3);
    }

    if (is_complex) {
	CVector v1, v2;
	ifs1 >> v1;
	ifs2 >> v2;
	if (!v1.Dim() || !v2.Dim())
	    print_error_and_exit (
		 "At least one vector is invalid or zero length", 3);
	if (v1.Dim() != v2.Dim())
	    print_error_and_exit ("Vectors differ in length", 2);
	n = v1.Dim();
	for (i = 0; i < n; i++) {
	    double diff = norm (v1[i]-v2[i]);
	    if (diff && is_reldiff)
	        diff /= max (norm(v1[i]), norm(v2[i]));
	    if (!vec_diff) vec_diff = (diff != 0);
	    if (diff > maxerr) {
		char cbuf[256];
		sprintf (cbuf, "Vectors differ in element %d", i+1);
		print_error_and_exit (cbuf, 1);
	    }
	}
    } else {
	RVector v1, v2;
	ifs1 >> v1;
	ifs2 >> v2;
	if (!v1.Dim() || !v2.Dim())
	    print_error_and_exit (
		 "At least one vector is invalid or zero length", 3);
	if (v1.Dim() != v2.Dim())
	    print_error_and_exit ("Vectors differ in length", 2);
	n = v1.Dim();
	for (i = 0; i < n; i++) {
	    double diff = fabs (v1[i]-v2[i]);
	    if (diff && is_reldiff)
	        diff /= max (fabs(v1[i]), fabs(v2[i]));
	    if (!vec_diff) vec_diff = (diff != 0);
	    if (diff > maxerr) {
		char cbuf[256];
		sprintf (cbuf, "Vectors differ in element %d", i+1);
		print_error_and_exit (cbuf, 1);
	    }
	}
    }

    if (vec_diff)
	print_error_and_exit ("Vectors differ, but are equal within limits", 0);

    return 0;
}

void print_usage_and_exit (int errcode)
{
    cout << "cmpvec: compare two toast vectors\n\n";
    cout << "syntax: cmpvec <file1> <file2> [-c] [-e=<maxerr>] [-r=<maxerr>] [-s] [-v]\n";
    cout << "        cmpvec -h\n";
    cout << "\n<file1> and <file2> must contain a single vector in toast format.\n";
    cout << "\nFlags:\n";
    cout << "c: Compare complex vectors. If not specified, real vectors are assumed.\n"; 
    cout << "e: <maxerr> defines the maximum difference between any two vector elements\n";
    cout << "   to be accepted as equal. If not specified, <maxerr>=0\n is assumed.\n";
    cout << "   This option cannot be combined with -r.\n";
    cout << "r: <maxerr> defines the maximum relative difference between any two vector\n";
    cout << "   elements to be accepted as equal. If not specified, <maxerr>=0\n is\n";
    cout << "   assumed. This option cannot be combined with -e.\n";
    cout << "s: Silent operation, even if vectors are not equal.\n";
    cout << "v: Verbose output\n";
    cout << "\nReturn values:\n";
    cout << "0: vectors are equal within the tolerance limit\n";
    cout << "1: at least one vector element pair is not equal\n";
    cout << "2: vectors are incompatible (not equal length)\n";
    cout << "3: at least one vector file does not contain a toast vector\n";
    cout << "4: sytax error\n";
    exit (errcode);
}

void print_error_and_exit (const char *msg, int errcode)
{
    if (!is_silent) {
        if (vecname[0] && vecname[1])
	    cout << vecname[0] << " and " << vecname[1] << ":\n";
	cout << msg << endl;
	if (errcode)
	    cout << "#### Comparison error! ####\n" << endl;
    }
    exit (errcode);
}
