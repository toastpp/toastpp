// ============================================================================
// opt2                           Martin Schweiger                       3.7.00
// console interface to meshopt routines
// opt2 acts as a filter: reads original mesh from stdin, modifies it according
// to switches, and writes new mesh to stdout
// ============================================================================

#include <iostream>
#include "felib.h"
#include "meshopt.h"

using namespace std;

void DisplayInfo (void);

int main (int argc, char *argv[])
{
    int optmode = 0;
    int i, res, len, *perm;
    Mesh mesh;

    CHECK_EXPIRED ();

    // command line parser
    for (i = 1; i < argc; i++) {
        xASSERT(argv[i][0] == '-', "Error parsing command line");
	switch (argv[i][1]) {
	case 'H':
	    DisplayInfo ();
	    exit (0);
	case 'm':
	    optmode = 0;
	    break;
	case 'b':
	    optmode = 1;
	    break;
	case 't':
	    optmode = 2;
	    break;
	}
    }

    cin >> mesh;
    len = mesh.nlen();
    perm = new int[len];
    for (i = 0; i < len; i++) perm[i] = i;

    switch (optmode) {
    case 0:
        res = Optimise_MMD (mesh, perm, 0, len);
	break;
    case 1:
        res = Optimise_MinBandwidth (mesh, perm, 0, len);
	break;
    case 2:
        res = Optimise_Tinney2 (mesh, perm, 0, len);
	break;
    }

    if (res) {
        cerr << "Optimisation failed. Aborting.\n";
	return 1;
    }

    mesh.Reorder (IVector (len, perm, SHALLOW_COPY));
    delete []perm;
    cout << mesh;
    return 0;
}

void DisplayInfo (void)
{
    cout << "\nProgram \033[1mopt2\033[0m\n";
    cout << VERSION_STRING << endl;
    cout << "Node order optimisation filter for TOAST FEM meshes\n";
    cout << "Reads a mesh from stdin, modifies it according to switches,\n";
    cout << "and writes the new mesh to stdout\n";
    cout << "\nSynopsis:\n";
    cout << "\topt2 -H\n";
    cout << "\topt2 [-m|b|t] < inmesh > outmesh\n";
    cout << "\nCmd line switches:\n";
    cout << "-H\tDisplay this help page and exit\n";
    cout << "-m\tMinimum degree optimisation (default)\n";
    cout << "-b\tMinimum bandwidth optimisation\n";
    cout << "-t\tTinney-2 optimisation\n";
}
