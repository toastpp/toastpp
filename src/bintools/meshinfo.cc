// meshinfo.cc
// Prints out information about a FEM mesh.
//
// Syntax: meshinfo meshfile

#include <mathlib.h>
#include <felib.h>
#include <fstream>

using namespace std;

void PrintUsage();
void ElSizeRange (const Mesh &mesh, double &emin, double &emax);

int main (int argc, char *argv[])
{
    char *meshfile;
    Mesh mesh;
    int i, j, cmd;
    double elsizemin, elsizemax;

    if (argc < 2) {
	PrintUsage ();
	exit (1);
    }
    meshfile = argv[1];

    cout << "Reading " << meshfile << " ..." << flush;
    ifstream ifs (meshfile);
    if (!ifs) {
	cerr << " FAILED!" << endl;
	exit (1);
    }
    ifs >> mesh;
    mesh.Setup();
    int nlen = mesh.nlen();
    int elen = mesh.elen();

    cout << " finished" << endl;
    cout << "Number of nodes: " << mesh.nlen() << endl;
    cout << "Number of elements: " << mesh.elen() << endl;
    cout << "Total semi-bandwidth: " << (mesh.MaxNodeDiff (BW_TOTAL) + 1)
	 << endl;
    cout << "Internal semi-bandwidth: " << (mesh.MaxNodeDiff (BW_INTERNAL) + 1)
	<< endl;

    // check element sizes
    cout << "Mesh size: " << mesh.FullSize() << endl;
    ElSizeRange (mesh, elsizemin, elsizemax);
    cout << "Element size range: " << elsizemin << " to " << elsizemax << endl;

    Point bbmin, bbmax;
    mesh.BoundingBox (bbmin, bbmax);
    cout << "Mesh bounding box: " << bbmin << " x " << bbmax << endl;

    if (elsizemin < 0) {
	cout << "Swap nodes in elements with negative sizes (1|0)? ";
	cin >> cmd;
	if (cmd) {
	    for (i = 0; i < elen; i++) {
		double elsize = mesh.ElSize (i);
		if (elsize < 0) {
		    cout << "Swapping nodes in element " << i << endl
			 << "\033[1A";
		    Element *pel = mesh.elist[i];
		    switch (pel->Type()) {
		    case ELID_TRI3:
		    case ELID_TRI3OLD:
			j = pel->Node[1];
			pel->Node[1] = pel->Node[2];
			pel->Node[2] = j;
			break;
		    case ELID_TET4:
			j = pel->Node[1];
			pel->Node[1] = pel->Node[2];
			pel->Node[2] = j;
			break;
		    default:
			cerr << "Element type not supported\n";
			exit (1);
		    }
		}
	    }
	    mesh.Setup();
	    ElSizeRange (mesh, elsizemin, elsizemax);
	    cout << "New element size range: " << elsizemin << " to "
		 << elsizemax << endl;
	    ofstream ofs ("meshinfo.msh");
	    ofs << mesh << endl;
	    cout << "Modified mesh written to meshinfo.msh" << endl;
	}
    }
}

void PrintUsage()
{
    cerr << "Usage: meshinfo meshfile" << endl;
}

void ElSizeRange (const Mesh &mesh, double &emin, double &emax)
{
    double elsize;
    for (int i = 0; i < mesh.elen(); i++) {
	elsize = mesh.ElSize (i);
	if (!i || elsize < emin) emin = elsize;
	if (!i || elsize > emax) emax = elsize;
    }
}
