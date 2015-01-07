// Take a mesh, a region index file, and an element-wise partial region
// file, and calculate the total size of each region, taking into account
// partial element volumes.

#include "mathlib.h"
#include "felib.h"
#include <iostream>
#include <fstream>

using namespace std;
void WriteEimHeader (const char *meshname, int imgsize, const char *eimname, 
    const char *type);
void WriteEim (const RVector &eim, int imgno, const char *eimname);

struct MatList {      // material properties
    double E;         // Young's modulus
    double nu;        // Poisson's ratio
    double dns;       // density
    double te;        // thermal expansion coefficient
};

int main (void) {
    char meshname[256], fname[256], cbuf[256];
    int i, j, el, nreg, elen, nlen;
    int *matidx;
    double perc, size;
    MatList *mat;
    Mesh mesh;

    cout << "Mesh file name: ";
    cin >> meshname;
    ifstream ifs (meshname);
    ifs >> mesh;
    mesh.Setup();
    elen = mesh.elen();
    nlen = mesh.nlen();
    ifs.close();

    cout << "Partial volume file name: ";
    cin >> fname;
    ifs.open (fname);
    ifs.getline (cbuf, 256); // header
    ifs >> i;
    if (i != elen) {
	cerr << "Partial volume list: invalid length" << endl;
	exit (1);
    }
    cout << "Number of regions: ";
    cin >> nreg;
    RVector regsize(nreg);
    RVector regeim(elen);

    for (i = 0; i < elen; i++) {
	ifs >> el;
	size = mesh.ElSize(el);
	cout << "element " << el << ", size = " << size << endl;
	for (j = 0; j < nreg; j++) {
	    ifs >> perc;
	    perc *= 0.01;
	    regsize[j] += size*perc;
	    regeim[el] += (j+1)*perc;
	}
    }

    cout << "Region sizes:" << endl;
    for (i = 0; i < nreg; i++)
	cout << i << '\t' << regsize[i] << endl;

    WriteEimHeader (meshname, elen, "region.eim", "N/A");
    WriteEim (regeim, 0, "region.eim");
    cout << "Region image written to region.eim" << endl;

    return 0;
}


void WriteEimHeader (const char *meshname, int imgsize, const char *eimname, 
    const char *type)
{
    ofstream ofs (eimname);
    ofs << "EIM" << endl;
    ofs << "Mesh = " << meshname << endl;
    ofs << "SolutionType = " << type << endl;
    ofs << "ImageSize = " << imgsize << endl;
    ofs << "EndHeader" << endl;
}

void WriteEim (const RVector &eim, int imgno, const char *eimname)
{
    ofstream ofs (eimname, ios::app);
    ofs << "Image " << imgno << endl;
    for (int i = 0; i < eim.Dim(); i++)
        ofs << eim[i] << ' ';
    ofs << endl;
}

