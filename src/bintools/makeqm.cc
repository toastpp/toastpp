// makeqm		Martin Schweiger		27.9.96
// generates a simple QM (source/detector) description file for
// use with toast.
// It first places the specified number of sources and detectors
// on a circle at equal angular spacing, then moves them radially
// until they sit on the boundary of the supplied mesh.
// The mesh should be centered at the origin. It is not
// necessarily circular, but should not have any concave boundary
// segments.

#include <mathlib.h>
#include <felib.h>
#include <iostream>
#include <fstream>

using namespace std;

int main (void)
{
    char meshname[256], qmname[256];
    int nqnm, i, j, cmd;
    int proftype, qproftype, mproftype;
    bool outline_from_mesh;
    double angle, rad, q0_offs, mesh_size;
    double width, sup, qwidth, qsup, mwidth, msup;
    Mesh mesh;
    Point zpt(2), opt(2), npt(2);
    ifstream ifs;

    cout << "(1) Load mesh for outline information" << endl;
    cout << "(2) Assume circular outline" << endl;
    cin >> cmd;

    switch (cmd) {
    case 1:
	cout << "Mesh file: ";
	cin >> meshname;
	ifs.open (meshname);
	if (!ifs) {
	    cerr << "makeqm: Could not open " << meshname << endl;
	    return 1;
	}
	ifs >> mesh;
	mesh.Setup();
	mesh_size = mesh.Size (&zpt);
	rad = 2.0 * mesh_size;
	outline_from_mesh = true;
	break;
    case 2:
	cout << "Object radius: ";
	cin >> rad;
	outline_from_mesh = false;
	break;
    default:
	xERROR("Invalid command.");
    }


    cout << "Output QM file: ";
    cin >> qmname;

    cout << "Number of sources (= number of detectors): ";
    cin >> nqnm;
    if (nqnm < 0) {
	cerr << "makeqm: Invalid argument" << endl;
	return 1;
    }

    cout << "Angular offset of first source [deg.]: ";
    cin >> q0_offs;
    q0_offs *= Pi/180.0;

    for (int which = 0; which < 2; which++) {
        if (which) cout << "MEASUREMENT profiles:" << endl;
	else       cout << "SOURCE profiles" << endl;
	cout << "\t(1) Point" << endl;
	cout << "\t(2) Gaussian" << endl;
	cout << "\t(3) Cosine" << endl;
	cout << "\t(4) Tophat" << endl;
	cin >> proftype;
	switch (proftype) {
	case 1:
	    break;
	case 2:
	    cout << "1/e half-width of Gaussian: ";
	    cin  >> width;
	    break;
	case 3:
	    cout << "Profile width parameter w [Q(x) = cos(pi/2 x/w)]: ";
	    cin  >> width;
	    cout << "Support radius a [0: a=w]: ";
	    cin  >> sup;
	    if (!sup) sup = width;
	    break;
	case 4:
	    cout << "Half-width of top-hat profile: ";
	    cin  >> width;
	    break;
	default:
	    xERROR("Invalid command");
	    break;
	}
	if (which) {
	    mproftype = proftype;
	    mwidth = width;
	    msup   = sup;
	} else {
	    qproftype = proftype;
	    qwidth = width;
	    qsup   = sup;
	}
    }

    ofstream ofs (qmname);
    ofs << "QM file" << endl;
    ofs << "Dimension 2" << endl;
    ofs << endl;
    ofs.precision (8);

    ofs << "SourceList " << nqnm << endl;
    for (i = 0; i < nqnm; i++) {
	angle = (double)i / (double)nqnm * 2.0 * Pi + q0_offs;
	opt[0] = rad * cos (angle);
	opt[1] = rad * sin (angle);
	if (outline_from_mesh) 
	    npt = mesh.BndIntersect (zpt, opt);
	else
	    npt = opt;
	for (j = 0; j < npt.Dim(); j++) ofs << npt[j] << ' ';
	switch (qproftype) {
	case 1:
	    ofs << "Point" << endl;
	    break;
	case 2:
	    ofs << "Gaussian " << qwidth << endl;
	    break;
	case 3:
	    ofs << "Cosine " << qwidth << ' ' << qsup << endl;
	    break;
	case 4:
	    ofs << "Tophat " << qwidth << endl;
	    break;
	}
    }
    ofs << endl;

    ofs << "MeasurementList " << nqnm << endl;
    for (i = 0; i < nqnm; i++) {
	angle = (i+0.5) / (double)nqnm * 2.0 * Pi + q0_offs;
	opt[0] = rad * cos (angle);
	opt[1] = rad * sin (angle);
	if (outline_from_mesh)
	    npt = mesh.BndIntersect (zpt, opt);
	else
	    npt = opt;
	for (j = 0; j < npt.Dim(); j++) ofs << npt[j] << ' ';
	switch (mproftype) {
	case 1:
	    ofs << "Point" << endl;
	    break;
	case 2:
	    ofs << "Gaussian " << mwidth << endl;
	    break;
	case 3:
	    ofs << "Cosine " << mwidth << ' ' << msup << endl;
	    break;
	case 4:
	    ofs << "Tophat " << mwidth << endl;
	    break;
	}
    }
    ofs << endl;

    ofs << "LinkList" << endl;
    for (i = 0; i < nqnm; i++) {
	ofs << nqnm << ':';
	for (j = 0; j < nqnm; j++) ofs << ' ' << j;
	ofs << endl;
    }
    ofs << endl;
      
    return 0;
}
