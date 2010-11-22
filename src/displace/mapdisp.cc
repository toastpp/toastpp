// This maps a displacement field, given by a base mesh and a distorted
// mesh with conforming structure, onto a third non-conforming mesh
// Output is a mesh conforming with the third mesh, distorted by the
// calculated displacement field.

#include <fstream.h>
#include <iostream.h>
#include "mathlib.h"
#include "felib.h"

int main (void)
{
    char fname[256];
    ifstream ifs;
    int i, j, n, el;
    Mesh mesh1_orig, mesh1_disp, mesh2_orig, mesh2_disp;

    cout << "Step 1: calculation of displacement field" << endl;
    cout << "Give the names of the original and the conforming distorted\n";
    cout << "mesh from which to calculate the displacements" << endl;
    cout << "Original mesh:  ";
    cin >> fname;
    ifs.open (fname);
    ifs >> mesh1_orig;
    ifs.close();

    cout << "Distorted mesh: ";
    cin >> fname;
    ifs.open (fname);
    ifs >> mesh1_disp;
    ifs.close();

    cout << "\nStep 2: mapping onto non-conforming mesh" << endl;
    cout << "Give the mesh name of the original non-conforming mesh onto\n";
    cout << "which to map the displacement field." << endl;
    cout << "Mesh name: ";
    cin >> fname;
    ifs.open (fname);
    ifs >> mesh2_orig;
    ifs.close();

    cout << "\nSetting up meshes" << endl;
    mesh1_orig.Setup();
    mesh1_disp.Setup();
    mesh2_orig.Setup();

    // now do the mapping
    mesh2_disp.Copy(mesh2_orig);
    n = mesh2_orig.nlen();

    for (i = 0; i < n; i++) {
	Point &pt = mesh2_orig.nlist[i];
	el = mesh1_orig.ElFind (pt);
	if (el < 0) { // problem here!
	    double d, mindist = 1e10;
	    for (j = 0; j < mesh1_orig.elen(); j++) {
		d = pt.Dist (mesh1_orig.ElCentre(j));
		if (d < mindist) el = j, mindist = d;
	    }
	}
	Element *pel = mesh1_orig.elist[el];
	int nnode = pel->nNode();
	int *node = pel->Node;
	RVector fun = pel->GlobalShapeF (mesh1_orig.nlist, pt);

	RVector disp(fun.Dim());
	for (j = 0; j < pel->nNode(); j++) {
	    RVector dn = mesh1_disp.nlist[node[j]]-mesh1_orig.nlist[node[j]];
	    disp += dn*fun[j];
	    mesh2_disp.nlist[i] += disp;
	}
    }
    
    ofstream ofs("dispmap.msh");
    ofs << mesh2_disp;
    cout << "Displacement mapped and applied to mesh.\n";
    cout << "Distorted mesh written to dispmap.msh" << endl;
    return 0;
}
