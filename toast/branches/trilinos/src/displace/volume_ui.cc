#include <iostream>
#include <fstream>
//#include <stdlib.h>
#include "mathlib.h"
#include "felib.h"

const char *WS = " \t";
const bool add_displacements = true;

using namespace std;
//RVector ReadDisp (char *disp_name, int nnodes);

void usage()
{
  cerr << "Usage: displace_ui [in_mesh.msh] [out_mesh] [disp.msh] \n";
  cerr << endl;
  cerr << "Computing volume of each tetrahedra before and after FEM solver." << endl;
  cerr << endl;
  exit(1);
}

int main (int argc, char **argv) {

    char *mesh_in_name = NULL;
    char *mesh_out_name = NULL;
    char *disp_name = NULL;
    char mesh_out_name_vol[255];

    double tol = 1e-12;
    char bname[256], dname[256], cbuf[256], *s;
    double pratio0, modulus0, te0, disp, sz, tot_sz;
    Mesh mesh;
    int nmat;
    struct MatList { double E, nu, dns, te; } *matlist;
    int ex, ey, ez;
    double dx, dy, dz;
    bool *egrid;
    int i, j, k, nd, res, cmd, el;
    int *rowptr, *colidx, nzero;
    int *drowptr, *dcolidx, dn;
    int encoding;
    bool bDisp, bBndf, bAuto;
    ifstream disp_f, bndf_f;
    double noise;
    int nblob;
    double sigma, fac1, fac2, h, w;
    double dT = 1.0;
    double vf, angle;
    double dsp, mind = 1e100, maxd = -1e100, mind2 = 1e100, maxd2 = -1e100;

    srand(1234567);

    if (argc < 4){
      usage();
    }
    mesh_in_name = argv[1];
    argc--;
    argv++;
    mesh_out_name = argv[1];
    argc--;
    argv++;

    ifstream mesh_f(mesh_in_name);
    mesh_f >> mesh;

    mesh.Setup ();

    int n = mesh.nlen();
    int elen = mesh.elen();
    int dim = mesh.Dimension();

    // record element volumes
    RVector size0(elen);
    double elsizemin, elsizemax, elsizeavg = 0.0;
    for (i = 0; i < elen; i++) {
	size0[i] = mesh.ElSize (i);
	if (size0[i] < elsizemin || !i) elsizemin = size0[i];
	if (size0[i] > elsizemax || !i) elsizemax = size0[i];
	elsizeavg += size0[i]/(double)elen;
    }

  
    //load displacements
    disp_name = argv[1];
    argc--;
    argv++;
    
    double x[n*3];
    
    ifstream ifs_disp (disp_name);
    
    for (i = 0; i < n; i++)
    {
      for (j = 0; j < 3; j++)
      {
        ifs_disp >> x[i*3+j];
	//cout << x[i*3+j] << " ";
      }
      //cout << endl;
    }
    // Add displacements to mesh
    if (add_displacements) {
        for (i = 0; i < n; i++) {
	    for (j = 0; j < 3; j++)
	        mesh.nlist[i][j] += x[i*3+j];
	}
	mesh.Setup();
    }

    for (i = 0; i < n; i++) {
        //cout << "node " << i << endl;
	dsp = fabs(x[i*3+2]);
	//cout << "dsp = " << dsp << ", ";
	if (dsp < mind) mind = dsp;
	if (dsp > maxd) maxd = dsp;
	dsp = fabs(x[i*3+0]);
	//cout << "dsp = " << dsp << ", ";
	if (dsp < mind) mind = dsp;
	if (dsp > maxd) maxd = dsp;
	dsp = fabs(x[i*3+1]);
	//cout << "dsp = " << dsp << endl;
	if (dsp < mind) mind = dsp;
	if (dsp > maxd) maxd = dsp;
	//cout << "mind = " << mind << ", maxd = " << maxd << endl;
    }
    cout << "Displacement range: " << mind << " to " << maxd << endl;
    
    int tmp = 0;
    
    strcpy(mesh_out_name_vol,mesh_out_name);    
    strcat(mesh_out_name_vol,"_vol.msh");
    ofstream ofs (mesh_out_name_vol);
    ofs << mesh.elen() << endl;
    
    for (el = 0; el < mesh.elen(); el++) {
	Element *pel = mesh.elist[el];
	ofs << size0[el] << " " << pel->Size() << endl;
      }
    cout << "Total mesh size after displacement: " << mesh.FullSize() << endl;

    return 0;
}                                                                              
