#include <fstream>
#include <stdlib.h>
#include "mathlib.h"
#include "felib.h"

using namespace std;

//#define OUTPUT_SYSMAT
//#define OUTPUT_SOLUTION
//#define OUTPUT_RHS
//#define ENCODE_DISPLACEMENT
//#define READ_VOXEL_MESH
//#define EXPORT_REGION_IMAGE
#define ENCODE_STRAINENERGY

#define SANJAY 1 // debugging. 1 for debug. 0 for no debug.

#define BIGSPRING 1e7
//#define SOLVE_CHOLESKY

const char *WS = " \t";
const bool add_displacements = true;

void WriteNimHeader (char *meshname, int imgsize, char *nimname, char *type);
void WriteNim (const RVector &nim, int imgno, char *nimname);
void WriteEimHeader (char *meshname, int imgsize, char *eimname, char *type);
void WriteEim (const RVector &eim, int imgno, char *eimname);
RVector ReadEim (char *eimname);

int main (int argc, char **argv) {
    double tol = 1e-12;
    char bname[256], dname[256], cbuf[256], *s;
    double pratio0, modulus0, te0, disp, sz, tot_sz;
    Mesh mesh;
    char *mesh_in_name = NULL;
    char *mesh_out_name = NULL;
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

    srand(1234567);

    cout << "Loading " << argv[1] << endl;  
    mesh_in_name = argv[1];

    argc--;
    argv++;

    //cout << "Loading " << argv[1] << endl;  
    mesh_out_name = argv[1];

    //cout << "Mesh file name: ";
    //cin >> bname;
    //ifstream mesh_f(bname);
    //mesh_f >> mesh;
    ifstream mesh_f(mesh_in_name);
    mesh_f >> mesh;

    mesh.Setup ();

    IVector matidx(mesh.elen());
    for (i=0; i < mesh.elen(); i++)
    	matidx[i]=0;
    
    int n = mesh.nlen();
    int elen = mesh.elen();
    int dim = mesh.Dimension();
    int it = 0;

    // record element volumes
    RVector size0(elen);
    double elsizemin, elsizemax, elsizeavg = 0.0;
    for (i = 0; i < elen; i++) {
	size0[i] = mesh.ElSize (i);
	if (size0[i] < elsizemin || !i) elsizemin = size0[i];
	if (size0[i] > elsizemax || !i) elsizemax = size0[i];
	elsizeavg += size0[i]/(double)elen;
    }

    cout << "Total mesh size before displacement: " << mesh.FullSize() << endl;
    cout << "Element sizes: ";
    cout << "min = " << elsizemin
	 << ", max = " << elsizemax
	 << ", avg = " << elsizeavg << endl;

    //Mesh mesh2(mesh);
    while ((elsizemin <= 0.0) && (it < 10))
      {
	cout << "Warning: negative element sizes:" << endl;
	cout << "iteration " << it << endl;
// 	for (i = 0; i < elen; i++)
// 	  {
// 	    if (size0[i] <= 0.0) 
// 	      {
// 		cout << "el=" << i << " (region " << matidx[i]
// 		     << "), size=" << size0[i]
// 		     << ", nodes:";
// 		for (j = 0; j < mesh.elist[i]->nNode(); j++)
// 		  cout << " " << mesh.elist[i]->Node[j];
// 		cout << endl;
// 	      }
// 	  }
	
	// try to relax the mesh around bad elements
	//Mesh mesh2(mesh);
	int nbadnd = 0, nrelaxed = 0;
	for (i = 0; i < elen; i++)
	  {
	    if (size0[i] <= 0.0) nbadnd += mesh.elist[i]->nNode();
	  }
	int *badnd = new int[nbadnd];
	for (i = nbadnd = 0; i < elen; i++)
	  {
	    if (size0[i] <= 0.0)
	      {
		for (j = 0; j < mesh.elist[i]->nNode(); j++)
		  badnd[nbadnd++] = mesh.elist[i]->Node[j];
	      }
	  }
	
	for (i = 0; i < nbadnd; i++) 
	  {
	    if (badnd[i] < 0) continue;
	    //cout << "Relaxing node " << badnd[i] << endl;
	    Point bc = mesh.NeighbourBarycentre(badnd[i]);
	    for (k = 0; k < bc.Dim(); k++)
	      mesh.nlist[badnd[i]][k] = bc[k];
	    for (j = i+1; j < nbadnd; j++)
	      if (badnd[j] == badnd[i]) badnd[j] = -1;
	    nrelaxed++;
	  }
	delete []badnd;
	mesh.Setup();
	it++;
	//cout << "iteration incremental, it = " << it << endl;
	
	// record element volumes
	//RVector size1(elen);
	elsizemin = 1000.0;
	elsizemax = 0.0;
	elsizeavg = 0.0;
	for (i = 0; i < elen; i++) 
	  {
	    size0[i] = mesh.ElSize (i);
	    if (size0[i] < elsizemin || !i) elsizemin = size0[i];
	    if (size0[i] > elsizemax || !i) elsizemax = size0[i];
	    elsizeavg += size0[i]/(double)elen;
	  }
	cout << "elsizemin = " << elsizemin << endl;
	
	
      }//for negative elements
    
    //delete []badnd;
    //mesh2.Setup();
    //ofstream ofs ("relaxed.msh");
    ofstream ofs (mesh_out_name);
    ofs << mesh << endl;
    cout << "Relaxed mesh written to " << mesh_out_name << endl;
    exit (0);

 

    return 0;
}                                                                              

void WriteNimHeader (char *meshname, int imgsize, char *nimname, char *type)
{
    ofstream ofs (nimname);
    ofs << "NIM" << endl;
    ofs << "Mesh = " << meshname << endl;
    ofs << "SolutionType = " << type << endl;
    ofs << "ImageSize = " << imgsize << endl;
    ofs << "EndHeader" << endl;
}

void WriteNim (const RVector &nim, int imgno, char *nimname)
{
    ofstream ofs (nimname, ios::app);
    ofs << "Image " << imgno << endl;
    for (int i = 0; i < nim.Dim(); i++)
        ofs << nim[i] << ' ';
    ofs << endl;
}

void WriteEimHeader (char *meshname, int imgsize, char *eimname, char *type)
{
    ofstream ofs (eimname);
    ofs << "EIM" << endl;
    ofs << "Mesh = " << meshname << endl;
    ofs << "SolutionType = " << type << endl;
    ofs << "ImageSize = " << imgsize << endl;
    ofs << "EndHeader" << endl;
}

void WriteEim (const RVector &eim, int imgno, char *eimname)
{
    ofstream ofs (eimname, ios::app);
    ofs << "Image " << imgno << endl;
    for (int i = 0; i < eim.Dim(); i++)
        ofs << eim[i] << ' ';
    ofs << endl;
}

RVector ReadEim (char *eimname)
{
    char cbuf[256];
    int i, imgsize = 0;
    ifstream ifs (eimname);
    ifs.getline (cbuf, 256);
    if (strcmp (cbuf, "EIM"))
	return  RVector(); // problem
    while (ifs.getline (cbuf, 256)) {
	if (!strcasecmp (cbuf, "EndHeader")) {
	    break;
	} else if (!strncasecmp (cbuf, "ImageSize", 9)) {
	    sscanf (cbuf+11, "%d", &imgsize);
	}
    }
    RVector eim(imgsize);
    do {
	ifs.getline (cbuf, 256);
    } while (strncasecmp (cbuf, "Image", 5));
    for (i = 0; i < imgsize; i++)
	ifs >> eim[i];
    return eim;
}
    
