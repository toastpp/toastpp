#include <fstream>
#include <stdlib.h>
#include "mathlib.h"
#include "felib.h"

using namespace std;

#define OUTPUT_SYSMAT 0
#define OUTPUT_SOLUTION 0
#define OUTPUT_RHS 0
//#define ENCODE_DISPLACEMENT
//#define READ_VOXEL_MESH
//#define EXPORT_REGION_IMAGE
#define ENCODE_STRAINENERGY
//#define false = 0
//#define true = 1

#define SANJAY 1 // debugging. 1 for debug. 0 for no debug.
#define BIGSPRING 1e7
//#define SOLVE_CHOLESKY

const char *WS = " \t";
const bool add_displacements = true;
const double DISP_NONE = 1e10; // flag for 'no displacement'

void WriteNimHeader (char *meshname, int imgsize, char *nimname, char *type);
void WriteNim (const RVector &nim, int imgno, char *nimname);
void WriteEimHeader (char *meshname, int imgsize, char *eimname, char *type);
void WriteEim (const RVector &eim, int imgno, char *eimname);
RVector ReadEim (char *eimname);
void usage()
{
  cerr << "Usage: displace_ui [in_mesh.msh] [out_mesh] <options> \n";
  cerr << "where <options> is one or more of the following:\n";
  cerr << "\t <-b boundary_displ_file> File with boundary displacements\n\n";
  cerr << "\t <-e mode>  element biomechanical properties.  :\n";
  cerr << "\t 		       1 (homogeneous) [Poisson's ratio (double) Young's modulus (double) Thermal expansion (double)] | :\n";
  cerr << "\t 		       2 (from file) [file name]  \n\n";
  cerr << "\t <-p perturbation_mode> Apply perturbations to thermal coefficients: \n";
  cerr << "\t 		       0 (No) | 1 (Gaussian element noise) [Noise level (double)] | \n";
  cerr << "\t 		       2 (Random Gaussian blobs) [Number of blobs (int) Blob level (double) Blob width (double)] | \n";
  cerr << "\t 		       3 (Read coefficient file) [coeff_file]:\n";
  cerr << "\t <-K temperature change(double)> Applied temperature change\n\n";
  cerr << "\t <-hvolf Force azimuth angle degrees (double) Force magnitude (double)> Add homogeneous volume forces\n\n";
  cerr << "\t <-disp> Mesh parameter encoding displacements \n\n";
  cerr << "\t <-strain> Mesh parameter encoding strain energy density \n\n";
  cerr << "\t <-elchange> Mesh parameter encoding element volume change \n\n";
  cerr << "\t <-vol> Mesh parameter encoding tetra volume before and after FEM solver \n\n";
  cerr << "\t <-disp_range min (double) max (double)> Displacement encoding range \n\n";
  cerr << endl;
  cerr << "Forward FEM solver." << endl;
  cerr << endl;
  exit(1);
}


int main (int argc, char **argv) {

    char *mesh_in_name = NULL;
    char *mesh_out_name = NULL;
    char *boundary_name = NULL;
    char *el_file_name = NULL;
    char *eim_name = NULL;
    char mesh_out_name_disp[255];
    char mesh_out_name_strain[255];
    char mesh_out_name_elvolchange[255];
    char mesh_out_name_vol[255];
    char mesh_out_name_realvolchange[255];

    double tol = 1e-6;
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

    bool el_homogeneous = false;
    bool el_file = false;
    bool boundary = false;
    bool gaussian_noise = false, gaussian_blobs = false, coefficient_file = false;
    bool hom_vol_forces = false, enc_vol = false, enc_disp = false, enc_strain =  false, enc_elvolchange =  false;
    bool disp_range = false; 
    bool ok =  false;

    srand(1234567);

    if (argc < 3){
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
    
  // Parse remaining arguments
  while (argc > 1) {
    ok = false;
    
    if ((ok == false) && (strcmp(argv[1], "-b") == 0)){
      argc--;
      argv++;
      boundary = true;
      bDisp = true; 
      ok = true;
      boundary_name = argv[1];
      argc--;
      argv++;

    }    
    if ((ok == false) && (strcmp(argv[1], "-e") == 0)){
      argc--;
      argv++;
      if (strcmp(argv[1], "1") == 0)
      {
        el_homogeneous = true;
	el_file = false;
	nmat = 1;
	matlist = new MatList[nmat];

	argc--;
	argv++;
	matlist[0].nu = atof(argv[1]);
	argc--;
	argv++;
	matlist[0].E = atof(argv[1]);
	argc--;
	argv++;
	matlist[0].te = atof(argv[1]);
     }
      else if (strcmp(argv[1], "2") == 0)
      {
	el_homogeneous = false;
	el_file = true;
        argc--;
        argv++;
	el_file_name = argv[1];
      }
      argc--;
      argv++;
      ok = true;
    }    
    if ((ok == false) && (strcmp(argv[1], "-p") == 0)){
      argc--;
      argv++;
      if (strcmp(argv[1], "0") == 0)
      {
        gaussian_noise = false;
	gaussian_blobs = false;
	coefficient_file = false;
      }
      else if (strcmp(argv[1], "1") == 0)
      {
        gaussian_noise = true;
	gaussian_blobs = false;
	coefficient_file = false;
	argc--;
	argv++;
	noise = atof(argv[1]);
      }
      else if (strcmp(argv[1], "2") == 0)
      {
        gaussian_noise = false;
	gaussian_blobs = true;
	coefficient_file = false;
	argc--;
	argv++;
	nblob = atoi(argv[1]);
	argc--;
	argv++;
	noise = atof(argv[1]);
	argc--;
	argv++;
	w = atof(argv[1]);
      }
      else if (strcmp(argv[1], "3") == 0)
      {
        gaussian_noise = false;
	gaussian_blobs = false;
	coefficient_file = true;
	argc--;
	argv++;
	eim_name = argv[1];
     }
      argc--;
      argv++;
      ok = true;
    } 
    if ((ok == false) && (strcmp(argv[1], "-K") == 0)){
      argc--;
      argv++;
      ok = true;
      dT = atof(argv[1]);
      argc--;
      argv++;

    }    
     if ((ok == false) && (strcmp(argv[1], "-hvolf") == 0)){
      argc--;
      argv++;
      hom_vol_forces = true; 
      ok = true;
      angle = atof(argv[1]);
      argc--;
      argv++;
      vf = atof(argv[1]);
      argc--;
      argv++; 
     }    
    if ((ok == false) && (strcmp(argv[1], "-disp") == 0)){
      argc--;
      argv++;
      enc_disp = true; 
      ok = true;
    }    
    if ((ok == false) && (strcmp(argv[1], "-vol") == 0)){
      argc--;
      argv++;
      enc_vol = true; 
      ok = true;
    }    
    if ((ok == false) && (strcmp(argv[1], "-strain") == 0)){
      argc--;
      argv++;
      enc_strain = true; 
      ok = true;
    }    
    if ((ok == false) && (strcmp(argv[1], "-elchange") == 0)){
      argc--;
      argv++;
      enc_elvolchange = true; 
      ok = true;
    }    
    if ((ok == false) && (strcmp(argv[1], "-disp_range") == 0)){
      argc--;
      argv++;
      disp_range = true; 
      ok = true;
      mind2 = atof(argv[1]);
      argc--;
      argv++;
      maxd2 = atof(argv[1]);
      argc--;
      argv++;
    }    
  
    if ( !ok ) {
      usage();
    }
  }
    if (boundary)
        disp_f.open (boundary_name);

    // Generate system stiffness matrix
    mesh.SparseRowStructure (rowptr, colidx, nzero);
    BlockExpand (rowptr, colidx, n, drowptr, dcolidx, dn, dim, dim);
    RCompRowMatrix F (dn, dn, drowptr, dcolidx);
    delete []rowptr;
    delete []colidx;
    delete []drowptr;
    delete []dcolidx;

    RVector rhs(dn);
    RVector x(dn);

    RVector pratio(mesh.elen());
    RVector modulus(mesh.elen());
    RVector te(mesh.elen());
    IVector matidx(mesh.elen());

    if (el_homogeneous)
      matidx = 0;

    if (el_file_name)
      {
        int nidx, idx;
	ifstream ifs (el_file_name);
	ifs >> nmat;
	matlist = new MatList[nmat];
	for (i = 0; i < nmat; i++) {
	  ifs >> matlist[i].E >> matlist[i].nu >> matlist[i].dns
	      >> matlist[i].te;
	  cout << "mat " << i << ", E = " << matlist[i].E << ", nu = " << matlist[i].nu << ", dns = " << matlist[i].dns
	  		<< ", te = " << matlist[i].te << endl;
	}
	ifs >> nidx;
	if (nidx != mesh.elen()) {
	  cerr << "Invalid length of material index list. Aborting." << endl;
	  exit(1);
	}
	for (i = 0; i < nidx; i++) {
	  ifs >> idx;
	  if (idx < 0 || idx >= nmat) {
	    cerr << "Invalid index in material index list. Aborting."
		 << endl;
	    exit (1);
	  }
	  matidx[i] = idx;
	}
      }
    
    // assign material properties
    for (i = 0; i < elen; i++) {
      pratio[i] = matlist[matidx[i]].nu;
      modulus[i] = matlist[matidx[i]].E;
      te[i] = matlist[matidx[i]].te;
    }
    
    //Apply perturbations to thermal coefficients
    if (gaussian_noise)
      {
	for (i = 0; i < elen; i++)
	  te[i] += gasdev (noise);
      }    
    else if (gaussian_blobs)
      {
	Point bbmin, bbmax, cnt(mesh.Dimension());
	mesh.BoundingBox (bbmin, bbmax);
	for (i = 0; i < nblob; i++) {
	  // find blob centre
	  for (j = 0; j < mesh.Dimension(); j++)
	    cnt[j] = bbmin[j] + (bbmax[j]-bbmin[j]) *
	      (double)rand()/(double)RAND_MAX;
	  // find blob sigma
	  sigma = (double)rand()/(double)RAND_MAX * w; sigma = 20.0;
	  // find blob height
	  h = (double)rand()/(double)RAND_MAX * noise;  //h = 10.0;
	  if ((double)rand() > 0.5*(double)RAND_MAX) h = -h;
	  // add blob
	  fac1 = h/sqrt(2.0*M_PI*sigma*sigma);
	  fac2 = -1.0/(2.0*sigma*sigma);
	  for (j = 0; j < elen; j++) {
	    Point elcnt = mesh.ElCentre(j);
	    double dst = elcnt.Dist (cnt);
	    te[j] += fac1 * exp(dst*dst*fac2);
	  }
	}
	
      }
    else if (coefficient_file)
      {
	te = ReadEim(eim_name);
      }
    
#ifdef EXPORT_REGION_IMAGE
    // temporary: map material index to a nim file
    RVector ndmat (n);
    IVector rndval(nmat);
    IVector ncount(n);
    for (i = 0; i < nmat; i++) {
	rndval[i] = i+1;
    }
    for (i = 0; i < elen; i++) {
	Element *pel = mesh.elist[i];
	for (j = 0; j < pel->nNode(); j++) {
	    ndmat[pel->Node[j]] += rndval[matidx[i]];
	    ncount[pel->Node[j]]++;
	}
    }
    for (i = 0; i < n; i++) ndmat[i] /= ncount[i];
    WriteNimHeader (mesh_in_name, n, "material.nim", "N/A");
    WriteNim (ndmat, 1, "material.nim");
#endif 

    cout << "Total mesh size before displacement: " << mesh.FullSize() << endl;
    cout << "Element sizes: ";
    cout << "min = " << elsizemin
	 << ", max = " << elsizemax
	 << ", avg = " << elsizeavg << endl;

    if (elsizemin <= 0.0) {
	cout << "Warning: negative element sizes:" << endl;
	for (i = 0; i < elen; i++)
	    if (size0[i] <= 0.0) {
		cout << "el=" << i << " (region " << matidx[i]
		     << "), size=" << size0[i]
		     << ", nodes:";
		for (j = 0; j < mesh.elist[i]->nNode(); j++)
		    cout << " " << mesh.elist[i]->Node[j];
		cout << endl;
	    }

	// try to relax the mesh around bad elements
	Mesh mesh2(mesh);
	int nbadnd = 0, nrelaxed = 0;
	for (i = 0; i < elen; i++)
	    if (size0[i] <= 0.0) nbadnd += mesh.elist[i]->nNode();
	int *badnd = new int[nbadnd];
	for (i = nbadnd = 0; i < elen; i++)
	    if (size0[i] <= 0.0)
		for (j = 0; j < mesh.elist[i]->nNode(); j++)
		    badnd[nbadnd++] = mesh.elist[i]->Node[j];
#ifdef UNDEF
	// check for nodes appearing more than once
	for (i = 0; i < nbadnd-1; i++)
	    for (j = i+1; j < nbadnd; j++)
		if (badnd[i] == badnd[j] && badnd[j] >= 0) {
		    cout << "Relaxing node " << badnd[i] << endl;
		    Point bc = mesh.NeighbourBarycentre(badnd[i]);
		    for (k = 0; k < bc.Dim(); k++)
			mesh2.nlist[badnd[i]][k] = bc[k];
		    badnd[i] = badnd[j] = -1;
		    nrelaxed++;
		}
#endif
	if (!nrelaxed) { // no duplicates - relax all nodes in list
	    for (i = 0; i < nbadnd; i++) {
		if (badnd[i] < 0) continue;
		cout << "Relaxing node " << badnd[i] << endl;
		Point bc = mesh.NeighbourBarycentre(badnd[i]);
		for (k = 0; k < bc.Dim(); k++)
		    mesh2.nlist[badnd[i]][k] = bc[k];
		for (j = i+1; j < nbadnd; j++)
		    if (badnd[j] == badnd[i]) badnd[j] = -1;
		nrelaxed++;
	    }
	}
	delete []badnd;
	mesh2.Setup();
	ofstream ofs ("relaxed.msh");
	ofs << mesh2 << endl;
	cout << "Relaxed mesh written to relaxed.msh" << endl;
	exit (0);
    }

    cout << "Assembling system matrix" << endl;
    AddToSysMatrix_elasticity (mesh, F, modulus, pratio);
    cout << "Finished!" << endl;

    // add thermal expansion
    AddToRHS_thermal_expansion (mesh, rhs, modulus, pratio, te, dT);

#ifdef UNDEF
    for (i = 0; i < mesh.nlen(); i++) {
      if (mesh.nlist[i][2] == 0) {
        double disp[3] = {0,0,0.01};
        for (j = 0; j < 3; j++) {
          F(i*3+j,i*3+j) *= BIGSPRING;
          rhs[i*3+j] = F(i*3+j,i*3+j)*disp[j];
        }
      } else if (mesh.nlist[i][2] == 3) {
        double disp[3] = {0,0,-0.01};
        for (j = 0; j < 3; j++) {
          F(i*3+j,i*3+j) *= BIGSPRING;
          rhs[i*3+j] = F(i*3+j,i*3+j)*disp[j];
        }
      }
    }
#endif

    // add explicit boundary displacements using Payne-Irons
    // Dirichlet condition
    RDenseMatrix Disp(n, 3);
    Disp = DISP_NONE;

    if (bDisp) {
        cout << "adding explicit boundary displacements" << endl;
	while (disp_f.getline (cbuf, 256)) {
	    s = strtok (cbuf, WS);
	    res = sscanf (s, "%d", &nd);
	    if (!res) continue;
	    for (i = 0; i < 3; i++) {
	        s = strtok (NULL, WS);
		res = sscanf (s, "%lf", &Disp(nd,i));
	    }
	}
    }
    
    for (i = 0; i < n; i++)
      for (j = 0; j < 3; j++)
      { 
        if (Disp(i,j) != DISP_NONE) {
	  F(i*3+j, i*3+j) *= BIGSPRING;
	  rhs[i*3+j] = F(i*3+j, i*3+j) * Disp(i,j);
        }
       }


    // add volume forces
    if (hom_vol_forces) {
	cout << "adding volume forces" << endl;
	RVector fvec(dn);
	RVector dir(dim);
	angle *= M_PI/180.0;
	dir[0] = cos (angle);
	dir[2] = sin (angle);
	for (i = 0; i < n; i++)
	    for (j = 0; j < dim; j++)
	        fvec[i*3+j] = vf*dir[j];
	AddToRHS_elasticity (mesh, rhs, &fvec, RHS_P);
    }

    // add boundary forces
    if (bBndf) {
        cout << "adding boundary forces" << endl;
        double bf;
	RVector fvec(dn);
	while (bndf_f.getline (cbuf, 256)) {
	    s = strtok (cbuf, WS);
	    res = sscanf (s, "%d", &nd);
	    if (!res) continue;
	    if (!mesh.nlist[nd].isBnd()) {
	        cerr << "Can't apply surface force to internal node\n";
		exit(1);
	    }
	    for (i = 0; i < 3; i++) {
	        s = strtok (NULL, WS);
		res = sscanf (s, "%lf", &bf);
		if (res) fvec[nd*3+i] = bf;
	    }
	}
	AddToRHS_elasticity (mesh, rhs, &fvec, RHS_BNDPF);
    }

#ifdef SOLVE_CHOLESKY
    cout << "Starting Cholesky solver" << endl;
    // calculate factorisation matrix
    RVector Fd(dn);
    F.SymbolicCholeskyFactorize (drowptr, dcolidx);
    cout << "Completed symbolic factorisation" << endl;
    RCompRowMatrix F_factor (dn, dn, drowptr, dcolidx);
    delete []drowptr;
    delete []dcolidx;

    // Factorize system
    CholeskyFactorize (F, F_factor, Fd, false);
    cout << "Completed factorisation" << endl;

    // Solve system
    CholeskySolve (F_factor, Fd, rhs, x);
    cout << "Completed solve" << endl;
#else
    cout << "Starting linear solver" << endl;
    RPreconditioner *precon = new RPrecon_Diag;
    precon->Reset (&F);
    int niter;
#ifdef OUTPUT_SYSMAT
    //ofstream sysm("sysmat.dat");
    //F.Print (sysm);
#endif
#ifdef OUTPUT_RHS
    //ofstream frhs("rhs.dat");
    //frhs << rhs << endl;
#endif
    niter = BiCGSTAB (F, rhs, x, tol, precon);
    //niter = GMRES (F, rhs, x, tol, precon);
    delete precon;
    cout << "converged after " << niter << " iterations" << endl;
#endif

#ifdef OUTPUT_SOLUTION
    //ofstream fsol("sol.dat");
    //fsol << x << endl;
#endif

    //cout << "after writing .dat files" << endl;
    RVector SED(mesh.elen()); // element-based strain energy density
    
    // strain energy of an element is defined as
    // 1/2 u^T K u
    // with u: nodal displacements and K stiffness matrix

    // calculate strain tensor components from displacement field
    for (el = 0; el < mesh.elen(); el++) {
        Element *pel = mesh.elist[el];
	int nd, nn = pel->nNode();

	// extract element displacements
	RVector u(nn*dim);
	RVector Ku(nn*dim);
	for (i = 0; i < nn; i++) {
	    nd = pel->Node[i];
	    for (j = 0; j < dim; j++)
	        u[i*dim+j] = x[nd*dim+j];
	}
	pel->ElasticityStiffnessMatrix (modulus[el], pratio[el]).Ax (u, Ku);
	SED[el] = 0.5 * (u & Ku); // element strain energy
    }


    // Add displacements to mesh
    if (add_displacements) {
        for (i = 0; i < n; i++) {
	    for (j = 0; j < 3; j++)
	    {
	        mesh.nlist[i][j] += x[i*3+j];
	    }
	}
	mesh.Setup();
    }

    for (i = 0; i < n; i++) {
        dsp = fabs(x[i*3+2]);
	if (dsp < mind) mind = dsp;
	if (dsp > maxd) maxd = dsp;
	dsp = fabs(x[i*3+0]);
	if (dsp < mind) mind = dsp;
	if (dsp > maxd) maxd = dsp;
	dsp = fabs(x[i*3+1]);
	if (dsp < mind) mind = dsp;
	if (dsp > maxd) maxd = dsp;
    }
    cout << "Displacement range: " << mind << " to " << maxd << endl;
    
    if (disp_range)
      {
	mind = mind2;
	maxd = maxd2;
      }
    	
    // Mesh parameter encoding 
    if (enc_disp)
      {
	strcpy(mesh_out_name_disp,mesh_out_name);
	strcat(mesh_out_name_disp,"_disp.msh");
	ofstream ofs (mesh_out_name_disp);
	ofs << mesh.nlen() << endl;
	for (i = 0; i < mesh.nlen(); i++) {
	  for (j = 0; j < 3; j++)
	    ofs << x[i*3+j] << ' ';
	  ofs << endl;
	}
	ofs.close();
      }

    int tmp = 0;    
    if (enc_strain) {
      
      strcpy(mesh_out_name_strain,mesh_out_name);
      strcat(mesh_out_name_strain,"_strain.msh");
      ofstream ofs (mesh_out_name_strain);
      ofs << mesh.elen() << endl;
      for (el = 0; el < mesh.elen(); el++) {
	Element *pel = mesh.elist[el];
	double sed = SED[el];
	ofs << sed << endl;
	//ofs << sed << " ";
	//tmp++;
	//if (tmp == 3)
	//  {
	//    ofs << endl;
	//    tmp = 0;
	//  }
      }
    } 
    tmp = 0;
    if (enc_elvolchange) {
      
      strcpy(mesh_out_name_elvolchange,mesh_out_name);
      strcat(mesh_out_name_elvolchange,"_elchange.msh");
      strcpy(mesh_out_name_realvolchange,mesh_out_name);
      strcat(mesh_out_name_realvolchange,"_realchange.msh");
      ofstream ofs (mesh_out_name_elvolchange);
      ofstream ofs_real (mesh_out_name_realvolchange);
      ofs << mesh.elen() << endl;
      ofs_real << mesh.elen() << endl;
      for (el = 0; el < mesh.elen(); el++) {
	Element *pel = mesh.elist[el];
	//int dsize = (int)( ( 1000.0 * (pel->Size()/size0[el]) )- 1000.0 );
	int dsize = (int)( 1000.0 * (pel->Size()/size0[el]) * 255 / 3000 );
	float dsize_real = (float)(( pel->Size()/size0[el] )-1);
	ofs << dsize << endl;
	ofs_real << dsize_real << endl;
	//ofs << dsize << " ";
	//tmp++;
	//if (tmp == 3)
	//  {
	//    ofs << endl;
	//    tmp = 0;
	// }
      }
    }
    if (enc_vol) {
      
      strcpy(mesh_out_name_vol,mesh_out_name);
      strcat(mesh_out_name_vol,"_vol.msh");
      ofstream ofs_vol (mesh_out_name_vol);
      ofs_vol << mesh.elen() << endl;
      for (el = 0; el < mesh.elen(); el++) {
	//Element *pel = mesh.elist[el];
	//ofs_vol << size0[el] << " " << pel->Size() << endl;
	ofs_vol << size0[el] << " " << mesh.ElSize(el) << endl;
      }
    }
    
    cout << "Total mesh size after displacement: " << mesh.FullSize() << endl;

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
    
