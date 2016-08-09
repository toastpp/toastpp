#include <fstream>
#include <stdlib.h>
#include "mathlib.h"
#include "felib.h"

using namespace std;

//#define OUTPUT_SYSMAT 1
//#define OUTPUT_SOLUTION 1
//#define OUTPUT_RHS 1
//#define ENCODE_DISPLACEMENT
//#define READ_VOXEL_MESH
//#define EXPORT_REGION_IMAGE
#define ENCODE_STRAINENERGY

#define SANJAY 1 // debugging. 1 for debug. 0 for no debug.

#define BIGSPRING 1e7
//#define SOLVE_CHOLESKY

const char *WS = " \t";
const bool add_displacements = true;

void AddPerturbations (const Mesh &mesh, const IVector &matidx, RVector &te);
void WriteNimHeader (char *meshname, int imgsize, char *nimname, char *type);
void WriteNim (const RVector &nim, int imgno, const char *nimname);
void WriteEimHeader (const char *meshname, int imgsize, const char *eimname,
    const char *type);
void WriteEim (const RVector &eim, int imgno, const char *eimname);
RVector ReadEim (const char *eimname);

int main (void) {
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

    srand(1234567);

#ifdef READ_VOXEL_MESH
    cout << "Binary file name: ";
    cin >> bname;

    cout << "Voxel dimensions (x y z): ";
    cin >> ex >> ey >> ez;

    cout << "Physical dimensions of a voxel (dx dy dz): ";
    cin >> dx >> dy >> dz;

    egrid = new bool[ex*ey*ez];
    FILE *f = fopen (bname, "rb");
    fread (egrid, sizeof (bool), ex*ey*ez, f);
    fclose (f);

    CreateVoxelMesh (ex, ey, ez, egrid, dx, dy, dz, mesh);
#else
    cout << "Mesh file name: ";
    cin >> bname;
    ifstream mesh_f(bname);
    mesh_f >> mesh;
#endif

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

    cout << "Displacement file ('-' for none): ";
    cin >> dname;
    if (bDisp = (strcmp (dname, "-") != 0))
        disp_f.open (dname);

    cout << "Boundary force file ('-' for none): ";
    cin >> dname;
    if (bBndf = (strcmp (dname, "-") != 0))
        bndf_f.open (dname);

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

    cout << "Element parameters  (1) homogeneous  (2) from file: ";
    cin >> cmd;
    if (cmd == 1) {
	nmat = 1;
	matlist = new MatList[nmat];
        cout << "Poisson's ratio: ";
	cin  >> matlist[0].nu;
	cout << "Young's modulus: ";
	cin  >> matlist[0].E;
	cout << "Thermal expansion: ";
	cin  >> matlist[0].te;
	matidx = 0;
    } else {
        int nidx, idx;
        cout << "Material file: ";
	cin >> cbuf;
	ifstream ifs (cbuf);
	ifs >> nmat;
	matlist = new MatList[nmat];
	for (i = 0; i < nmat; i++) {
	    ifs >> matlist[i].E >> matlist[i].nu >> matlist[i].dns
		>> matlist[i].te;
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

    cout << "Apply perturbations to thermal coefficients" << endl;
    cout << "0: No" << endl;
    cout << "1: Gaussian element noise" << endl;
    cout << "2: Random Gaussian blobs" << endl;
    cout << "3: Localised perturbations" << endl;
    cout << "4: Read coefficient file" << endl;
    cin >> cmd;
    switch (cmd) {
    case 1: {
	double noise;
	cout << "Noise level: ";
	cin >> noise;
	for (i = 0; i < elen; i++)
	    te[i] += gasdev (noise);
        } break;
    case 2: {
	int nblob;
	double noise, sigma, fac1, fac2, h, w;
	Point bbmin, bbmax, cnt(mesh.Dimension());
	mesh.BoundingBox (bbmin, bbmax);
	cout << "Number of blobs: ";
	cin >> nblob;
	cout << "Blob level: ";
	cin >> noise;
	cout << "Blob width: ";
	cin >> w;
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
        } break;
    case 3:
	AddPerturbations (mesh, matidx, te);
        break;
    case 4: {
	char eimname[256];
	cout << "thermal coefficient file (EIM): ";
	cin >> eimname;
	te = ReadEim (eimname);
	} break;
    }
    WriteEimHeader (bname, elen, "thcoeff.eim", "N/A");
    WriteEim (te, 0, "thcoeff.eim");
    ofstream ofs7 ("tmp2.txt");   
    
    ofs7 << "Thermal coefficient distribution written to thcoeff.eim"
	 << endl;
    
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
    WriteNimHeader (bname, n, "material.nim", "N/A");
    WriteNim (ndmat, 1, "material.nim");
#endif 

    ofs7 << "Total mesh size before displacement: " << mesh.FullSize() << endl;
    ofs7 << "Element sizes: ";
    ofs7 << "min = " << elsizemin
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
    ofs7 << "Finished!" << endl;

    // add explicit boundary displacements using Payne-Irons
    // Dirichlet condition
    if (bDisp) {
        ofs7 << "adding explicit boundary displacements" << endl;
        while (disp_f.getline (cbuf, 256)) {
	    s = strtok (cbuf, WS);
	    ofs7 << "s= " << s << endl;
	    res = sscanf (s, "%d", &nd);
	    ofs7 << "res = " << res << endl;
	    if (!res) continue;
	    for (i = 0; i < 3; i++) {
	        s = strtok (NULL, WS);
		res = sscanf (s, "%lf", &disp);
		ofs7 << " i = " << i << ", disp = " << disp << endl;
		if (res) {
		    ofs7 << "F(nd*3+i,nd*3+i) before = " << F(nd*3+i,nd*3+i) << endl;
		    F(nd*3+i,nd*3+i) *= BIGSPRING;
		    ofs7 << "F(nd*3+i,nd*3+i) after product with BIGSPRING = " << F(nd*3+i,nd*3+i) << endl;
		    rhs[nd*3+i] = F(nd*3+i,nd*3+i)*disp;
		    ofs7 << "rhs(nd*3+i) = " << rhs[nd*3+i] << endl;
		}
	    }
	}
    }

    // add thermal expansion
    double dT;
    cout << "Applied temperature change [K]: ";
    cin >> dT;
    if (dT) 
	AddToRHS_thermal_expansion (mesh, rhs, modulus, pratio, te, dT);

    // add volume forces
    cout << "Add homog. volume forces (1/0)? ";
    cin >> cmd;
    if (cmd) {
	ofs7 << "adding volume forces" << endl;
        double vf, angle;
	RVector fvec(dn);
	RVector dir(dim);
	cout << "Force azimuth angle [deg]:";
	cin >> angle;
	angle *= M_PI/180.0;
	dir[0] = cos (angle);
	dir[2] = sin (angle);
        cout << "Force magnitude: ";
	cin >> vf;
	for (i = 0; i < n; i++)
	    for (j = 0; j < dim; j++)
	        fvec[i*3+j] = vf*dir[j];
	AddToRHS_elasticity (mesh, rhs, &fvec, RHS_P);
    }

    // add boundary forces
    if (bBndf) {
        ofs7 << "adding boundary forces" << endl;
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
		ofs7 << "i = " << i << ", res = " << res << endl;
		if (res) fvec[nd*3+i] = bf;
		ofs7 << "nd*3+i = " << nd*3+i << ", bf=fvec[nd*3+i]= " << bf << endl;
	    }
	}
	AddToRHS_elasticity (mesh, rhs, &fvec, RHS_BNDPF);
    }

#ifdef SOLVE_CHOLESKY
    ofs7 << "Starting Cholesky solver" << endl;
    // calculate factorisation matrix
    RVector Fd(dn);
    F.SymbolicCholeskyFactorize (drowptr, dcolidx);
    ofs7 << "Completed symbolic factorisation" << endl;
    RCompRowMatrix F_factor (dn, dn, drowptr, dcolidx);
    delete []drowptr;
    delete []dcolidx;

    // Factorize system
    CholeskyFactorize (F, F_factor, Fd, false);
    ofs7 << "Completed factorisation" << endl;

    // Solve system
    CholeskySolve (F_factor, Fd, rhs, x);
    ofs7 << "Completed solve" << endl;
#else
    ofs7 << "Starting linear solver" << endl;
    RPreconditioner *precon = new RPrecon_Diag;
    precon->Reset (&F);
    int niter;
#ifdef OUTPUT_SYSMAT
    ofstream sysm("sysmat.dat");
    F.Print (sysm);
#endif
#ifdef OUTPUT_RHS
    ofstream frhs("rhs.dat");
    frhs << rhs << endl;
#endif
    niter = BiCGSTAB (F, rhs, x, tol, precon);
    //niter = GMRES (F, rhs, x, tol, precon);
    delete precon;
    ofs7 << "converged after " << niter << " iterations" << endl;
#endif

#ifdef OUTPUT_SOLUTION
    ofstream fsol("sol.dat");
    fsol << x << endl;
#endif

    ofs7 << "after writing .dat files" << endl;
    
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
	    ofs7 << "node " << i << ", disp = ";
	    for (j = 0; j < 3; j++)
	        mesh.nlist[i][j] += x[i*3+j];
		cout << x[i*3+j] << ", "; 
	}
	mesh.Setup();
    }

    double dsp, mind = 1e100, maxd = -1e100;
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
    ofs7 << "Displacement range: " << mind << " to " << maxd << endl;

    cout << "Displacement encoding range [min max]\n";
    cout << "(- for auto)" << endl;
    cin.getline (cbuf, 256);
    cin.getline (cbuf, 256);
    bAuto = (strcmp (cbuf, "-") == 0);
    if (!bAuto) {
        sscanf (cbuf, "%lf%lf", &mind, &maxd);
    }

    cout << "Write out displacements?\n";
    cout << "(0) no\n";
    cout << "(1) yes\n";
    cin  >> encoding;

    if (encoding) {
	ofstream ofs("displacement.dat");
	for (i = 0; i < n; i++) {
	    ofs << x[i*3+0] << '\t' << x[i*3+1] << '\t' << x[i*3+2] << endl;
	}
    }

    cout << "Write out strain enery densities?\n";
    cout << "(0) no\n";
    cout << "(1) yes\n";
    cin >> encoding;

    if (encoding) {
        RVector nSED(mesh.nlen());
	IVector count(mesh.nlen());
	for (el = 0; el < mesh.elen(); el++) {
	    Element *pel = mesh.elist[el];
	    double sed = SED[el];
	    for (i = 0; i < pel->nNode(); i++) {
	        nd = pel->Node[i];
		nSED[nd] += sed;
		count[nd]++;
	    }
	}
	ofstream ofs("SED.dat");
	for (nd = 0; nd < mesh.nlen(); nd++)
	    ofs << nSED[nd]/count[nd] << endl;
    }

    cout << "Write out element volume change?\n";
    cout << "(0) no\n";
    cout << "(1) yes\n";
    cin >> encoding;

    if (encoding) {
	RVector ds(mesh.nlen());
	RVector count(mesh.nlen());
	for (el = 0; el < mesh.elen(); el++) {
	    Element *pel = mesh.elist[el];
	    double dsize = pel->Size()/size0[el];
	    for (i = 0; i < pel->nNode(); i++) {
	        nd = pel->Node[i];
		ds[nd] += dsize;
		count[nd]++;
	    }
	}
	ofstream ofs("volchange.dat");
	for (nd = 0; nd < mesh.nlen(); nd++)
	    ofs << ds[nd] / count[nd] << endl;
    }

    ofs7 << "Total mesh size after displacement: " << mesh.FullSize() << endl;

    // Write out mesh
    ofstream ofs ("displace.msh");
    ofs << mesh << endl;

    return 0;
}                                                                              

void AddPerturbations (const Mesh &mesh, const IVector &matidx, RVector &te)
{
    int i, cmd, reg;
    int dim = mesh.Dimension();
    Point cnt(dim);
    double val, rad;
    do {
	cout << "Thermal coefficient perturbation value: ";
	cin >> val;
	cout << "Perturbation centre ";
	cout << (dim == 3 ? "[x y z]: " : "[x y]: ");
	cin >> cnt[0] >> cnt[1];
	if (dim == 3) cin >> cnt[2];
	cout << "Perturbation radius: ";
	cin >> rad;
	cout << "Region affected (-1 for all): ";
	cin >> reg;
	
	for (i = 0; i < mesh.elen(); i++) {
	    if (reg < 0 || matidx[i] == reg) {
		if (cnt.Dist (mesh.ElCentre(i)) < rad) {
		    te[i] += val;
		}
	    }
	}

	cout << "(1) continue" << endl;
	cout << "(0) finished" << endl;
	cout << "[1|0] >> ";
	cin >> cmd;
    } while (cmd == 1);
}

void WriteNimHeader (const char *meshname, int imgsize, const char *nimname,
    const char *type)
{
    ofstream ofs (nimname);
    ofs << "NIM" << endl;
    ofs << "Mesh = " << meshname << endl;
    ofs << "SolutionType = " << type << endl;
    ofs << "ImageSize = " << imgsize << endl;
    ofs << "EndHeader" << endl;
}

void WriteNim (const RVector &nim, int imgno, const char *nimname)
{
    ofstream ofs (nimname, ios::app);
    ofs << "Image " << imgno << endl;
    for (int i = 0; i < nim.Dim(); i++)
        ofs << nim[i] << ' ';
    ofs << endl;
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

RVector ReadEim (const char *eimname)
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
    
