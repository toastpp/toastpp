#include <fstream.h>
#include "mathlib.h"
#include "felib.h"

//#define OUTPUT_SYSMAT
//#define OUTPUT_SOLUTION
//#define OUTPUT_RHS
//#define ENCODE_DISPLACEMENT
#define ENCODE_STRAINENERGY

#define SANJAY 1 // debugging. 1 for debug. 0 for no debug.

#define BIGSPRING 1e7
//#define SOLVE_CHOLESKY

const char *WS = " \t";
const bool add_displacements = false;

int main (void) {
    double tol = 1e-12;
    char bname[256], dname[256], cbuf[256], *s;
    double pratio0, modulus0, disp, sz, tot_sz;
    Mesh mesh;
    int ex, ey, ez;
    double dx, dy, dz;
    bool *egrid;
    int i, j, k, nd, res, cmd, el;
    int *rowptr, *colidx, nzero;
    int *drowptr, *dcolidx, dn;
    bool bDisp, bBndf, bAuto;
    ifstream disp_f, bndf_f;

    cout << "Binary file name: ";
    cin >> bname;
 
   // ifstream mesh_f (bname); // modified.sk. as next line only.
    ifstream mesh_f(bname, ios::in|ios::binary);
    
    cout << "Voxel dimensions (x y z): ";
    // cin >> ex >> ey >> ez; // temporarily using next line:
    //ex = 154; ey = 155; ez = 142;
    ex = 154; ey = 155; ez = 100;
    cout << "Physical dimensions of a voxel (dx dy dz): ";
    // cin >> dx >> dy >> dz; // temporarily using next lines:
	dx = 0.056; dy = 0.056; dz = 0.06;

    egrid = new bool[ex*ey*ez];
    mesh_f.read ((bool*) egrid, ex*ey*ez*sizeof(bool));

    CreateVoxelMesh (ex, ey, ez, egrid, dx, dy, dz, mesh);

    ofstream fofs ("mesh_test.msh");
	if(SANJAY==2){	cout << mesh;
	cout << endl;}

    fofs << mesh;

    // mesh_f >> mesh; // this was when mesh were read in from a data file.
    mesh.Setup ();
	if(SANJAY==1) cout << "check " << endl;

    int n = mesh.nlen();
    int dim = mesh.Dimension();
    cout << "Total mesh size before displacement: " << mesh.FullSize() << endl;

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
    cout << "Element parameters  (1) homogeneous  (2) from file: ";
    cin >> cmd;
    if (cmd == 1) {
        cout << "Poisson's ratio: ";
	cin  >> pratio0;
	cout << "Young's modulus: ";
	cin  >> modulus0;
	pratio = pratio0;
	modulus = modulus0;
    } else {
        int nmat, nidx, idx;
        cout << "Material file: ";
	cin >> cbuf;
	ifstream ifs (cbuf);
	ifs >> nmat;
	struct MatList { double E, nu, dns; } *matlist = new MatList[nmat];
	for (i = 0; i < nmat; i++) {
	    ifs >> matlist[i].E >> matlist[i].nu >> matlist[i].dns;
	}
	ifs >> nidx;
	if (nidx != mesh.elen()) {
	    cerr << "Invalid length of material index list. Aborting." << endl;
	    exit(1);
	}
	for (i = 0; i < nidx; i++) {
	    ifs >> idx;
	    if (idx < 0 || idx >=nmat) {
	        cerr << "Invalid index in material index list. Aborting." << endl;
		exit (1);
	    }
	    pratio[i] = matlist[idx].nu;
	    modulus[i] = matlist[idx].E;
	}
    }
    cout << "Assembling system matrix" << endl;
    AddToSysMatrix_elasticity (mesh, F, modulus, pratio);
    cout << "Finished!" << endl;

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

#ifdef UNDEF
    // add explicit boundary displacements using Payne-Irons
    // Dirichlet condition
    if (bDisp) {
        while (disp_f.getline (cbuf, 256)) {
	    s = strtok (cbuf, WS);
	    res = sscanf (s, "%d", &nd);
	    if (!res) continue;
	    for (i = 0; i < 3; i++) {
	        s = strtok (NULL, WS);
		res = sscanf (s, "%lf", &disp);
		if (res) {
		    F(nd*3+i,nd*3+i) *= BIGSPRING;
		    rhs[nd*3+i] = F(nd*3+i,nd*3+i)*disp;
		}
	    }
	}
    }
#endif

    // add volume forces
    cout << "Add homog. volume forces (1/0)? ";
    cin >> cmd;
    if (cmd) {
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
    cout << "converged after " << niter << " iterations" << endl;
#endif

#ifdef OUTPUT_SOLUTION
    ofstream fsol("sol.dat");
    fsol << x << endl;
#endif

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
        for (i = 0; i < n; i++)
	    for (j = 0; j < 3; j++)
	        mesh.nlist[i][j] += x[i*3+j];
	//mesh.Setup();
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
    cout << "Displacement range: " << mind << " to " << maxd << endl;

    cout << "Displacement encoding range [min max]\n";
    cout << "(- for auto)" << endl;
    cin.getline (cbuf, 256);
    cin.getline (cbuf, 256);
    bAuto = (strcmp (cbuf, "-") == 0);
    if (!bAuto) {
        sscanf (cbuf, "%lf%lf", &mind, &maxd);
    }

#ifdef ENCODE_DISPLACEMENT
    // Encode z-displacement magnitude in mua parameter
    for (i = 0; i < n; i++) {
        dsp = fabs(x[i*3+2]);
	mesh.plist[i].SetMua (dsp);
        //mesh.plist[i].SetMua ((dsp-mind)/(maxd-mind));
    }

    // Encode x-displacement component in mus parameter
    for (i = 0; i < n; i++) {
        dsp = fabs(x[i*3+0]);
	mesh.plist[i].SetMus(dsp);
        //mesh.plist[i].SetMus ((dsp-mind)/(maxd-mind));
    }

    // Encode y-displacement component in n parameter
    for (i = 0; i < n; i++) {
        dsp = fabs(x[i*3+1]);
	mesh.plist[i].SetN (dsp);
        //mesh.plist[i].SetN ((dsp-mind)/(maxd-mind));
    }
#endif

#ifdef ENCODE_STRAINENERGY
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
    for (nd = 0; nd < mesh.nlen(); nd++)
        mesh.plist[nd].SetMua (nSED[nd]/count[nd]);
#endif
    cout << "Total mesh size after displacement: " << mesh.FullSize() << endl;

    // Write out mesh
    ofstream ofs ("displace.msh");
    ofs << mesh << endl;

    return 0;
}                                                                              
