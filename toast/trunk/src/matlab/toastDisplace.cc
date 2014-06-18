//#include <fstream.h>
#include <math.h>
#include <stdlib.h>
#include "mathlib.h"
#include "felib.h"
#include "toastmex.h"
#include "mex.h"
#include "util.h"
using namespace std;

const double M_PI=3.141592653589793;
//#define OUTPUT_SYSMAT 1
//#define OUTPUT_SOLUTION 1
//#define OUTPUT_RHS 1
//#define ENCODE_DISPLACEMENT
//#define READ_VOXEL_MESH
//#define EXPORT_REGION_IMAGE
#define ENCODE_STRAINENERGY
//#define OUTPUT 1
//#define SOLVEIT 1
#define SANJAY 1 // debugging. 1 for debug. 0 for no debug.

#define BIGSPRING 1e7
//#define SOLVE_CHOLESKY

const char *WS = " \t";
const bool add_displacements = true;


void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	double tol = 1e-12;
	char bname[256], dname[256], cbuf[256], *s;
	double pratio0, modulus0, te0, disp, sz, tot_sz;
	//Mesh mesh;
	int nmat;
	//struct MatList { double E, nu, dns, te; } *matlist;
	int ex, ey, ez;
	double dx, dy, dz;
	bool *egrid;
	int i, j, k, nd, res, cmd, el;
	int *rowptr, *colidx, nzero;
	int *drowptr, *dcolidx, dn;
	int encoding;
	bool   bAuto;
	srand(1234567);


	//mesh handle
	Mesh *mesh = (Mesh*)Handle2Ptr (mxGetScalar (prhs[0]));

	//mesh->Setup ();

	int n = mesh->nlen();
	int elen = mesh->elen();
	int dim = mesh->Dimension();
	//cerr<<n<<"  "<<elen<<"  "<<dim<<endl;

	// record element volumes
	RVector size0(elen);
	double elsizemin= 0.0, elsizemax= 0.0, elsizeavg = 0.0;
	for (i = 0; i < elen; i++) {
		size0[i] = mesh->ElSize (i);
		if (size0[i] < elsizemin || !i) elsizemin = size0[i];
		if (size0[i] > elsizemax || !i) elsizemax = size0[i];
		elsizeavg += size0[i]/(double)elen;
	}

	// Generate system stiffness matrix
	mesh->SparseRowStructure (rowptr, colidx, nzero);
	BlockExpand (rowptr, colidx, n, drowptr, dcolidx, dn, dim, dim);
	RCompRowMatrix F (dn, dn, drowptr, dcolidx);
	delete []rowptr;
	delete []colidx;
	delete []drowptr;
	delete []dcolidx;

	RVector rhs(dn);
	RVector x(dn);

	RVector pratio(mesh->elen());
	RVector modulus(mesh->elen());
	RVector te(mesh->elen());
	IVector matidx(mesh->elen());


	//Element parameters 

	CopyVector (modulus, prhs[1]);
	CopyVector (pratio, prhs[2]);

	if (pratio.Dim() !=mesh->elen())
		mexErrMsgTxt("Poisson Ratio vector should have dimension of number of elements\n");
	if (modulus.Dim() !=mesh->elen())
		mexErrMsgTxt( "Youngs modulus vector should have dimension of number of elements\n");
	// assign material properties



	cout << "Total mesh size before displacement: " << mesh->FullSize() << endl;
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
				for (j = 0; j < mesh->elist[i]->nNode(); j++)
					cout << " " << mesh->elist[i]->Node[j];
				cout << endl;
			}

			// try to relax the mesh around bad elements
			Mesh mesh2(*mesh);
			int nbadnd = 0, nrelaxed = 0;
			for (i = 0; i < elen; i++)
				if (size0[i] <= 0.0) nbadnd += mesh->elist[i]->nNode();
			int *badnd = new int[nbadnd];
			for (i = nbadnd = 0; i < elen; i++)
				if (size0[i] <= 0.0)
					for (j = 0; j < mesh->elist[i]->nNode(); j++)
						badnd[nbadnd++] = mesh->elist[i]->Node[j];
#ifdef UNDEF
			// check for nodes appearing more than once
			for (i = 0; i < nbadnd-1; i++)
				for (j = i+1; j < nbadnd; j++)
					if (badnd[i] == badnd[j] && badnd[j] >= 0) {
						cout << "Relaxing node " << badnd[i] << endl;
						Point bc = mesh->NeighbourBarycentre(badnd[i]);
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
							Point bc = mesh->NeighbourBarycentre(badnd[i]);
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
//cout<<F(0,0)<<" M "<<modulus[0]<<" P "<<pratio[0]<<endl;
	AddToSysMatrix_elasticity (*mesh, F, modulus, pratio);
	cout << "Finished!" << endl;
	CopyMatrix (&plhs[2],F);
//cout<<F(0,0)<<" M "<<modulus[0]<<" P "<<pratio[0]<<endl;
	// add explicit boundary displacements using Payne-Irons
	// Dirichlet condition

	if (!mxIsChar(prhs[3]))
	{
		cout << "adding explicit boundary displacements" << endl;
		RDenseMatrix dbmat;
		CopyMatrix (dbmat, prhs[3]);
		for (int j = 0; j<dbmat.nRows(); j++){
			int nd = (int) (dbmat(j,0)-.5) ;
			for (i = 0; i < 3; i++) {
				double  disp = dbmat(j,i+1);
				//cout << "F(nd*3+i,nd*3+i) before = " << F(nd*3+i,nd*3+i) << endl;
				F(nd*3+i,nd*3+i) *= BIGSPRING;
				//cout << "F(nd*3+i,nd*3+i) after product with BIGSPRING = " << F(nd*3+i,nd*3+i) << endl;
				rhs[nd*3+i] = F(nd*3+i,nd*3+i)*disp;
				//cout << "rhs(nd*3+i) = " << rhs[nd*3+i] << endl;
			}
		}
	}



	// add volume forces

	if (!mxIsChar(prhs[5]))
	{
		
		RVector vfam;
		CopyVector (vfam, prhs[5]);
		RVector fvec(dn);
		RVector dir(dim);
		double vf, angle;
		angle = vfam[0];
		angle *= M_PI/180.0;
		dir[0] = cos (angle);
		dir[2] = sin (angle);
		vf=vfam[1];
		cout << "Add homog. volume forces " <<endl;
		cout << "Force azimuth angle [deg]: " << angle << endl;
		cout << "Force magnitude: " << vf << endl;
		for (i = 0; i < n; i++)
			for (j = 0; j < dim; j++)
				fvec[i*3+j] = vf*dir[j];
		AddToRHS_elasticity (*mesh, rhs, &fvec, RHS_P);
	}

	// add boundary forces



	if (!mxIsChar(prhs[4])){
		cout << "adding boundary forces" << endl;
		RDenseMatrix bfmat;
		RVector fvec(dn);
		CopyMatrix (bfmat, prhs[4]);
		for (int j = 0; j<bfmat.nRows(); j++){
			int nd = (int) (bfmat(j,0)-.5) ;
			if (!mesh->nlist[nd].isBnd()) {
				mexErrMsgTxt("Can't apply surface force to internal node\n");
			}

			for (i = 0; i < 3; i++) {
				double  bf = bfmat(j,i+1);
				// cout << "i = " << i << ", res = " << res << endl;
				fvec[nd*3+i] = bf;
				//cout << "nd*3+i = " << nd*3+i << ", bf=fvec[nd*3+i]= " << bf << endl;

			}
		}
		AddToRHS_elasticity (*mesh, rhs, &fvec, RHS_BNDPF);
	}

#ifdef SOLVEIT

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
	niter = BiCGSTAB (F, rhs, x, tol, precon);
	//niter = GMRES (F, rhs, x, tol, precon);
	delete precon;
	cout << "converged after " << niter << " iterations" << endl;
#endif

	CopyVector (&plhs[0],x);
	CopyMatrix (&plhs[1],F);
	CopyVector (&plhs[2],rhs);
#else 
	CopyMatrix (&plhs[0],F);
	CopyVector (&plhs[1],rhs);
#endif



#ifdef OUTPUT








	RVector SED(mesh->elen()); // element-based strain energy density

	// strain energy of an element is defined as
	// 1/2 u^T K u
	// with u: nodal displacements and K stiffness matrix

	// calculate strain tensor components from displacement field
	for (el = 0; el < mesh->elen(); el++) {
		Element *pel = mesh->elist[el];
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
			cout << "node " << i << ", disp = ";
			for (j = 0; j < 3; j++)
				mesh->nlist[i][j] += x[i*3+j];
			cout << x[i*3+j] << ", "; 
		}
		mesh->Setup();
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

	cout << "Mesh parameter encoding:\n";
	cout << "(0) none\n";
	cout << "(1) displacements\n";
	cout << "(2) strain energy density\n";
	cout << "(3) element volume change\n";
	cin  >> encoding;

	if (encoding == 1) {

		// Encode z-displacement magnitude in mua parameter
		cout << "z-displacement" << endl;
		for (i = 0; i < n; i++) {
			dsp = fabs(x[i*3+2]);
			cout << "i = " << i << ", i*3+2 = " << i*3+2 << ", x[i*3+2] = " << x[i*3+2] << ", dsp = " << dsp << endl;
			mesh->plist[i].SetMua (dsp);
			//mesh->plist[i].SetMua ((dsp-mind)/(maxd-mind));
		}

		// Encode x-displacement component in mus parameter
		cout << "x-displacement" << endl;
		for (i = 0; i < n; i++) {
			dsp = fabs(x[i*3+0]);
			cout << "i = " << i << ", i*3+2 = " << i*3+2 << ", x[i*3+2] = " << x[i*3+2] << ", dsp = " << dsp << endl;
			mesh->plist[i].SetMus(dsp);
			//mesh->plist[i].SetMus ((dsp-mind)/(maxd-mind));
		}

		// Encode y-displacement component in n parameter
		cout << "y-displacement" << endl;
		for (i = 0; i < n; i++) {
			dsp = fabs(x[i*3+1]);
			cout << "i = " << i << ", i*3+2 = " << i*3+2 << ", x[i*3+2] = " << x[i*3+2] << ", dsp = " << dsp << endl;
			mesh->plist[i].SetN (dsp);
			//mesh->plist[i].SetN ((dsp-mind)/(maxd-mind));
		}

	} else if (encoding == 2) {

		RVector nSED(mesh->nlen());
		IVector count(mesh->nlen());
		for (el = 0; el < mesh->elen(); el++) {
			Element *pel = mesh->elist[el];
			double sed = SED[el];
			for (i = 0; i < pel->nNode(); i++) {
				nd = pel->Node[i];
				nSED[nd] += sed;
				count[nd]++;
			}
		}
		for (nd = 0; nd < mesh->nlen(); nd++)
			mesh->plist[nd].SetMua (nSED[nd]/count[nd]);

	} else if (encoding == 3) {

		RVector ds(mesh->nlen());
		RVector count(mesh->nlen());
		for (el = 0; el < mesh->elen(); el++) {
			Element *pel = mesh->elist[el];
			double dsize = pel->Size()/size0[el];
			for (i = 0; i < pel->nNode(); i++) {
				nd = pel->Node[i];
				ds[nd] += dsize;
				count[nd]++;
			}
		}
		for (nd = 0; nd < mesh->nlen(); nd++)
			mesh->plist[nd].SetMua (ds[nd]/count[nd]);

	}

	cout << "Total mesh size after displacement: " << mesh->FullSize() << endl;

	// Write out mesh
	ofstream ofs ("displace.msh");
	ofs << mesh << endl;
#endif
	//   return 0; // end of the mex function
}                                                                              
