// reconstruct thermal expansion parameters on a region basis
// given an original and a distorted mesh
// We use element size difference as objective function

#include <fstream>
#include <time.h>
#include "mathlib.h"
#include "felib.h"
#include <iostream>

#define ENCODE_STRAINENERGY
#define MAXREGION 100
#define BIGSPRING 1e7

#define METHOD 1
// method 1 uses region sizes as data
// method 2 uses nodal displacements
// method 3 uses similarity measure on sampled grids

using namespace std;
struct MatList {      // material properties
    double E;         // Young's modulus
    double nu;        // Poisson's ratio
    double dns;       // density
    double te;        // thermal expansion coefficient
};

struct Displacement { // nodal displacement record
    int nd;           // node
    double d[3];      // displacement vector
};

const char *WS = " \t";
const bool add_displacements = true;
const double DISP_NONE = 1e10; // flag for 'no displacement'

void WriteNimHeader (char *meshname, int imgsize, char *nimname, char *type);
void WriteNim (const RVector &nim, int imgno, char *nimname);
void WriteEimHeader (char *meshname, int imgsize, char *eimname, char *type);
void WriteEim (const RVector &eim, int imgno, char *eimname);


double of1 (const Mesh &mesh1, const Mesh &mesh2)
{
    // objective function version 1:
    // root mean square difference of element sizes
    
    double els1, els2, ds, sum = 0.0;

    for (int i = 0; i < mesh1.elen(); i++) {
	els1 = mesh1.elist[i]->Size();
	els2 = mesh2.elist[i]->Size();
	ds   = els1-els2;
	sum += ds*ds;
    }
    return sqrt(sum);
}

double of2 (const RVector &rsize1, const RVector &rsize2)
{
    // objective function version 2
    // relative region size differences

    int i, nreg = rsize1.Dim();
    double dr, sum = 0.0;

    for (i = 0; i < nreg; i++) {
	dr = (rsize1[i] - rsize2[i]) / rsize1[i];
	sum += dr*dr;
    }
    return sum;
}

double of3 (const Mesh &mesh1, const Mesh &mesh2)
{
    // objective function version 3
    // root sum square of node displacements
    int n = mesh1.nlen();
    int dim = mesh1.Dimension();
    int i, j;

    double d, of = 0.0;
    for (i = 0; i < n; i++) {
	for (j = 0; j < dim; j++) {
	    d = mesh1.nlist[i][j] - mesh2.nlist[i][j];
	    of += d*d;
	}
    }
    return sqrt(of);
}

double of4 (const int *el1, const int *el2, int n)
{
    // objective function version 4
    // correspondence of pixel values for sampled meshes

    double sum = 0.0;
    for (int i = 0; i < n; i++)
	if (el1[i] != el2[i]) sum += 1.0;
    return sum/n;
}

double RegCorr (const int *el1, const int *el2, int n, int reg)
{
    // returns a measure for misalignment of region 'reg' between the
    // two sampling grids
    double sum = 0.0;
    for (int i = 0; i < n; i++) {
	if ((el1[i] == reg && el2[i] != reg) ||
	    (el2[i] == reg && el1[i] != reg))
	    sum += 1.0;
    }
    return sum/n;
}

void Regsize (const Mesh &mesh, int *matidx, RVector &grpsize)
{
    // returns the size of each region in the mesh, as defined by the
    // matidx index list

    grpsize = 0.0; // reset
    int i, reg, elen = mesh.elen();
    for (i = 0; i < elen; i++) {
	reg = matidx[i];
	if (reg >= grpsize.Dim()) {
	    cerr << "Index out of range" << endl;
	    exit(0);
	}
	grpsize[reg] += mesh.ElSize(i);
    }
}

void RegsizePerc (const Mesh &mesh, RDenseMatrix &percs, int *matidx, RVector &grpsize)
{
    // returns the size of each region in the mesh, as defined by the
    // matidx index list

    grpsize = 0.0; // reset
    int i, j, reg, elen = mesh.elen();
    for (i = 0; i < elen; i++) {
	//reg = matidx[i];
	//if (reg >= grpsize.Dim()) {
	//    cerr << "Index out of range" << endl;
	 //   exit(0);
	//}
	//cout << "element " << i << ", region " << reg << ", size tetra = " << mesh.ElSize(i) << endl;
	//cout << "volume total region = " << grpsize[reg] << endl;
	for (j = 0; j < grpsize.Dim(); j++)
	{
	  //grpsize[reg] += mesh.ElSize(i)*percs(i,reg);
	  grpsize[j] += mesh.ElSize(i)*percs(i,j);
	  //cout << "percs = " << percs(i,reg) << ", volume region after = " << grpsize[reg] << endl; 
	}
    }
}

void Displace (Mesh &mesh, const RVector &x)
{
    // apply displacements to mesh nodes
    int n = mesh.nlen(), dim = mesh.Dimension();
    for (int i = 0; i < n; i++)
	for (int j = 0; j < dim; j++)
	    mesh.nlist[i][j] += x[i*dim+j];
    mesh.Setup();
}

void SetupSysmat (const Mesh &mesh, RCompRowMatrix &F)
{
    int *rowptr, *colidx, *drowptr, *dcolidx, nzero, dn; 
    int n = mesh.nlen();
    int dim = mesh.Dimension();

    mesh.SparseRowStructure (rowptr, colidx, nzero);
    BlockExpand (rowptr, colidx, n, drowptr, dcolidx, dn, dim, dim);
    F.New (dn, dn);
    F.Initialise (drowptr, dcolidx);
}


void Solve (const Mesh &mesh, RCompRowMatrix &F, RDenseMatrix &disp,
    MatList *mat, int *matidx, int nreg, RVector &x)
{
    const double dT = 1.0; // we assume a temperature diff of 1K

    int elen = mesh.elen(), nlen = mesh.nlen(), dim = mesh.Dimension();
    int i, j, reg, niter, dn = nlen*dim;
    double fd;
    double tol = 1e-6;
    RVector E(elen), nu(elen), te(elen);
    RVector rhs(dn);

    for (i = 0; i < elen; i++) {
	reg = matidx[i];
	E[i] = mat[reg].E;
	nu[i] = mat[reg].nu;
	te[i] = mat[reg].te;
    }
    F.Zero();

    // add elasticity terms
    AddToSysMatrix_elasticity (mesh, F, E, nu);

    // add thermal expansion
    AddToRHS_thermal_expansion (mesh, rhs, E, nu, te, dT);

    // add explicit boundary displacements using Payne-Irons
    // Dirichlet condition
    for (i = 0; i < nlen; i++)
	for (j = 0; j < dim; j++)
	    if (disp(i,j) != DISP_NONE) {
		fd = F(i*dim+j, i*dim+j) *= BIGSPRING;
		rhs[i*3+j] = fd * disp(i,j);
	    }

    // solve linear system
    RPrecon_Diag precon;
    //RPrecon_IC precon;
    //RPrecon_Null precon;
    precon.Reset (&F);
    x.New(dn);
    //niter = GMRES (F, rhs, x, tol, &precon);
    niter = BiCGSTAB (F, rhs, x, tol, &precon, dn/10);
    //niter = PCG (F, rhs, x, tol, &precon);
}

int main (void) {

    int i, j, k, nlen, elen, el, dim, nreg, nreg1, nreg2, *matidx1, *matidx2, *nelreg1, *nelreg2, itmax = 4;
    bool bDisp, diffMat, perc1, perc2;
    double *regsize1[2], *regsize2[2];
    char fname[256], mname1[256], mname2[256], mname_perc1[256], mname_perc2[256],dname[256], cbuf[256];
    char meshname[256];
    double perc;
    ifstream ifs;
    Mesh mesh1, mesh2;
    MatList *mat1, *mat2;
    //double **percs1, **percs2;

    cout << "Original mesh name: ";
    cin >> meshname;
    ifs.open (meshname);
    ifs >> mesh1;
    ifs.close();
    mesh1.Setup();

    cout << "Distorted mesh name: ";
    cin >> fname;
    ifs.open (fname);
    ifs >> mesh2;
    ifs.close();
    mesh2.Setup();

    cout << "Material/region original mesh name: ";
    cin >> mname1;

    cout << "Percentage file Material/region original mesh name: ";
    cin >> mname_perc1;
    perc1 = (strcmp (mname_perc1, "-") != 0);

    cout << "Material/region distorted mesh name (- for none): ";
    cin >> mname2;
    diffMat = (strcmp (mname2, "-") != 0);

    cout << "Percentage file Material/region distorted mesh name: ";
    cin >> mname_perc2;
    perc2 = (strcmp (mname_perc2, "-") != 0);

    cout << "Displacement file name (- for none): ";
    cin >> dname;
    bDisp = (strcmp (dname, "-") != 0);
    cout << endl;

    // sanity checks
    nlen = mesh1.nlen();
    elen = mesh1.elen();
    dim  = mesh1.Dimension();
    if (nlen == mesh2.nlen() && elen == mesh2.elen()) {
	cout << "Nodes: " << nlen << endl;
	cout << "Elements: " << elen << endl;
    } else {
	cerr << "Meshes incompatible" << endl;
	exit (1);
    }
    cout << "Mesh 1 size: " << mesh1.FullSize() << endl;
    cout << "Mesh 2 size: " << mesh2.FullSize() << endl;
    cout << "Element size difference measure: " << of1 (mesh1, mesh2) << endl;
    
    // scan original mesh regions
    ifs.open (mname1);
    ifs >> nreg1;
    cout << "original mesh material file, nreg = " << nreg1 << endl;
    mat1 = new MatList[nreg1];
    nelreg1 = new int[nreg1];
    regsize1[0] = new double[nreg1];
    regsize1[1] = new double[nreg1];
    RVector rgsize1(nreg1);
    for (i = 0; i < nreg1; i++) {
	ifs >> mat1[i].E >> mat1[i].nu >> mat1[i].dns >> mat1[i].te;
	nelreg1[i] = 0;
	for (j = 0; j < 2; j++)
	    regsize1[j][i] = 0.0;
    }
    ifs >> i;
    if (i != elen) {
	cerr << "Material index list: invalid length" << endl;
	exit (1);
    }
    matidx1 = new int[elen];
    for (i = 0; i < elen; i++) {
	ifs >> matidx1[i];
	if (matidx1[i] < 0 || matidx1[i] >= nreg1) {
	    cerr << "Material index list: index out of range" << endl;
	    exit (1);
	}
	nelreg1[matidx1[i]]++;
	regsize1[0][matidx1[i]] += mesh1.ElSize (i);
	regsize1[1][matidx1[i]] += mesh2.ElSize (i);
    }
    cout << "Found " << nreg1 << " regions:\n";
    ifs.close();
    
    // scan distorted mesh regions
    if (diffMat)
    {
      ifs.open (mname2);
      ifs >> nreg2;
      cout << "distorted mesh material file, nreg = " << nreg2 << endl;
      if (nreg1 != nreg2){
	    cerr << "Material files with different number of regions" << endl;
	    exit (1);
	}
      nreg = nreg2;
      mat2 = new MatList[nreg];
      nelreg2 = new int[nreg];
      regsize2[0] = new double[nreg];
      regsize2[1] = new double[nreg];
      RVector rgsize2(nreg);
      for (i = 0; i < nreg; i++) {
	  ifs >> mat2[i].E >> mat2[i].nu >> mat2[i].dns >> mat2[i].te;
	 // cout << "region " << i << ", (E,nu,dns,te)=(" << 
	 // 	mat2[i].E << "," << mat2[i].nu << "," << mat2[i].dns << "," << mat2[i].te
	 // 	<< ")" << endl;
	  nelreg2[i] = 0;
	  for (j = 0; j < 2; j++)
	      regsize2[j][i] = 0.0;
      }
      ifs >> i;
     // cout << "ifs -> i = " << i << endl;
      if (i != elen) {
	  cerr << "Material index list: invalid length" << endl;
	  exit (1);
      }
      matidx2 = new int[elen];
      for (i = 0; i < elen; i++) {
	  ifs >> matidx2[i];
	 // cout << "elen = " << i << ", matidx2 = " << matidx2[i] << endl;
	  if (matidx2[i] < 0 || matidx2[i] >= nreg) {
	      cerr << "Material index list: index out of range" << endl;
	      exit (1);
	  }
	  nelreg2[matidx2[i]]++;
	  regsize2[0][matidx1[i]] += mesh1.ElSize (i);
	  regsize2[1][matidx2[i]] += mesh2.ElSize (i);
      }
      cout << "Found " << nreg << " regions:\n";
      ifs.close();
    }//end scan distorted mesh
    
    RDenseMatrix percs1(elen,nreg1);
    
    if (perc1)
    {
      ifs.open (mname_perc1);
      
      ifs >> i;
      if (i != elen) {
	  cerr << "Partial volume list: invalid length" << endl;
	  exit (1);
      }
      for (i = 0; i < elen; i++) {
	  ifs >> el;
	  for (j = 0; j < 7; j++) {
	      ifs >> perc;
	      if (j < nreg1)
	      {
		perc *= 0.01;
	        percs1(i,j) = perc;
	      }
	  }
      }
      ifs.close();
    }
    
    RDenseMatrix percs2(elen,nreg);

    if (perc2)
    {
      ifs.open (mname_perc2);
      
      ifs >> i;
      if (i != elen) {
	  cerr << "Partial volume list: invalid length" << endl;
	  exit (1);
      }
      for (i = 0; i < elen; i++) {
	  ifs >> el;
	  for (j = 0; j < 7; j++) {
	      ifs >> perc;
	      if (j < nreg)
	      {
	        perc *= 0.01;
	        percs2(i,j) = perc;
	      }
	  }
      }
      ifs.close();
    }

    // scan displacement file
    RDenseMatrix Disp(nlen, 3);
    Disp = DISP_NONE;
    if (bDisp) {
	char cbuf[256], *s;
	int res, nd;
	double disp;
	ifs.open (dname);
        while (ifs.getline (cbuf, 256)) {
	    s = strtok (cbuf, WS);
	    res = sscanf (s, "%d", &nd);
	    if (!res) continue;
	    for (i = 0; i < 3; i++) {
	        s = strtok (NULL, WS);
		res = sscanf (s, "%lf", &disp);
		if (res)
		    Disp(nd,i) = disp;
	    }
	}

	ifs.close();
    }

    // Initialise system matrix
    RCompRowMatrix F;
    SetupSysmat (mesh1, F);
    RVector rsize1(nreg), rsize2(nreg);
    
    if (perc1)
    {
       RegsizePerc (mesh1, percs1, matidx1, rsize1);
       cout << "perc1" << endl;
    }
    else
    {
       Regsize (mesh1, matidx1, rsize1);  // region sizes of original mesh
    }
    cout << "initial rsize1 = " << rsize1 << endl;
    
    if (diffMat)
    {
      if (perc2)
      {
        RegsizePerc (mesh2, percs2, matidx2, rsize2);
	cout << "diffMat, perc2" << endl;
      }
      else
      {
        Regsize (mesh2, matidx2, rsize2);  // region sizes of distorted mesh
	cout << "diffMat, not perc2" << endl;
      }
    }
    else
    {
      if (perc1)
      {
        RegsizePerc (mesh2, percs1, matidx1, rsize2);
	cout << "no diffMat, perc1" << endl;
      }
      else
      {
        Regsize (mesh2, matidx1, rsize2);
      }
    }
    cout << "initial rsize2 = " << rsize2 << endl;  


#if METHOD == 1
    // This method uses the region size as data

    // Now build the Jacobian regsize_i / thermal coeff_j
    double of, dth = 1e-6;
    time_t t_start,t_end;
    RDenseMatrix J (nreg, nreg);
    RVector x, rsize(nreg);
    RVector te_tot(nreg);
    for (i = 0; i < nreg; i++) {
	te_tot[i] = mat1[i].te;
	mat1[i].te = 0.0;
    }
    cout << "initial region size difference = " << rsize2-rsize1 << endl;
    double error_tmp = 0.0, total_error = 0.0, min_error = 1e-3;
    for (i = 0; i < nreg; i++)
	{
	  error_tmp = rsize2[i]-rsize1[i];
	  total_error += error_tmp*error_tmp;
	}
    total_error = sqrt(total_error);
    cout << "total_error = " << total_error << endl;
    
    for (int iter = 0; iter < itmax; iter++) {

	for (i = 0; i < nreg; i++) {
	    mat1[i].te += dth;
	    cout << "region " << i << ", te = " << mat1[i].te << endl;
	    (void) time(&t_start);
	    Solve (mesh1, F, Disp, mat1, matidx1, nreg, x);
	    (void) time(&t_end);
	    cout << "time for solver = " << (int)(t_end-t_start) << endl;
	    Displace (mesh1, x);
	    
	    if (perc1)
	      RegsizePerc (mesh1, percs1, matidx1, rsize);
	    else
	      Regsize (mesh1, matidx1, rsize);
	    
	    cout << "rsize = " << rsize << endl;
	    for (j = 0; j < nreg; j++)
	    {
		J(j,i) = (rsize[j]-rsize1[j])/dth;
		cout << "region " << j << ", rsize-rsize1 = " << rsize[j]-rsize1[j] 
			<< ", J(" << j << "," << i << ") = " << J(j,i) << endl;
	    }
	    Displace (mesh1, -x); // undo displacements
	    mat1[i].te -= dth;
	}

	RVector dx(rsize2-rsize1); // region size difference
	cout << "region size difference (rsize2-rsize1) = " << rsize2-rsize1 << endl;

	// Lets see if this Jacobian is invertible
	IVector indx;
	double d;
	LUFactorize (J, indx, d);
	cout << "after LUFactorize" << endl;
	cout << "J = " << J << endl;
	cout << "indx = " << indx << endl;
	cout << "d = " << d << endl;
	LUSolve (J, indx, dx);
	cout << "after LUSolve" << endl;
	cout << "J = " << J << endl;
	cout << "indx = " << indx << endl;
	cout << "dx = " << d << endl;

	// update thermal coefficients
	for (i = 0; i < nreg; i++) {
	    te_tot[i] += dx[i];
	    cout << "dx = " << dx[i] << ", accumulated te_tot = " << te_tot[i] << endl;
	    //mat[i].te += dx[i];
	}
	// update mesh
	for (i = 0; i < nreg; i++)
	    mat1[i].te = dx[i];
	Solve (mesh1, F, Disp, mat1, matidx1, nreg, x);
	Displace (mesh1, x);
	mesh1.Setup();
	
	if (perc1)
	  RegsizePerc (mesh1, percs1, matidx1, rsize1);
	else
	  Regsize (mesh1, matidx1, rsize1);
	
	cout << "after update mesh, rsize1 = " << rsize1 << endl;
	for (i = 0; i < nreg; i++)
	    mat1[i].te = 0;
	    
	cout << "Iteration " << iter+1 << endl;
	error_tmp = 0.0;
	total_error = 0.0;
	for (i = 0; i < nreg; i++)
	{
	    cout << "Reg " << i << "\t te = " << te_tot[i]
		 << "\t size = " << rsize1[i] << ", diff = " << rsize2[i]-rsize1[i]  
		 << endl;
	  error_tmp = rsize2[i]-rsize1[i];
	  total_error += error_tmp*error_tmp;
	}
	total_error = sqrt(total_error);
	cout << "total_error = " << total_error << endl;
     
        if (total_error < min_error)
        {
           cout << "total error < 1e-3" << endl;
           break;
        }
   }
  cout << "normalized te's" << endl;
  cout.setf(ios::fixed);
  cout.precision(7);
  for (i = 0; i < nreg; i++)
  {
    cout << "region " << i << ", te = " << (mat1[i].te - (mat1[0].te/mat1[i].E)) << endl;
  }

    

#elif METHOD == 2
    // This method uses the individual node displacements as data

    // form Jacobian by explicit perturbation of region thermal coeffs
    const double dth = 1e-6;
    cerr << "Allocating Jacobian " << nlen*dim << " x " << nreg << endl;
    RDenseMatrix J (nlen*dim, nreg);
    cerr << "Finished" << endl;
    RVector x;
    RVector te_tot(nreg);
    for (i = 0; i < nreg; i++) {
	te_tot[i] = mat[i].te;
	mat[i].te = 0.0;
    }

    cout << "Initial error: " << of3 (mesh1, mesh2) << endl;
    WriteEimHeader (meshname, mesh1.elen(), "recon_te.eim", "N/A");

    for (int iter = 0; iter < 4; iter++) {

	cout << "Generating Jacobian:" << endl;
	for (i = 0; i < nreg; i++) {
	    cout << "region " << i << endl;
	    mat[i].te += dth;
	    Solve (mesh1, F, Disp, mat, matidx, nreg, x);
	    for (j = 0; j < nlen; j++) {
		for (k = 0; k < dim; k++)
		    J(j*dim+k,i) = x[j*dim+k]/dth;
	    }
	    mat[i].te -= dth;
	    //cout << "Reg. " << i+1 << " of " << nreg << endl
	//	 << "\033[1A";
	}
	// the data: nodal displacements
	RVector b(nlen*dim);
	for (j = 0; j < nlen; j++)
	    for (k = 0; k < dim; k++)
		b[j*dim+k] = mesh2.nlist[j][k] - mesh1.nlist[j][k];

	// gradient
	RVector g(nreg);
	g = ATx (J, b);

	// overdetermined solver	
	double lambda = 1e-3;
	RSymMatrix JTJ = ATA(J);
	RVector diag = JTJ.Diag();
	double mn = mean(diag);
	for (i = 0; i < JTJ.nRows(); i++)
	    JTJ(i,i) += lambda*mn;
	CHdecomp (JTJ, true);
	RVector dx = CHsubst (JTJ, g);

	// TEMPORARY: random scaling
	//dx *= 10.0;
	
	// update thermal coefficients
	for (i = 0; i < nreg; i++) {
	    te_tot[i] += dx[i];
	    //mat[i].te += dx[i];
	}
	// update mesh
	for (i = 0; i < nreg; i++)
	    mat[i].te = dx[i];
	Solve (mesh1, F, Disp, mat, matidx, nreg, x);
	for (j = 0; j < nlen; j++)
	    for (k = 0; k < dim; k++)
		mesh1.nlist[j][k] += x[j*dim+k];
	mesh1.Setup();
	Regsize (mesh1, matidx, rgsize);
	for (i = 0; i < nreg; i++)
	    mat[i].te = 0;

	cout << "Iteration " << iter+1 << endl;
	for (i = 0; i < nreg; i++)
	    cout << "Reg " << i << "\t te = " << te_tot[i]
		 << "\t size = " << rgsize[i]
		 << endl;
	cout << "Error: " << of3 (mesh1, mesh2) << endl;
	cout << endl;

	// write out the elementwise thermal coeffs
	RVector eim(mesh1.elen());
	for (i = 0; i < mesh1.elen(); i++)
	    eim[i] = te_tot[matidx[i]];
	WriteEim (eim, iter, "recon_te.eim");
    }

#elif METHOD == 3

    // this method uses correspondence of sampled images for
    // objective function

    const int gridsize = 256;
    const double dth = 1e-1;

    Point bbmin, bbmax;
    mesh2.BoundingBox (bbmin, bbmax);
    // note that the bounding box should be identical for both meshes
    // if we keep the boundary nodes fixed.
    RVector bbsize = bbmax - bbmin;
    double maxsize = vmax(bbsize);
    IVector gdim(bbsize.Dim());
    int glen = 1;
    for (i = 0; i < gdim.Dim(); i++) {
	gdim[i] = (int)(bbsize[i] * gridsize / maxsize);
	glen *= gdim[i];
    }
    cerr << "Mesh sampling grid: " << gdim << endl;

    // raster the target mesh (mesh2)
    int *elref2 = GenerateElementPixelRef (mesh2, gdim, &bbmin, &bbmax);

#ifdef DEBUG_OUTPUT
    // DEBUG OUTPUT
    ofstream dofs ("dbg.dat");
    int val, ofs = gdim[0]*gdim[1]*(gdim[2]/2);
    for (i = 0; i < gdim[0]; i++)
	for (j = 0; j < gdim[1]; j++) {
	    int el = elref2[i + j*gdim[0] + ofs];
	    if (el >= 0) val = matidx[el];
	    else val = -1;
	    dofs << val << endl;
	}
#endif

    cout << "after the methods, it continues to do a lot of things ... " << endl;
    
    // form Jacobian by explicit perturbation of region thermal coeffs
    RDenseMatrix J (nreg, nreg);
    RVector x;
    RVector te_tot(nreg);
    for (i = 0; i < nreg; i++) {
	te_tot[i] = mat[i].te;
	mat[i].te = 0.0;
    }

    int *elref1orig =  GenerateElementPixelRef (mesh1, gdim, &bbmin, &bbmax);
    cout << "Initial error: " << of4 (elref1orig, elref2, glen) << endl;
    double *corr = new double[nreg];
    double *corr2 = new double[nreg];
    for (i = 0; i < nreg; i++)
	corr[i] = RegCorr (elref1orig, elref2, glen, i);

    for (int iter = 0; iter < 10; iter++) {

	for (i = 0; i < nreg; i++) {
	    mat[i].te += dth;
	    Solve (mesh1, F, Disp, mat, matidx, nreg, x);
	    Displace (mesh1, x);
	    int *elref1 = GenerateElementPixelRef (mesh1, gdim, &bbmin,&bbmax);

	    // verbose output
	    int md;
	    for (j = md = 0; j < glen; j++)
		if (elref1[j] != elref1orig[j]) md++;
	    cerr << md << " of " << glen << " voxels reassigned" << endl;

	    for (j = 0; j < nreg; j++) {
		corr2[j] = RegCorr (elref1, elref2, glen, j);
		cerr << "reg " << j << ": correspondence orig=" << corr[j]
		     << ", pert=" << corr2[j] << endl;
	    }
	    delete []elref1;
	    for (j = 0; j < nreg; j++)
		J(j,i) = (corr[j]-corr2[j])/dth;
	    Displace (mesh1, -x);
	    mat[i].te -= dth;
	}

	cerr << "Jacobian:\n" << J << endl;

	// the data: pixel correspondences
	RVector b(nreg);
	for (i = 0; i < nreg; i++) b[i] = corr[i];

	// gradient
	RVector g(nreg);
	g = ATx (J, b);

	// overdetermined solver	
	double lambda = 1e-5;
	RSymMatrix JTJ = ATA(J);
	RVector diag = JTJ.Diag();
	double mn = mean(diag);
	for (i = 0; i < JTJ.nRows(); i++)
	    JTJ(i,i) += lambda*mn;
	CHdecomp (JTJ, true);
	RVector dx = CHsubst (JTJ, g);

	// update thermal coefficients
	for (i = 0; i < nreg; i++) {
	    te_tot[i] += dx[i];
	    //mat[i].te += dx[i];
	}
	// update mesh
	for (i = 0; i < nreg; i++)
	    mat[i].te = dx[i];
	Solve (mesh1, F, Disp, mat, matidx, nreg, x);
	for (j = 0; j < nlen; j++)
	    for (k = 0; k < dim; k++)
		mesh1.nlist[j][k] += x[j*dim+k];
	mesh1.Setup();
	for (i = 0; i < nreg; i++)
	    mat[i].te = 0;

	cout << "Iteration " << iter+1 << endl;
	for (i = 0; i < nreg; i++)
	    cout << "Reg " << i << "\t te =" << te_tot[i] << endl;

	int *elref1 = GenerateElementPixelRef (mesh1, gdim, &bbmin,&bbmax);
	cout << "Error: " << of4 (elref1, elref2, glen) << endl;
	delete []elref1;
	cout << endl;
    }
    delete []elref1orig;
    delete []elref2;

#endif


    delete []mat1;
    delete []mat2;
    delete []matidx1;
    delete []matidx2;
    delete []nelreg1;
    delete []nelreg2;
    for (j = 0; j < 2; j++)
    {
	delete []regsize1[j];
	delete []regsize2[j];
    }
    return 0;



#ifdef UNDEF
    double tol = 1e-12;
    char bname[256], dname[256], cbuf[256], *s;
    double pratio0, modulus0, te0, disp, sz, tot_sz;
    Mesh mesh;
    int ex, ey, ez;
    double dx, dy, dz;
    bool *egrid;
    int i, j, k, nd, res, cmd, el;
    int *rowptr, *colidx, nzero;
    int *drowptr, *dcolidx, dn;
    int encoding;
    bool bDisp, bBndf, bAuto;
    ifstream disp_f, bndf_f;

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
    int dim = mesh.Dimension();
    cout << "Total mesh size before displacement: " << mesh.FullSize() << endl;

    // record element volumes
    RVector size0(mesh.elen());
    for (i = 0; i < mesh.elen(); i++) {
	size0[i] = mesh.elist[i]->Size();
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
    cout << "Element parameters  (1) homogeneous  (2) from file: ";
    cin >> cmd;
    if (cmd == 1) {
        cout << "Poisson's ratio: ";
	cin  >> pratio0;
	cout << "Young's modulus: ";
	cin  >> modulus0;
	cout << "Thermal expansion: ";
	cin  >> te0;
	pratio = pratio0;
	modulus = modulus0;
	te = te0;
    } else {
        int nmat, nidx, idx;
        cout << "Material file: ";
	cin >> cbuf;
	ifstream ifs (cbuf);
	ifs >> nmat;
	struct MatList { double E, nu, dns, te; } *matlist = new MatList[nmat];
	for (i = 0; i < nmat; i++) {
	    ifs >> matlist[i].E >> matlist[i].nu >> matlist[i].dns
		>> matlist[i].te;
	}
	ifs >> nidx;
	if (nidx != mesh.elen()) {
	    cerr << "Invalid length of material index list. Aborting." << endl;
	    exit(1);
	}

	// temporary: map material index to a nim file
	RVector ndmat(mesh.nlen());
	IVector nndmat(mesh.nlen());

	for (i = 0; i < nidx; i++) {
	    ifs >> idx;
	    if (idx < 0 || idx >=nmat) {
	        cerr << "Invalid index in material index list. Aborting." << endl;
		exit (1);
	    }
	    pratio[i] = matlist[idx].nu;
	    modulus[i] = matlist[idx].E;
	    te[i] = matlist[idx].te;

	    // map material index from element to node
	    Element *pel = mesh.elist[i];
	    for (j = 0; j < pel->nNode(); j++) {
		ndmat[pel->Node[j]] += idx;
		nndmat[pel->Node[j]]++;
	    }

	}
	for (i = 0; i < ndmat.Dim(); i++)
	    ndmat[i] /= nndmat[i];
	WriteNimHeader (bname, mesh.nlen(), "material.nim", "N/A");
	WriteNim (ndmat, 1, "material.nim");
    }
    cout << "Assembling system matrix" << endl;
    AddToSysMatrix_elasticity (mesh, F, modulus, pratio);
    cout << "Finished!" << endl;

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

    // add thermal expansion
    double dT;
    cout << "Applied temperature change [K]: ";
    cin >> dT;
    AddToRHS_thermal_expansion (mesh, rhs, modulus, pratio, te, dT);

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
        for (i = 0; i < n; i++) {
	    for (j = 0; j < 3; j++)
	        mesh.nlist[i][j] += x[i*3+j];
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

    } else if (encoding == 2) {

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

    } else if (encoding == 3) {
	
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
	for (nd = 0; nd < mesh.nlen(); nd++)
	    mesh.plist[nd].SetMua (ds[nd]/count[nd]);

    }

    cout << "Total mesh size after displacement: " << mesh.FullSize() << endl;

    // Write out mesh
    ofstream ofs ("displace.msh");
    ofs << mesh << endl;

    return 0;

#endif // UNDEF
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

#ifdef UNDEF
double LineSearch (const Mesh &mesh, const Mesh &mesh_tgt, RVector &p,
    RVector &x, MatList *mat, RCompRowMatrix &F, RDenseMatrix &Disp,
    int *matidx)
{
    Mesh mesh_tmp(mesh);
    int i, j, nlen = mesh.nlen(), nreg = p.Dim();
    int dim = mesh.Dimension();
    double te_store[nreg];
    for (i = 0; i < nreg; i++) te_store[i] = mat[i].te;

    double f0 = of3 (mesh, mesh_tgt);

    double alpha = 1.0;
    for (i = 0; i < nreg; i++)
	mat[i].te = p*alpha;
    Solve(mesh_tmp, F, Disp, mat, matidx, nreg, x);
    for (j = 0; j < nlen; j++)
	for (k = 0; k < dim; k++)
	    mesh_tmp.nlist[j][k] += 
#endif
