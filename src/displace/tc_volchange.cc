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

void Regsize (const Mesh &mesh, int *matidx, RVector &grpsize, int *volchanges_state, int nreg)
{
  // returns the size of each region in the mesh, as defined by the
  // matidx index list
  
  grpsize = 0.0; // reset
  int i, count = 0, reg, elen = mesh.elen();
  RVector grpsize_tmp(nreg);
  for (i = 0; i < elen; i++) {
    reg = matidx[i];
    if (reg >= nreg) {
      cerr << "Index out of range" << endl;
      exit(0);
    }
    grpsize_tmp[reg] += mesh.ElSize(i);
  }
  for(i = 0; i < nreg; i++)
    {
      if (volchanges_state[i] != 0)
	{
	  grpsize[count] = grpsize_tmp[i];
	  count++;
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
    {
      if (disp(i,j) != DISP_NONE) {
	fd = F(i*dim+j, i*dim+j) *= BIGSPRING;
	rhs[i*3+j] = fd * disp(i,j);
      }
     }
  
  
  RPrecon_Diag precon;
  //RPrecon_IC precon;
  //RPrecon_Null precon;
  precon.Reset (&F);
  x.New(dn);
  //niter = GMRES (F, rhs, x, tol, &precon);
  //niter = BiCGSTAB (F, rhs, x, tol, &precon, dn/10);
  niter = BiCGSTAB (F, rhs, x, tol, &precon, dn/10);
  //niter = PCG (F, rhs, x, tol, &precon);
}

int main (void) {
  
  int i, j, k, nlen, elen, el, dim, nreg1, nreg2, *matidx1, res, nd, itmax = 4;
  double *volchanges;
  char *s;
  int *volchanges_state;
  char fname[256], mname1[256], volchangename[256], matvolchangename[256], boundary_name[256], dname[256], cbuf[256];
  char meshname[256];
  ifstream ifs;
  ofstream ofs;
  Mesh mesh1;
  int *rowptr, *colidx, nzero;
  int *drowptr, *dcolidx, dn;
  
  MatList *mat1;
  int count = 0;
  ifstream disp_f, bndf_f;

  cout << "Original mesh name: ";
  cin >> meshname;
  ifs.open (meshname);
  ifs >> mesh1;
  ifs.close();
  mesh1.Setup();
  
  cout << "boundary condition file name: ";
  cin >> boundary_name;
  disp_f.open (boundary_name);
  
  cout << "Material/region original mesh name: ";
  cin >> mname1;
  
  cout << "Desired volume change file name: ";
  cin >> volchangename;
  
  cout << "Material/region final mesh name: ";
  cin >> matvolchangename;
  
  // sanity checks
  nlen = mesh1.nlen();
  elen = mesh1.elen();
  dim  = mesh1.Dimension();
  
  cout << "Mesh 1 size: " << mesh1.FullSize() << endl;
  
  RVector size0(elen);
  for (i = 0; i < elen; i++)
    size0[i] = mesh1.ElSize(i);
  
  // Initialise system matrix
  RCompRowMatrix F;
  SetupSysmat (mesh1, F);
  
  // scan original mesh regions
  ifs.open (mname1);
  ifs >> nreg1;
  cout << "original mesh material file, nreg = " << nreg1 << endl;
  mat1 = new MatList[nreg1];
  for (i = 0; i < nreg1; i++) {
    ifs >> mat1[i].E >> mat1[i].nu >> mat1[i].dns >> mat1[i].te;
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
  }
  cout << "Found " << nreg1 << " regions:\n";
  ifs.close();
  
  
  // scan displacement file
  RDenseMatrix Disp(nlen, 3);
  Disp = DISP_NONE;
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
  
  // scan desired volume changes for each region
  ifs.open (volchangename);
  ifs >> nreg2;
  cout << "desired volchange file, nreg = " << nreg2 << endl;
  if (nreg1 != nreg2)
    {
      cout << "regions of material and desired volchange files are different" << endl;
      exit(1);
    }
  
  volchanges =  new double[nreg1];
  volchanges_state = new int[nreg1];
  int nactive_materials = 0;
  for (i = 0; i < nreg1; i++)
    {
      volchanges[i] = 0;
      volchanges_state[i] = 0;
    }
  
  for (i = 0; i < nreg2; i++) {
    ifs >> volchanges[i];
    if (volchanges[i] != 999)
      {
	volchanges_state[i] = 1;
	nactive_materials++;
	cout << "region " << i << " is active, volchange = " << volchanges[i] << endl;
      }
  }
  
  RVector rsize1(nactive_materials), rsize2(nactive_materials);
  
  Regsize (mesh1, matidx1, rsize1, volchanges_state, nreg1);  // region sizes of original mesh
  
  count = 0;
  for (i = 0; i < nreg1; i++)
    {
      if (volchanges_state[i] == 1)
	{
	  rsize2[count] = rsize1[count] + (rsize1[count] * (volchanges[i]/100.0));
	  cout << "region " << i << ", rsize1 = " << rsize1[count] << ", volchange = " << volchanges[i] 
	       << ", rsize2 = " << rsize2[count] << endl;
	  count++;
	}
    }
  ifs.close();
  
  cout << "initial rsize1 = " << rsize1 << endl;
  cout << "final rsize2 = " << rsize2 << endl;
  
  // Now build the Jacobian regsize_i / thermal coeff_j
  double of, dth = 1e-6;
  time_t t_start,t_end;
  RDenseMatrix J (nactive_materials, nactive_materials);
  RVector x, rsize(nactive_materials);
  RVector te_tot(nreg1);
  for (i = 0; i < nreg1; i++) {
    te_tot[i] = mat1[i].te;
    mat1[i].te = 0.0;
  }
  cout << "initial region size difference = " << rsize2-rsize1 << endl;
  double error_tmp = 0.0, total_error = 0.0, min_error = 1e-3;
  for (i = 0; i < nactive_materials; i++)
    {
      error_tmp = rsize2[i]-rsize1[i];
      total_error += error_tmp*error_tmp;
    }
  total_error = sqrt(total_error);
  cout << "total_error = " << total_error << endl;
  
  for (int iter = 0; iter < itmax; iter++) {
    
    count = 0;
    for (i = 0; i < nreg1; i++) {
      if (volchanges_state[i] != 0)
	{
	  
	  mat1[i].te += dth;
	  cout << "region " << i << ", te = " << mat1[i].te << endl;
	  (void) time(&t_start);
	  Solve (mesh1, F, Disp, mat1, matidx1, nreg1, x);
	  (void) time(&t_end);
	  cout << "time for solver = " << (int)(t_end-t_start) << endl;
	  Displace (mesh1, x);
	  
	  Regsize (mesh1, matidx1, rsize, volchanges_state, nreg1);
	  
	  cout << "rsize = " << rsize << endl;
	  for (j = 0; j < nactive_materials; j++)
	    {
	      J(j,count) = (rsize[j]-rsize1[j])/dth;
	      cout << "region " << j << ", rsize-rsize1 = " << rsize[j]-rsize1[j] 
		   << ", J(" << j << "," << count << ") = " << J(j,count) << endl;
	    }
	  Displace (mesh1, -x); // undo displacements
	  mat1[i].te -= dth;
	  count++;
	}
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
    cout << "dx = " << dx << endl;
    
    // update thermal coefficients
    count = 0;
    for (i = 0; i < nreg1; i++) {
      cout << "reg = " << i << ", te_tot = " << te_tot[i] << endl;
      if (volchanges_state[i] != 0)
	{
	  te_tot[i] += dx[count];
	  cout << "dx = " << dx[count] << ", accumulated te_tot = " << te_tot[i] << endl;
	  count++;
	}
    }
    // update mesh
    cout << "updating mesh" << endl;
    count = 0;
    for (i = 0; i < nreg1; i++)
      {
	    if (volchanges_state[i] != 0)
	    {
	      mat1[i].te = dx[count];
	      count++;
	    }
	    cout << "reg " << i << ", mat1.te = " << mat1[i].te << endl;
      }
    Solve (mesh1, F, Disp, mat1, matidx1, nreg1, x);
    Displace (mesh1, x);
    mesh1.Setup();
    
    Regsize (mesh1, matidx1, rsize1, volchanges_state, nreg1);
    
    cout << "after update mesh, rsize1 = " << rsize1 << endl;
    for (i = 0; i < nreg1; i++)
      mat1[i].te = 0;
    
    cout << "Iteration " << iter+1 << endl;
    error_tmp = 0.0;
    total_error = 0.0;
    count = 0;
    for (i = 0; i < nreg1; i++)
      {
	cout << "Reg " << i << "\t te = " << te_tot[i] << endl;
	if (volchanges_state[i] == 1)
	  {
	    cout << "\t size = " << rsize1[count] << ", diff = " << rsize2[count]-rsize1[count]  
		 << endl;
	    error_tmp = rsize2[count]-rsize1[count];
	    total_error += error_tmp*error_tmp;
	    count++;
	  }
      }
    total_error = sqrt(total_error);
    cout << "total_error = " << total_error << endl;
    if (total_error < min_error)
    {
      cout << "total error < 1e-3" << endl;
      break;
    }
  }//end of iterations
  
  ofs.open (matvolchangename);
  
  ofs << nreg1 << endl;
  for (i = 0; i < nreg1; i++) 
    ofs << mat1[i].E << " " << mat1[i].nu << " "  << mat1[i].dns << " "  << te_tot[i] << endl;;
  
  ofs << elen << endl;
  
  for (i = 0; i < elen; i++) 
    ofs << matidx1[i] << endl;
  
  ofs.close(); 
  
  delete []mat1;
  delete []matidx1;
  
  return 0;
}                                                                              
