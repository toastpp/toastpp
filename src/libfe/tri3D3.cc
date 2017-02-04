// ==========================================================================
// Module libfe
// File tri3D3.cc
// Definition of class Triangle3D3
// ==========================================================================

#define FELIB_IMPLEMENTATION

#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "mathlib.h"
#include "felib.h"
#include "lin_qr.h"
#include "tri_qr.h"

using namespace std;

//#define DEBUG_TRI3D3
/* utility for crossproducts */
FELIB RVector cross3(RVector& a, RVector& b){
    RVector c(3); 
    c[0]  = (a[1]*b[2] - a[2]*b[1]);
    c[1]  = (a[2]*b[0] - a[0]*b[2]);
    c[2]  = (a[0]*b[1] - a[1]*b[0]);
    return c;
}

Triangle3D3::Triangle3D3 () : Triangle3()
{
}

Triangle3D3::Triangle3D3 (const Triangle3D3& el): Triangle3(el)
{
}

Triangle3D3::~Triangle3D3 ()
{
  delete [] V2D;
}

void Triangle3D3::Initialise (const NodeList& nlist)
{
    //    cout << "Initialising Triangle3D3 element \n";
    V2D = new RVector [3];
    for (int iv = 0; iv < 3 ; iv ++)
      V2D[iv].New(2);
    // create  mapped points    
    V2D[0][0] = 0; V2D[0][1] = 0;

    RVector r10 = nlist[Node[1]] - nlist[Node[0]];
    RVector r20 = nlist[Node[2]] - nlist[Node[0]];
    double d10 = l2norm(r10);
    double d20 = l2norm(r20);
    double cs  = (r10 & r20)/(d10*d20);
    double ss  = sqrt(1.0 - (cs*cs));

    V2D[1][0] = d10; V2D[1][1] = 0;
    V2D[2][0] = d20*cs; V2D[2][1] = d20*ss;
#ifdef DEBUG_TRI3D3 
    cout << "length of edges " << d10 << "," << d20 << endl;
    cout << "angle 2-1-3 : cos " << cs << " sin " << ss << endl;
    cout << "2D points : " << V2D[0] << " " << V2D[1] << " " << V2D[2] << endl;
#endif 
    // set up mappings 
    //

    //    Map2Dto3D.New(3,2); // not needed for now...
    Map3Dto2D.New(2,3);
    MapLocaltoGlobal.New(3,3);

    RVector nu(3);
    nu = cross3(r10,r20); // normal to plane
    nu /= l2norm(nu);       // normalise 
    //    cout << r10 << " " << r20 << " " << nu << endl;
    // local to global map - trivial 

    for(int i = 0; i < 3 ; i++) {
      MapLocaltoGlobal(i,0) = r10[i];
      MapLocaltoGlobal(i,1) = r20[i];
      MapLocaltoGlobal(i,2) = nu[i];
    }

    //    cout << "local to global \n" <<  MapLocaltoGlobal << endl;
    // other mappings

    RVector r1p(3), r2p(3); // normals IN plane
    r1p = cross3(r10,nu);
    r2p = cross3(r20,nu);

    RDenseMatrix B(2,3); // temporary. Actually it is global to local map.
    double xi = r2p & r10;
    double eta = r1p & r20;
    B(0,0) = r2p[0]/xi;  B(0,1) = r2p[1]/xi;  B(0,2) = r2p[2]/xi;
    B(1,0) = r1p[0]/eta; B(1,1) = r1p[1]/eta; B(1,2) = r1p[2]/eta;
    RDenseMatrix M(2,2);
    M(0,0) = V2D[1][0]; M(1,0) = V2D[1][1];
    M(0,1) = V2D[2][0]; M(1,1) = V2D[2][1];

    //    cout << "local coordinate normalisations : " << xi << "," << eta << endl;
    //    cout << "M : \n" << M << endl;
    //    cout << "B : \n" << B << endl;

    Map3Dto2D = M * B;
#ifdef DEBUG_TRI3D3 
    cout << "3D to 2D mapping : \n" << Map3Dto2D << endl;
#endif

    // Set up triangle geometry parameters
    double x0 = V2D[0][0], y0 = V2D[0][1];
    double x1 = V2D[1][0], y1 = V2D[1][1];
    double x2 = V2D[2][0], y2 = V2D[2][1];

    a0 = x1*y2 - x2*y1;  b0 = y1-y2;  c0 = x2-x1;
    a1 = x2*y0 - x0*y2;  b1 = y2-y0;  c1 = x0-x2;
    a2 = x0*y1 - x1*y0;  b2 = y0-y1;  c2 = x1-x0;

    // I really don't know if this next line works ...
    Element_Unstructured_2D::Initialise (nlist);
    dASSERT(size > 0, "Element size not positive");

#ifdef TRI3_STORE_INTFF
    intff.New(3);
    intff = sym_intff * size;
#endif

#ifdef TRI3_STORE_INTDD
    intdd.New(3);
    intdd = ComputeIntDD (nlist);
#endif
    /* don't know about this either */
    //    if (!subsampling_initialised) {
    //    Point loc(2);
    //	double bloc;
    //	int i, j;
    //    for (i = 0; i < nsample_lin; i++) {
    //	    bloc = (i+0.5)/(double)nsample_lin;
    //	    for (j = 0; j < 3; j++)
    //		absc_bndsample[j][i].New(2);
    //	    absc_bndsample[0][i][0] = bloc;
    //	    absc_bndsample[0][i][1] = 0.0;
    //	    absc_bndsample[1][i][0] = bloc;
    //	    absc_bndsample[1][i][1] = bloc;
    //	    absc_bndsample[2][i][0] = 0.0;
    //	    absc_bndsample[2][i][1] = bloc;
    //	}
    //	subsampling_initialised = true;
    //    }

  // lastly, push nodes to unit sphere. Have to add three edge points
  for (int i = 0; i < 3 ; i++) {
    uspts[i] = nlist[Node[i]];
    uspts[i+3] = 0.5*(nlist[Node[i]]+nlist[Node[(i+1)%3]]); // "edge points"
    uspts[i] /= l2norm(uspts[i]);
    uspts[i+3] /= l2norm(uspts[i+3]);
  }
}


Point Triangle3D3::Local (const NodeList &nlist, const Point& glob3D) const
{
    dASSERT(glob3D.Dim() == 3, "Argument 1 invalid dimension");
    RVector glob(2);
    glob = Map3Dto2D * (glob3D - nlist[Node[0]]);

    Point loc(2);
    double scale = 1.0/(a0+a1+a2);
    loc[0] = (a1 + b1*glob[0] + c1*glob[1]) * scale;
    loc[1] = (a2 + b2*glob[0] + c2*glob[1]) * scale;

    return loc;
}

Point Triangle3D3::SurfToLocal (int side, const Point &p) const
{
  cerr << "Triangle3D3::SurfToLocal not implemented \n";
  /* to be done */
  //    return Triangle_SurfToLocal (side, p);
  return Point();
}

RVector Triangle3D3::DirectionCosine (int side, RDenseMatrix &jacin)
{
  cerr << "Triangle3D3::DirectionCosine not implemented \n";
  /* to be done */
  //    RVector cosin(2);
  //  switch(side) {
  //  case 0: cosin[0] = -b2, cosin[1] = -c2; break;
  //  case 1: cosin[0] = -b0, cosin[1] = -c0; break;
  //  case 2: cosin[0] = -b1, cosin[1] = -c1; break;
  //  default: xERROR(Side index out of range);
  //  }
  //  return cosin/length(cosin);
	return RVector();
};



bool Triangle3D3::GContains (const Point& glob3D, const NodeList& nlist) const
{
  /* done */
    dASSERT(glob3D.Dim() == 3, "Argument 1 invalid dimension");
    RVector glob(2);
    glob = Map3Dto2D * (glob3D - nlist[Node[0]]);

    double xx = glob[0], yy = glob[1];

    // check bounding box
    if (xx < bbmin[0] || xx > bbmax[0] || yy < bbmin[1] || yy > bbmax[1])
        return false;

    double x0, x1, x2, y0, y1, y2, y0r, yyr, fac;
    const double EPS = 1e-10;

    for (int i = 0; i < 3; i++) {
        x0 = V2D[i][0],       y0 = V2D[i][1];
	x1 = V2D[(i+1)%3][0], y1 = V2D[(i+1)%3][1];
	x2 = V2D[(i+2)%3][0], y2 = V2D[(i+2)%3][1];
	if (fabs (x1-x2) < EPS) {
	    if ((x0 < x1 && xx > x1) || (x0 > x1 && xx < x1)) return false;
	} else {
	    fac = (y2-y1)/(x2-x1);
	    y0r = (x0-x1)*fac + y1;
	    yyr = (xx-x1)*fac + y1;
	    if ((y0 < y0r && yy > yyr) || (y0 > y0r && yy < yyr)) return false;
	}
    }
    return true;
}

RVector Triangle3D3::GlobalShapeF (const NodeList& nlist, const Point& glob3D)
    const
{
    dASSERT(glob3D.Dim() == 3, "Invalid point dimension");
    RVector glob(2);
    glob = Map3Dto2D * (glob3D - nlist[Node[0]]);
#ifdef DEBUG_TRI3D3 
    cout << "input point " << glob3D << " locally mapped point : " << glob << endl;
#endif
    RVector fun(3);
    double scale = 1.0/(2.0*size);
    fun[0] = scale * (a0 + b0*glob[0] + c0*glob[1]);
    fun[1] = scale * (a1 + b1*glob[0] + c1*glob[1]);
    fun[2] = scale * (a2 + b2*glob[0] + c2*glob[1]);
    return fun;
}

RDenseMatrix Triangle3D3::GlobalShapeD (const NodeList &nlist,
    const Point &glob3D) const
{
    dASSERT(glob3D.Dim() == 3, "Argument 1 invalid dimension");
    RVector glob(2);
    glob = Map3Dto2D * (glob3D - nlist[Node[0]]);

    // mapping seems uneccessary, since function is always constant
    RDenseMatrix der (2,3);
    double scale = 1.0/(2.0*size);
    der(0,0) = b0*scale;
    der(0,1) = b1*scale;
    der(0,2) = b2*scale;
    der(1,0) = c0*scale;
    der(1,1) = c1*scale;
    der(1,2) = c2*scale;
    return der;
}

RDenseMatrix  Triangle3D3::FTAMat () const 
{
    RDenseMatrix A(2,3);
    for(int i = 0; i < 2; i++)
      for(int j = 0; j < 3; j++)
	A(i,j) = MapLocaltoGlobal(j,i);
    return A;  
}


//------------- integrate functions on unit sphere --------------

RVector Triangle3D3::LocalShapeQF (const Point &loc) const
{
    dASSERT(loc.Dim() == 2, "Invalid point dimension");
    RVector fun(6);
    double L0 = 1.0-loc[0]-loc[1];
    double L1 = loc[0];
    double L2 = loc[1];
    fun[0] = L0 * (2.0*L0 - 1.0);
    fun[1] = L1 * (2.0*L1 - 1.0);
    fun[2] = L2 * (2.0*L2 - 1.0);
    fun[3] = 4.0*L0*L1;
    fun[4] = 4.0*L1*L2;
    fun[5] = 4.0*L0*L2;
    return fun;
}

RDenseMatrix Triangle3D3::LocalShapeQD (const Point &loc) const
{
    dASSERT(loc.Dim() == 2, "Invalid point dimension");
    RDenseMatrix der(2,6);
    double lx = loc[0], ly = loc[1];
    der(0,0) = der(1,0) = 4.0*(lx+ly) - 3.0;
    der(0,1) = 4.0*lx - 1.0;
    der(1,1) = 0.0;
    der(0,2) = 0.0;
    der(1,2) = 4.0*ly - 1.0;
    der(0,3) = 4.0*(1.0-2.0*lx-ly);
    der(1,3) = -4.0*lx;
    der(0,4) = 4.0*ly;
    der(1,4) = 4.0*lx;
    der(0,5) = -4.0*ly;
    der(1,5) = 4.0*(1.0-lx-2.0*ly);
    return der;
}
double Triangle3D3::IntUnitSphereFF (const NodeList& nlist, const int si, const int sj) const   
{ 
  // Returns integral of u_si(r)*u_sj(r) over unitsphere patch 
  // Since patch is quadratic, we have to invoke quadratic interpolation,
  // hence the above two functions (QF, QD).
  // The function f(r) is sampled at nodes into vector P. This can be three
  // noded or 6 noded.

  double result = 0.0;
  // quadrature integration
  for(int i = 0; i < 7 ; i++){
    Point2D eta;
    eta[0] = xg_1[i]; eta[1] =xg_2[i];
    RVector fun3 = LocalShapeF(eta);
    RDenseMatrix fd = LocalShapeQD(eta);
    // find Jacobian (length of cross product of x and y derivatives)
    RVector a(3), b(3);
    for (int j = 0; j < 6 ; j++) {
      for(int k = 0; k < 3; k++){
	a[k] += uspts[j][k] * fd(0,j);
	b[k] += uspts[j][k] * fd(1,j);
      }
    }
    double jac = l2norm(cross3(a,b));
    result += jac*fun3[si]*fun3[sj]*wg[i];
  }

  return result;
}
double Triangle3D3::IntUnitSpherePFF (const NodeList& nlist, const int si, const int sj, const RVector& P) const   
{ 
  // Returns integral of f(r)*u_si(r)*u_sj(r) over unitsphere patch 
  // Since patch is quadratic, we have to invoke quadratic interpolation,
  // hence the above two functions (QF, QD).
  // The function f(r) is sampled at nodes into vector P. This can be three
  // noded or 6 noded.
  // CORRECTION : 6-11-07. Global P samples, resampled onto Plocal. Can only be size 3.

  double result = 0.0;
  RVector Plocal(3);
  for (int i = 0; i < 3 ; i++)
    Plocal[i] = P[Node[i]];

  int pn = Plocal.Dim();
  // quadrature integration
  for(int i = 0; i < 7 ; i++){
    Point2D eta;
    eta[0] = xg_1[i]; eta[1] =xg_2[i];
    RVector fun3 = LocalShapeF(eta);
    RVector fun6 = LocalShapeQF(eta);
    RDenseMatrix fd = LocalShapeQD(eta);
    double pval = Plocal & (pn ==3 ? fun3 : fun6);   // interpolated value of f(eta)
    // find Jacobian (length of cross product of x and y derivatives)
    RVector a(3), b(3);
    for (int j = 0; j < 6 ; j++) {
      for(int k = 0; k < 3; k++){
	a[k] += uspts[j][k] * fd(0,j);
	b[k] += uspts[j][k] * fd(1,j);
      }
    }
    double jac = l2norm(cross3(a,b));
    result += jac*pval*fun3[si]*fun3[sj]*wg[i];
  }

  return result;
}
double Triangle3D3::IntUnitSphereP (const NodeList& nlist, const RVector& P) const   
{ 
  // Returns integral of f(r) over unitsphere patch 
  // Since patch is quadratic, we have to invoke quadratic interpolation,
  // hence the above two functions (QF, QD).
  // The function f(r) is sampled at nodes into vector P. This can be three
  // noded or 6 noded.
  // CORRECTION : 6-11-07. Global P samples, resampled onto Plocal. Can only be size 3.


  double result = 0.0;
  RVector Plocal(3);
  for (int i = 0; i < 3 ; i++)
    Plocal[i] = P[Node[i]];

  int pn = Plocal.Dim();
  // quadrature integration
  for(int i = 0; i < 7 ; i++){
    Point2D eta;
    eta[0] = xg_1[i]; eta[1] =xg_2[i];
    RVector fun = (pn == 3 ? LocalShapeF(eta) : LocalShapeQF(eta));
    RDenseMatrix fd = LocalShapeQD(eta);
    double pval = Plocal & fun;   // interpolated value of f(eta)
    // find Jacobian (length of cross product of x and y derivatives)
    RVector a(3), b(3);
    for (int j = 0; j < 6 ; j++) {
      for(int k = 0; k < 3; k++){
	a[k] += uspts[j][k] * fd(0,j);
	b[k] += uspts[j][k] * fd(1,j);
      }
    }
    double jac = l2norm(cross3(a,b));
    result += jac*pval*wg[i];
  }

  return result;
}









