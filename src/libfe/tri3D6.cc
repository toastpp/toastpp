// ==========================================================================
// Module libfe
// File tri3D6.cc
// Definition of class Triangle3D6
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

//#define DEBUG_TRI3D6
/* utility for crossproducts */
//RVector cross3(RVector& a, RVector& b){
//    RVector c(3); 
//    c[0]  = (a[1]*b[2] - a[2]*b[1]);
//    c[1]  = (a[2]*b[0] - a[0]*b[2]);
//    c[2]  = (a[0]*b[1] - a[1]*b[0]);
//    return c;
//}

Triangle3D6::Triangle3D6 () : Triangle6()
{
}

Triangle3D6::Triangle3D6 (const Triangle3D6& el): Triangle6(el)
{
}

Triangle3D6::~Triangle3D6 ()
{
  delete [] V2D;
}

void Triangle3D6::Initialise (const NodeList& nlist)
{
    //    cout << "Initialising Triangle3D6 element \n";
    V2D = new RVector [6];
    for (int iv = 0; iv < 6 ; iv ++)
      V2D[iv].New(2);
    // create  mapped points    
    // points 0 1 2 are identical to Triangle3D3
    V2D[0][0] = 0; V2D[0][1] = 0;

    RVector r10 = nlist[Node[1]] - nlist[Node[0]];
    RVector r20 = nlist[Node[2]] - nlist[Node[0]];
    double d10 = l2norm(r10);
    double d20 = l2norm(r20);
    double cs  = (r10 & r20)/(d10*d20);
    double ss  = sqrt(1.0 - (cs*cs));

    V2D[1][0] = d10; V2D[1][1] = 0;
    V2D[2][0] = d20*cs; V2D[2][1] = d20*ss;
#ifdef DEBUG_TRI3D6 
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

    // local to global map - trivial 

    for(int i = 0; i < 3 ; i++) {
      MapLocaltoGlobal(i,0) = r10[i];
      MapLocaltoGlobal(i,1) = r20[i];
      MapLocaltoGlobal(i,2) = nu[i];
    }

    // other mappings

    RVector r1p(3), r2p(3); // normals IN plane
    r1p = cross3(r10,nu);
    r2p = cross3(r20,nu);
    //    r1p[0] = (r10[1]*nu[2] - r10[2]*nu[1]);
    //    r1p[1] = (r10[2]*nu[0] - r10[0]*nu[2]);
    //    r1p[2] = (r10[0]*nu[1] - r10[1]*nu[0]);

    //    r2p[0] = (r20[1]*nu[2] - r20[2]*nu[1]);
    //    r2p[1] = (r20[2]*nu[0] - r20[0]*nu[2]);
    //    r2p[2] = (r20[0]*nu[1] - r20[1]*nu[0]);

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
#ifdef DEBUG_TRI3D6 
    cout << "3D to 2D mapping : \n" << Map3Dto2D << endl;
#endif
    // now get the other 2D vectors using mapping.
    for (int iv = 3 ; iv < 6; iv++)
      V2D[iv] = Map3Dto2D*(nlist[Node[iv]] - nlist[Node[0]]);

#ifdef DEBUG_TRI3D6 
    cout << "2D points : " << V2D[3] << " " << V2D[4] << " " << V2D[5] << endl;
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

#ifdef TRI6_STORE_INTFF
    intff.New(6);
    intff = sym_intff * size;
#endif

    //    if (!subsampling_initialised) {
    //        Point loc(2);
    //	int i, j, idx;
    //        for (i = idx = 0; i < nsample_lin; i++) {
    //	    loc[0] = (double)i/(double)(nsample_lin-1);
    //	    for (j = 0; j < nsample_lin-i; j++) {
    //	        loc[1] = (double)j/(double)(nsample_lin-1);
    //		absc_sample[idx].New(2);
    //		absc_sample[idx] = loc;
    //		idx++;
    //	    }
    //	}
    //	subsampling_initialised = true;
    //    }

  // lastly push nodes to unit sphere, for use by unit sphere integration.
    //    cout << "push to unit sphere" << endl;
  for (int i = 0; i < 6 ; i++) {
    uspts[i] = nlist[Node[i]] / l2norm(nlist[Node[i]]);
  }
}


Point Triangle3D6::Local (const NodeList &nlist, const Point& glob3D) const
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

Point Triangle3D6::SurfToLocal (int side, const Point &p) const
{
  cerr << "Triangle3D6::SurfToLocal not implemented \n";
  /* to be done */
  //    return Triangle_SurfToLocal (side, p);
	return Point();
}

RVector Triangle3D6::DirectionCosine (int side, RDenseMatrix &jacin)
{
  cerr << "Triangle3D6::DirectionCosine not implemented \n";
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



bool Triangle3D6::GContains (const Point& glob3D, const NodeList& nlist) const
{
  /* done */
    dASSERT(glob3D.Dim() == 3, "Argument 1 invalid dimension");
    RVector glob(2);
    glob = Map3Dto2D * (glob3D - nlist[Node[0]]);

    double xx = glob[0], yy = glob[1];

    // check bounding box
    //if (xx < bbmin[0] || xx > bbmax[0] || yy < bbmin[1] || yy > bbmax[1])
    //    return false;

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

RVector Triangle3D6::GlobalShapeF (const NodeList& nlist, const Point& glob3D)
    const
{
    dASSERT(glob3D.Dim() == 3, "Invalid point dimension");
    RVector glob(2);
    glob = Map3Dto2D * (glob3D - nlist[Node[0]]);
#ifdef DEBUG_TRI3D6 
    cout << "input point " << glob3D << " locally mapped point : " << glob << endl;
#endif

    RVector fun(6);
    double scale = 1.0/(2.0*size);
    double L0 = scale * (a0 + b0*glob[0] + c0*glob[1]);
    double L1 = scale * (a1 + b1*glob[0] + c1*glob[1]);
    double L2 = scale * (a2 + b2*glob[0] + c2*glob[1]);
    fun[0] = L0 * (2.0*L0 - 1.0);
    fun[1] = L1 * (2.0*L1 - 1.0);
    fun[2] = L2 * (2.0*L2 - 1.0);
    fun[3] = 4.0*L0*L1;
    fun[4] = 4.0*L1*L2;
    fun[5] = 4.0*L0*L2;
    return fun;

}

RDenseMatrix Triangle3D6::GlobalShapeD (const NodeList &nlist,
    const Point &glob3D) const
{
    dASSERT(glob3D.Dim() == 3, "Argument 1 invalid dimension");
    RVector glob(2);
    glob = Map3Dto2D * (glob3D - nlist[Node[0]]);

    RDenseMatrix der(2,6);
    double scale = 1.0/(2.0*size);
    double L0 = scale * (a0 + b0*glob[0] + c0*glob[1]);
    double L1 = scale * (a1 + b1*glob[0] + c1*glob[1]);
    double L2 = scale * (a2 + b2*glob[0] + c2*glob[1]);
    der(0,0) = b0*scale * (4.0*L0-1.0);
    der(1,0) = c0*scale * (4.0*L0-1.0);
    der(0,1) = b1*scale * (4.0*L1-1.0);
    der(1,1) = c1*scale * (4.0*L1-1.0);
    der(0,2) = b2*scale * (4.0*L2-1.0);
    der(1,2) = c2*scale * (4.0*L2-1.0);
    scale *= 4.0;
    der(0,3) = scale * (b0*L1 + b1*L0);
    der(1,3) = scale * (c0*L1 + c1*L0);
    der(0,4) = scale * (b1*L2 + b2*L1);
    der(1,4) = scale * (c1*L2 + c2*L1);
    der(0,5) = scale * (b0*L2 + b2*L0);
    der(1,5) = scale * (c0*L2 + c2*L0);
    return der;
}


RDenseMatrix  Triangle3D6::FTAMat () const 
{
    RDenseMatrix A(2,3);
    for(int i = 0; i < 2; i++)
      for(int j = 0; j < 3; j++)
	A(i,j) = MapLocaltoGlobal(j,i);
    return A;  
}
//------------- integrate functions on unit sphere --------------
double Triangle3D6::IntUnitSphereP (const NodeList& nlist, const RVector& P) const   
{ 
    // Returns integral of f(r) over unitsphere patch 

  double result = 0.0;
  RVector Plocal(6);
  for (int i = 0; i < 6 ; i++)
    Plocal[i] = P[Node[i]];
  // quadrature integration
  for(int i = 0; i < 7 ; i++){
    Point2D eta;
    eta[0] = xg_1[i]; eta[1] =xg_2[i];
    RVector fun = LocalShapeF(eta);
    RDenseMatrix fd = LocalShapeD(eta);
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
    //    cout << "eta : " << eta << " jac " << jac << endl;
    result += jac*pval*wg[i];
  }

  return result;
}
double Triangle3D6::IntUnitSphereFF (const NodeList& nlist, const int si, const int sj) const   
{ 
    // Returns integral of u_si(r) u_sj(r) over unitsphere patch 

  double result = 0.0;
  // quadrature integration
  for(int i = 0; i < 7 ; i++){
    Point2D eta;
    eta[0] = xg_1[i]; eta[1] =xg_2[i];
    RVector fun = LocalShapeF(eta);
    RDenseMatrix fd = LocalShapeD(eta);

    // find Jacobian (length of cross product of x and y derivatives)

    RVector a(3), b(3);
    for (int j = 0; j < 6 ; j++) {
      for(int k = 0; k < 3; k++){
	a[k] += uspts[j][k] * fd(0,j);
	b[k] += uspts[j][k] * fd(1,j);
      }
    }
    double jac = l2norm(cross3(a,b));
    //    cout << "eta : " << eta << " jac " << jac << endl;
    result += jac*fun[si]*fun[sj]*wg[i];
  }

  return result;
}
double Triangle3D6::IntUnitSpherePFF (const NodeList& nlist, const int si, const int sj, const RVector& P) const   
{ 
    // Returns integral of f(r) u_si(r) u_sj(r) over unitsphere patch 

  double result = 0.0;
  RVector Plocal(6);
  for (int i = 0; i < 6 ; i++)
    Plocal[i] = P[Node[i]];
  // quadrature integration
  for(int i = 0; i < 7 ; i++){
    Point2D eta;
    eta[0] = xg_1[i]; eta[1] =xg_2[i];
    RVector fun = LocalShapeF(eta);
    RDenseMatrix fd = LocalShapeD(eta);
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
    //    cout << "eta : " << eta << " jac " << jac << endl;
    result += jac*pval*fun[si]*fun[sj]*wg[i];
  }

  return result;
}




