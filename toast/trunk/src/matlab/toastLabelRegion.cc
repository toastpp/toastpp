// =========================================================================
// toastLabelRegion
// Given a mesh and a closed surface composed of triangles, this function
// assigns a region label to each mesh node inside the surface. Labels of
// external nodes remain unchanged.
//
// RH parameters:
//     1: mesh handle (pointer)
//     2: array of surface node coordinates (real n x 3)
//     3: triangle index list, 1-based (int m x 3)
//     4: internal label value (int)
// LH parameters:
//     1: array of mesh node labels
//
// Notes:
// * this currently only works for 3D meshes
// * the surface elements must be 3-noded triangles
// =========================================================================

#include "mex.h"
#include "mathlib.h"
#include "felib.h"
#include "toastmex.h"
#include "util.h"

// surface triangle parameters
typedef struct {
	double x1, y1, z1; // vtx 1
	double x2, y2, z2; // vtx 2
	double x3, y3, z3; // vtx 3
	double a, b, c, d; // plane params
	double d1, d2, d3; // dist of vtx i from opposite edge
} TriParam;

// =========================================================================
// set up the surface triangle parameters

TriParam *SetupSurface (const RDenseMatrix &vtx, const IDenseMatrix &idx)
{
    int i, j;
    int ntri = idx.nRows();
    TriParam *prm = new TriParam[ntri];

    // generate plane parameters E: ax+by+cz+d=0 for all triangles
    for (i = 0; i < ntri; i++) {
	TriParam &pp = prm[i];
	j = idx(i,0);
	pp.x1 = vtx(j,0), pp.y1 = vtx(j,1), pp.z1 = vtx(j,2);
	j = idx(i,1);
	pp.x2 = vtx(j,0), pp.y2 = vtx(j,1), pp.z2 = vtx(j,2);
	j = idx(i,2);
	pp.x3 = vtx(j,0), pp.y3 = vtx(j,1), pp.z3 = vtx(j,2);
	pp.a = pp.y1*(pp.z2-pp.z3) - pp.z1*(pp.y2-pp.y3) + 
	    (pp.y2*pp.z3 - pp.y3*pp.z2);
	pp.b = -(pp.x1*(pp.z2-pp.z3) - pp.z1*(pp.x2-pp.x3) +
	    (pp.x2*pp.z3 - pp.x3*pp.z2));
	pp.c = pp.x1*(pp.y2-pp.y3) - pp.y1*(pp.x2-pp.x3) +
	    (pp.x2*pp.y3 - pp.x3*pp.y2);
	pp.d = -(pp.x1*(pp.y2*pp.z3-pp.y3*pp.z2) - 
	    pp.y1*(pp.x2*pp.z3-pp.x3*pp.z2) + pp.z1*(pp.x2*pp.y3-pp.x3*pp.y2));

	RVector dr(3), a(3);
	dr[0] = pp.x1-pp.x2, dr[1] = pp.y1-pp.y2, dr[2] = pp.z1-pp.z2;
	a[0]  = pp.x3-pp.x2, a[1] =  pp.y3-pp.y2, a[2]  = pp.z3-pp.z2;
	pp.d1 = dot (dr, a)/l2normsq(a);
	dr[0] = pp.x2-pp.x3, dr[1] = pp.y2-pp.y3, dr[2] = pp.z2-pp.z3;
	a[0]  = pp.x1-pp.x3, a[1]  = pp.y1-pp.y3, a[2]  = pp.z1-pp.z3;
	pp.d2 = dot (dr, a)/l2normsq(a);
	dr[0] = pp.x3-pp.x1, dr[1] = pp.y3-pp.y1, dr[2] = pp.z3-pp.z1;
	a[0]  = pp.x2-pp.x1, a[1]  = pp.y2-pp.y1, a[2]  = pp.z2-pp.z1;
	pp.d3 = dot (dr, a)/l2normsq(a);
    }
    return prm;
}

// =========================================================================
// calculate the intersection of a line with the plane of a triangle

bool intersect (const RVector &pos, const RVector &dir, const TriParam &pp,
    RVector &X, double &t, bool half)
{
    const double EPS = 1e-8f;
    double den = pp.a*dir[0] + pp.b*dir[1] + pp.c*dir[2];
    if (fabs (den) < EPS) return false;
    t = -(pp.a*pos[0] + pp.b*pos[1] + pp.c*pos[2] + pp.d) / den;

    if (half && t < 0.0) return false;
    X[0] = pos[0] + t*dir[0];
    X[1] = pos[1] + t*dir[1];
    X[2] = pos[2] + t*dir[2];
    return true;
}

// =========================================================================
// check if an intersection point is inside the area of a triangle

bool in_tri (const RVector &X, const TriParam &pp)
{
    double a, b, den;
    double x1, x2, x3, xp, y1, y2, y3, yp;

    if (fabs(pp.a) > fabs(pp.b)) {
	if (fabs(pp.a) > fabs(pp.c)) { // project to yz axis
	    x1 = pp.y1, x2 = pp.y2, x3 = pp.y3, xp = X[1];
	    y1 = pp.z1, y2 = pp.z2, y3 = pp.z3, yp = X[2];
	} else {                       // project to xy axis
	    x1 = pp.x1, x2 = pp.x2, x3 = pp.x3, xp = X[0];
	    y1 = pp.y1, y2 = pp.y2, y3 = pp.y3, yp = X[1];
	}
    } else {
	if (fabs(pp.b) > fabs(pp.c)) { // project to xz axis
	    x1 = pp.x1, x2 = pp.x2, x3 = pp.x3, xp = X[0];
	    y1 = pp.z1, y2 = pp.z2, y3 = pp.z3, yp = X[2];
	} else {                       // project to xy axis
	    x1 = pp.x1, x2 = pp.x2, x3 = pp.x3, xp = X[0];
	    y1 = pp.y1, y2 = pp.y2, y3 = pp.y3, yp = X[1];
	}
    }

    den = (yp-y1)*(x3-x2) + (x1-xp)*(y3-y2);
    if (den == 0.0) return false;
    a = ((y2-y1)*(x3-x2) + (x1-x2)*(y3-y2))/den;
    if (a < 1.0) return false;
    if (den = x3-x2) {
	b = (x1 - x2 + a*(xp-x1))/den;
    } else if (den = y3-y2) {
	b = (y1 - y2 + a*(yp-y1))/den;
    } else if (den = x3-x2) {
	b = (x1 - x2 + a*(xp-x1))/den;
    } else return false;
    if (b < 0.0 || b > 1.0) return false;
    
    return true;
}

// =========================================================================
// check if point p is inside the surface described by prm

bool inside (int ntri, const TriParam *prm, const Point &pos)
{
    int i, nx = 0;
    RVector dir(3), X(3);
    double t;

    // find a search direction - FIX THIS!
    //double phi = rand()*Pi*2.0;
    //double tht = 0.5*Pi-acos(rand()*2.0-1.0);
    //dir[0] = cos(tht)*cos(phi);
    //dir[1] = cos(tht)*sin(phi);
    //dir[2] = sin(tht);
    dir[0] = 1, dir[1] = dir[2] = 0;

    for (i = 0; i < ntri; i++) {
	if (!intersect (pos, dir, prm[i], X, t, true)) continue;
	if (in_tri (X, prm[i])) nx++;
    }
    return nx & 1; // odd number of intersections
}

// =========================================================================
// go through mesh nodes and check for internal ones

void AssignInternalLabels (Mesh *mesh, int ntri, const TriParam *prm,
    int label)
{
    int i;
    int nlen = mesh->nlen();

    for (i = 0; i < nlen; i++) {
	if (inside (ntri, prm, mesh->nlist[i]))
	    mesh->nlist[i].SetRegion (label);
    }
}

void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    int i;

    // mesh handle
    //int hMesh = (int)mxGetScalar (prhs[0]);
    //Mesh *mesh = (Mesh*)hMesh;
    Mesh *mesh = (Mesh*)Handle2Ptr (mxGetScalar (prhs[0]));
    int nlen = mesh->nlen();
    int dim  = mesh->Dimension();

    // surface node list
    RDenseMatrix svtx;
    CopyMatrix (svtx, prhs[1]);
    if (svtx.nCols() != dim) {
	mexPrintf ("toastLabelRegion: invalid surface vertex dimensions\n");
	return;
    }

    // surface triangle list
    IDenseMatrix sidx;
    CopyMatrix (sidx, prhs[2]);
    int ntri = sidx.nRows();
    int nnd  = sidx.nCols();
    if (nnd != dim) {
	mexPrintf ("toastLabelRegion: invalid surface index dimensions\n");
	return;
    }
    // make zero-based
    int *v = sidx.ValPtr();
    for (i = 0; i < ntri*nnd; i++) v[i]--;

    // internal node label
    int ilabel = (int)mxGetScalar (prhs[3]);

    // set up surface parameters
    TriParam *prm = SetupSurface (svtx, sidx);

    // now assign labels
    AssignInternalLabels (mesh, ntri, prm, ilabel);

    // return the new node label list
    mxArray *label = mxCreateDoubleMatrix (nlen, 1, mxREAL);
    double *pr = mxGetPr (label);
    for (i = 0; i < nlen; i++)
	*pr++ = mesh->nlist[i].Region();
    plhs[0] = label;
}
