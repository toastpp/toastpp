// toastPointCloud
//
// Fill a volume bounded by a surface mesh with a point cloud
// This can subsequently be used to do a Delaunay triangulation
// for generating a volume mesh.
//
// Input meshes must be closed surfaces consisting of 3-noded triangles
 
#include "mex.h"
#include "felib.h"
#include "toastmex.h"

using namespace std;

// Element parameter structure
struct TriParam {
	float x1, y1, z1; // vtx 1
	float x2, y2, z2; // vtx 2
	float x3, y3, z3; // vtx 3
	float a, b, c, d; // plane params
	float d1, d2, d3; // dist of vtx i from opposite edge
};

// ========================================================================
// Prototypes

bool inside (int ntri, const TriParam *pp, const Point3D &pos);

bool intersect (const Point3D &pos, const Point3D &dir, const TriParam &pp,
		Point3D &X, double t, bool half);

bool in_tri (const Point3D &X, const TriParam &pp);

TriParam *SetupElParams (const RDenseMatrix &vtx, const IDenseMatrix &idx);

// ========================================================================
// Check for point inside surface shell

bool inside (int ntri, const TriParam *pp, const Point3D &pos)
{
    int i, nx = 0;
    Point3D X;
    double t;
    
    // find a search direction
    double phi = drand48()*Pi*2.0;
    double tht = 0.5*Pi-acos(drand48()*2.0-1.0);
    Point3D dir (cos(tht)*cos(phi), cos(tht)*sin(phi), sin(tht));

    for (i = 0; i < ntri; i++) {
	if (!intersect (pos, dir, pp[i], X, t, true)) continue;
	if (in_tri (X, pp[i])) nx++;
    }
    return nx & 1; // odd number of intersections
}

// ========================================================================

bool intersect (const Point3D &pos, const Point3D &dir, const TriParam &pp,
		Point3D &X, double t, bool half)
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

// ========================================================================

bool in_tri (const Point3D &X, const TriParam &pp)
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

// ========================================================================

TriParam *SetupElParams (const RDenseMatrix &vtx, const IDenseMatrix &idx)
{
    int i, n;

    int nvtx = vtx.nRows();
    int ntri = idx.nRows();

    TriParam *pp = new TriParam[ntri];
    for (i = 0; i < ntri; i++) {
	n = idx(i,0);
	pp[i].x1 = vtx(n,0); pp[i].y1 = vtx(n,1); pp[i].z1 = vtx(n,2);
	n = idx(i,1);
	pp[i].x2 = vtx(n,0); pp[i].y2 = vtx(n,1); pp[i].z2 = vtx(n,2);
	n = idx(i,2);
	pp[i].x3 = vtx(n,0); pp[i].y3 = vtx(n,1); pp[i].z3 = vtx(n,2);

	pp[i].a = pp[i].y1*(pp[i].z2-pp[i].z3) - pp[i].z1*(pp[i].y2-pp[i].y3) + (pp[i].y2*pp[i].z3 - pp[i].y3*pp[i].z2);
	pp[i].b = -(pp[i].x1*(pp[i].z2-pp[i].z3) - pp[i].z1*(pp[i].x2-pp[i].x3) + (pp[i].x2*pp[i].z3 - pp[i].x3*pp[i].z2));
	pp[i].c = pp[i].x1*(pp[i].y2-pp[i].y3) - pp[i].y1*(pp[i].x2-pp[i].x3) + (pp[i].x2*pp[i].y3 - pp[i].x3*pp[i].y2);
	pp[i].d = -(pp[i].x1*(pp[i].y2*pp[i].z3-pp[i].y3*pp[i].z2) - pp[i].y1*(pp[i].x2*pp[i].z3-pp[i].x3*pp[i].z2) + pp[i].z1*(pp[i].x2*pp[i].y3-pp[i].x3*pp[i].y2));

	RVector dr(3), a(3);
	dr[0] = pp[i].x1-pp[i].x2;
	dr[1] = pp[i].y1-pp[i].y2;
	dr[2] = pp[i].z1-pp[i].z2;
	a[0]  = pp[i].x3-pp[i].x2;
	a[1]  = pp[i].y3-pp[i].y2;
	a[2]  = pp[i].z3-pp[i].z2;
	pp[i].d1 = dot (dr, a)/l2normsq(a);
	dr[0] = pp[i].x2-pp[i].x3;
	dr[1] = pp[i].y2-pp[i].y3;
	dr[2] = pp[i].z2-pp[i].z3;
	a[0]  = pp[i].x1-pp[i].x3;
	a[1]  = pp[i].y1-pp[i].y3;
	a[2]  = pp[i].z1-pp[i].z3;
	pp[i].d2 = dot (dr, a)/l2normsq(a);
	dr[0] = pp[i].x3-pp[i].x1;
	dr[1] = pp[i].y3-pp[i].y1;
	dr[2] = pp[i].z3-pp[i].z1;
	a[0]  = pp[i].x2-pp[i].x1;
	a[1]  = pp[i].y2-pp[i].y1;
	a[2]  = pp[i].z2-pp[i].z1;
	pp[i].d3 = dot (dr, a)/l2normsq(a);
    }
    return pp;
}

// ========================================================================

RDenseMatrix *CreateCloud (const RDenseMatrix &vtx, const IDenseMatrix &idx,
    int npt)
{
    int i, j, n, nin;
    int nvtx = vtx.nRows();
    int ntri = idx.nRows();

    // set up element parameters
    TriParam *pp = SetupElParams (vtx, idx);
    
    // mesh bounding box
    Point3D bbmin (vtx(0,0), vtx(0,1), vtx(0,2));
    Point3D bbmax (vtx(0,0), vtx(0,1), vtx(0,2));
    for (i = 1; i < nvtx; i++) {
	for (j = 0; j < 3; j++) {
	    if      (vtx(i,j) < bbmin[j]) bbmin[j] = vtx(i,j);
	    else if (vtx(i,j) > bbmax[j]) bbmax[j] = vtx(i,j);
	}
    }
    cout << "Mesh bounding box:\n";
    cout << bbmin << " x " << bbmax << endl;

    if (!npt) {
	// estimate target number of points
	double bbvol = 1;
	for (j = 0; j < 3; j++) bbvol *= (bbmax[j]-bbmin[j]);
	double bbrad = pow(bbvol, 1.0/3.0);
	double A = 4.0*Pi*bbrad*bbrad; // approx surface area
	double Ap = A/nvtx;            // area per vertex
	double dst = sqrt(Ap/Pi);      // vertex distance
	double V = bbrad*bbrad*bbrad;  // approx volume
	double Vp = dst*dst*dst;       // volume ver volume vertex
	npt = (int)(V/Vp);
    }
    cout << "Generating " << npt << " internal vertices\n";

    RDenseMatrix *cloud =  new RDenseMatrix(npt,3);

    // do a Monte-Carlo integration
    Point3D pos;
    for (n = nin = 0; nin < npt; n++) {
	for (j = 0; j < 3; j++)
	    pos[j] = drand48()*(bbmax[j]-bbmin[j]) + bbmin[j];
	if (inside (ntri, pp, pos)) {
	    for (j = 0; j < 3; j++)
		(*cloud)(nin,j) = pos[j];
	    nin++;
	}
    }
    double volfrac = (double)nin/(double)n;
    double vol = volfrac;
    for (j = 0; j < 3; j++) vol *= (bbmax[j]-bbmin[j]);
    cout << "Mesh volume: " << vol << endl;
    delete pp;

    return cloud;
}

// ========================================================================

void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // vertex list
    RDenseMatrix vtx;
    CopyMatrix (vtx, prhs[0]);
    if (vtx.nCols() != 3)
	mexErrMsgTxt ("Invalid vertex dimension (should be nvtx x 3)");

    // index list
    RDenseMatrix ridx;
    CopyMatrix (ridx, prhs[1]);
    if (ridx.nCols() != 3)
	mexErrMsgTxt ("Invalid index dimension (should be ntri x 3)");

    // target number of points
    int npt = 0;
    if (nrhs > 2)
	npt = (int)(mxGetScalar (prhs[2])+0.5);

    // copy to integer matrix
    IDenseMatrix idx(ridx.nRows(), ridx.nCols());
    int i, n = ridx.nRows() * ridx.nCols();
    int minidx, maxidx;
    double *rptr = ridx.ValPtr();
    int    *iptr = idx.ValPtr();
    for (i = 0; i < n; i++) {
	iptr[i] = (int)(rptr[i]-0.5);
	if (!i || iptr[i] < minidx) minidx = iptr[i];
	if (!i || iptr[i] > maxidx) maxidx = iptr[i];
    }

    // sanity checks
    if (minidx < 0 || maxidx >= vtx.nRows())
	mexErrMsgTxt ("Element index out of range");

    RDenseMatrix *cloud = CreateCloud (vtx, idx, npt);
    if (nlhs > 0)
	CopyMatrix (&plhs[0], *cloud);
    delete cloud;
}
