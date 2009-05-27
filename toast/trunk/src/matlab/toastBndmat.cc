// =========================================================================
// toastSysmat
// Generate a boundary system matrix from a mesh by integrating over the
// boundary
//
// RH parameters:
//     1: mesh handle (double)
// LH parameters:
//     1: system matrix (sparse double matrix)
//     2: (optional) boundary part of system matrix (sparse)
//     3: (optional) boundary matrix prefactors (alpha)
// =========================================================================


#include "mex.h"
#include "stoastlib.h"
#include "toastmex.h"
#include "fwdsolver.h"
#include "util.h"

#define EQ 0
#define GT 1
#define GE 2
#define LT 3
#define LE 4

char *Relstr[5] = {"EQ", "GT", "GE", "LT", "LE"};
char *Dimstr[3] = {"X", "Y", "Z"};

using namespace std;

// =========================================================================
// Implementation

void Assert (bool cond, const char *msg)
{
    if (!cond) {
	char cbuf[256] = "toastSysmat: ";
	strcat (cbuf, msg);
	mexErrMsgTxt (cbuf);
    }
}

bool Inside (Mesh *mesh, int nd, int dim, int rel, double v)
{
    const double EPS = 1e-8;
    double nv = mesh->nlist[nd][dim];
    switch (rel) {
    case EQ:
	return (fabs(nv-v) < EPS);
    case GT:
	return (nv > v);
    case GE:
	return (nv > v-EPS);
    case LT:
	return (nv < v);
    case LE:
	return (nv < v+EPS);
    default:
	return false;
    }
}


void CalcBndmat (Mesh *mesh, char *intstr, int dim, int rel, double v,
    mxArray **res)
{
    // sysmatrix structure
    int el, i, j, is, js, nside, nnode;
    int n = mesh->nlen();
    int nel = mesh->elen();
    int *rowptr, *colidx, nzero;
    double val;
    bool subreg = (dim >= 0);
    mesh->SparseRowStructure (rowptr, colidx, nzero);
    RCompRowMatrix F (n, n, rowptr, colidx);
    delete []rowptr;
    delete []colidx;

    if (!strcasecmp (intstr, "FF")) {
	for (el = 0; el < nel; el++) {
	    Element *pel = mesh->elist[el];
	    nside = pel->nSide();
	    nnode = pel->nNode();
	    for (i = 0; i < nnode; i++) {
		is = pel->Node[i];
		if (subreg && !Inside (mesh, is, dim, rel, v)) continue;
		for (j = 0; j < nnode; j++) {
		    js = pel->Node[j];
		    if (subreg && !Inside (mesh, js, dim, rel, v)) continue;
		    val = pel->BndIntFF (i, j);
		    F.Add (is, js, val);
		}
	    }
	}


#ifdef UNDEF
	RVector dummy(n);
	dummy = 1.0;
	AddToSysMatrix (*mesh, F, &dummy, ASSEMBLE_BNDPFF);
#endif
    }
    CopyMatrix (res, F);
}

// ============================================================================
// MATLAB entry point

void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    QMMesh *mesh = (QMMesh*)Handle2Ptr (mxGetScalar (prhs[0]));

    char intstr[256];
    RVector prm;
    int n, i;
    bool elbasis = false;
    int rel = -1;
    int dim = -1;
    double val;

    mxGetString (prhs[1], intstr, 256);
    
    if (nrhs > 2) {
	char regionstr[256];
	char dimstr[32], relstr[32];
	mxGetString (prhs[2], regionstr, 256);
	sscanf (regionstr, "%s %s %lf", dimstr, relstr, &val);
    cerr << relstr << endl;
	switch (toupper(dimstr[0])) {
	case 'X': dim = 0; break;
	case 'Y': dim = 1; break;
	case 'Z': dim = 2; break;
	}
	for (i = 0; i <= 4; i++) {
	    if (!strcasecmp (relstr, Relstr[i])) {
		rel = i; break;
	    }
	}
    }

    if (dim >= 0)
	cerr << "toastBndmat: region limit: " << Dimstr[dim] << ' '
	     << Relstr[rel] << ' ' << val << endl;

    CalcBndmat (mesh, intstr, dim, rel, val, &plhs[0]);
}
