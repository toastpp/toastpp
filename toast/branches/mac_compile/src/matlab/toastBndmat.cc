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
char *Dimstr[4] = {"X", "Y", "Z", "R"};

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
    double nv;
    if (dim < 3)
	nv = mesh->nlist[nd][dim];
    else
	nv = length(mesh->nlist[nd]);

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

bool Contained (int idx, const IVector &vtxlist)
{
    int i, len = vtxlist.Dim();
    for (i = 0; i < len; i++)
	if (vtxlist[i] == idx) return true;
    return false;
}

void CalcBndmat (Mesh *mesh, char *intstr, int dim, int rel, double v,
    mxArray **res)
{
    // sysmatrix structure
    int el, sd, i, j, is, js, ir, jr, nside, nnode;
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
	    for (sd = 0; sd < nside; sd++) {
		if (!pel->IsBoundarySide(sd)) continue;
		nnode = pel->nSideNode(sd);
		for (i = 0; i < nnode; i++) {
		    ir = pel->SideNode(sd,i);
		    is = pel->Node[ir];
		    if (subreg && !Inside (mesh, is, dim, rel, v)) continue;
		    for (j = 0; j < nnode; j++) {
			jr = pel->SideNode(sd,j);
			js = pel->Node[jr];
			if (subreg && !Inside (mesh, js, dim, rel,v)) continue;
			val = pel->BndIntFFSide (ir, jr, sd);
			F.Add (is, js, val);
		    }
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

void CalcBndmat (Mesh *mesh, char *intstr, const IVector &vtxlist,
    mxArray **res)
{
    // sysmatrix structure
    int el, sd, i, j, is, js, ir, jr, nside, nnode;
    int n = mesh->nlen();
    int nel = mesh->elen();
    int *rowptr, *colidx, nzero;
    double val;
    mesh->SparseRowStructure (rowptr, colidx, nzero);
    RCompRowMatrix F (n, n, rowptr, colidx);
    delete []rowptr;
    delete []colidx;

    if (!strcasecmp (intstr, "FF")) {
	for (el = 0; el < nel; el++) {
	    Element *pel = mesh->elist[el];
	    nside = pel->nSide();
	    for (sd = 0; sd < nside; sd++) {
		if (!pel->IsBoundarySide(sd)) continue;
		nnode = pel->nSideNode(sd);
		for (i = 0; i < nnode; i++) {
		    ir = pel->SideNode(sd,i);
		    is = pel->Node[ir];
		    if (!Contained (is,vtxlist)) continue;
		    for (j = 0; j < nnode; j++) {
			jr = pel->SideNode(sd,j);
			js = pel->Node[jr];
			if (!Contained (js,vtxlist)) continue;
			val = pel->BndIntFFSide (ir, jr, sd);
			F.Add (is, js, val);
		    }
		}
	    }
	}
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
	if (mxIsChar (prhs[2])) {
	    char regionstr[256];
	    char dimstr[32], relstr[32];
	    mxGetString (prhs[2], regionstr, 256);
	    sscanf (regionstr, "%s %s %lf", dimstr, relstr, &val);

	    switch (toupper(dimstr[0])) {
	    case 'X': dim = 0; break;
	    case 'Y': dim = 1; break;
	    case 'Z': dim = 2; break;
	    case 'R': dim = 3; break;
	    }
	    for (i = 0; i <= 4; i++) {
		if (!strcasecmp (relstr, Relstr[i])) {
		    rel = i; break;
		}
	    }
	    CalcBndmat (mesh, intstr, dim, rel, val, &plhs[0]);
	} else {
	    RVector tmp;
	    CopyVector(tmp,prhs[2]);
	    IVector vtxlist (tmp.Dim());
	    for (i = 0; i < tmp.Dim(); i++)
		vtxlist[i] = (int)(tmp[i]-0.5); // switch to 0-based
	    CalcBndmat (mesh, intstr, vtxlist, &plhs[0]);
	}
    } else {
	CalcBndmat (mesh, intstr, -1, rel, val, &plhs[0]);
    }

}
