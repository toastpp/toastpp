// ========================================================================
// Implementation of class MatlabToast
// element-related methods
// ========================================================================

#include "matlabtoast.h"
#include "toastmex.h"

using namespace std;
using namespace toast;

// =========================================================================

void MatlabToast::ElDof (int nlhs, mxArray *plhs[], int nrhs,
    const mxArray *prhs[])
{
    Mesh *mesh = GETMESH_SAFE(0);

    int idx = (int)mxGetScalar(prhs[1]) - 1; // convert to zero-based
    if (idx < 0 || idx >= mesh->elen())
        mexErrMsgTxt ("Eldof: element index out of range");

    Element *pel = mesh->elist[idx];
    int nnd = pel->nNode();

    mxArray *dof = mxCreateDoubleMatrix (1, nnd, mxREAL);
    double *pr = mxGetPr (dof);
    for (int i = 0; i < nnd; i++)
        *pr++ = pel->Node[i]+1;

    plhs[0] = dof;
}

// =========================================================================

void MatlabToast::ElSize (int nlhs, mxArray *plhs[], int nrhs,
    const mxArray *prhs[])
{
    Mesh *mesh = GETMESH_SAFE(0);

    if (nrhs > 1) {
        int idx = (int)mxGetScalar(prhs[1]) - 1; // convert to zero-based
	if (idx < 0 || idx >= mesh->elen())
	    mexErrMsgTxt ("ElementSize:: element index out of range");
	plhs[0] = mxCreateDoubleScalar (mesh->ElSize(idx));
    } else {
        int n = mesh->elen();
	plhs[0] = mxCreateDoubleMatrix (n, 1, mxREAL);
	double *pr = mxGetPr (plhs[0]);
	for (int i = 0; i < n; i++)
	    pr[i] = mesh->ElSize(i);
    }
}

// =========================================================================

void MatlabToast::ElData (int nlhs, mxArray *plhs[], int nrhs,
    const mxArray *prhs[])
{
    int i, j;
    double *pr;

    Mesh *mesh = GETMESH_SAFE(0);

    // element index
    int el = (int)mxGetScalar (prhs[1]) - 1;
    if (el < 0 || el >= mesh->elen()) { // invalid index
	plhs[0] = mxCreateDoubleMatrix (0, 0, mxREAL);
	return;
    }

    Element *pel = mesh->elist[el];
    int dim = pel->Dimension();
    int nnd = pel->nNode();
    int nsd = pel->nSide();
    int nsn = pel->nSideNode(0);
    for (i = 1; i < nsd; i++)
	nsn = ::max (nsn, pel->nSideNode(i));

    mxArray *vtx = mxCreateDoubleMatrix (nnd, dim, mxREAL);
    pr = mxGetPr (vtx);
    for (i = 0; i < dim; i++)
	for (j = 0; j < nnd; j++)
	    *pr++ = mesh->nlist[pel->Node[j]][i];

    mxArray *idx = mxCreateDoubleMatrix (nsd, nsn, mxREAL);
    pr = mxGetPr (idx);
    for (i = 0; i < nsn; i++)
	for (j = 0; j < nsd; j++)
	    if (i < pel->nSideNode(j))
		*pr++ = pel->Node[pel->SideNode(j,i)]+1;
	    else
		*pr++ = 0;
    
    plhs[0] = vtx;
    plhs[1] = idx;
    if (nlhs > 2)
        plhs[2] = mxCreateDoubleScalar(pel->Type());
}

// =========================================================================

void MatlabToast::ElMat (int nlhs, mxArray *plhs[], int nrhs,
    const mxArray *prhs[])
{
    char cbuf[256];
    mxArray *elmat;
    double *pr;

    Mesh *mesh = GETMESH_SAFE(0);

    int idx = (int)mxGetScalar(prhs[1]) - 1; // convert to zero-based
    Element *pel = mesh->elist[idx];
    int i, j, k, l, nnd = pel->nNode();
    int dim = mesh->Dimension();
    int nprm = 3;

    mxGetString (prhs[2], cbuf, 256);

    if (!strcmp(cbuf, "F")) {
	elmat = mxCreateDoubleMatrix (nnd, 1, mxREAL);
	pr = mxGetPr(elmat);
	for (i = 0; i < nnd; i++)
	    pr[i] = pel->IntF(i);
    } else if (!strcmp(cbuf, "FF")) {
	elmat = mxCreateDoubleMatrix (nnd, nnd, mxREAL);
	pr = mxGetPr(elmat);
	for (i = 0; i < nnd; i++) {
	    pr[i*nnd+i] = pel->IntFF(i,i);
	    for (j = 0; j < i; j++)
		pr[i*nnd+j] = pr[j*nnd+i] = pel->IntFF(i,j);
	}
    } else if (!strcmp(cbuf, "FFF")) {
	mwSize dims[3] = {nnd,nnd,nnd};
	elmat = mxCreateNumericArray (3, dims, mxDOUBLE_CLASS, mxREAL);
	pr = mxGetPr(elmat);
	for (i = 0; i < nnd; i++) {
	    for (j = 0; j < nnd; j++) {
		for (k = 0; k < nnd; k++) {
		    pr[(i*nnd+j)*nnd+k] = pel->IntFFF(i,j,k);
		}
	    }
	}
    } else if (!strcmp(cbuf, "DD")) {
	elmat = mxCreateDoubleMatrix (nnd, nnd, mxREAL);
	pr = mxGetPr(elmat);
	for (i = 0; i < nnd; i++) {
	    pr[i*nnd+i] = pel->IntDD(i,i);
	    for (j = 0; j < i; j++)
		pr[i*nnd+j] = pr[j*nnd+i] = pel->IntDD(i,j);
	}
    } else if (!strcmp(cbuf, "FD")) {
	mwSize dims[3] = {nnd, nnd, dim};
	elmat = mxCreateNumericArray (3, dims, mxDOUBLE_CLASS, mxREAL);
	pr = mxGetPr(elmat);
	for (i = 0; i < nnd; i++)
	    for (j = 0; j < nnd; j++) {
		RVector fd = pel->IntFD(i,j);
		if (fd.Dim() == 0) {
		    plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
		    return;
		}
		for (k = 0; k < dim; k++)
		    pr[i + nnd*(j + k*nnd)] = fd[k];
	    }
    } else if (!strcmp(cbuf, "FDD")) {
	mwSize dims[3] = {nnd, nnd, nnd};
	elmat = mxCreateNumericArray (3, dims, mxDOUBLE_CLASS, mxREAL);
	pr = mxGetPr(elmat);
	for (i = 0; i < nnd; i++)
	    for (j = 0; j < nnd; j++)
		for (k = 0; k < nnd; k++)
		    pr[i+nnd*(j+nnd*k)] = pel->IntFDD(i,j,k);
    } else if (!strcmp(cbuf, "dd")) {
	mwSize dims[4] = {nnd, dim, nnd, dim};
	elmat = mxCreateNumericArray (4, dims, mxDOUBLE_CLASS, mxREAL);
	RSymMatrix intdd = pel->Intdd();
	pr = mxGetPr(elmat);
	for (l = 0; l < dim; l++)
	    for (k = 0; k < nnd; k++)
		for (j = 0; j < dim; j++)
		    for (i = 0; i < nnd; i++)
			*pr++ = intdd(i*dim+j,k*dim+l);

    } else if (!strcmp(cbuf, "PFF")) {
	elmat = mxCreateDoubleMatrix (nnd, nnd, mxREAL);
	RVector prm;
	CopyVector (prm, prhs[3]);
	RSymMatrix intPFF = pel->IntPFF(prm);
	pr = mxGetPr(elmat);
	for (i = 0; i < nnd; i++) {
	    pr[i*nnd+i] = intPFF(i,i);
	    for (j = 0; j < i; j++)
	        pr[i*nnd+j] = pr[j*nnd+i] = intPFF(i,j);
	}
	nprm++;
    } else if (!strcmp(cbuf, "PDD")) {
        elmat = mxCreateDoubleMatrix (nnd, nnd, mxREAL);
	RVector prm;
	CopyVector (prm, prhs[3]);
	RSymMatrix intPDD = pel->IntPDD(prm);
	pr = mxGetPr(elmat);
	for (i = 0; i < nnd; i++) {
	    pr[i*nnd+i] = intPDD(i,i);
	    for (j = 0; j < i; j++)
	        pr[i*nnd+j] = pr[j*nnd+i] = intPDD(i,j);
	}
	nprm++;
    } else if (!strcmp(cbuf, "BndF")) {
	// for now, only integrals over all boundary sides are supported
	elmat = mxCreateDoubleMatrix (nnd, 1, mxREAL);
	pr = mxGetPr(elmat);
	RVector bndintf = pel->BndIntF();
	for (i = 0; i < pel->nNode(); i++)
	    pr[i] = bndintf[i];

    } else if (!strcmp(cbuf, "BndFF")) {
	int ii, jj, sd;
	if (nrhs > nprm) sd = (int)mxGetScalar(prhs[nprm+1]) - 1; // side index
	else sd = -1;
	elmat = mxCreateDoubleMatrix (nnd, nnd, mxREAL);
	pr = mxGetPr(elmat);

	if (sd >= 0) { // integral over a single side
	    for (ii = 0; ii < pel->nSideNode(sd); ii++) {
		i = pel->SideNode(sd,ii);
		pr[i*nnd+i] = pel->BndIntFFSide(i,i,sd);
		for (jj = 0; jj < ii; jj++) {
		    j = pel->SideNode(sd,jj);
		    pr[i*nnd+j] = pr[j*nnd+i] = pel->BndIntFFSide(i,j,sd);
		}
	    }
	} else { // integral over all boundary sides
	    for (sd = 0; sd < pel->nSide(); sd++) {
		if (!pel->IsBoundarySide (sd)) continue;
		for (ii = 0; ii < pel->nSideNode(sd); ii++) {
		    i = pel->SideNode(sd,ii);
		    pr[i*nnd+i] += pel->BndIntFFSide(i,i,sd);
		    for (jj = 0; jj < ii; jj++) {
			j = pel->SideNode(sd,jj);
			pr[i*nnd+j] = pr[j*nnd+i] += pel->BndIntFFSide(i,j,sd);
		    }
		}
	    }
	}
    } else if (!strcmp(cbuf, "BndPFF")) {
        int sd;
        if (nrhs > nprm) sd =  (int)mxGetScalar(prhs[nprm+1]) - 1; // side index
	else sd = -1;
	elmat = mxCreateDoubleMatrix (nnd, nnd, mxREAL);
	pr = mxGetPr(elmat);
	RVector prm;
	CopyVector (prm, prhs[3]);

	if (sd >= 0) { // integral over a single side
	    xERROR("Not implemented yet!");
	} else { // integral over all boundary sides
	    for (i = 0; i < nnd; i++) {
	        for (j = 0; j < nnd; j++)
		    pr[i*nnd+j] = pel->BndIntPFF(i,j,prm);
	    }
	}
    } else {
	mexErrMsgTxt ("Elmat: Integral type string not recognised");
    }
    plhs[0] = elmat;
}

// =========================================================================

void MatlabToast::ShapeFunc (int nlhs, mxArray *plhs[], int nrhs,
    const mxArray *prhs[])
{
    int i, j;

    Mesh *mesh = GETMESH_SAFE(0);
    int nlen = mesh->nlen();
    int elen = mesh->elen();
    int dim  = mesh->Dimension();
    bool isglobal = false;

    // element index
    int idx = (int)(floor(mxGetScalar (prhs[1])-0.5));
    ASSERTARG(idx >= 0 && idx < elen, 2, "Invalid element index");
    Element *pel = mesh->elist[idx];
    int nn = pel->nNode();

    double *ppt = mxGetPr (prhs[2]);

    // point given in global coordinates?
    if (nrhs > 3 && mxIsChar(prhs[3])) {
          char cbuf[256];
	  mxGetString (prhs[3], cbuf, 256);
	  isglobal = !strcasecmp (cbuf, "global");
    }

    // evaluation points
    const mwSize *gdim = mxGetDimensions (prhs[2]);
    int npoint = gdim[1];
    ASSERTARG(gdim[0] == dim, 1, "Invalid point dimensions");

    mxArray *sf = mxCreateDoubleMatrix (nn, npoint, mxREAL);
    double *pfunc = mxGetPr (sf);

    Point pt(dim);
    for (j = 0; j < npoint; j++) {      
        for (i = 0; i < dim; i++)
	    pt[i] = ppt[i+j*dim];

	// calculate shape functions
	RVector fun = (isglobal ?
		       pel->GlobalShapeF (mesh->nlist, pt) :
		       pel->LocalShapeF (pt));

	// store shape functions
	for (i = 0; i < nn; i++)
	    pfunc[i+j*nn] = fun[i];
    }
    plhs[0] = sf;    
}

// =========================================================================

void MatlabToast::ShapeGrad (int nlhs, mxArray *plhs[], int nrhs,
    const mxArray *prhs[])
{
  int i, j, k;

    Mesh *mesh = GETMESH_SAFE(0);
    int nlen = mesh->nlen();
    int elen = mesh->elen();
    int dim  = mesh->Dimension();
    bool isglobal = false;

    // element index
    int idx = (int)(floor(mxGetScalar (prhs[1])-0.5));
    ASSERTARG(idx >= 0 && idx < elen, 2, "Invalid element index");
    Element *pel = mesh->elist[idx];
    int nn = pel->nNode();

    double *ppt = mxGetPr (prhs[2]);

    // point given in global coordinates?
    if (nrhs > 3 && mxIsChar(prhs[3])) {
          char cbuf[256];
	  mxGetString (prhs[3], cbuf, 256);
	  isglobal = !strcasecmp (cbuf, "global");
    }

    // evaluation points
    const mwSize *gdim = mxGetDimensions (prhs[2]);
    int npoint = gdim[1];
    ASSERTARG(gdim[0] == dim, 1, "Invalid point dimensions");

    mwSize dims[3] = {dim, nn, npoint};
    mxArray *sg = mxCreateNumericArray (3, dims, mxDOUBLE_CLASS, mxREAL);
    double *pgrad = mxGetPr (sg);

    Point pt(dim);
    for (j = 0; j < npoint; j++) {
        for (i = 0; i < dim; i++)
	    pt[i] = ppt[i+j*dim];

	// calculate shape function derivatives
	RDenseMatrix fgrad = (isglobal ?
			      pel->GlobalShapeD (mesh->nlist, pt) :
			      pel->LocalShapeD (pt));
	double *pg = fgrad.ValPtr();

	// store shape function derivatives
	for (i = 0; i < nn; i++)
	    for (k = 0; k < dim; k++)
	        pgrad[k + (i + j*nn)*dim] = pg[i + k*nn];
    }
    plhs[0] = sg;
}

