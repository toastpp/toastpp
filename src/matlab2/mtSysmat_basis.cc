#include "matlabtoast.h"
#include "mexutil.h"
#include "util.h"

static void CalcSysmat (const Raster2 *raster, RVector &mua, RVector &mus,
		 RVector &ref, double freq, mxArray **res)
{
    int i, nz;
    idxtype *rowptr, *colidx;

    const Mesh &mesh = raster->mesh();
    int n = mesh.nlen();
    int nel = mesh.elen();
    int blen = raster->BLen();
    mesh.SparseRowStructure (rowptr, colidx, nz);
    RCompRowMatrix K(n,n,rowptr,colidx);
    delete []rowptr;
    delete []colidx;

    // parameter transformations
    RVector cmua = mua*c0/ref;
    RVector ckap = c0/(3.0*ref*(mua+mus));
    RVector c2a(blen);
    for (int i = 0; i < blen; i++)
	c2a[i] = (ref[i] ? c0/(2.0*ref[i]*A_Keijzer (ref[i])) : 0.0);

    // temporary: map c2a parameter to mesh basis for ASSEMBLE_BNDPFF
    RVector hc2a(n);
    raster->Map_BasisToMesh(c2a, hc2a);

    for (i = 0; i < nel; i++) {
	raster->AddToElMatrix (i, K, &cmua, ASSEMBLE_PFF);
	raster->AddToElMatrix (i, K, &ckap, ASSEMBLE_PDD);
	raster->AddToElMatrix (i, K, &hc2a, ASSEMBLE_BNDPFF);
    }
    CopyMatrix (res, K);
}

void MatlabToast::Sysmat_basis (int nlhs, mxArray *plhs[], int nrhs,
    const mxArray *prhs[])
{
    Raster2 *raster = (Raster2*)GETBASIS_SAFE(0);

    RVector mua, mus, ref;
    double freq;
    int blen;

    CopyVector (mua, prhs[1]);
    CopyVector (mus, prhs[2]);
    CopyVector (ref, prhs[3]);
    freq = mxGetScalar (prhs[4]);

    blen = raster->BLen();
    ASSERTARG (blen == mua.Dim(), 2, "wrong size");
    ASSERTARG (blen == mus.Dim(), 3, "wrong size");
    ASSERTARG (blen == ref.Dim(), 4, "wrong size");

    CalcSysmat (raster, mua, mus, ref, freq, &plhs[0]);
}

