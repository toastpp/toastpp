// ==========================================================================
// SolverLIN: Linear solver for difference imaging
// ==========================================================================

#include "stoastlib.h"
#include "pscaler.h"
#include "solverlin.h"
#include "fwdsolver.h"
#include "jacobian.h"
#include "util.h"

using namespace std;

// ==========================================================================
// external references

extern int bWriteJ;
extern char g_meshname[256];

// ==========================================================================

static double lambda = 1.0;

SolverLIN::SolverLIN (ParamParser *_pp): Solver (_pp)
{
}

void SolverLIN::Solve (CFwdSolver &FWS, const Raster &raster,
    const Scaler *pscaler, const ObjectiveFunction &OF, const RVector &data,
    const RVector &sd, Solution &bsol, Solution &msol,
    const CCompRowMatrix &qvec, const CCompRowMatrix &mvec, double omega)
{
    cout << "Starting linear solver ..." << endl;

    int i, niter;
    const QMMesh *mesh = FWS.MeshPtr();
    int ndat = data.Dim();
    int nlen = mesh->nlen();
    int slen = raster.SLen();
    int n    = slen*2;

    RVector proj(ndat);
    RVector kap = bsol.GetParam (OT_CKAPPA);

    // allocate field vectors
    CVector *dphi = new CVector[mesh->nQ];
    for (i = 0; i < mesh->nQ; i++) dphi[i].New (nlen);
    CVector *aphi = new CVector[mesh->nM];
    for (i = 0; i < mesh->nM; i++) aphi[i].New (nlen);

    LOGOUT("Resetting forward solver");
    FWS.Reset (msol, omega);
    FWS.CalcFields (qvec, dphi);
    FWS.CalcFields (mvec, aphi);
    proj = FWS.ProjectAll_real (mvec, dphi);

    // subtract baseline from data
    RVector b = (data-proj)/sd;

    // generate Jacobian
    LOGOUT("Calculating Jacobian ...");
    RDenseMatrix *J = new RDenseMatrix (ndat, n);
    GenerateJacobian (&raster, mesh, mvec, dphi, aphi,
        FWS.GetDataScaling(), *J);

    // apply data space rescaling T
    J->RowScale (inv(sd));

    // apply solution space rescaling S
    pscaler->ScaleJacobian (bsol.GetActiveParams(), *J);

    if (bWriteJ) {
	WriteNimHeader (g_meshname, nlen, "pmdf_mua.nim", "MUA");
	WriteNimHeader (g_meshname, nlen, "pmdf_kappa.nim", "MUS");
	for (i = 0; i < J->nRows(); i++) {
	    RVector row = J->Row (i);
	    RVector jmua (row, 0, slen);
	    RVector jkappa (row, slen, slen);
	    RVector mjmua(n);
	    RVector mjkappa (n);
	    raster.Map_SolToMesh (jmua, mjmua);
	    raster.Map_SolToMesh (jkappa, mjkappa);
	    WriteImage (mjmua, i, "pmdf_mua.nim");
	    WriteImage (mjkappa, i, "pmdf_kappa.nim");
	}
    }

    // now apply solver dx = J^T (JJT)^-1 dy
    RSymMatrix JJT = AAT(*J);

    RVector diag = JJT.Diag();
    double mn = mean (diag);
    cout << "Mean of diagonal: " << mn << endl;

    for (i = 0; i < JJT.nRows(); i++)
	JJT(i,i) += lambda*mn;

    RVector JTdy(ndat), dx(n);
    if (FWS.LinSolver() == LSOLVER_DIRECT) {
	LOGOUT ("Entering Cholesky solver ...");
	CHdecomp (JJT, true);
	JTdy = CHsubst (JJT, b);
	LOGOUT ("Finished");
    } else {
	LOGOUT ("Entering PCG solver ...");
	double tol = FWS.GetLinSolverTol();
	RPrecon_Diag precon;
	precon.Reset (&JJT);
	niter = PCG (JJT, b, JTdy, tol, &precon);
	LOGOUT("Converged to tolerance %f after %d iterations", tol, niter);
    }
    dx = ATx (*J, JTdy);
    RVector x = dx; // need to add baseline params here

    bsol.SetActiveParams (bsol.GetActiveParams() + pscaler->Unscale(x));
    raster.Map_SolToMesh (bsol, msol);
    msol.WriteImg_mua (0, "recon_mua.nim");
    msol.WriteImg_mus (0, "recon_mus.nim");
    RVector img(raster.GLen());

    const double c = 0.3/1.4;
    RVector cmua   = bsol.GetParam(OT_CMUA);
    RVector ckappa = bsol.GetParam(OT_CKAPPA);
    
    raster.Map_SolToGrid (cmua, img);
    RVector mua = img/c;
    raster.Map_SolToGrid (ckappa, img);
    RVector kappa = img/c;

    cerr << "mua range: " << vmin (mua) << ' ' << vmax(mua) << endl;
    cerr << "kappa range: " << vmin (kappa) << ' ' << vmax(kappa) << endl;

    WritePPM (mua, raster.GDim(), 0, 0, "recon_mua.ppm");
    WritePPM (kappa, raster.GDim(), 0, 0, "recon_kappa.ppm");
}

void SolverLIN::ReadParams (ParamParser &pp)
{
}

void SolverLIN::WriteParams (ParamParser &pp)
{
    pp.PutString ("SOLVER", "LINEAR");
}
