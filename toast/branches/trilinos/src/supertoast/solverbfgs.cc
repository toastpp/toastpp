// ==========================================================================
// SolverBFGS: BFGS method (Broyden-Fletcher-Goldfarb-Shanno)
// ==========================================================================

#include "stoastlib.h"
#include "solverbfgs.h"
#include "fwdsolver.h"
#include "jacobian.h"
#include "util.h"

using namespace std;

// ==========================================================================
// external references

//extern clock_t clock0;
extern int g_imgfmt;

// ==========================================================================
// local prototypes

static bool LineSearch (CFwdSolver &FWS, const Raster &raster,
    const Scaler *pscaler, const ObjectiveFunction &OF,
    const CCompRowMatrix &qvec, const CCompRowMatrix &mvec,
    const RVector &data, const RVector &sd, double omega,
    const RVector &p, const RVector &grad, double f0, double &fmin,
    double &lambda, Solution &meshsol, RVector &proj, bool &proj_valid);

// ==========================================================================

SolverBFGS::SolverBFGS (ParamParser *_pp): Solver (_pp)
{
    itmax = 50;
}

// ==========================================================================

void SolverBFGS::Solve (CFwdSolver &FWS, const Raster &raster,
    const Scaler *pscaler, const ObjectiveFunction &OF, const RVector &data,
    const RVector &sd, Solution &bsol, Solution &msol,
    const CCompRowMatrix &qvec, const CCompRowMatrix &mvec, double omega)
{
    int i, j, k;
    const QMMesh *mesh = FWS.MeshPtr();
    int ndat = data.Dim();
    int nlen = mesh->nlen();
    int slen = raster.SLen();
    int blen = raster.BLen();
    int glen = raster.GLen();
    int n    = slen*2;
    bool pvalid;
    double of, fmin, alpha, sum, vt;
    RVector proj(ndat);
    RVector g(n);
    RVector p(n);
    RVector s(n), y(n);

    CVector *dphi = new CVector[mesh->nQ];
    for (i = 0; i < mesh->nQ; i++) dphi[i].New (nlen);

    int i_count = 0;

    // find initial guess for minimizer
    RVector x0 = pscaler->Scale (bsol.GetActiveParams());
    RVector x(x0);

    // find initial objective function
    FWS.Reset (msol, omega);
    FWS.CalcFields (qvec, dphi);
    proj = FWS.ProjectAll_real (mvec, dphi);
    of = OF.get_value (data, proj, sd);

    cout << "Iteration 0: OF: " << of << endl;

    cout << "  Allocating inverse Hessian:  " << n << "x" << n << " ("
	 << (4*n*(n+1))/1048576 << " Mbytes)" << endl;
    RSymMatrix HI(n,n);

#if BFGS_HESSIAN == BFGS_HESSIAN_IDENTITY
    // initial guess for Hessian: identity
    for (i = 0; i < n; i++) HI(i,i) = 1.0;
    cout << "  Initial Hessian: Identity" << endl;
#elif BFGS_HESSIAN == BFGS_HESSIAN_DIAG
    CVector *aphi = new CVector[mesh.nM];
    for (i = 0; i < mesh.nM; i++) aphi[i].New (nlen);
    FWS.CalcFields (mesh, mesh.nM, mvec, aphi);
    cout << "  Allocating Jacobian: " << ndat << "x" << n << " ("
	 << (8*ndat*n)/1048576 << " Mbytes)" << endl;
    RDenseMatrix J (ndat, n);
    RVector M(n);
    RVector kap0(n/2);
    for(i = 0; i < n/2; i++) kap0[i] = x0[i+n/2];
    GenerateJacobian (raster, dphi, aphi, FWS.datatype, J, sd, kap0);
    pscaler->ScaleJacobian (bsol.GetActiveParams(), J);
    ATA_diag (J, M);
    J.New(0,0);
    for (i = 0; i < n; i++) HI(i,i) = 1.0/M[i];
    cout << "  Initial Hessian: diag(JTJ)" << endl;
#else
    cerr << "  ** Invalid Hessian initialisation. Aborting." << endl;
    exit (1);
#endif

    // initial gradient
    OF.add_gradient_data (g, raster, FWS, proj, dphi, mvec);

    // initial guess of step length
    alpha = of / l2norm (g);

    // begin Quasi-Newton iterations
    while (i_count < itmax) {

        HI.Ax (-g, p);

	// line search
	LineSearch (FWS, raster, pscaler, OF, qvec, mvec, data, sd, omega,
		    x, p, of, fmin, alpha, msol, proj, pvalid);

	// update approximate solution
	x += p*alpha;
	bsol.SetActiveParams (pscaler->Unscale(x));
	raster.Map_SolToMesh (bsol, msol, true);

	switch (g_imgfmt) {
	case IMGFMT_NIM:
	    msol.WriteImg_mua (i_count+1, "recon_mua.nim");
	    msol.WriteImg_mus (i_count+1, "recon_mus.nim");
	    break;
	case IMGFMT_RAW:
	    Solution rsol (OT_NPARAM, blen);
	    raster.Map_SolToBasis (bsol, rsol, true);
	    rsol.WriteImg_mua (i_count+1, "recon_mua.raw");
	    rsol.WriteImg_mus (i_count+1, "recon_mus.raw");
	    break;
	}

	// new gradient
	FWS.Reset (msol, omega);
	FWS.CalcFields (qvec, dphi);
	proj = FWS.ProjectAll_real (mvec, dphi);
	RVector g1(n);
	OF.add_gradient_data (g1, raster, FWS, proj, dphi, mvec);
	of = OF.get_value (data, proj, sd);

	RVector x1(bsol.GetActiveParams());
	s = x1-x0;
	y = g1-g;

	// update inverse Hessian
	// Byrd et. al. 1996 (p.2)

	cout << "  Updating inverse Hessian" << endl;
	ofstream ofs("HI.dat"); // temporary storage of new inverse Hessian
	RVector Vj(n), HVj(n), VTHVj(n);
	double rho = 1.0 / (y & s);
	for (j = 0; j < n; j++) {
  	    cout << j << endl;
  	    // build j-th column of V = y s^T
	    for (i = 0; i < n; i++)
	        Vj[i] = -rho*y[i]*s[j];
	    Vj[j] += 1.0;
	    // build j-th column of HV
	    for (i = 0; i < n; i++) {
	        for (sum = 0.0, k = 0; k < n; k++)
		    sum += HI(i,k)*Vj[k];
		HVj[i] = sum;
	    }
	    // build j-th column of V^T HV
	    for (i = 0; i < n; i++) {
	        for (sum = 0.0, k = 0; k < n; k++) {
		    vt = -y[k]*s[i];
		    if (i == k) vt += 1.0;
		    sum += vt*HVj[k];
		}
		VTHVj[i] = sum;
	    }
	    // add rho * s s^T
	    for (i = 0; i < n; i++)
	        VTHVj[i] += rho * s[i]*s[j];
	    // write column j to disk
	    ofs << VTHVj << endl;
	}
	ofs.close();
	// read new inverse Hessian back into H
	ifstream ifs("HI.dat");
	for (j = 0; j < n; j++) {
	    ifs >> VTHVj;
	    for (i = j; i < n; i++) // only need lower triangle
	        HI(i,j) = VTHVj[i];
	}
	ifs.close();
	//RVector h1(n), h2(n);
	//H.Ax(s, h1);
	//double den1 = 1.0/(h1 & s);
	//double den2 = 1.0/(y & s);
	//for (i = 0; i < n; i++) {
	//    double h = 0.0;
	//    for (j = 0; j < n; j++) h += s[j] * H(j,i);
	//}
	//for (i = 0; i < n; i++)
  	//    for (j = 0; j <= i; j++) // symmetric, so only need lower triangle
	//        H(i,j) += -h1[i]*h2[j]*den1 + y[i]*y[j]*den2;

	i_count++;
	x0 = x1;
	g  = g1;
	cout << "Iteration: " << i_count << "  OF: " << of << endl;
    }
}

// ==========================================================================

void SolverBFGS::ReadParams (ParamParser &pp)
{
    // === MAX ITERATION COUNT ===

    if (!(pp.GetInt ("BFGS_ITMAX", itmax) || pp.GetInt ("NONLIN_ITMAX", itmax))
	  || itmax <= 0) {
	do {
	    cout << "\nMax number of BFGS iterations (>0):\n";
	    cout << ">> ";
	    cin >> itmax;
	} while (itmax <= 0);
    }
}

// ==========================================================================

void SolverBFGS::WriteParams (ParamParser &pp)
{
    pp.PutString ("SOLVER", "BFGS");
    pp.PutInt ("BFGS_ITMAX", itmax);
}

// ==========================================================================

bool LineSearch (CFwdSolver &FWS, const Raster &raster, const Scaler *pscaler,
    const ObjectiveFunction &OF, const CCompRowMatrix &qvec,
    const CCompRowMatrix &mvec, const RVector &data, const RVector &sd,
    double omega, const RVector &p, const RVector &grad, double f0,
    double &fmin, double &lambda, Solution &meshsol, RVector &proj,
    bool &proj_valid)
{
    LOGOUT ("Linesearch Start");
    const int MAXIT = 16;
    double x0 = 0.0, xm, fm, fb, xb = lambda;
    proj_valid = false;
    
    RVector p1 = p+grad*xb;
    raster.Map_ActiveSolToMesh (pscaler->Unscale(p1), meshsol);
    while (!meshsol.Valid()) {
	LOGOUT ("** Invalid nodal parameters. Reducing step size");
	xb *= 0.5;
	p1 = p+grad*xb;
	raster.Map_ActiveSolToMesh (pscaler->Unscale(p1), meshsol);
    }
    proj = FWS.ProjectAll_real (qvec, mvec, meshsol, omega);
    fb = ObjectiveFunction::get_value (data, proj, sd);
    LOGOUT("Step  %f  OF %f", xb, fb);

    if (fb < f0) { // increase interval
        xm = xb; fm = fb;
	xb *= 2.0;
	p1 = p+grad*xb;
	raster.Map_ActiveSolToMesh (pscaler->Unscale(p1), meshsol);
	if (!meshsol.Valid()) {
	    LOGOUT ("** Invalid nodal parameters. Truncating step size.");
	    fb = fm*10.0; // invalidate this step
	} else {
	    proj = FWS.ProjectAll_real (qvec, mvec, meshsol, omega);
	    fb = ObjectiveFunction::get_value (data, proj, sd);
	    LOGOUT("Step  %f  OF %f", xb, fb);
	}
	while (fb < fm) {
	    x0 = xm; f0 = fm;
	    xm = xb; fm = fb;
	    xb *= 2.0;
	    p1 = p+grad*xb;
	    raster.Map_ActiveSolToMesh (pscaler->Unscale(p1), meshsol);
	    proj = FWS.ProjectAll_real (qvec, mvec, meshsol, omega);
	    fb = ObjectiveFunction::get_value (data, proj, sd);
	    LOGOUT("Step  %f  OF %f", xb, fb);
	}
    } else { // decrease interval
        xm = 0.5*xb;
	p1 = p+grad*xm;
	raster.Map_ActiveSolToMesh (pscaler->Unscale(p1), meshsol);
	proj = FWS.ProjectAll_real (qvec, mvec, meshsol, omega);
	fm = ObjectiveFunction::get_value (data, proj, sd);
	LOGOUT("Step  %f  OF %f", xm, fm);
	int itcount = 0;
	while (fm > f0) {
  	    if (++itcount > MAXIT) return false;
	    xb = xm; fb = fm;
	    xm = 0.5*xb;
	    p1 = p+grad*xm;
	    raster.Map_ActiveSolToMesh (pscaler->Unscale(p1), meshsol);
	    proj = FWS.ProjectAll_real (qvec, mvec, meshsol, omega);
	    fm = ObjectiveFunction::get_value (data, proj, sd);
	    LOGOUT("Step  %f  OF %f", xm, fm);
	}
    }
    // quadratic interpolation
    double a = ((f0-fb)/(x0-xb) - (f0-fm)/(x0-xm)) / (xb-xm);
    double b = (f0-fb)/(x0-xb) - a*(x0+xb);
    lambda = -b/(2.0*a);
    p1 = p+grad*lambda;
    raster.Map_ActiveSolToMesh (pscaler->Unscale(p1), meshsol);
    proj = FWS.ProjectAll_real (qvec, mvec, meshsol, omega);
    fmin = ObjectiveFunction::get_value (data, proj, sd);
    if (fmin > fm) {  // interpolation didn't give improvement
        lambda = xm, fmin = fm;
    } else proj_valid = true;
    LOGOUT("Final %f  OF %f", lambda, fmin);
    LOGOUT ("Linesearch End");
    return true;
}
