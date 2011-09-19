#include <iostream>
#include "mathlib.h"
#include "felib.h"
#include "stoastlib.h"
#include "timing.h"
#include "toastcuda.h"

#ifdef USE_CUDA_FLOAT
#include <thrust/extrema.h>
#include <cusp/csr_matrix.h>
#include "vector_cusp.h"
#include "crmatrix_cusp.h"
#endif

using namespace std;

static double refind = 1.4;

// =================================================================

#ifdef UNDEF
void Test_Vector ()
{
    const int n = 100;
    int i;
    TVector<float> v1(n);
    for (i = 0; i < n; i++)
        v1[i] = i;
    cuspTVector<float> dv1;
    dv1.Set (v1);
    cuspTVector<float> dv2;
    dv2 = dv1;
    TVector<float> v2;
    v2 = dv2.Get();
    cout << "v2 = " << v2 << endl;
}
#endif

// =================================================================

void Test_Ax_float (QMMesh &mesh)
{
    cout << "CUDA test: Ax (float)" << endl;

    int i, j, n = mesh.nlen();
    int nQ = mesh.nQ;
    double qwidth = 2;
    SourceMode srctp = SRCMODE_NEUMANN;
    double t0, dt;
    double tol = 1e-8;

    // allocate a system matrix
    int *rowptr, *colidx, nzero;
    mesh.SparseRowStructure (rowptr, colidx, nzero);
    RCompRowMatrix F (n, n, rowptr, colidx);
  
    // assemble it
    double mua0 = 0.01;
    double mus0 = 1.0;
    double kap0 = 1.0/(3.0*(mua0+mus0));
    double A = A_Keijzer(refind);

    RVector cmua(n); cmua = c0/refind * mua0;
    RVector ckap(n); ckap = c0/refind * kap0;
    RVector c2a(n); c2a = c0/(2.0*refind*A);
    AddToSysMatrix (mesh, F, &cmua, ASSEMBLE_PFF);
    AddToSysMatrix (mesh, F, &ckap, ASSEMBLE_PDD);
    AddToSysMatrix (mesh, F, &c2a, ASSEMBLE_BNDPFF);

    // convert to single precision
    FCompRowMatrix FF (n, n, rowptr, colidx);
    double *pf = F.ValPtr();
    float *pff = FF.ValPtr();
    for (i = 0; i < nzero; i++)
        *pff++ = (float)*pf++;

    // build source vectors
    FVector *qvec = new FVector[nQ];
    for (i = 0; i < nQ; i++) {
	RVector rq = QVec_Gaussian (mesh, mesh.Q[i], qwidth, srctp);
	qvec[i].New(n);
	for (j = 0; j < n; j++) qvec[i][j] = (float)rq[j];
    }

    // calculate Ax (serial)
    FVector *bs = new FVector[nQ];
    for (i = 0; i < nQ; i++)
	bs[i].New(n);

    t0 = clock();
    for (i = 0; i < nQ; i++)
	FF.Ax (qvec[i], bs[i]);
    dt = (double)(clock()-t0)/(double)CLOCKS_PER_SEC;
    cout << "Time(Ax) [serial]: " << dt << endl;

    // calculate Ax (parallel)
    FVector *bp = new FVector[nQ];
    for (i = 0; i < nQ; i++)
	bp[i].New(n);

    t0 = clock();
    FF.Ax (qvec, bp, nQ);
    dt = (double)(clock()-t0)/(double)CLOCKS_PER_SEC;
    cout << "Time(Ax) [parallel]: " << dt << endl;

    // calculate by hand for validation
    int *rp = FF.rowptr;
    int *ci = FF.colidx;
    float *val = FF.ValPtr();
    FVector q = qvec[0];
    float *x = q.data_buffer();
    FVector b2(n);
    for (i = 0; i < n; i++) {
	double sum = 0.0;
	for (j = rp[i]; j < rp[i+1]; j++)
	    sum += val[j] * x[ci[j]];
	b2[i] = sum;
    }

    cout << "Validation: ||Ax_serial||   = " << l2norm(bs[0]) << endl
	 << "            ||Ax_parallel|| = " << l2norm(bp[0]) << endl
	 << "            ||Ax_target||   = " << l2norm(b2)    << endl;

    delete []qvec;
    delete []bs;
    delete []bp;
}

// =================================================================

void Test_Ax_double (QMMesh &mesh)
{
    cout << "CUDA test: Ax (double)" << endl;

    int i, j, n = mesh.nlen();
    int nQ = mesh.nQ;
    double qwidth = 2;
    SourceMode srctp = SRCMODE_NEUMANN;
    double t0, dt;
    double tol = 1e-8;

    // allocate a system matrix
    int *rowptr, *colidx, nzero;
    mesh.SparseRowStructure (rowptr, colidx, nzero);
    RCompRowMatrix F (n, n, rowptr, colidx);
  
    // assemble it
    double mua0 = 0.01;
    double mus0 = 1.0;
    double kap0 = 1.0/(3.0*(mua0+mus0));
    double A = A_Keijzer(refind);

    RVector cmua(n); cmua = c0/refind * mua0;
    RVector ckap(n); ckap = c0/refind * kap0;
    RVector c2a(n); c2a = c0/(2.0*refind*A);
    AddToSysMatrix (mesh, F, &cmua, ASSEMBLE_PFF);
    AddToSysMatrix (mesh, F, &ckap, ASSEMBLE_PDD);
    AddToSysMatrix (mesh, F, &c2a, ASSEMBLE_BNDPFF);

    // build source vectors
    RVector *qvec = new RVector[nQ];
    for (i = 0; i < nQ; i++) {
	qvec[i] = QVec_Gaussian (mesh, mesh.Q[i], qwidth, srctp);
    }

    // calculate Ax (serial)
    RVector *bs = new RVector[nQ];
    for (i = 0; i < nQ; i++)
	bs[i].New(n);

    t0 = clock();
    for (i = 0; i < nQ; i++)
	F.Ax (qvec[i], bs[i]);
    dt = (double)(clock()-t0)/(double)CLOCKS_PER_SEC;
    cout << "Time(Ax) [serial]: " << dt << endl;

    // calculate Ax (parallel)
    RVector *bp = new RVector[nQ];
    for (i = 0; i < nQ; i++)
	bp[i].New(n);

    t0 = clock();
    F.Ax (qvec, bp, nQ);
    dt = (double)(clock()-t0)/(double)CLOCKS_PER_SEC;
    cout << "Time(Ax) [parallel]: " << dt << endl;

    // calculate by hand for validation
    int *rp = F.rowptr;
    int *ci = F.colidx;
    double *val = F.ValPtr();
    RVector q = qvec[0];
    double *x = q.data_buffer();
    RVector b2(n);
    for (i = 0; i < n; i++) {
	double sum = 0.0;
	for (j = rp[i]; j < rp[i+1]; j++)
	    sum += val[j] * x[ci[j]];
	b2[i] = sum;
    }

    cout << "Validation: ||Ax_serial||   = " << l2norm(bs[0]) << endl
	 << "            ||Ax_parallel|| = " << l2norm(bp[0]) << endl
	 << "            ||Ax_target||   = " << l2norm(b2)    << endl;

    delete []qvec;
    delete []bs;
    delete []bp;
}

// =================================================================

void Test_BiCGSTAB (QMMesh &mesh, int iq)
{
    cout << "CUDA test: BiCGSTAB" << endl;

    double qwidth = 2;
    SourceMode srctp = SRCMODE_NEUMANN;

    int i, n = mesh.nlen();
    double t0, dt;
    double tol = 1e-8;
    int itmax = 200;

    // allocate a system matrix
    int *rowptr, *colidx, nzero;
    mesh.SparseRowStructure (rowptr, colidx, nzero);
    RCompRowMatrix F (n, n, rowptr, colidx);
  
    // assemble it
    double mua0 = 0.01;
    double mus0 = 1.0;
    double kap0 = 1.0/(3.0*(mua0+mus0));
    double A = A_Keijzer(refind);

    RVector cmua(n); cmua = c0/refind * mua0;
    RVector ckap(n); ckap = c0/refind * kap0;
    RVector c2a(n); c2a = c0/(2.0*refind*A);
    AddToSysMatrix (mesh, F, &cmua, ASSEMBLE_PFF);
    AddToSysMatrix (mesh, F, &ckap, ASSEMBLE_PDD);
    AddToSysMatrix (mesh, F, &c2a, ASSEMBLE_BNDPFF);

    // convert to single precision
    FCompRowMatrix FF (n, n, rowptr, colidx);
    double *pf = F.ValPtr();
    float *pff = FF.ValPtr();
    for (i = 0; i < nzero; i++)
        *pff++ = (float)*pf++;

    // make a source vector
    FVector b(n);
    RVector rq = QVec_Gaussian (mesh, mesh.Q[iq], qwidth, srctp);
    for (i = 0; i < n; i++) b[i] = (float)rq[i];

    // Solve Ax=b
    FVector x(n);
    FPreconditioner *fp = 0;
    double tol0 = tol;
    t0 = clock();
    int it = BiCGSTAB (FF, b, x, tol, fp, itmax);
    dt = (double)(clock()-t0)/(double)CLOCKS_PER_SEC;

    // Verify
    FVector b2(n);
    b2 = Ax(FF,x);

    // Results
    cout << "Target:  ||b||  = " << l2norm(b) << endl;
    cout << "Result:  ||Ax|| = " << l2norm(b2) << endl;
    cout << "Iterations:     " << it << " (max: " << itmax << ")" << endl;
    cout << "Error:          " << tol << " (target: " << tol0 << ")" << endl;
    cout << "Time(BiCGSTAB): " << dt << endl;

    delete []rowptr;
    delete []colidx;
}

// =================================================================

void Test_FwdSolver (QMMesh &mesh)
{
    cout << "CUDA test: FwdSolver" << endl;

    int i, j;
    double tol = 1e-8;
    double c = c0/refind;
    double qwidth = 2;
    double mwidth = 2;
    SourceMode srctp = SRCMODE_NEUMANN;
    int n = mesh.nlen();
    int nQ = mesh.nQ;
    int nM = mesh.nM;

    FFwdSolver fws("BICGSTAB", tol);
    fws.SetLinSolverMaxit (200);
    FCompRowMatrix qvec, mvec;
    FVector *dphi;
    double t0, dt;

    double mua0 = 0.01;
    double mus0 = 1.0;
    double kap0 = 1.0/(3.0*(mua0+mus0));
    double A = A_Keijzer (refind);

    Solution sol(OT_NPARAM, n);
    RVector cmua(n); cmua = c*mua0;
    RVector ckap(n); ckap = c*kap0;
    RVector ref(n);  ref  = refind;
    RVector c2a(n);  c2a  = c0/(2.0*refind*A);
    sol.SetParam (OT_CMUA, cmua);
    sol.SetParam (OT_CKAPPA, ckap);
    sol.SetParam (OT_N, ref);
    sol.SetParam (OT_C2A, c2a);

    // build the field vectors
    dphi = new FVector[nQ];
    for (i = 0; i < nQ; i++) dphi[i].New(n);

    // build source vectors
    qvec.New (nQ, n);
    for (i = 0; i < nQ; i++) {
	RVector rq = QVec_Gaussian (mesh, mesh.Q[i], qwidth, srctp);
	FVector fq(n);
	for (j = 0; j < n; j++) fq[j] = (float)rq[j];
	qvec.SetRow (i, fq);
    }

    // build measurement vectors
    mvec.New (nM, n);
    for (i = 0; i < nM; i++) {
        RVector rm = QVec_Gaussian (mesh, mesh.M[i], mwidth, SRCMODE_NEUMANN);
	FVector fm(n);
	for (j = 0; j < n; j++) fm[j] = (float)rm[j];
	for (j = 0; j < n; j++) fm[j] *= mesh.plist[j].C2A();
	mvec.SetRow (i, fm);
    }

    fws.Allocate (mesh);
    fws.Reset (sol, 0);
    IterativeSolverResult res;
    t0 = clock();
    fws.CalcFields (qvec, dphi, &res);
    dt = (double)(clock()-t0)/(double)CLOCKS_PER_SEC;

    FVector proj = fws.ProjectSingle (0, mvec, dphi[0]);
    cout << "Max. iterations:  " << res.it_count << endl;
    cout << "Max. rel. error:  " << res.rel_err << endl;
    cout << "Time(CalcFields): " << dt << endl;

    delete []dphi;
}

// =================================================================

void Test_FwdSolver_double (QMMesh &mesh)
{
    cout << "CUDA test: FwdSolver" << endl;

    int i, j;
    double tol = 1e-8;
    double c = c0/refind;
    double qwidth = 2;
    double mwidth = 2;
    SourceMode srctp = SRCMODE_NEUMANN;
    int n = mesh.nlen();
    int nQ = mesh.nQ;
    int nM = mesh.nM;

    RFwdSolver fws("BICGSTAB", tol);
    fws.SetLinSolverMaxit (200);
    RCompRowMatrix qvec, mvec;
    RVector *dphi;
    double t0, dt;

    Solution sol(OT_NPARAM, n);
    RVector cmua(n); cmua = c*0.01;
    RVector ckap(n); ckap = c*(1.0/(3.0*(0.01+1)));
    RVector ref(n);  ref  = refind;
    RVector c2a(n);  c2a  = c0/(2*refind*A_Keijzer(refind));
    sol.SetParam (OT_CMUA, cmua);
    sol.SetParam (OT_CKAPPA, ckap);
    sol.SetParam (OT_N, ref);
    sol.SetParam (OT_C2A, c2a);

    // build the field vectors
    dphi = new RVector[nQ];
    for (i = 0; i < nQ; i++) dphi[i].New(n);

    // build source vectors
    qvec.New (nQ, n);
    for (i = 0; i < nQ; i++) {
	RVector rq = QVec_Gaussian (mesh, mesh.Q[i], qwidth, srctp);
	qvec.SetRow (i, rq);
    }

    // build measurement vectors
    mvec.New (nM, n);
    for (i = 0; i < nM; i++) {
        RVector rm = QVec_Gaussian (mesh, mesh.M[i], mwidth, SRCMODE_NEUMANN);
	for (j = 0; j < n; j++) rm[j] *= mesh.plist[j].C2A();
	mvec.SetRow (i, rm);
    }

    fws.Allocate (mesh);
    fws.Reset (sol, 0);
    t0 = clock();
    fws.CalcFields (qvec,dphi);
    dt = (double)(clock()-t0)/(double)CLOCKS_PER_SEC;

    RVector proj = fws.ProjectSingle (0, mvec, dphi[0]);
    cout << proj << endl;

    //for (i = 0; i < nQ; i++)
    //  cout << "||phi" << i << "||: " << l2norm(dphi[i]) << endl;

    cout << "Time(CalcFields): " << dt << endl;

    delete []dphi;
}

// =================================================================

void Test_FwdSolver_scomplex (QMMesh &mesh)
{
    cout << "CUDA test: FwdSolver" << endl;

    int i, j;
    double tol = 1e-8;
    double c = c0/refind;
    double qwidth = 2;
    double mwidth = 2;
    SourceMode srctp = SRCMODE_NEUMANN;
    int n = mesh.nlen();
    int nQ = mesh.nQ;
    int nM = mesh.nM;

    double freq = 100;
    double omega = freq * 2.0*Pi*1e-6;

    SCFwdSolver fws("BICGSTAB", tol);
    fws.SetLinSolverMaxit (200);
    SCCompRowMatrix qvec, mvec;
    SCVector *dphi;
    double t0, dt;

    Solution sol(OT_NPARAM, n);
    RVector cmua(n); cmua = c*0.01;
    RVector ckap(n); ckap = c*(1.0/(3.0*(0.01+1)));
    RVector ref(n);  ref  = refind;
    RVector c2a(n);  c2a  = c0/(2*refind*A_Keijzer(refind));
    sol.SetParam (OT_CMUA, cmua);
    sol.SetParam (OT_CKAPPA, ckap);
    sol.SetParam (OT_N, ref);
    sol.SetParam (OT_C2A, c2a);

    // build the field vectors
    dphi = new SCVector[nQ];
    for (i = 0; i < nQ; i++) dphi[i].New(n);

    // build source vectors
    qvec.New (nQ, n);
    for (i = 0; i < nQ; i++) {
	RVector rq = QVec_Gaussian (mesh, mesh.Q[i], qwidth, srctp);
	SCVector fq(n);
	for (j = 0; j < n; j++) fq[j].re = (float)rq[j];
	qvec.SetRow (i, fq);
    }

    // build measurement vectors
    mvec.New (nM, n);
    for (i = 0; i < nM; i++) {
        RVector rm = QVec_Gaussian (mesh, mesh.M[i], mwidth, SRCMODE_NEUMANN);
	SCVector fm(n);
	for (j = 0; j < n; j++) fm[j].re = (float)rm[j];
	for (j = 0; j < n; j++) fm[j].re *= mesh.plist[j].C2A();
	mvec.SetRow (i, fm);
    }

    fws.Allocate (mesh);
    fws.Reset (sol, omega);
    t0 = clock();
    fws.CalcFields (qvec,dphi);
    dt = (double)(clock()-t0)/(double)CLOCKS_PER_SEC;

    SCVector proj = fws.ProjectSingle (0, mvec, dphi[0]);
    cout << proj << endl;

    //for (i = 0; i < nQ; i++)
    //  cout << "||phi" << i << "||: " << l2norm(dphi[i]) << endl;

    cout << "Time(CalcFields): " << dt << endl;

    delete []dphi;
}

// =================================================================

int main (int argc, char *argv[])
{
#ifdef USE_CUDA_FLOAT
    cuda_EchoDeviceProperties();
#endif

    QMMesh mesh;
    ifstream ifs (argv[1]);
    ifs >> mesh;
    mesh.Setup();

    ifstream qmf (argv[2]);
    mesh.LoadQM (qmf);

    int iq = 0;
    if (argc > 3) sscanf (argv[3], "%d", &iq);

    //Test_Vector();
    Test_Ax_float (mesh);
    Test_Ax_double (mesh);
    //Test_BiCGSTAB (mesh, iq);
    //Test_FwdSolver (mesh);
    //Test_FwdSolver_double (mesh);
    //Test_FwdSolver_scomplex (mesh);

    return 0;
}
