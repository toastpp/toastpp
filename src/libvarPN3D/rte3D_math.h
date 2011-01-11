#include<mathlib.h>
double factorial(int n);

double doublefactorial(int n);

toast::complex cpowi(toast::complex &c, const int m);

int sign(int m);

int signum(int m);

toast::complex wigCoeff(const int l1, const int l2, const int m1, const int m2);

toast::complex kronD(const IVector &a, const IVector &b);

//double plm(const int l, const int m, const double x);
RVector plm(const int l, const int m, const RVector& vec);

void LegendreTable(const int l, const int numpts, const RDenseMatrix &pt, RDenseMatrix* &LT);

void sphericalHarmonics(const int order, const int numpts, const RDenseMatrix& pt, RDenseMatrix* &Ylm);

//CDenseMatrix sphericalHarmonicsConj(const int order, const int numpts, const RDenseMatrix& pt);

CCompRowMatrix kronsd(const CCompRowMatrix &A, const CDenseMatrix &B);
void kronplus(const int spatrow, const int spatcol, const IVector& node_angN, const IVector& offset, const double Aval, const RCompRowMatrix &B, RCompRowMatrix& C);
//void BIntUnitSphere(const int sphOrder1, const int size, const RDenseMatrix& ptsPlus, const CVector& wtsPlus, const RVector& bnormal, CCompRowMatrix &Aintsc, CCompRowMatrix &Aintss, CCompRowMatrix &Aintc, CCompRowMatrix& intSdotnPlusHemisphere, CCompRowMatrix& intSdotnMinusHemisphere);
void BIntUnitSphere(const int size1, const int size2, const int sphOrder1, const int sphOrder1, const RDenseMatrix& ptsPlus, const RVector& wtsPlus, const RVector& bnormal, RDenseMatrix* &Ylm, RCompRowMatrix& intSdotnPlusHemisphere, RCompRowMatrix& intSdotnMinusHemisphere);


/*void BIntUnitSphere(const int sphOrder1, const RDenseMatrix& ptsPlus, const CVector& wtsPlus, const RVector& bnormal, CDenseMatrix& bintplus, CDenseMatrix& bintminus);

void BIntUnitSphere_plus(const int size, const int sphOrder, const RDenseMatrix& pts, const CVector &wts, const RVector& bnormal,  CCompRowMatrix& bintplus);

void BIntUnitSphere_minus(const int size, const int sphOrder, const RDenseMatrix& pts, const CVector &wts, const RVector& bnormal, const bool is_isotropic, const RVector& dirVec, CCompRowMatrix& bintminus);

void genmat_angint_3D_test(const int size, const int sphOrder,  const RDenseMatrix& pts, const CVector& wts, CDenseMatrix& Aintssss, CDenseMatrix& Aintscc, CDenseMatrix& Aintssc, CDenseMatrix& Aintcc);

void getQuadPointsWts(const int order, RDenseMatrix& pts, CVector& wts);
*/
int getPos(const int l, const int m);

RCompRowMatrix shrink(const RDenseMatrix &dnsmat);

void tabulatePtsWtsForSphere();
void testPtsWtsForSphere(const int sphOrder, const  int angN, const RDenseMatrix& pts, const CVector &wts, CDenseMatrix* &Yl2m2, CDenseMatrix* &Yl1m1);
/*void testQuadHemisphere();
void testQuadHemisphere2();

void BIntUnitSphere_test(const int size, const int sphOrder);
*/



