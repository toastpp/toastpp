//#include <complex.h>

#define SQR(X) ( (X)*(X))
#define lmYindex(l,m) ( abs(m) <= l ? l*(l+1) + m : -1 )
#define lind(k) (int(floor(sqrt(double(k)))))
#define mind(k) (int(k - SQR(lind(k)) - lind(k) )) 
#define CGzm(l,m)  (sqrt((double((l+m)*(l-m)))/(double((2*l+1)*(2*l-1))) ))
#define CGzp(l,m)  (sqrt((double((l+m+1)*(l-m+1)))/(double((2*l+1)*(2*l+3))) ))
#define CGpem(l,m) (sqrt((double((l-m)*(l-m-1)))/(double((2*l+1)*(2*l-1))) ))
#define CGpep(l,m) (-sqrt((double((l+m+1)*(l+m+2)))/(double((2*l+1)*(2*l+3)))))
#define CGmem(l,m) (-sqrt((double((l+m)*(l+m-1)))/(double((2*l+1)*(2*l-1))) ))
#define CGmep(l,m) (sqrt((double((l-m+1)*(l-m+2)))/(double((2*l+1)*(2*l+3))) ))

/* now the raising and lowering second order coefficient */
#define CGzzm(l,m) (CGzm(l,m))*(CGzm((l-1),m))
#define CGzz0(l,m) (CGzm(l,m))*(CGzp((l-1),m))
#define CGzz1(l,m) (CGzp(l,m))*(CGzm((l+1),m))
#define CGzzp(l,m) (CGzp(l,m))*(CGzp((l+1),m))

#define CGaa(l,m) (CGpep(l,m))*(CGpep((l+1),(m+1)))
#define CGab(l,m) (CGpep(l,m))*(CGpem((l+1),(m+1)))
#define CGac(l,m) (CGpep(l,m))*(CGmep((l+1),(m+1)))
#define CGad(l,m) (CGpep(l,m))*(CGmem((l+1),(m+1)))
//
#define CGba(l,m) (CGpem(l,m))*(CGpep((l-1),(m+1)))
#define CGbb(l,m) (CGpem(l,m))*(CGpem((l-1),(m+1)))
#define CGbc(l,m) (CGpem(l,m))*(CGmep((l-1),(m+1)))
#define CGbd(l,m) (CGpem(l,m))*(CGmem((l-1),(m+1)))
//
#define CGca(l,m) (CGmep(l,m))*(CGpep((l-1),(m+1)))
#define CGcb(l,m) (CGmep(l,m))*(CGpem((l-1),(m+1)))
#define CGcc(l,m) (CGmep(l,m))*(CGmep((l-1),(m+1)))
#define CGcd(l,m) (CGmep(l,m))*(CGmem((l-1),(m+1)))
//
#define CGda(l,m) (CGmep(l,m))*(CGpep((l-1),(m-1)))
#define CGdb(l,m) (CGmep(l,m))*(CGpem((l-1),(m-1)))
#define CGdc(l,m) (CGmep(l,m))*(CGmep((l-1),(m-1)))
#define CGdd(l,m) (CGmep(l,m))*(CGmem((l-1),(m-1)))



using namespace toast;

const double irtpi = 1/sqrt(M_PI);
const double irt2 = sqrt(0.5);
const toast::complex iirt2 = sqrt( toast::complex(-0.5,0));
const toast::complex myi =  toast::complex(0.0,1.0);

RDenseMatrix Chop(const RDenseMatrix& ) ;
RDenseMatrix RealMat(const CDenseMatrix& ) ;
RDenseMatrix ImagMat(const CDenseMatrix& ) ;
RDenseMatrix Sparse2Dense(const RCompRowMatrix& );
CDenseMatrix Sparse2Dense(const CCompRowMatrix& );

RCompRowMatrix Chop(RCompRowMatrix& ); // TOAST complains if argument is const
RCompRowMatrix RealMat(CCompRowMatrix& );
RCompRowMatrix ImagMat(CCompRowMatrix& );

void genmat_angint_3D_PN(RCompRowMatrix&, RCompRowMatrix&, RCompRowMatrix&, RCompRowMatrix&, RCompRowMatrix&, CCompRowMatrix&, const int);
void genmat_angint_sdm_3D_PN(RCompRowMatrix&, RCompRowMatrix&, RCompRowMatrix&,  RCompRowMatrix&, RCompRowMatrix&, RCompRowMatrix&, RCompRowMatrix&, RCompRowMatrix&, RCompRowMatrix&, CCompRowMatrix&, const int);
