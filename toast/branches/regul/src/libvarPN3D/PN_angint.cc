#include <mathlib.h>
#include <stdio.h>
#include <felib.h>
#include "PN_incl.h"
#include "PN_angint.h"
#include "sphericalHarmonic_algebra.h"
#include "toast.h"
using namespace toast;

double intSinCosY(const int l, const int m)
{	
	toast::complex result(0, 0);
	CVector a(2), b(2), c(2), d(2);
	IDenseMatrix ac(2, 2), bc(2, 2), cc(2, 2), dc(2, 2); 
	sincosY(l, m, a, b, c, d, ac, bc, cc, dc);
	if(ac(0, 0) == 0 && ac(0, 1) == 0) result += a[0];
	if(ac(1, 0) == 0 && ac(1, 1) == 0) result += a[1]; 
	if(bc(0, 0) == 0 && bc(0, 1) == 0) result += b[0];
	if(bc(1, 0) == 0 && bc(1, 1) == 0) result += b[1]; 
	if(cc(0, 0) == 0 && cc(0, 1) == 0) result += c[0];
	if(cc(1, 0) == 0 && cc(1, 1) == 0) result += c[1];  
	if(dc(0, 0) == 0 && dc(0, 1) == 0) result += d[0];
	if(dc(1, 0) == 0 && dc(1, 1) == 0) result += d[1]; 
	return(result.re); 
}

double intSinSinY(const int l, const int m)
{	
	toast::complex result(0.0, 0.0);
	CVector a(2), b(2), c(2), d(2);
	IDenseMatrix ac(2, 2), bc(2, 2), cc(2, 2), dc(2, 2); 
	sinsinY(l, m, a, b, c, d, ac, bc, cc, dc);
	if(ac(0, 0) == 0 && ac(0, 1) == 0) result += a[0];
	if(ac(1, 0) == 0 && ac(1, 1) == 0) result += a[1]; 
	if(bc(0, 0) == 0 && bc(0, 1) == 0) result += b[0];
	if(bc(1, 0) == 0 && bc(1, 1) == 0) result += b[1]; 
	if(cc(0, 0) == 0 && cc(0, 1) == 0) result += c[0];
	if(cc(1, 0) == 0 && cc(1, 1) == 0) result += c[1];  
	if(dc(0, 0) == 0 && dc(0, 1) == 0) result += d[0];
	if(dc(1, 0) == 0 && dc(1, 1) == 0) result += d[1]; 
	return(result.re); 
}

double intCosY(const int l, const int m)
{	
	toast::complex result(0, 0);
	CVector e(2), f(2);
	IDenseMatrix ec(2, 2), fc(2, 2); 
	cosY(l, m, e, f, ec, fc);
	if(ec(0, 0) == 0 && ec(0, 1) == 0) result += e[0];
	if(ec(1, 0) == 0 && ec(1, 1) == 0) result += e[1]; 
	if(fc(0, 0) == 0 && fc(0, 1) == 0) result += f[0];
	if(fc(1, 0) == 0 && fc(1, 1) == 0) result += f[1]; 
	return(result.re); 
}

/**
Computes all the angular integrals required by variable order PN approximation
**/
void genmat_angint(const int sphOrder, const int angN, RCompRowMatrix& Aint, RCompRowMatrix& Aintsc, RCompRowMatrix& Aintss, RCompRowMatrix& Aintc, RCompRowMatrix& Aintscsc, RCompRowMatrix& Aintscss, RCompRowMatrix& Aintscc, RCompRowMatrix& Aintssss, RCompRowMatrix& Aintssc, RCompRowMatrix& Aintcc)
{
	RDenseMatrix dnsAint(angN, angN), dnsAintsc(angN, angN), dnsAintss(angN, angN), dnsAintc(angN, angN);
	RDenseMatrix dnsAintscsc(angN, angN), dnsAintscss(angN, angN), dnsAintcc(angN, angN);
	RDenseMatrix dnsAintscc(angN, angN), dnsAintssss(angN, angN), dnsAintssc(angN, angN);

	CVector a1(2), b1(2), c1(2), d1(2), e1(2), f1(2), p1(2), a2(2), b2(2), c2(2), d2(2), e2(2), f2(2), p2(2);
        IDenseMatrix a1c(2, 2), b1c(2, 2), c1c(2, 2), d1c(2, 2), e1c(2, 2), f1c(2, 2), p1c(2, 2);
	IDenseMatrix a2c(2, 2), b2c(2, 2), c2c(2, 2), d2c(2, 2), e2c(2, 2), f2c(2, 2), p2c(2, 2);
	 

	int is, js, indl1, indl2;

	for(int l1 = 0; l1 <= sphOrder; l1++){
		indl1 = l1*l1;
		for(int m1 = -1*l1; m1 <= l1; m1++){
			is = indl1 + l1 + m1;    	
			for(int l2 = 0; l2 <= sphOrder; l2++){
				indl2 = l2*l2;
				for(int m2 = -1*l2; m2 <= l2; m2++){
					
					js = indl2 + l2 + m2;
				
					sphY(l1, m1, p1, p1c);
					sphY(l2, m2, p2, p2c);
					dnsAint(is, js) = Integrate2(p1, p2, p1c, p2c);

					sincosY(l2, m2, a2, b2, c2, d2, a2c, b2c, c2c, d2c);
					dnsAintsc(is, js) = Integrate2(p1, a2, p1c, a2c) + Integrate2(p1, b2, p1c, b2c);
					dnsAintsc(is, js) += Integrate2(p1, c2, p1c, c2c) + Integrate2(p1, d2, p1c, d2c);

					sinsinY(l2, m2, a2, b2, c2, d2, a2c, b2c, c2c, d2c);
					dnsAintss(is, js) = Integrate2(p1, a2, p1c, a2c) + Integrate2(p1, b2, p1c, b2c);
					dnsAintss(is, js) += Integrate2(p1, c2, p1c, c2c) + Integrate2(p1, d2, p1c, d2c);

					cosY(l2, m2, e2, f2, e2c, f2c);
					dnsAintc(is, js) = Integrate2(p1, e2, p1c, e2c) + Integrate2(p1, f2, p1c, f2c);

					sincosY(l1, m1, a1, b1, c1, d1, a1c, b1c, c1c, d1c);
					sincosY(l2, m2, a2, b2, c2, d2, a2c, b2c, c2c, d2c);
					dnsAintscsc(is, js) = Integrate2(a1, a2, a1c, a2c) + Integrate2(a1, b2, a1c, b2c);
					dnsAintscsc(is, js) += Integrate2(a1, c2, a1c, c2c) + Integrate2(a1, d2, a1c, d2c);
					dnsAintscsc(is, js) += Integrate2(b1, a2, b1c, a2c) + Integrate2(b1, b2, b1c, b2c);
					dnsAintscsc(is, js) += Integrate2(b1, c2, b1c, c2c) + Integrate2(b1, d2, b1c, d2c);
					dnsAintscsc(is, js) += Integrate2(c1, a2, c1c, a2c) + Integrate2(c1, b2, c1c, b2c);
					dnsAintscsc(is, js) += Integrate2(c1, c2, c1c, c2c) + Integrate2(c1, d2, c1c, d2c);
					dnsAintscsc(is, js) += Integrate2(d1, a2, d1c, a2c) + Integrate2(d1, b2, d1c, b2c);
					dnsAintscsc(is, js) += Integrate2(d1, c2, d1c, c2c) + Integrate2(d1, d2, d1c, d2c);

					
					sinsinY(l2, m2, a2, b2, c2, d2, a2c, b2c, c2c, d2c);
					dnsAintscss(is, js) = Integrate2(a1, a2, a1c, a2c) + Integrate2(a1, b2, a1c, b2c);
					dnsAintscss(is, js) += Integrate2(a1, c2, a1c, c2c) + Integrate2(a1, d2, a1c, d2c);
					dnsAintscss(is, js) += Integrate2(b1, a2, b1c, a2c) + Integrate2(b1, b2, b1c, b2c);
					dnsAintscss(is, js) += Integrate2(b1, c2, b1c, c2c) + Integrate2(b1, d2, b1c, d2c);
					dnsAintscss(is, js) += Integrate2(c1, a2, c1c, a2c) + Integrate2(c1, b2, c1c, b2c);
					dnsAintscss(is, js) += Integrate2(c1, c2, c1c, c2c) + Integrate2(c1, d2, c1c, d2c);
					dnsAintscss(is, js) += Integrate2(d1, a2, d1c, a2c) + Integrate2(d1, b2, d1c, b2c);
					dnsAintscss(is, js) += Integrate2(d1, c2, d1c, c2c) + Integrate2(d1, d2, d1c, d2c);

					cosY(l2, m2, e2, f2, e2c, f2c);
					dnsAintscc(is, js) = Integrate2(a1, e2, a1c, e2c) + Integrate2(a1, f2, a1c, f2c);
					dnsAintscc(is, js) += Integrate2(b1, e2, b1c, e2c) + Integrate2(b1, f2, b1c, f2c);
					dnsAintscc(is, js) += Integrate2(c1, e2, c1c, e2c) + Integrate2(c1, f2, c1c, f2c);
					dnsAintscc(is, js) += Integrate2(d1, e2, d1c, e2c) + Integrate2(d1, f2, d1c, f2c);

					cosY(l1, m1, e1, f1, e1c, f1c);
					dnsAintcc(is, js) = Integrate2(e1, e2, e1c, e2c) + Integrate2(e1, f2, e1c, f2c);
					dnsAintcc(is, js) += Integrate2(f1, e2, f1c, e2c) + Integrate2(f1, f2, f1c, f2c);


					sinsinY(l1, m1, a1, b1, c1, d1, a1c, b1c, c1c, d1c);
					dnsAintssc(is, js) = Integrate2(a1, e2, a1c, e2c) + Integrate2(a1, f2, a1c, f2c);
					dnsAintssc(is, js) += Integrate2(b1, e2, b1c, e2c) + Integrate2(b1, f2, b1c, f2c);
					dnsAintssc(is, js) += Integrate2(c1, e2, c1c, e2c) + Integrate2(c1, f2, c1c, f2c);
					dnsAintssc(is, js) += Integrate2(d1, e2, d1c, e2c) + Integrate2(d1, f2, d1c, f2c);

					sinsinY(l2, m2, a2, b2, c2, d2, a2c, b2c, c2c, d2c);
					dnsAintssss(is, js) = Integrate2(a1, a2, a1c, a2c) + Integrate2(a1, b2, a1c, b2c);
					dnsAintssss(is, js) += Integrate2(a1, c2, a1c, c2c) + Integrate2(a1, d2, a1c, d2c);
					dnsAintssss(is, js) += Integrate2(b1, a2, b1c, a2c) + Integrate2(b1, b2, b1c, b2c);
					dnsAintssss(is, js) += Integrate2(b1, c2, b1c, c2c) + Integrate2(b1, d2, b1c, d2c);
					dnsAintssss(is, js) += Integrate2(c1, a2, c1c, a2c) + Integrate2(c1, b2, c1c, b2c);
					dnsAintssss(is, js) += Integrate2(c1, c2, c1c, c2c) + Integrate2(c1, d2, c1c, d2c);
					dnsAintssss(is, js) += Integrate2(d1, a2, d1c, a2c) + Integrate2(d1, b2, d1c, b2c);
					dnsAintssss(is, js) += Integrate2(d1, c2, d1c, c2c) + Integrate2(d1, d2, d1c, d2c);
					}
			}
		}
	}
	Aint = shrink(dnsAint);Aintsc = shrink(dnsAintsc); Aintss = shrink(dnsAintss);Aintc = shrink(dnsAintc);
	Aintscsc = shrink(dnsAintscsc);Aintscss = shrink(dnsAintscss);Aintscc = shrink(dnsAintscc); Aintssc = shrink(dnsAintssc);
	Aintssss = shrink(dnsAintssss);Aintcc = shrink(dnsAintcc);
}

/**
Phase function discretization
NOTE!! Only Henyey-Greenstein phase function has been implemented currently
**/
RVector phaseFuncDisc(int sphOrder, toast::complex (*phaseFunc)(const double g, const double costheta), const double g)
{
 RVector phaseFn(sphOrder+1);
 for(int l=0; l <= sphOrder; l++)
		phaseFn[l] = pow(g, (double)l);
 return(phaseFn);	
}

/** Phase integrals 
**/
void genmat_apu(toast::complex (*phaseFunc)(const double g, const double costheta), const double g, const int angN, const int sphOrder, RCompRowMatrix& apu1, RCompRowMatrix& apu1sc, RCompRowMatrix& apu1ss, RCompRowMatrix& apu1c)
{
 RDenseMatrix dnsapu1(angN, angN), dnsapu1sc(angN, angN), dnsapu1ss(angN, angN), dnsapu1c(angN, angN);
     RVector phaseFn;
    phaseFn = phaseFuncDisc(sphOrder, phaseFunc, g);
    std::cout<<"****************** Discretized phase function ****************"<<std::endl;
    std::cout<<phaseFn<<std::endl;
    std::cout<<"**************************************************************"<<std::endl;
    int indl, indm;
    int indl1, indm1, indl2, indm2, is, js;
    CVector a1(2), b1(2), c1(2), d1(2), e1(2), f1(2), p1(2), p2(2), p(2);
    IDenseMatrix a1c(2, 2), b1c(2, 2), c1c(2, 2), d1c(2, 2), e1c(2, 2), f1c(2, 2), p1c(2, 2), p2c(2, 2), pc(2, 2);
 
    for(int l1=0; l1 <= sphOrder; l1++){
    	indl1 =  l1*l1;
	for(int m1 = -1*l1; m1 <= l1; m1++){
		indm1 = l1 + m1;
		is = indl1 + indm1;
	        for(int l2=0; l2<= sphOrder; l2++){
			indl2 = l2*l2;
			for(int m2 = -1*l2; m2<=l2; m2++){
				indm2 = l2 + m2;
				js = indl2 + indm2;
				
				sphY(l1, m1, p1, p1c);
				sphY(l2, m2, p2, p2c);
				dnsapu1(is, js) += phaseFn[l2]*Integrate2(p2, p1, p2c, p1c);
				
				sincosY(l1, m1, a1, b1, c1, d1, a1c, b1c, c1c, d1c);
				dnsapu1sc(is, js) += phaseFn[l2]*(Integrate2(a1, p2, a1c, p2c) + Integrate2(b1, p2, b1c, p2c) + Integrate2(c1, p2, c1c, p2c) + Integrate2(d1, p2, d1c, p2c));
				
				sinsinY(l1, m1, a1, b1, c1, d1, a1c, b1c, c1c, d1c);
				dnsapu1ss(is, js) += phaseFn[l2]*(Integrate2(a1, p2, a1c, p2c) + Integrate2(b1, p2, b1c, p2c) + Integrate2(c1, p2, c1c, p2c) + Integrate2(d1, p2, d1c, p2c));

				cosY(l1, m1, e1, f1, e1c, f1c);
				dnsapu1c(is, js) += phaseFn[l2]*(Integrate2(e1, p2, e1c, p2c) + Integrate2(f1, p2, f1c, p2c));

			
			}
		}
	}
   }

  apu1 = shrink(dnsapu1); apu1sc = shrink(dnsapu1sc); apu1ss = shrink(dnsapu1ss);
  apu1c = shrink(dnsapu1c);
}

