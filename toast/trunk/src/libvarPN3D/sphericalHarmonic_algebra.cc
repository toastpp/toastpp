#include "sphericalHarmonic_algebra.h"
using namespace toast;
double fac[] = {1.000000e+000, 1.000000e+000, 2.000000e+000, 6.000000e+000, 2.400000e+001, 1.200000e+002, 7.200000e+002, 5.040000e+003, 4.032000e+004, 3.628800e+005, 3.628800e+006, 3.991680e+007, 4.790016e+008, 6.227021e+009, 8.717829e+010, 1.307674e+012, 2.092279e+013, 3.556874e+014, 6.402374e+015, 1.216451e+017, 2.432902e+018, 5.109094e+019, 1.124001e+021, 2.585202e+022, 6.204484e+023, 1.551121e+025, 4.032915e+026, 1.088887e+028, 3.048883e+029, 8.841762e+030, 2.652529e+032, 8.222839e+033, 2.631308e+035, 8.683318e+036, 2.952328e+038, 1.033315e+040, 3.719933e+041, 1.376375e+043, 5.230226e+044, 2.039788e+046, 8.159153e+047, 3.345253e+049, 1.405006e+051, 6.041526e+052, 2.658272e+054, 1.196222e+056, 5.502622e+057, 2.586232e+059, 1.241392e+061, 6.082819e+062, 3.041409e+064, 1.551119e+066, 8.065818e+067, 4.274883e+069, 2.308437e+071, 1.269640e+073, 7.109986e+074, 4.052692e+076, 2.350561e+078, 1.386831e+080, 8.320987e+081, 5.075802e+083, 3.146997e+085, 1.982608e+087, 1.268869e+089, 8.247651e+090, 5.443449e+092, 3.647111e+094, 2.480036e+096, 1.711225e+098, 1.197857e+100, 8.504786e+101, 6.123446e+103, 4.470115e+105, 3.307885e+107, 2.480914e+109, 1.885495e+111, 1.451831e+113, 1.132428e+115, 8.946182e+116, 7.156946e+118, 5.797126e+120, 4.753643e+122, 3.945524e+124, 3.314240e+126, 2.817104e+128, 2.422710e+130, 2.107757e+132, 1.854826e+134, 1.650796e+136, 1.485716e+138, 1.352002e+140, 1.243841e+142, 1.156773e+144, 1.087366e+146, 1.032998e+148, 9.916779e+149, 9.619276e+151, 9.426890e+153, 9.332622e+155, 9.332622e+157};

double dfactorial[] = {1.000000e+000, 1.000000e+000, 2.000000e+000, 3.000000e+000, 8.000000e+000, 1.500000e+001, 4.800000e+001, 1.050000e+002, 3.840000e+002, 9.450000e+002, 3.840000e+003, 1.039500e+004, 4.608000e+004, 1.351350e+005, 6.451200e+005, 2.027025e+006, 1.032192e+007, 3.445943e+007, 1.857946e+008, 6.547291e+008, 3.715891e+009, 1.374931e+010, 8.174961e+010, 3.162341e+011, 1.961991e+012, 7.905854e+012, 5.101175e+013, 2.134580e+014, 1.428329e+015, 6.190283e+015, 4.284987e+016, 1.918988e+017, 1.371196e+018, 6.332660e+018, 4.662066e+019, 2.216431e+020, 1.678344e+021, 8.200795e+021, 6.377707e+022, 3.198310e+023, 2.551083e+024, 1.311307e+025, 1.071455e+026, 5.638620e+026, 4.714401e+027, 2.537379e+028, 2.168624e+029, 1.192568e+030, 1.040940e+031, 5.843584e+031, 5.204698e+032, 2.980228e+033, 2.706443e+034, 1.579521e+035, 1.461479e+036, 8.687364e+036, 8.184284e+037, 4.951798e+038, 4.746885e+039, 2.921561e+040, 2.848131e+041, 1.782152e+042, 1.765841e+043, 1.122756e+044, 1.130138e+045, 7.297912e+045, 7.458913e+046, 4.889601e+047, 5.072061e+048, 3.373825e+049, 3.550443e+050, 2.395416e+051, 2.556319e+052, 1.748653e+053, 1.891676e+054, 1.311490e+055, 1.437674e+056, 1.009847e+057, 1.121385e+058, 7.977794e+058, 8.971083e+059, 6.462013e+060, 7.356288e+061, 5.363471e+062, 6.179282e+063, 4.558950e+064, 5.314183e+065, 3.966287e+066, 4.676481e+067, 3.529995e+068, 4.208833e+069, 3.212296e+070, 3.872126e+071, 2.987435e+072, 3.639799e+073, 2.838063e+074, 3.494207e+075, 2.752921e+076, 3.424322e+077, 2.725392e+078, 3.424322e+079};

/** factorial
**/
double factorial(int n){
     xASSERT(n >= 0, "Factorial is defined only for non-negative integers");
     if(n <= 100) return(fac[n]);
     else return(n * factorial(n-1));
}

/** double factorial
**/
double doublefactorial(int n){
	if(n >= 0 & n<=100) return(dfactorial[n]);
	else if(n == -1) return(1);
	else return(n*doublefactorial(n-2));
}
/** Computes C^{n} where C is a 'complex' number and 'n' is a positive integer
**/
toast::complex cpowi(toast::complex &c, const int m){
        xASSERT(m >= 0, "cpowi is defined only for non-negative exponents");
	toast::complex ans(1, 0);
	for(int i = 0; i < m; i++) ans = ans*c;

	return(ans);
}

/** Computes (-1)^{m}
**/
int sign(int m){
	if(abs(m) % 2 == 0) return(1);
	else return(-1);
}	

/** Returns 0, 1, and 2 respectively based on whether m=0, m<0 and m>0.
**/
int signum(int m){
	if(m == 0)
		return(0);
	else if(m < 0)
		return(1);
	else
		return(2);
}

/* evaluates \int_{S^{2}} Y_{\ell_{1}, m_{1}} Y_{\ell_{2}, m_{2}} where both Y's are complex spherical harmonics
*/
toast::complex kronD(const IVector &a, const IVector &b)
{
        int l1 = a[0], l2 = b[0], m1 = a[1], m2 = b[1];
	if(l1 >= 0 && l2>=0 && abs(m1) <= l1 && abs(m2) <= l2 && l1==l2 && m1 == -1*m2) return(toast::complex(sign(m2), 0.0));
	else return(toast::complex(0.0, 0.0));
}

/** Computes associated Legendre polynomials of a given order on a set of points 
*	l -> Maximum order
*	numpts -> number of points
*	pt -> the three dimensional point set
*	LT[] -> Legendre polynomials in a table form
**/
void LegendreTable(const int l, const int numpts, const RDenseMatrix &pt, RDenseMatrix* &LT){
	int j, k, r, s;
    
	RVector x = pt.Col(2);
        RVector ons(x.Dim(), 1.0);
	for(int i=0; i <= l; i++)
	{
		LT[i].SetRow(2*i, pow(ons- x*x, i/2.0)*sign(i)*doublefactorial(2*i-1)); //m= l
	 	if(i>0)
			LT[i].SetRow(2*i-1, x*pow(ons - x*x, (i-1)/2.0)*(2*i-1)*sign(i-1)*doublefactorial(2*i-3)); //m=l-1
		
		k = 2*i-2;r = 2*(i-1)-1; s = 2*(i-2);
		for(int m = i-2; m >= 0; m--)
		{
			LT[i].SetRow(k, (x*(2*i-1)*LT[i-1].Row(r) - LT[i-2].Row(s)*(i+m-1))/(double)(i-m));
	  		k--; r--; s--;
		}
		j=0; k = 2*i;
		for(int m = -1*i; m < 0; m++)
		{
			LT[i].SetRow(j, LT[i].Row(k)*sign(m)*factorial(i-abs(m))/(double)factorial(i+abs(m)));
			j++; k--; 
		}
	}
}

/**Real spherical harmonics computed on a set of points
* order -> Maximum order of spherical harmonics
* numpts ->  number of points on which to evaluate
* pt -> 3D point set
* sphHarm[l] -> spherical harmonics over point set of order 'l'
**/
void sphericalHarmonics(const int order, const int numpts, const RDenseMatrix& pt, RDenseMatrix* &sphHarm){
	int m, sgm;
	double dtmp1, dtmp2, normfac;
	CDenseMatrix *cplxSphHarm;
	cplxSphHarm = new CDenseMatrix[order+1];
	for(int l=0; l<=order; l++)
	{
		cplxSphHarm[l].New(2*l+1, numpts);
	}
	toast::complex ctmp;
        RDenseMatrix *LT;
	LT = new RDenseMatrix[order+1];
	for(int l=0; l<=order; l++)
		LT[l].New(2*l+1, numpts);

	LegendreTable(order, numpts, pt, LT);
	for(int l=0; l<=order; l++)
	{
		int k = 0;
		for(m = -1*l; m <= l; m++){
	        	sgm = signum(m);
			dtmp1 =  sqrt((2*l+1)*factorial(l - m)/(double)(4*M_PI*factorial(l+m)));
			for(int i=0; i<numpts; i++){
				dtmp2 = dtmp1*LT[l].Get(k, i);
				normfac = sqrt(1-pt.Get(i, 2)*pt.Get(i, 2));
				switch(sgm){
					case 0: 
						cplxSphHarm[l](k, i).re = dtmp2; cplxSphHarm[l](k, i).im = 0.0;
						break;
					case 1:
						if(normfac!=0){
							ctmp.re = pt.Get(i, 0)/normfac; ctmp.im = -1*pt.Get(i, 1)/normfac;
							ctmp = cpowi(ctmp, abs(m));
							cplxSphHarm[l](k, i) = ctmp*toast::complex(dtmp2, 0);
						}
						break;
	       				case 2:
						if(normfac!=0){
							ctmp.re = pt.Get(i, 0)/normfac; ctmp.im = pt.Get(i, 1)/normfac;
							ctmp = cpowi(ctmp, abs(m));
							cplxSphHarm[l](k, i) = ctmp*toast::complex(dtmp2, 0);
						}
						break;
				}
	    		}	
		k++;	
		}
	}
	toast::complex temp; 
        int indm;
	for(int l=0; l<=order; l++)
	{
		for(m=-1*l; m<= l; m++){
			indm = l + m;
			if(m>0)
			{
				for(int i=0; i<numpts; i++)
				{
					temp = (cplxSphHarm[l](m+l, i) + toast::complex(sign(m), 0)*cplxSphHarm[l](l-m, i))/toast::complex(sqrt(2), 0);
					sphHarm[l](indm, i) = temp.re; 
				}
			}
			else if(m<0)
			{
				for(int i=0; i<numpts; i++)
				{
					temp = (cplxSphHarm[l](l+abs(m), i) - toast::complex(sign(m), 0)*cplxSphHarm[l](l-abs(m), i))/toast::complex(0, sqrt(2));
					sphHarm[l](indm, i) = temp.re;	 
				}
			}
			else
			{
				for(int i=0; i<numpts; i++)
				{
					temp = cplxSphHarm[l](indm, i);
					sphHarm[l](indm, i) = temp.re;
				}
			}
		}
	}

	delete []cplxSphHarm;
	delete []LT;
}

/** Y_{l, m}^{R} = a_{m}Y_{l, m} + b_{m}Y_{l, -m} where the supercript 'R' denotes the real-valued spherical harmonics.
* This function gives the value of a_{m}
**/
toast::complex am(int m)
{
	if(m==0) return(toast::complex(1.0, 0));
	else if(m>0) return(toast::complex(1.0/sqrt(2), 0));
	else return(toast::complex(0, sign(m)/sqrt(2)));
}

/** Y_{l, m}^{R} = a_{m}Y_{l, m} + b_{m}Y_{l, -m} where the supercript 'R' denotes the real-valued spherical harmonics.
* This function gives the value of b_{m}
**/
toast::complex bm(int m)
{
	if(m==0) return(toast::complex(0.0, 0));
	else if(m>0) return(toast::complex(sign(m)/sqrt(2), 0));
	else return(toast::complex(0, -1.0/sqrt(2)));
}


/* Sin(\theta)Cos(\phi)Y_{l, m}^{R} = Sin(\theta)Cos(\phi)(a_{m}Y_{l, m} + b_{m}Y^{*}_{l, m}) where the superscript 'R' denote the 
real-valued spherical harmonics.
a_{m} = 1, b_{m} = 0, if m=0
a_{m} = b_{m} = 1/\sqrt{2}, if m>0
a_{m} = 1/\sqrt{-2}, b_{m} = -1/\sqrt{-2}, otherwise

This routine gives the coefficients of terms Y_{l-1, m+1}, Y_{l+1, m+1}, Y_{l-1, m-1}, Y_{l+1, m-1}, Y_{l-1, -m+1}, 
Y_{l+1, -m+1}, Y_{l-1, -m-1} and Y_{l+1, -m-1}.
*/
void sincosY(const int l, const int m, CVector& a, CVector& b, CVector& c, CVector& d, IDenseMatrix& a1c, IDenseMatrix& b1c, IDenseMatrix& c1c, IDenseMatrix& d1c)
{

	toast::complex alpha = toast::complex(0.5, 0);
	a[0] = am(m)*Blm(l, -m)*alpha; b[0] = -am(m)*Blm(l+1, m+1)*alpha; 
	c[0] = -am(m)*Blm(l, m)*alpha; d[0] = am(m)*Blm(l+1, -m+1)*alpha;

	a[1] = bm(m)*Blm(l, m)*alpha; b[1] = -bm(m)*Blm(l+1, -m+1)*alpha; 
	c[1] = -bm(m)*Blm(l, -m)*alpha; d[1] = bm(m)*Blm(l+1, m+1)*alpha;

	a1c(0, 0) = l-1; a1c(0, 1) = m+1; b1c(0, 0) = l+1; b1c(0, 1) = m+1;
	c1c(0, 0) = l-1; c1c(0, 1) = m-1; d1c(0, 0) = l+1; d1c(0, 1) = m-1;

	a1c(1, 0) = l-1; a1c(1, 1) = -m+1; b1c(1, 0) = l+1; b1c(1, 1) = -m+1;
	c1c(1, 0) = l-1; c1c(1, 1) = -m-1; d1c(1, 0) = l+1; d1c(1, 1) = -m-1; 

}
/* Sin(\theta)Sin(\phi)Y_{l, m}^{R} = Sin(\theta)Sin(\phi)(a_{m}Y_{l, m} + b_{m}Y^{*}_{l, m}) where the superscript 'R' denote the 
real-valued spherical harmonics.
a_{m} = 1, b_{m} = 0, if m=0
a_{m} = b_{m} = 1/\sqrt{2}, if m>0
a_{m} = 1/\sqrt{-2}, b_{m} = -1/\sqrt{-2}, otherwise

This routine gives the coefficients of terms Y_{l-1, m+1}, Y_{l+1, m+1}, Y_{l-1, m-1}, Y_{l+1, m-1}, Y_{l-1, -m+1}, 
Y_{l+1, -m+1}, Y_{l-1, -m-1} and Y_{l+1, -m-1}. 
*/
void sinsinY(const int l, const int m, CVector& a, CVector& b, CVector& c, CVector& d, IDenseMatrix& a1c, IDenseMatrix& b1c, IDenseMatrix& c1c, IDenseMatrix& d1c)
{

	toast::complex alpha = toast::complex(0, -0.5);
	a[0] = am(m)*Blm(l, -m)*alpha; b[0] = -am(m)*Blm(l+1, m+1)*alpha; 
	c[0] = am(m)*Blm(l, m)*alpha; d[0] = -am(m)*Blm(l+1, -m+1)*alpha;

	a[1] = bm(m)*Blm(l, m)*alpha; b[1] = -bm(m)*Blm(l+1, -m+1)*alpha; 
	c[1] = bm(m)*Blm(l, -m)*alpha; d[1] = -bm(m)*Blm(l+1, m+1)*alpha;

	a1c(0, 0) = l-1; a1c(0, 1) = m+1; b1c(0, 0) = l+1; b1c(0, 1) = m+1;
	c1c(0, 0) = l-1; c1c(0, 1) = m-1; d1c(0, 0) = l+1; d1c(0, 1) = m-1;

	a1c(1, 0) = l-1; a1c(1, 1) = -m+1; b1c(1, 0) = l+1; b1c(1, 1) = -m+1;
	c1c(1, 0) = l-1; c1c(1, 1) = -m-1; d1c(1, 0) = l+1; d1c(1, 1) = -m-1; 
}

/* Cos(\theta)Y_{l, m}^{R} = Cos(\theta)(a_{m}Y_{l, m} + b_{m}Y^{*}_{l, m}) where the superscript 'R' denote the 
real-valued spherical harmonics.
a_{m} = 1, b_{m} = 0, if m=0
a_{m} = b_{m} = 1/\sqrt{2}, if m>0
a_{m} = 1/\sqrt{-2}, b_{m} = -1/\sqrt{-2}, otherwise

This routine gives the coefficients of terms Y_{l-1, m}, Y_{l+1, m}, Y_{l-1, -m}, Y_{l+1, -m}. 
*/
void cosY(const int l, const int m, CVector& e, CVector& f, IDenseMatrix& e1c, IDenseMatrix& f1c)
{

	e[0] = am(m)*Alm(l, m); f[0] = am(m)*Alm(l+1, m);
	e[1] = bm(m)*Alm(l, -m); f[1] = bm(m)*Alm(l+1, -m);

	e1c(0, 0) = l-1; e1c(0, 1) = m;
	f1c(0, 0) = l+1; f1c(0, 1) = m;

	e1c(1, 0) = l-1; e1c(1, 1) = -m;
	f1c(1, 0) = l+1; f1c(1, 1) = -m;
}

/* Y_{l, m}^{R} = (a_{m}Y_{l, m} + b_{m}Y^{*}_{l, m}) where the superscript 'R' denote the 
real-valued spherical harmonics.
a_{m} = 1, b_{m} = 0, if m=0
a_{m} = b_{m} = 1/\sqrt{2}, if m>0
a_{m} = 1/\sqrt{-2}, b_{m} = -1/\sqrt{-2}, otherwise

This routine gives the coefficients of terms Y_{l, m}, Y_{l, -m}. 
*/
void sphY(const int l, const int m, CVector& p, IDenseMatrix& p1c)
{

	p[0] = am(m); p[1] = bm(m);
	
	p1c(0, 0) = l; p1c(0, 1) = m;
	p1c(1, 0) = l; p1c(1, 1) = -m;

}

/*
Computes integral of the form
\int (a_{0}Y_{l1, m1} + a_{1}Y_{l2, m2})(b_{0}Y_{l3, m3} + b_{1}Y_{l4, m4})
*/
double Integrate2(CVector &a, CVector &b, IDenseMatrix &a1c, IDenseMatrix &b1c)
{
   toast::complex temp;
   temp =  a[0]*b[0]*kronD(a1c.Row(0), b1c.Row(0)) + a[0]*b[1]*kronD(a1c.Row(0), b1c.Row(1)) + a[1]*b[0]*kronD(a1c.Row(1), b1c.Row(0)) + a[1]*b[1]*kronD(a1c.Row(1), b1c.Row(1));
   return(temp.re);
}


