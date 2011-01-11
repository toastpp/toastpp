#include <sstream>
#include <unistd.h>
#include <time.h>
#include <fstream.h>
#include <iomanip.h>
#include <stdio.h>
#include <stdlib.h>
#include<stdio.h>
#include <mathlib.h>
#include <felib.h>
#include "rte3D_math.h"
using namespace toast;
double fac[] = {1.000000e+000, 1.000000e+000, 2.000000e+000, 6.000000e+000, 2.400000e+001, 1.200000e+002, 7.200000e+002, 5.040000e+003, 4.032000e+004, 3.628800e+005, 3.628800e+006, 3.991680e+007, 4.790016e+008, 6.227021e+009, 8.717829e+010, 1.307674e+012, 2.092279e+013, 3.556874e+014, 6.402374e+015, 1.216451e+017, 2.432902e+018, 5.109094e+019, 1.124001e+021, 2.585202e+022, 6.204484e+023, 1.551121e+025, 4.032915e+026, 1.088887e+028, 3.048883e+029, 8.841762e+030, 2.652529e+032, 8.222839e+033, 2.631308e+035, 8.683318e+036, 2.952328e+038, 1.033315e+040, 3.719933e+041, 1.376375e+043, 5.230226e+044, 2.039788e+046, 8.159153e+047, 3.345253e+049, 1.405006e+051, 6.041526e+052, 2.658272e+054, 1.196222e+056, 5.502622e+057, 2.586232e+059, 1.241392e+061, 6.082819e+062, 3.041409e+064, 1.551119e+066, 8.065818e+067, 4.274883e+069, 2.308437e+071, 1.269640e+073, 7.109986e+074, 4.052692e+076, 2.350561e+078, 1.386831e+080, 8.320987e+081, 5.075802e+083, 3.146997e+085, 1.982608e+087, 1.268869e+089, 8.247651e+090, 5.443449e+092, 3.647111e+094, 2.480036e+096, 1.711225e+098, 1.197857e+100, 8.504786e+101, 6.123446e+103, 4.470115e+105, 3.307885e+107, 2.480914e+109, 1.885495e+111, 1.451831e+113, 1.132428e+115, 8.946182e+116, 7.156946e+118, 5.797126e+120, 4.753643e+122, 3.945524e+124, 3.314240e+126, 2.817104e+128, 2.422710e+130, 2.107757e+132, 1.854826e+134, 1.650796e+136, 1.485716e+138, 1.352002e+140, 1.243841e+142, 1.156773e+144, 1.087366e+146, 1.032998e+148, 9.916779e+149, 9.619276e+151, 9.426890e+153, 9.332622e+155, 9.332622e+157};

double dfactorial[] = {1.000000e+000, 1.000000e+000, 2.000000e+000, 3.000000e+000, 8.000000e+000, 1.500000e+001, 4.800000e+001, 1.050000e+002, 3.840000e+002, 9.450000e+002, 3.840000e+003, 1.039500e+004, 4.608000e+004, 1.351350e+005, 6.451200e+005, 2.027025e+006, 1.032192e+007, 3.445943e+007, 1.857946e+008, 6.547291e+008, 3.715891e+009, 1.374931e+010, 8.174961e+010, 3.162341e+011, 1.961991e+012, 7.905854e+012, 5.101175e+013, 2.134580e+014, 1.428329e+015, 6.190283e+015, 4.284987e+016, 1.918988e+017, 1.371196e+018, 6.332660e+018, 4.662066e+019, 2.216431e+020, 1.678344e+021, 8.200795e+021, 6.377707e+022, 3.198310e+023, 2.551083e+024, 1.311307e+025, 1.071455e+026, 5.638620e+026, 4.714401e+027, 2.537379e+028, 2.168624e+029, 1.192568e+030, 1.040940e+031, 5.843584e+031, 5.204698e+032, 2.980228e+033, 2.706443e+034, 1.579521e+035, 1.461479e+036, 8.687364e+036, 8.184284e+037, 4.951798e+038, 4.746885e+039, 2.921561e+040, 2.848131e+041, 1.782152e+042, 1.765841e+043, 1.122756e+044, 1.130138e+045, 7.297912e+045, 7.458913e+046, 4.889601e+047, 5.072061e+048, 3.373825e+049, 3.550443e+050, 2.395416e+051, 2.556319e+052, 1.748653e+053, 1.891676e+054, 1.311490e+055, 1.437674e+056, 1.009847e+057, 1.121385e+058, 7.977794e+058, 8.971083e+059, 6.462013e+060, 7.356288e+061, 5.363471e+062, 6.179282e+063, 4.558950e+064, 5.314183e+065, 3.966287e+066, 4.676481e+067, 3.529995e+068, 4.208833e+069, 3.212296e+070, 3.872126e+071, 2.987435e+072, 3.639799e+073, 2.838063e+074, 3.494207e+075, 2.752921e+076, 3.424322e+077, 2.725392e+078, 3.424322e+079};


double factorial(int n){
     xASSERT(n >= 0, Factorial is defined only for non-negative integers);
     if(n <= 100) return(fac[n]);
     else return(n * factorial(n-1));
}

double doublefactorial(int n){
	if(n >= 0 & n<=100) return(dfactorial[n]);
	else if(n == -1) return(1);
	else return(n*doublefactorial(n-2));
}

toast::complex cpowi(toast::complex &c, const int m){
        xASSERT(m >= 0, cpowi is defined only for non-negative exponents);
	toast::complex ans(1, 0);
	for(int i = 0; i < m; i++) ans = ans*c;

	return(ans);
}

int sign(int m){
	if(abs(m) % 2 == 0) return(1);
	else return(-1);
}	

int signum(int m){
	if(m == 0)
		return(0);
	else if(m < 0)
		return(1);
	else
		return(2);
}

complex wigCoeff(const int l1, const int l2, const int l3, const int m1, const int m2, const int m3)
{
	/* computes /int Y_{l1, m1} Y_{l2, m2} Y^{*}_{l3, m3}
		inputs:
			l1, l2, l3, m1, m2, m3
	 */
	
	//if(l1>=0 && l2 >=0 && l3 >=0 && l3 >= abs(l1-l2) && l3 <= (l1+l2) && m1 + m2 + m3 == 0) 
	if(m1+m2+m3 == 0)
	{ 
	  double coeff1=0;
	 
	  if(l1+l2-l3 >=0 && l1-l2+l3>=0 && -l1+l2+l3 >=0 )
	  	coeff1 = sign(l1-l2-m3)*sqrt((factorial(l1+l2-l3)*factorial(l1-l2+l3)*factorial(-l1+l2+l3))/(double)factorial(l1+l2+l3+1));
	
	  double coeff2= sqrt(factorial(l1+m1)*factorial(l1-m1)*factorial(l2+m2)*factorial(l2-m2)*factorial(l3+m3)*factorial(l3-m3));
	  double coeff3=0;
          int mint = max(l2-l3-m1, l1-l3+m2);
	  mint = max(mint, 0);
	  int maxt = min(l1+l2-l3, l1-m1);
	  maxt = min(maxt,  l2+m2);
	  for(int t = mint; t <= maxt;t++)
	  { 
		if(t<0 || l3-l2+t+m1 < 0 || l3-l1+t-m2 < 0 || l1+l2-l3-t <0 || l1-t-m1<0 || l2-t+m2<0 ) continue;
		coeff3 += sign(t)/(double)(factorial(t)*factorial(l3-l2+t+m1)*factorial(l3-l1+t-m2)*factorial(l1+l2-l3-t)*factorial(l1-t-m1)*factorial(l2-t+m2));
	  }
	  return(complex(coeff1*coeff2*coeff3, 0)); 
        }
	else return(complex(0, 0));
	
}
/* evaluates \int_{S^{2}} Y_{\ell_{1}, m_{1}} Y_{\ell_{2}, m_{2}} where both Y's are complex spherical harmonics
*/
complex kronD(const IVector &a, const IVector &b)
{
        int l1 = a[0], l2 = b[0], m1 = a[1], m2 = b[1];
	if(l1 >= 0 && l2>=0 && abs(m1) <= l1 && abs(m2) <= l2 && l1==l2 && m1 == -1*m2) return(complex(sign(m2), 0.0));
	else return(complex(0.0, 0.0));
}

/*double plm(const int l, const int m, const double x)
{	
	xASSERT(l>=0 && m >= 0 && l >= m , plm is defined only for non-negative l and m and l>=m );
	if(l==m) return(sign(l)*doublefactorial(2*l-1)*pow(1-x*x, l/2.0));
	if (m == l-1) return(x*(2*l-1)*sign(l-1)*doublefactorial(2*l-3)*pow(1-x*x, (l-1)/2.0));
	return((x*(2*l-1)*plm(l-1, m, x) - (l+m-1)*plm(l-2, m, x))/(double)(l-m));
}

RDenseMatrix LegendreTable(const int l, const int numpts, const RDenseMatrix &pt){
	RDenseMatrix LT;
	LT.New(2*l+1, numpts);
	double x;
	int j, k;

	k=2*l;
	for(int m = l; m >= 0; m--)
	{
	  for(int i=0; i<numpts; i++)
	  {
		double x = pt.Get(i, 2);
		LT(k, i) = plm(l, m, x);
	  }
	  k--; 
	}
	j=0; k = 2*l;
	for(int m = -1*l; m < 0; m++)
	{
		for(int i=0; i<numpts; i++)
			LT(j, i) = LT.Get(k, i)*sign(m)*factorial(l-abs(m))/(double)factorial(l+abs(m));
		j++; k--; 
	}
	return(LT);
}
*/
/*RVector plm(const int l, const int m, const RVector& vec)
{
	xASSERT(l>=0 && m >= 0 && l >= m , plm is defined only for non-negative l and m and l>=m );
	if(l==m){
		RVector ans(vec.Dim()), ons(vec.Dim(), 1.0);
		//for(int i=0; i < vec.Dim(); i++)
			ans =pow(ons- vec*vec, l/2.0)*sign(l)*doublefactorial(2*l-1);
		return(ans);
	} //return(sign(l)*doublefactorial(2*l-1)*pow(1-x*x, l/2.0));
	if (m == l-1){
		RVector ans(vec.Dim()),  ons(vec.Dim(), 1.0);
		for(int i=0; i < vec.Dim(); i++)
			ans = vec*pow(ons-vec*vec, (l-1)/2.0)*(2*l-1)*sign(l-1)*doublefactorial(2*l-3);
		return(ans);
	} //return(x*(2*l-1)*sign(l-1)*doublefactorial(2*l-3)*pow(1-x*x, (l-1)/2.0));
	
	return((vec*(2*l-1)*plm(l-1, m, vec) - plm(l-2, m, vec)*(l+m-1))/(double)(l-m));

}
RDenseMatrix LegendreTable(const int l, const int numpts, const RDenseMatrix &pt){
	RDenseMatrix LT;
	LT.New(2*l+1, numpts);
	double x;
	int j, k;

	k=2*l;
	for(int m = l; m >= 0; m--)
	{
	  //for(int i=0; i<numpts; i++)
	  ///{
		RVector x = pt.Col(2);
		LT.SetRow(k, plm(l, m, x));
	  //}
	  k--; 
	}
	j=0; k = 2*l;
	for(int m = -1*l; m < 0; m++)
	{
		//for(int i=0; i<numpts; i++)
			LT.SetRow(j, LT.Row(k)*sign(m)*factorial(l-abs(m))/(double)factorial(l+abs(m)));
		j++; k--; 
	}
	return(LT);
}
*/
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
	//return(LT);
}

void sphericalHarmonics(const int order, const int numpts, const RDenseMatrix& pt, RDenseMatrix* &sphHarm){
	int m, sgm;
	double dtmp1, dtmp2, normfac;
	CDenseMatrix *cplxSphHarm;
	cplxSphHarm = new CDenseMatrix[order+1];
	//sphHarm = new sphHarm[order+1];
	for(int l=0; l<=order; l++)
	{
		cplxSphHarm[l].New(2*l+1, numpts);
	//	sphHarm[l].New(2*l+1, numpts);
	}
	complex ctmp;
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
							cplxSphHarm[l](k, i) = ctmp*complex(dtmp2, 0);
						}
						break;
	       				case 2:
						if(normfac!=0){
							ctmp.re = pt.Get(i, 0)/normfac; ctmp.im = pt.Get(i, 1)/normfac;
							ctmp = cpowi(ctmp, abs(m));
							cplxSphHarm[l](k, i) = ctmp*complex(dtmp2, 0);
						}
						break;
				}
	    		}	
		k++;	
		}
	}
	complex temp; 
        int indm;
	for(int l=0; l<=order; l++)
	{
		for(m=-1*l; m<= l; m++){
			indm = l + m;
			if(m>0)
			{
				for(int i=0; i<numpts; i++)
				{
					temp = (cplxSphHarm[l](m+l, i) + complex(sign(m), 0)*cplxSphHarm[l](l-m, i))/complex(sqrt(2), 0);
					sphHarm[l](indm, i) = temp.re; 
				}
			}
			else if(m<0)
			{
				for(int i=0; i<numpts; i++)
				{
					temp = (cplxSphHarm[l](l+abs(m), i) - complex(sign(m), 0)*cplxSphHarm[l](l-abs(m), i))/complex(0, sqrt(2));
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
	//cout<<"nRows: "<<sphHarm.nRows()<<" nCols:" << sphHarm.nCols()<<endl;
	//return(sphHarm);	
}
/*RDenseMatrix sphericalHarmonics(const int order, const int numpts, const RDenseMatrix& pt){
	int m, sgm;
	double dtmp1, dtmp2, normfac;
	RDenseMatrix sphHarm;
	CDenseMatrix cplxSphHarm;
	complex ctmp;
        RDenseMatrix LT;
	LT = LegendreTable(order, numpts, pt);
	sphHarm.New(LT.nRows(), LT.nCols());
	cplxSphHarm.New(LT.nRows(), LT.nCols());

	int k = 0;
	for(m = -1*order; m <= order; m++){
	        sgm = signum(m);
		dtmp1 =  sqrt((2*order+1)*factorial(order - m)/(double)(4*M_PI*factorial(order+m)));
		for(int i=0; i<numpts; i++){
			dtmp2 = dtmp1*LT.Get(k, i);
			normfac = sqrt(1-pt.Get(i, 2)*pt.Get(i, 2));
			switch(sgm){
				case 0: 
					cplxSphHarm(k, i).re = dtmp2; cplxSphHarm(k, i).im = 0.0;
					break;
				case 1:
					if(normfac!=0){
						ctmp.re = pt.Get(i, 0)/normfac; ctmp.im = -1*pt.Get(i, 1)/normfac;
						ctmp = cpowi(ctmp, abs(m));
						cplxSphHarm(k, i) = ctmp*complex(dtmp2, 0);
					}
					break;
	       			case 2:
					if(normfac!=0){
						ctmp.re = pt.Get(i, 0)/normfac; ctmp.im = pt.Get(i, 1)/normfac;
						ctmp = cpowi(ctmp, abs(m));
						cplxSphHarm(k, i) = ctmp*complex(dtmp2, 0);
					}
					break;
			}
	    	}	
		k++;	
	}
	complex temp; 
        int indm;
	for(m=-1*order; m<= order; m++){
		indm = order + m;
		if(m>0)
		{
			for(int i=0; i<numpts; i++)
			{
				temp = (cplxSphHarm(m+order, i) + complex(sign(m), 0)*cplxSphHarm(order-m, i))/complex(sqrt(2), 0);
				sphHarm(indm, i) = temp.re; 
			}
		}
		else if(m<0)
		{
			for(int i=0; i<numpts; i++)
			{
				temp = (cplxSphHarm(order-m, i) - complex(sign(m), 0)*cplxSphHarm(order+m, i))/complex(0, sqrt(2));
				sphHarm(indm, i) = temp.re;	 
			}
		}
		else
		{
			for(int i=0; i<numpts; i++)
			{
				temp = cplxSphHarm(indm, i);
				sphHarm(indm, i) = temp.re;
			}
		}
	}
	//cout<<"nRows: "<<sphHarm.nRows()<<" nCols:" << sphHarm.nCols()<<endl;
	return(sphHarm);	
}
*/
/*CDenseMatrix sphericalHarmonicsConj(const int order, const int numpts, const RDenseMatrix& pt)
{
   CDenseMatrix Ylm;
   Ylm = sphericalHarmonics(order, numpts, pt);
   CDenseMatrix Ystarlm(Ylm.nRows(), Ylm.nCols());
   int j=0, k = 2*order;
   for(int m = -1*order; m <= order; m++){
		Ystarlm.SetRow(j, complex(sign(m),0)*Ylm.Row(k)); j++; k--;
   }
          
   return(Ystarlm);	
}*/

/*CDenseMatrix sphericalHarmonicsConj2(const int order, const int numpts, const RDenseMatrix& pt, const CDenseMatrix& Ylm)
{
   CDenseMatrix Ystarlm(Ylm.nRows(), Ylm.nCols());
   int j=0, k = 2*order;
   for(int m = -1*order; m <= order; m++){
		Ystarlm.SetRow(j, complex(sign(m),0)*Ylm.Row(k)); j++; k--;
   }
          
   return(Ystarlm);	
}*/

CCompRowMatrix kronsd(const CCompRowMatrix &A, const CDenseMatrix &B)
{
    int ia, ib, ka, kb, ja, jb, i, idx;
    int va_i, vb_i, v_i;
    int na = A.nRows(), ma = A.nCols(), va = A.nVal();
    int nb = B.nRows(), mb = B.nCols(), vb = nb*mb;
    int n = na*nb, m = ma*mb, v = va*vb;
    complex a_ij, b_ij;
    const complex *Aval = A.ValPtr();
    int *rowptr = new int[n+1];
    int *colidx = new int[v];
    complex *val = new complex[v];

    rowptr[0] = 0;
    i = idx = 0;
    for (ia = 0; ia < na; ia++) {
	va_i = A.rowptr[ia+1] - A.rowptr[ia]; // nonzeros in row ia of A
	for (ib = 0; ib < nb; ib++) {
	    vb_i = mb; // nonzeros in row ib of B
	    v_i = va_i * vb_i;
	    rowptr[i+1] = rowptr[i] + v_i;
	    for (ka = A.rowptr[ia]; ka < A.rowptr[ia+1]; ka++) {
		ja = A.colidx[ka];
		a_ij = Aval[ka];
		for (kb = 0; kb < mb; kb++) {
		    jb = kb;
		    b_ij = B.Get(ib, kb);
		    colidx[idx] = ja*mb + jb;
		    val[idx] = a_ij * b_ij;
		    idx++;
		}
	    }
	    i++;
	}
    }
    CCompRowMatrix C (n, m, rowptr, colidx, val);
    delete []rowptr;
    delete []colidx;
    delete []val;
    return C;
}

void kronplus(const int spatrow, const int spatcol, const IVector& node_angN, const IVector& offset, const double a_ij, const RCompRowMatrix &B, RCompRowMatrix& C)
{
    double *Cval = C.ValPtr();
    int row_offset = 0, col_offset = 0;
    const int *rowptr, *colidx;
    B.GetSparseStructure(&rowptr, &colidx);
    const double *bval = B.ValPtr(); 
    int nr = B.nRows();
    int nc = B.nCols();
 
    row_offset = offset[spatrow];
	
    col_offset = offset[spatcol];
  
    /* for(int ib = 0; ib < node_angN[spatrow]; ib++)
     {
	 for(int jb=0; jb < node_angN[spatcol]; jb++)
	 {
		//cout<<"i: "<<row_offset+ib<<" j: "<<col_offset+jb<<endl;
		C(row_offset+ib , col_offset+jb) = C.Get(row_offset+ib , col_offset+jb) + complex(a_ij, 0)*B.Get(ib, jb); 	
	}
     }*/
    //cout<<"Inside kronplus"<<endl;
   // cout<<"spatrowcols: "<<spatrow<<"  "<<spatcol<<endl;
    for(int ib=0; ib < nr; ib++)
    {
	for(int j=rowptr[ib]; j < rowptr[ib+1]; j++)
	{
		int jb = colidx[j];
		if(ib<node_angN[spatrow] && jb <node_angN[spatcol]){
	        //cout<<"ib: "<<ib<<" jb: "<<jb<<endl;
		C(row_offset+ib , col_offset+jb) = C.Get(row_offset+ib , col_offset+jb) + a_ij*bval[j]; 
		}
	}
    }
    
}

RCompRowMatrix shrink(const RDenseMatrix &dnsmat)
{
    int i, j;
    int *rowptr, *colidx;
    double *val;
    int m =  dnsmat.nRows(), n= dnsmat.nCols();

    rowptr = new int[m+1];
    rowptr[0] = 0;
    for (i = 0; i < m; i++) {
	int nz=0;
	for (j = 0; j < n; j++)
		if(!(fabs(dnsmat.Get(i, j)) < 1e-15))
			nz++;
	rowptr[i+1] = rowptr[i] + nz;
    }
    
    colidx = new int[rowptr[m]];
    val = new double[rowptr[m]];
    int k=0;
    for(i=0; i < m; i++){
	for(j=0; j<n; j++){
		if(!(fabs(dnsmat.Get(i, j)) < 1e-15)){
			colidx[k] = j; 
			val[k] = dnsmat.Get(i, j);
			k++;
		}
	}
    }
    RCompRowMatrix C (m, n, rowptr, colidx, val);

    delete []rowptr;
    delete []colidx;
    delete []val;	
    return C;
}

inline CDenseMatrix sdmatmult(const CCompRowMatrix &A, const CDenseMatrix &B)
{
    dASSERT(A.nCols() == B.nRows(), Invalid sizes of matrices);
    int i, j, k, m, ra, ra1, ra2;
    int nr = A.nRows();
    int nc = B.nCols();
    CDenseMatrix C(nr, nc);
    const complex *val = A.ValPtr();
    complex *valc = C.data_buffer();
    for (i = 0; i < nr; i++) {
	ra1 = A.rowptr[i];
	ra2 = A.rowptr[i+1];
	for (ra = ra1; ra < ra2; ra++) {
	    k = A.colidx[ra];
	    for(j = 0; j < nc; j++){
		valc[i*nc+j] += val[ra]*B.Get(k, j); 
		}
     }
	}
  return C;	
}

void BIntUnitSphere(const int size1, const int size2, const int sphOrder1, const int sphOrder2,  const RDenseMatrix& pts, const RVector& wts, const RVector& bnormal, RDenseMatrix* &Ylm, RCompRowMatrix& bintplus, RCompRowMatrix& bintminus)
{
	/*Computes the integral on the sphere of the form 
	 * \int_{S^{n-1}} (s.n)_{plusminus} \psi_{i} \psi_{j}
	 * Inputs:
	 * 	order -> Order of spherical harmonics to be used
	 * 	int_order -> Integration order which determines the number of Gauss Points
	 * 	bnormal -> unit normal to the boundary face/element 
	 * 	isplus -> '+' side or '-' side of the s.n 
	 * 			TRUE denotes '+' , FALSE denotes '-'
	 * Output:
	 * 	bint -> boundary integral*/

   int numpts = pts.nRows(); 
    RDenseMatrix dnssdotnplus(1, numpts), dnssdotnminus(1, numpts);
//     CVector sdotnplus(numpts);
     for(int i=0; i < numpts; i++){
	     double tmp = bnormal[0]*pts.Get(i, 0) + bnormal[1]*pts.Get(i, 1) + bnormal[2]*pts.Get(i, 2);
	     if(tmp>0){
	     	dnssdotnplus(0, i) = tmp; 
		//sdotnplus[i] = tmp;
	     }
	     else if(tmp<0){
		dnssdotnminus(0, i) = -tmp;
		}
     }
    RCompRowMatrix sdotnplus = shrink(dnssdotnplus);
    RCompRowMatrix sdotnminus = shrink(dnssdotnminus);
    const int *rowptr1, *colidx1, *rowptr2, *colidx2;
    sdotnplus.GetSparseStructure(&rowptr1, &colidx1);
    sdotnminus.GetSparseStructure(&rowptr2, &colidx2);
       
     int is, js, indl1, indm1, indl2, indm2;
     RDenseMatrix dnsbintplus(size1, size2), dnsbintminus(size1, size2);
     for(int l1 = 0; l1 <= sphOrder1; l1++){
	indl1 = getPos(l1, -1*l1);
	for(int m1 = -1*l1; m1 <= l1; m1++){
		indm1 = l1 + m1;
		is = indl1 + indm1;
	        for(int l2 = 0; l2 <= sphOrder2; l2++){
			indl2 = getPos(l2, -l2);
			for(int m2 = -l2; m2 <= l2; m2++){
				indm2 = l2 + m2;
				js = indl2 + indm2;
				for(int i=0; i<sdotnplus.nRows(); i++)
				{
					for(int j=rowptr1[i]; j < rowptr1[i+1]; j++)
						dnsbintplus(is, js) += sdotnplus(i, colidx1[j])*Ylm[l1](indm1, colidx1[j])*Ylm[l2](indm2, colidx1[j])*wts[colidx1[j]]*4*M_PI;
				}
				for(int i=0; i<sdotnminus.nRows(); i++)
				{
					for(int j=rowptr2[i]; j < rowptr2[i+1]; j++)
						dnsbintminus(is, js) += sdotnminus(i, colidx2[j])*Ylm[l1](indm1, colidx2[j])*Ylm[l2](indm2, colidx2[j])*wts[colidx2[j]]*4*M_PI;
				}
				//if (((dnsbintplus.Get(is, js)-dnsbintminus.Get(is, js)) -(Aintsc.Get(is, js)*bnormal[0] + Aintss.Get(is, js)*bnormal[1] + Aintc.Get(is, js)*bnormal[2])).re > 1e-12 && ((dnsbintplus.Get(is, js)-dnsbintminus.Get(is, js)) -(Aintsc.Get(is, js)*bnormal[0] + Aintss.Get(is, js)*bnormal[1] + Aintc.Get(is, js)*bnormal[2])).im > 1e-12)
			     //cout<<l1<<" "<<m1<<" "<<l2<<" "<<m2<<" "<<dnsbintplus.Get(is, js)<<endl;//<<" "<<dnsbintminus.Get(is, js)+Aintss.Get(is, js)<<endl;					
			}
		}
	
	     }
	}
			//CDenseMatrix dnsbintplus = sdmatmult(dotspdns(sdotnplus, Yl1m1), Yl2m2_modi);
	
         bintplus = shrink(dnsbintplus);
	 bintminus = shrink(dnsbintminus);	

} 
/*	int angN = bintplus.nRows();
	CDenseMatrix dnsAint(angN, angN), dnsAintsc(angN, angN), dnsAintss(angN, angN), dnsAintc(angN, angN);
	
	complex a1, b1, c1, d1, e1, f1, a2, b2, c2, d2, e2, f2;
	complex coeff_sc, coeff_ss, phi_coeff;
	coeff_sc = complex(0.5, 0); coeff_ss = complex(0, -0.5);

        int a1c[2], b1c[2], c1c[2], d1c[2], e1c[2], f1c[2], p1c[2];
	int a2c[2], b2c[2], c2c[2], d2c[2], e2c[2], f2c[2], p2c[2];

	int is, js, indl1, indl2;

	for(int l1 = 0; l1 <= sphOrder; l1++){
		indl1 = getPos(l1, -1*l1);
		for(int m1 = -1*l1; m1 <= l1; m1++){
			a1 = sqrt(complex((l1 + m1)*(l1+m1-1)/(double)((2*l1+1)*(2*l1-1)), 0));
			b1 = sqrt(complex((l1 - m1 +1)*(l1-m1+2)/(double)((2*l1+1)*(2*l1+3)), 0));
			c1 = sqrt(complex((l1 - m1)*(l1-m1-1)/(double)((2*l1+1)*(2*l1-1)), 0)); 
			d1 = sqrt(complex((l1 + m1 +1)*(l1+m1+2)/(double)((2*l1+1)*(2*l1+3)), 0));
			e1 = sqrt(complex((l1 + m1)*(l1-m1)/(double)((2*l1+1)*(2*l1-1)), 0)); 
			f1 = sqrt(complex((l1 + m1 +1)*(l1- m1+1)/(double)((2*l1+1)*(2*l1+3)), 0));
			a1c[0] = l1-1; a1c[1] = m1-1;
			b1c[0] = l1+1; b1c[1] = m1-1;
			c1c[0] = l1-1; c1c[1] = m1+1;
			d1c[0] = l1+1; d1c[1] = m1+1;
			e1c[0] = l1-1; e1c[1] = m1;
			f1c[0] = l1+1; f1c[1] = m1;
                        p1c[0] = l1; p1c[1] = m1;
			
			is = indl1 + l1 + m1;    	
			for(int l2 = 0; l2 <= sphOrder; l2++){
				indl2 = getPos(l2, -1*l2);
				phi_coeff = PHI_COEFF(l2);
				for(int m2 = -1*l2; m2 <= l2; m2++){
					a2 = sqrt(complex((l2 - m2)*(l2 - m2 - 1)/(double)((2*l2+1)*(2*l2-1)), 0));
					b2 = sqrt(complex((l2 + m2 +1)*(l2 + m2 + 2)/(double)((2*l2+1)*(2*l2+3)), 0));
					c2 = sqrt(complex((l2 + m2)*(l2 + m2 - 1)/(double)((2*l2+1)*(2*l2-1)), 0)); 
					d2 = sqrt(complex((l2 - m2 + 1)*(l2 - m2 + 2)/(double)((2*l2+1)*(2*l2+3)), 0));
					e2 = sqrt(complex((l2 + m2)*(l2-m2)/(double)((2*l2+1)*(2*l2-1)), 0)); 
					f2 = sqrt(complex((l2 + m2 +1)*(l2- m2+1)/(double)((2*l2+1)*(2*l2+3)), 0));
					a2c[0] = l2-1; a2c[1] = m2+1;
					b2c[0] = l2+1; b2c[1] = m2+1;
					c2c[0] = l2-1; c2c[1] = m2-1;
					d2c[0] = l2+1; d2c[1] = m2-1;
					e2c[0] = l2-1; e2c[1] = m2;
					f2c[0] = l2+1; f2c[1] = m2;
					p2c[0] = l2; p2c[1] = m2;

					js = indl2 + l2 + m2;
				     
					dnsAint(is, js) = phi_coeff *kronD(p1c, p2c);
					
				     dnsAintsc(is, js) = (-a1*kronD(a1c, p2c) + c1*kronD(c1c, p2c) + b1*kronD(b1c, p2c) - d1*kronD(d1c, p2c));
				     dnsAintsc(is, js) *= coeff_sc*phi_coeff ;
				     
				     dnsAintss(is, js) = (-a1*kronD(a1c, p2c) + b1*kronD(b1c, p2c));
				     dnsAintss(is, js) += (-c1*kronD(c1c, p2c) + d1*kronD(d1c, p2c));
				     dnsAintss(is, js) *= coeff_ss*phi_coeff ;

				     dnsAintc(is, js) = phi_coeff *(e1*kronD(e1c, p2c) + f1*kronD(f1c,p2c));
				}
			}
		}
	     }

	    bintplus = dnsAintsc*complex(bnormal[0], 0) + dnsAintss*complex(bnormal[1], 0) + dnsAintc*complex(bnormal[2], 0);
	    bintminus = dnsAintsc*complex(bnormal[0], 0) + dnsAintss*complex(bnormal[1], 0) + dnsAintc*complex(bnormal[2], 0);

}
*/

/*void BIntUnitSphere(const int sphOrder1, const int size, const RDenseMatrix& ptsPlus, const CVector& wtsPlus, const RVector& bnormal, CCompRowMatrix& Aintsc, CCompRowMatrix &Aintss, CCompRowMatrix &Aintc, CCompRowMatrix& intSdotnPlusHemisphere, CCompRowMatrix& intSdotnMinusHemisphere)
{


	double alpha, beta;
	double z = (fabs(bnormal[2]) > 1e-06) ? bnormal[2] : 0;
	alpha = acos(z);
	if(sin(alpha) == 0) beta =0; //beta doesn't matter in this case
	else
	{
		double beta1 = acos(bnormal[0]/sin(alpha));
        	double beta2 = asin(bnormal[1]/sin(alpha));
		double cos1 = (fabs(cos(beta1)) > 1e-6) ? cos(beta1) : 0, sin2 = (fabs(sin(beta2)) > 1e-6)? sin(beta2) : 0;

		if (cos1 >=0 && sin2>=0) beta = beta1; // 0 <= beta <= Pi/2
		else if(cos1<=0 && sin2>=0) beta = beta1; //Pi/2 <= beta <= Pi
		else if (cos1 <=0 && sin2<=0) beta = acos(fabs(bnormal[0]/sin(alpha))) + M_PI; //Pi <= beta <=3Pi/2
		else if (cos1 >=0 && sin2 <=0) beta = beta2; //3*Pi/2 <= beta <= Pi
	}

	if(fabs(bnormal[0]- sin(alpha)*cos(beta)) > 1e-12 || fabs(bnormal[1] - sin(alpha)*sin(beta)) > 1e-12 || fabs(bnormal[2] - cos(alpha)) > 1e-12)
		cout<<"Warning: Failure bnormal"<<bnormal<<"******** ("<<sin(alpha)*cos(beta)<<", "<<sin(alpha)*sin(beta)<<", "<<cos(alpha)<<")"<<endl;
	
	double  cosalpha = cos(alpha), sinalpha = sin(alpha), cosbeta = cos(beta), sinbeta = sin(beta);
	int numPts = ptsPlus.nRows();	
	RDenseMatrix Amat(3, 3), Bmat(3, 3), rotPts(numPts, ptsPlus.nCols());
	
	Amat(0, 0) = cosalpha; Amat(0, 1) = 0.0; Amat(0, 2) = sinalpha;
	Amat(1, 0) = 0; Amat(1, 1) = 1; Amat(1, 2) = 0;
	Amat(2, 0) = -sinalpha; Amat(2, 1) = 0; Amat(2, 2) = cosalpha;

	Bmat(0, 0) = cosbeta; Bmat(0, 1) = -sinbeta; Bmat(0, 2) = 0.0;
	Bmat(1, 0) = sinbeta; Bmat(1, 1) = cosbeta; Bmat(1, 2) = 0.0;
	Bmat(2, 0) = 0.0; Bmat(2, 1) = 0.0; Bmat(2, 2) = 1.0;

	for(int i=0; i < ptsPlus.nRows(); i++)
	  rotPts.SetRow(i, Bmat*Amat*ptsPlus.Row(i));

	
	CDenseMatrix dnsintSdotnSphere(size, size), dnsintSdotnPlusHemisphere(size, size), dnsintSdotnMinusHemisphere(size, size);

	int is, js, indl1, indl2, indm1, indm2;
	CDenseMatrix *Ystarl1m1, *Yl2m2;
	Ystarl1m1 = new CDenseMatrix[sphOrder1+1];
	Yl2m2 = new CDenseMatrix[sphOrder1+1];

	for(int l2=0; l2<=sphOrder1; l2++)
		Yl2m2[l2] = sphericalHarmonics(l2, numPts, rotPts);
	for(int l1=0; l1 <= sphOrder1; l1++)
		Ystarl1m1[l1] = sphericalHarmonicsConj2(l1, numPts, rotPts, Yl2m2[l1]);

	for(int l1 = 0; l1 <= sphOrder1; l1++){
		indl1 = getPos(l1, -1*l1);
		for(int m1 = -1*l1; m1 <= l1; m1++){
			indm1 = l1+m1;
			is = indl1 + indm1;    	
			for(int l2 = 0; l2 <= sphOrder1; l2++){
				indl2 = getPos(l2, -1*l2);

				for(int m2 = -1*l2; m2 <= l2; m2++){
					indm2 = l2 + m2;
					js = indl2 + l2 + m2;
				     
					
				     dnsintSdotnSphere(is, js) = Aintsc.Get(is, js)*bnormal[0] + Aintss.Get(is, js)*bnormal[1] + Aintc.Get(is, js)*bnormal[2];
				     	for(int i=0; i<numPts; i++)
					{
						dnsintSdotnPlusHemisphere(is, js) += (Ystarl1m1[l1](indm1, i)*Yl2m2[l2](indm2, i)*(rotPts.Get(i, 0)*bnormal[0] + rotPts.Get(i, 1)*bnormal[1] + rotPts.Get(i, 2)*bnormal[2]))*wtsPlus[i];
						if((rotPts.Get(i, 0)*sinalpha*cosbeta + rotPts.Get(i, 1)*sinalpha*sinbeta + rotPts.Get(i, 2)*cosalpha)<-1e-12) 
						{
							cout<<"warning: Sdotn negative!!!"<<endl;
							cout<<"bnormal: "<<bnormal<<" computed vector: ["<<sin(alpha)*cos(beta)<<" "<<sin(alpha)*sin(beta)<<" "<<cos(alpha)<<endl;
							cout<<"Magnitude SDOTN: "<<(rotPts.Get(i, 0)*sinalpha*cosbeta + rotPts.Get(i, 1)*sinalpha*sinbeta + rotPts.Get(i, 2)*cosalpha)<<endl;
						}
					}
					//cout<<l1<<" "<<m1<<" "<<l2<<" "<<m2<<" "<<dnsintSdotnPlusHemisphere(is, js)<<endl;
				}
			}
		}
	}

	dnsintSdotnMinusHemisphere = dnsintSdotnPlusHemisphere- dnsintSdotnSphere;

	intSdotnPlusHemisphere = shrink(dnsintSdotnPlusHemisphere);
	intSdotnMinusHemisphere = shrink(dnsintSdotnMinusHemisphere);
	delete []Ystarl1m1;
	delete []Yl2m2;

}
*/

/*void BIntUnitSphere(const int sphOrder1, const RDenseMatrix& pts, const CVector& wtsPlus, const RVector& bnormal, CDenseMatrix& intSdotnPlusHemiSphere, CDenseMatrix& intSdotnMinusHemiSphere)
{
	cout<<"Enter BIntUnitSphere ..."<<endl;
	int numPluspts =0;
	for(int i=0; i < pts.nRows(); i++)
		if(pts(i, 0)*bnormal[0] + pts(i, 1)*bnormal[1] + pts(i, 2)*bnormal[2] >= 0)
			numPluspts++;
	//cout<<"Number of points selected: " <<numPluspts<<" wts available: "<<wtsPlus.Dim()<<endl;	
	cout<<"Number of points computed ..."<<endl;
	RDenseMatrix ptsPlus(numPluspts, 3);
	int k=0;
	for(int i=0; i < pts.nRows(); i++)
	{
		if(pts(i, 0)*bnormal[0] + pts(i, 1)*bnormal[1] + pts(i, 2)*bnormal[2] >= 0)
		{
			ptsPlus.SetRow(k, pts.Row(i));
			k++;
		}

	}
	//cout<<"Points collated ...Numpts: "<<numPluspts<<" k: "<<k<<endl;
	//CDenseMatrix Amat(numPluspts, numPluspts);
	//CVector bvec(numPluspts);
	//cout<<"Forming the system matrix ..."<<endl;
	
	//k=0;
	//for(int l=0; l<= 18; l++)
	//{
	//	cout<<l<<endl;
	//	CDenseMatrix Ylm = sphericalHarmonics(l, numPluspts, ptsPlus);
	//	for(int m = -l; m <= l; m++)
	//	{
	//		if (k>=numPluspts) break; 
	//
	//		int indm = l+m;
	//		for(int i=0; i < numPluspts; i++)
	//				Amat(k, i) = Ylm(indm, i);
				
	//		if(l==0 && m==0) bvec[k] = sqrt(M_PI);
	//		else if(l%2 == 1 && m == 0)
	//		{ 
	//		  RVector res(1), inp(1);
	//		  inp[0] = 0;
	//		  res = plm(l, 1, inp);
	//		  bvec[k] = sqrt((2*l+1)*M_PI)*res[0]/((double)l*(l+1));
	//		}
	//		else bvec[k] = 0; 
	//		k++;
			
	//	}
	//}
	//cout<<endl;
	//cout<<"System of equations formed...."<<endl;	
	//CVector wts= inverse(Amat)*bvec;
	//CVector wts(numPluspts), cvec(numPluspts), dvec(numPluspts);
	//cout<<"Solving ..."<<endl;
	//int isnotsingular = QRFactorize(Amat, cvec, dvec);
	//if(!isnotsingular) cout<< "Singular system: "<<isnotsingular<<endl;

	//QRSolve(Amat, cvec, dvec, bvec, wts);
	//cout<<sum(wts)<<endl;
	CDenseMatrix intSdotnSphere(intSdotnPlusHemiSphere.nRows(), intSdotnPlusHemiSphere.nCols());
	complex a1, b1, c1, d1, e1, f1, a2, b2, c2, d2, e2, f2;
	complex coeff_sc, coeff_ss;
	coeff_sc = complex(0.5, 0); coeff_ss = complex(0, -0.5);

        int a1c[2], b1c[2], c1c[2], d1c[2], e1c[2], f1c[2], p1c[2];
	int a2c[2], b2c[2], c2c[2], d2c[2], e2c[2], f2c[2], p2c[2];

	int is, js, indl1, indl2, indm1, indm2;
	CDenseMatrix *Ystarl1m1, *Yl2m2;
	Ystarl1m1 = new CDenseMatrix[sphOrder1+1];
	Yl2m2 = new CDenseMatrix[sphOrder1+1];

	for(int l2=0; l2<=sphOrder1; l2++)
		Yl2m2[l2] = sphericalHarmonics(l2, numPluspts, ptsPlus);
	for(int l1=0; l1 <= sphOrder1; l1++)
		Ystarl1m1[l1] = sphericalHarmonicsConj2(l1, numPluspts, ptsPlus, Yl2m2[l1]);

	for(int l1 = 0; l1 <= sphOrder1; l1++){
		indl1 = getPos(l1, -1*l1);
		for(int m1 = -1*l1; m1 <= l1; m1++){
			a1 = sqrt(complex((l1 + m1)*(l1+m1-1)/(double)((2*l1+1)*(2*l1-1)), 0));
			b1 = sqrt(complex((l1 - m1 +1)*(l1-m1+2)/(double)((2*l1+1)*(2*l1+3)), 0));
			c1 = sqrt(complex((l1 - m1)*(l1-m1-1)/(double)((2*l1+1)*(2*l1-1)), 0)); 
			d1 = sqrt(complex((l1 + m1 +1)*(l1+m1+2)/(double)((2*l1+1)*(2*l1+3)), 0));
			e1 = sqrt(complex((l1 + m1)*(l1-m1)/(double)((2*l1+1)*(2*l1-1)), 0)); 
			f1 = sqrt(complex((l1 + m1 +1)*(l1- m1+1)/(double)((2*l1+1)*(2*l1+3)), 0));
			a1c[0] = l1-1; a1c[1] = m1-1;
			b1c[0] = l1+1; b1c[1] = m1-1;
			c1c[0] = l1-1; c1c[1] = m1+1;
			d1c[0] = l1+1; d1c[1] = m1+1;
			e1c[0] = l1-1; e1c[1] = m1;
			f1c[0] = l1+1; f1c[1] = m1;
                        p1c[0] = l1; p1c[1] = m1;
			
			indm1 = l1+m1;
			is = indl1 + indm1;    	
			for(int l2 = 0; l2 <= sphOrder1; l2++){
				indl2 = getPos(l2, -1*l2);
				complex phicoeff = PHI_COEFF(l2);
				if(phicoeff != complex(1.0, 0))
					cout<<phicoeff<<endl;

				for(int m2 = -1*l2; m2 <= l2; m2++){
					p2c[0] = l2; p2c[1] = m2;

					indm2 = l2 + m2;
					js = indl2 + l2 + m2;
				     
					
				     intSdotnSphere(is, js) = (-a1*kronD(a1c, p2c) + c1*kronD(c1c, p2c) + b1*kronD(b1c, p2c) - d1*kronD(d1c, p2c))*coeff_sc*bnormal[0];
				     
				     intSdotnSphere(is, js) += ((-a1*kronD(a1c, p2c) + b1*kronD(b1c, p2c)) + (-c1*kronD(c1c, p2c) + d1*kronD(d1c, p2c)))*coeff_ss*bnormal[1];

				     intSdotnSphere(is, js) += ((e1*kronD(e1c, p2c) + f1*kronD(f1c,p2c)))*bnormal[2];

				     intSdotnSphere(is, js) *= phicoeff;
				
				     	for(int i=0; i<numPluspts; i++)
					{
						intSdotnPlusHemiSphere(is, js) += (Ystarl1m1[l1](indm1, i)*Yl2m2[l2](indm2, i)*(ptsPlus(i, 0)*bnormal[0] + ptsPlus(i, 1)*bnormal[1] + ptsPlus(i, 2)*bnormal[2]));
					}
						intSdotnPlusHemiSphere(is, js) *= phicoeff;
					
				}
			}
		}
	}

	intSdotnMinusHemiSphere = intSdotnPlusHemiSphere - intSdotnSphere;


}
*/

int getPos(const int l, const int m){
	int indm = 0;
//	for(int i=0; i < l;  i++)
//		indm += (2*i+1);
	return(l*l + (l+m));
}

/*void tabulatePtsWtsForSphere()
{


 double pts17[110][3] = {{0, 0, 1.0}, {0, 1.0, 0}, {1.0, 0, 0}, {0, 0, -1.0}, {0, -1.0, 0}, {-1.0, 0, 0}, 

{.57735026918962576, .57735026918962576, .57735026918962576}, {.57735026918962576, .57735026918962576, -0.57735026918962576}, {.57735026918962576, -0.57735026918962576, .57735026918962576}, {-0.57735026918962576, .57735026918962576, .57735026918962576}, {-0.57735026918962576, -0.57735026918962576, .57735026918962576}, {-0.57735026918962576, .57735026918962576, -0.57735026918962576}, {.57735026918962576, -0.57735026918962576, -0.57735026918962576}, {-0.57735026918962576, -0.57735026918962576, -0.57735026918962576}, 

{0.18511563534473617, 0.18511563534473617, 0.9651240350865941}, {0.18511563534473617, 0.18511563534473617, -0.9651240350865941}, {0.18511563534473617, -0.18511563534473617, 0.9651240350865941}, {-0.18511563534473617, 0.18511563534473617, 0.9651240350865941}, {-0.18511563534473617, -0.18511563534473617, 0.9651240350865941}, {-0.18511563534473617, 0.18511563534473617, -0.9651240350865941}, {0.18511563534473617, -0.18511563534473617, -0.9651240350865941}, {-0.18511563534473617, -0.18511563534473617, -0.9651240350865941}, 

{0.18511563534473617, 0.9651240350865941, 0.18511563534473617}, {0.18511563534473617, 0.9651240350865941, -0.18511563534473617}, {0.18511563534473617, -0.9651240350865941, 0.18511563534473617}, {-0.18511563534473617, 0.9651240350865941, 0.18511563534473617}, {-0.18511563534473617, -0.9651240350865941, 0.18511563534473617}, {-0.18511563534473617, 0.9651240350865941, -0.18511563534473617}, {0.18511563534473617, -0.9651240350865941, -0.18511563534473617}, {-0.18511563534473617, -0.9651240350865941, -0.18511563534473617}, 

{0.9651240350865941, 0.18511563534473617, 0.18511563534473617}, {0.9651240350865941, 0.18511563534473617, -0.18511563534473617}, {0.9651240350865941, -0.18511563534473617, 0.18511563534473617}, {-0.9651240350865941, 0.18511563534473617, 0.18511563534473617}, {-0.9651240350865941, -0.18511563534473617, 0.18511563534473617}, {-0.9651240350865941, 0.18511563534473617, -0.18511563534473617}, {0.9651240350865941, -0.18511563534473617, -0.18511563534473617}, {-0.9651240350865941, -0.18511563534473617, -0.18511563534473617}, 

{0.39568947305594191, 0.39568947305594191, 0.8287699812525922}, {0.39568947305594191, 0.39568947305594191, -0.8287699812525922}, {0.39568947305594191, -0.39568947305594191, 0.8287699812525922}, {-0.39568947305594191, 0.39568947305594191, 0.8287699812525922}, {-0.39568947305594191, -0.39568947305594191, 0.8287699812525922}, {-0.39568947305594191, 0.39568947305594191, -0.8287699812525922}, {0.39568947305594191, -0.39568947305594191, -0.8287699812525922}, {-0.39568947305594191, -0.39568947305594191, -0.8287699812525922}, 

{0.39568947305594191, 0.8287699812525922, 0.39568947305594191}, {0.39568947305594191, 0.8287699812525922, -0.39568947305594191}, {0.39568947305594191, -0.8287699812525922, 0.39568947305594191}, {-0.39568947305594191, 0.8287699812525922, 0.39568947305594191}, {-0.39568947305594191, -0.8287699812525922, 0.39568947305594191}, {-0.39568947305594191, 0.8287699812525922, -0.39568947305594191}, {0.39568947305594191, -0.8287699812525922, -0.39568947305594191}, {-0.39568947305594191, -0.8287699812525922, -0.39568947305594191}, 

{0.8287699812525922, 0.39568947305594191, 0.39568947305594191}, {0.8287699812525922, 0.39568947305594191, -0.39568947305594191}, {0.8287699812525922, -0.39568947305594191, 0.39568947305594191}, {-0.8287699812525922, 0.39568947305594191, 0.39568947305594191}, {-0.8287699812525922, -0.39568947305594191, 0.39568947305594191}, {-0.8287699812525922, 0.39568947305594191, -0.39568947305594191}, {0.8287699812525922, -0.39568947305594191, -0.39568947305594191}, {-0.8287699812525922, -0.39568947305594191, -0.39568947305594191}, 

{0.69042104838229218, 0.69042104838229218, 0.21595729184584883}, {0.69042104838229218, 0.69042104838229218, -0.21595729184584883}, {0.69042104838229218, -0.69042104838229218, 0.21595729184584883}, {-0.69042104838229218, 0.69042104838229218, 0.21595729184584883}, {-0.69042104838229218, -0.69042104838229218, 0.21595729184584883}, {-0.69042104838229218, 0.69042104838229218, -0.21595729184584883}, {0.69042104838229218, -0.69042104838229218, -0.21595729184584883}, {-0.69042104838229218, -0.69042104838229218, -0.21595729184584883},

{0.69042104838229218, 0.21595729184584883, 0.69042104838229218}, {0.69042104838229218, 0.21595729184584883, -0.69042104838229218}, {0.69042104838229218, -0.21595729184584883, 0.69042104838229218}, {-0.69042104838229218, 0.21595729184584883, 0.69042104838229218}, {-0.69042104838229218, -0.21595729184584883, 0.69042104838229218}, {-0.69042104838229218, 0.21595729184584883, -0.69042104838229218}, {0.69042104838229218, -0.21595729184584883, -0.69042104838229218}, {-0.69042104838229218, -0.21595729184584883, -0.69042104838229218},

{0.21595729184584883, 0.69042104838229218, 0.69042104838229218}, {0.21595729184584883, 0.69042104838229218, -0.69042104838229218}, {0.21595729184584883, -0.69042104838229218, 0.69042104838229218}, {-0.21595729184584883, 0.69042104838229218, 0.69042104838229218}, {-0.21595729184584883, -0.69042104838229218, 0.69042104838229218}, {-0.21595729184584883, 0.69042104838229218, -0.69042104838229218}, {0.21595729184584883, -0.69042104838229218, -0.69042104838229218}, {-0.21595729184584883, -0.69042104838229218, -0.69042104838229218},

{0.47836902881215020, 0, 0.8781589106040662}, {0.47836902881215020, 0, -0.8781589106040662}, {-0.47836902881215020, 0, 0.8781589106040662}, {-0.47836902881215020, 0, -0.8781589106040662}, {0.8781589106040662, 0, 0.47836902881215020}, {0.8781589106040662, 0, -0.47836902881215020}, {-0.8781589106040662, 0, 0.47836902881215020}, {-0.8781589106040662, 0, -0.47836902881215020}, 

{0, 0.8781589106040662, 0.47836902881215020}, {0, 0.8781589106040662, -0.47836902881215020}, {0, -0.8781589106040662, 0.47836902881215020}, {0, -0.8781589106040662, -0.47836902881215020}, {0, 0.47836902881215020, 0.8781589106040662}, {0, 0.47836902881215020, -0.8781589106040662}, {0, -0.47836902881215020, 0.8781589106040662}, {0, -0.47836902881215020, -0.8781589106040662}, 

{0.47836902881215020, 0.8781589106040662, 0}, {0.47836902881215020, -0.8781589106040662, 0}, {-0.47836902881215020, 0.8781589106040662, 0}, {-0.47836902881215020, -0.8781589106040662, 0}, {0.8781589106040662, 0.47836902881215020, 0}, {0.8781589106040662, -0.47836902881215020, 0}, {-0.8781589106040662, 0.47836902881215020, 0}, {-0.8781589106040662, -0.47836902881215020, 0}};

 
	
	
	double x, y, z;
	int numpts = 110;

	RDenseMatrix pts(numpts, 3);
	for(int i=0; i < numpts; i++)
		for(int j=0; j< 3; j++)
			pts(i, j) = pts17[i][j];
 	CDenseMatrix A1(numpts, numpts);
	CVector b1(numpts);
	FILE *fid;
	
	int k=0;
	for(int l=0; l<= 25; l++)
	{
		RDenseMatrix Ylm = sphericalHarmonics(l, numpts, pts); 
		for(int m = -l; m <= l; m++)
		{
			int indm = l+m;
			if(k < numpts)
			{
				for(int i=0; i < numpts; i++)
					A1(k, i) = Ylm(indm, i);
				
				if(l==0 && m==0) b1[k] = sqrt(4*M_PI);
				else b1[k] = 0; 
				k++;
				cout<<"k: "<<k<<endl;
			}
		}
	}

	fid = fopen("SphereQuad_sysmat.txt", "w");
	fprintf(fid, "[");
	for(int i=0; i < numpts; i++)
	{
		for(int j=0; j < numpts; j++)
	  		fprintf(fid, "%f, ", A1.Get(i, j).re);
		cout<<"Writing file 1: "<< i<<endl;
		fprintf(fid, ";\n");
	 }
	 fprintf(fid, "];");
	 fclose(fid);

	 fid = fopen("SphereQuad_rhs.txt", "w");
	fprintf(fid, "[");
	for(int i=0; i < numpts; i++)
	{
	 	fprintf(fid, "%f;\n", b1[i].re);
		cout<<"Writing file 2: "<< i<<endl;

	}
	fprintf(fid, "];");
	fclose(fid);
}
*/
/*void testPtsWtsForSphere(const int sphOrder, const int angN, const RDenseMatrix& pts, const CVector &wts, CDenseMatrix* &Yl2m2, CDenseMatrix* &Yl1m1)
{

  int numpts = pts.nRows();
  CDenseMatrix dnsAint(angN, angN);

  int is, js, indl1, indl2;

  for(int l1 = 0; l1 <= 3; l1++){
	indl1 = getPos(l1, -1*l1);
	for(int m1 = -1*l1; m1 <= l1; m1++){
	        			
		is = indl1 + l1 + m1;    	
		for(int l2 = 0; l2 <= 3; l2++){
			indl2 = getPos(l2, -1*l2);
			for(int m2 = -1*l2; m2 <= l2; m2++){
				js = indl2 + l2 + m2;

				for(int i=0; i<numpts; i++)
					dnsAint(is, js) += Yl1m1[l1](l1+m1, i)*Yl2m2[l2](l2+m2, i)*wts[i]*4*M_PI;
				 
				cout<<l1<<" "<<m1<<" "<<l2<<" "<<m2<<" "<<dnsAint(is, js)<<endl; 

			}
		}
	}
   }

  				

}*/

