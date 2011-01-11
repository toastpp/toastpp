/***************************************************************************
 * forward_rte1_3D.cc             Surya Mohan           05/05/2010         *
 *                                                                         *
 ***************************************************************************/

#ifdef __BORLANDC__
#include <strstrea.h>
#include <conio.h>
#include <process.h>

typedef unsigned pid_t;
#else
#include <sstream>
#include <unistd.h>

#endif
#include <iostream>
#include <fstream>
#include <string.h>
#include <algorithm>
#include <vector>
#include<set>
#include <time.h>
#include <iomanip.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mathlib.h>
#include "matrix.h"
#include <felib.h>
#include "source.h"
#include "pparse.h"
#include "rte3D_math.h"
#include "phaseFunc.h"
#define USE_INTONSPHERE
#define USE_SUBSAMPLING
#define VARYDELTA

using namespace toast;

#define MIN(A,B) ( (A) < (B) ? (A) : (B))
#define MAX(A,B) ( (A) > (B) ? (A) : (B))
#define Alm(l, m) (sqrt(complex((((l)+(m))*((l)-(m)))/((double)(2*(l)+1)*(2*(l)-1)), 0)))
#define Blm(l, m) (sqrt(complex((((l)+(m))*((l)+(m)-1))/((double)(2*(l)+1)*(2*(l)-1)), 0)))


class MyDataContext{
public:
RCompRowMatrix Sint, Sdx, Sdy, Sdz;
RCompRowMatrix Sgrad, Sx, Sy, Sz;
RCompRowMatrix Sdxx, Sdxy, Sdyx, Sdyy,  Sdxz, Sdzx,  Sdyz, Sdzy, Sdzz;
RCompRowMatrix Aint, Aintsc, Aintss, Aintc;
RCompRowMatrix Aintscsc,  Aintscss, Aintscc, Aintssss,  Aintssc, Aintcc;
RCompRowMatrix SPS, SPSdx, SPSdy, SPSdz;
RCompRowMatrix spatA3_rte, spatA3_sdmx, spatA3_sdmy, spatA3_sdmz;
RCompRowMatrix apu1, apu1sc, apu1ss, apu1c;
RCompRowMatrix A2, b1;
CVector A2x;
CDenseMatrix Aintx, Aintscx, Aintssx, Aintcx;
CDenseMatrix apu1x, apu1scx, apu1ssx, apu1cx;
CDenseMatrix Aintscscx, Aintscssx, Aintssssx, Aintsccx, Aintsscx, Aintccx;
int mfc_count;
int spatN, sphOrder, maxAngN;
IVector nodal_sphOrder, node_angN, offset;
const double *sintval, *sdxval, *sdyval, *sdzval, *sxval, *syval, *szval; 
const double *sdxxval, *sdxyval, *sdyxval, *sdyyval, *sdxzval, *sdzxval, *sdyzval, *sdzyval, *sdzzval; 
const double *spsval, *spsdxval, *spsdyval, *spsdzval, *spata3_rteval, *spata3_sdmxval, *spata3_sdmyval, *spata3_sdmzval; 

const double *aintval, *aintscval, *aintssval, *aintcval, *apu1val, *apu1scval, *apu1ssval, *apu1cval;
const double *aintscscval, *aintscssval, *aintssssval, *aintsccval, *aintsscval, *aintccval;
CCompRowMatrix Xmat;
RDenseMatrix *Ylm;
complex *xmatval;
double w, c;
MyDataContext(QMMesh &spatMesh, const IVector& nodal_sphOrder, RVector &delta, RVector &mua, RVector &mus, RVector &ref, const double g, complex (*phaseFunc)(const double g, const double costheta), double w, double c, const RDenseMatrix& pts, const RVector &wts)
{ 
	this->w = w;
	this->c = c;
	this->nodal_sphOrder = nodal_sphOrder;
	sphOrder = vmax(nodal_sphOrder);

	Ylm = new RDenseMatrix[sphOrder +1];
        	
	for(int l=0; l<=sphOrder; l++)
		Ylm[l].New(2*l+1, pts.nRows());
	sphericalHarmonics(sphOrder, pts.nRows(), pts, Ylm);


	/*RDenseMatrix *temp;
        RDenseMatrix tempMat(1, 3);
	tempMat(0, 0) = 0; tempMat(0, 1) = 1; tempMat(0, 2) = 0;
	cout<<"******Input**********"<<endl;
	cout<<tempMat<<endl;
	temp =  new RDenseMatrix[4];
        for(int l=0; l<=3; l++)
		temp[l].New(2*l+1, 1);
	sphericalHarmonics(3, 1, tempMat, temp);
	for(int l=0; l<=3; l++)
	{
		cout<<"***********Order "<<l<<"****************"<<endl;
		cout<<temp[l]<<endl;
	}*/
	  	
	spatN = spatMesh.nlen();
	        	
	node_angN.New(spatN);
	for(int i=0; i < spatN; i++)
	{
		node_angN[i] = 0;
		for(int l=0; l <= nodal_sphOrder[i]; l++)
			node_angN[i] += (2*l+1);
	}
	maxAngN = vmax(node_angN);

	offset.New(spatN);
	offset[0] = 0;
	for(int i=1; i < spatN; i++)
		offset[i] = offset[i-1] + node_angN[i-1]; 

	Aintx.New(spatN, maxAngN); Aintscx.New(spatN, maxAngN); Aintssx.New(spatN, maxAngN); Aintcx.New(spatN, maxAngN);
	apu1x.New(spatN, maxAngN); apu1scx.New(spatN, maxAngN); apu1ssx.New(spatN, maxAngN); apu1cx.New(spatN, maxAngN);
	Aintscscx.New(spatN, maxAngN); Aintscssx.New(spatN, maxAngN); Aintssssx.New(spatN, maxAngN); Aintsccx.New(spatN, maxAngN);
	Aintsscx.New(spatN, maxAngN); Aintccx.New(spatN, maxAngN);
	
	cout<<"Generating boundary integrals ..."<<endl;
       // genmat_boundint_3D(spatMesh, nodal_sphOrder, node_angN, offset, pts, wts, Aintsc, Aintss, Aintc, A2, b1);
	genmat_boundint_3D(spatMesh, nodal_sphOrder, node_angN, offset, pts, wts, Ylm, A2, b1);
	cout<<"Done computing boundary integrals ..."<<endl;	

	gen_spatint_3D(spatMesh, mua, mus, ref, delta, w, c, Sint, Sdx, Sdy, Sdz, Sx, Sy, Sz, Sdxx, Sdxy, Sdyx, Sdyy, Sdxz, Sdzx, Sdyz, Sdzy, Sdzz, spatA3_rte, spatA3_sdmx, spatA3_sdmy, spatA3_sdmz, SPS, SPSdx, SPSdy, SPSdz);
	cout<<"Generated Spatial Integrals ..."<<endl;
	int angN = vmax(node_angN);
      	genmat_angint_3D(sphOrder, angN, Aint, Aintsc, Aintss, Aintc, Aintscsc, Aintscss, Aintscc,  Aintssss, Aintssc, Aintcc);
	cout<<"Generated angular integrals ..."<<endl;
	cout<<"Generating phase integrals ..."<<endl;	
	genmat_apu(phaseFunc, g, angN, sphOrder, apu1, apu1sc, apu1ss, apu1c);
	
	//fwrite(sphOrder, angN, Aint, Aintsc,  Aintss, Aintc, Aintscsc, Aintscss, Aintscc,  Aintssss, Aintssc, Aintcc, apu1, apu1sc, apu1ss, apu1c);

	/*cout<<"A: "<<Aint.rowptr[angN]<<" Ax: "<<Aintsc.rowptr[angN]<<" Ay: "<<Aintss.rowptr[angN]<<" Az: "<<Aintc.rowptr[angN]<<endl;
	cout<<"Axx: "<<Aintscsc.rowptr[angN]<<" Axy: "<<Aintscss.rowptr[angN]<<" Axz: "<<Aintscc.rowptr[angN]<<" Ayy: "<<Aintssss.rowptr[angN]<<endl;
	cout<<"Ayz: "<<Aintssc.rowptr[angN]<<" Azz: "<<Aintcc.rowptr[angN]<<endl;

	cout<<"B: "<<apu1.rowptr[angN]<<" Bx: "<<apu1sc.rowptr[angN]<<" By: "<<apu1ss.rowptr[angN]<<" Bz: "<<apu1c.rowptr[angN]<<endl;*/


	
	//cout<<apu1<<"\n******************\n"<<apu1sc<<"\n******************\n"<<apu1ss<<"\n******************\n"<<apu1c<<endl;

/*	CCompRowMatrix apu1_test, apu1sc_test, apu1ss_test, apu1c_test;	
	genmat_apu_test(phaseFunc, g, angN, sphOrder, pts, wts, apu1_test, apu1sc_test, apu1ss_test, apu1c_test);

	cout<<apu1<<"\n******************\n"<<apu1_test<<"\n******************\n"<<apu1sc<<"\n******************\n"<<apu1sc_test<<endl;
	cout<<apu1ss<<"\n******************\n"<<apu1ss_test<<"\n******************\n"<<apu1c<<"\n******************\n"<<apu1c_test<<endl;
*/

        apu1.Transpone(); apu1sc.Transpone(); apu1ss.Transpone(); apu1c.Transpone(); 	
	Aint.Transpone(); Aintsc.Transpone(); Aintss.Transpone(); Aintc.Transpone();   
	Aintscsc.Transpone(); Aintscss.Transpone(); Aintscc.Transpone(); Aintssc.Transpone(); 
	Aintcc.Transpone(); Aintssss.Transpone();

        
	aintval = Aint.ValPtr(); aintscval = Aintsc.ValPtr(); aintssval = Aintss.ValPtr(); aintcval = Aintc.ValPtr();
	apu1val = apu1.ValPtr(); apu1scval = apu1sc.ValPtr(); apu1ssval = apu1ss.ValPtr(); apu1cval = apu1c.ValPtr();
	aintscscval = Aintscsc.ValPtr(); aintscssval = Aintscss.ValPtr(); aintssssval = Aintssss.ValPtr(); aintsccval = Aintscc.ValPtr();
	aintsscval = Aintssc.ValPtr(); aintccval = Aintcc.ValPtr();


	sintval = Sint.ValPtr(); sdxval = Sdx.ValPtr(); sdyval = Sdy.ValPtr();
	sdzval = Sdz.ValPtr(); sxval = Sx.ValPtr(); syval = Sy.ValPtr(); 
	szval = Sz.ValPtr(); sdxxval = Sdxx.ValPtr(); sdxyval = Sdxy.ValPtr(); 
	sdyxval = Sdyx.ValPtr(); sdyyval = Sdyy.ValPtr(); sdxzval = Sdxz.ValPtr(); 
	sdzxval = Sdzx.ValPtr(); sdyzval = Sdyz.ValPtr(); sdzyval = Sdzy.ValPtr(); 
	sdzzval = Sdzz.ValPtr(); spsval = SPS.ValPtr(); spsdxval = SPSdx.ValPtr(); 
	spsdyval = SPSdy.ValPtr(); spsdzval = SPSdz.ValPtr(); spata3_rteval = spatA3_rte.ValPtr(); 
	spata3_sdmxval = spatA3_sdmx.ValPtr(); spata3_sdmyval = spatA3_sdmy.ValPtr(); 
	spata3_sdmzval = spatA3_sdmz.ValPtr();      
	


	A2x.New(sum(node_angN));

	Xmat.New(spatN, maxAngN);
	int *xrowptr = new int[spatN + 1];
	
	xrowptr[0]=0;
	for(int i=1;i < spatN+1; i++)
		xrowptr[i] = xrowptr[i-1] + node_angN[i-1];
	
	int *xcolidx = new int[xrowptr[spatN]];
	int k=0;
	for(int i = 0; i < spatN; i++)
	{
		for(int j=0; j < node_angN[i]; j++)
		{
			xcolidx[k] = j;
			k++;
		}
	 }

	Xmat.Initialise(xrowptr, xcolidx);
	xmatval = Xmat.ValPtr();

}

~MyDataContext()
{
  delete []Ylm;
};

/*void fwrite(const int sphOrder, const int angN, CCompRowMatrix& Aint, CCompRowMatrix& Aintsc, CCompRowMatrix& Aintss, CCompRowMatrix& Aintc, CCompRowMatrix& Aintscsc, CCompRowMatrix& Aintscss, CCompRowMatrix& Aintscc, CCompRowMatrix& Aintssss, CCompRowMatrix& Aintssc, CCompRowMatrix& Aintcc, CCompRowMatrix& apu1, CCompRowMatrix& apu1sc, CCompRowMatrix& apu1ss, CCompRowMatrix& apu1c)
{

 writeOneFile(angN, "Aint.dat", Aint); writeOneFile(angN, "Aintsc.dat", Aintsc); writeOneFile(angN, "Aintss.dat", Aintss); 
 writeOneFile(angN, "Aintc.dat", Aintc); writeOneFile(angN, "Aintscsc.dat", Aintscsc); writeOneFile(angN, "Aintscss.dat", Aintscss); 
 writeOneFile(angN, "Aintssss.dat", Aintssss); writeOneFile(angN, "Aintscc.dat", Aintscc); writeOneFile(angN, "Aintssc.dat", Aintssc);
 writeOneFile(angN, "Aintcc.dat", Aintcc); writeOneFile(angN, "apu1.dat", apu1); writeOneFile(angN, "apu1sc.dat", apu1sc);
 writeOneFile(angN, "apu1ss.dat", apu1ss); writeOneFile(angN, "apu1c.dat", apu1c); 
 
 



  cout<<"FILES WRITTEN"<<endl;
}

void writeOneFile(const int angN, const char *fname, CCompRowMatrix& A)
{
	FILE *fid;
	fid = fopen(fname, "w");
	complex *aval = A.ValPtr();
  	for(int i=0; i<angN; i++)
  	{
		for(int na = A.rowptr[i]; na < A.rowptr[i+1]; na++)
		{	
			int col = A.colidx[na];
			double re = aval[na].re;
			fprintf(fid, "%d %d %f\n", i+1, col+1, re);
		}
  	}  
  	fclose(fid);

}*/
void gen_spatint_3D(const QMMesh& mesh, const RVector& muabs, const RVector& muscat, const RVector& ref, const RVector& delta, double w, double c, RCompRowMatrix& Sint, RCompRowMatrix& Sdx, RCompRowMatrix& Sdy, RCompRowMatrix& Sdz, RCompRowMatrix& Sx, RCompRowMatrix& Sy, RCompRowMatrix& Sz, RCompRowMatrix& Sdxx, RCompRowMatrix& Sdxy, RCompRowMatrix& Sdyx, RCompRowMatrix& Sdyy, RCompRowMatrix& Sdxz, RCompRowMatrix& Sdzx, RCompRowMatrix& Sdyz, RCompRowMatrix& Sdzy, RCompRowMatrix& Sdzz, RCompRowMatrix& spatA3_rte, RCompRowMatrix& spatA3_sdmx, RCompRowMatrix& spatA3_sdmy, RCompRowMatrix& spatA3_sdmz, RCompRowMatrix& SPS, RCompRowMatrix& SPSdx, RCompRowMatrix& SPSdy, RCompRowMatrix& SPSdz)
{
   int sysdim = mesh.nlen();       // dimensions are size of nodes.

   int *rowptr, *colidx, nzero;
   mesh.SparseRowStructure (rowptr, colidx, nzero);
   
   Sint.New (sysdim, sysdim);
   Sint.Initialise (rowptr, colidx);
   Sx.New (sysdim, sysdim);
   Sx.Initialise (rowptr, colidx);
   Sy.New (sysdim, sysdim);
   Sy.Initialise (rowptr, colidx);
   Sz.New (sysdim, sysdim);
   Sz.Initialise (rowptr, colidx);
   SPS.New (sysdim, sysdim);
   SPS.Initialise (rowptr, colidx); 
   spatA3_rte.New (sysdim, sysdim);
   spatA3_rte.Initialise (rowptr, colidx);
   Sdx.New (sysdim, sysdim);
   Sdx.Initialise (rowptr, colidx);
   Sdy.New (sysdim, sysdim);
   Sdy.Initialise (rowptr, colidx);
   Sdz.New (sysdim, sysdim);
   Sdz.Initialise (rowptr, colidx);
   Sdxx.New (sysdim, sysdim);
   Sdxx.Initialise (rowptr, colidx);
   Sdxy.New(sysdim, sysdim);
   Sdxy.Initialise (rowptr, colidx);
   Sdyx.New (sysdim, sysdim);
   Sdyx.Initialise (rowptr, colidx);
   Sdyy.New (sysdim, sysdim);
   Sdyy.Initialise (rowptr, colidx);
   Sdxz.New (sysdim, sysdim);
   Sdxz.Initialise (rowptr, colidx);
   Sdzx.New(sysdim, sysdim);
   Sdzx.Initialise(rowptr, colidx);
   Sdyz.New (sysdim, sysdim);
   Sdyz.Initialise (rowptr, colidx);
   Sdzy.New(sysdim, sysdim);
   Sdzy.Initialise(rowptr, colidx);
   Sdzz.New (sysdim, sysdim);
   Sdzz.Initialise (rowptr, colidx);
   SPSdx.New (sysdim, sysdim);
   SPSdx.Initialise (rowptr, colidx); 
   spatA3_sdmx.New (sysdim, sysdim);
   spatA3_sdmx.Initialise (rowptr, colidx);
   SPSdy.New (sysdim, sysdim);
   SPSdy.Initialise (rowptr, colidx); 
   spatA3_sdmy.New (sysdim, sysdim);
   spatA3_sdmy.Initialise (rowptr, colidx);
   SPSdz.New (sysdim, sysdim);
   SPSdz.Initialise (rowptr, colidx); 
   spatA3_sdmz.New (sysdim, sysdim);
   spatA3_sdmz.Initialise (rowptr, colidx);

   
   double elk_ij, elb_ij, elsx_ij, elsy_ij, elsz_ij;
   int el, nodel, i, j, k,is, js;


   
   for (el = 0; el < mesh.elen(); el++) {
	nodel = mesh.elist[el]->nNode();
	
	double dss = delta[el]; //streamline diffusion value for this element.
	double sigmatot = muabs[el] + muscat[el];
	
	RSymMatrix eldd = mesh.elist[el]->Intdd(); // "all at once!""
	for (i = 0; i < nodel; i++) {
	    if ((is = mesh.elist[el]->Node[i]) >= sysdim) continue;
	    for (j = 0; j < nodel; j++) {
		if ((js = mesh.elist[el]->Node[j]) >= sysdim) continue;
		
		elb_ij = mesh.elist[el]->IntFF (i, j);
		Sint(is,js) +=elb_ij*w/c; 
		SPS(is, js) += elb_ij*muscat[el];
		spatA3_rte(is, js) += elb_ij*sigmatot;


		elsx_ij = mesh.elist[el]->IntFd (j,i,0);
		Sx(is,js) += elsx_ij;
		Sdx(is,js) += dss*elsx_ij*w/c;
		SPSdx(is, js) += dss*elsx_ij*muscat[el];
		spatA3_sdmx(is, js) += dss*elsx_ij*sigmatot;

		elsy_ij = mesh.elist[el]->IntFd (j,i,1);
		Sy(is,js) += elsy_ij;
		Sdy(is,js) += dss*elsy_ij*w/c;
		SPSdy(is, js) += dss*elsy_ij*muscat[el];
		spatA3_sdmy(is, js) += dss*elsy_ij*sigmatot;

		elsz_ij = mesh.elist[el]->IntFd (j,i,2);
		Sz(is,js) += elsz_ij;
  		Sdz(is,js) += dss*elsz_ij*w/c;
		SPSdz(is, js) += dss*elsz_ij*muscat[el];
		spatA3_sdmz(is, js) += dss*elsz_ij*sigmatot;

		Sdxx(is,js) += dss * eldd(i*3,j*3);
	       	Sdxy(is,js) += dss * eldd(i*3,j*3+1);
     		Sdyx(is,js) += dss * eldd(i*3+1,j*3);
       		Sdyy(is,js) += dss * eldd(i*3+1,j*3+1);	
	       	Sdxz(is,js) += dss * eldd(i*3,j*3+2);
     		Sdzx(is,js) += dss * eldd(i*3+2,j*3);
	       	Sdyz(is,js) += dss * eldd(i*3+1,j*3+2);
     		Sdzy(is,js) += dss * eldd(i*3+2,j*3+1);
       		Sdzz(is,js) += dss * eldd(i*3+2,j*3+2);	
	    }
	}
   }

  delete []rowptr;
  delete []colidx;
}
void sincosY(const int l, const int m, CVector& a, CVector& b, CVector& c, CVector& d, IDenseMatrix& a1c, IDenseMatrix& b1c, IDenseMatrix& c1c, IDenseMatrix& d1c)
{

	if(m==0)
	{
		a[0] = Blm(l, 0)/complex(2.0, 0); b[0] = complex(-1, 0)*Blm(l+1, 1)/complex(2.0, 0); 
		c[0] = complex(-1, 0)*Blm(l, 0)/complex(2.0, 0); d[0] = Blm(l+1, 1)/complex(2.0, 0);
		
		a[1] = 0; b[1] = 0; 
		c[1] = 0; d[1] = 0;
		
		a1c(0, 0) = l-1; a1c(0, 1) = 1; b1c(0, 0) = l+1; b1c(0, 1) = 1;
		c1c(0, 0) = l-1; c1c(0, 1) = -1; d1c(0, 0) = l+1; d1c(0, 1) = -1; 
		
		a1c(1, 0) = -1; a1c(1, 1) = -1; b1c(1, 0) = -1; b1c(1, 1) = -1;
		c1c(1, 0) = -1; c1c(1, 1) = -1; d1c(1, 0) = -1; d1c(1, 1) = -1; 				
	}
	else if(m>0)
	{
	    	a[0] = Blm(l, -m)/(complex(2.0*sqrt(2), 0)); b[0] =  complex(-1, 0)*Blm(l+1, m+1)/(complex(2.0*sqrt(2), 0)); 
	    	c[0] = complex(-1, 0)*Blm(l, m)/(complex(2.0*sqrt(2), 0)); d[0] = Blm(l+1, -m+1)/(complex(2.0*sqrt(2), 0));
	    	a[1] = complex(sign(m), 0)*Blm(l, m)/(complex(2.0*sqrt(2), 0)); b[1] = complex(-1*sign(m), 0)*Blm(l+1, -m+1)/(complex(2.0*sqrt(2), 0));
		c[1] =  complex(-1*sign(m), 0)*Blm(l, -m)/(complex(2.0*sqrt(2), 0)); d[1] = complex(sign(m), 0)*Blm(l+1, m+1)/(complex(2.0*sqrt(2), 0));

	   	a1c(0, 0) = l-1; a1c(0, 1) = m+1; b1c(0, 0) = l+1; b1c(0, 1) = m+1;
	   	c1c(0, 0) = l-1; c1c(0, 1) = m-1; d1c(0, 0) = l+1; d1c(0, 1) = m-1;

	   	a1c(1, 0) = l-1; a1c(1, 1) = -m+1; b1c(1, 0) = l+1; b1c(1, 1) = -m+1;
	   	c1c(1, 0) = l-1; c1c(1, 1) = -m-1; d1c(1, 0) = l+1; d1c(1, 1) = -m-1; 
 
	}
	else
	{
	   	a[0] =  Blm(l, m)/(complex(0, 2.0*sqrt(2))); b[0] =  complex(-1, 0)*Blm(l+1, -m+1)/(complex(0, 2.0*sqrt(2))); 
	    	c[0] =   complex(-1, 0)*Blm(l, -m)/(complex(0, 2.0*sqrt(2))); d[0] = Blm(l+1, m+1)/(complex(0, 2.0*sqrt(2)));
	    	a[1] = complex(sign(m+1), 0)*Blm(l, -m)/(complex(0, 2.0*sqrt(2))); b[1] = complex(-1*sign(m+1), 0)*Blm(l+1, m+1)/(complex(0, 2.0*sqrt(2)));
		c[1] =  complex(-1*sign(m+1), 0)*Blm(l, m)/(complex(0, 2.0*sqrt(2))); d[1] =  complex(sign(m+1), 0)*Blm(l+1, -m+1)/(complex(0, 2.0*sqrt(2)));


	   	a1c(0, 0) = l-1; a1c(0, 1) = -m+1; b1c(0, 0) = l+1; b1c(0, 1) = -m+1;
	   	c1c(0, 0) = l-1; c1c(0, 1) = -m-1; d1c(0, 0) = l+1; d1c(0, 1) = -m-1;

	   	a1c(1, 0) = l-1; a1c(1, 1) = m+1; b1c(1, 0) = l+1; b1c(1, 1) = m+1;
	   	c1c(1, 0) = l-1; c1c(1, 1) = m-1; d1c(1, 0) = l+1; d1c(1, 1) = m-1; 
	}

}


void sinsinY(const int l, const int m, CVector& a, CVector& b, CVector& c, CVector& d, IDenseMatrix& a1c, IDenseMatrix& b1c, IDenseMatrix& c1c, IDenseMatrix& d1c)
{

	if(m==0)
	{
		a[0] = Blm(l, 0)/complex(0, 2.0); b[0] = complex(-1, 0)*Blm(l+1, 1)/complex(0, 2.0); 
		c[0] = Blm(l, 0)/complex(0, 2.0); d[0] = complex(-1, 0)*Blm(l+1, 1)/complex(0, 2.0);
		
		a[1] = 0; b[1] = 0; 
		c[1] = 0; d[1] = 0;
		
		a1c(0, 0) = l-1; a1c(0, 1) = 1; b1c(0, 0) = l+1; b1c(0, 1) = 1;
		c1c(0, 0) = l-1; c1c(0, 1) = -1; d1c(0, 0) = l+1; d1c(0, 1) = -1; 
		
		a1c(1, 0) = -1; a1c(1, 1) = -1; b1c(1, 0) = -1; b1c(1, 1) = -1;
		c1c(1, 0) = -1; c1c(1, 1) = -1; d1c(1, 0) = -1; d1c(1, 1) = -1; 				
	}
	else if(m>0)
	{
	    	a[0] = Blm(l, -m)/(complex(0, 2.0*sqrt(2))); b[0] =  complex(-1, 0)*Blm(l+1, m+1)/(complex(0, 2.0*sqrt(2))); 
	    	c[0] = Blm(l, m)/(complex(0, 2.0*sqrt(2))); d[0] = complex(-1, 0)*Blm(l+1, -m+1)/(complex(0, 2.0*sqrt(2)));
	    	a[1] = complex(sign(m), 0)*Blm(l, m)/(complex(0, 2.0*sqrt(2))); b[1] = complex(-1*sign(m), 0)*Blm(l+1, -m+1)/(complex(0, 2.0*sqrt(2)));
		c[1] =  complex(sign(m), 0)*Blm(l, -m)/(complex(0, 2.0*sqrt(2))); d[1] = complex(-1*sign(m), 0)*Blm(l+1, m+1)/(complex(0, 2.0*sqrt(2)));

	   	a1c(0, 0) = l-1; a1c(0, 1) = m+1; b1c(0, 0) = l+1; b1c(0, 1) = m+1;
	   	c1c(0, 0) = l-1; c1c(0, 1) = m-1; d1c(0, 0) = l+1; d1c(0, 1) = m-1;

	   	a1c(1, 0) = l-1; a1c(1, 1) = -m+1; b1c(1, 0) = l+1; b1c(1, 1) = -m+1;
	   	c1c(1, 0) = l-1; c1c(1, 1) = -m-1; d1c(1, 0) = l+1; d1c(1, 1) = -m-1; 
 
	}
	else
	{
	   	a[0] =  complex(-1, 0)*Blm(l, m)/(complex(2.0*sqrt(2), 0)); b[0] =  Blm(l+1, -m+1)/(complex(2.0*sqrt(2), 0)); 
	    	c[0] =   complex(-1, 0)*Blm(l, -m)/(complex(2.0*sqrt(2), 0)); d[0] = Blm(l+1, m+1)/(complex(2.0*sqrt(2), 0));
	    	a[1] = complex(sign(m), 0)*Blm(l, -m)/(complex(2.0*sqrt(2), 0)); b[1] = complex(-1*sign(m), 0)*Blm(l+1, m+1)/(complex(2.0*sqrt(2), 0));
		c[1] =  complex(sign(m), 0)*Blm(l, m)/(complex(2.0*sqrt(2), 0)); d[1] =  complex(-1*sign(m), 0)*Blm(l+1, -m+1)/(complex(2.0*sqrt(2), 0));


	   	a1c(0, 0) = l-1; a1c(0, 1) = -m+1; b1c(0, 0) = l+1; b1c(0, 1) = -m+1;
	   	c1c(0, 0) = l-1; c1c(0, 1) = -m-1; d1c(0, 0) = l+1; d1c(0, 1) = -m-1;

	   	a1c(1, 0) = l-1; a1c(1, 1) = m+1; b1c(1, 0) = l+1; b1c(1, 1) = m+1;
	   	c1c(1, 0) = l-1; c1c(1, 1) = m-1; d1c(1, 0) = l+1; d1c(1, 1) = m-1; 
	}
}

void cosY(const int l, const int m, CVector& e, CVector& f, IDenseMatrix& e1c, IDenseMatrix& f1c)
{
	if(m==0)
	{
		e[0] = Alm(l, 0); f[0] = Alm(l+1, 0);
		e[1] = 0; f[1] = 0;
		
		e1c(0, 0) = l-1; e1c(0, 1) = 0;
		f1c(0, 0) = l+1; f1c(0, 1) = 0;
		
		e1c(1, 0) = -1; e1c(1, 1) = -1;
		f1c(1, 0) = -1; f1c(1, 1) = -1;

	}	
	else if(m>0)
	{
		e[0] = Alm(l, m)/complex(sqrt(2), 0); f[0] = Alm(l+1, m)/complex(sqrt(2), 0);
		e[1] = complex(sign(m), 0)*Alm(l, -m)/complex(sqrt(2), 0); f[1] = complex(sign(m), 0)*Alm(l+1, -m)/complex(sqrt(2), 0);
	
		e1c(0, 0) = l-1; e1c(0, 1) = m;
		f1c(0, 0) = l+1; f1c(0, 1) = m;

		e1c(1, 0) = l-1; e1c(1, 1) = -m;
		f1c(1, 0) = l+1; f1c(1, 1) = -m;
	}
	else
	{
		e[0] = Alm(l, -m)/complex(0, sqrt(2)); f[0] = Alm(l+1, -m)/complex(0, sqrt(2));
		e[1] = complex(-1*sign(m), 0)*Alm(l, m)/complex(0, sqrt(2)); f[1] = complex(-1*sign(m), 0)*Alm(l+1, m)/complex(0, sqrt(2));
	
		e1c(0, 0) = l-1; e1c(0, 1) = -m;
		f1c(0, 0) = l+1; f1c(0, 1) = -m;

		e1c(1, 0) = l-1; e1c(1, 1) = m;
		f1c(1, 0) = l+1; f1c(1, 1) = m;

	}

}

void sphY(const int l, const int m, CVector& p, IDenseMatrix& p1c)
{
	if(m==0)
	{
		p[0] = 1.0; p[1] = 0;
		
		p1c(0, 0) = l; p1c(0, 1) = m;
		p1c(1, 0) = -1; p1c(1, 1) = -3;
		
	}	
	else if(m>0)
	{
		p[0] = complex(1.0, 0)/complex(sqrt(2), 0); p[1] = complex(sign(m), 0)/complex(sqrt(2), 0);
	
		p1c(0, 0) = l; p1c(0, 1) = m;
		p1c(1, 0) = l; p1c(1, 1) = -m;

	}
	else
	{
		p[0] = complex(1.0, 0)/complex(0, sqrt(2)); p[1] = complex(-1*sign(m), 0)/complex(0, sqrt(2));
	
		p1c(0, 0) = l; p1c(0, 1) = -m;
		p1c(1, 0) = l; p1c(1, 1) = m;

	}
}

double Integrate2(CVector &a, CVector &b, IDenseMatrix &a1c, IDenseMatrix &b1c)
{
   complex temp;
   temp =  a[0]*b[0]*kronD(a1c.Row(0), b1c.Row(0)) + a[0]*b[1]*kronD(a1c.Row(0), b1c.Row(1)) + a[1]*b[0]*kronD(a1c.Row(1), b1c.Row(0)) + a[1]*b[1]*kronD(a1c.Row(1), b1c.Row(1));
   return(temp.re);
}

void genmat_angint_3D(const int sphOrder, const int angN, RCompRowMatrix& Aint, RCompRowMatrix& Aintsc, RCompRowMatrix& Aintss, RCompRowMatrix& Aintc, RCompRowMatrix& Aintscsc, RCompRowMatrix& Aintscss, RCompRowMatrix& Aintscc, RCompRowMatrix& Aintssss, RCompRowMatrix& Aintssc, RCompRowMatrix& Aintcc)
{
	RDenseMatrix dnsAint(angN, angN), dnsAintsc(angN, angN), dnsAintss(angN, angN), dnsAintc(angN, angN);
	RDenseMatrix dnsAintscsc(angN, angN), dnsAintscss(angN, angN), dnsAintcc(angN, angN);
	RDenseMatrix dnsAintscc(angN, angN), dnsAintssss(angN, angN), dnsAintssc(angN, angN);

	CVector a1(2), b1(2), c1(2), d1(2), e1(2), f1(2), p1(2), a2(2), b2(2), c2(2), d2(2), e2(2), f2(2), p2(2);
        IDenseMatrix a1c(2, 2), b1c(2, 2), c1c(2, 2), d1c(2, 2), e1c(2, 2), f1c(2, 2), p1c(2, 2);
	IDenseMatrix a2c(2, 2), b2c(2, 2), c2c(2, 2), d2c(2, 2), e2c(2, 2), f2c(2, 2), p2c(2, 2);
	 

	int is, js, indl1, indl2;

	for(int l1 = 0; l1 <= sphOrder; l1++){
		indl1 = getPos(l1, -1*l1);
		cout<<"angint: order: "<<l1<<endl;
		for(int m1 = -1*l1; m1 <= l1; m1++){
		        			
			is = indl1 + l1 + m1;    	
			for(int l2 = 0; l2 <= sphOrder; l2++){
				indl2 = getPos(l2, -1*l2);
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

RVector phaseFuncDisc(int sphOrder, complex (*phaseFunc)(const double g, const double costheta), const double g)
{
 RVector phaseFn(sphOrder+1);
 /*int numpts = pts.nRows();
 RDenseMatrix LT;
 RVector vec;
 for(int l=0; l <= sphOrder; l++)
 {	
	LT = LegendreTable(l, numpts, pts);
	//vec = LT.ColSum();
	for(int i = 0; i < numpts; i++)
		phaseFn[l] += wts[i]*complex(LT.Get(l, i), 0.0)*phaseFunc(g, pts.Get(i, 2));
  }
	return(phaseFn); 	*/
 for(int l=0; l <= sphOrder; l++)
		phaseFn[l] = pow(g, (double)l);
 return(phaseFn);	
}

complex am(int m)
{
	if(m==0) return(complex(1.0, 0.0));
	else if (m > 0) return(complex(1/sqrt(2.0), 0.0));
	else return(complex(0.0, sign(m)/sqrt(2.0))); 
}

complex bm(int m)
{
	if(m==0) return(complex(0.0, 0.0));
	else if (m > 0) return(complex(1/sqrt(2.0), 0.0));
	else return(complex(0.0, -1*sign(m)/sqrt(2.0))); 
}

complex delta(int a, int b)
{
	return(a == b ? complex(1.0, 0.0) : complex(0.0, 0.0)); 
}

/*void genmat_apu(complex (*phaseFunc)(const double g, const double costheta), const double g, const int angN, const int sphOrder, RCompRowMatrix& apu1, RCompRowMatrix& apu1sc, RCompRowMatrix& apu1ss, RCompRowMatrix& apu1c)
{
    RDenseMatrix dnsapu1(angN, angN), dnsapu1sc(angN, angN), dnsapu1ss(angN, angN), dnsapu1c(angN, angN);
    RVector phaseFn;
    phaseFn = phaseFuncDisc(sphOrder, phaseFunc, g);
    cout<<"****************** Discretized phase function ****************"<<endl;
    cout<<phaseFn<<endl;
    cout<<"**************************************************************"<<endl;
    int indl, indm;
    int indl1, indm1, indl2, indm2, is, js;
    complex temp;

    for(int l1=0; l1 <= sphOrder; l1++){
    	indl1 =  getPos(l1, -1*l1);
        cout<<"apu: order: "<<l1<<endl;
	for(int m1 = -1*l1; m1 <= l1; m1++){
		indm1 = l1 + m1;
		is = indl1 + indm1;
	        for(int l2=0; l2<= sphOrder; l2++){
			indl2 = getPos(l2, -1*l2);

			for(int m2 = -1*l2; m2<=l2; m2++){
				indm2 = l2 + m2;
				js = indl2 + indm2;
				
				
				for(int l=0; l <= sphOrder; l++)
    				{
					for(int m=-1*l; m <=l; m++)
					{
						complex I1 = (am(m1)*delta(-m, m1)*sign(m) + bm(m1)*delta(m, m1));
						complex I2 = delta(l, l2)*(am(m2)*delta(m, m2) + bm(m2)*sign(m)*delta(-m, m2));
						 
					 	temp = I2*delta(l, l1)*I1*phaseFn[l];
						dnsapu1(is, js) += temp.re;

						complex sincosyy = (Blm(l1, -m1)*delta(l, l1-1) - Blm(l1+1, m1+1)*delta(l, l1+1))*(am(m1)*delta(-m, m1+1)*sign(m) + bm(m1)*delta(m, m1+1));
						sincosyy += (-Blm(l1, m1)*delta(l, l1-1) + Blm(l1+1, -m1+1)*delta(l, l1+1))*(am(m1)*delta(-m, m1-1)*sign(m) + bm(m1)*delta(m, m1-1));
						temp = I2*0.5*sincosyy*phaseFn[l]; 
						dnsapu1sc(is, js) += temp.re;

						complex sinsinyy = (Blm(l1, -m1)*delta(l, l1-1) - Blm(l1+1, m1+1)*delta(l, l1+1))*(am(m1)*delta(-m, m1+1)*sign(m) - bm(m1)*delta(m, m1+1));
						sinsinyy += (Blm(l1, m1)*delta(l, l1-1) - Blm(l1+1, -m1+1)*delta(l, l1+1))*(am(m1)*delta(-m, m1-1)*sign(m) - bm(m1)*delta(m, m1-1));
						temp = complex(0, -0.5)*phaseFn[l]*I2*sinsinyy;
						dnsapu1ss(is, js) += temp.re;

						complex cosyy = (Alm(l1, m1)*delta(l, l1-1) + Alm(l1+1, m1)*delta(l, l1+1))*I1;
						
						temp = I2*cosyy*phaseFn[l];
						dnsapu1c(is, js) += temp.re; 
					}
				}
				//cout<<l1<<" "<<m1<<" "<<l2<<" "<<m2<<" "<<dnsapu1c(is, js)<<endl;
			}
		}
	}
   }

   
   apu1 = shrink(dnsapu1); apu1sc = shrink(dnsapu1sc); apu1ss = shrink(dnsapu1ss); apu1c = shrink(dnsapu1c); 	

}*/

void genmat_apu(complex (*phaseFunc)(const double g, const double costheta), const double g, const int angN, const int sphOrder, RCompRowMatrix& apu1, RCompRowMatrix& apu1sc, RCompRowMatrix& apu1ss, RCompRowMatrix& apu1c)
{
 RDenseMatrix dnsapu1(angN, angN), dnsapu1sc(angN, angN), dnsapu1ss(angN, angN), dnsapu1c(angN, angN);
     RVector phaseFn;
    phaseFn = phaseFuncDisc(sphOrder, phaseFunc, g);
    cout<<"****************** Discretized phase function ****************"<<endl;
    cout<<phaseFn<<endl;
    cout<<"**************************************************************"<<endl;
    int indl, indm;
    int indl1, indm1, indl2, indm2, is, js;
    CVector a1(2), b1(2), c1(2), d1(2), e1(2), f1(2), p1(2), p2(2), p(2);
    IDenseMatrix a1c(2, 2), b1c(2, 2), c1c(2, 2), d1c(2, 2), e1c(2, 2), f1c(2, 2), p1c(2, 2), p2c(2, 2), pc(2, 2);
 
    for(int l1=0; l1 <= sphOrder; l1++){
    	indl1 =  getPos(l1, -1*l1);
	cout<<"apu: order: "<<l1<<endl;
	for(int m1 = -1*l1; m1 <= l1; m1++){
		indm1 = l1 + m1;
		is = indl1 + indm1;
	        for(int l2=0; l2<= sphOrder; l2++){
			indl2 = getPos(l2, -1*l2);

			for(int m2 = -1*l2; m2<=l2; m2++){
				indm2 = l2 + m2;
				js = indl2 + indm2;
				

				//for(int l=0; l <= sphOrder; l++)
    				//{
				//	for(int m=-1*l; m <=l; m++)
				//	{
				//		sphY(l, m, p, pc);
				//		sphY(l1, m1, p1, p1c);
				//		sphY(l2, m2, p2, p2c);
				//		dnsapu1(is, js) += phaseFn[l]*Integrate2(p, p1, pc, p1c)*Integrate2(p, p2, pc, p2c);
				//		sincosY(l1, m1, a1, b1, c1, d1, a1c, b1c, c1c, d1c);
				//		dnsapu1sc(is, js) += phaseFn[l]*Integrate2(p, p2, pc, p2c)*(Integrate2(a1, p, a1c, pc) + Integrate2(b1, p, b1c, pc) + Integrate2(c1, p, c1c, pc) + Integrate2(d1, p, d1c, pc));
				//		sinsinY(l1, m1, a1, b1, c1, d1, a1c, b1c, c1c, d1c);
				//		dnsapu1ss(is, js) += phaseFn[l]*Integrate2(p, p2, pc, p2c)*(Integrate2(a1, p, a1c, pc) + Integrate2(b1, p, b1c, pc) + Integrate2(c1, p, c1c, pc) + Integrate2(d1, p, d1c, pc));

				//		cosY(l1, m1, e1, f1, e1c, f1c);
				//		dnsapu1c(is, js) += phaseFn[l]*Integrate2(p, p2, pc, p2c)*(Integrate2(e1, p, e1c, pc) + Integrate2(f1, p, f1c, pc));

				//	}
				//}

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


void initialiseA2b1(const Mesh &mesh, const IVector &node_angN, const IVector &offset, RCompRowMatrix& A2, RCompRowMatrix& b1)
{
   int el, nodel, i, j, k,is, js;
   int *crrowptr, *crcolidx, nzero;
   int sysdim = mesh.nlen();
   mesh.SparseRowStructure(crrowptr, crcolidx, nzero);
   
   int *status = new int[crrowptr[sysdim]];
   for(i=0; i<nzero; i++)
	status[i] = 0; // 1 implies nonzero 0 denotes zero;
   
   for (el = 0; el < mesh.elen(); el++) {
        if(!mesh.elist[el]->HasBoundarySide ()) continue;
	nodel = mesh.elist[el]->nNode();
	// now determine the element integrals
	for(int sd = 0; sd <  mesh.elist[el]->nSide(); sd++)  {
	  // if sd is not a boundary side. skip 
	  if(!mesh.elist[el]->IsBoundarySide (sd)) continue;
	  for (int i = 0; i < nodel; i++) {
	    if ((is = mesh.elist[el]->Node[i]) >= sysdim) continue;
	    for (int j = 0; j < nodel; j++) {
		if ((js = mesh.elist[el]->Node[j]) >= sysdim) continue;
		for (int rp = crrowptr[is]; rp < crrowptr[is+1]; rp++){
        		if (crcolidx[rp] == js) status[rp] = 1;
	    }
	  } 
	 }
	  /*for(int nd1 = 0; nd1 < mesh.elist[el]->nSideNode(sd); nd1++) {
	    i = mesh.elist[el]->SideNode(sd,nd1);
	    is = mesh.elist[el]->Node[i];
	    for(int nd2 = 0; nd2 < mesh.elist[el]->nSideNode(sd); nd2++){
		j = mesh.elist[el]->SideNode(sd, nd2);
		js = mesh.elist[el]->Node[j];
		for (int rp = crrowptr[is]; rp < crrowptr[is+1]; rp++){
        		if (crcolidx[rp] == js) status[rp] = 1;
		}
	    }
	  }*/
         }
    }
    
    int tnzero=0;
    for(i=0; i < nzero; i++)
	if(status[i]) tnzero++; 
   int *spatrowptr, *spatcolidx;
   spatrowptr = new int[sysdim + 1];
   spatcolidx = new int[tnzero];
   spatrowptr[0] = 0;
   j=0;
   for(i = 0; i < sysdim; i++)
   {
	int rp1 = crrowptr[i];
	int rp2 = crrowptr[i+1];
        k=0;
	for(int rp = rp1; rp < rp2; rp++)
	{
		if(status[rp]){ 
			k++;
			spatcolidx[j] = crcolidx[rp];
			j++;			
		} 
	}
	spatrowptr[i+1] = spatrowptr[i] + k;
   } 
   delete []status;
    int ia, ib, ka, kb, ja, jb, idx;
    int va_i, vb_i, v_i;
    int na = sysdim, ma = sysdim, va = spatrowptr[sysdim], col;
    int *rowptr = new int[sum(node_angN)+1];
    for(i = 0; i < sum(node_angN) + 1; i++)
    	rowptr[i] = 0;
    
    k=1; 
    for (ia = 0; ia < na; ia++)
    {
	for(j = 0; j < node_angN[ia]; j++)
	{
		for(i = spatrowptr[ia]; i < spatrowptr[ia+1]; i++)
		{
			col = spatcolidx[i];
			rowptr[k] += node_angN[col];
		 }

		k++;
        } 
    }
  
    for(i = 1; i < sum(node_angN)+1; i++)
	rowptr[i] += rowptr[i-1];
	
    
   int *colidx = new int[rowptr[sum(node_angN)]];
   int col_offset;
   k=0;
   for(ia = 0; ia < na; ia++)
   {
	for(j = 0; j < node_angN[ia]; j++)
	{	
		for(i = spatrowptr[ia]; i < spatrowptr[ia+1]; i++)
		{
			col = spatcolidx[i];

			col_offset = offset[col];

			/*if(ia == 1669)
				cout<<col<< "  "<<node_angN[col]<<"  "<<node_angN[1669]<<endl;*/

			for(int l = 0; l < node_angN[col]; l++)
			{
				colidx[k] = col_offset + l;
				k++;
			}
		}
	}
   }
   A2.Initialise(rowptr, colidx);
   b1.Initialise(rowptr, colidx);
   A2.Zero(); 
   b1.Zero();

   delete []rowptr;
   delete []colidx;
   delete []spatrowptr;
   delete []spatcolidx;

}

void findMaxLocalSphOrder(const Mesh &mesh, const IVector& sphOrder, const IVector& node_angN, const int el, int &maxSphOrder, int &maxAngN)
{
  int nodel = mesh.elist[el]->nNode();
  maxAngN = 0; maxSphOrder =0; 
  for(int i=0; i<nodel; i++)
  {
	int is = mesh.elist[el]->Node[i];
	maxSphOrder = MAX(sphOrder[is], maxSphOrder);
	maxAngN = MAX(node_angN[is], maxAngN);

  }

}

//void genmat_boundint_3D(const Mesh& mesh,  const IVector& sphOrder, const IVector& node_angN, const IVector& offset, const RDenseMatrix& ptsPlus, const CVector& wtsPlus, CCompRowMatrix &Aintsc, CCompRowMatrix &Aintss, CCompRowMatrix &Aintc, CCompRowMatrix& A2, CCompRowMatrix& b1)
void genmat_boundint_3D(const Mesh& mesh,  const IVector& sphOrder, const IVector& node_angN, const IVector& offset, const RDenseMatrix& ptsPlus, const RVector& wtsPlus, RDenseMatrix* &Ylm, RCompRowMatrix& A2, RCompRowMatrix& b1)

  /*      
     produces complete matrix of integrals in space and angle for boundary term
  */
{
  // S2Mesh must be 3D - should check for that...

   const int sysdim = mesh.nlen();       // dimensions are size of nodes.
   const int fullsysdim = sum(node_angN);     // full size of angles X space nodes
   
   A2.New(fullsysdim, fullsysdim);
   b1.New(fullsysdim, fullsysdim);
   
   initialiseA2b1(mesh, node_angN, offset, A2, b1);
      
   double ela_ij;
   int el, nodel, is, js;
   // now create matrix, by looping over elements that have a boundary
   for (el = 0; el < mesh.elen(); el++) {
	//cout<<"el: "<<el<<endl;
	if(!(el*100%mesh.elen()))
		cout<<el*100/mesh.elen() <<"% work done ..."<<endl;
        if(!mesh.elist[el]->HasBoundarySide ()) continue;

	nodel = mesh.elist[el]->nNode();
       
	// now determine the element integrals
	for(int sd = 0; sd <  mesh.elist[el]->nSide(); sd++)  {
	  // if sd is not a boundary side. skip 
	  if(!mesh.elist[el]->IsBoundarySide (sd)) continue;
            
	  // get boundary normal...
	  RVector nhat = mesh.ElDirectionCosine(el,sd);
	  
	  RCompRowMatrix  Angbintplus, Angbintminus;
	  int lmaxAngN, lmaxSphOrder;
          findMaxLocalSphOrder(mesh, sphOrder, node_angN, el, lmaxSphOrder, lmaxAngN);
	  //cout<<"maxes: "<<lmaxSphOrder<<"  "<<lmaxAngN<<endl;
	  BIntUnitSphere(lmaxAngN, lmaxAngN, lmaxSphOrder, lmaxSphOrder, ptsPlus, wtsPlus, nhat, Ylm, Angbintplus, Angbintminus);
          //cout<<"rmatsize: "<<Angbintminus.nRows()<<"  "<<Angbintminus.nCols()<<"   "<<Angbintplus.nRows()<<"  "<<Angbintplus.nCols()<<endl;

	  for (int i = 0; i < nodel; i++) {
	    if ((is = mesh.elist[el]->Node[i]) >= sysdim) continue;
	    for (int j = 0; j < nodel; j++) {
		if ((js = mesh.elist[el]->Node[j]) >= sysdim) continue;
		ela_ij = mesh.elist[el]->BndIntFFSide (i, j,sd);
		kronplus(is, js, node_angN, offset, ela_ij, Angbintplus, A2);
	 	kronplus(is, js, node_angN, offset, ela_ij, Angbintminus, b1);	
	    }
	  } 
	
        } // end loop on sides
   	
    
   } // end loop on elements
   A2.Shrink();
   b1.Shrink();
}

};// end MyDataContext class

// =========================================================================
// global parameters

SourceMode srctp = SRCMODE_NEUMANN;   // source type
ParamParser pp;
QMMesh qmmesh;
NodeList &nlist=qmmesh.nlist;
ElementList &elist=qmmesh.elist;
//double wts17[] = {-14.776315, 4.131531, 3.053487, -4.592317, 4.131531, 3.053487, 0.436314, -1.996058, -0.292499, -1.939062, -1.210250, 0.379318, -1.267245, -0.349494, 7.892854, 1.489065, 6.799725, 6.420396, 7.513526, 2.961523, 2.582194, 1.868393, 1.534250, -4.590531, 0.441119, -3.644677, -2.551546, 0.588395, -3.497401, -0.504735, 0.166289, -2.252024, -0.926841, -1.306169, -0.213039, -0.779565, -1.158893, -1.872696, -2.679047, 1.500590, -1.576812, -0.155850, -1.258085, -1.022607, 0.398355, 0.079628, -1.829717, 3.526964, -0.727481, 2.736109, 1.633874, -1.038863, 2.424728, 0.063373, -0.920597, 2.393454, 0.181638, 1.602599, 0.500364, -0.129743, 1.291218, 0.972493, 1.106781, -1.442977, 0.522757, -0.895899, -0.311875, 0.559702, -0.858954, -0.024322, 2.447745, 0.171931, 1.863723, 1.920962, 2.504984, 0.698714, 0.755953, 0.114692, 3.070372, -0.681341, 2.486350, 1.067690, 1.651712, 1.321340, -0.097319, 0.737318, -3.722966, -1.578274, -4.327823, -0.973417, -0.973235, -0.080292, -1.578086, 0.524559, -1.020065, 0.477729, -1.020065, 0.477729, -3.949549, -1.200000, -3.949549, -1.200000, 0.625000, 0.625000, 0.625000, 0.625000, 0.445251, 0.445251, 0.445251, 0.445251};
double wts17[] = {.0038282704949371616, .0038282704949371616, .0038282704949371616, .0038282704949371616, .0038282704949371616, .0038282704949371616, .0097937375124875125, .0097937375124875125, .0097937375124875125, .0097937375124875125, .0097937375124875125, .0097937375124875125, .0097937375124875125, .0097937375124875125, .0082117372831911110, .0082117372831911110, .0082117372831911110, .0082117372831911110, .0082117372831911110, .0082117372831911110, .0082117372831911110, .0082117372831911110, .0082117372831911110, .0082117372831911110, .0082117372831911110, .0082117372831911110, .0082117372831911110, .0082117372831911110, .0082117372831911110, .0082117372831911110, .0082117372831911110, .0082117372831911110, .0082117372831911110, .0082117372831911110, .0082117372831911110, .0082117372831911110, .0082117372831911110, .0082117372831911110, .0095954713360709628, .0095954713360709628, .0095954713360709628, .0095954713360709628, .0095954713360709628, .0095954713360709628, .0095954713360709628, .0095954713360709628, .0095954713360709628, .0095954713360709628, .0095954713360709628, .0095954713360709628, .0095954713360709628, .0095954713360709628, .0095954713360709628, .0095954713360709628, .0095954713360709628, .0095954713360709628, .0095954713360709628, .0095954713360709628, .0095954713360709628, .0095954713360709628, .0095954713360709628, .0095954713360709628, .0099428148911781033, .0099428148911781033, .0099428148911781033, .0099428148911781033, .0099428148911781033, .0099428148911781033, .0099428148911781033, .0099428148911781033, .0099428148911781033, .0099428148911781033, .0099428148911781033, .0099428148911781033, .0099428148911781033, .0099428148911781033, .0099428148911781033, .0099428148911781033, .0099428148911781033, .0099428148911781033, .0099428148911781033, .0099428148911781033, .0099428148911781033, .0099428148911781033, .0099428148911781033, .0099428148911781033, .0096949963616630283, .0096949963616630283, .0096949963616630283, .0096949963616630283, .0096949963616630283, .0096949963616630283, .0096949963616630283, .0096949963616630283, .0096949963616630283, .0096949963616630283, .0096949963616630283, .0096949963616630283, .0096949963616630283, .0096949963616630283, .0096949963616630283, .0096949963616630283, .0096949963616630283, .0096949963616630283, .0096949963616630283, .0096949963616630283, .0096949963616630283, .0096949963616630283, .0096949963616630283, .0096949963616630283}; 



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

// =========================================================================
// local prototypes

inline CVector matrixFreeCaller(const CVector& x, void * context);
void SelectSourceProfile (int &qtype, double &qwidth, SourceMode &srctp);
void SelectMeasurementProfile (ParamParser &pp, int &mtype, double &mwidth);
void genmat_source_3D(const IVector& sphOrder, const IVector& node_angN, const IVector& offset, RCompRowMatrix* & Source, const Mesh& mesh,  const int Nsource,const int ns, const RVector& dirVec, const bool is_isotropic, const RDenseMatrix& pts, const RVector& wts, const RCompRowMatrix& b1, RDenseMatrix* &Ylm);
void genmat_sourcevalvector_cos_3D(const IVector& sphOrder, const IVector& node_angN, const IVector& offset, RCompRowMatrix& Svec, const Mesh& mesh, const int Nsource, const RVector& dirVec, const bool is_isotropic, const RDenseMatrix& pts, const RVector& wts, RDenseMatrix* &Ylm);
void genmat_toastsourcevalvector3D_cos(const IVector& sphOrder, const IVector& node_angN, const IVector& offset, RCompRowMatrix& Svec, const Mesh& mesh, const RCompRowMatrix qvec, const int iq, const RVector& dirVec, const bool is_isotropic, const RDenseMatrix& pts, const RVector& wts, RDenseMatrix* &Ylm);
void genmat_toastsource3D(const IVector& sphOrder, RCompRowMatrix* & Source, const Mesh& mesh, const RCompRowMatrix qvec, const int ns, const RVector& dirVec, const bool is_isotropic, const RDenseMatrix& pts, const RVector& wts, const RCompRowMatrix& b1, RDenseMatrix* &Ylm);
void WriteData (const RVector &data, char *fname);
void WriteDataBlock (const QMMesh &mesh, const RVector &data, char *fname);
void OpenNIM (const char *nimname, const char *meshname, int size);
void WriteNIM (const char *nimname, const RVector &img, int size, int no);
bool ReadNim (char *nimname, RVector &img);
void WritePGM (const RVector &img, const IVector &gdim, char *fname);
void WritePPM (const RVector &img, const IVector &gdim,
double *scalemin, double *scalemax, char *fname);
CVector getDiag(void * context);

// error handler for FE library routines *************************************
void LocalErrorhandler (char *msg)
{
    cerr << "\nread_toastbemmesh (PID " << getpid() << ")\n" << msg << endl << flush;
    cerr << "Aborted.\n";
    //    logfile << msg << endl << "Aborted." << endl;
    exit (1);
}

void genmat_toastsource3D(const IVector& sphOrder, const IVector& node_angN, const IVector& offset, RCompRowMatrix* & Source, const Mesh& mesh, const RCompRowMatrix qvec, const int ns, const RVector& dirVec, const bool is_isotropic, const RDenseMatrix& pts, const RVector& wts, const RCompRowMatrix& b1, RDenseMatrix* &Ylm)
   /*      
    Function generates the source values vector for FEM of the radiative 
    transfer equation
  */
{
   int sysdim = mesh.nlen();       // dimensions are size of nodes.
   int fullsysdim = sum(node_angN);   // full size of angles X space nodes
   RCompRowMatrix Svec;
   for (int i = 0; i < ns; i++) {
     Source[i].New(fullsysdim,1);
     genmat_toastsourcevalvector3D_cos(sphOrder, node_angN, offset, Svec, mesh, qvec,i, dirVec, is_isotropic,  pts, wts, Ylm);
     b1.AB(Svec, Source[i]);
     cout<<"Prepared source number "<<i<<endl;
   }
}


void genmat_toastsourcevalvector3D_cos(const IVector& sphOrder, const IVector& node_angN, const IVector& offset, RCompRowMatrix& Svec, const Mesh& mesh, const RCompRowMatrix qvec, const int iq, const RVector& dirVec, const bool is_isotropic,  const RDenseMatrix& pts, const RVector& wts, RDenseMatrix* &Ylm)
  /*      
    Function generates the source values vector for FEM of the radiative 
    transfer equation
  */
{
   int el, nodel, i, j, k,is, js;
   int sysdim = mesh.nlen();       // dimensions are size of nodes.
   int fullsysdim = sum(node_angN);   // full size of angles X space nodes
   
   int *rowptr, *colidx;
   rowptr = new int[fullsysdim+1];
   colidx = new int[fullsysdim];
   rowptr[0] = 0;
   for(int i=1;  i <= fullsysdim; i++)
	rowptr[i] = i;
   for(int i=0; i<fullsysdim; i++)
	colidx[i] = 0;
   Svec.New (fullsysdim,1);
   Svec.Initialise(rowptr, colidx);

   int row_offset = 0;
   RDenseMatrix Ystarlm;
   RDenseMatrix dirMat(1, 3);

   // now create vector, by looping over elements that have a boundary
   for (int jq = qvec.rowptr[iq]; jq < qvec.rowptr[iq+1]; jq++) {
   	int Nsource = qvec.colidx[jq];
   	double sweight = norm(qvec.Get(iq,Nsource)); 
   	for (el = 0; el < mesh.elen(); el++) {
        	if(!mesh.elist[el]->IsNode(Nsource)) continue; // source not in this el
        	if(!mesh.elist[el]->HasBoundarySide()) continue;
		for(int sd = 0; sd <  mesh.elist[el]->nSide(); sd++) {
	  		// if sd is not a boundary side. skip 
	  		if(!mesh.elist[el]->IsBoundarySide (sd)) continue;
 			RVector nhat = mesh.ElDirectionCosine(el,sd);
			dirMat(0, 0) = dirVec[0]; dirMat(0, 1) = dirVec[1]; dirMat(0, 2) = dirVec[2]; 
	  		for(int nd = 0; nd < mesh.elist[el]->nSideNode(sd); nd++) {
	    			is = mesh.elist[el]->SideNode(sd,nd);
	    			js = mesh.elist[el]->Node[is];
				double ela_i =  mesh.elist[el]->IntF(is);
	    			if(is_isotropic)
				{
					RDenseMatrix Angsvec(node_angN[js], 1);
					for(int l=0; l<= sphOrder[js]; l++){
						int indl = getPos(l, -1*l);
						//Ystarlm = sphericalHarmonics(l, pts.nRows(), pts);
						for(int m=-l; m<=l; m++){
							int indm = l + m;
							for(int i=0; i < pts.nRows(); i++)
							{
								double sdotn = nhat[0]*pts.Get(i, 0) + nhat[1]*pts.Get(i, 1) + nhat[2]*pts.Get(i, 2);
								if(sdotn<0)
									Angsvec(indl+indm, 0) += Ylm[l](indm, i)*(-sdotn)*wts[i]*4*M_PI;
							}
						}
					}
					for(int i = 0; i < node_angN[js]; i++)
						Svec(offset[js] + i, 0) += Angsvec(i, 0)*sweight*ela_i;

				}
				else{
					RDenseMatrix Angsvec(node_angN[js], 1);
					RDenseMatrix *Ystarlm;
					Ystarlm = new RDenseMatrix[sphOrder[js]+1];
					for(int l=0; l<= sphOrder[js]; l++)
						Ystarlm[l].New(2*l+1, 1);
					sphericalHarmonics(sphOrder[js], 1, dirMat, Ystarlm);
					for(int l=0; l<= sphOrder[js]; l++){
						int indl = getPos(l, -1*l);
						//Ystarlm = sphericalHarmonics(l, 1, dirMat);
						for(int m=-l; m<=l; m++){
							int indm = l + m;
							Angsvec(indl+indm, 0) = Ystarlm[l](indm, 0);
						}
					}
					for(int i = 0; i < node_angN[js]; i++)
						Svec(offset[js] + i, 0) += Angsvec(i, 0)*sweight;

					delete []Ystarlm;
				}
	    		}

		} // end loop on element sides

   	} // end loop on elements
  }
   delete []rowptr;
   delete []colidx;
}  


void genmat_source_3D(const IVector& sphOrder, const IVector& node_angN, const IVector& offset, RCompRowMatrix* & Source, const Mesh& mesh,  const int* Nsource, const int ns, const RVector& dirVec, const bool is_isotropic, const RDenseMatrix& pts, const RVector& wts, const RCompRowMatrix& b1, RDenseMatrix* &Ylm)
  /*      
    Function generates the source values vector for FEM of the radiative 
    transfer equation
  */
{
   int sysdim = mesh.nlen();         // dimensions are size of nodes.
   int fullsysdim = sum(node_angN);       // full size of angles X space nodes

   RCompRowMatrix Svec;
   for (int i = 0; i < ns; i++) {
     Source[i].New(fullsysdim,1);
     genmat_sourcevalvector_cos_3D(sphOrder, node_angN, offset, Svec, mesh, Nsource[i], dirVec, is_isotropic, pts, wts, Ylm);
     b1.AB(Svec, Source[i]);
   }
}
void genmat_sourcevalvector_cos_3D(const IVector& sphOrder, const IVector& node_angN, const IVector& offset, RCompRowMatrix& Svec, const Mesh& mesh, const int Nsource, const RVector& dirVec, const bool is_isotropic, const RDenseMatrix& pts, const RVector& wts, RDenseMatrix* &Ylm)
  /*      
    Function generates the source values vector for FEM of the radiative 
    transfer equation
  */
{
   int el, nodel, i, j, k,is, js;
   int sysdim = mesh.nlen();       // dimensions are size of nodes.
   int fullsysdim = sum(node_angN);   // full size of angles X space nodes
   
   int *rowptr, *colidx;
   rowptr = new int[fullsysdim+1];
   colidx = new int[fullsysdim];
   rowptr[0] = 0;
   for(int i=1;  i <= fullsysdim; i++)
	rowptr[i] = i;
   for(int i=0; i<fullsysdim; i++)
	colidx[i] = 0;
   Svec.New (fullsysdim,1);
   Svec.Initialise(rowptr, colidx);

   int row_offset = 0;
   RDenseMatrix dirMat(1, 3);
   RDenseMatrix Ystarlm;

   // now create vector, by looping over elements that have a boundary
   for (el = 0; el < mesh.elen(); el++) {
        if(!mesh.elist[el]->IsNode(Nsource)) continue; // source not in this el
        if(!mesh.elist[el]->HasBoundarySide()) continue;
	for(int sd = 0; sd <  mesh.elist[el]->nSide(); sd++) {
	  // if sd is not a boundary side. skip 
	  if(!mesh.elist[el]->IsBoundarySide (sd)) continue;
	  
	  RVector nhat = mesh.ElDirectionCosine(el,sd);
	  dirMat(0, 0) = dirVec[0]; dirMat(0, 1) = dirVec[1]; dirMat(0, 2) = dirVec[2]; 
	 //cout<<"direction vector: "<<dirMat<<endl;
	  for(int nd = 0; nd < mesh.elist[el]->nSideNode(sd); nd++) {
	    is = mesh.elist[el]->SideNode(sd,nd);
	    js = mesh.elist[el]->Node[is];
	   if(js != Nsource) continue;
	    double ela_i =  mesh.elist[el]->IntF(is);
	    if(is_isotropic)
	    {
			RDenseMatrix Angsvec(node_angN[js], 1);
			for(int l=0; l<= sphOrder[js]; l++){
				int indl = getPos(l, -1*l);
				//Ystarlm = sphericalHarmonics(l, pts.nRows(), pts);
				for(int m=-l; m<=l; m++){
					int indm = l + m;
					for(int i=0; i < pts.nRows(); i++)
					{
						double sdotn = nhat[0]*pts.Get(i, 0) + nhat[1]*pts.Get(i, 1) + nhat[2]*pts.Get(i, 2);
						if(sdotn<0)
							Angsvec(indl+indm, 0) += Ylm[l](indm, i)*(-sdotn)*wts[i]*4*M_PI;
					}
				}
			}
			for(int i = 0; i < node_angN[js]; i++)
				Svec(offset[js] + i, 0) = Angsvec(i, 0)*ela_i;
			
		}
	    	else{
			RDenseMatrix Angsvec(node_angN[js], 1);
			RDenseMatrix *Ystarlm;
			Ystarlm = new RDenseMatrix[sphOrder[js]+1];
			for(int l=0; l<= sphOrder[js]; l++)
			 Ystarlm[l].New(2*l+1, 1);
			sphericalHarmonics(sphOrder[js], 1, dirMat, Ystarlm);
			for(int l=0; l<= sphOrder[js]; l++){
				int indl = getPos(l, -1*l);
				//Ystarlm = sphericalHarmonics(l, 1, dirMat);
				for(int m=-l; m<=l; m++){
					int indm = l + m;
					Angsvec(indl+indm, 0) = Ystarlm[l](indm, 0);
					cout<<"l: "<<l<<" m: "<<m<<" "<<Angsvec(indl+indm, 0)<<endl;
				}
			}
			for(int i = 0; i < node_angN[js]; i++)
				Svec(offset[js] + i, 0) = Angsvec(i, 0)*ela_i;
			delete []Ystarlm;

			/*RDenseMatrix Angsvec(node_angN[js], 1);
			for(int l=0; l<= sphOrder[js]; l++){
				int indl = getPos(l, -1*l);
				//Ystarlm = sphericalHarmonics(l, pts.nRows(), pts);
				for(int m=-l; m<=l; m++){
					int indm = l + m;
					for(int i=0; i < pts.nRows(); i++)
					{
						double sdotn = dirVec[0]*pts.Get(i, 0) + dirVec[1]*pts.Get(i, 1) + dirVec[2]*pts.Get(i, 2);
						Angsvec(indl+indm, 0) += Ylm[l](indm, i)*exp(-(1-sdotn)*(1-sdotn)/(2*0.01*0.01))*wts[i]*4*M_PI;
					}
					cout<<"l: "<<l<<" m: "<<m<<" "<<Angsvec(indl+indm, 0)<<endl;
				}
			}
				
			for(int i = 0; i < node_angN[js]; i++)
				Svec(offset[js] + i, 0) += Angsvec(i, 0)*ela_i;*/

		}
	}
    } 


   } // end loop on elements
   delete []rowptr;
   delete []colidx;
}
void genmat_toastintsourcevalvector_cos_3D(const IVector& sphOrder, const IVector& node_angN, const IVector& offset, RVector& Svec, const Mesh& mesh, const RCompRowMatrix qvec, const int iq, const RVector& dirVec, const bool is_isotropic,  const RDenseMatrix& pts, const RVector& wts, RDenseMatrix* &Ylm)
  /*      
    Function generates the source values vector for FEM of the radiative 
    transfer equation
  */
{
   int el, nodel, i, j, k,is, js;
   int sysdim = mesh.nlen();       // dimensions are size of nodes.
   int fullsysdim = sum(node_angN);   // full size of angles X space nodes
   
   int row_offset = 0;
   RDenseMatrix Ystarlm;
   RDenseMatrix dirMat(1, 3);

   // now create vector, by looping over elements that have a boundary
   for (int jq = qvec.rowptr[iq]; jq < qvec.rowptr[iq+1]; jq++) {
   	int Nsource = qvec.colidx[jq];
   	double sweight = norm(qvec.Get(iq,Nsource)); 
   	for (el = 0; el < mesh.elen(); el++) {
		for(int sd = 0; sd <  mesh.elist[el]->nSide(); sd++) {
	  		// if sd is not a boundary side. skip 
 			RVector nhat = mesh.ElDirectionCosine(el,sd);
			dirMat(0, 0) = -1*nhat[0]; dirMat(0, 1) = -1*nhat[1]; dirMat(0, 2) = -1*nhat[2]; 
	  		for(int nd = 0; nd < mesh.elist[el]->nSideNode(sd); nd++) {
	    			is = mesh.elist[el]->SideNode(sd,nd);
	    			js = mesh.elist[el]->Node[is];
				double ela_i =  mesh.elist[el]->IntF(is);
	    			if(is_isotropic)
				{
					RDenseMatrix Angsvec(node_angN[js], 1);
					for(int l=0; l<= sphOrder[js]; l++){
						int indl = getPos(l, -1*l);
						//Ystarlm = sphericalHarmonics(l, pts.nRows(), pts);
						for(int m=-l; m<=l; m++){
							int indm = l + m;
							for(int i=0; i < pts.nRows(); i++)
							{
								double sdotn = nhat[0]*pts.Get(i, 0) + nhat[1]*pts.Get(i, 1) + nhat[2]*pts.Get(i, 2);
								if(sdotn<0)
									Angsvec(indl+indm, 0) += Ylm[l](indm, i)*(-sdotn)*wts[i]*4*M_PI;
							}
						}
					}
					for(int i = 0; i < node_angN[js]; i++)
						Svec[offset[js] + i] += Angsvec(i, 0)*sweight*ela_i;

				}
				else{
					RDenseMatrix Angsvec(node_angN[js], 1);
					RDenseMatrix *Ystarlm;
					Ystarlm = new RDenseMatrix[sphOrder[js]+1];
					for(int l=0; l<= sphOrder[js]; l++)
						Ystarlm[l].New(2*l+1, 1);
					sphericalHarmonics(sphOrder[js], 1, dirMat, Ystarlm);
					for(int l=0; l<= sphOrder[js]; l++){
						int indl = getPos(l, -1*l);
						//Ystarlm = sphericalHarmonics(l, 1, dirMat);
						for(int m=-l; m<=l; m++){
							int indm = l + m;
							Angsvec(indl+indm, 0) = Ystarlm[l](indm, 0);
						}
					}
					for(int i = 0; i < node_angN[js]; i++)
						Svec[offset[js] + i] += Angsvec(i, 0)*sweight;

					delete []Ystarlm;
				}
	    		}

		} // end loop on element sides

   	} // end loop on elements
  }
}  

void genmat_intsourcevalvector_cos_3D(const IVector& sphOrder, const IVector& node_angN, const IVector& offset, RVector& Svec, const Mesh& mesh, const int Nsource, const RVector& dirVec, const bool is_isotropic, const RDenseMatrix& pts, const RVector& wts, RDenseMatrix* &Ylm)
  /*      
    Function generates the source values vector for FEM of the radiative 
    transfer equation
  */
{
   int el, nodel, i, j, k,is, js;
   int sysdim = mesh.nlen();       // dimensions are size of nodes.
   int fullsysdim = sum(node_angN);   // full size of angles X space nodes
   
   int row_offset = 0;
   RDenseMatrix dirMat(1, 3);
   RDenseMatrix Ystarlm;

   // now create vector, by looping over elements that have a boundary
   for (el = 0; el < mesh.elen(); el++) {
        if(!mesh.elist[el]->IsNode(Nsource)) continue; // source not in this el
	for(int sd = 0; sd <  mesh.elist[el]->nSide(); sd++) {
	  RVector nhat = mesh.ElDirectionCosine(el,sd);
	  dirMat(0, 0) = dirVec[0]; dirMat(0, 1) = dirVec[1]; dirMat(0, 2) = dirVec[2]; 
	  //cout<<"direction vector: "<<dirMat<<endl;
	  for(int nd = 0; nd < mesh.elist[el]->nSideNode(sd); nd++) {
	    is = mesh.elist[el]->SideNode(sd,nd);
	    js = mesh.elist[el]->Node[is];
	   if(js != Nsource) continue;
	    double ela_i =  mesh.elist[el]->IntF(is);
	    if(is_isotropic)
	    {
			RDenseMatrix Angsvec(node_angN[js], 1);
			for(int l=0; l<= sphOrder[js]; l++){
				int indl = getPos(l, -1*l);
				//Ystarlm = sphericalHarmonics(l, pts.nRows(), pts);
				for(int m=-l; m<=l; m++){
					int indm = l + m;
					for(int i=0; i < pts.nRows(); i++)
					{
						double sdotn = nhat[0]*pts.Get(i, 0) + nhat[1]*pts.Get(i, 1) + nhat[2]*pts.Get(i, 2);
						if(sdotn<0)
							Angsvec(indl+indm, 0) += Ylm[l](indm, i)*(-sdotn)*wts[i]*4*M_PI;
					}
				}
			}
			for(int i = 0; i < node_angN[js]; i++)
				Svec[offset[js] + i] += Angsvec(i, 0)*ela_i;
			
		}
	    	else{
			RDenseMatrix Angsvec(node_angN[js], 1);
			RDenseMatrix *Ystarlm;
			Ystarlm = new RDenseMatrix[sphOrder[js]+1];
			for(int l=0; l<= sphOrder[js]; l++)
			 Ystarlm[l].New(2*l+1, 1);
			sphericalHarmonics(sphOrder[js], 1, dirMat, Ystarlm);
			for(int l=0; l<= sphOrder[js]; l++){
				int indl = getPos(l, -1*l);
				//Ystarlm = sphericalHarmonics(l, 1, dirMat);
				for(int m=-l; m<=l; m++){
					int indm = l + m;
					Angsvec(indl+indm, 0) = Ystarlm[l](indm, 0)*(4*M_PI)/(2*l+1);
					cout<<"l: "<<l<<" m: "<<m<<" "<<Angsvec(indl+indm, 0)<<endl;

				}
			}
			for(int i = 0; i < node_angN[js]; i++)
				Svec[offset[js] + i] += Angsvec(i, 0);
			delete []Ystarlm;

		
		}

		}
	    }



   } // end loop on elements
}
// main routine **************************************************************

int main (int argc, char *argv[])
{
    char cbuf[200];
    int el;
    //xASSERT(argc >= 3, Current implementation only when QM file is specified);
    //    logfile << "Reading mesh" << endl;
    cout << "Reading mesh " << argv[1] << endl;
    ifstream ifs;
    ifs.open (argv[1]);
    xASSERT (ifs.is_open(), Mesh file not found.);
    ifs >> qmmesh;
    xASSERT (ifs.good(), Problem reading mesh.);
    ifs.close ();
    cout << "* " << qmmesh.elen() << " elements, " << qmmesh.nlen()
	 << " nodes\n";
    int dimension = nlist[0].Dim();
    for (int i = 1; i < nlist.Len(); i++)
	xASSERT(nlist[i].Dim() == dimension, Inconsistent node dimensions.);
    xASSERT(dimension >= 2 && dimension <= 3, Mesh dimension must be 2 or 3.);
    qmmesh.Setup();

    // set up angular "mesh"

    //******** optical parameters

    //****** sources.
    //****** should read in QM file, but just set something for now.
    cout << "Forming the source\n";
    int ns = 1, nM;
    int *Nsource = new int [ns];
    int    qprof, mprof;   // source/measurement profile (0=Gaussian, 1=Cosine)
    double qwidth, mwidth; // source/measurement support radius [mm]
    RCompRowMatrix qvec, mvec;
    char file_extn[200];

    cin>>file_extn;
    cout<<"File name prefix: "<<file_extn<<endl;
    //    Measurement datatype;
    if(argc < 3) {
         cin >> ns;
	 cout<< "Number of sources: "<<ns<<endl;
         Nsource = new int[ns];
	 cout << "source node coords: "<<endl;
	 for(int i=0; i<ns; i++){
		cin >>Nsource[i];
     	 cout<<nlist[Nsource[i]]<<"   "<<endl;}
    }
    else {
    cout << "QM file " << argv[2] << endl;
    ifs.open (argv[2]);
    xASSERT (ifs.is_open(), QM file not found.);
    qmmesh.LoadQM (ifs);
    xASSERT (ifs.good(), Problem reading QM.);
    ifs.close ();
    ns = qmmesh.nQ;
    nM = qmmesh.nM;
    cout << ns << " sources\n";
    SelectSourceProfile (qprof, qwidth, srctp);
    SelectMeasurementProfile (pp, mprof, mwidth);
    // build the source vectors
    qvec.New (ns, qmmesh.nlen());
    for (int i = 0; i < ns; i++) {
	RVector q(qmmesh.nlen());
	switch (qprof) {
	case 0:
	    //SetReal (q, QVec_Point (qmmesh, qmmesh.Q[i], srctp));
	    q = QVec_Point (qmmesh, qmmesh.Q[i], srctp);
	    break;
	case 1:
	    //SetReal (q, QVec_Gaussian (qmmesh, qmmesh.Q[i], qwidth, srctp));
	    q = QVec_Gaussian (qmmesh, qmmesh.Q[i], qwidth, srctp);
	    break;
	case 2:
	    //SetReal (q, QVec_Cosine (qmmesh, qmmesh.Q[i], qwidth, srctp));
	    q = QVec_Cosine (qmmesh, qmmesh.Q[i], qwidth, srctp);
	    break;
	}
	qvec.SetRow (i, q);
    }
    cout << "Sources set "<<endl;
    // build the measurement vectors
    mvec.New (nM, qmmesh.nlen());
    LOGOUT1_INIT_PROGRESSBAR ("Meas. vectors", 50, nM);
    for (int i = 0; i < nM; i++) {
	RVector m(qmmesh.nlen());
	switch (mprof) {
	case 0:
	    //SetReal (m, QVec_Point (qmmesh, qmmesh.M[i], SRCMODE_NEUMANN));
	    m = QVec_Point (qmmesh, qmmesh.M[i], SRCMODE_NEUMANN);
	    break;
	case 1:
	    //SetReal (m, QVec_Gaussian (qmmesh, qmmesh.M[i], mwidth, SRCMODE_NEUMANN));
	    m = QVec_Gaussian (qmmesh, qmmesh.M[i], mwidth, SRCMODE_NEUMANN);
	    break;
	case 2:
	    //SetReal (m, QVec_Cosine (qmmesh, qmmesh.M[i], mwidth, SRCMODE_NEUMANN));
	    m = QVec_Cosine (qmmesh, qmmesh.M[i], mwidth, SRCMODE_NEUMANN);
	    break;
	}
	for (int j = 0; j < qmmesh.nlen(); j++) 
	  m[j] *= qmmesh.plist[j].C2A();
	mvec.SetRow (i, m);
    }
    }
    //***** parameters 
    double freq = 0, g;
    int phaseExpOrder, is_iso;
    bool is_isotropic;
    RVector dirVec(3);
    cin>>g;
    cout<< "value of g (WARNING !! This g is considered constant throughout the domain): "<< g<<endl;
    cin >> freq;
    cout << "value for frequency (MHz) : "<<freq<<endl;
    cin>> is_iso;
    is_isotropic = is_iso>0 ? true : false;
    cout << "The source is isotropic ? "<<is_isotropic<<endl;
    cin>> dirVec[0]; cin>> dirVec[1]; cin>>dirVec[2];
    dirVec = dirVec*1/(double)length(dirVec); // normalize it just in case
    cout<< "The direction vector for the source (if not isotropic) is "<<dirVec<<endl; 

    double w = freq * 2.0*M_PI*1e-6;
    double c = 0.3;
     int numpts;
     numpts = 110;
     RDenseMatrix pts(numpts, 3); 
     RVector wts(numpts);
     for(int i=0; i < numpts; i++){
	wts[i] = wts17[i]; 
	for(int j=0; j < 3; j++)
		pts(i, j) = pts17[i][j];
	}	

    RVector muscat(qmmesh.elen());
    RVector muabs(qmmesh.elen());
    RVector ref(qmmesh.elen());

    cout<< "Reading mua, mus and refractive index values for all the elements from file: "<<endl;
    for(el = 0; el <  qmmesh.elen(); el++){ // assign something to parameters
      cin >> muabs[el];
      cin >> muscat[el]; 
      cin >> ref[el]; 
    }
    cout << "mua " << muabs[0] << " mus " << muscat[0] <<" ref "<<ref[0]<<endl;
    // "smoothing parameter delta of streamline diffusion modification.
    // Just set to constant for now.
    IVector sphOrder(qmmesh.nlen());
    cout<<"Reading spherical harmonic expansion order for each node: "<<endl;
    for(int i=0; i < qmmesh.nlen(); i++)
	cin >> sphOrder[i]; 

    RVector delta(qmmesh.elen());
    double min = 1e20, max = -1e20;
    cout << "delta " ;
    for(el = 0; el <  qmmesh.elen(); el++){ // assign something to delta 
#ifdef VARYDELTA 
	 // this requires Nsource to be set...
         Point centre = Point3D(0, 0, 0);
	 int nNode =  elist[el]->nNode();
	 for(int i=0; i < nNode; i++) 
		centre += nlist[elist[el]->Node[i]];
	 centre = centre/nNode;
	 double sval =  elist[el]->Size()/((muscat[el]));
	  // double sval = 1/(5.0*(muabs[el] + muscat[el]));
	  // double dist = l2norm(centre - nlist[1669]);
	 //  double rval = 5.0*sval/dist; 
          // delta[el] = MIN(sval,rval);
	 delta[el] = sval;
#else
	 min = max = delta[0];
#endif 
	 //	 cout << delta[el] << " ";
    }
    cout << "min : " << vmin(delta) << "max : " << vmax(delta) << endl;
   MyDataContext ctxt(qmmesh, sphOrder, delta, muabs, muscat, ref, g, &phaseFunc, w, c, pts, wts);

     
    CVector proj(ns*nM); // projection data
    RCompRowMatrix* Source;
    RVector b2(sum(ctxt.node_angN));
    if( !(Source = new  RCompRowMatrix [ns]))
          cerr << "Memory Allocation error Source = new  RCompRowMatrix\n";
    cout<<"calling computing the source vector ..."<<endl;
   if(argc<3)	
     {
      genmat_source_3D(ctxt.nodal_sphOrder, ctxt.node_angN, ctxt.offset, Source, qmmesh, Nsource, ns, dirVec, is_isotropic, pts, wts, ctxt.b1, ctxt.Ylm);
      genmat_intsourcevalvector_cos_3D(ctxt.nodal_sphOrder, ctxt.node_angN, ctxt.offset, b2, qmmesh, Nsource[0], dirVec, is_isotropic, pts, wts, ctxt.Ylm);
     }
    else 
    {
      genmat_toastsource3D(ctxt.nodal_sphOrder,  ctxt.node_angN, ctxt.offset, Source, qmmesh, qvec, ns, dirVec, is_isotropic, pts, wts, ctxt.b1, ctxt.Ylm);
      genmat_toastintsourcevalvector_cos_3D(ctxt.nodal_sphOrder,  ctxt.node_angN, ctxt.offset, b2, qmmesh, qvec, 0, dirVec, is_isotropic, pts, wts, ctxt.Ylm);
    }

    cout << "calculating the radiance\n";
    int sysdim = qmmesh.nlen();
    int fullsysdim = sum(ctxt.node_angN);
    CVector RHS(fullsysdim);
    CVector * Phi = new CVector [ns];
    CVector * Phisum = new CVector [ns];
    CPrecon_Diag * AACP = new  CPrecon_Diag;
    CVector idiag = getDiag(&ctxt);
    AACP->ResetFromDiagonal(idiag);
    double tol = 1e-9;
    ofstream osPhi("Phi.sol");
    ofstream osPhisum("Phisum.sol");
    ofstream osRHS("RHS.vec");
    double res;
    int iter;
    clock_t start = clock();
    int row_offset = 0;
    for (int j = 0; j < ns ; j++) {   
      cout << "start source " << j << endl;

      Phi[j].New(fullsysdim);
      Phisum[j].New(sysdim);
       for(int i = 0; i < fullsysdim; i++)
	RHS[i] = Source[j].Get(i,0) + b2[i];

      GMRES(&matrixFreeCaller, &ctxt, RHS, Phi[j], tol, AACP, 100);
      osRHS << RHS << endl;
      cout << "finished source " << j << endl;
      osPhi << "Phi " << j << "\n" << Phi[j] << endl;
      for (int k = 0; k < sysdim; k++)
      {
	  Phisum[j][k] += Phi[j][ctxt.offset[k]]*sqrt(4*M_PI);
      }
      
      osPhisum << "Phisum " << j << "\n" << Phisum[j] << endl;
      for (int im = 0; im < nM; im++) {
	for(int in = 0; in < sysdim; in++)
	  proj[j*nM + im] += Phisum[j][in]*mvec.Get(im,in) ;
      }
    }
    clock_t end = clock();
    osPhi.close();
    osPhisum.close();
    osRHS.close();
    char flnmod[300], farg[300], ftime[300];
    strcpy(flnmod, file_extn); strcpy(farg, file_extn);strcpy(ftime, file_extn);
    strcat(flnmod, "_lnmod.nim"); strcat(farg, "_arg.nim"); strcat(ftime, "_time.txt");
    OpenNIM (flnmod, argv[1], sysdim);
    OpenNIM (farg, argv[1], sysdim);
    for (int i = 0; i < ns; i++) {
	   WriteNIM (flnmod, LogMod(Phisum[i]), sysdim, i);
	   WriteNIM (farg, Arg(Phisum[i]), sysdim, i);
    }
    cout << "  Log Mod field written to "<< flnmod << endl;
    cout << "  Arg field written to "<<farg << endl;


    // output data files


	WriteData (LogMod(proj), "fmod_hr.fem");
	cout<<"fmod_hr.fem written ..."<<endl;
        WriteData (Arg(proj), "farg_hr.fem");
	cout<<"farg_hr.fem written ..."<<endl;
        WriteDataBlock (qmmesh, LogMod(proj), "fmod_hr.dat");
	cout<<"fmod_hr.dat written ..."<<endl;
	WriteDataBlock (qmmesh, Arg(proj), "farg_hr.dat");
	cout<<"farg_hr.dat written ..."<<endl;

    FILE *fid;
    fid = fopen(ftime, "w");
    fprintf(fid, "Time taken by solver: %f\n", (double)(end-start)/CLOCKS_PER_SEC);
    fclose(fid);
    cout<<"The solver took "<<(double)(end-start)/CLOCKS_PER_SEC<<" seconds"<<endl;

	delete []Phi;
	delete []Phisum;
	delete []Nsource;
	delete []Source;

}

inline void RCAx(const RCompRowMatrix &A, const CVector& x, CVector &res)
{
    dASSERT_2PRM(x.Dim() == A.nCols(),
    "Parameter 1 invalid size (expected %d, actual %d)",
    A.nCols(), x.Dim());
    if (res.Dim() != A.nRows()) res.New(A.nRows());

    int r, i, i2;
    complex br;
    const double *aval;
    aval = A.ValPtr();

    for (r = i = 0; r < A.nRows();) {
	i2 = A.rowptr[r+1];
	for (br = complex(0, 0); i < i2; i++)
	    br += x[A.colidx[i]]*aval[i];
	res[r++] = br;
    }

}

inline CVector matrixFreeCaller(const CVector& x, void * context)
{
	MyDataContext *ctxt = (MyDataContext*) context;
	int dof = sum(ctxt->node_angN);
	CVector result(dof);
	int spatN = ctxt->spatN;
	int maxAngN = ctxt->maxAngN;

	memcpy (ctxt->xmatval, x.data_buffer(), dof*sizeof(complex));

	int i, j, k, m, ra, ra1, ra2, rb, rb1, rb2;
    	int nr = ctxt->Xmat.nRows();
    	int nc = ctxt->Aint.nCols();
    	
	ctxt->Aintx.Zero(); ctxt->Aintscx.Zero(); ctxt->Aintssx.Zero(); ctxt->Aintcx.Zero();
	ctxt->apu1x.Zero(); ctxt->apu1scx.Zero(); ctxt->apu1ssx.Zero(); ctxt->apu1cx.Zero();
	ctxt->Aintscscx.Zero(); ctxt->Aintscssx.Zero(); ctxt->Aintssssx.Zero(); ctxt->Aintsccx.Zero();
	ctxt->Aintsscx.Zero(); ctxt->Aintccx.Zero();
	
	complex *aintxval = ctxt->Aintx.data_buffer(); complex *aintscxval = ctxt->Aintscx.data_buffer(); complex *aintssxval = ctxt->Aintssx.data_buffer();
	complex *aintcxval = ctxt->Aintcx.data_buffer(); complex *apu1xval = ctxt->apu1x.data_buffer(); complex *apu1scxval = ctxt->apu1scx.data_buffer();  
	complex *apu1ssxval = ctxt->apu1ssx.data_buffer(); complex *apu1cxval = ctxt->apu1cx.data_buffer(); 
	complex *aintscscxval = ctxt->Aintscscx.data_buffer(); complex *aintscssxval = ctxt->Aintscssx.data_buffer();  
	complex *aintssssxval = ctxt->Aintssssx.data_buffer(); complex *aintsccxval = ctxt->Aintsccx.data_buffer();  
	complex *aintsscxval = ctxt->Aintsscx.data_buffer();  complex *aintccxval = ctxt->Aintccx.data_buffer();  


	complex xval;
    	for (i = 0; i < nr; i++) {
    		ra1 = ctxt->Xmat.rowptr[i];
		ra2 = ctxt->Xmat.rowptr[i+1];
		for (ra = ra1; ra < ra2; ra++) {
			j = ctxt->Xmat.colidx[ra];
			xval = ctxt->xmatval[ra];

	    		rb1 = ctxt->Aint.rowptr[j];
	    		rb2 = ctxt->Aint.rowptr[j+1];
	    		for (rb = rb1; rb < rb2; rb++) {
	        		k = ctxt->Aint.colidx[rb];
				aintxval[i*maxAngN + k] += xval*ctxt->aintval[rb];		
	    		}

			rb1 = ctxt->Aintsc.rowptr[j];
	    		rb2 = ctxt->Aintsc.rowptr[j+1];
	    		for (rb = rb1; rb < rb2; rb++) {
	        		k = ctxt->Aintsc.colidx[rb];
				aintscxval[i*maxAngN + k] += xval*ctxt->aintscval[rb];		
	    		}

			rb1 = ctxt->Aintss.rowptr[j];
	    		rb2 = ctxt->Aintss.rowptr[j+1];
	    		for (rb = rb1; rb < rb2; rb++) {
	        		k = ctxt->Aintss.colidx[rb];
				aintssxval[i*maxAngN + k] += xval*ctxt->aintssval[rb];		
	    		}
			
			rb1 = ctxt->Aintc.rowptr[j];
	    		rb2 = ctxt->Aintc.rowptr[j+1];
	    		for (rb = rb1; rb < rb2; rb++) {
	        		k = ctxt->Aintc.colidx[rb];
				aintcxval[i*maxAngN + k] += xval*ctxt->aintcval[rb];		
	    		}

			rb1 = ctxt->apu1.rowptr[j];
	    		rb2 = ctxt->apu1.rowptr[j+1];
	    		for (rb = rb1; rb < rb2; rb++) {
	        		k = ctxt->apu1.colidx[rb];
				apu1xval[i*maxAngN + k] += xval*ctxt->apu1val[rb];		
	    		}

			rb1 = ctxt->apu1sc.rowptr[j];
	    		rb2 = ctxt->apu1sc.rowptr[j+1];
	    		for (rb = rb1; rb < rb2; rb++) {
	        		k = ctxt->apu1sc.colidx[rb];
				apu1scxval[i*maxAngN + k] += xval*ctxt->apu1scval[rb];		
	    		}

			rb1 = ctxt->apu1ss.rowptr[j];
	    		rb2 = ctxt->apu1ss.rowptr[j+1];
	    		for (rb = rb1; rb < rb2; rb++) {
	        		k = ctxt->apu1ss.colidx[rb];
				apu1ssxval[i*maxAngN + k] += xval*ctxt->apu1ssval[rb];		
	    		}
			
			rb1 = ctxt->apu1c.rowptr[j];
	    		rb2 = ctxt->apu1c.rowptr[j+1];
	    		for (rb = rb1; rb < rb2; rb++) {
	        		k = ctxt->apu1c.colidx[rb];
				apu1cxval[i*maxAngN + k] += xval*ctxt->apu1cval[rb];		
	    		}

			rb1 = ctxt->Aintscsc.rowptr[j];
	    		rb2 =ctxt->Aintscsc.rowptr[j+1];
	    		for (rb = rb1; rb < rb2; rb++) {
	        		k = ctxt->Aintscsc.colidx[rb];
				aintscscxval[i*maxAngN + k] += xval*ctxt->aintscscval[rb];		
	    		}

			rb1 = ctxt->Aintscss.rowptr[j];
	    		rb2 = ctxt->Aintscss.rowptr[j+1];
	    		for (rb = rb1; rb < rb2; rb++) {
	        		k = ctxt->Aintscss.colidx[rb];
				aintscssxval[i*maxAngN + k] += xval*ctxt->aintscssval[rb];		
	    		}

			rb1 = ctxt->Aintssss.rowptr[j];
	    		rb2 = ctxt->Aintssss.rowptr[j+1];
	    		for (rb = rb1; rb < rb2; rb++) {
	        		k = ctxt->Aintssss.colidx[rb];
				aintssssxval[i*maxAngN + k] += xval*ctxt->aintssssval[rb];		
	    		}
			
			rb1 = ctxt->Aintscc.rowptr[j];
	    		rb2 = ctxt->Aintscc.rowptr[j+1];
	    		for (rb = rb1; rb < rb2; rb++) {
	        		k = ctxt->Aintscc.colidx[rb];
				aintsccxval[i*maxAngN + k] += xval*ctxt->aintsccval[rb];		
	    		}
			
			rb1 = ctxt->Aintssc.rowptr[j];
	    		rb2 = ctxt->Aintssc.rowptr[j+1];
	    		for (rb = rb1; rb < rb2; rb++) {
	        		k = ctxt->Aintssc.colidx[rb];
				aintsscxval[i*maxAngN + k] += xval*ctxt->aintsscval[rb];		
	    		}
			
			rb1 = ctxt->Aintcc.rowptr[j];
	    		rb2 = ctxt->Aintcc.rowptr[j+1];
	    		for (rb = rb1; rb < rb2; rb++) {
	        		k = ctxt->Aintcc.colidx[rb];
				aintccxval[i*maxAngN + k] += xval*ctxt->aintccval[rb];		
	    		}


		}
    }

    int scol;
    for(int is = 0; is < spatN; is++)
    {
	for(int ia = 0; ia < ctxt->node_angN[is]; ia++)
	{
		complex temp(0, 0);
		for(int js = ctxt->Sint.rowptr[is]; js < ctxt->Sint.rowptr[is+1]; js++)
		{
			scol = ctxt->Sint.colidx[js];
			temp += complex(0, ctxt->sintval[js])*aintxval[scol*maxAngN + ia];
			temp += complex(0, ctxt->sdxval[js])*aintscxval[scol*maxAngN + ia];
			temp += complex(0, ctxt->sdyval[js])*aintssxval[scol*maxAngN + ia];
			temp += complex(0, ctxt->sdzval[js])*aintcxval[scol*maxAngN +  ia];
			temp += aintxval[scol*maxAngN + ia]*ctxt->spata3_rteval[js];
			temp += aintscxval[scol*maxAngN + ia]*(ctxt->spata3_sdmxval[js] - ctxt->sxval[js]);
			temp += aintssxval[scol*maxAngN + ia]*(ctxt->spata3_sdmyval[js] - ctxt->syval[js]);
			temp += aintcxval[scol*maxAngN +  ia]*(ctxt->spata3_sdmzval[js] - ctxt->szval[js]);
			temp += aintscscxval[scol*maxAngN + ia]*ctxt->sdxxval[js];
			temp += aintscssxval[scol*maxAngN + ia]*(ctxt->sdxyval[js] + ctxt->sdyxval[js]);
			temp += aintssssxval[scol*maxAngN + ia]*ctxt->sdyyval[js];
			temp += aintsccxval[scol*maxAngN + ia]*(ctxt->sdxzval[js] + ctxt->sdzxval[js]);
			temp += aintsscxval[scol*maxAngN + ia]*(ctxt->sdyzval[js] + ctxt->sdzyval[js]);
			temp += aintccxval[scol*maxAngN + ia]*ctxt->sdzzval[js];
			temp -= apu1xval[scol*maxAngN + ia]*ctxt->spsval[js];
			temp -= apu1scxval[scol*maxAngN + ia]*ctxt->spsdxval[js];
			temp -= apu1ssxval[scol*maxAngN +  ia]*ctxt->spsdyval[js];
			temp -= apu1cxval[scol*maxAngN +  ia]*ctxt->spsdzval[js];
		}
		result[ctxt->offset[is] + ia]  = temp;
	} 
	}


 	RCAx(ctxt->A2, x, ctxt->A2x);
    	complex *res  = result.data_buffer();
    	complex *arg1 = ctxt->A2x.data_buffer();
    	for (int i=0; i < dof; i++)
    		*res++ += *arg1++;

    	return result;
}

CVector getDiag(void * context)
{
    MyDataContext *ctxt = (MyDataContext*) context;
    int spatN = ctxt->spatN; 
    int nDim = sum(ctxt->node_angN);
    CVector result(nDim);
    int arow, brow;
    complex a0_rte, a0_sdm, a1_rte, a1_sdm;
    complex a0, a1, a2, a3, a4;
    double  coeff = ctxt->w/ctxt->c;

    for(int i = 0; i < spatN; i++)
    {
	for(int j = 0; j < ctxt->node_angN[i]; j++)
	{
		a0_rte = complex(0, ctxt->Aint.Get(j, j) * ctxt->Sint.Get(i, i));
		a0_sdm = complex(0, ctxt->Aintsc.Get(j, j)*ctxt->Sdx.Get(i, i));
		a0_sdm += complex(0, ctxt->Aintss.Get(j, j) * ctxt->Sdy.Get(i, i));
  		a0_sdm += complex(0, ctxt->Aintc.Get(j, j)*ctxt->Sdz.Get(i, i));

		a1_rte = ctxt->Aintsc.Get(j, j)*ctxt->Sx.Get(i, i);	
		a1_rte += ctxt->Aintss.Get(j, j)*ctxt->Sy.Get(i, i);
		a1_rte += ctxt->Aintc.Get(j, j)*ctxt->Sz.Get(i, i);

		a1_sdm = ctxt->Aintscsc.Get(j, j)*ctxt->Sdxx.Get(i, i);
		a1_sdm +=  ctxt->Aintscss.Get(j, j)*ctxt->Sdxy.Get(i, i);
        	a1_sdm += ctxt->Aintscss.Get(j, j)*ctxt->Sdyx.Get(i, i);
		a1_sdm += ctxt->Aintssss.Get(j, j)*ctxt->Sdyy.Get(i, i);	
		a1_sdm += ctxt->Aintscc.Get(j, j)*ctxt->Sdxz.Get(i, i);
		a1_sdm +=  ctxt->Aintscc.Get(j, j)*ctxt->Sdzx.Get(i, i);	
		a1_sdm += ctxt->Aintssc.Get(j, j)*ctxt->Sdyz.Get(i, i);
		a1_sdm += ctxt->Aintssc.Get(j, j)*ctxt->Sdzy.Get(i, i);	
		a1_sdm += ctxt->Aintcc.Get(j, j)*ctxt->Sdzz.Get(i, i);
	
		a3 = ctxt->Aint.Get(j, j)*ctxt->spatA3_rte.Get(i, i) + ctxt->Aintsc.Get(j, j)*ctxt->spatA3_sdmx.Get(i, i) + ctxt->Aintss.Get(j, j)*ctxt->spatA3_sdmy.Get(i, i) + ctxt->Aintc.Get(j, j)*ctxt->spatA3_sdmz.Get(i, i);


		a4 = ctxt->apu1.Get(j, j)*ctxt->SPS.Get(i, i) + ctxt->apu1sc.Get(j, j)*ctxt->SPSdx.Get(i, i) + ctxt->apu1ss.Get(j, j)*ctxt->SPSdy.Get(i, i) + ctxt->apu1c.Get(j, j)*ctxt->SPSdz.Get(i, i);
	
		a2 = ctxt->A2.Get(ctxt->offset[i] + j, ctxt->offset[i] + j);

		a0 = a0_rte + a0_sdm;
	
		a1 = a1_sdm - a1_rte;
	
		result[ctxt->offset[i] + j] = a0+a1+a2+a3-a4;
	}
     }
     return result;
}


//=========================================================================


void SelectSourceProfile (int &qtype, double &qwidth, SourceMode &srctp)
{
    char cbuf[256];
    int cmd;

    bool typeok = false;
    if (pp.GetString ("SOURCETYPE", cbuf)) {
	if (!strcasecmp (cbuf, "NEUMANN")) {
	    srctp = SRCMODE_NEUMANN;
	    typeok = true;
	} else if (!strcasecmp (cbuf, "ISOTROPIC")) {
	    srctp = SRCMODE_ISOTROPIC;
	    typeok = true;
	}
    }
    while (!typeok) {
	cout << "\nSource type:\n";
	cout << "(1) Neumann boundary source\n";
	cout << "(2) Isotropic point source\n";
	cout << "[1|2] >> ";
	cin  >> cmd;
	cout<<cmd<<endl;
	switch (cmd) {
	    case 1: srctp = SRCMODE_NEUMANN;   typeok = true; break;
	    case 2: srctp = SRCMODE_ISOTROPIC; typeok = true; break;
	}
    }
    pp.PutString ("SOURCETYPE",
        srctp == SRCMODE_NEUMANN ? "NEUMANN" : "ISOTROPIC");

    qtype = -1;
    if (pp.GetString ("SOURCEPROFILE", cbuf)) {
	if (!strcasecmp (cbuf, "POINT")) {
	    qtype = 0;
	} else if (!strcasecmp (cbuf, "GAUSSIAN")) {
	    qtype = 1;
	} else if (!strcasecmp (cbuf, "COSINE")) {
	    qtype = 2;
	}
    }
    while (qtype < 0) {
	cout << "\nSource profile:\n";
	cout << "(1) Point\n";
	cout << "(2) Gaussian\n";
	cout << "(3) Cosine\n";
	cout << "[1|2|3] >> ";
	cin  >> qtype;
	cout<<qtype<<endl;
	qtype -= 1;
    }
    if (qtype > 0 && !pp.GetReal ("SOURCEWIDTH", qwidth)) {
	switch (qtype) {
	case 1:
	    cout << "\nSource 1/e radius [mm]:\n>> ";
	    break;
	case 2:
	    cout << "\nSource support radius [mm]:\n>> ";
	    break;
	}
	cin >> qwidth;
	cout<<qwidth<<endl;
    }
    switch (qtype) {
    case 0:
	pp.PutString ("SOURCEPROFILE", "POINT");
	break;
    case 1:
	pp.PutString ("SOURCEPROFILE", "GAUSSIAN");
	pp.PutReal ("SOURCEWIDTH", qwidth);
	break;
    case 2:
	pp.PutString ("SOURCEPROFILE", "COSINE");
	pp.PutReal ("SOURCEWIDTH", qwidth);
	break;
    }
}
// ============================================================================

void SelectMeasurementProfile (ParamParser &pp, int &mtype, double &mwidth)
{
    char cbuf[256];
    mtype = -1;
    if (pp.GetString ("MEASUREMENTPROFILE", cbuf)) {
	if (!strcasecmp (cbuf, "POINT")) {
	    mtype = 0;
	} else if (!strcasecmp (cbuf, "GAUSSIAN")) {
	    mtype = 1;
	} else if (!strcasecmp (cbuf, "COSINE")) {
	    mtype = 2;
	}
    }
    while (mtype < 0) {
	cout << "\nMeasurement profile:\n";
	cout << "(1) Point\n";
	cout << "(2) Gaussian\n";
	cout << "(3) Cosine\n";
	cout << "[1|2|3] >> ";
	cin  >> mtype;
	cout << mtype<<endl;
	mtype -= 1;
    }
    if (mtype > 0 && !pp.GetReal ("MEASUREMENTWIDTH", mwidth)) {
	switch (mtype) {
	case 1:
	    cout << "\nMeasurement 1/e radius [mm]:\n>> ";
	    break;
	case 2:
	    cout << "\nMeasurement support radius [mm]:\n>> ";
	    break;
	}
	cin >> mwidth;
	cout << mwidth << endl;
    }
    switch (mtype) {
    case 0:
	pp.PutString ("MEASUREMENTPROFILE", "POINT");
	break;
    case 1:
	pp.PutString ("MEASUREMENTPROFILE", "GAUSSIAN");
	pp.PutReal ("MEASUREMENTWIDTH", mwidth);
	break;
    case 2:
	pp.PutString ("MEASUREMENTPROFILE", "COSINE");
	pp.PutReal ("MEASUREMENTWIDTH", mwidth);
	break;
    }
}
void WriteData (const RVector &data, char *fname)
{
    ofstream ofs (fname);
    ofs << setprecision(14);
    ofs << data << endl;
}

void WriteDataBlock (const QMMesh &mesh, const RVector &data, char *fname)
{
    int q, m, i;
    ofstream ofs (fname);
    for (m = i = 0; m < mesh.nM; m++) {
	for (q = 0; q < mesh.nQ; q++) {
	    if (mesh.Connected (q,m)) ofs << data[i++];
	    else                      ofs << '-';
	    ofs << (q == mesh.nQ-1 ? '\n' : '\t');
	}
    }   
}

void OpenNIM (const char *nimname, const char *meshname, int size)
{
    ofstream ofs (nimname);
    ofs << "NIM" << endl;
    ofs << "Mesh = " << meshname << endl;
    ofs << "SolutionType = N/A" << endl;
    ofs << "ImageSize = " << size << endl;
    ofs << "EndHeader" << endl;
}

void WriteNIM (const char *nimname, const RVector &img, int size, int no)
{
    ofstream ofs (nimname, ios::app);
    ofs << "Image " << no << endl;
    for (int i = 0; i < size; i++)
        ofs << img[i] << ' ';
    ofs << endl;
}

bool ReadNim (char *nimname, RVector &img)
{
    char cbuf[256];
    int i, imgsize = 0;

    ifstream ifs (nimname);
    if (!ifs.getline (cbuf, 256)) return false;
    if (strcmp (cbuf, "NIM")) return false;
    do {
        ifs.getline (cbuf, 256);
	if (!strncasecmp (cbuf, "ImageSize", 9))
	    sscanf (cbuf+11, "%d", &imgsize);
    } while (strcasecmp (cbuf, "EndHeader"));
    if (!imgsize) return false;
    img.New(imgsize);
    do {
        ifs.getline (cbuf, 256);
    } while (strncasecmp (cbuf, "Image", 5));
    for (i = 0; i < imgsize; i++)
        ifs >> img[i];
    return true;
}

#ifdef WRITEPPM
void WritePGM (const RVector &img, const IVector &gdim, char *fname)
{
    int i, ii, dim = gdim[0]*gdim[1];
    double imgmin = 1e100, imgmax = -1e100;
    unsigned char *pixmap = new unsigned char[dim];

    if (gdim.Dim() > 2)
        ii = (gdim[2]/2)*gdim[0]*gdim[1]; // simply plot middle layer
    else
        ii = 0;

    for (i = 0; i < dim; i++) {
        if (img[i+ii] < imgmin) imgmin = img[i+ii];
	if (img[i+ii] > imgmax) imgmax = img[i+ii];
    }
    double scale = 256.0/(imgmax-imgmin);
    for (i = 0; i < dim; i++) {
        int v = (int)((img[i+ii]-imgmin)*scale);
	if      (v < 0  ) v = 0;
	else if (v > 255) v = 255;
	pixmap[i] = (unsigned char)v;
    }
    ofstream ofs(fname);
    ofs << "P5" << endl;
    ofs << "# CREATOR: TOAST (test_grid)" << endl;
    ofs << gdim[0] << ' ' << gdim[1] << endl;
    ofs << 255 << endl;
    for (i = 0; i < dim; i++) ofs << pixmap[i];
    ofs << endl;
    delete []pixmap;
}

void WritePPM (const RVector &img, const IVector &gdim,
    double *scalemin, double *scalemax, char *fname)
{
    typedef struct {
        unsigned char r,g,b;
    } RGB;

    int i, ii, dim = gdim[0]*gdim[1];
    unsigned char *pixmap = new unsigned char[dim];
    double imgmin, imgmax, scale;

    static RGB colmap[256];
    static bool have_colmap = false;

    if (!have_colmap) {
        int r, g, b;
        ifstream ifs (colormap);
	for (i = 0; i < 256; i++) {
	    ifs >> r >> g >> b;
	    colmap[i].r = (unsigned char)r;
	    colmap[i].g = (unsigned char)g;
	    colmap[i].b = (unsigned char)b;
	}
	have_colmap = true;
    }

    if (gdim.Dim() > 2)
        ii = (gdim[2]/2)*gdim[0]*gdim[1]; // simply plot middle layer
    else
        ii = 0;

    // rescale image
    if (scalemin) {
        imgmin = *scalemin;
    } else {
	for (i = 0, imgmin = 1e100; i < dim; i++)
	    if (img[i+ii] < imgmin) imgmin = img[i+ii];
    }
    if (scalemax) {
        imgmax = *scalemax;
    } else {
      for (i = 0, imgmax = -1e100; i < dim; i++)
	  if (img[i+ii] > imgmax) imgmax = img[i+ii];
    }
    scale = 256.0/(imgmax-imgmin);

    for (i = 0; i < dim; i++) {
        int v = (int)((img[i+ii]-imgmin)*scale);
	if      (v < 0  ) v = 0;
	else if (v > 255) v = 255;
	pixmap[i] = (unsigned char)v;
    }
    ofstream ofs(fname);
    ofs << "P6" << endl;
    ofs << "# CREATOR: TOAST (test_grid)" << endl;
    ofs << gdim[0] << ' ' << gdim[1] << endl;
    ofs << 255 << endl;
    for (i = 0; i < dim; i++)
        ofs << colmap[pixmap[i]].r
	    << colmap[pixmap[i]].g
	    << colmap[pixmap[i]].b;
    ofs << endl;

    delete []pixmap;
}

#endif
