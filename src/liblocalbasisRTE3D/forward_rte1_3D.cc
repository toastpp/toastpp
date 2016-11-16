/***************************************************************************
 * forward_rte1_3D.cc             Simon Arridge           27.11.06         *
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
#include <time.h>
#include <fstream>
#include <iomanip>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <mathlib.h>
#include "matrix.h"
#include <felib.h>
#include "source.h"
#include "pparse.h"
#include "toast.h"
#include "rte3D.h"
#define VARYDELTA
#define USE_INTONSPHERE

#define MIN(A,B) ( (A) < (B) ? (A) : (B))

class MyDataContext{
public:
double sigma0_0_0; 
ScatKernType sktyp;  
CCompRowMatrix Sint, Sdx, Sdy, Sdz;
RCompRowMatrix Sgrad, Sx, Sy, Sz;
RCompRowMatrix Sdxx, Sdxy, Sdyx, Sdyy,  Sdxz, Sdzx,  Sdyz, Sdzy, Sdzz;
RCompRowMatrix Aint, Aintsc, Aintss, Aintc, Anvec, Anvec_sc, Anvec_ss,  Anvec_c;
RCompRowMatrix Aintscsc,  Aintscss, Aintscc,  Aintssss,  Aintssc,  Aintcc;
RCompRowMatrix SPS, SPSdx, SPSdy, SPSdz;
RCompRowMatrix spatA3_rte, spatA3_sdmx, spatA3_sdmy, spatA3_sdmz;
RCompRowMatrix apu1, apu1sc, apu1ss, apu1c;
RCompRowMatrix A2;
CCompRowMatrix b1;
CDenseMatrix Aintx, Aintscx, Aintssx, Aintcx;
CDenseMatrix Aintscscx, Aintscssx, Aintssssx, Aintsccx, Aintsscx, Aintccx;
CDenseMatrix apu1x, apu1scx, apu1ssx, apu1cx, Xmat;
CVector A2x;

const std::complex<double> *sintval, *sdxval, *sdyval, *sdzval;
const double *sxval, *syval, *szval; 
const double *sdxxval, *sdxyval, *sdyxval, *sdyyval, *sdxzval, *sdzxval, *sdyzval, *sdzyval, *sdzzval; 
const double *spsval, *spsdxval, *spsdyval, *spsdzval, *spata3_rteval, *spata3_sdmxval, *spata3_sdmyval, *spata3_sdmzval; 

const double *aintval, *aintscval, *aintssval, *aintcval, *apu1val, *apu1scval, *apu1ssval, *apu1cval;
const double *aintscscval, *aintscssval, *aintssssval, *aintsccval, *aintsscval, *aintccval;
std::complex<double> *xmatval;

idxtype *browptr, *bcolidx;
int nzb;
int spatN, angN;

double w, c;

MyDataContext(QMMesh &spatMesh, Mesh &angMesh, RVector &delta, RVector &muabs, RVector &muscat, ScatKernType &sktyp, double w, double c)
{ 
	this->sktyp = sktyp;
	this->w = w;
	this->c = c;
	spatN = spatMesh.nlen();
	angN = angMesh.nlen();

	cout<<"Generating spatial integrals ..."<<endl;
	genmat_spatint_nobf_3D(spatMesh, muabs, muscat, Sint, Sgrad, Sx, Sy, Sz, SPS, spatA3_rte);
	genmat_spatint_sdm_nobf_3D(spatMesh, delta, muabs, muscat, Sdx, Sdy, Sdz, Sdxx, Sdxy, Sdyx, Sdyy,  Sdxz, Sdzx, Sdyz, Sdzy, Sdzz, SPSdx, SPSdy, SPSdz, spatA3_sdmx, spatA3_sdmy, spatA3_sdmz);

	cout<<"Generating angular integrals ..."<<endl;
	genmat_angint_3D(Aint, Aintsc, Aintss, Aintc, Anvec, angMesh);
    	genmat_angint_sdm_3D(Aintscsc,  Aintscss,   Aintscc,  Aintssss,  Aintssc,  Aintcc,  Anvec_sc, Anvec_ss,  Anvec_c, angMesh);

	cout<<"Generating phase integrals ..."<<endl;
	genmat_apu(angMesh, Anvec, Anvec_sc, Anvec_ss, Anvec_c, apu1, apu1sc, apu1ss, apu1c, sktyp);

	cout<<"Generating boundary integrals ..."<<endl;
  	genmat_boundint_3D(spatMesh,  angMesh, A2, b1);
	

	Sint = Sint*std::complex<double>(0, w/c);
   	Sdx = Sdx*std::complex<double>(0, w/c);
  	Sdy = Sdy*std::complex<double>(0, w/c);
   	Sdz = Sdz*std::complex<double>(0, w/c);

	int angN = angMesh.nlen();
	angMesh.SparseRowStructure (browptr, bcolidx, nzb);
   	
	/*Allocating memory for A2x where A2 is the matrix corresponding to boundary*/
	A2x.New(spatN*angN);
	CDenseMatrix temp(angN, spatN), temp2(spatN, angN); 
        
	
	Xmat = temp2;

	/*Preparing angular integrals for computing Kronecker products implicitly*/	
    	Aint.Transpone(); Aintsc.Transpone(); Aintss.Transpone(); Aintc.Transpone(); 
    	Aintscsc.Transpone(); Aintscss.Transpone(); Aintssss.Transpone(); Aintscc.Transpone(); 
    	Aintssc.Transpone(); Aintcc.Transpone(); apu1.Transpone();  apu1sc.Transpone();  apu1ss.Transpone();  
    	apu1c.Transpone(); 	

	/*Dereferencing the val pointers of angular matrices*/ 
	aintval = Aint.ValPtr(); aintscval = Aintsc.ValPtr(); aintssval = Aintss.ValPtr(); aintcval = Aintc.ValPtr();
	apu1val = apu1.ValPtr(); apu1scval = apu1sc.ValPtr(); apu1ssval = apu1ss.ValPtr(); apu1cval = apu1c.ValPtr();
	aintscscval = Aintscsc.ValPtr(); aintscssval = Aintscss.ValPtr(); aintssssval = Aintssss.ValPtr(); aintsccval = Aintscc.ValPtr();
	aintsscval = Aintssc.ValPtr(); aintccval = Aintcc.ValPtr();

	/*Dereferencing the val pointers of spatial matrices*/
	Aintx.New(spatN, angN); Aintscx.New(spatN, angN); Aintssx.New(spatN, angN); Aintcx.New(spatN, angN);
	apu1x.New(spatN, angN); apu1scx.New(spatN, angN); apu1ssx.New(spatN, angN); apu1cx.New(spatN, angN);
	Aintscscx.New(spatN, angN); Aintscssx.New(spatN, angN); Aintssssx.New(spatN, angN); Aintsccx.New(spatN, angN);
	Aintsscx.New(spatN, angN); Aintccx.New(spatN, angN);

	Xmat.New(spatN, angN);
	xmatval = Xmat.data_buffer();

	sintval = Sint.ValPtr(); sdxval = Sdx.ValPtr(); sdyval = Sdy.ValPtr();
	sdzval = Sdz.ValPtr(); sxval = Sx.ValPtr(); syval = Sy.ValPtr(); 
	szval = Sz.ValPtr(); sdxxval = Sdxx.ValPtr(); sdxyval = Sdxy.ValPtr(); 
	sdyxval = Sdyx.ValPtr(); sdyyval = Sdyy.ValPtr(); sdxzval = Sdxz.ValPtr(); 
	sdzxval = Sdzx.ValPtr(); sdyzval = Sdyz.ValPtr(); sdzyval = Sdzy.ValPtr(); 
	sdzzval = Sdzz.ValPtr(); spsval = SPS.ValPtr(); spsdxval = SPSdx.ValPtr(); 
	spsdyval = SPSdy.ValPtr(); spsdzval = SPSdz.ValPtr(); spata3_rteval = spatA3_rte.ValPtr(); 
	spata3_sdmxval = spatA3_sdmx.ValPtr(); spata3_sdmyval = spatA3_sdmy.ValPtr(); 
	spata3_sdmzval = spatA3_sdmz.ValPtr();      
	
}

~MyDataContext()
{
  delete []browptr;
  delete []bcolidx;
};

private:

/*Generating the spatial matrices required by RTE*/
void genmat_spatint_nobf_3D(const Mesh& mesh, const RVector &muabs, const RVector &muscat, CCompRowMatrix& Sint, RCompRowMatrix& Sgrad, RCompRowMatrix& Sx, RCompRowMatrix& Sy, RCompRowMatrix& Sz, RCompRowMatrix &SPS, RCompRowMatrix &spatA3_rte)
{
   int sysdim = mesh.nlen();       // dimensions are size of nodes.

   idxtype *rowptr, *colidx;
   int nzero;
   int *elrowptr, *elcolidx;
   mesh.SparseRowStructure (rowptr, colidx, nzero);
   Sint.New (sysdim, sysdim);
   Sint.Initialise (rowptr, colidx);
   Sgrad.New (sysdim, sysdim);
   Sgrad.Initialise (rowptr, colidx);
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
   
   double elk_ij, elb_ij, elsx_ij, elsy_ij, elsz_ij;
   int el, nodel, i, j, k,is, js;

   for (el = 0; el < mesh.elen(); el++) {
	nodel = mesh.elist[el]->nNode();

	for (i = 0; i < nodel; i++) {
	    if ((is = mesh.elist[el]->Node[i]) >= sysdim) continue;
	    for (j = 0; j < nodel; j++) {
		if ((js = mesh.elist[el]->Node[j]) >= sysdim) continue;
		elb_ij = mesh.elist[el]->IntFF (i, j);
		Sint(is,js) += elb_ij;
		SPS(is, js) += elb_ij*muscat[el];
		spatA3_rte(is, js) += elb_ij*(muabs[el] + muscat[el]);
		elk_ij = mesh.elist[el]->IntDD (i, j);
		Sgrad(is,js) += elk_ij;
		elsx_ij = mesh.elist[el]->IntFd (j,i,0);
		Sx(is,js) += elsx_ij;
		elsy_ij = mesh.elist[el]->IntFd (j,i,1);
		Sy(is,js) += elsy_ij;
		elsz_ij = mesh.elist[el]->IntFd (j,i,2);
		Sz(is,js) += elsz_ij;
	    }
	}
   }

   delete []rowptr;
   delete []colidx;


}

/*Generating the spatial matrices required by RTE with SDM*/
void genmat_spatint_sdm_nobf_3D(const Mesh& mesh,  const RVector& delta, const RVector &muabs, const RVector &muscat, CCompRowMatrix& Sdx, CCompRowMatrix& Sdy,  CCompRowMatrix& Sdz, RCompRowMatrix& Sdxx, RCompRowMatrix& Sdxy, RCompRowMatrix& Sdyx, RCompRowMatrix& Sdyy, RCompRowMatrix& Sdxz, RCompRowMatrix& Sdzx, RCompRowMatrix& Sdyz, RCompRowMatrix& Sdzy, RCompRowMatrix& Sdzz, RCompRowMatrix &SPSdx, RCompRowMatrix &SPSdy, RCompRowMatrix &SPSdz, RCompRowMatrix &spatA3_sdmx, RCompRowMatrix &spatA3_sdmy, RCompRowMatrix &spatA3_sdmz)
{
   int sysdim = mesh.nlen();       // dimensions are size of nodes.

   idxtype *rowptr, *colidx;
   int nzero;
   int *elrowptr, *elcolidx;
   mesh.SparseRowStructure (rowptr, colidx, nzero);

   Sdx.New (sysdim, sysdim);
   Sdx.Initialise (rowptr, colidx);
   Sdy.New (sysdim, sysdim);
   Sdy.Initialise (rowptr, colidx);
   Sdz.New (sysdim, sysdim);
   Sdz.Initialise (rowptr, colidx);
   Sdxx.New (sysdim, sysdim);
   Sdxx.Initialise (rowptr, colidx);
   Sdxy.New (sysdim, sysdim);
   Sdxy.Initialise (rowptr, colidx);
   Sdyx.New (sysdim, sysdim);
   Sdyx.Initialise (rowptr, colidx);
   Sdyy.New (sysdim, sysdim);
   Sdyy.Initialise (rowptr, colidx);
   Sdxz.New (sysdim, sysdim);
   Sdxz.Initialise (rowptr, colidx);
   Sdzx.New (sysdim, sysdim);
   Sdzx.Initialise (rowptr, colidx);
   Sdyz.New (sysdim, sysdim);
   Sdyz.Initialise (rowptr, colidx);
   Sdzy.New (sysdim, sysdim);
   Sdzy.Initialise (rowptr, colidx);
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
	
	RSymMatrix eldd = mesh.elist[el]->Intdd(); // "all at once!""

	for (i = 0; i < nodel; i++) {
	    if ((is = mesh.elist[el]->Node[i]) >= sysdim) continue;
	    for (j = 0; j < nodel; j++) {
		if ((js = mesh.elist[el]->Node[j]) >= sysdim) continue;

		elsx_ij = mesh.elist[el]->IntFd (j,i,0);
       		Sdx(is,js) += dss*elsx_ij;
		elsy_ij = mesh.elist[el]->IntFd (j,i,1);
	      	Sdy(is,js) += dss*elsy_ij;
		elsz_ij = mesh.elist[el]->IntFd (j,i,2);
		Sdz(is,js) += dss*elsz_ij;
	       		
		SPSdx(is, js) += dss*elsx_ij*muscat[el];
		SPSdy(is, js) += dss*elsy_ij*muscat[el];
		SPSdz(is, js) += dss*elsz_ij*muscat[el];

		spatA3_sdmx(is, js) += dss*elsx_ij*(muabs[el] + muscat[el]);
		spatA3_sdmy(is, js) += dss*elsy_ij*(muabs[el] + muscat[el]);
		spatA3_sdmz(is, js) += dss*elsz_ij*(muabs[el] + muscat[el]);

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

/** Thresholds and shrinks a real dense matrix to give a real sparse matrix
* NOTE!! The threshold is set to 1e-15
**/
RCompRowMatrix shrink(const RDenseMatrix &dnsmat)
{
    int i, j;
    idxtype *rowptr, *colidx;
    double *val;
    int m =  dnsmat.nRows(), n= dnsmat.nCols();

    rowptr = new idxtype[m+1];
    rowptr[0] = 0;
    for (i = 0; i < m; i++) {
	int nz=0;
	for (j = 0; j < n; j++)
		if(!(fabs(dnsmat.Get(i, j)) < 1e-15))
			nz++;
	rowptr[i+1] = rowptr[i] + nz;
    }
    
    colidx = new idxtype[rowptr[m]];
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

/**
Computes all the angular integrals required by RTE
**/
void genmat_angint_3D(RCompRowMatrix& Aint, RCompRowMatrix& Aintsc, RCompRowMatrix& Aintss, RCompRowMatrix& Aintc, RCompRowMatrix& Anvec, const Mesh& S2mesh)
{
  const int& SN =  S2mesh.nlen();       // dimensions are size of nodes.
  const int& SE =  S2mesh.elen();       // number of spherical elements.

  idxtype *angrowptr, *angcolidx;
  int nzero;
  S2mesh.SparseRowStructure (angrowptr, angcolidx, nzero);
  Anvec.New(SN,SN);  // and this one...
  Anvec.Initialise(angrowptr, angcolidx);

  RDenseMatrix dnsAint(SN, SN), dnsAintsc(SN, SN), dnsAintss(SN, SN), dnsAintc(SN, SN);


  int el, i,j,is,js;


  RVector sx(SN);    // sample Sin t Cos p on nodes
  RVector sy(SN);    // sample Sin t Sin p on nodes
  RVector sz(SN);    // sample Cos t on nodes

  for (i = 0; i < SN; i++) {  // create the samples on sphere
    const Node& np = S2mesh.nlist[i] ; // do we really need to copy it ?
    double ilen = 1.0/length(np); // normalise by length, just to be sure
    sx[i] = np[0]*ilen;
    sy[i] = np[1]*ilen;
    sz[i] = np[2]*ilen;
  }

  for(el = 0; el < SE ; el++){
    for(i = 0; i < S2mesh.elist[el]->nNode(); i++) {
       if ((is = S2mesh.elist[el]->Node[i]) >= SN) continue;
       for(j = 0; j < S2mesh.elist[el]->nNode(); j++) {
	 if ((js = S2mesh.elist[el]->Node[j]) >= SN) continue;
#ifdef USE_INTONSPHERE
	 dnsAint(is,js) += S2mesh.elist[el]->IntUnitSphereFF (S2mesh.nlist,i, j);
	 dnsAintsc(is,js) += S2mesh.elist[el]->IntUnitSpherePFF (S2mesh.nlist,i, j, sx);
	 dnsAintss(is,js) += S2mesh.elist[el]->IntUnitSpherePFF (S2mesh.nlist,i, j, sy);
	 dnsAintc(is,js)  += S2mesh.elist[el]->IntUnitSpherePFF (S2mesh.nlist,i, j, sz);
#else
	 dnsAint(is,js) += S2mesh.elist[el]->IntFF (i, j);
	 dnsAintsc(is,js) += S2mesh.elist[el]->IntPFF (i, j, sx);
	 dnsAintss(is,js) += S2mesh.elist[el]->IntPFF (i, j, sy);
	 dnsAintc(is,js)  += S2mesh.elist[el]->IntPFF (i, j, sz);
#endif
       }
    }
  }

 Aint = shrink(dnsAint); Aintsc = shrink(dnsAintsc); Aintss = shrink(dnsAintss); Aintc = shrink(dnsAintc);

 for(int jj = 0; jj < SN; jj++) {  // loop on "directions
   RVector Anvt(SN);
   int Jind;

    for(int el = 0; el < SE; el ++) { // find elements containing node jj
      if(!S2mesh.elist[el]->IsNode(jj)) continue ;
      // node is in this element. Find local number
      Jind = -1; // used to break search
      for( i = 0; i < S2mesh.elist[el]->nNode() && Jind < 0; i++) {
	if( S2mesh.elist[el]->Node[i] == jj)
	  Jind = i;
      }
      // now integrate shape function of node jj with this element
      for(i = 0; i < S2mesh.elist[el]->nNode(); i++) {
        if ((is = S2mesh.elist[el]->Node[i]) >= SN) continue;  
#ifdef USE_INTONSPHERE
	Anvt[is] += S2mesh.elist[el]->IntUnitSphereFF (S2mesh.nlist, i, Jind);
#else
	Anvt[is] += S2mesh.elist[el]->IntFF (i, Jind);
#endif
      } 
    } // end loop on elements
    // and finally, insert this vector of integrals into global matrix
    for(j = 0; j < SN; j++) {    // insert into "global matrices"
      if(Anvt[j] == 0.0) continue;      // this could be done better...
      Anvec(j,jj) = Anvt[j];
    }
  } // end loop on "directions"
 delete []angrowptr;
 delete []angcolidx;

}

/**
Computes all the angular integrals required by RTE with SDM
**/
void genmat_angint_sdm_3D(RCompRowMatrix& Aintscsc,  RCompRowMatrix& Aintscss,  RCompRowMatrix& Aintscc,  RCompRowMatrix& Aintssss,  RCompRowMatrix& Aintssc,  RCompRowMatrix& Aintcc, RCompRowMatrix& Anvec_sc, RCompRowMatrix& Anvec_ss, RCompRowMatrix& Anvec_c, const Mesh& S2mesh)
{
  const int& SN =  S2mesh.nlen();       // dimensions are size of nodes.
  const int& SE =  S2mesh.elen();       // number of spherical elements.

  idxtype *angrowptr, *angcolidx;
  int nzero;
  S2mesh.SparseRowStructure (angrowptr, angcolidx, nzero);

  RDenseMatrix dnsAintscsc(SN, SN), dnsAintscss(SN, SN), dnsAintcc(SN, SN);
  RDenseMatrix dnsAintscc(SN, SN), dnsAintssss(SN, SN), dnsAintssc(SN, SN);

         // and these three
  Anvec_sc.New(SN,SN);  
  Anvec_sc.Initialise(angrowptr, angcolidx);
  Anvec_ss.New(SN,SN);  
  Anvec_ss.Initialise(angrowptr, angcolidx);
  Anvec_c.New(SN,SN);  
  Anvec_c.Initialise(angrowptr, angcolidx);

  int el, i,j,is,js;


  RVector sx(SN);     // sample Sin t Cos p on nodes
  RVector sy(SN);     // sample Sin t Sin p on nodes
  RVector sz(SN);     // sample Cos t on nodes

  RVector sxx(SN);    // sample (Sin t Cos p)^2 on nodes
  RVector sxy(SN);    // sample (Sin t Cos p)(Sin t Sin p) on nodes
  RVector sxz(SN);    // sample (Sin t Cos p) Cos t on nodes
  RVector syy(SN);    // sample (Sin t Sin p)^2 on nodes
  RVector syz(SN);    // sample (Sin t Sin p) Cos t  on nodes
  RVector szz(SN);    // sample ( Cos t)^2 on nodes

  for (i = 0; i < SN; i++) {  // create the samples on sphere
    const Node& np = S2mesh.nlist[i] ; // do we really need to copy it ?
    double ilen = 1.0/length(np); // normalise by length, just to be sure
    sx[i] = np[0]*ilen;
    sy[i] = np[1]*ilen;
    sz[i] = np[2]*ilen;

    sxx[i] = np[0]*ilen*np[0]*ilen;
    sxy[i] = np[0]*ilen*np[1]*ilen;
    sxz[i] = np[0]*ilen*np[2]*ilen;
    syy[i] = np[1]*ilen*np[1]*ilen;
    syz[i] = np[1]*ilen*np[2]*ilen;
    szz[i] = np[2]*ilen*np[2]*ilen;
  }

  for(el = 0; el < SE ; el++){
    for(i = 0; i < S2mesh.elist[el]->nNode(); i++) {
       if ((is = S2mesh.elist[el]->Node[i]) >= SN) continue;
       for(j = 0; j < S2mesh.elist[el]->nNode(); j++) {
	 if ((js = S2mesh.elist[el]->Node[j]) >= SN) continue;
#ifdef USE_INTONSPHERE
	 dnsAintscsc(is,js) += S2mesh.elist[el]->IntUnitSpherePFF (S2mesh.nlist, i, j, sxx);
	 dnsAintscss(is,js) += S2mesh.elist[el]->IntUnitSpherePFF (S2mesh.nlist, i, j, sxy);
	 dnsAintscc(is,js)  += S2mesh.elist[el]->IntUnitSpherePFF (S2mesh.nlist, i, j, sxz);
	 dnsAintssss(is,js) += S2mesh.elist[el]->IntUnitSpherePFF (S2mesh.nlist, i, j, syy);
	 dnsAintssc(is,js)  += S2mesh.elist[el]->IntUnitSpherePFF (S2mesh.nlist, i, j, syz);
	 dnsAintcc(is,js)   += S2mesh.elist[el]->IntUnitSpherePFF (S2mesh.nlist, i, j, szz);
#else
	 dnsAintscsc(is,js) += S2mesh.elist[el]->IntPFF (i, j, sxx);
	 dnsAintscss(is,js) += S2mesh.elist[el]->IntPFF (i, j, sxy);
	 dnsAintscc(is,js)  += S2mesh.elist[el]->IntPFF (i, j, sxz);
	 dnsAintssss(is,js) += S2mesh.elist[el]->IntPFF (i, j, syy);
	 dnsAintssc(is,js)  += S2mesh.elist[el]->IntPFF (i, j, syz);
	 dnsAintcc(is,js)   += S2mesh.elist[el]->IntPFF (i, j, szz);
#endif
       }
    }
  }

  Aintscsc = shrink(dnsAintscsc);Aintscss = shrink(dnsAintscss);Aintscc = shrink(dnsAintscc); Aintssc = shrink(dnsAintssc);
  Aintssss = shrink(dnsAintssss);Aintcc = shrink(dnsAintcc);

  for(int jj = 0; jj < SN; jj++) {  // loop on "directions
    RVector Anvtsc(SN);
    RVector Anvtss(SN);
    RVector Anvtc(SN);
    int Jind;

    for(int el = 0; el < SE; el ++) { // find elements containing node jj
      if(!S2mesh.elist[el]->IsNode(jj)) continue ;
      // node is in this element. Find local number
      Jind = -1; // used to break search
      for( i = 0; i < S2mesh.elist[el]->nNode() && Jind < 0; i++) {
	if( S2mesh.elist[el]->Node[i] == jj)
	  Jind = i;
      }
      // now integrate shape function of node jj with this element
      for(i = 0; i < S2mesh.elist[el]->nNode(); i++) {
        if ((is = S2mesh.elist[el]->Node[i]) >= SN) continue;  
#ifdef USE_INTONSPHERE
	Anvtsc[is] += S2mesh.elist[el]->IntUnitSpherePFF (S2mesh.nlist,i, Jind,sx);
	Anvtss[is] += S2mesh.elist[el]->IntUnitSpherePFF (S2mesh.nlist,i, Jind,sy);
	Anvtc[is]  += S2mesh.elist[el]->IntUnitSpherePFF (S2mesh.nlist,i, Jind,sz);
#else
	Anvtsc[is] += S2mesh.elist[el]->IntPFF (i, Jind,sx);
	Anvtss[is] += S2mesh.elist[el]->IntPFF (i, Jind,sy);
	Anvtc[is]  += S2mesh.elist[el]->IntPFF (i, Jind,sz);
#endif
      } 
    } // end loop on elements
    // and finally, insert this vector of integrals into global matrix
    for(j = 0; j < SN; j++) {    // insert into "global matrices"
      if(Anvtsc[j] == 0.0) continue;      // this could be done better...
      Anvec_sc(j,jj) = Anvtsc[j];
      if(Anvtss[j] == 0.0) continue;      // this could be done better...
      Anvec_ss(j,jj) = Anvtss[j];
      if(Anvtc[j] == 0.0) continue;       // this could be done better...
      Anvec_c(j,jj) = Anvtc[j];
    }
  } // end loop on "directions"
 delete []angrowptr;
   delete []angcolidx;

}

/** Phase integrals 
*  NOTE!! This function has been specialized for g=0
**/
void genmat_apu(const Mesh &S2mesh, RCompRowMatrix &Anvec, RCompRowMatrix &Anvec_sc, RCompRowMatrix &Anvec_ss, RCompRowMatrix &Anvec_c, RCompRowMatrix& apu1, RCompRowMatrix& apu1sc, RCompRowMatrix& apu1ss, RCompRowMatrix& apu1c, ScatKernType sktyp)
{
	const int& SN = S2mesh.nlen();
	apu1.New(SN,SN);
        apu1sc.New(SN,SN);
        apu1ss.New(SN,SN);
        apu1c.New(SN,SN);
        switch(sktyp)
	{
	 	case  MUSHOMOG_G0 :
       	 	{
			 for(int bb = 0; bb < SN; bb++){ // not efficient !
	 			RCompRowMatrix acol = Anvec.Subcols(bb,bb+1);
	 			RCompRowMatrix acolsc = Anvec_sc.Subcols(bb,bb+1);
	 			RCompRowMatrix acolss = Anvec_ss.Subcols(bb,bb+1);
	 			RCompRowMatrix acolc  = Anvec_c.Subcols(bb,bb+1);
	 	 		for(int cc = 0; cc < SN; cc++){
	   				RCompRowMatrix bcol = Anvec.Subcols(cc,cc+1);
	   				bcol.Transpone();// row vector of column

	   				apu1 +=  kron(acol,bcol);
	   				apu1sc += kron(acolsc,bcol);
	   				apu1ss += kron(acolss,bcol);
	   				apu1c  += kron(acolc,bcol);
	 			}
       			}
			apu1 *= 1/(4*M_PI);//-1*sigma[0];
			apu1sc *= 1/(4*M_PI);//-1*sigma[0];
			apu1ss *= 1/(4*M_PI);//-1*sigma[0];
			apu1c *= 1/(4*M_PI);//-1*sigma[0];
     
			break;
		 }
		 default :
                       cout<<"Unknown g type (g should be zero for this case) "<<endl;
		}
} 

/** Preallocating memory for boundary integral terms.
The matrix is considered to be sparse both spatially and angularly
which is eventually shrunk after computation.
**/
void initialiseA2b1(const Mesh &mesh, const Mesh &S2mesh, RCompRowMatrix& A2, CCompRowMatrix& b1)
{
    int ia, ib, ka, kb, ja, jb, i, idx;
    int va_i, vb_i, v_i;
    int snzero, anzero;
    idxtype *srowptr, *scolidx, *arowptr, *acolidx; 
    mesh.SparseRowStructure(srowptr, scolidx, snzero);
    S2mesh.SparseRowStructure(arowptr, acolidx, anzero);

   int sysdim = mesh.nlen();
   int el, nodel, j, k,is, js;
 
   int *status = new int[srowptr[sysdim]];
   for(i=0; i<snzero; i++)
	status[i] = 0; // 1 implies nonzero 0 denotes zero;
  /*'status' variable updates the nodes on the boundary*/ 
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
		for (int rp = srowptr[is]; rp < srowptr[is+1]; rp++){
        		if (scolidx[rp] == js) status[rp] = 1;
	    }
	  } 
	 }
       }
    }
    
    int tnzero=0;
    for(i=0; i < snzero; i++)
	if(status[i]) tnzero++;
   /*The new spatial sparsity pattern limited just to the boundary is computed*/ 
   int *spatrowptr, *spatcolidx;
   spatrowptr = new int[sysdim + 1];
   spatcolidx = new int[tnzero];
   spatrowptr[0] = 0;
   j=0;
   for(i = 0; i < sysdim; i++)
   {
	int rp1 = srowptr[i];
	int rp2 = srowptr[i+1];
        k=0;
	for(int rp = rp1; rp < rp2; rp++)
	{
		if(status[rp]){ 
			k++;
			spatcolidx[j] = scolidx[rp];
			j++;			
		} 
	}
	spatrowptr[i+1] = spatrowptr[i] + k;
   } 
   delete []status;




    int na = mesh.nlen(), ma = mesh.nlen();
    int nb = S2mesh.nlen(), mb = S2mesh.nlen();
    int n = na*nb, m = ma*mb, v = snzero*anzero;

    idxtype *rowptr = new idxtype[n+1];
    idxtype *colidx = new idxtype[v];

    rowptr[0] = 0;
    i = idx = 0;
    /*The overall sparsity pattern is computed based on spatial and angular sparsity patterns*/
    for (ia = 0; ia < na; ia++) {
	va_i = spatrowptr[ia+1] - spatrowptr[ia]; // nonzeros in row ia of A
	for (ib = 0; ib < nb; ib++) {
	    vb_i = arowptr[ib+1] - arowptr[ib]; // nonzeros in row ib of B
	    v_i = va_i * vb_i;
	    rowptr[i+1] = rowptr[i] + v_i;
	    for (ka = spatrowptr[ia]; ka < spatrowptr[ia+1]; ka++) {
		ja = spatcolidx[ka];
		for (kb = arowptr[ib]; kb < arowptr[ib+1]; kb++) {
		    jb = acolidx[kb];
		    colidx[idx] = ja*mb + jb;
		    idx++;
		}
	    }
	    i++;
	}
    }
   A2.Initialise(rowptr, colidx);
   b1.Initialise(rowptr, colidx);
   A2.Zero(); 
   b1.Zero();

   delete []rowptr;
   delete []colidx;
   delete []srowptr;
   delete []scolidx;
   delete []spatrowptr;
   delete []spatcolidx;
   delete []arowptr;
   delete []acolidx;

}

/** Adds aB to its appropriate place in the system matrix
* spatrow -> spatial row where 'a' is drawn from
* spatcol -> spatial column where 'a' is drawn from
* node_angN -> number of angular degrees of freedom for all the spatial nodes
* offset -> starting location in the system matrix for each spatial node
* a_ij -> 'a'
* B -> B
* C -> output (System matrix) 
**/
template<class MT>
void kronsdplus(const int spatrow, const int spatcol, const int angN, const double a_ij, const RCompRowMatrix &B, TCompRowMatrix<MT>& C)
{
    MT *Cval = C.ValPtr();
    int jb, ib;
    int row_offset = 0, col_offset = 0;
   
    for(int i=0; i < spatrow; i++) row_offset += angN;
    for(int i=0; i < spatcol; i++) col_offset += angN; 
 
    int *browptr, *bcolidx;
    int col;

    for(int ib=0; ib < angN; ib++)
    {
     for(int jb = B.rowptr[ib]; jb<B.rowptr[ib+1]; jb++)
     {
        col = B.colidx[jb];

	C(row_offset+ib , col_offset+col) = C.Get(row_offset+ib , col_offset+col) + B.Get(ib, col)*a_ij; 	
     }
    }  
}

/**Compute the boundary integral terms
**/
void genmat_boundint_3D(const Mesh& mesh,  const Mesh& S2mesh, RCompRowMatrix& A2, CCompRowMatrix& b1)
{
   const int& SN =  S2mesh.nlen();       // dimensions are size of nodes.
   const int& SE =  S2mesh.elen();       // number of spherical elements.
   const int sysdim = mesh.nlen();       // dimensions are size of nodes.
   const int fullsysdim = sysdim*SN;     // full size of angles X space nodes

   A2.New (fullsysdim, fullsysdim);
   b1.New(fullsysdim, fullsysdim);
   initialiseA2b1(mesh, S2mesh, A2, b1);

   double ela_ij;
   int el, nodel, i, j, k,is, js;

   // first create structure for angular integrals 
   idxtype *angrowptr, *angcolidx;
   int nzero;
   S2mesh.SparseRowStructure (angrowptr, angcolidx, nzero);

   RVector sx(SN);    // sample Sin t Cos p on nodes
   RVector sy(SN);    // sample Sin t Sin p on nodes
   RVector sz(SN);    // sample Cos t on nodes

   for (i = 0; i < SN; i++) {  // create the samples on sphere
    const Node& np = S2mesh.nlist[i] ; // do we really need to copy it ?
    double ilen = 1.0/length(np); // normalise by length, just to be sure
    sx[i] = np[0]*ilen;
    sy[i] = np[1]*ilen;
    sz[i] = np[2]*ilen;
   }
  
   RCompRowMatrix  Angbintplus(SN,SN), Angbintminus(SN, SN);
   Angbintplus.Initialise(angrowptr, angcolidx);
   Angbintminus.Initialise(angrowptr, angcolidx);
   RVector f1(SN), f2(SN); 
   // now create matrix, by looping over elements that have a boundary
   for (el = 0; el < mesh.elen(); el++) {
	if(!(el*100%mesh.elen()))
		cout<<el*100/mesh.elen() <<"% work done ..."<<endl;

        if(!mesh.elist[el]->HasBoundarySide ()) continue;

	nodel = mesh.elist[el]->nNode();
 
	// now determine the element integrals

	for(int sd = 0; sd <  mesh.elist[el]->nSide(); sd++) {
	  // if sd is not a boundary side. skip 
	  if(!mesh.elist[el]->IsBoundarySide (sd)) continue;
 	    
	  // get boundary normal...
	  RVector nhat = mesh.ElDirectionCosine(el,sd);
	  Angbintplus.Zero(); Angbintminus.Zero();
	  int ia,ja,isa,jsa,ela;
	  for( j = 0; j < SN; j++) {
	  	double tmp = nhat[0]*sx[j] + nhat[1]*sy[j] + nhat[2]*sz[j]; 
	  	f1[j] =  (tmp > 0.0 ? tmp : 0.0);
	  	f2[j] = (tmp < 0.0 ? -tmp : 0.0);
	   }
	   for(ela = 0; ela < SE ; ela++){
	   	for(ia = 0; ia < S2mesh.elist[ela]->nNode(); ia++) {
	   		if ((isa = S2mesh.elist[ela]->Node[ia]) >= SN) continue;
	   		for(ja = 0; ja < S2mesh.elist[ela]->nNode(); ja++) {
				if ((jsa = S2mesh.elist[ela]->Node[ja]) >= SN) continue;
#ifdef USE_INTONSPHERE
				Angbintplus(isa,jsa) += S2mesh.elist[ela]->IntUnitSpherePFF (S2mesh.nlist, ia, ja, f1);
				Angbintminus(isa,jsa) += S2mesh.elist[ela]->IntUnitSpherePFF (S2mesh.nlist, ia, ja, f2);

#else
				Angbintplus(isa,jsa) += S2mesh.elist[ela]->IntPFF (ia, ja, f1);
				Angbintminus(isa,jsa) += S2mesh.elist[ela]->IntPFF (ia, ja, f2);

#endif

  			}
	
	    	}
	  }
 
	  for (i = 0; i < nodel; i++) {
	    if ((is = mesh.elist[el]->Node[i]) >= sysdim) continue;
	    for (j = 0; j < nodel; j++) {
		if ((js = mesh.elist[el]->Node[j]) >= sysdim) continue;
		ela_ij = mesh.elist[el]->SurfIntFF (i, j,sd);
		kronsdplus(is, js, SN, ela_ij, Angbintplus, A2);
	 	kronsdplus(is, js, SN, ela_ij, Angbintminus, b1);	
	     }
	   }
	} // end loop on element sides

   } // end loop on elements
   A2.Shrink();
   b1.Shrink();
 delete []angrowptr;
   delete []angcolidx;

}
};

/**Computes source vector for a single point source in the interior of the domain
**/
void genmat_b2(const Mesh &mesh, const Mesh& S2mesh, CCompRowMatrix& Svec, const int Nsource, const RVector& dirVec, const bool is_cosine, const int angMesh_node)
 {
   int el, nodel, i, j, k,is, js;
   int sysdim = mesh.nlen();       // dimensions are size of nodes.
   int SN = S2mesh.nlen();
   int fullsysdim = sysdim*SN;   // full size of angles X space nodes
   
   idxtype *rowptr, *colidx;
   rowptr = new idxtype[fullsysdim+1];
   colidx = new idxtype[fullsysdim];
   rowptr[0] = 0;
   for(int i=1;  i <= fullsysdim; i++)
	rowptr[i] = i;
   for(int i=0; i<fullsysdim; i++)
	colidx[i] = 0;
   Svec.New (fullsysdim,1);
   Svec.Initialise(rowptr, colidx);

   int offset = 0;
   RDenseMatrix dirMat(1, 3);

   RVector sx(SN);    // sample Sin t Cos p on nodes
   RVector sy(SN);    // sample Sin t Sin p on nodes
   RVector sz(SN);    // sample Cos t on nodes

   for (i = 0; i < SN; i++) {  // create the samples on sphere
    const Node& np = S2mesh.nlist[i] ; // do we really need to copy it ?
    double ilen = 1.0/length(np); // normalise by length, just to be sure
    sx[i] = np[0]*ilen;
    sy[i] = np[1]*ilen;
    sz[i] = np[2]*ilen;
   }
  

   // now create vector, by looping over elements that have a boundary
   for (el = 0; el < mesh.elen(); el++) {
        if(!mesh.elist[el]->IsNode(Nsource)) continue; // source not in this el
	for(int sd = 0; sd <  mesh.elist[el]->nSide(); sd++) {
	   dirMat(0, 0) = dirVec[0]; dirMat(0, 1) = dirVec[1]; dirMat(0, 2) = dirVec[2];

	  for(int nd = 0; nd < mesh.elist[el]->nSideNode(sd); nd++) {
	    is = mesh.elist[el]->SideNode(sd,nd);
	    js = mesh.elist[el]->Node[is];
	    if(js == Nsource)
	    {
		offset = 0;
		for(int i=0; i < js; i++) offset += SN; 

	    	if(is_cosine)
		{ 
			for(int i=0; i < SN; i++) Svec(offset + i, 0) = 1.0;
		}
	    	else{
			
			Svec(offset + angMesh_node, 0) = 1.0;
		}
	    }

	  }
	} // end loop on element sides

   } // end loop on elements

   delete []rowptr;
   delete []colidx;
}

// =========================================================================
// global parameters

SourceMode srctp = SRCMODE_NEUMANN;   // source type
ParamParser pp;
QMMesh qmmesh;
NodeList &nlist=qmmesh.nlist;
ElementList &elist=qmmesh.elist;
inline CVector matrixFreeCaller(const CVector& x, void * context);
void SelectSourceProfile (int &qtype, double &qwidth, SourceMode &srctp);
void SelectMeasurementProfile (ParamParser &pp, int &mtype, double &mwidth);
void genmat_toastsource_3D(CCompRowMatrix* & Source, const Mesh& mesh, const Mesh& S2mesh, const CCompRowMatrix qvec, const int ns, const CCompRowMatrix& b1, const bool is_cosine, const int angMesh_node);
void genmat_toastsourcevalvector_3D(CCompRowMatrix& Svec, const Mesh& mesh, const Mesh& S2mesh,  const CCompRowMatrix qvec, const int iq, const bool is_cosine, const int angMesh_node);
void genmat_source_3D(CCompRowMatrix* & Source, const Mesh& mesh,  const Mesh& S2mesh, const int* Nsource, const int ns, const CCompRowMatrix& b1, const bool is_cosine, const int angMesh_node);
void genmat_sourcevalvector_3D(CCompRowMatrix& Svec, const Mesh& mesh, const Mesh& S2mesh, const int Nsource, const bool is_cosine, const int angMesh_node);
void WriteData (const RVector &data, char *fname);
void WriteDataBlock (const QMMesh &mesh, const RVector &data, char *fname);
void OpenNIM (const char *nimname, const char *meshname, int size);
void WriteNIM (const char *nimname, const RVector &img, int size, int no);
bool ReadNim (char *nimname, RVector &img);
void WritePGM (const RVector &img, const IVector &gdim, char *fname);
void WritePPM (const RVector &img, const IVector &gdim,
double *scalemin, double *scalemax, char *fname);
CVector getDiag(void * context);

// main routine **************************************************************

int main (int argc, char *argv[])
{
    char cbuf[200];
    int el;
    
    cout << "Reading mesh " << argv[1] << endl;
    ifstream ifs;
    ifs.open (argv[1]);
    xASSERT (ifs.is_open(), "Mesh file not found.");
    ifs >> qmmesh;
    xASSERT (ifs.good(), "Problem reading mesh.");
    ifs.close ();
    cout << "Spatial mesh has " << qmmesh.elen() << " elements, " << qmmesh.nlen()
	 << " nodes\n";
    int dimension = nlist[0].Dim();
    for (int i = 1; i < nlist.Len(); i++)
	xASSERT(nlist[i].Dim() == dimension, "Inconsistent node dimensions.");
    xASSERT(dimension >= 2 && dimension <= 3, "Mesh dimension must be 2 or 3.");
    qmmesh.Setup();


    Mesh S2Mesh;

    ifstream ifm(argv[2]);
    ifm >> S2Mesh;
    ifm.close();
    cout << "Angular mesh has " << S2Mesh.elen() << " elements, " << S2Mesh.nlen()
         << " nodes\n";
    S2Mesh.Setup();
    const int& SN =  S2Mesh.nlen();       // dimensions are size of nodes.

    char file_extn[200];
    cin>>file_extn;
    cout<<"File name prefix: "<<file_extn<<endl;// prefix for the output files

    int ns, nM;
    int *Nsource;
    int    qprof, mprof;   // source/measurement profile (0=Gaussian, 1=Cosine)
    double qwidth, mwidth; // source/measurement support radius [mm]
    CCompRowMatrix qvec, mvec;
    if(argc < 4) { //point source
    	cin>>ns;
    	Nsource = new int[ns];
    	for(int i=0; i<ns; i++)
    	{
    		cin>>Nsource[i];
    		cout << i<<"th source node coords: "<<nlist[Nsource[i]]<<endl;
	}
    }
    else {
    cout << "QM file " << argv[3] << endl;
    ifs.open (argv[3]);
    xASSERT (ifs.is_open(), "QM file not found.");
    qmmesh.LoadQM (ifs);
    xASSERT (ifs.good(), "Problem reading QM.");
    ifs.close ();
    ns = qmmesh.nQ;
    nM = qmmesh.nM;
    cout << ns << " sources\n";
    SelectSourceProfile (qprof, qwidth, srctp);
    SelectMeasurementProfile (pp, mprof, mwidth);
    // build the source vectors
    qvec.New (ns, qmmesh.nlen());
    for (int i = 0; i < ns; i++) {
	CVector q(qmmesh.nlen());
	switch (qprof) {
	case 0:
	    SetReal (q, QVec_Point (qmmesh, qmmesh.Q[i], srctp));
	    break;
	case 1:
	    SetReal (q, QVec_Gaussian (qmmesh, qmmesh.Q[i], qwidth, srctp));
	    break;
	case 2:
	    SetReal (q, QVec_Cosine (qmmesh, qmmesh.Q[i], qwidth, srctp));
	    break;
	}
	qvec.SetRow (i, q);
    }
    // build the measurement vectors
    mvec.New (nM, qmmesh.nlen());
    LOGOUT1_INIT_PROGRESSBAR ("Meas. vectors", 50, nM);
    for (int i = 0; i < nM; i++) {
	CVector m(qmmesh.nlen());
	switch (mprof) {
	case 0:
	    SetReal (m, QVec_Point (qmmesh, qmmesh.M[i], SRCMODE_NEUMANN));
	    break;
	case 1:
	    SetReal (m, QVec_Gaussian (qmmesh, qmmesh.M[i], mwidth,
				       SRCMODE_NEUMANN));
	    break;
	case 2:
	    SetReal (m, QVec_Cosine (qmmesh, qmmesh.M[i], mwidth,
				     SRCMODE_NEUMANN));
	    break;
	}

	// **** WARNING: This doesn't work anymore!
	// **** Need alternative way to retrieve C2A
	//for (int j = 0; j < qmmesh.nlen(); j++) 
	//  m[j] *= qmmesh.plist[j].C2A();
	mvec.SetRow (i, m);
    }
    }
    //******** optical parameters
    double freq = 0, hg;
    int is_iso;
    bool is_cosine;
    RVector dirVec(3);
    
    cin >> hg;  
    cout << "value of g (WARNING !! This g is considered constant throughout the domain): "<<hg<<endl;
    ScatKernType sktyp = (hg == 0.0 ?   MUSHOMOG_G0 :   MUSHOMOG_GCONST);
    cin >> freq; 
    cout << "value for frequency (MHz) : "<<freq<<endl;
    cin>> is_iso;
    is_cosine = is_iso>0 ? true : false;
    cout<<"The source is cosine or directed (1. Cosine 0. Directed): "<<is_cosine<<endl;
    cin>>dirVec[0]; cin>>dirVec[1]; cin>>dirVec[2];
    dirVec = dirVec*1.0/length(dirVec);// normalize the direction vector just in case
    cout<<"The direction vector for the source (if it is directed): "<<dirVec[0]<<", "<<dirVec[1]<<", "<<dirVec[2]<<endl;
    double w = freq * 2.0*M_PI*1e-6;
    double c = 0.3;


    RVector muscat(qmmesh.elen());
    RVector muabs(qmmesh.elen());
    RVector ref(qmmesh.elen());
    RVector g(qmmesh.elen());

    cout<< "Reading mua, mus and refractive index values for all the elements from file: "<<endl;
    for(el = 0; el <  qmmesh.elen(); el++){ // assign something to parameters
      cin >> muabs[el];
      cin >> muscat[el];
      cin >> ref[el]; 
      g[el] = hg; 
    }

    RVector sigmatot(qmmesh.elen());
    RVector sigma(qmmesh.elen()); // to be assigned
    RVector intst;
    
    calc_paramdistr_nobf_new_3D(sigma, sigmatot, intst, qmmesh, muscat, muabs,  g, S2Mesh);

    // "smoothing parameter delta of streamline diffusion modification.
    // Just set to constant for now.
    RVector delta(qmmesh.elen());
    double min = 1e20, max = -1e20;
    for(el = 0; el <  qmmesh.elen(); el++){ // assign something to delta 
#ifdef VARYDELTA 
        double sval =  elist[el]->Size()/((muscat[el]));
	delta[el] = sval;
#else
	 min = max = delta[0];
#endif 
    }

    MyDataContext ctxt(qmmesh, S2Mesh, delta, muabs, muscat, sktyp, w, c);

    RVector norm(S2Mesh.nlen());
    for(int i=0; i < S2Mesh.nlen(); i++)
    {
	Node &np = S2Mesh.nlist[i];
	double ilen = 1.0/length(np); // normalise by length, just to be sure
    	np[0] = np[0]*ilen;
    	np[1] = np[1]*ilen;
    	np[2] = np[2]*ilen;

	norm[i] = sqrt((np[0] - dirVec[0])*(np[0] - dirVec[0]) + (np[1] - dirVec[1])*(np[1] - dirVec[1])+(np[2] - dirVec[2])*(np[2] - dirVec[2]));
    }

    double minnorm=norm[0];
    int angMesh_node=0;
    for(int i=1; i < SN; i++){
	if(norm[i] < minnorm) {minnorm = norm[i];angMesh_node = i;}
    }
  
    cout<<"The source is being directed along: "<<S2Mesh.nlist[angMesh_node]; 
    CVector proj(ns*nM); // projection data
    CCompRowMatrix* Source;
    if( !(Source = new  CCompRowMatrix [ns]))
          cerr << "Memory Allocation error Source = new  RCompRowMatrix\n";
    if(argc<4)
      genmat_source_3D(Source, qmmesh, S2Mesh , Nsource, ns, ctxt.b1, is_cosine, angMesh_node);
    else
      genmat_toastsource_3D(Source, qmmesh, S2Mesh, qvec, ns, ctxt.b1, is_cosine, angMesh_node);

   CCompRowMatrix b2;
   genmat_b2(qmmesh, S2Mesh, b2, Nsource[0], dirVec, is_cosine, angMesh_node);
 
    cout << "calculating the radiance\n";
    int sysdim = qmmesh.nlen();
    int fullsysdim = sysdim * SN;
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
    for (int j = 0; j < ns ; j++) {   
      cout << "Radiance with the source number:  " << j << endl;

      Phi[j].New(fullsysdim);
      Phisum[j].New(sysdim);
      for(int i = 0; i < fullsysdim; i++){
	RHS[i] = Source[j].Get(i,0) + b2.Get(i, 0);
	}
	
      GMRES(&matrixFreeCaller, &ctxt, RHS, Phi[j], tol, AACP, 100);
      osRHS << RHS << endl;
      osPhi << "Phi " << j << "\n" << Phi[j] << endl;
      for (int k = 0; k < sysdim; k++) {
	for (int t = 0; t < SN; t++){
	  Phisum[j][k] += Phi[j][k*SN+t]*intst[t];
	}
      }
      osPhisum << "Phisum " << j << "\n" << Phisum[j] << endl;
      for (int im = 0; im < nM; im++) {
	for(int in = 0; in < sysdim; in++)
	  proj[j*nM + im] += mvec.Get(im,in) * Phisum[j][in];
      }
    }
    clock_t end = clock();
    osPhi.close();
    osPhisum.close();
    osRHS.close();
    
    char flnmod[300], farg[300], ftime[300];
    strcpy(flnmod, file_extn); strcpy(farg, file_extn);strcpy(ftime, file_extn);
    strcat(flnmod, "_lnmod.nim"); strcat(farg, "_arg.nim");strcat(ftime, "_time.txt");
    OpenNIM (flnmod, argv[1], sysdim);
    OpenNIM (farg, argv[1], sysdim);
    for (int i = 0; i < ns; i++) {
	   WriteNIM (flnmod, LogMod(Phisum[i]), sysdim, i);
	   WriteNIM (farg, Arg(Phisum[i]), sysdim, i);
    }
    cout << "  Log Mod field written to "<< flnmod << endl;
    cout << "  Arg field written to "<<farg << endl;
 
    FILE *fid;
    fid = fopen(ftime, "w");
    fprintf(fid, "Time taken by solver: %f\n", (double)(end-start)/CLOCKS_PER_SEC);
    fclose(fid);

    cout<<"The solver took "<<(double)(end-start)/CLOCKS_PER_SEC<<" seconds"<<endl;

	delete []Phi;
	delete []Phisum;
	delete []Nsource;

}

/** Computes Ax 
	A -> Real sparse matrix
	x -> Complex vector
     NOTE!! It's rightful place in crmatrix class using templates
**/
inline void RCAx(const RCompRowMatrix &A, const CVector& x, CVector &res)
{
    dASSERT(x.Dim() == A.nCols(),
    "Parameter 1 invalid size (expected %d, actual %d)",
    A.nCols(), x.Dim());
    if (res.Dim() != A.nRows()) res.New(A.nRows());

    int r, i, i2;
    std::complex<double> br;
    const double *aval;
    aval = A.ValPtr();

    for (r = i = 0; r < A.nRows();) {
	i2 = A.rowptr[r+1];
	for (br = std::complex<double>(0, 0); i < i2; i++)
	    br += x[A.colidx[i]]*aval[i];
	res[r++] = br;
    }

}

/** Computes Sx required by the GMRES solver
*	S -> System matrix
*	x -> current guess of the solution
**/
inline CVector matrixFreeCaller(const CVector& x, void * context)
{
    MyDataContext *ctxt = (MyDataContext*) context;
    	
	int spatN = ctxt->spatN;
	int angN = ctxt->angN;
	int dof = spatN*angN;
	CVector result(dof);
	
	/*Implict Kronecker product implementation
	*	(S \circplus A)x = Sx_{r}A^{T}
	* where 'x_{r}' is a matrix resulting form reshaping of 'x'. 
	*/
	
	/*Reshaping 'x' to 'x_{r}'*/
	memcpy (ctxt->xmatval, x.data_buffer(), dof*
		sizeof(std::complex<double>));

	int i, j, k, m, ra, ra1, ra2, rb, rb1, rb2;
    	int nr = ctxt->Xmat.nRows();
    	int nc = ctxt->Xmat.nCols();
    	
	/*Intialize Ax's to zero where A is the angular matrix*/	
	ctxt->Aintx.Zero(); ctxt->Aintscx.Zero(); ctxt->Aintssx.Zero(); ctxt->Aintcx.Zero();
	ctxt->apu1x.Zero(); ctxt->apu1scx.Zero(); ctxt->apu1ssx.Zero(); ctxt->apu1cx.Zero();
	ctxt->Aintscscx.Zero(); ctxt->Aintscssx.Zero(); ctxt->Aintssssx.Zero(); ctxt->Aintsccx.Zero();
	ctxt->Aintsscx.Zero(); ctxt->Aintccx.Zero();

	/*Dereference the val pointers of Ax's*/
	std::complex<double> *aintxval = ctxt->Aintx.data_buffer();
	std::complex<double> *aintscxval = ctxt->Aintscx.data_buffer();
	std::complex<double> *aintssxval = ctxt->Aintssx.data_buffer();
	std::complex<double> *aintcxval = ctxt->Aintcx.data_buffer();
	std::complex<double> *apu1xval = ctxt->apu1x.data_buffer();
	std::complex<double> *apu1scxval = ctxt->apu1scx.data_buffer();  
	std::complex<double> *apu1ssxval = ctxt->apu1ssx.data_buffer();
	std::complex<double> *apu1cxval = ctxt->apu1cx.data_buffer(); 
	std::complex<double> *aintscscxval = ctxt->Aintscscx.data_buffer();
	std::complex<double> *aintscssxval = ctxt->Aintscssx.data_buffer();  
	std::complex<double> *aintssssxval = ctxt->Aintssssx.data_buffer();
	std::complex<double> *aintsccxval = ctxt->Aintsccx.data_buffer();  
	std::complex<double> *aintsscxval = ctxt->Aintsscx.data_buffer();
	std::complex<double> *aintccxval = ctxt->Aintccx.data_buffer();  

	/*Computing x_{r}A^{T}*/
	std::complex<double> xval;
    	for (i = 0; i < nr; i++) {
		for (j = 0; j < nc; j++) {
			xval = ctxt->Xmat.Get(i, j);

	    		rb1 = ctxt->Aint.rowptr[j];
	    		rb2 = ctxt->Aint.rowptr[j+1];
	    		for (rb = rb1; rb < rb2; rb++) {
	        		k = ctxt->Aint.colidx[rb];
				aintxval[i*angN + k] += xval*ctxt->aintval[rb];		
	    		}

			rb1 = ctxt->Aintsc.rowptr[j];
	    		rb2 = ctxt->Aintsc.rowptr[j+1];
	    		for (rb = rb1; rb < rb2; rb++) {
	        		k = ctxt->Aintsc.colidx[rb];
				aintscxval[i*angN + k] += xval*ctxt->aintscval[rb];		
	    		}

			rb1 = ctxt->Aintss.rowptr[j];
	    		rb2 = ctxt->Aintss.rowptr[j+1];
	    		for (rb = rb1; rb < rb2; rb++) {
	        		k = ctxt->Aintss.colidx[rb];
				aintssxval[i*angN + k] += xval*ctxt->aintssval[rb];		
	    		}
			
			rb1 = ctxt->Aintc.rowptr[j];
	    		rb2 = ctxt->Aintc.rowptr[j+1];
	    		for (rb = rb1; rb < rb2; rb++) {
	        		k = ctxt->Aintc.colidx[rb];
				aintcxval[i*angN + k] += xval*ctxt->aintcval[rb];		
	    		}

			rb1 = ctxt->apu1.rowptr[j];
	    		rb2 = ctxt->apu1.rowptr[j+1];
	    		for (rb = rb1; rb < rb2; rb++) {
	        		k = ctxt->apu1.colidx[rb];
				apu1xval[i*angN + k] += xval*ctxt->apu1val[rb];		
	    		}

			rb1 = ctxt->apu1sc.rowptr[j];
	    		rb2 = ctxt->apu1sc.rowptr[j+1];
	    		for (rb = rb1; rb < rb2; rb++) {
	        		k = ctxt->apu1sc.colidx[rb];
				apu1scxval[i*angN + k] += xval*ctxt->apu1scval[rb];		
	    		}

			rb1 = ctxt->apu1ss.rowptr[j];
	    		rb2 = ctxt->apu1ss.rowptr[j+1];
	    		for (rb = rb1; rb < rb2; rb++) {
	        		k = ctxt->apu1ss.colidx[rb];
				apu1ssxval[i*angN + k] += xval*ctxt->apu1ssval[rb];		
	    		}
			
			rb1 = ctxt->apu1c.rowptr[j];
	    		rb2 = ctxt->apu1c.rowptr[j+1];
	    		for (rb = rb1; rb < rb2; rb++) {
	        		k = ctxt->apu1c.colidx[rb];
				apu1cxval[i*angN + k] += xval*ctxt->apu1cval[rb];		
	    		}

			rb1 = ctxt->Aintscsc.rowptr[j];
	    		rb2 =ctxt->Aintscsc.rowptr[j+1];
	    		for (rb = rb1; rb < rb2; rb++) {
	        		k = ctxt->Aintscsc.colidx[rb];
				aintscscxval[i*angN + k] += xval*ctxt->aintscscval[rb];		
	    		}

			rb1 = ctxt->Aintscss.rowptr[j];
	    		rb2 = ctxt->Aintscss.rowptr[j+1];
	    		for (rb = rb1; rb < rb2; rb++) {
	        		k = ctxt->Aintscss.colidx[rb];
				aintscssxval[i*angN + k] += xval*ctxt->aintscssval[rb];		
	    		}

			rb1 = ctxt->Aintssss.rowptr[j];
	    		rb2 = ctxt->Aintssss.rowptr[j+1];
	    		for (rb = rb1; rb < rb2; rb++) {
	        		k = ctxt->Aintssss.colidx[rb];
				aintssssxval[i*angN + k] += xval*ctxt->aintssssval[rb];		
	    		}
			
			rb1 = ctxt->Aintscc.rowptr[j];
	    		rb2 = ctxt->Aintscc.rowptr[j+1];
	    		for (rb = rb1; rb < rb2; rb++) {
	        		k = ctxt->Aintscc.colidx[rb];
				aintsccxval[i*angN + k] += xval*ctxt->aintsccval[rb];		
	    		}
			
			rb1 = ctxt->Aintssc.rowptr[j];
	    		rb2 = ctxt->Aintssc.rowptr[j+1];
	    		for (rb = rb1; rb < rb2; rb++) {
	        		k = ctxt->Aintssc.colidx[rb];
				aintsscxval[i*angN + k] += xval*ctxt->aintsscval[rb];		
	    		}
			
			rb1 = ctxt->Aintcc.rowptr[j];
	    		rb2 = ctxt->Aintcc.rowptr[j+1];
	    		for (rb = rb1; rb < rb2; rb++) {
	        		k = ctxt->Aintcc.colidx[rb];
				aintccxval[i*angN + k] += xval*ctxt->aintccval[rb];		
	    		}


		}
    }
   
    /*Computing S(x_{r}A^{T})*/
    int scol;
    for(int is = 0; is < spatN; is++)
    {
	for(int ia = 0; ia < angN; ia++)
	{
	    std::complex<double> temp(0, 0);
	    for(int js = ctxt->Sint.rowptr[is]; js < ctxt->Sint.rowptr[is+1]; js++)
	        {
			scol = ctxt->Sint.colidx[js];
			temp += aintxval[scol*angN +  ia]*(ctxt->sintval[js] + ctxt->spata3_rteval[js]);
			temp += aintscxval[scol*angN +  ia]*(ctxt->sdxval[js] + ctxt->spata3_sdmxval[js] - ctxt->sxval[js]);
			temp +=aintssxval[scol*angN +  ia]*(ctxt->sdyval[js] + ctxt->spata3_sdmyval[js] - ctxt->syval[js]);
			temp += aintcxval[scol*angN +  ia]*(ctxt->sdzval[js] + ctxt->spata3_sdmzval[js] - ctxt->szval[js]);
			temp += aintscscxval[scol*angN +  ia]*ctxt->sdxxval[js];
			temp += aintscssxval[scol*angN +  ia]*(ctxt->sdxyval[js] + ctxt->sdyxval[js]);
			temp += aintssssxval[scol*angN +  ia]*ctxt->sdyyval[js];
			temp += aintsccxval[scol*angN +  ia]*(ctxt->sdxzval[js] + ctxt->sdzxval[js]);
			temp += aintsscxval[scol*angN +  ia]*(ctxt->sdyzval[js] + ctxt->sdzyval[js]);
			temp += aintccxval[scol*angN +  ia]*ctxt->sdzzval[js];
			temp -= apu1xval[scol*angN +  ia]*ctxt->spsval[js];
			temp -= apu1scxval[scol*angN +  ia]*ctxt->spsdxval[js];
			temp -= apu1ssxval[scol*angN +  ia]*ctxt->spsdyval[js];
			temp -= apu1cxval[scol*angN +  ia]*ctxt->spsdzval[js];
		}
		result[is*angN + ia]  = temp;
	} 
    }

    /*Computing A_{2}x explicitly where A_{2} is the matrix resulting from boundary*/
    RCAx(ctxt->A2, x, ctxt->A2x);
    std::complex<double> *res  = result.data_buffer();
    std::complex<double> *arg1 = ctxt->A2x.data_buffer();
    for (int i=0; i < dof; i++)
    	*res++ += *arg1++;

    return result;
   
   }

/** Computes the diagonal of the system matrix which is required for preconditioning
**/
CVector getDiag(void * context)
{
    MyDataContext *ctxt = (MyDataContext*) context;
    int spatN = ctxt->spatN; 
    int angN = ctxt->angN;
    int nDim = spatN*angN;
    CVector result(nDim);
    int arow, brow;
    std::complex<double> a0_rte, a0_sdm, a1_rte, a1_sdm;
    std::complex<double> a0, a1, a2, a3, a4;
    double  coeff = ctxt->w/ctxt->c;
	
    for(int j=0; j < nDim; j++)
    {
	arow = j/angN; brow = j%angN;
	a0_rte = ctxt->Sint.Get(arow, arow)*ctxt->Aint.Get(brow, brow);
	a0_sdm = ctxt->Sdx.Get(arow, arow)*ctxt->Aintsc.Get(brow, brow);
	a0_sdm += ctxt->Sdy.Get(arow, arow)*ctxt->Aintss.Get(brow, brow);
  	a0_sdm += ctxt->Sdz.Get(arow, arow)*ctxt->Aintc.Get(brow, brow);

	a1_rte = ctxt->Aintsc.Get(brow, brow)*ctxt->Sx.Get(arow, arow);	
	a1_rte += ctxt->Aintss.Get(brow, brow)*ctxt->Sy.Get(arow, arow);
	a1_rte += ctxt->Aintc.Get(brow, brow)*ctxt->Sz.Get(arow, arow);

	a1_sdm = ctxt->Aintscsc.Get(brow, brow)*ctxt->Sdxx.Get(arow, arow);
	a1_sdm +=  ctxt->Aintscss.Get(brow, brow)*ctxt->Sdxy.Get(arow, arow);
        a1_sdm += ctxt->Aintscss.Get(brow, brow)*ctxt->Sdyx.Get(arow, arow);
	a1_sdm += ctxt->Aintssss.Get(brow, brow)*ctxt->Sdyy.Get(arow, arow);	
	a1_sdm += ctxt->Aintscc.Get(brow, brow)*ctxt->Sdxz.Get(arow, arow);
	a1_sdm +=  ctxt->Aintscc.Get(brow, brow)*ctxt->Sdzx.Get(arow, arow);	
	a1_sdm += ctxt->Aintssc.Get(brow, brow)*ctxt->Sdyz.Get(arow, arow);
	a1_sdm += ctxt->Aintssc.Get(brow, brow)*ctxt->Sdzy.Get(arow, arow);	
	a1_sdm += ctxt->Aintcc.Get(brow, brow)*ctxt->Sdzz.Get(arow, arow);
	
	a3 = ctxt->Aint.Get(brow, brow)*ctxt->spatA3_rte.Get(arow, arow) + ctxt->Aintsc.Get(brow, brow)*ctxt->spatA3_sdmx.Get(arow, arow) + ctxt->Aintss.Get(brow, brow)*ctxt->spatA3_sdmy.Get(arow, arow) + ctxt->Aintc.Get(brow, brow)*ctxt->spatA3_sdmz.Get(arow, arow);


	a4 = ctxt->apu1.Get(brow, brow)*ctxt->SPS.Get(arow, arow) + ctxt->apu1sc.Get(brow, brow)*ctxt->SPSdx.Get(arow, arow) + ctxt->apu1ss.Get(brow, brow)*ctxt->SPSdy.Get(arow, arow) + ctxt->apu1c.Get(brow, brow)*ctxt->SPSdz.Get(arow, arow);
	
	a2 = ctxt->A2.Get(j, j);

        a0 = a0_rte + a0_sdm;
	
	a1 = a1_sdm - a1_rte;
	
	result[j] = a0;
	result[j] += a1;
	result[j] += a2;
	result[j] += a3;
	result[j] += a4;
	}
 
        return result;
}


//=========================================================================

/** Computes source vectors for point sources on the boundary 
**/
void genmat_source_3D(CCompRowMatrix* & Source, const Mesh& mesh,  const Mesh& S2mesh, const int* Nsource, const int ns, const CCompRowMatrix& b1, const bool is_cosine, const int angMesh_node)
{
   const int& SN =  S2mesh.nlen();   // dimensions are size of nodes.
   int sysdim = mesh.nlen();         // dimensions are size of nodes.
   int fullsysdim = sysdim*SN;       // full size of angles X space nodes

   if( !(Source = new  CCompRowMatrix [ns]))
          cerr << "Memory Allocation error Source = new  RCompRowMatrix\n";

   CCompRowMatrix Svec;
   for (int i = 0; i < ns; i++) {
     genmat_sourcevalvector_3D(Svec, mesh, S2mesh, Nsource[i], is_cosine, angMesh_node);
     Source[i].New(fullsysdim,1);
     b1.AB(Svec,Source[i]);        // Toast weirdness
   }
}


/**Computes source vector for a single point source on the boundary
**/
void genmat_sourcevalvector_3D(CCompRowMatrix& Svec, const Mesh& mesh, const Mesh& S2mesh, const int Nsource, const bool is_cosine, const int angMesh_node)
{
   int el, nodel, i, j, k,is, js;
   const int& SN =  S2mesh.nlen();       // dimensions are size of nodes.
   const int& SE =  S2mesh.elen();       // number of spherical elements.
   int sysdim = mesh.nlen();       // dimensions are size of nodes.
   int fullsysdim = sysdim*SN;   // full size of angles X space nodes
   idxtype* srp;
   if( !(srp = new idxtype [sysdim+1]))
      cerr << "Memory Allocation error srp = new int\n";
   idxtype* sci;
   if( !(sci = new idxtype [sysdim]))   // structure of boundary element matrix
      cerr << "Memory Allocation error sci = new int\n";
   for(i = 0; i < sysdim; i++) // make it dense for simplicity
      sci[i] = 0;

   Svec.New (fullsysdim,1);

   idxtype *arp;
   if( !(arp = new idxtype [SN+1]))
      cerr << "Memory Allocation error arp = new int\n";
   for(i = 0; i <= SN; i++) // make it dense for simplicity
      arp[i] = i;

   idxtype *aci;
   if( !(aci = new idxtype [SN]))   // structure of boundary element matrix
      cerr << "Memory Allocation error sci = new int\n";
   for(i = 0; i < SN; i++) // make it dense for simplicity
      aci[i] = 0;



   idxtype *angrowptr, *angcolidx;
   int nzero;
   S2mesh.SparseRowStructure (angrowptr, angcolidx, nzero);

   RVector sx(SN);    // sample Sin t Cos p on nodes
   RVector sy(SN);    // sample Sin t Sin p on nodes
   RVector sz(SN);    // sample Cos t on nodes

   for (i = 0; i < SN; i++) {  // create the samples on sphere
    const Node& np = S2mesh.nlist[i] ; // do we really need to copy it ?
    double ilen = 1.0/length(np); // normalise by length, just to be sure
    sx[i] = np[0]*ilen;
    sy[i] = np[1]*ilen;
    sz[i] = np[2]*ilen;
   }

   // now create vector, by looping over elements that have a boundary
   for (el = 0; el < mesh.elen(); el++) {
        if(!mesh.elist[el]->IsNode(Nsource)) continue; // source not in this el
        if(!mesh.elist[el]->HasBoundarySide()) continue;

	for(int sd = 0; sd <  mesh.elist[el]->nSide(); sd++) {
	  // if sd is not a boundary side. skip 
	  if(!mesh.elist[el]->IsBoundarySide (sd)) continue;

	  for(i = 0; i < sysdim+1; i++) // make it dense for simplicity
	    srp[i] = i;
	  CCompRowMatrix sint(sysdim,1);
	  sint.Initialise(srp,sci);

	  for(int nd = 0; nd < mesh.elist[el]->nSideNode(sd); nd++) {
	    is = mesh.elist[el]->SideNode(sd,nd);
	    js = mesh.elist[el]->Node[is];
	    sint(js,0) = std::complex<double>(1.0, 0);
	  }
	  // get boundary normal...
	  RVector nhat = mesh.ElDirectionCosine(el,sd);
	  RVector shat = -nhat; // inward directed vector

	  // now do angular integrals
	  CCompRowMatrix  Angsvec(SN,1);
	  Angsvec.Initialise(arp, aci);
	  RVector f(SN);    // sample (shat . nhat) on nodes
	  int ia,ja,isa,jsa,ela;
	  for( j = 0; j < SN; j++) {
	    double tmp = nhat[0]*sx[j] + nhat[1]*sy[j] + nhat[2]*sz[j]; 
	    f[j] =  (tmp < 0.0 ? -1*tmp : 0.0);
	  }
	  for(ela = 0; ela < SE ; ela++){
	    for(ia = 0; ia < S2mesh.elist[ela]->nNode(); ia++) {
	      if ((isa = S2mesh.elist[ela]->Node[ia]) >= SN) continue;
	      for(ja = 0; ja < S2mesh.elist[ela]->nNode(); ja++) {
		if ((jsa = S2mesh.elist[ela]->Node[ja]) >= SN) continue;
#ifdef USE_INTONSPHERE
		if(is_cosine)
			Angsvec(isa,0) += S2mesh.elist[ela]->IntUnitSpherePFF ( S2mesh.nlist, ia, ja, f);
#else
		if(is_cosine)
			Angsvec(isa,0) += S2mesh.elist[ela]->IntPFF (ia, ja, f);
#endif
	      }
	    }
	  }
	  if(!is_cosine)
		Angsvec(angMesh_node, 0) = 1.0;
  	  Svec += kron(sint,Angsvec);		
	} // end loop on element sides
   } // end loop on elements
   delete [] sci;
   delete [] srp;
   delete [] aci;
   delete [] arp;
}


/* Computes the source vectors for all the boundary sources when a QM file has been specified
*/
void genmat_toastsource_3D(CCompRowMatrix* & Source, const Mesh& mesh, const Mesh& S2mesh,  const CCompRowMatrix qvec, const int ns, const CCompRowMatrix& b1, const bool is_cosine, const int angMesh_node)
{
   const int& SN =  S2mesh.nlen();       // dimensions are size of nodes.
   const int& SE =  S2mesh.elen();       // number of spherical elements.
   int sysdim = mesh.nlen();       // dimensions are size of nodes.
   int fullsysdim = sysdim*SN;   // full size of angles X space nodes
   if( !(Source = new  CCompRowMatrix [ns]))
          cerr << "Memory Allocation error Source = new  RCompRowMatrix\n";

   CCompRowMatrix Svec;
   for (int i = 0; i < ns; i++) {
     genmat_toastsourcevalvector_3D(Svec, mesh,  S2mesh, qvec,i, is_cosine, angMesh_node);
     Source[i].New(fullsysdim,1);
     b1.AB(Svec,Source[i]);        // Toast weirdness
   }
}

/* Computes the source vector per a boundary source when a QM file has been specified
*/
void genmat_toastsourcevalvector_3D(CCompRowMatrix& Svec, const Mesh& mesh, 
const Mesh& S2mesh,  const CCompRowMatrix qvec, const int iq, const bool is_cosine, const int angMesh_node)
{
   int el, nodel, i, j, k,is, js;
   const int& SN =  S2mesh.nlen();       // dimensions are size of nodes.
   const int& SE =  S2mesh.elen();       // number of spherical elements.
   int sysdim = mesh.nlen();       // dimensions are size of nodes.
   int fullsysdim = sysdim*SN;   // full size of angles X space nodes
   idxtype* srp;
   if( !(srp = new idxtype [sysdim+1]))
      cerr << "Memory Allocation error srp = new int\n";
   idxtype* sci;
   if( !(sci = new idxtype [sysdim]))   // structure of boundary element matrix
      cerr << "Memory Allocation error sci = new int\n";
   for(i = 0; i < sysdim; i++) // make it dense for simplicity
      sci[i] = 0;

   Svec.New (fullsysdim,1);

   idxtype *arp;
   if( !(arp = new idxtype [SN+1]))
      cerr << "Memory Allocation error arp = new int\n";
   for(i = 0; i <= SN; i++) // make it dense for simplicity
      arp[i] = i;

   idxtype *aci;
   if( !(aci = new idxtype [SN]))   // structure of boundary element matrix
      cerr << "Memory Allocation error sci = new int\n";
   for(i = 0; i < SN; i++) // make it dense for simplicity
      aci[i] = 0;

   idxtype *angrowptr, *angcolidx;
   int nzero;
   S2mesh.SparseRowStructure (angrowptr, angcolidx, nzero);

   RVector sx(SN);    // sample Sin t Cos p on nodes
   RVector sy(SN);    // sample Sin t Sin p on nodes
   RVector sz(SN);    // sample Cos t on nodes


   for (i = 0; i < SN; i++) {  // create the samples on sphere
    const Node& np = S2mesh.nlist[i] ; // do we really need to copy it ?
    double ilen = 1.0/length(np); // normalise by length, just to be sure
    sx[i] = np[0]*ilen;
    sy[i] = np[1]*ilen;
    sz[i] = np[2]*ilen;
   }
   
   for(i = 0; i < sysdim+1; i++) // make it dense for simplicity
	    srp[i] = i;
   CCompRowMatrix sint(sysdim,1);
   sint.Initialise(srp,sci);

   CCompRowMatrix  Angsvec(SN,1);
   Angsvec.Initialise(arp, aci);
   RVector f(SN);    // sample (shat . nhat) on nodes
 
   // now create vector, by looping over elements that have a boundary
   for (int jq = qvec.rowptr[iq]; jq < qvec.rowptr[iq+1]; jq++) {
     int Nsource = qvec.colidx[jq];
     double sweight = norm(qvec.Get(iq,Nsource)); // we have to assume source is real

     for (el = 0; el < mesh.elen(); el++) {
        if(!mesh.elist[el]->IsNode(Nsource)) continue; // source not in this el
        if(!mesh.elist[el]->HasBoundarySide()) continue;

	for(int sd = 0; sd <  mesh.elist[el]->nSide(); sd++) {
	  // if sd is not a boundary side. skip 
	  if(!mesh.elist[el]->IsBoundarySide (sd)) continue;

          sint.Zero();
	  
	  for(int nd = 0; nd < mesh.elist[el]->nSideNode(sd); nd++) {
	    is = mesh.elist[el]->SideNode(sd,nd);
	    js = mesh.elist[el]->Node[is];
	    sint(js,0) = std::complex<double>(1.0, 0);
	  }
	  // get boundary normal...
	  RVector nhat = mesh.ElDirectionCosine(el,sd);
	  RVector shat = -nhat; // inward directed vector

	  // now do angular integrals
	  Angsvec.Zero();
	  int ia,ja,isa,jsa,ela;
	  //	  cout << "(s.n)- : ";
	  for( j = 0; j < SN; j++) {
	    double tmp = nhat[0]*sx[j] + nhat[1]*sy[j] + nhat[2]*sz[j]; 
	    f[j] =  (tmp < 0.0 ? tmp : 0.0);
	  }
	  for(ela = 0; ela < SE ; ela++){
	    for(ia = 0; ia < S2mesh.elist[ela]->nNode(); ia++) {
	      if ((isa = S2mesh.elist[ela]->Node[ia]) >= SN) continue;
	      for(ja = 0; ja < S2mesh.elist[ela]->nNode(); ja++) {
		if ((jsa = S2mesh.elist[ela]->Node[ja]) >= SN) continue;
#ifdef USE_INTONSPHERE
		if(is_cosine)
			Angsvec(isa,0) += S2mesh.elist[ela]->IntUnitSpherePFF ( S2mesh.nlist, ia, ja, f);
#else
		if(is_cosine)
			Angsvec(isa,0) += S2mesh.elist[ela]->IntPFF (ia, ja, f);
#endif
	      }
	    }
	  }
	  if(!is_cosine)
		Angsvec(angMesh_node, 0) = 1.0;
  	  Svec += kron(sint,Angsvec)*std::complex<double>(sweight, 0);		
	} // end loop on element sides
     } // end loop on elements
   }
   delete [] sci;
   delete [] srp;
   delete [] aci;
   delete [] arp;
}

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
