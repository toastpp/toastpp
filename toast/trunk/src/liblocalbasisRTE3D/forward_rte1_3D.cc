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
#include <fstream.h>
#include <iomanip.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <mathlib.h>
#include "matrix.h"
#include <felib.h>
#include "source.h"
#include "pparse.h"
#include "rte3D.h"
#define VARYDELTA
#define USE_INTONSPHERE

using namespace toast;

#define MIN(A,B) ( (A) < (B) ? (A) : (B))

class MyDataContext{
public:
double sigma0_0_0; 
ScatKernType sktyp;  
CCompRowMatrix Sint, Sdx, Sdy, Sdz;
CCompRowMatrix Sgrad, Sx, Sy, Sz;
CCompRowMatrix Sdxx, Sdxy, Sdyx, Sdyy,  Sdxz, Sdzx,  Sdyz, Sdzy, Sdzz;
CCompRowMatrix Aint, Aintsc, Aintss, Aintc, Anvec, Anvec_sc, Anvec_ss,  Anvec_c;
CCompRowMatrix Aintscsc,  Aintscss, Aintscc,  Aintssss,  Aintssc,  Aintcc;
CCompRowMatrix SPS, SPSdx, SPSdy, SPSdz;
CCompRowMatrix spatA3_rte, spatA3_sdmx, spatA3_sdmy, spatA3_sdmz;
CCompRowMatrix apu1, apu1sc, apu1ss, apu1c;
CCompRowMatrix A2, b1;
CDenseMatrix Aintx, Aintscx, Aintssx, Aintcx;
CDenseMatrix Aintscscx, Aintscssx, Aintssssx, Aintsccx, Aintsscx, Aintccx;
CDenseMatrix apu1x, apu1scx, apu1ssx, apu1cx;
CVector A2x;
CDenseMatrix A0_RTE, A0_SDM, A1_RTE, A1_SDM, A3, A4, A, Xmat, A0;

const complex *sintval, *sdxval, *sdyval, *sdzval, *sxval, *syval, *szval; 
const complex *sdxxval, *sdxyval, *sdyxval, *sdyyval, *sdxzval, *sdzxval, *sdyzval, *sdzyval, *sdzzval; 
const complex *spsval, *spsdxval, *spsdyval, *spsdzval, *spata3_rteval, *spata3_sdmxval, *spata3_sdmyval, *spata3_sdmzval; 

const complex *aintval, *aintscval, *aintssval, *aintcval, *apu1val, *apu1scval, *apu1ssval, *apu1cval;
const complex *aintscscval, *aintscssval, *aintssssval, *aintsccval, *aintsscval, *aintccval;
complex *xmatval;

int *browptr, *bcolidx, nzb;
int mfc_count;
int spatN, angN;

double w, c;

MyDataContext(QMMesh &spatMesh, Mesh &angMesh, RVector &delta, RVector &muabs, RVector &muscat, ScatKernType &sktyp, double w, double c)
{ 
	this->sktyp = sktyp;
	this->w = w;
	this->c = c;
	spatN = spatMesh.nlen();
	angN = angMesh.nlen();
	
	genmat_spatint_nobf_3D(spatMesh, muabs, muscat, Sint, Sgrad, Sx, Sy, Sz, SPS, spatA3_rte);
        cout<<"Done nobf ..."<<endl;
	genmat_spatint_sdm_nobf_3D(spatMesh, delta, muabs, muscat, Sdx, Sdy, Sdz, Sdxx, Sdxy, Sdyx, Sdyy,  Sdxz, Sdzx, Sdyz, Sdzy, Sdzz, SPSdx, SPSdy, SPSdz, spatA3_sdmx, spatA3_sdmy, spatA3_sdmz);
	genmat_angint_3D(Aint, Aintsc, Aintss, Aintc, Anvec, angMesh);
        cout<<"Calling genmat_angint_sdm_3D ..."<<endl;
    	genmat_angint_sdm_3D(Aintscsc,  Aintscss,   Aintscc,  Aintssss,  Aintssc,  Aintcc,  Anvec_sc, Anvec_ss,  Anvec_c, angMesh);
	cout<<"Calling genmat_apu ..."<<endl;
	genmat_apu(angMesh, Anvec, Anvec_sc, Anvec_ss, Anvec_c, apu1, apu1sc, apu1ss, apu1c, sktyp);
	cout<<"Calling genmat_boundint_3D ..."<<endl;
  	genmat_boundint_3D(spatMesh,  angMesh, A2, b1);
	

	Sint = Sint*complex(0, w/c);
   	Sdx = Sdx*complex(0, w/c);
  	Sdy = Sdy*complex(0, w/c);
   	Sdz = Sdz*complex(0, w/c);

	int angN = angMesh.nlen();
	angMesh.SparseRowStructure (browptr, bcolidx, nzb);
    
	A2x.New(spatN*angN);
	CDenseMatrix temp(angN, spatN), temp2(spatN, angN); 
        
	
	A0 = temp2; A0_RTE = temp2; A0_SDM = temp2; A1_RTE = temp2; A1_SDM = temp2; A3 = temp2; A4 = temp2; A = temp2; Xmat = temp2;


    	Aint.Transpone(); Aintsc.Transpone(); Aintss.Transpone(); Aintc.Transpone(); 
    	Aintscsc.Transpone(); Aintscss.Transpone(); Aintssss.Transpone(); Aintscc.Transpone(); 
    	Aintssc.Transpone(); Aintcc.Transpone(); apu1.Transpone();  apu1sc.Transpone();  apu1ss.Transpone();  
    	apu1c.Transpone(); 	

	Aintx.New(spatN, angN); Aintscx.New(spatN, angN); Aintssx.New(spatN, angN); Aintcx.New(spatN, angN);
	apu1x.New(spatN, angN); apu1scx.New(spatN, angN); apu1ssx.New(spatN, angN); apu1cx.New(spatN, angN);
	Aintscscx.New(spatN, angN); Aintscssx.New(spatN, angN); Aintssssx.New(spatN, angN); Aintsccx.New(spatN, angN);
	Aintsscx.New(spatN, angN); Aintccx.New(spatN, angN);

	Xmat.New(spatN, angN);

	sintval = Sint.ValPtr(); sdxval = Sdx.ValPtr(); sdyval = Sdy.ValPtr();
	sdzval = Sdz.ValPtr(); sxval = Sx.ValPtr(); syval = Sy.ValPtr(); 
	szval = Sz.ValPtr(); sdxxval = Sdxx.ValPtr(); sdxyval = Sdxy.ValPtr(); 
	sdyxval = Sdyx.ValPtr(); sdyyval = Sdyy.ValPtr(); sdxzval = Sdxz.ValPtr(); 
	sdzxval = Sdzx.ValPtr(); sdyzval = Sdyz.ValPtr(); sdzyval = Sdzy.ValPtr(); 
	sdzzval = Sdzz.ValPtr(); spsval = SPS.ValPtr(); spsdxval = SPSdx.ValPtr(); 
	spsdyval = SPSdy.ValPtr(); spsdzval = SPSdz.ValPtr(); spata3_rteval = spatA3_rte.ValPtr(); 
	spata3_sdmxval = spatA3_sdmx.ValPtr(); spata3_sdmyval = spatA3_sdmy.ValPtr(); 
	spata3_sdmzval = spatA3_sdmz.ValPtr();      
	
/*	int *xrowptr = new int[spatN + 1];
	
	xrowptr[0]=0;
	for(int i=1;i < spatN+1; i++)
		xrowptr[i] = xrowptr[i-1] + angN;
	
	int *xcolidx = new int[xrowptr[spatN]];
	int k=0;
	for(int i = 0; i < spatN; i++)
	{
		for(int j=0; j < angN; j++)
		{
			xcolidx[k] = j;
			k++;
		}
	 }

	Xmat.Initialise(xrowptr, xcolidx);
	*/
 
	mfc_count =0;

}

~MyDataContext()
{
  delete []browptr;
  delete []bcolidx;
};

private:

void genmat_spatint_nobf_3D(const Mesh& mesh, const RVector &muabs, const RVector &muscat, CCompRowMatrix& Sint, CCompRowMatrix& Sgrad, CCompRowMatrix& Sx, CCompRowMatrix& Sy, CCompRowMatrix& Sz, CCompRowMatrix &SPS, CCompRowMatrix &spatA3_rte)
{
   int sysdim = mesh.nlen();       // dimensions are size of nodes.

   int *rowptr, *colidx, nzero;
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

void genmat_spatint_sdm_nobf_3D(const Mesh& mesh,  const RVector& delta, const RVector &muabs, const RVector &muscat, CCompRowMatrix& Sdx, CCompRowMatrix& Sdy,  CCompRowMatrix& Sdz, CCompRowMatrix& Sdxx, CCompRowMatrix& Sdxy, CCompRowMatrix& Sdyx, CCompRowMatrix& Sdyy, CCompRowMatrix& Sdxz, CCompRowMatrix& Sdzx, CCompRowMatrix& Sdyz, CCompRowMatrix& Sdzy, CCompRowMatrix& Sdzz, CCompRowMatrix &SPSdx, CCompRowMatrix &SPSdy, CCompRowMatrix &SPSdz, CCompRowMatrix &spatA3_sdmx, CCompRowMatrix &spatA3_sdmy, CCompRowMatrix &spatA3_sdmz)
{
   int sysdim = mesh.nlen();       // dimensions are size of nodes.

   int *rowptr, *colidx, nzero;
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

void genmat_angint_3D(CCompRowMatrix& Aint, CCompRowMatrix& Aintsc, CCompRowMatrix& Aintss, CCompRowMatrix& Aintc, CCompRowMatrix& Anvec, const Mesh& S2mesh)
{
  // this function returns SN X SN matrices of angular shape integrals.

  // S2Mesh must be 3D - should check for that...

  const int& SN =  S2mesh.nlen();       // dimensions are size of nodes.
  const int& SE =  S2mesh.elen();       // number of spherical elements.

  int *angrowptr, *angcolidx, nzero;
  S2mesh.SparseRowStructure (angrowptr, angcolidx, nzero);
  cout<<"Number of non-zeros in the angular mesh: "<<nzero<<endl;
  Aint.New(SN,SN);   // integrals of products of shape functions 
  Aint.Initialise(angrowptr, angcolidx);
  Aintsc.New(SN,SN); // integrals of products of shape functions X Sin t Cos p
  Aintsc.Initialise(angrowptr, angcolidx);
  Aintss.New(SN,SN); // integrals of products of shape functions X Sin t Sin p
  Aintss.Initialise(angrowptr, angcolidx);
  Aintc.New(SN,SN);  // integrals of products of shape functions X Cos t
  Aintc.Initialise(angrowptr, angcolidx);
  Anvec.New(SN,SN);  // and this one...
  Anvec.Initialise(angrowptr, angcolidx);

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
	 Aint(is,js) += S2mesh.elist[el]->IntUnitSphereFF (S2mesh.nlist,i, j);
	 Aintsc(is,js) += S2mesh.elist[el]->IntUnitSpherePFF (S2mesh.nlist,i, j, sx);
	 Aintss(is,js) += S2mesh.elist[el]->IntUnitSpherePFF (S2mesh.nlist,i, j, sy);
	 Aintc(is,js)  += S2mesh.elist[el]->IntUnitSpherePFF (S2mesh.nlist,i, j, sz);
#else
	 Aint(is,js) += S2mesh.elist[el]->IntFF (i, j);
	 Aintsc(is,js) += S2mesh.elist[el]->IntPFF (i, j, sx);
	 Aintss(is,js) += S2mesh.elist[el]->IntPFF (i, j, sy);
	 Aintc(is,js)  += S2mesh.elist[el]->IntPFF (i, j, sz);
#endif
       }
    }
  }

 //cout<<Aintsc<<endl;
  // next logic is all Tanja - I don't completely get it...
  // CHECK carefully with Tanja!!

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
//    cout << "Node " << jj << " is " << Jind << " in element " << el << endl;
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

void genmat_angint_sdm_3D(CCompRowMatrix& Aintscsc,  CCompRowMatrix& Aintscss,  CCompRowMatrix& Aintscc,  CCompRowMatrix& Aintssss,  CCompRowMatrix& Aintssc,  CCompRowMatrix& Aintcc, CCompRowMatrix& Anvec_sc, CCompRowMatrix& Anvec_ss, CCompRowMatrix& Anvec_c, const Mesh& S2mesh)
{
  // this function returns SN X SN matrices of angular shape integrals.

  // S2Mesh must be 3D - should check for that...

  const int& SN =  S2mesh.nlen();       // dimensions are size of nodes.
  const int& SE =  S2mesh.elen();       // number of spherical elements.

  int *angrowptr, *angcolidx, nzero;
  S2mesh.SparseRowStructure (angrowptr, angcolidx, nzero);

  Aintscsc.New(SN,SN); // products of shape functions (Sin t Cos p)^2
  Aintscsc.Initialise(angrowptr, angcolidx);
  Aintscss.New(SN,SN); // products of shape functions (Sin t Cos p)(Sin t Sint )
  Aintscss.Initialise(angrowptr, angcolidx);
  Aintscc.New(SN,SN);  // products of shape functions (Sin t Cos p) Cos t
  Aintscc.Initialise(angrowptr, angcolidx);
  Aintssss.New(SN,SN); // products of shape functions (Sin t Sin p)^2
  Aintssss.Initialise(angrowptr, angcolidx);
  Aintssc.New(SN,SN); //  products of shape functions (Sin t Sin p) Cos t
  Aintssc.Initialise(angrowptr, angcolidx);
  Aintcc.New(SN,SN);  //  products of shape functions ( Cos t )^2
  Aintcc.Initialise(angrowptr, angcolidx);
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
	 Aintscsc(is,js) += S2mesh.elist[el]->IntUnitSpherePFF (S2mesh.nlist, i, j, sxx);
	 Aintscss(is,js) += S2mesh.elist[el]->IntUnitSpherePFF (S2mesh.nlist, i, j, sxy);
	 Aintscc(is,js)  += S2mesh.elist[el]->IntUnitSpherePFF (S2mesh.nlist, i, j, sxz);
	 Aintssss(is,js) += S2mesh.elist[el]->IntUnitSpherePFF (S2mesh.nlist, i, j, syy);
	 Aintssc(is,js)  += S2mesh.elist[el]->IntUnitSpherePFF (S2mesh.nlist, i, j, syz);
	 Aintcc(is,js)   += S2mesh.elist[el]->IntUnitSpherePFF (S2mesh.nlist, i, j, szz);
#else
	 Aintscsc(is,js) += S2mesh.elist[el]->IntPFF (i, j, sxx);
	 Aintscss(is,js) += S2mesh.elist[el]->IntPFF (i, j, sxy);
	 Aintscc(is,js)  += S2mesh.elist[el]->IntPFF (i, j, sxz);
	 Aintssss(is,js) += S2mesh.elist[el]->IntPFF (i, j, syy);
	 Aintssc(is,js)  += S2mesh.elist[el]->IntPFF (i, j, syz);
	 Aintcc(is,js)   += S2mesh.elist[el]->IntPFF (i, j, szz);
#endif
       }
    }
  }



  // next logic is all Tanja - I don't completely get it...
  // CHECK carefully with Tanja!!

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
//    cout << "Node " << jj << " is " << Jind << " in element " << el << endl;
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

void genmat_apu(const Mesh &S2mesh, CCompRowMatrix &Anvec, CCompRowMatrix &Anvec_sc, CCompRowMatrix &Anvec_ss, CCompRowMatrix &Anvec_c, CCompRowMatrix& apu1, CCompRowMatrix& apu1sc, CCompRowMatrix& apu1ss, CCompRowMatrix& apu1c, ScatKernType sktyp)
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
	 			CCompRowMatrix acol = Anvec.Subcols(bb,bb+1);
	 			CCompRowMatrix acolsc = Anvec_sc.Subcols(bb,bb+1);
	 			CCompRowMatrix acolss = Anvec_ss.Subcols(bb,bb+1);
	 			CCompRowMatrix acolc  = Anvec_c.Subcols(bb,bb+1);
	 	 		for(int cc = 0; cc < SN; cc++){
	   				CCompRowMatrix bcol = Anvec.Subcols(cc,cc+1);
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
		/*case MUSINHOMOG_G0 :
			break;
     		case MUSHOMOG_GCONST :
      		{
			CVector apujjvec(SN*SN), apuscjjvec(SN*SN), apussjjvec(SN*SN), apucjjvec(SN*SN);
			for(int bb = 0; bb < SN; bb++){ // not efficient !
	 			CCompRowMatrix acol = Anvec.Subcols(bb,bb+1);
	 			CCompRowMatrix acolsc = Anvec_sc.Subcols(bb,bb+1);
	 			CCompRowMatrix acolss = Anvec_ss.Subcols(bb,bb+1);
	 			CCompRowMatrix acolc  = Anvec_c.Subcols(bb,bb+1);
	 			
				CCompRowMatrix apuii = kron(acol,Anvec);
	 			CCompRowMatrix apuscii = kron(acolsc,Anvec);
	 			CCompRowMatrix apussii = kron(acolss,Anvec);
        		        CCompRowMatrix apucii  = kron(acolc,Anvec); // *********** Surya 08/01/2010 AN ERROR MAY BE: TO REVERT REPLACE WITH THE LINE BELOW ************
	 			//RCompRowMatrix apucii  = kron(acolss,Anvec); 
 	 			CVector sigcol(SN); // assign a column of sigma[el];
     	 			for (int ii = 0; ii < SN; ii++)
	   				sigcol[ii] = complex(sigma[0](ii,bb), 0);
	 			apujjvec +=  apuii*sigcol;
	 			apuscjjvec += apuscii*sigcol;
	 			apussjjvec += apussii*sigcol;
	 			apucjjvec  += apucii*sigcol;
				 int* arowptr;
	 			if( !(arowptr = new int [SN+1]))
	   				cerr << "Memory Allocation error arowptr = new int\n";
	 			int* acolidx;
	 			if( !(acolidx = new int [SN*SN]))
	   				cerr << "Memory Allocation error acolidx = new int\n";
	 			int rp = 0;
	 			for(int ii = 0, cp = 0; ii < SN; ii++) {
	   				arowptr[ii] = rp; rp += SN;
	   				for(int kk = 0; kk < SN; kk++)
	     					acolidx[cp++] = kk;
	 			}
	 			arowptr[SN] = rp;
	 			apu1.Initialise(arowptr,acolidx);
	 			apu1sc = apu1;
	 			apu1ss = apu1;
	 			apu1c  = apu1;
	 			for (int ii = 0; ii < SN; ii++){
	   				for(int kk = 0; kk < SN; kk++) {
	     					apu1(ii,kk) = apujjvec[ii*SN + kk];
	     					apu1sc(ii,kk) = apuscjjvec[ii*SN + kk];
	     					apu1ss(ii,kk) = apussjjvec[ii*SN + kk];
	     					apu1c(ii,kk)  = apucjjvec[ii*SN + kk];
	   				}
	 			}
       				delete []arowptr;
       				delete []acolidx;
       			}
			apu1 *= -1;
			apu1sc *= -1;
			apu1ss *= -1;
			apu1c *= -1;

	 		break;
		}*/	
			
	}
} 

void initialiseA2b1(const Mesh &mesh, const Mesh &S2mesh, CCompRowMatrix& A2, CCompRowMatrix& b1)
{
   /*int el, nodel, i, j, k,is, js;
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
	 }*/
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
/*         }
    }
    
   cout<<"computed 'status' ..."<<endl;
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
   cout<<"spatrowptr and spatcolidx computed ..."<<endl; 
    int ia, ib, ka, kb, ja, jb, idx;
    int va_i, vb_i, v_i;
    int na = sysdim, ma = sysdim, va = spatrowptr[sysdim], col;
    int nb = S2mesh.nlen();
    int *rowptr = new int[na*nb+1];
    for(i = 0; i < na*nb + 1; i++)
    	rowptr[i] = 0;
    
    k=1; 
    for (ia = 0; ia < na; ia++)
    {
	for(j = 0; j < nb; j++)
	{
		for(i = spatrowptr[ia]; i < spatrowptr[ia+1]; i++)
		{
			col = spatcolidx[i];
			rowptr[k] += nb;
		 }

		k++;
        } 
    }
  
    for(i = 1; i < na*nb+1; i++)
	rowptr[i] += rowptr[i-1];
 
    
   int *colidx = new int[rowptr[na*nb]];
   int col_offset;
   k=0;
   for(ia = 0; ia < na; ia++)
   {
	for(j = 0; j < nb; j++)
	{	
		for(i = spatrowptr[ia]; i < spatrowptr[ia+1]; i++)
		{
			col = spatcolidx[i];

			col_offset =0;
			for(int m = 0; m < col; m++) col_offset += nb;

			for(int l = 0; l < nb; l++)
			{
				colidx[k] = col_offset + l;
				k++;
			}
		}
	}
   }
*/

    int ia, ib, ka, kb, ja, jb, i, idx;
    int va_i, vb_i, v_i;
    int snzero, anzero;
    int *srowptr, *scolidx, *arowptr, *acolidx; 
    mesh.SparseRowStructure(srowptr, scolidx, snzero);
    S2mesh.SparseRowStructure(arowptr, acolidx, anzero);

   int sysdim = mesh.nlen();
   int el, nodel, j, k,is, js;
 
   int *status = new int[srowptr[sysdim]];
   for(i=0; i<snzero; i++)
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
		for (int rp = srowptr[is]; rp < srowptr[is+1]; rp++){
        		if (scolidx[rp] == js) status[rp] = 1;
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
    
   cout<<"computed 'status' ..."<<endl;
    int tnzero=0;
    for(i=0; i < snzero; i++)
	if(status[i]) tnzero++; 
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

    int *rowptr = new int[n+1];
    int *colidx = new int[v];

    rowptr[0] = 0;
    i = idx = 0;
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
    /*cout<<"Number of columns in 102420: "<<rowptr[102421] - rowptr[102420]<<endl;
    for(int i=rowptr[102420] ; i< rowptr[102421]; i++)
	cout<<colidx[i]<<" "<<endl;*/
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

void kronsdplus(const int spatrow, const int spatcol, const int angN, const double a_ij, const CCompRowMatrix &B, CCompRowMatrix& C)
{
    complex *Cval = C.ValPtr();
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
	//cout<<"i: "<<row_offset+ib<<" j: "<<col_offset+col<<endl;

	C(row_offset+ib , col_offset+col) = C.Get(row_offset+ib , col_offset+col) + complex(a_ij, 0)*B.Get(ib, col); 	
     }
    }  

    /*for(int ib = 0; ib < angN; ib++)
    {
	 for(int jb=0; jb < angN; jb++)
	 {
		cout<<"i: "<<row_offset+ib<<" j: "<<col_offset+jb<<endl;
		C(row_offset+ib , col_offset+jb) = C.Get(row_offset+ib , col_offset+jb) + complex(a_ij, 0)*B.Get(ib, jb); 	
	}
     }*/	
}

void genmat_boundint_3D(const Mesh& mesh,  const Mesh& S2mesh, CCompRowMatrix& A2, CCompRowMatrix& b1)
  /*      
     produces complete matrix of integrals in space and angle for boundary term
  */
{
  // S2Mesh must be 3D - should check for that...

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
   int *angrowptr, *angcolidx, nzero;
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
  
   CCompRowMatrix  Angbintplus(SN,SN), Angbintminus(SN, SN);
   Angbintplus.Initialise(angrowptr, angcolidx);
   Angbintminus.Initialise(angrowptr, angcolidx);
   RVector f1(SN), f2(SN); 
   // now create matrix, by looping over elements that have a boundary
   for (el = 0; el < mesh.elen(); el++) {
        //cout<<"element number: "<<el<<endl;
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
		ela_ij = mesh.elist[el]->BndIntFFSide (i, j,sd);
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
void genmat_b2_cos(const Mesh &mesh, const Mesh& S2mesh, CCompRowMatrix& Svec, const int Nsource, const RVector& dirVec, const bool is_isotropic)
  /*      
    Function generates the source values vector for FEM of the radiative 
    transfer equation
  */
{
   int el, nodel, i, j, k,is, js;
   int sysdim = mesh.nlen();       // dimensions are size of nodes.
   int SN = S2mesh.nlen();
   int fullsysdim = sysdim*SN;   // full size of angles X space nodes
   
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
        //if(!mesh.elist[el]->HasBoundarySide()) continue;
	for(int sd = 0; sd <  mesh.elist[el]->nSide(); sd++) {
	  // if sd is not a boundary side. skip 
	  //if(!mesh.elist[el]->IsBoundarySide (sd)) continue;
	  
	  //RVector nhat = mesh.ElDirectionCosine(el,sd);
	  //dirMat(0, 0) = -1*nhat[0]; dirMat(0, 1) = -1*nhat[1]; dirMat(0, 2) = -1*nhat[2]; 
	   dirMat(0, 0) = dirVec[0]; dirMat(0, 1) = dirVec[1]; dirMat(0, 2) = dirVec[2];

	  for(int nd = 0; nd < mesh.elist[el]->nSideNode(sd); nd++) {
	    is = mesh.elist[el]->SideNode(sd,nd);
	    js = mesh.elist[el]->Node[is];
	    if(js == Nsource)
	    {
		offset = 0;
		for(int i=0; i < js; i++) offset += SN; 

	    	if(is_isotropic)
		{ 
			for(int i=0; i < SN; i++) Svec(offset + i, 0) = 1.0;
		}
	    	else{
			//cout<<"Source is not isotropic ..."<<endl;
			for(int i=0; i<SN; i++)
			{
			  if((sx[i] == 0) && (sy[i] == -1) && (sz[i] == 0))
			  {
				//cout<<"Imposed nonisotropic source at node: "<<i<<endl;
			  	Svec(offset + i, 0) = 1.0;
			   }
			  
			}		
		}
	    }

	  }
	} // end loop on element sides

   } // end loop on elements

  /*for(int i =0 ; i < sysdim*SN; i++)
	{
	  if (Svec(i, 0) != 0)
		cout<<"i: "<<i<<" "<<Svec(i, 0)<<endl;
	}
*/

   delete []rowptr;
   delete []colidx;
}

// =========================================================================
// global parameters

SourceMode srctp = SRCMODE_NEUMANN;   // source type
double avg_cmua = 1.0, avg_ckappa = 1.0;
ParamParser pp;
QMMesh qmmesh;
NodeList &nlist=qmmesh.nlist;
ElementList &elist=qmmesh.elist;

// =========================================================================
// local prototypes

inline CVector matrixFreeCaller(const CVector& x, void * context);
inline void testspeeds(void* context);

void SelectSourceProfile (int &qtype, double &qwidth, SourceMode &srctp);
void SelectMeasurementProfile (ParamParser &pp, int &mtype, double &mwidth);
void genmat_toastsource3D(CCompRowMatrix* & Source, const Mesh& mesh, const Mesh& S2mesh, const CCompRowMatrix qvec, const int ns, const CCompRowMatrix& b1, const bool is_isotropic, const int angMesh_node);
void genmat_toastsourcevalvector3D_cos(CCompRowMatrix& Svec, const Mesh& mesh, const Mesh& S2mesh,  const CCompRowMatrix qvec, const int iq, const bool is_isotropic, const int angMesh_node);
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

// main routine **************************************************************

int main (int argc, char *argv[])
{
    char cbuf[200];
    int el;

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
    //    elist[0]->Initialise(nlist);
    qmmesh.Setup();

    // set up angular "mesh"

    Mesh S2Mesh;

    ifstream ifm(argv[2]);
    ifm >> S2Mesh;
    ifm.close();
    cout << "Angular  " << S2Mesh.elen() << " elements, " << S2Mesh.nlen()
         << " nodes\n";
    S2Mesh.Setup();
    const int& SN =  S2Mesh.nlen();       // dimensions are size of nodes.

    char file_extn[200];
    cin>>file_extn;
    cout<<"File name prefix: "<<file_extn<<endl;

        //****** sources.
    //****** should read in QM file, but just set something for now.
    cout << "Forming the source\n";
    int ns = 1, nM;
    int *Nsource = new int [ns];
    int    qprof, mprof;   // source/measurement profile (0=Gaussian, 1=Cosine)
    double qwidth, mwidth; // source/measurement support radius [mm]
    CCompRowMatrix qvec, mvec;
    //    Measurement datatype;
    if(argc < 4) {
    //****** should read in QM file, but just set something for now.
    cin>>Nsource[0];
    cout << "source node coords: "<<nlist[Nsource[0]]<<endl;
    }
    else {

    //    logfile << "Reading mesh" << endl;
    cout << "QM file " << argv[3] << endl;
    ifs.open (argv[3]);
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
	for (int j = 0; j < qmmesh.nlen(); j++) 
	  m[j] *= qmmesh.plist[j].C2A();
	mvec.SetRow (i, m);
    }
    cout<<"set measurement vector ..."<<endl;	
    }
    //******** optical parameters
    double freq = 0, hg;
    int is_iso;
    bool is_isotropic;
    RVector dirVec(3);
    
    cin >> hg;  
    cout << "value of g (WARNING !! This g is considered constant throughout the domain): "<<hg<<endl;
    ScatKernType sktyp = (hg == 0.0 ?   MUSHOMOG_G0 :   MUSHOMOG_GCONST);
    cin >> freq; 
    cout << "value for frequency (MHz) : "<<freq<<endl;
    cin>> is_iso;
    is_isotropic = is_iso>0 ? true : false;
    cout<<"Is the source isotropic or not: "<<is_isotropic<<endl;
    cin>>dirVec[0]; cin>>dirVec[1]; cin>>dirVec[2];
    dirVec = dirVec*1.0/length(dirVec);
    cout<<"Source direction vector is: ["<<dirVec[0]<<", "<<dirVec[1]<<", "<<dirVec[2]<<"]"<<endl;
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
    cout << "mua " << muabs[0] << " mus " << muscat[0] <<endl;

    RVector sigmatot(qmmesh.elen());
    RVector sigma(qmmesh.elen()); // to be assigned
    RVector intst;
    //***** parameters 
    cout << "Frequency " << freq << " MHz\n";
    
    //sigma = new RDenseMatrix [qmmesh.elen() ];
    calc_paramdistr_nobf_new_3D(sigma, sigmatot, intst, qmmesh, muscat, muabs,  g, S2Mesh);

    cout << "intst " << intst << endl;
    // "smoothing parameter delta of streamline diffusion modification.
    // Just set to constant for now.

    cout<<"Source node is: "<<nlist[1669]<<endl;
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
	 //double sval =  elist[el]->Size()/((muscat[el]));
         double sval = 1/(5.0*(muabs[el] + muscat[el]));
	 double dist = l2norm(centre - nlist[1669]);
	 double rval = 5.0*sval/dist; 
         delta[el] = MIN(sval,rval);
	 //delta[el] = sval;
#else
	 min = max = delta[0];
#endif 
	 //	 cout << delta[el] << " ";
    }
    cout << "min : " << vmin(delta) << "max : " << vmax(delta) << endl;

    //************ system matrices
    cout << "calculating finite element integrals\n";
    /**/

    MyDataContext ctxt(qmmesh, S2Mesh, delta, muabs, muscat, sktyp, w, c);
    testspeeds(&ctxt);

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
      genmat_source_3D(Source, qmmesh, S2Mesh , Nsource, ns, ctxt.b1);
    else
      genmat_toastsource3D(Source, qmmesh, S2Mesh, qvec, ns, ctxt.b1, is_isotropic, angMesh_node);

  /* CCompRowMatrix b2;
   RVector dirVec(3); 
   dirVec[0] = 0; dirVec[1] = 1.0; dirVec[2] = 0; 
   genmat_b2_cos(qmmesh, S2Mesh, b2, 4630, dirVec, is_isotropic);*/
 
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
      cout << "start source " << j << endl;

      Phi[j].New(fullsysdim);
      Phisum[j].New(sysdim);
      for(int i = 0; i < fullsysdim; i++)
	RHS[i] = Source[j].Get(i,0);// + b2.Get(i, 0);
	//RHS[i] = b2.Get(i, 0);
	
      GMRES(&matrixFreeCaller, &ctxt, RHS, Phi[j], tol, AACP, 100);
      // BiCGSTAB(&matrixFreeCaller, &ctxt, RHS, Phi[j], tol, AACP, 100);
      osRHS << RHS << endl;
      cout << "finished source " << j << endl;
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
    
    char flnmod[300], farg[300];
    strcpy(flnmod, file_extn); strcpy(farg, file_extn);
    strcat(flnmod, "_lnmod.nim"); strcat(farg, "_arg.nim"); 
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

    cout<<"The solver took "<<(double)(end-start)/CLOCKS_PER_SEC<<" seconds"<<endl;

	delete []Phi;
	delete []Phisum;
	delete []Nsource;
	//delete []Source;

}

void setrow(CCompRowMatrix &mat, int r, CVector &rv){
	complex *val = mat.ValPtr();
	int ncols = rv.Dim();
	for(int j=0; j < ncols; j++)
		val[r*ncols + j] = rv[j];
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

inline void dsmatmult(const CDenseMatrix &A, const CCompRowMatrix &B, CDenseMatrix &C)
{
    dASSERT(A.nCols() == B.nRows(), Invalid sizes of matrices);
    int i, j, k, m, ra, ra1, ra2;
    int nar = A.nRows(), nac = A.nCols();
    int angN = B.nCols();
    //CDenseMatrix C(nar, angN);
    const complex *val = B.ValPtr();
    complex *valc = C.data_buffer();
    for(i=0; i < nac; i++){
	ra1 = B.rowptr[i];
        ra2 = B.rowptr[i+1];
	for(ra = ra1; ra<ra2; ra++){
		j = B.colidx[ra];
		for(k=0; k<nar; k++)	 
			valc[k*angN + j] += A.Get(k, i)*val[ra];}}
  //return C;	
}
inline void testspeeds(void * context)
{
   
      MyDataContext *ctxt = (MyDataContext*) context;
    int spatN = ctxt->spatN; 
    int angN = ctxt->angN;
    CVector result(spatN*angN);
    FILE *fid;
    fid = fopen("MatVecMult_speeds.txt", "w");
	
    CVector x(spatN*angN);
    for(int i=0; i < spatN*angN; i++)
	x[i] = 1.0;
    
    complex *val = ctxt->Xmat.data_buffer();
    memcpy (val, x.data_buffer(), angN*spatN*sizeof(complex));

    clock_t startA0 = clock();
    ctxt->Aintx.Zero(); ctxt->Aintscx.Zero(); ctxt->Aintssx.Zero(); ctxt->Aintcx.Zero();
    dsmatmult(ctxt->Xmat, ctxt->Aint, ctxt->Aintx); dsmatmult(ctxt->Xmat, ctxt->Aintsc, ctxt->Aintscx); 
    dsmatmult(ctxt->Xmat, ctxt->Aintss, ctxt->Aintssx); dsmatmult(ctxt->Xmat, ctxt->Aintc, ctxt->Aintcx); 
    ctxt->A0_RTE.Zero(); ctxt->A0_SDM.Zero(); 
    int i, j,k, ra, ra1, ra2;
    for (i = 0; i < ctxt->Sint.nRows(); i++) {
	ra1 = ctxt->Sint.rowptr[i];
	ra2 = ctxt->Sint.rowptr[i+1];
	for (ra = ra1; ra < ra2; ra++) {
	    k = ctxt->Sint.colidx[ra];
	    for(j = 0; j < angN; j++){
		ctxt->A0_RTE(i, j) = ctxt->A0_RTE.Get(i, j) +  ctxt->sintval[ra]*ctxt->Aintx.Get(k, j); 
		ctxt->A0_SDM(i, j) = ctxt->A0_SDM.Get(i, j) + ctxt->sdxval[ra]*ctxt->Aintscx.Get(k, j) + ctxt->sdyval[ra]*ctxt->Aintssx.Get(k, j) + ctxt->sdzval[ra]*ctxt->Aintcx.Get(k, j);
		}
     }
   }
 
    ctxt->A = (ctxt->A0_RTE + ctxt->A0_SDM); 
    clock_t endA0 = clock();
    fprintf(fid, "Time taken by A0 : %f\n", (double)(endA0-startA0)/CLOCKS_PER_SEC);

    clock_t startA1 = clock();
    ctxt->Aintscx.Zero(); ctxt->Aintssx.Zero(); ctxt->Aintcx.Zero();
    ctxt->Aintscscx.Zero(); ctxt->Aintscssx.Zero(); ctxt->Aintssssx.Zero(); ctxt->Aintsccx.Zero();
    ctxt->Aintsscx.Zero(); ctxt->Aintccx.Zero();
    dsmatmult(ctxt->Xmat, ctxt->Aintsc, ctxt->Aintscx); 
    dsmatmult(ctxt->Xmat, ctxt->Aintss, ctxt->Aintssx); dsmatmult(ctxt->Xmat, ctxt->Aintc, ctxt->Aintcx); 
    dsmatmult(ctxt->Xmat, ctxt->Aintscsc, ctxt->Aintscscx); dsmatmult(ctxt->Xmat, ctxt->Aintscss, ctxt->Aintscssx);  
    dsmatmult(ctxt->Xmat, ctxt->Aintssss, ctxt->Aintssssx); dsmatmult(ctxt->Xmat, ctxt->Aintscc, ctxt->Aintsccx);  
    dsmatmult(ctxt->Xmat, ctxt->Aintssc, ctxt->Aintsscx); dsmatmult(ctxt->Xmat, ctxt->Aintcc, ctxt->Aintccx); 
    ctxt->A1_RTE.Zero(); ctxt->A1_SDM.Zero(); 
    for (i = 0; i < ctxt->Sint.nRows(); i++) {
	ra1 = ctxt->Sint.rowptr[i];
	ra2 = ctxt->Sint.rowptr[i+1];
	for (ra = ra1; ra < ra2; ra++) {
	    k = ctxt->Sint.colidx[ra];
	    for(j = 0; j < angN; j++){
		ctxt->A1_RTE(i, j) = ctxt->A1_RTE.Get(i, j) + ctxt->sxval[ra]*ctxt->Aintscx.Get(k, j) + ctxt->syval[ra]*ctxt->Aintssx.Get(k, j) + ctxt->szval[ra]*ctxt->Aintcx.Get(k, j);
		
		ctxt->A1_SDM(i, j) = ctxt->A1_SDM.Get(i, j) + ctxt->sdxxval[ra]*ctxt->Aintscscx.Get(k, j) + (ctxt->sdxyval[ra]+ctxt->sdyxval[ra])*ctxt->Aintscssx.Get(k, j);
		ctxt->A1_SDM(i, j) = ctxt->A1_SDM.Get(i, j) + ctxt->sdyyval[ra]*ctxt->Aintssssx.Get(k, j) + (ctxt->sdyzval[ra]+ctxt->sdzyval[ra])*ctxt->Aintsscx.Get(k, j);
		ctxt->A1_SDM(i, j) = ctxt->A1_SDM.Get(i, j) + ctxt->sdzzval[ra]*ctxt->Aintccx.Get(k, j) + (ctxt->sdxzval[ra]+ctxt->sdzxval[ra])*ctxt->Aintsccx.Get(k, j);
		}
     }
   }
 
    ctxt->A = (ctxt->A1_RTE + ctxt->A1_SDM); 
    clock_t endA1 = clock();
    fprintf(fid, "Time taken by A1 : %f\n", (double)(endA1-startA1)/CLOCKS_PER_SEC);

    clock_t startA2 = clock();
    ctxt->A2.Ax(x, ctxt->A2x);
    clock_t endA2 = clock();
    fprintf(fid, "Time taken by A1 : %f\n", (double)(endA2-startA2)/CLOCKS_PER_SEC);

    clock_t startA3 = clock();
    ctxt->Aintx.Zero(); ctxt->Aintscx.Zero(); ctxt->Aintssx.Zero(); ctxt->Aintcx.Zero();
    dsmatmult(ctxt->Xmat, ctxt->Aint, ctxt->Aintx); dsmatmult(ctxt->Xmat, ctxt->Aintsc, ctxt->Aintscx); 
    dsmatmult(ctxt->Xmat, ctxt->Aintss, ctxt->Aintssx); dsmatmult(ctxt->Xmat, ctxt->Aintc, ctxt->Aintcx); 
    ctxt->A3.Zero(); 
    for (i = 0; i < ctxt->Sint.nRows(); i++) {
	ra1 = ctxt->Sint.rowptr[i];
	ra2 = ctxt->Sint.rowptr[i+1];
	for (ra = ra1; ra < ra2; ra++) {
	    k = ctxt->Sint.colidx[ra];
	    for(j = 0; j < angN; j++){
		ctxt->A3(i, j) = ctxt->A3.Get(i, j) + ctxt->spata3_rteval[ra]*ctxt->Aintx.Get(k, j) + ctxt->spata3_sdmxval[ra]*ctxt->Aintscx.Get(k, j) + ctxt->spata3_sdmyval[ra]*ctxt->Aintssx.Get(k, j) + ctxt->spata3_sdmzval[ra]*ctxt->Aintcx.Get(k, j);
		}
     }
   }
 
    ctxt->A = (ctxt->A3); 
    clock_t endA3 = clock();
    fprintf(fid, "Time taken by A3 : %f\n", (double)(endA3-startA3)/CLOCKS_PER_SEC);

    clock_t startA4 = clock();
    ctxt->apu1x.Zero(); ctxt->apu1scx.Zero(); ctxt->apu1ssx.Zero(); ctxt->apu1cx.Zero();
    dsmatmult(ctxt->Xmat, ctxt->apu1, ctxt->apu1x); dsmatmult(ctxt->Xmat, ctxt->apu1sc, ctxt->apu1scx);  
    dsmatmult(ctxt->Xmat, ctxt->apu1ss, ctxt->apu1ssx); dsmatmult(ctxt->Xmat, ctxt->apu1c, ctxt->apu1cx); 
    ctxt->A4.Zero(); 
    for (i = 0; i < ctxt->Sint.nRows(); i++) {
	ra1 = ctxt->Sint.rowptr[i];
	ra2 = ctxt->Sint.rowptr[i+1];
	for (ra = ra1; ra < ra2; ra++) {
	    k = ctxt->Sint.colidx[ra];
	    for(j = 0; j < angN; j++){
	 ctxt->A4(i, j) = ctxt->A4.Get(i, j) + ctxt->spsval[ra]*ctxt->apu1x.Get(k, j) + ctxt->spsdxval[ra]*ctxt->apu1scx.Get(k, j) + ctxt->spsdyval[ra]*ctxt->apu1ssx.Get(k, j) + ctxt->spsdzval[ra]*ctxt->apu1cx.Get(k, j); 
		}
     }
   }

    ctxt->A = (ctxt->A4); 
    clock_t endA4 = clock();
    fprintf(fid, "Time taken by A4 : %f\n", (double)(endA4-startA4)/CLOCKS_PER_SEC);

   /*     */
       /* 
    
    */

       //ctxt->A3.Zero(); ctxt->A4.Zero();
   
	        /*
	
	       */
//+ (ctxt->A1_SDM - ctxt->A1_RTE) + ctxt->A3 - ctxt->A4;

    /*complex *res  = result.data_buffer();
    complex *arg1 = ctxt->A2x.data_buffer();
    complex *arg2 = ctxt->A.data_buffer();
    for (int i=0; i < spatN*angN; i++)
    	*res++ += *arg1++ + *arg2++;*/ 
   	/*ctxt->A2.Ax(x, ctxt->A2x);
    	complex *res  = result.data_buffer();
    	complex *arg1 = ctxt->A2x.data_buffer();
    	for (int i=0; i < dof; i++)
    		*res++ += *arg1++;*/

    fclose(fid);
    cout<<"Speed check file written"<<endl;

   }


inline CVector matrixFreeCaller(const CVector& x, void * context)
{
   /*MyDataContext *ctxt = (MyDataContext*) context;
	int spatN = ctxt->spatN;
	int angN = ctxt->angN;
	int dof = spatN*angN;
	CVector result(dof);


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

    int scol;
    for(int is = 0; is < spatN; is++)
    {
	for(int ia = 0; ia < ctxt->angN; ia++)
	{
		complex temp(0, 0);
		for(int js = ctxt->Sint.rowptr[is]; js < ctxt->Sint.rowptr[is+1]; js++)
		{
			scol = ctxt->Sint.colidx[js];
			temp += (ctxt->sintval[js] + ctxt->spata3_rteval[js])*aintxval[scol*angN +  ia];
			temp += (ctxt->sdxval[js] + ctxt->spata3_sdmxval[js] - ctxt->sxval[js])*aintscxval[scol*angN +  ia];
			temp += (ctxt->sdyval[js] + ctxt->spata3_sdmyval[js] - ctxt->syval[js])*aintssxval[scol*angN +  ia];
			temp += (ctxt->sdzval[js] + ctxt->spata3_sdmzval[js] - ctxt->szval[js])*aintcxval[scol*angN +  ia];
			temp += ctxt->sdxxval[js]*aintscscxval[scol*angN +  ia];
			temp += (ctxt->sdxyval[js] + ctxt->sdyxval[js])*aintscssxval[scol*angN +  ia];
			temp += ctxt->sdyyval[js]*aintssssxval[scol*angN +  ia];
			temp += (ctxt->sdxzval[js] + ctxt->sdzxval[js])*aintsccxval[scol*angN +  ia];
			temp += (ctxt->sdyzval[js] + ctxt->sdzyval[js])*aintsscxval[scol*angN +  ia];
			temp += ctxt->sdzzval[js]*aintccxval[scol*angN +  ia];
			temp -= ctxt->spsval[js]*apu1xval[scol*angN +  ia];
			temp -= ctxt->spsdxval[js]*apu1scxval[scol*angN +  ia];
			temp -= ctxt->spsdyval[js]*apu1ssxval[scol*angN +  ia];
			temp -= ctxt->spsdzval[js]*apu1cxval[scol*angN +  ia];
		}
		result[is*angN + ia]  = temp;
	} 
	}


 
	ctxt->A2.Ax(x, ctxt->A2x);
    	complex *res  = result.data_buffer();
    	complex *arg1 = ctxt->A2x.data_buffer();
    	for (int i=0; i < dof; i++)
    		*res++ += *arg1++;

    	return result;*/
      MyDataContext *ctxt = (MyDataContext*) context;
    int spatN = ctxt->spatN; 
    int angN = ctxt->angN;
    CVector result(spatN*angN);
  
    complex *val = ctxt->Xmat.data_buffer();
    memcpy (val, x.data_buffer(), angN*spatN*sizeof(complex));
    
    ctxt->Aintx.Zero(); ctxt->Aintscx.Zero(); ctxt->Aintssx.Zero(); ctxt->Aintcx.Zero();
    ctxt->apu1x.Zero(); ctxt->apu1scx.Zero(); ctxt->apu1ssx.Zero(); ctxt->apu1cx.Zero();
    ctxt->Aintscscx.Zero(); ctxt->Aintscssx.Zero(); ctxt->Aintssssx.Zero(); ctxt->Aintsccx.Zero();
    ctxt->Aintsscx.Zero(); ctxt->Aintccx.Zero();
    dsmatmult(ctxt->Xmat, ctxt->Aint, ctxt->Aintx); dsmatmult(ctxt->Xmat, ctxt->Aintsc, ctxt->Aintscx); 
    dsmatmult(ctxt->Xmat, ctxt->Aintss, ctxt->Aintssx); dsmatmult(ctxt->Xmat, ctxt->Aintc, ctxt->Aintcx); 
    dsmatmult(ctxt->Xmat, ctxt->Aintscsc, ctxt->Aintscscx); dsmatmult(ctxt->Xmat, ctxt->Aintscss, ctxt->Aintscssx);  
    dsmatmult(ctxt->Xmat, ctxt->Aintssss, ctxt->Aintssssx); dsmatmult(ctxt->Xmat, ctxt->Aintscc, ctxt->Aintsccx);  
    dsmatmult(ctxt->Xmat, ctxt->Aintssc, ctxt->Aintsscx); dsmatmult(ctxt->Xmat, ctxt->Aintcc, ctxt->Aintccx);  
    dsmatmult(ctxt->Xmat, ctxt->apu1, ctxt->apu1x); dsmatmult(ctxt->Xmat, ctxt->apu1sc, ctxt->apu1scx);  
    dsmatmult(ctxt->Xmat, ctxt->apu1ss, ctxt->apu1ssx); dsmatmult(ctxt->Xmat, ctxt->apu1c, ctxt->apu1cx); 

    ctxt->A2.Ax(x, ctxt->A2x); 

    ctxt->A0_RTE.Zero(); ctxt->A0_SDM.Zero(); ctxt->A1_RTE.Zero(); ctxt->A1_SDM.Zero(); ctxt->A3.Zero(); ctxt->A4.Zero();
    int i, j,k, ra, ra1, ra2;
    for (i = 0; i < ctxt->Sint.nRows(); i++) {
	ra1 = ctxt->Sint.rowptr[i];
	ra2 = ctxt->Sint.rowptr[i+1];
	for (ra = ra1; ra < ra2; ra++) {
	    k = ctxt->Sint.colidx[ra];
	    for(j = 0; j < angN; j++){
		ctxt->A0_RTE(i, j) = ctxt->A0_RTE.Get(i, j) +  ctxt->sintval[ra]*ctxt->Aintx.Get(k, j); 
		ctxt->A0_SDM(i, j) = ctxt->A0_SDM.Get(i, j) + ctxt->sdxval[ra]*ctxt->Aintscx.Get(k, j) + ctxt->sdyval[ra]*ctxt->Aintssx.Get(k, j) + ctxt->sdzval[ra]*ctxt->Aintcx.Get(k, j);

	        ctxt->A1_RTE(i, j) = ctxt->A1_RTE.Get(i, j) + ctxt->sxval[ra]*ctxt->Aintscx.Get(k, j) + ctxt->syval[ra]*ctxt->Aintssx.Get(k, j) + ctxt->szval[ra]*ctxt->Aintcx.Get(k, j);
		
		ctxt->A1_SDM(i, j) = ctxt->A1_SDM.Get(i, j) + ctxt->sdxxval[ra]*ctxt->Aintscscx.Get(k, j) + (ctxt->sdxyval[ra]+ctxt->sdyxval[ra])*ctxt->Aintscssx.Get(k, j);
		ctxt->A1_SDM(i, j) = ctxt->A1_SDM.Get(i, j) + ctxt->sdyyval[ra]*ctxt->Aintssssx.Get(k, j) + (ctxt->sdyzval[ra]+ctxt->sdzyval[ra])*ctxt->Aintsscx.Get(k, j);
		ctxt->A1_SDM(i, j) = ctxt->A1_SDM.Get(i, j) + ctxt->sdzzval[ra]*ctxt->Aintccx.Get(k, j) + (ctxt->sdxzval[ra]+ctxt->sdzxval[ra])*ctxt->Aintsccx.Get(k, j);

		ctxt->A3(i, j) = ctxt->A3.Get(i, j) + ctxt->spata3_rteval[ra]*ctxt->Aintx.Get(k, j) + ctxt->spata3_sdmxval[ra]*ctxt->Aintscx.Get(k, j) + ctxt->spata3_sdmyval[ra]*ctxt->Aintssx.Get(k, j) + ctxt->spata3_sdmzval[ra]*ctxt->Aintcx.Get(k, j);

	        ctxt->A4(i, j) = ctxt->A4.Get(i, j) + ctxt->spsval[ra]*ctxt->apu1x.Get(k, j) + ctxt->spsdxval[ra]*ctxt->apu1scx.Get(k, j) + ctxt->spsdyval[ra]*ctxt->apu1ssx.Get(k, j) + ctxt->spsdzval[ra]*ctxt->apu1cx.Get(k, j); 
		}
     }
   }
 
    ctxt->A = (ctxt->A0_RTE + ctxt->A0_SDM) + (ctxt->A1_SDM - ctxt->A1_RTE) + ctxt->A3 - ctxt->A4;

    complex *res  = result.data_buffer();
    complex *arg1 = ctxt->A2x.data_buffer();
    complex *arg2 = ctxt->A.data_buffer();
    for (int i=0; i < spatN*angN; i++)
    	*res++ += *arg1++ + *arg2++;
    return result;

   
   }

CVector getDiag(void * context)
{
    MyDataContext *ctxt = (MyDataContext*) context;
    int spatN = ctxt->spatN; 
    int angN = ctxt->angN;
    int nDim = spatN*angN;
    CVector result(nDim);
    int arow, brow;
    complex a0_rte, a0_sdm, a1_rte, a1_sdm;
    complex a0, a1, a2, a3, a4;
    double  coeff = ctxt->w/ctxt->c;
	
    for(int j=0; j < nDim; j++)
    {
	arow = j/angN; brow = j%angN;
	a0_rte = ctxt->Aint.Get(brow, brow) * ctxt->Sint.Get(arow, arow);
	a0_sdm =  ctxt->Aintsc.Get(brow, brow)*ctxt->Sdx.Get(arow, arow);
	a0_sdm += ctxt->Aintss.Get(brow, brow) * ctxt->Sdy.Get(arow, arow);
  	a0_sdm += ctxt->Aintc.Get(brow, brow)*ctxt->Sdz.Get(arow, arow);

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

void genmat_toastsource3D(CCompRowMatrix* & Source, const Mesh& mesh, const Mesh& S2mesh,  const CCompRowMatrix qvec, const int ns, const CCompRowMatrix& b1, const bool is_isotropic, const int angMesh_node)
  /*      
    Function generates the source values vector for FEM of the radiative 
    transfer equation
  */
{
   const int& SN =  S2mesh.nlen();       // dimensions are size of nodes.
   const int& SE =  S2mesh.elen();       // number of spherical elements.
   int sysdim = mesh.nlen();       // dimensions are size of nodes.
   int fullsysdim = sysdim*SN;   // full size of angles X space nodes
   if( !(Source = new  CCompRowMatrix [ns]))
          cerr << "Memory Allocation error Source = new  RCompRowMatrix\n";

   CCompRowMatrix Svec;
   for (int i = 0; i < ns; i++) {
     genmat_toastsourcevalvector3D_cos(Svec, mesh,  S2mesh, qvec,i, is_isotropic, angMesh_node);
     Source[i].New(fullsysdim,1);
     b1.AB(Svec,Source[i]);        // Toast weirdness
   }
}


void genmat_toastsourcevalvector3D_cos(CCompRowMatrix& Svec, const Mesh& mesh, 
const Mesh& S2mesh,  const CCompRowMatrix qvec, const int iq, const bool is_isotropic, const int angMesh_node)
  /*      
    Function generates the source values vector for FEM of the radiative 
    transfer equation
  */
{
   int el, nodel, i, j, k,is, js;
   const int& SN =  S2mesh.nlen();       // dimensions are size of nodes.
   const int& SE =  S2mesh.elen();       // number of spherical elements.
   int sysdim = mesh.nlen();       // dimensions are size of nodes.
   int fullsysdim = sysdim*SN;   // full size of angles X space nodes
   int* srp;
   if( !(srp = new int [sysdim+1]))
      cerr << "Memory Allocation error srp = new int\n";
   int* sci;
   if( !(sci = new int [sysdim]))   // structure of boundary element matrix
      cerr << "Memory Allocation error sci = new int\n";
   for(i = 0; i < sysdim; i++) // make it dense for simplicity
      sci[i] = 0;

   Svec.New (fullsysdim,1);

   //   static int srp[3] = {0,1,2};
   //   static int sci[2] = {0,0};
   //   RCompRowMatrix sint(2,1);       // From Tanya's code. Very strange !
   //   sint.Initialise(srp,sci);
   //   sint(0,0) = 1; sint(1,0) = 1;

   int *arp;
   if( !(arp = new int [SN+1]))
      cerr << "Memory Allocation error arp = new int\n";
   for(i = 0; i <= SN; i++) // make it dense for simplicity
      arp[i] = i;

   int *aci;
   if( !(aci = new int [SN]))   // structure of boundary element matrix
      cerr << "Memory Allocation error sci = new int\n";
   for(i = 0; i < SN; i++) // make it dense for simplicity
      aci[i] = 0;

   int *angrowptr, *angcolidx, nzero;
   S2mesh.SparseRowStructure (angrowptr, angcolidx, nzero);

   //cout << "Angles : " << SN << " nodes " << SE << " elements " << nzero << " nonzeros\n";
   RVector sx(SN);    // sample Sin t Cos p on nodes
   RVector sy(SN);    // sample Sin t Sin p on nodes
   RVector sz(SN);    // sample Cos t on nodes


   for (i = 0; i < SN; i++) {  // create the samples on sphere
    const Node& np = S2mesh.nlist[i] ; // do we really need to copy it ?
    double ilen = 1.0/length(np); // normalise by length, just to be sure
    sx[i] = np[0]*ilen;
    sy[i] = np[1]*ilen;
    sz[i] = np[2]*ilen;
    //cout << sx[i] << " " << sy[i] << " " << sz[i] << endl;
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

	  //srp = new int [sysdim+1];   // structure of boundary element matrix
	 
	  //	  cout << "genmat_sourcevalvector_cos : el " << el << " sd " << sd << endl;
          sint.Zero();
	  
	  for(int nd = 0; nd < mesh.elist[el]->nSideNode(sd); nd++) {
	    is = mesh.elist[el]->SideNode(sd,nd);
	    js = mesh.elist[el]->Node[is];
	    sint(js,0) = complex(1.0, 0);
	  }
	  //	  cout << sint << endl;
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
	    //	    cout << f[j] << " ";
	  }
	  //	  cout << endl;
	  for(ela = 0; ela < SE ; ela++){
	    for(ia = 0; ia < S2mesh.elist[ela]->nNode(); ia++) {
	      if ((isa = S2mesh.elist[ela]->Node[ia]) >= SN) continue;
	      for(ja = 0; ja < S2mesh.elist[ela]->nNode(); ja++) {
		if ((jsa = S2mesh.elist[ela]->Node[ja]) >= SN) continue;
		//XS		cout << "updating Angsvec(" << isa << ",0)\n";
#ifdef USE_INTONSPHERE
		if(is_isotropic)
			Angsvec(isa,0) += S2mesh.elist[ela]->IntUnitSpherePFF ( S2mesh.nlist, ia, ja, f);
#else
		if(is_isotropic)
			Angsvec(isa,0) += S2mesh.elist[ela]->IntPFF (ia, ja, f);
#endif
	      }
	    }
	  }
	  if(!is_isotropic)
		Angsvec(angMesh_node, 0) = 1.0;
	  //	  cout << "Angsvec : " << Angsvec << endl;
  	  Svec += kron(sint,Angsvec)*complex(sweight, 0);		
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
