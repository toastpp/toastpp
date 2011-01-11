/***************************************************************************
 * genmat_angint_sdm_3D.cc         Simon Arridge               27.11.06    *
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
#include <stream.h>
#include <iomanip.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <mathlib.h>
#include <felib.h>
#include "rte3D.h"

#define USE_INTONSPHERE
    using namespace toast;


 
void genmat_angint_sdm_3D(RCompRowMatrix& Aintscsc,  RCompRowMatrix& Aintscss,  RCompRowMatrix& Aintscc,  RCompRowMatrix& Aintssss,  RCompRowMatrix& Aintssc,  RCompRowMatrix& Aintcc, RCompRowMatrix& Anvec_sc, RCompRowMatrix& Anvec_ss, RCompRowMatrix& Anvec_c, const Mesh& S2mesh)
{
  // this function returns SN X SN matrices of angular shape integrals.

  // S2Mesh must be 3D - should check for that...

  const int& SN =  S2mesh.nlen();       // dimensions are size of nodes.
  const int& SE =  S2mesh.elen();       // number of spherical elements.

  int *angrowptr, *angcolidx, nzero;
  S2mesh.SparseRowStructure (angrowptr, angcolidx, nzero);

  cout << SN << " nodes " << SE << " elements " << nzero << " nonzeros\n";
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
    cout << sx[i] << " " << sy[i] << " " << sz[i] << " " << sxx[i] << " " << sxy[i] << " " << sxz[i] << " " << syy[i] << " " << syz[i] << " " << szz[i] << endl;
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
