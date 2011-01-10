/***************************************************************************
 * genmat_angint_3D.cc            Simon Arridge                27.11.06    *
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
#include <toast.h>
#include <rte3D.h>

#define USE_INTONSPHERE
    using namespace toast;


 
void genmat_angint_3D(RCompRowMatrix& Aint, RCompRowMatrix& Aintsc, RCompRowMatrix& Aintss, RCompRowMatrix& Aintc, RCompRowMatrix& Anvec, const Mesh& S2mesh)
{
  // this function returns SN X SN matrices of angular shape integrals.

  // S2Mesh must be 3D - should check for that...

  const int& SN =  S2mesh.nlen();       // dimensions are size of nodes.
  const int& SE =  S2mesh.elen();       // number of spherical elements.

  int *angrowptr, *angcolidx, nzero;
  S2mesh.SparseRowStructure (angrowptr, angcolidx, nzero);

  cout << SN << " nodes " << SE << " elements " << nzero << " nonzeros\n";
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
    //    cout << sx[i] << " " << sy[i] << " " << sz[i] << endl;
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
