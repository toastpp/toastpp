//==========================================================================
// genmat_sourcevalvector_cos_2p5D.cc       S.Arridge               21.11.06
//
//==========================================================================

#include <stdio.h>
#include <stdlib.h>
#include <stream.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <mathlib.h>
#include <felib.h>
#include "toast.h"
#include "rte3D.h"
#define USE_INTONSPHERE

    using namespace toast;


void genmat_sourcevalvector_cos_3D(CCompRowMatrix& Svec, const Mesh& mesh, const Mesh& S2mesh, const int Nsource)
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

   // now create vector, by looping over elements that have a boundary
   for (el = 0; el < mesh.elen(); el++) {
        if(!mesh.elist[el]->IsNode(Nsource)) continue; // source not in this el
        if(!mesh.elist[el]->HasBoundarySide()) continue;

	for(int sd = 0; sd <  mesh.elist[el]->nSide(); sd++) {
	  // if sd is not a boundary side. skip 
	  if(!mesh.elist[el]->IsBoundarySide (sd)) continue;

	  //srp = new int [sysdim+1];   // structure of boundary element matrix
	  for(i = 0; i < sysdim+1; i++) // make it dense for simplicity
	    srp[i] = i;
	  CCompRowMatrix sint(sysdim,1);
	  sint.Initialise(srp,sci);
	  //cout << "genmat_sourcevalvector_cos : el " << el << " sd " << sd << endl;

	  for(int nd = 0; nd < mesh.elist[el]->nSideNode(sd); nd++) {
	    is = mesh.elist[el]->SideNode(sd,nd);
	    js = mesh.elist[el]->Node[is];
	    sint(js,0) = complex(1.0, 0);
	  }
	  //	  cout << sint << endl;
	  // get boundary normal...
	  RVector nhat = mesh.ElDirectionCosine(el,sd);
	  RVector shat = -nhat; // inward directed vector
	  //	  double ang = atan2(shat[1],shat[0]); // inward angle
	  //	  cout << "el " << el << " sd " << sd << " shat " << shat<< endl;

	  // now do angular integrals
	  CCompRowMatrix  Angsvec(SN,1);
	  Angsvec.Initialise(arp, aci);
	  RVector f(SN);    // sample (shat . nhat) on nodes
	  int ia,ja,isa,jsa,ela;
	  //	  cout << "(s.n)- : ";
	  for( j = 0; j < SN; j++) {
	    double tmp = nhat[0]*sx[j] + nhat[1]*sy[j] + nhat[2]*sz[j]; 
	    f[j] =  (tmp < 0.0 ? -1*tmp : 0.0);
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
		Angsvec(isa,0) += S2mesh.elist[ela]->IntUnitSpherePFF ( S2mesh.nlist, ia, ja, f);
#else
		Angsvec(isa,0) += S2mesh.elist[ela]->IntPFF (ia, ja, f);
#endif
	      }
	    }
	  }
	  //	  cout << "Angsvec : " << Angsvec << endl;
  	  Svec += kron(sint,Angsvec);		
	} // end loop on element sides
   } // end loop on elements
   delete [] sci;
   delete [] srp;
   delete [] aci;
   delete [] arp;
}


