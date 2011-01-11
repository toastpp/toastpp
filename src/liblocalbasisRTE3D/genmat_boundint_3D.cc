//==========================================================================
// genmat_boundint_3D.cc               S.Arridge                    27.11.06
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


void genmat_boundint_3D(RCompRowMatrix& Bint, const Mesh& mesh,  const Mesh& S2mesh)
  /*      
     produces complete matrix of integrals in space and angle for boundary term
  */
{
  // S2Mesh must be 3D - should check for that...

   const int& SN =  S2mesh.nlen();       // dimensions are size of nodes.
   const int& SE =  S2mesh.elen();       // number of spherical elements.
   const int sysdim = mesh.nlen();       // dimensions are size of nodes.
   const int fullsysdim = sysdim*SN;     // full size of angles X space nodes

   Bint.New (fullsysdim, fullsysdim);
   int *elrowptr, *elcolidx;
   if(!(elrowptr = new int [sysdim +1]))// holds sparse struture of element matrix
     cerr << "Memory Allocation error elrowptr = new int\n";

   double ela_ij;
   int el, nodel, i, j, k,is, js;

   // first create structure for angular integrals 
   int *angrowptr, *angcolidx, nzero;
   S2mesh.SparseRowStructure (angrowptr, angcolidx, nzero);

   cout << "Angles : " << SN << " nodes " << SE << " elements " << nzero << " nonzeros\n";
   RVector sx(SN);    // sample Sin t Cos p on nodes
   RVector sy(SN);    // sample Sin t Sin p on nodes
   RVector sz(SN);    // sample Cos t on nodes

   for (i = 0; i < SN; i++) {  // create the samples on sphere
    const Node& np = S2mesh.nlist[i] ; // do we really need to copy it ?
    double ilen = 1.0/length(np); // normalise by length, just to be sure
    sx[i] = np[0]*ilen;
    sy[i] = np[1]*ilen;
    sz[i] = np[2]*ilen;
    cout << sx[i] << " " << sy[i] << " " << sz[i] << endl;
   }

   // now create matrix, by looping over elements that have a boundary
   for (el = 0; el < mesh.elen(); el++) {
        if(!mesh.elist[el]->HasBoundarySide ()) continue;

	nodel = mesh.elist[el]->nNode();
 
	// create structure of element matrix
        if(!(elcolidx = new int [nodel * nodel]))
	    cerr << "Memory Allocation error elcolidx = new int\n";
	int *elnlist;
	if(!(elnlist= new int[nodel]))
	    cerr << "Memory Allocation error elnlist = new int\n";	   
	for(i = 0; i < nodel; i++) {
	  if ((is = mesh.elist[el]->Node[i]) >= sysdim) continue;
	  // use insertion sort to put the indices in order
	  int insertflag = 0;
	  for(j = 0; j < i; j++) {
	    if(is < elnlist[j]) {
	      //shift, insert, skip
	      for(k = i; k > j; k--)
		elnlist[k] = elnlist[k-1];
	      elnlist[j] = is;
	      insertflag = 1;
	      break;
	    }
	  }
	  if(!insertflag) elnlist[i] = is;
	}
	int rp = 0, cp = 0;
	for(i = 0, k=0; i < sysdim; i++){
	  elrowptr[i] = rp;
	  if(k< nodel && i == elnlist[k]) {
	    rp = rp + nodel;
	    k++;
	    for(j = 0; j < nodel; j++)
	      elcolidx[cp++] = elnlist[j];
	  }
	}
	elrowptr[sysdim] = rp;
	/*
	cout << "elnlist : ";
	for(i = 0; i < nodel; i++)
	  cout << elnlist[i] << " ";
	cout << endl;
	cout << "elcolidx : ";
	for(i = 0; i < nodel*nodel; i++)
	  cout << elcolidx[i] << " ";
	cout << endl;
	cout << "elrowptr : ";
	for(i = 0; i < sysdim+1; i++)
	  cout << elrowptr[i] << " ";
	cout << endl;
	*/
	delete [] elnlist;

	// now determine the element integrals

	for(int sd = 0; sd <  mesh.elist[el]->nSide(); sd++) {
	  // if sd is not a boundary side. skip 
	  if(!mesh.elist[el]->IsBoundarySide (sd)) continue;
 	  RCompRowMatrix elBind(sysdim,sysdim); // local "boundary element"
	  elBind.Initialise (elrowptr, elcolidx); 

  
	  // get boundary normal...
	  RVector nhat = mesh.ElDirectionCosine(el,sd);

	  for (i = 0; i < nodel; i++) {
	    if ((is = mesh.elist[el]->Node[i]) >= sysdim) continue;
	    for (j = 0; j < nodel; j++) {
		if ((js = mesh.elist[el]->Node[j]) >= sysdim) continue;
		ela_ij = mesh.elist[el]->BndIntFFSide (i, j,sd);
		elBind(is,js) += ela_ij;
	    }
	  } // end loop on spatial integrals
	  //	  cout << "element " << el << " side " << sd << "\n" << elBind <<endl;
	  // now do angular integrals
	  RCompRowMatrix  Angbint(SN,SN);
	  Angbint.Initialise(angrowptr, angcolidx);
	  RVector f(SN);    // sample (shat . nhat) on nodes
	  int ia,ja,isa,jsa,ela;
//	  cout << "(s.n)+ : ";
	  for( j = 0; j < SN; j++) {
	    double tmp = nhat[0]*sx[j] + nhat[1]*sy[j] + nhat[2]*sz[j]; 
	    f[j] =  (tmp > 0.0 ? tmp : 0.0);
//	    cout << f[j] << " ";
	  }
//	  cout << endl;
	  for(ela = 0; ela < SE ; ela++){
	    for(ia = 0; ia < S2mesh.elist[ela]->nNode(); ia++) {
	      if ((isa = S2mesh.elist[ela]->Node[ia]) >= SN) continue;
	      for(ja = 0; ja < S2mesh.elist[ela]->nNode(); ja++) {
		if ((jsa = S2mesh.elist[ela]->Node[ja]) >= SN) continue;
#ifdef USE_INTONSPHERE
		Angbint(isa,jsa) += S2mesh.elist[ela]->IntUnitSpherePFF (S2mesh.nlist, ia, ja, f);
#else
		Angbint(isa,jsa) += S2mesh.elist[ela]->IntPFF (ia, ja, f);
#endif
	      }
	    }
	  }

//      	  cout << "\nelement " << el << " side " << sd << "\n:" <<Angbint << endl;
  	  Bint = Bint + kron(elBind,Angbint);		
	  //  	  cout << "Full matrix " << Bint << endl;
	} // end loop on element sides
//	cout  << "deleting colidx\n";
	delete []elcolidx; 

   } // end loop on elements

   delete []elrowptr;
   delete []angrowptr;
   delete []angcolidx;
}


