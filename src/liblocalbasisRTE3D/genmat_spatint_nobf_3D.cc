//==========================================================================
// genmat_spatint_nobf_3D.cc         S.Arridge                      27.11.06
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
    using namespace toast;


void genmat_spatint_nobf_3D(const Mesh& mesh, RCompRowMatrix& Sint, RCompRowMatrix& Sgrad, RCompRowMatrix& Sx, RCompRowMatrix& Sy, RCompRowMatrix& Sz, RCompRowMatrix* &SP)
{
   int sysdim = mesh.nlen();       // dimensions are size of nodes.

// reset storage allocation from mesh neighbour list
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
   
   if( !(SP = new RCompRowMatrix [mesh.elen() ]))
//       cerr <<  "Memory Allocation error SP = new RCompRowMatrix\n";
	;
   if( !(elrowptr = new int [sysdim +1])) // holds sparse struture of element matrix
//     cerr << "Memory Allocation error elrowptr = new int\n";
	;

   double elk_ij, elb_ij, elsx_ij, elsy_ij, elsz_ij;
   int el, nodel, i, j, k,is, js;

   for (el = 0; el < mesh.elen(); el++) {
     //     cout << el << " ";
	nodel = mesh.elist[el]->nNode();
 
	// create structure of element matrix
        if( !(elcolidx = new int [nodel * nodel]))
//	  cerr << "Memory Allocation error elcolidx = new int\n";
	;
	//	int *elnlist = new int [nodel];
	IVector elnlist(nodel);
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
	  if(k < nodel && i == elnlist[k]) {
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
	//	delete [] elnlist;

	SP[el].New(sysdim,sysdim);
	SP[el].Initialise (elrowptr, elcolidx);
	// now determine the element integrals
	for (i = 0; i < nodel; i++) {
	    if ((is = mesh.elist[el]->Node[i]) >= sysdim) continue;
	    for (j = 0; j < nodel; j++) {
		if ((js = mesh.elist[el]->Node[j]) >= sysdim) continue;
		elb_ij = mesh.elist[el]->IntFF (i, j);
		Sint(is,js) += elb_ij;
      		SP[el](is,js) += elb_ij;
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
	delete []elcolidx; 
   }


   delete []elrowptr;

   delete []rowptr;
   delete []colidx;

}
