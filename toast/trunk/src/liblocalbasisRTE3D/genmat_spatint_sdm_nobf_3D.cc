//==========================================================================
// genmat_spatint_sdm_nobf_3D.cc         S.Arridge                 27.11.06
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


// generates some sparse matrices of "system" size.
// delta is an input array of parameters pertinent to the streamline diffusion
// functions.

void genmat_spatint_sdm_nobf_3D(const Mesh& mesh,  const RVector& delta, RCompRowMatrix& Sdx, RCompRowMatrix& Sdy,  RCompRowMatrix& Sdz, RCompRowMatrix& Sdxx, RCompRowMatrix& Sdxy, RCompRowMatrix& Sdyx, RCompRowMatrix& Sdyy, RCompRowMatrix& Sdxz, RCompRowMatrix& Sdzx, RCompRowMatrix& Sdyz, RCompRowMatrix& Sdzy, RCompRowMatrix& Sdzz, RCompRowMatrix* &SPdx, RCompRowMatrix* &SPdy, RCompRowMatrix* &SPdz)
{
   cout << "Starting genmat_spatint_sdm_nobf_3D" << endl;
   int sysdim = mesh.nlen();       // dimensions are size of nodes.

// reset storage allocation from mesh neighbour list
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
      
   if( !(SPdx = new RCompRowMatrix [mesh.elen() ]))
     cerr << "Memory Allocation error SPdx = new RCompRowMatrix\n";
   if( !(SPdy = new RCompRowMatrix [mesh.elen() ]))
     cerr << "Memory Allocation error SPdy = new RCompRowMatrix\n";
   if( !(SPdz = new RCompRowMatrix [mesh.elen() ]))
     cerr << "Memory Allocation error SPdz = new RCompRowMatrix\n";
   if( !(elrowptr = new int [sysdim +1])) // holds sparse struture of element matrix
     cerr << "Memory Allocation error elrowptr = new int\n";

   cout << "Memory allocated OK" << endl;
   double elk_ij, elb_ij, elsx_ij, elsy_ij, elsz_ij;
   int el, nodel, i, j, k,is, js;

   for (el = 0; el < mesh.elen(); el++) {
	nodel = mesh.elist[el]->nNode();
  	double dss = delta[el]; //streamline diffusion value for this element.
	// create structure of element matrix
        if(!(elcolidx = new int [nodel * nodel]))
	  cerr << "Memory Allocation error elcolidx = new int\n";
	int *elnlist;
	if( !(elnlist = new int[nodel]))
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
	  if(k < nodel && i == elnlist[k]) {
	    rp = rp + nodel;
	    k++;
	    for(j = 0; j < nodel; j++)
	      elcolidx[cp++] = elnlist[j];
	  }
	}
 
	elrowptr[sysdim] = rp;
	/*
	//	cout << "elnlist : ";
	for(i = 0; i < nodel; i++)
	  cout << elnlist[i] << " ";
	cout << endl;
	//	cout << "elcolidx : ";
	for(i = 0; i < nodel*nodel; i++)
	  cout << elcolidx[i] << " ";
	cout << endl;
	//	cout << "elrowptr : ";
	for(i = 0; i < sysdim+1; i++)
	  cout << elrowptr[i] << " ";
	cout << endl;
	*/
	
	delete [] elnlist;

	SPdx[el].New(sysdim,sysdim);
	SPdx[el].Initialise (elrowptr, elcolidx);
	SPdy[el].New(sysdim,sysdim);
	SPdy[el].Initialise (elrowptr, elcolidx);
	SPdz[el].New(sysdim,sysdim);
	SPdz[el].Initialise (elrowptr, elcolidx);
	// now determine the element matrices
	RSymMatrix eldd = mesh.elist[el]->Intdd(); // "all at once!""
	//      	cout << "eldd \n" << eldd << endl;

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
	       		
      		SPdx[el](is,js) += dss*elsx_ij;
       		SPdy[el](is,js) += dss*elsy_ij;
       		SPdz[el](is,js) += dss*elsz_ij;
		
		//	      	cout << "el " << el << "(" << is << "," << js << ") x:" << dss*elsx_ij << " y:" << dss*elsy_ij  << " z:" << dss*elsz_ij << endl; 
	       
		// the mixed derivative terms are retrived from eldd matrix
	       
		//		cout << eldd(i*3,j*3) << " " << eldd(i*3,j*3+1) << " " << eldd(i*3+1,j*3) << " " << eldd(i*3+1,j*3+1)<< " " << eldd(i*3,j*3+2) << " " << eldd(i*3+2,j*3) << " " << eldd(i*3+1,j*3+2)<< " " << eldd(i*3+2,j*3+1)<< " " << eldd(i*3+2,j*3+2) << " " << endl;		
		       
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
	delete []elcolidx; 
   }

   delete []elrowptr;

   delete []rowptr;
   delete []colidx;

}
