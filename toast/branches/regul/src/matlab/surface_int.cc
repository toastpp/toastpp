/*****************************************************************************
 * surf_shape_int_test.cc              Simon Arridge                07.09.05 *
 *                                                                           *
 * integrate the shape function and forward and adjoint solutions on a mesh  *
 * surface. Assume for now that these are read in                            *
 *****************************************************************************/

//#include <time.h>
//#include <fstream>
//#include <iomanip.h>
//#include <stdio.h>
//#include <stdlib.h>
//#include <math.h>
//#include <string.h>
//#include "mathlib.h"
#include "felib.h"
#include "mex.h"
#include "toastmex.h"
#include "stoastlib.h"
#include "fwdsolver.h"
#include "util.h"
using namespace std;

#define MAXNB 1000


//int dimension;			// mesh dimension (2 or 3)
//Environment env;

//QMMesh qmmesh;
//NodeList &nlist=qmmesh.nlist;
//ElementList &elist=qmmesh.elist;

void shape_deriv_mat(QMMesh &qmmesh, RVector& nodalfunc,  
		     RDenseMatrix& shapmatFF, RDenseMatrix& shapmatDD);

// error handler for FE library routines *************************************
#ifdef UNDEF
void LocalErrorhandler (char *msg)
{
    cerr << "\nread_toastbemmesh (PID " << getpid() << ")\n" << msg << endl << flush;
    cerr << "Aborted.\n";
    //    logfile << msg << endl << "Aborted." << endl;
    exit (1);
}
#endif

// main routine **************************************************************
#ifdef UNDEF
int main (int argc, char *argv[])
{
    char cbuf[200];
    char pmdfname[256];

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
    dimension = nlist[0].Dim();
    for (int i = 1; i < nlist.Len(); i++)
	xASSERT(nlist[i].Dim() == dimension, Inconsistent node dimensions.);
    xASSERT(dimension >= 2 && dimension <= 3, Mesh dimension must be 2 or 3.);
    //    elist[0]->Initialise(nlist);
    qmmesh.Setup();

    RVector nodalfunc(qmmesh.nlen()); //  function sampled on nodes.
    RDenseMatrix shapmatFF(qmmesh.nlen(),qmmesh.nlen()); // should be sparse.
    RDenseMatrix shapmatDD(qmmesh.nlen(),qmmesh.nlen()); // should be sparse.

    nodalfunc = 1; // test with a constant function
    shape_deriv_mat(qmmesh, nodalfunc, shapmatFF, shapmatDD);
    
    //    cout << "shape matrix FF : \n" << shapmatFF << endl;
    //    cout << "shape matrix DD : \n" << shapmatDD << endl;

    /* invent "forward and adjoint fields", for testing */
    RVector phi(qmmesh.nlen()), adjphi(qmmesh.nlen());
    RVector Jphi(qmmesh.nlen()), Jadjphi(qmmesh.nlen());

    phi = 1; adjphi = 1; Jphi = 1; Jadjphi = 1; // not realistic, just tests.
    double areaint = adjphi & (shapmatFF*phi);
    double Jareaint = Jadjphi & (shapmatDD*Jphi);

//    cout << "nodal integrals :\n" << nodalint << endl;
    double rad = sqrt(areaint/(4*M_PI));
    cout << "total integral is " << areaint << "\tradius is " << rad << endl;
    double Jrad = sqrt(Jareaint/(4*M_PI));
    cout << "J total integral is "  << Jareaint << endl;
}
#endif

void shape_deriv_mat(QMMesh &qmmesh, RVector& nodalfunc,  
		     RDenseMatrix& shapmatFF,RDenseMatrix& shapmatDD)
/*
     this function takes a QMMesh, and a function defined on the nodes of
     the mesh, and returns the matrices of integrals of products of shape 
     functions, and scale products of gradients of shape functions, over
     the surface of the mesh
*/
{
   NodeList &nlist=qmmesh.nlist;
   ElementList &elist=qmmesh.elist;

   int el, nodel, i, j, is, js;

   for (el  = 0; el < qmmesh.elen(); el++){
	nodel = qmmesh.elist[el]->nNode();


	for (i = 0; i < nodel; i++) {
	  is = qmmesh.elist[el]->Node[i];  // global node number
	  for (j = 0; j < nodel; j++) { // not efficient - use FF matrix
	    js = qmmesh.elist[el]->Node[j];  // global node number

	    shapmatFF(is,js) += qmmesh.elist[el]->IntPFF(i,j,nodalfunc);
	    shapmatDD(is,js) += qmmesh.elist[el]->IntPDD(i,j,nodalfunc);
	  }
	}
    }

}

void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{  
   // int hMesh = (int)mxGetScalar (prhs[0]);
 // QMMesh *qmmesh = (QMMesh*)hMesh;
    QMMesh *qmmesh = (QMMesh*)Handle2Ptr (mxGetScalar (prhs[0]));
   

    RVector nodalfunc;
    CopyVector (nodalfunc, prhs[1]);

    RDenseMatrix shapmatFF(qmmesh->nlen(),qmmesh->nlen()); // should be sparse.
    RDenseMatrix shapmatDD(qmmesh->nlen(),qmmesh->nlen()); // should be sparse.

    shape_deriv_mat(*qmmesh, nodalfunc, shapmatFF,shapmatDD);

    CopyMatrix (&plhs[0], shapmatFF);
    CopyMatrix (&plhs[1], shapmatDD);
}
