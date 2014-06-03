#include "mex.h"
#include "mathlib.h"
#include "felib.h"
#include "toastmex.h"
#include "util.h"
#include "sharmlib.h"
//using namespace toast;
using namespace std;


void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  //SHOW OUTPUT
  double sout = mxGetScalar (prhs[0]);

cerr << "# matrix s_nodes: " <<sout << endl;

  // copy node positions
  RDenseMatrix s_nodes;
  CopyMatrix (s_nodes, prhs[1]);
  if (sout>0.5)
    {    cerr << "# matrix s_nodes: " << s_nodes.nRows() << " x "
	      << s_nodes.nCols() << endl;}
 
  // copy element connectivity matrix
  IDenseMatrix s_elements;
  CopyMatrix (s_elements, prhs[2]);
  if (sout >0.5){
    cerr << "# matrix s_elements: " << s_elements.nRows() << " x "
	 << s_elements.nCols() << endl;}
  
 
  diffusphere duck2010; 
  
  // call the single region calculation
  int n = s_nodes.nRows();
  RVector a(n), b(n);
  duck2010.DoDif (s_elements, s_nodes, a, b);
 
  if (sout >0.5){
    cerr << "= Finished single region calculation" << endl;}
	//CopyVector (a, prhs[1]);
//	CopyVector (b, prhs[2]);
  // copy a and b matrices back to matlab
 CopyVector (&plhs[0], a);
 // if (sout >0.5){
 //   cerr << "= Returning matrix a: " << a.nRows() << " x " << a.nCols() << endl;}

  CopyVector (&plhs[1], b);
 // if (sout >0.5){
 //   cerr << "= Returning matrix b: " << b.nRows() << " x " << b.nCols() << endl;}

}


