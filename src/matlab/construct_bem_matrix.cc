#include "mex.h"
#include "mathlib.h"
#include "felib.h"
#include "toastmex.h"
#include "bemdifflib.h"

using namespace std;

void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  //SHOW OUTPUT
  double sout = mxGetScalar (prhs[0]);


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
  
  // absorption coeff
  double mua = mxGetScalar (prhs[3]);
  if (sout >0.5){
    cerr << "# mua = " << mua << endl;}

  // scattering coeff
  double mus = mxGetScalar (prhs[4]);
  if (sout >0.5){
    cerr << "# mus = " << mus << endl;}

  //frequency
  double freq = mxGetScalar (prhs[5]);
  if (sout >0.5){
    cerr << "# frequency = " << freq << endl;}

  //Refrective index
  double nu = mxGetScalar (prhs[6]);
  if (sout >0.5){
    cerr << "# refracive idx = " << nu << endl;}

  bemdiff jan10; 
  
  // call the single region calculation
  CDenseMatrix a, b;
  jan10.single_region_refidx_3D (s_nodes, s_elements, mua, mus, freq, nu, a, b);
  if (sout >0.5){
    cerr << "= Finished single region calculation" << endl;}

  // copy a and b matrices back to matlab
  CopyMatrix (&plhs[0], a);
  if (sout >0.5){
    cerr << "= Returning matrix a: " << a.nRows() << " x " << a.nCols() << endl;}

  CopyMatrix (&plhs[1], b);
  if (sout >0.5){
    cerr << "= Returning matrix b: " << b.nRows() << " x " << b.nCols() << endl;}

}
