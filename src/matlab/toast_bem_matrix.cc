#include "mex.h"
#include "mathlib.h"
#include "felib.h"
#include "toastmex.h"
#include "util.h"
#include "bemlib.h"

using namespace std;

void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	//SHOW OUTPUT
	double sout = mxGetScalar (prhs[0]);
	BEM_Surface *bemsurf;
	int i, j, idx;

	if (mxGetM(prhs[1]) > 1) { // pass explicit node and element list

		cerr << "Using node/element matrices to initialise surface" << endl;
		// copy node positions
		RDenseMatrix s_nodes;
		CopyMatrix (s_nodes, prhs[1]);
		if (sout>0.5)
		{    cerr << "# matrix s_nodes: " << s_nodes.nRows() << " x "
		<< s_nodes.nCols() << endl;}

		// copy element connectivity matrix
		IDenseMatrix s_elements;
		CopyMatrix (s_elements, prhs[2]);
		int idxmin = 1000000;
		for (i = 0; i < s_elements.nRows(); i++)
			for (j = 0; j < s_elements.nCols(); j++)
				if (s_elements(i,j) < idxmin) idxmin = s_elements(i,j);
		if (idxmin > 0) {
			cerr << "Warning: index list looks 1-based. Shifting indices to 0" << endl;
			for (i = 0; i < s_elements.nRows(); i++)
				for (j = 0; j < s_elements.nCols(); j++)
					s_elements(i,j) -= 1;
		}

		if (sout >0.5){
			cerr << "# matrix s_elements: " << s_elements.nRows() << " x "
				<< s_elements.nCols() << endl;}


		bemsurf = new BEM_Surface(s_nodes, s_elements);
		idx = 3;
		cerr << "Finished!" << endl;

	} else { // pass FEM volume mesh

		cerr << "Using mesh to initialise surface" << endl;
		Mesh *mesh = (Mesh*)Handle2Ptr (mxGetScalar (prhs[1]));
		bemsurf = new BEM_Surface (*mesh);
		idx = 2;
		cerr << "Finished!" << endl;

	}

	BEM_Region bemreg(bemsurf);
	BEM_Kernel *kernel = new BEM_Kernel_Helmholtz;



	// absorption coeff
	double mua = mxGetScalar (prhs[idx++]);
	if (sout >0.5){
		cerr << "# mua = " << mua << endl;}

	// scattering coeff
	double mus = mxGetScalar (prhs[idx++]);
	if (sout >0.5){
		cerr << "# mus = " << mus << endl;}

	//frequency
	double freq = mxGetScalar (prhs[idx++]);
	if (sout >0.5){
		cerr << "# frequency = " << freq << endl;}

	//Refrective index
	double nu = mxGetScalar (prhs[idx++]);
	if (sout >0.5){
		cerr << "# refracive idx = " << nu << endl;}


	bemreg.SetMua(mua);
	bemreg.SetMus(mus);
	bemreg.SetRef(nu);
	bemreg.SetFreq(freq);
	bemreg.SetKernel(kernel);

	cerr << "finished setup kernel" << endl;

	//bemlib jan10; 

	// call the single region calculation
	CDenseMatrix a, b;
	bemreg.ConstructRegionMatrix (a, b);
	//jan10.single_region_refidx_3D (s_nodes, s_elements, mua, mus, freq, nu, a, b);
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
