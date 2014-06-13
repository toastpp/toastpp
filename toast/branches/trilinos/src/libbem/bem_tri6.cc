#define BEMLIB_IMPLEMENTATION

#include "bemlib.h"
#include "bem_tri6.h"

// Gaussian quadrature points

const int n_gauss6=7;

const double shapf0[n_gauss6] ={-0.11111111111111, -0.02807494322308, -0.02807494322308, -0.05258390110255, -0.08076859419189, -0.08076859419189, 0.47435260858554};
const double shapf1[n_gauss6] ={ 0.44444444444444,  0.11229977289232,  0.88413424176407,  0.11229977289232,  0.32307437676755,  0.04103582626314, 0.32307437676755};
const double shapf2[n_gauss6] ={-0.11111111111111, -0.05258390110255, -0.02807494322308, -0.02807494322308,  0.47435260858554, -0.08076859419189, -0.08076859419189};
const double shapf3[n_gauss6] ={ 0.44444444444444,  0.11229977289232,  0.11229977289232,  0.88413424176407,  0.32307437676755,  0.32307437676755,  0.04103582626314};
const double shapf4[n_gauss6] ={-0.11111111111111, -0.02807494322308, -0.05258390110255, -0.02807494322308, -0.08076859419189,  0.47435260858554, -0.08076859419189};
const double shapf5[n_gauss6] ={0.44444444444444,   0.88413424176407,  0.11229977289232,  0.11229977289232,  0.04103582626314,  0.32307437676755,  0.32307437676755};

const double shapd0[n_gauss6] ={ -0.333333333333, -0.880568256420, -0.880568256420,  0.761136512841,  0.594853970706,  0.594853970706, -2.189707941412 };
const double shapd1[n_gauss6] ={  0.000000000000,  1.641704769261,  0.000000000000, -1.641704769261, -2.784561912119,  0.000000000000,  2.784561912119 };
const double shapd2[n_gauss6] ={  0.333333333333, -0.761136512841,  0.880568256420,  0.880568256420,  2.189707941412, -0.594853970706, -0.594853970706 };
const double shapd3[n_gauss6] ={  1.333333333333,  1.880568256420,  0.238863487159,  1.880568256420,  0.405146029294,  3.189707941412,  0.405146029294 };
const double shapd4[n_gauss6] ={  0.000000000000,  0.000000000000,  0.000000000000,  0.000000000000,  0.000000000000,  0.000000000000,  0.000000000000 };
const double shapd5[n_gauss6] ={ -1.333333333333, -1.880568256420, -0.238863487159, -1.880568256420, -0.405146029294, -3.189707941412, -0.405146029294 };
const double shapd6[n_gauss6] ={ -0.333333333333, -0.880568256420, -0.880568256420,  0.761136512841,  0.594853970706,  0.594853970706, -2.189707941412 };
const double shapd7[n_gauss6] ={ -1.333333333333, -0.238863487159, -1.880568256420, -1.880568256420, -3.189707941412, -0.405146029294, -0.405146029294 };
const double shapd8[n_gauss6] ={  0.000000000000,  0.000000000000,  0.000000000000,  0.000000000000,  0.000000000000,  0.000000000000,  0.000000000000};
const double shapd9[n_gauss6] ={  1.333333333333,  0.238863487159,  1.880568256420,  1.880568256420,  3.189707941412,  0.405146029294,  0.405146029294};
const double shapd10[n_gauss6] ={ 0.333333333333,  0.880568256420, -0.761136512841,  0.880568256420, -0.594853970706,  2.189707941412, -0.594853970706 };
const double shapd11[n_gauss6] ={ 0.000000000000,  0.000000000000,  1.641704769261, -1.641704769261,  0.000000000000, -2.784561912119,  2.784561912119 };

const double wg6[n_gauss6]={0.225/2.,0.13239415278850619/2.,0.13239415278850619/2.,
                          0.13239415278850619/2.,
		          0.12593918054482713/2.,0.12593918054482713/2.,
		          0.12593918054482713/2.};

const int n_gauss_o=6;
const double gpts2=0.93246951420315202781;
const double gpts1=0.66120938646626451366;
const double gpts0=0.23861918608319690863;

const double wgpts2=0.17132449237917034504;
const double wgpts1=0.36076157304813860757;
const double wgpts0=0.46791393457269104739;

const double xg_k_1[6]={-gpts2,-gpts1,-gpts0, gpts0, gpts1, gpts2};
const double xg_k_2[6]={-gpts2,-gpts1,-gpts0, gpts0, gpts1, gpts2};

const double wg_k[6]={wgpts2,wgpts1,wgpts0,wgpts0,wgpts1,wgpts2};


BEM_Triangle6::BEM_Triangle6 (BEM_Surface *s, IVector &ndidx): BEM_Element (s, ndidx)
{
	Initialise (ndidx);
	qj = 0;
}

BEM_Triangle6::~BEM_Triangle6()
{
	if (qj) delete []qj;
}

CVector BEM_Triangle6::Integrate_Singular (BEM_Kernel *kernel, int nodep, const Point3D &load, bool invert)
{
	const int dimension = 6;
	const int dim2 = dimension*2;
	const int i_kern = 4;

	int icd; 
	double det_J_r;
	std::complex<double> kern1, kern2;
	Point2D loc;

	static RVector jac_tab(dimension+1), shapf(dimension);
	static CVector kern_tab(i_kern);
	CVector a_aux(dim2);
	//==================================================================
	//  purpose: Singular integrals in triangular element
	//  element is divided by 2 sub-triangles
	//==================================================================

	for (int ic = 0; ic < dimension; ic++) {

		icd = ic + dimension;

		//  node 0------------------------------------------------
		if (nodep == node[0]) { 
			for(int ig=0; ig<n_gauss_o; ig++){

				det_J_r=(1.+xg_k_1[ig])/8.;

				for(int jg=0; jg<n_gauss_o; jg++){

					// Integral t1

					loc[0]=(1.+xg_k_1[ig])*(1.-xg_k_2[jg])/4.; 
					loc[1]=(1.+xg_k_1[ig])*(1.+xg_k_2[jg])/4.;
					shapf = ShapeF (loc);

					jac_tab=Jacobi(loc);//, node);

					//	lomef<< jac_tab<<" : "<<ielemq<<" "<< ksi1<<" "<< ksi2<<endl;
					kern_tab = kernel->Calculate (this, loc, load);

					kern1 = kern_tab[0]*jac_tab[1]+kern_tab[1]*jac_tab[2]+
						kern_tab[2]*jac_tab[3];//Simon!

					a_aux[ic] += kern1*wg_k[ig]*wg_k[jg]*shapf[ic]*
						jac_tab[0]*det_J_r;//Simon!

					a_aux[icd]+= kern_tab[3]*wg_k[ig]*wg_k[jg]*shapf[ic]*
						jac_tab[0]*det_J_r;

					//a_aux[icd]+= wg_k[ig]*wg_k[jg]*jac_tab[0]*det_J_r;

					//fileout<<"ielemq= "<<ielemq<<"\tic= "<<ic<<"\ta_aux[icd]= "<<
					//  a_aux[icd]<<endl;
				}
			}
		}

		//node 1-------------------------------------------------
		if(nodep==node[1]){ 
			for(int ig=0; ig<n_gauss_o; ig++){

				det_J_r=(1.+xg_k_1[ig])/16.;

				for(int jg=0; jg<n_gauss_o; jg++){

					// Integral t1

					loc[0]=(1.-xg_k_1[ig])/4.;
					loc[1]=(1.+xg_k_1[ig])*(1.-xg_k_2[jg])/4.;
					shapf=ShapeF(loc);    

					jac_tab=Jacobi(loc);//, node);

					// differenct kernels for greens / ma and ms 
					kern_tab = kernel->Calculate (this, loc, load);

					kern1 = kern_tab[0]*jac_tab[1]+kern_tab[1]*jac_tab[2]+
						kern_tab[2]*jac_tab[3];//Simon!

					a_aux[ic] += kern1*wg_k[ig]*wg_k[jg]*shapf[ic]*
						jac_tab[0]*det_J_r;//Simon!

					a_aux[icd]+= kern_tab[3]*wg_k[ig]*wg_k[jg]*shapf[ic]*
						jac_tab[0]*det_J_r;
					//a_aux[icd]+= wg_k[ig]*wg_k[jg]*jac_tab[0]*det_J_r;


					// Integral t2    
					loc[0]=(2.-(1.+xg_k_1[ig])*xg_k_2[jg])/4.;
					loc[1]=(1.+xg_k_1[ig])*(1.+xg_k_2[jg])/4.;
					shapf = ShapeF (loc);

					jac_tab=Jacobi(loc);//, node);
					kern_tab = kernel->Calculate (this, loc, load);

					kern1 = kern_tab[0]*jac_tab[1]+kern_tab[1]*jac_tab[2]+
						kern_tab[2]*jac_tab[3];//Simon!

					a_aux[ic] += kern1*wg_k[ig]*wg_k[jg]*shapf[ic]*
						jac_tab[0]*det_J_r;//Simon!

					a_aux[icd]+= kern_tab[3]*wg_k[ig]*wg_k[jg]*shapf[ic]*
						jac_tab[0]*det_J_r;
					//a_aux[icd]+= wg_k[ig]*wg_k[jg]*jac_tab[0]*det_J_r;


					//fileout<<"ielemq= "<<ielemq<<"\tic= "<<ic<<"\ta_aux[icd]= "<<
					//  a_aux[icd]<<endl;

				}
			}
		}
		//node 2-------------------------------------------------
		if(nodep==node[2]){  
			for(int ig=0; ig<n_gauss_o; ig++){

				det_J_r=(1.+xg_k_1[ig])/8.;

				for(int jg=0; jg<n_gauss_o; jg++){

					// Integral t1

					loc[0]=(1.-xg_k_1[ig])/2.; 
					loc[1]=(1.+xg_k_1[ig])*(1.-xg_k_2[jg])/4.; 
					shapf = ShapeF (loc);

					jac_tab=Jacobi(loc);//, node);
					kern_tab = kernel->Calculate (this, loc, load);

					kern1 = kern_tab[0]*jac_tab[1]+kern_tab[1]*jac_tab[2]+
						kern_tab[2]*jac_tab[3];//Simon!

					a_aux[ic] += kern1*wg_k[ig]*wg_k[jg]*shapf[ic]*
						jac_tab[0]*det_J_r;//Simon!

					a_aux[icd]+= kern_tab[3]*wg_k[ig]*wg_k[jg]*shapf[ic]*
						jac_tab[0]*det_J_r;
					//a_aux[icd]+= wg_k[ig]*wg_k[jg]*jac_tab[0]*det_J_r;

					//fileout<<"ielemq= "<<ielemq<<"\tic= "<<ic<<"\ta_aux[icd]= "<<
					// a_aux[icd]<<endl;
				}
			}
		}

		//node 3-------------------------------------------------
		if(nodep==node[3]){ 
			for(int ig=0; ig<n_gauss_o; ig++){

				det_J_r=(1.+xg_k_1[ig])/16.;

				for(int jg=0; jg<n_gauss_o; jg++){

					// Integral t1

					loc[0]=(2.+(1.+xg_k_1[ig])*xg_k_2[jg])/4.; 
					loc[1]=(1.-xg_k_1[ig])/4.; 
					shapf = ShapeF (loc);

					jac_tab=Jacobi(loc);//, node);
					kern_tab = kernel->Calculate (this, loc, load);

					kern1 = kern_tab[0]*jac_tab[1]+kern_tab[1]*jac_tab[2]+
						kern_tab[2]*jac_tab[3];//Simon!

					a_aux[ic] += kern1*wg_k[ig]*wg_k[jg]*shapf[ic]*
						jac_tab[0]*det_J_r;//Simon!

					a_aux[icd]+= kern_tab[3]*wg_k[ig]*wg_k[jg]*shapf[ic]*
						jac_tab[0]*det_J_r;
					//a_aux[icd]+= wg_k[ig]*wg_k[jg]*jac_tab[0]*det_J_r;

					// Integral t2    

					loc[0]=(1.-xg_k_1[ig])/4.; 
					loc[1]=(2.-(1.+xg_k_1[ig])*xg_k_2[jg])/4.; 
					shapf = ShapeF (loc);

					jac_tab=Jacobi(loc);//, node);
					kern_tab = kernel->Calculate (this, loc, load);

					kern1 = kern_tab[0]*jac_tab[1]+kern_tab[1]*jac_tab[2]+
						kern_tab[2]*jac_tab[3];//Simon!

					a_aux[ic] += kern1*wg_k[ig]*wg_k[jg]*shapf[ic]*
						jac_tab[0]*det_J_r;//Simon!

					a_aux[icd]+= kern_tab[3]*wg_k[ig]*wg_k[jg]*shapf[ic]*
						jac_tab[0]*det_J_r;
					//a_aux[icd]+= wg_k[ig]*wg_k[jg]*jac_tab[0]*det_J_r;


					//fileout<<"ielemq= "<<ielemq<<"\tic= "<<ic<<"\ta_aux[icd]= "<<
					//  a_aux[icd]<<endl;
				}
			}
		}

		//node 4-------------------------------------------------
		if(nodep==node[4]){  
			for(int ig=0; ig<n_gauss_o; ig++){

				det_J_r=(1.+xg_k_1[ig])/8.;

				for(int jg=0; jg<n_gauss_o; jg++){

					// Integral t1

					loc[0]=(1.+xg_k_1[ig])*(1.+xg_k_2[jg])/4.;
					loc[1]=(1.-xg_k_1[ig])/2.; 
					shapf = ShapeF (loc);

					jac_tab=Jacobi(loc);//, node);
					kern_tab = kernel->Calculate (this, loc, load);

					kern1 = kern_tab[0]*jac_tab[1]+kern_tab[1]*jac_tab[2]+
						kern_tab[2]*jac_tab[3];//Simon!

					a_aux[ic] += kern1*wg_k[ig]*wg_k[jg]*shapf[ic]*
						jac_tab[0]*det_J_r;//Simon!

					a_aux[icd]+= kern_tab[3]*wg_k[ig]*wg_k[jg]*shapf[ic]*
						jac_tab[0]*det_J_r;
					//a_aux[icd]+= wg_k[ig]*wg_k[jg]*jac_tab[0]*det_J_r;

					//fileout<<"ielemq= "<<ielemq<<"\tic= "<<ic<<"\ta_aux[icd]= "<<
					// a_aux[icd]<<endl;
				}
			}
		}

		//node 5-------------------------------------------------
		if(nodep==node[5]){  
			for(int ig=0; ig<n_gauss_o; ig++){

				det_J_r=(1.+xg_k_1[ig])/16.;

				for(int jg=0; jg<n_gauss_o; jg++){

					// Integral t1

					loc[0]=(1.+xg_k_1[ig])*(1.+xg_k_2[jg])/4.;
					loc[1]=(1.-xg_k_1[ig])/4.; 
					shapf = ShapeF (loc);

					jac_tab=Jacobi(loc);//, node);
					kern_tab = kernel->Calculate (this, loc, load);

					kern1 = kern_tab[0]*jac_tab[1]+kern_tab[1]*jac_tab[2]+
						kern_tab[2]*jac_tab[3];//Simon!

					a_aux[ic] += kern1*wg_k[ig]*wg_k[jg]*shapf[ic]*
						jac_tab[0]*det_J_r;//Simon!

					a_aux[icd]+= kern_tab[3]*wg_k[ig]*wg_k[jg]*shapf[ic]*
						jac_tab[0]*det_J_r;
					//a_aux[icd]+= wg_k[ig]*wg_k[jg]*jac_tab[0]*det_J_r;

					// Integral t2 

					loc[0]=(1.+xg_k_1[ig])*(1.-xg_k_2[jg])/4.;
					loc[1]=(2.+(1.+xg_k_1[ig])*xg_k_2[jg])/4.; 
					shapf = ShapeF (loc);

					jac_tab=Jacobi(loc);//, node);
					kern_tab = kernel->Calculate (this, loc, load);

					kern1 = kern_tab[0]*jac_tab[1]+kern_tab[1]*jac_tab[2]+
						kern_tab[2]*jac_tab[3];//Simon!

					a_aux[ic] += kern1*wg_k[ig]*wg_k[jg]*shapf[ic]*
						jac_tab[0]*det_J_r;//Simon!

					a_aux[icd]+= kern_tab[3]*wg_k[ig]*wg_k[jg]*shapf[ic]*
						jac_tab[0]*det_J_r;
					//a_aux[icd]+= wg_k[ig]*wg_k[jg]*jac_tab[0]*det_J_r;

					//fileout<<"ielemq= "<<ielemq<<"\tic= "<<ic<<"\ta_aux[icd]= "<<
					//  a_aux[icd]<<endl;

				}
			}
		}
	}

	//-------------------------------------------------
	return a_aux;
}

CVector BEM_Triangle6::Integrate_Nonsingular (BEM_Kernel *kernel, const Point3D &load, bool invert)
{
	const int dimension = 6;
	const int dim2 = dimension*2;
	const int i_kern = 4;
	int ic, icd;
	double  a_wsj;
	std::complex<double> kern1;
	static RVector jac_tab(dimension+1);//, shapf(dimension);
	static CVector kern_tab(i_kern);
	CVector a_aux(dim2);
	//==================================================================
	// purpose: to calculate nonsingular integrals 
	//  node is in local numerating system
	//==================================================================

	// loop over quadrature points
	for(int ig=0; ig<n_gauss6; ++ig) {

		// weights for quadrature point ig (for all nodes)
		RVector &qshapf = QuadratureShapeF (ig);
		kern_tab = kernel->Calculate (this, QuadraturePoint(ig), load);
		RVector jac_tab = QuadratureJacobi(ig);
		// Calculate the kernel with the precalculated shapef
		double det_jac1 = wg6[ig]*jac_tab[0];
		double x_nm1 = jac_tab[1];
		double y_nm1 = jac_tab[2];
		double z_nm1 = jac_tab[3];
		kern1 = kern_tab[0]*x_nm1 + kern_tab[1]*y_nm1 + kern_tab[2]*z_nm1;

		// loop over nodes
		for (ic = 0; ic < dimension; ic++) {

			icd=ic+dimension;

			//for the left hand side
			a_wsj=qshapf[ic]*det_jac1;//already multiplied by wg[ig]
			a_aux[ic] += kern1*a_wsj;

			//for the right hand side
			a_aux[icd]+= kern_tab[3]*a_wsj;
		}
	}
	return a_aux;
}

RVector BEM_Triangle6::ShapeF (Point2D &loc) const
{
	// WARNING: This assumes Jan's node ordering!
	// Change for toast node order

	double ksi1, ksi2, ksi3;
	static RVector shapf(6);
	ksi1 = loc[0], ksi2 = loc[1];
	ksi3=1.-ksi1-ksi2;
	// calculate the quadratic shape functions
	shapf[0]=-ksi3*(1.-2.*ksi3);  shapf[1]=4.*ksi1*ksi3;
	shapf[2]=-ksi1*(1.-2.*ksi1);  shapf[3]=4.*ksi1*ksi2;
	shapf[4]=-ksi2*(1.-2.*ksi2);  shapf[5]=4.*ksi2*ksi3;
  
    return shapf;
}


RVector BEM_Triangle6::ShapeD (Point2D &loc) const
{
	// WARNING: This assumes Jan's node ordering!
	// Change for toast node order

	double ksi1, ksi2, ksi3, ksi1_4, ksi2_4, ksi3_4;
    static RVector shapd(6*2);

	ksi1 = loc[0], ksi2 = loc[1];
    ksi3=1.-ksi1-ksi2; ksi3_4=ksi3*4.;
    ksi1_4=ksi1*4.;  ksi2_4=ksi2*4.;
  
    //calculate the derivatives of the shape functions with respect to ksi1
    shapd[0] = 1.-ksi3_4; shapd[1] = ksi3_4-ksi1_4; shapd[2] =-1.+ksi1_4;
    shapd[3] = ksi2_4;    shapd[4] = 0.;            shapd[5] =-ksi2_4;
  
    //calculate the derivatives of the shape functions with respect to ksi2
    shapd[6] = 1.-ksi3_4; shapd[7] =-ksi1_4;        shapd[8] = 0.;
    shapd[9] = ksi1_4;    shapd[10]=-1.+ksi2_4;     shapd[11]= ksi3_4-ksi2_4;
  
    return shapd;
}


Point2D &BEM_Triangle6::QuadraturePoint (int i) const
{
	// returns coordinates of quadrature point i

	const double beta1=0.47014206410511508;
	const double beta2=0.10128650732345633;
	const double alfa1=0.0597158717897698;
	const double alfa2=0.7974269853530873;

	const double xg_1[7]={1./3.,alfa1,beta1,beta1,alfa2,beta2,beta2};
	const double xg_2[7]={1./3.,beta1,alfa1,beta1,beta2,alfa2,beta2};
	//const double xg_3[7]={1./3.,beta1,beta1,alfa1,beta2,beta2,alfa2};

	static bool need_setup = true;
	static Point2D qp[n_gauss6];
	if (need_setup) {
		for (int j = 0; j < n_gauss6; j++) {
			qp[j][0] = xg_1[j];
			qp[j][1] = xg_2[j];
		}
		need_setup = false;
	}
	return qp[i];
}

RVector &BEM_Triangle6::QuadratureShapeF (int i) const
{
	// returns values of shape functions at quadrature point i
	static RVector qf[n_gauss6];
	static bool need_setup = true;
	if (need_setup) {
		for (int j = 0; j < n_gauss6; j++)
			qf[j] = ShapeF (QuadraturePoint (j));
		need_setup = false;
	}
	return qf[i];
}

RVector &BEM_Triangle6::QuadratureShapeD (int i) const
{
	// returns values of shape derivatives at quadrature point i
	static RVector qd[n_gauss6];
	static bool need_setup = true;
	if (need_setup) {
		for (int j = 0; j < n_gauss6; j++)
			qd[j] = ShapeD (QuadraturePoint (j));
		need_setup = false;
	}
	return qd[i];
}

RVector &BEM_Triangle6::QuadratureJacobi (int i) const
{
	// returns values of Jacobi matrix at quadrature point i
	if (!qj) {
		qj = new RVector[n_gauss6];
		for (int j = 0; j < n_gauss6; j++)
			qj[j] = Jacobi (QuadraturePoint (j));
	}
	return qj[i];
}
