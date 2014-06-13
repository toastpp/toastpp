#define BEMDIFFLIB_IMPLEMENTATION

#include <iostream>
#include <fstream>
#include <math.h>
#include <iomanip>
#include "mathlib.h"
#include "bem_lib_3D.h"
#include "bemdiff.h"

bemdiff::bemdiff()
{
}

bemdiff::~bemdiff()
{
}


void bemdiff::single_region_ms(RDenseMatrix &S_nodes, IDenseMatrix &S_elements, 
			       double mua, double mus, double freq, CDenseMatrix &a, CDenseMatrix &b)
{

  // #################  Optical Properties ##############################################;
  
  //Initialise optical parameters vectors;
  D.New(1);                                      // Diffusion Coefficient
  k.New(1);                                      // Wave Number
  //used for the dk/dmua
  ma.New(1);
  ms.New(1);   
  iwc.New(1);                              
  const double nu=1.;                            // Refractive Index (to be included)
  const double c=0.3*1.e+12/nu;                  // Speed of light in the medium in mm/ps
  double omega=2.*pi*freq;                       // Angular speed
  D[0] = 1./(3.*(mua+mus));                      // Diffusion coefficient
  k[0] = sqrt(std::complex<double>(mua*c, -omega)/(c*D[0]));  // Wave number
  ma[0]=mua;
  ms[0]=mus;
  iwc[0]=std::complex<double>(0,omega/c);
  
  // ############## The Geometrical  definitions ########################################;  
 
  
  n_elements = S_elements.nRows();               // Number of Elements in the  region
  n_nodes    = S_nodes.nRows();                  // Number of Nodes in the  region
  node = S_elements;                             // node - integer matrix stored as an n_elements by dimension
   
  x=S_nodes.Col(0);                              // The Coordinates of the nodes 
  y=S_nodes.Col(1);
  z=S_nodes.Col(2); 

  // ################ The Jan definitions ###############################################;  
  dimension  = 6;                                // Dimension of the Triangles 
  dim2       = 12;                               // number of derivatives of Shape functions 
  i_kern     = 4;                                // length of the kernel function return
  int i_col=0;
  CVector a_aux(dim2);
  int icd,  nodeq, norder[6];
  
 
  double xp, yp, zp;
 
  a.New(n_nodes,n_nodes);                        // a matrix for the return    
  b.New(n_nodes,n_nodes);                        // b matrix for the return
 
  // ############### Precalculate the Jacobian for nonsin #################################;
 
  RVector jac_tab(dimension+1);
 
  det_jac1.New(n_elements,n_gauss);             // Determinant 
  x_nm1.New(n_elements,n_gauss);
  y_nm1.New(n_elements,n_gauss);
  z_nm1.New(n_elements,n_gauss);

  for(int ielem=0; ielem<n_elements; ++ielem)
    {
      for(int ig=0; ig<n_gauss; ++ig)
	{
	  jac_tab= Jacobi_mat(  ielem, ig );    // Call the jacobi_mat function
	  det_jac1(ielem,ig)= wg[ig]*jac_tab[0];// multiplied by wg[ig] to save time in nonsin function
	  
	  x_nm1(ielem,ig) = jac_tab[1];         // he normal compoments 
	  y_nm1(ielem,ig) = jac_tab[2]; 
	  z_nm1(ielem,ig) = jac_tab[3];
	}
    }
  
  //############### The main Loop ###########################################################;
     
  for(int nodep=0; nodep<n_nodes; nodep++)
      {
	xp = x[nodep]; yp = y[nodep]; zp = z[nodep];	  // define the coordinates of the load point "nodep"
	
	for(int ielemq=0; ielemq<n_elements; ielemq++)    // now take each element in turn as the field element "ielemq"
	  {

	      for(int ic=0; ic<dimension; ++ic)
	      {
		norder[ic]=node(ielemq,ic);
	      }
	    for(int ic=0; ic<dimension; ++ic)             // and take each node of the ielemq as nodeq
	      {	
		nodeq=  node(ielemq,ic);
		icd=ic+dimension;
		//                                            use this to treat all nodes in the element as singular 	
		if(nodep!=norder[0]&&nodep!=norder[1]&&nodep!=norder[2]&& nodep!=norder[3]&&nodep!=norder[4]&&nodep!=norder[5])
		
		//                                            and this to treat only the same node as singular 
		//if(nodep!=nodeq)                          
		  {				  
		    //non-singular integrations kernel2() and kernel1() calculation;
		    a_aux=nonsin_mat( ic, ielemq, xp, yp, zp, 2);
		    a(nodep,nodeq) += a_aux[ic];
		    b(nodep,nodeq) += a_aux[icd]; 
		  }
		else
		  {
		    //singular integrations
		    //a_aux=singu4(0, ic, nodep, ielemq, xp, yp, zp); //use this for stricter treatment of singularity
		    a_aux=singu2(0, ic, nodep, ielemq, xp, yp, zp, 2);
		    a(nodep,nodeq) += a_aux[ic]; //This was missing from Jans Code 	
		    b(nodep,nodeq) += a_aux[icd];
		  }	
	      }//nodeq
	  }//ielemq
	// so far, all coefficients of A and B matrices, have been calculated;
	// except for the diagonal terms of the A-matrix;
	// The diagonal terms of matrix A can be determined as the 0.5 ;
	// a(nodep,nodep) = complex(0.50,0.);	// Done in Matlab code !!!;
      }//nodep
 
}// end of bem_ms

//#######################################################################################################
//'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

void bemdiff::single_region_ma(RDenseMatrix &S_nodes, IDenseMatrix &S_elements, 
			       double mua, double mus, double freq, CDenseMatrix &a, CDenseMatrix &b)
{

  // #################  Optical Properties ##############################################;
  
  //Initialise optical parameters vectors;
  D.New(1);                                      // Diffusion Coefficient
  k.New(1);                                      // Wave Number
  //used for the dk/dmua
  ma.New(1);
  ms.New(1);   
  iwc.New(1);                              
  const double nu=1.;                            // Refractive Index (to be included)
  const double c=0.3*1.e+12/nu;                  // Speed of light in the medium in mm/ps
  double omega=2.*pi*freq;                       // Angular speed
  D[0] = 1./(3.*(mua+mus));                      // Diffusion coefficient
  k[0] = sqrt(std::complex<double>(mua*c, -omega)/(c*D[0]));  // Wave number
  ma[0]=mua;
  ms[0]=mus;
  iwc[0]=std::complex<double>(0,omega/c);
  
  // ############## The Geometrical  definitions ########################################;  
 
  
  n_elements = S_elements.nRows();               // Number of Elements in the  region
  n_nodes    = S_nodes.nRows();                  // Number of Nodes in the  region
  node = S_elements;                             // node - integer matrix stored as an n_elements by dimension
   
  x=S_nodes.Col(0);                              // The Coordinates of the nodes 
  y=S_nodes.Col(1);
  z=S_nodes.Col(2); 

  // ################ The Jan definitions ###############################################;  
  dimension  = 6;                                // Dimension of the Triangles 
  dim2       = 12;                               // number of derivatives of Shape functions 
  i_kern     = 4;                                // length of the kernel function return
  int i_col=0;
  CVector a_aux(dim2);
  int icd, nodeq,  norder[6];

 
  double xp, yp, zp;
 
  a.New(n_nodes,n_nodes);                        // a matrix for the return    
  b.New(n_nodes,n_nodes);                        // b matrix for the return
 
  // ############### Precalculate the Jacobian for nonsin #################################;
 
  RVector jac_tab(dimension+1);
 
  det_jac1.New(n_elements,n_gauss);             // Determinant 
  x_nm1.New(n_elements,n_gauss);
  y_nm1.New(n_elements,n_gauss);
  z_nm1.New(n_elements,n_gauss);

  for(int ielem=0; ielem<n_elements; ++ielem)
    {
      for(int ig=0; ig<n_gauss; ++ig)
	{
	  jac_tab= Jacobi_mat(  ielem, ig );    // Call the jacobi_mat function
	  det_jac1(ielem,ig)= wg[ig]*jac_tab[0];// multiplied by wg[ig] to save time in nonsin function
	  
	  x_nm1(ielem,ig) = jac_tab[1];         // he normal compoments 
	  y_nm1(ielem,ig) = jac_tab[2]; 
	  z_nm1(ielem,ig) = jac_tab[3];
	}
    }
  
  //############### The main Loop ###########################################################;
     
  for(int nodep=0; nodep<n_nodes; nodep++)
      {
	xp = x[nodep]; yp = y[nodep]; zp = z[nodep];	  // define the coordinates of the load point "nodep"
	
	for(int ielemq=0; ielemq<n_elements; ielemq++)    // now take each element in turn as the field element "ielemq"
	  {
	      for(int ic=0; ic<dimension; ++ic)
	      {
		norder[ic]=node(ielemq,ic);
	      }
	    for(int ic=0; ic<dimension; ++ic)             // and take each node of the ielemq as nodeq
	      {	
		nodeq=  node(ielemq,ic);
		icd=ic+dimension;
		//                                            use this to treat all nodes in the element as singular 	
		if(nodep!=norder[0]&&nodep!=norder[1]&&nodep!=norder[2]&& nodep!=norder[3]&&nodep!=norder[4]&&nodep!=norder[5])
		
		//                                            and this to treat only the same node as singular 
		//if(nodep!=nodeq)                          
		  {				  
		    //non-singular integrations kernel2() and kernel1() calculation;
		    a_aux=nonsin_mat( ic, ielemq, xp, yp, zp, 1);
		    a(nodep,nodeq) += a_aux[ic];
		    b(nodep,nodeq) += a_aux[icd]; 
		  }
		else
		  {
		    //singular integrations
		    //a_aux=singu4(0, ic, nodep, ielemq, xp, yp, zp); //use this for stricter treatment of singularity
		    a_aux=singu2(0, ic, nodep, ielemq, xp, yp, zp, 1);
		    a(nodep,nodeq) += a_aux[ic]; //This was missing from Jans Code 	
		    b(nodep,nodeq) += a_aux[icd];
		  }	
	      }//nodeq
	  }//ielemq
	// so far, all coefficients of A and B matrices, have been calculated;
	// except for the diagonal terms of the A-matrix;
	// The diagonal terms of matrix A can be determined as the 0.5 ;
	// a(nodep,nodep) = complex(0.50,0.);	// Done in Matlab code !!!;
      }//nodep
 
}// end of Bem_ma




void bemdiff::single_region_3D(RDenseMatrix &S_nodes, IDenseMatrix &S_elements, 
			       double mua, double mus, double freq, CDenseMatrix &a, CDenseMatrix &b)
{



  //cerr<<" mua "<<mua<<"mus "<<mus<<"freq"<<freq<<endl;
  // #################  Optical Properties ##############################################;
  
  //Initialise optical parameters vectors;
  D.New(1);                                      // Diffusion Coefficient
  k.New(1);                                      // Wave Number
  //used for the dk/dmua
  ma.New(1);
  ms.New(1);   
  iwc.New(1);                              
  const double nu=1.;                            // Refractive Index (to be included)
  const double c=0.3*1.e+12/nu;                  // Speed of light in the medium in mm/ps
  double omega=2.*pi*freq;                       // Angular speed
  D[0] = 1./(3.*(mua+mus));                      // Diffusion coefficient
  k[0] = sqrt(std::complex<double>(mua*c, -omega)/(c*D[0]));  // Wave number
  ma[0]=mua;
  ms[0]=mus;
  iwc[0]=std::complex<double>(0,omega/c);
  
  // ############## The Geometrical  definitions ########################################;  
 
  
  n_elements = S_elements.nRows();               // Number of Elements in the  region
  n_nodes    = S_nodes.nRows();                  // Number of Nodes in the  region
  node = S_elements;                             // node - integer matrix stored as an n_elements by dimension
   
  x=S_nodes.Col(0);                              // The Coordinates of the nodes 
  y=S_nodes.Col(1);
  z=S_nodes.Col(2); 

  // ################ The Jan definitions ###############################################;  
  dimension  = 6;                                // Dimension of the Triangles 
  dim2       = 12;                               // number of derivatives of Shape functions 
  i_kern     = 4;                                // length of the kernel function return
  int i_col=0;
  CVector a_aux(dim2);
  int icd, nodeq, norder[6];
  

  double xp, yp, zp;
 
  a.New(n_nodes,n_nodes);                        // a matrix for the return    
  b.New(n_nodes,n_nodes);                        // b matrix for the return
 
  // ############### Precalculate the Jacobian for nonsin #################################;
 
  RVector jac_tab(dimension+1);
 
  det_jac1.New(n_elements,n_gauss);             // Determinant 
  x_nm1.New(n_elements,n_gauss);
  y_nm1.New(n_elements,n_gauss);
  z_nm1.New(n_elements,n_gauss);

  for(int ielem=0; ielem<n_elements; ++ielem)
    {
      for(int ig=0; ig<n_gauss; ++ig)
	{
	  jac_tab= Jacobi_mat(  ielem, ig );    // Call the jacobi_mat function
	  det_jac1(ielem,ig)= wg[ig]*jac_tab[0];// multiplied by wg[ig] to save time in nonsin function
	  
	  x_nm1(ielem,ig) = jac_tab[1];         // he normal compoments 
	  y_nm1(ielem,ig) = jac_tab[2]; 
	  z_nm1(ielem,ig) = jac_tab[3];
	}
    }
  
  //############### The main Loop ###########################################################;
     
  for(int nodep=0; nodep<n_nodes; nodep++)
      {
	xp = x[nodep]; yp = y[nodep]; zp = z[nodep];	  // define the coordinates of the load point "nodep"
	
	for(int ielemq=0; ielemq<n_elements; ielemq++)    // now take each element in turn as the field element "ielemq"
	  {

	    for(int ic=0; ic<dimension; ++ic)
	      {
		norder[ic]=node(ielemq,ic);
	      }
	    for(int ic=0; ic<dimension; ++ic)             // and take each node of the ielemq as nodeq
	      {	
		nodeq=  node(ielemq,ic);
		icd=ic+dimension;
		//                                            use this to treat all nodes in the element as singular 	
		if(nodep!=norder[0]&&nodep!=norder[1]&&nodep!=norder[2]&& nodep!=norder[3]&&nodep!=norder[4]&&nodep!=norder[5])
		
		//                                            and this to treat only the same node as singular 
		//if(nodep!=nodeq)                          
		  {				  
		    //non-singular integrations kernel2() and kernel1() calculation;
		    a_aux=nonsin_mat( ic, ielemq, xp, yp, zp, 0);
                    //a_aux=nonsin( ic, ielemq, xp, yp, zp);
		    a(nodep,nodeq) += a_aux[ic];
		    b(nodep,nodeq) += a_aux[icd]; 
		  }
		else
		  {
		    //singular integrations
		    //a_aux=singu4(0, ic, nodep, ielemq, xp, yp, zp); //use this for stricter treatment of singularity
		    a_aux=singu2(0, ic, nodep, ielemq, xp, yp, zp, 0);
		    //   lomef<< a_aux[ic]<<" "<<ic<<" "<< nodep<<" "<< ielemq<<" "<< xp<<" "<< yp<<" "<< zp<<endl;
		    a(nodep,nodeq) += a_aux[ic]; //This was missing from Jans Code 	
		    b(nodep,nodeq) += a_aux[icd];
		  }	
	      }//nodeq
	  }//ielemq
	// so far, all coefficients of A and B matrices, have been calculated;
	// except for the diagonal terms of the A-matrix;
	// The diagonal terms of matrix A can be determined as the 0.5 ;
	// a(nodep,nodep) = complex(0.50,0.);	// Done in Matlab code !!!;
      }//nodep
  // lomef<<"region"<<endl;
}// end of BemSingleRegion

//################ REFRACTIVE INDEX ###############################################


void bemdiff::single_region_refidx_3D(RDenseMatrix &S_nodes, IDenseMatrix &S_elements, 
			       double mua, double mus, double freq, double nu, CDenseMatrix &a, CDenseMatrix &b)
{



  //cerr<<" mua "<<mua<<"mus "<<mus<<"freq"<<freq<<endl;
  // #################  Optical Properties ##############################################;
  
  //Initialise optical parameters vectors;
  D.New(1);                                      // Diffusion Coefficient
  k.New(1);                                      // Wave Number
  //used for the dk/dmua
  ma.New(1);
  ms.New(1);   
  iwc.New(1);                              
  // const double nu=1.;                            // Refractive Index (to be included)
  const double c=0.3*1.e+12/nu;                  // Speed of light in the medium in mm/ps
  double omega=2.*pi*freq;                       // Angular speed
  D[0] = 1./(3.*(mua+mus));                      // Diffusion coefficient
  k[0] = sqrt(std::complex<double>(mua*c, -omega)/(c*D[0]));  // Wave number
  ma[0]=mua;
  ms[0]=mus;
  iwc[0]=std::complex<double>(0,omega/c);
  
  // ############## The Geometrical  definitions ########################################;  
 
  
  n_elements = S_elements.nRows();               // Number of Elements in the  region
  n_nodes    = S_nodes.nRows();                  // Number of Nodes in the  region
  node = S_elements;                             // node - integer matrix stored as an n_elements by dimension
   
  x=S_nodes.Col(0);                              // The Coordinates of the nodes 
  y=S_nodes.Col(1);
  z=S_nodes.Col(2); 

  // ################ The Jan definitions ###############################################;  
  dimension  = 6;                                // Dimension of the Triangles 
  dim2       = 12;                               // number of derivatives of Shape functions 
  i_kern     = 4;                                // length of the kernel function return
  int i_col=0;
  CVector a_aux(dim2);
  int icd, nodeq, norder[6];

  double xp, yp, zp;
 
  a.New(n_nodes,n_nodes);                        // a matrix for the return    
  b.New(n_nodes,n_nodes);                        // b matrix for the return
 
  // ############### Precalculate the Jacobian for nonsin #################################;
 
  RVector jac_tab(dimension+1);
 
  det_jac1.New(n_elements,n_gauss);             // Determinant 
  x_nm1.New(n_elements,n_gauss);
  y_nm1.New(n_elements,n_gauss);
  z_nm1.New(n_elements,n_gauss);

  for(int ielem=0; ielem<n_elements; ++ielem)
    {
      for(int ig=0; ig<n_gauss; ++ig)
	{
	  jac_tab= Jacobi_mat(  ielem, ig );    // Call the jacobi_mat function
	  det_jac1(ielem,ig)= wg[ig]*jac_tab[0];// multiplied by wg[ig] to save time in nonsin function
	  
	  x_nm1(ielem,ig) = jac_tab[1];         // he normal compoments 
	  y_nm1(ielem,ig) = jac_tab[2]; 
	  z_nm1(ielem,ig) = jac_tab[3];
	}
    }
  
  //############### The main Loop ###########################################################;
     
  for(int nodep=0; nodep<n_nodes; nodep++)
      {
	xp = x[nodep]; yp = y[nodep]; zp = z[nodep];	  // define the coordinates of the load point "nodep"
	
	for(int ielemq=0; ielemq<n_elements; ielemq++)    // now take each element in turn as the field element "ielemq"
	  {

	    for(int ic=0; ic<dimension; ++ic)
	      {
		norder[ic]=node(ielemq,ic);
	      }
	    for(int ic=0; ic<dimension; ++ic)             // and take each node of the ielemq as nodeq
	      {	
		nodeq=  node(ielemq,ic);
		icd=ic+dimension;
		//                                            use this to treat all nodes in the element as singular 	
		if(nodep!=norder[0]&&nodep!=norder[1]&&nodep!=norder[2]&& nodep!=norder[3]&&nodep!=norder[4]&&nodep!=norder[5])
		
		//                                            and this to treat only the same node as singular 
		//if(nodep!=nodeq)                          
		  {				  
		    //non-singular integrations kernel2() and kernel1() calculation;
		    a_aux=nonsin_mat( ic, ielemq, xp, yp, zp, 0);
                    //a_aux=nonsin( ic, ielemq, xp, yp, zp);
		    a(nodep,nodeq) += a_aux[ic];
		    b(nodep,nodeq) += a_aux[icd]; 
		  }
		else
		  {
		    //singular integrations
		    //a_aux=singu4(0, ic, nodep, ielemq, xp, yp, zp); //use this for stricter treatment of singularity
		    a_aux=singu2(0, ic, nodep, ielemq, xp, yp, zp, 0);
		    //   lomef<< a_aux[ic]<<" "<<ic<<" "<< nodep<<" "<< ielemq<<" "<< xp<<" "<< yp<<" "<< zp<<endl;
		    a(nodep,nodeq) += a_aux[ic]; //This was missing from Jans Code 	
		    b(nodep,nodeq) += a_aux[icd];
		  }	
	      }//nodeq
	  }//ielemq
	// so far, all coefficients of A and B matrices, have been calculated;
	// except for the diagonal terms of the A-matrix;
	// The diagonal terms of matrix A can be determined as the 0.5 ;
	// a(nodep,nodep) = complex(0.50,0.);	// Done in Matlab code !!!;
      }//nodep
  // lomef<<"region"<<endl;
}// end of BemSingleRegion






//########### OTHER FUNCTIONS ######################################################


inline RVector bemdiff::shape0f(double ksi1, double ksi2)
{
  static RVector shapf(dimension);

// zero order shape functions in respect to ksi1 ksi2
  shapf[0]=1.-ksi1-ksi2; shapf[1]=ksi1; shapf[2]=ksi2;

  return shapf;
}


inline RVector bemdiff::shape0d() 
{
  static RVector shapd(dim2);
 
  //calculate the derivatives of the zeroth order shape functions
  shapd[0]=-1.; shapd[1]= 1.; shapd[2]= 0.;

  shapd[3]=-1.; shapd[4]= 0.; shapd[5]= 1.;

  return shapd;
}

inline RVector bemdiff::shapef(double ksi1, double ksi2) 
{

  //==================================================================
  //  purpose: Calculate the  shape functions in respect to ksi1 ksi2
  //  
  //==================================================================
  double ksi3;
  static RVector shapf(dimension);
  ksi3=1.-ksi1-ksi2;
  // calculate the quadratic shape functions
  shapf[0]=-ksi3*(1.-2.*ksi3);  shapf[1]=4.*ksi1*ksi3;
  shapf[2]=-ksi1*(1.-2.*ksi1);  shapf[3]=4.*ksi1*ksi2;
  shapf[4]=-ksi2*(1.-2.*ksi2);  shapf[5]=4.*ksi2*ksi3;
  
  return shapf;
}


//####################################################################################

inline RVector bemdiff::shaped(double ksi1, double ksi2) 
{
 
//==================================================================
//  purpose: Calculate derivative of shape functions 
//           in respect to ksi1 ksi2 
//==================================================================
  double ksi3, ksi1_4, ksi2_4, ksi3_4;
  static RVector shapd(dim2);

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


//####################################################################################

RVector bemdiff::Jacobi(int i_region, int ielem, double ksi1, double ksi2)
{
 
  //==================================================================
  //  purpose: Calculate the Jacobian matrix entries
  //  
  //==================================================================

  int icd, iic;
  double dx_d_ksi1, dy_d_ksi1, dz_d_ksi1, dx_d_ksi2, dy_d_ksi2, dz_d_ksi2, det_jacob, nx, ny, nz;
  static RVector jac_tab(dimension+1), shapd(dim2);

  // local coordinate system
  dx_d_ksi1 = 0.; dy_d_ksi1 = 0.; dz_d_ksi1 = 0.; 
  dx_d_ksi2 = 0.; dy_d_ksi2 = 0.; dz_d_ksi2 = 0.;

  shapd=shaped(ksi1,ksi2);
  
  for(int ic=0; ic<dimension;++ic)
    {
      iic=node(ielem,ic); icd=ic+dimension;
      dx_d_ksi1 += shapd[ic]*x[iic];
      dy_d_ksi1 += shapd[ic]*y[iic];
      dz_d_ksi1 += shapd[ic]*z[iic];
      
      dx_d_ksi2 += shapd[icd]*x[iic];
      dy_d_ksi2 += shapd[icd]*y[iic];
      dz_d_ksi2 += shapd[icd]*z[iic];
    }   
  
  // calculate the unit normal components
  nx = dy_d_ksi1*dz_d_ksi2-dz_d_ksi1*dy_d_ksi2;
  ny = dz_d_ksi1*dx_d_ksi2-dx_d_ksi1*dz_d_ksi2;
  nz = dx_d_ksi1*dy_d_ksi2-dy_d_ksi1*dx_d_ksi2;
   
  // calculate the determinant of Jacobian as a function of point "ksi1, ksi2"
  det_jacob = 1./sqrt(nx*nx + ny*ny + nz*nz);
  
  // jac_tab table collect the Jacobian and directional components 
  // of the outward normal vector to the boundary (numbering system-anticlocwise) 
  jac_tab[0] = 1./det_jacob;
  jac_tab[1] = -nx*det_jacob; 
  jac_tab[2] = -ny*det_jacob; 
  jac_tab[3] = -nz*det_jacob;

  return jac_tab;
}


//####################################################################################

RVector bemdiff::Jacobi_mat( int ielem,int ig )
{
  //==================================================================
  //  purpose: Calculate the Jacobi matrix entries
  //           with precalculated shaped for 7 gaussian points
  //==================================================================
  int icd, iic;
  double dx_d_ksi1, dy_d_ksi1, dz_d_ksi1, dx_d_ksi2, dy_d_ksi2, dz_d_ksi2, det_jacob, nx, ny, nz;
  static RVector jac_tab(dimension+1);
 
  // local coordinate system
  dx_d_ksi1 = 0.; dy_d_ksi1 = 0.; dz_d_ksi1 = 0.; 
  dx_d_ksi2 = 0.; dy_d_ksi2 = 0.; dz_d_ksi2 = 0.;

  double  shapd[12]={shapd0[ig],shapd1[ig],shapd2[ig],shapd3[ig],shapd4[ig],shapd5[ig]
		     ,shapd6[ig],shapd7[ig],shapd8[ig],shapd9[ig],shapd10[ig],shapd11[ig]};

  for(int ic=0; ic<dimension;++ic)
    {
      iic=node(ielem,ic); icd=ic+dimension;
      dx_d_ksi1 += shapd[ic]*x[iic];
      dy_d_ksi1 += shapd[ic]*y[iic];
      dz_d_ksi1 += shapd[ic]*z[iic];
      
      dx_d_ksi2 += shapd[icd]*x[iic];
      dy_d_ksi2 += shapd[icd]*y[iic];
      dz_d_ksi2 += shapd[icd]*z[iic];
    }   
  
  // calculate the unit normal components
  nx = dy_d_ksi1*dz_d_ksi2-dz_d_ksi1*dy_d_ksi2;
  ny = dz_d_ksi1*dx_d_ksi2-dx_d_ksi1*dz_d_ksi2;
  nz = dx_d_ksi1*dy_d_ksi2-dy_d_ksi1*dx_d_ksi2;
    
  // calculate the determinant of Jacobian as a function of point "ksi1, ksi2"
  det_jacob = 1./sqrt(nx*nx + ny*ny + nz*nz);
  
  //jac_tab table collect the Jacobian and directional components 
  //of the outward normal vector to the boundary (numbering system-anticlocwise) 
  jac_tab[0] = 1./det_jacob;
  jac_tab[1] = -nx*det_jacob; 
  jac_tab[2] = -ny*det_jacob; 
  jac_tab[3] = -nz*det_jacob;

  return jac_tab;
}


//####################################################################################

CVector bemdiff::kernel(int i_region, int ielem, double ksi1, double ksi2, double xp, double yp, double zp)
{
  
  //==================================================================
  //  purpose: Calculate the kernals (Green's Function)
  //  Diffusion equation
  //==================================================================
  int iic;  
  double xq, yq, zq, rpq1, rpq2, xpxq, ypyq, zpzq;
  std::complex<double> krpq1, k_aux, e_aux;
  static RVector shapf(dimension);
  static CVector kern_tab(i_kern);

  // local coordinate system
  xq = 0.;    yq = 0.;    zq = 0.;
  shapf=shapef(ksi1, ksi2);
  
  for(int ic=0; ic<dimension;++ic){
    iic=node(ielem,ic);
    xq += shapf[ic]*x[iic];
    yq += shapf[ic]*y[iic];
    zq += shapf[ic]*z[iic];
  }

  xpxq=xq-xp; ypyq=yq-yp; zpzq=zq-zp;
  rpq2=1./(xpxq*xpxq+ypyq*ypyq+zpzq*zpzq);
  rpq1=sqrt(rpq2);
  
  krpq1= -k[i_region]/rpq1;
  e_aux=exp(krpq1);
  k_aux=-(k[i_region]+rpq1)*e_aux*pi4*rpq2;

  // kernel1 ->  kern_tab[0]; kern_tab[1]; kern_tab[2]; 
  kern_tab[0] = k_aux*xpxq;
  kern_tab[1] = k_aux*ypyq;
  kern_tab[2] = k_aux*zpzq;
  
  // the kernel2 -> kern_tab[3] 
  kern_tab[3] = e_aux*rpq1*pi4;
  
  return kern_tab;
}


//####################################################################################

CVector bemdiff::kernel_mat( int ielem, int ig , double xp, double yp, double zp)
{
  
  //==================================================================
  //  purpose: Calculate the kernals (Green's Function)
  //           Diffusion equation 
  //           only for the n_gauss : 7 gaussian points precalculated
  //==================================================================
  int iic;  
  double xq, yq, zq, rpq1, rpq2, xpxq, ypyq, zpzq;
  std::complex<double> krpq1, k_aux, e_aux;
  static CVector kern_tab(i_kern);
    
  // local coordinate system
  xq = 0.;    yq = 0.;    zq = 0.;

  double  shapf[6]={shapf0[ig],shapf1[ig],shapf2[ig],shapf3[ig],shapf4[ig],shapf5[ig]};
  for(int ic=0; ic<dimension;++ic){
    
    iic=node(ielem,ic);
    xq += shapf[ic]*x[iic];
    yq += shapf[ic]*y[iic];
    zq += shapf[ic]*z[iic];
  }
  xpxq=xq-xp; ypyq=yq-yp; zpzq=zq-zp;
  rpq2=1./(xpxq*xpxq+ypyq*ypyq+zpzq*zpzq);
  rpq1=sqrt(rpq2);
  
  krpq1= -k[0]/rpq1;
  e_aux=exp(krpq1);
  k_aux=-(k[0]+rpq1)*e_aux*pi4*rpq2;
  
  // kernel1 ->  kern_tab[0]; kern_tab[1]; kern_tab[2]; 
  kern_tab[0] = k_aux*xpxq;
  kern_tab[1] = k_aux*ypyq;
  kern_tab[2] = k_aux*zpzq;
  
  // the kernel2 -> kern_tab[3] 
  kern_tab[3] = e_aux*rpq1*pi4;
  
  return kern_tab;
}

//####################################################################################

CVector bemdiff::Dkernel_dms_mat( int ielem, int ig , double xp, double yp, double zp)
{
  
  //==================================================================
  //  purpose: Calculate the kernals (Green's Function)
  //           Diffusion equation 
  //           only for the n_gauss : 7 gaussian points precalculated
  //==================================================================
  int iic;  
  double xq, yq, zq, xpxq, ypyq, zpzq, rpq;
  std::complex<double> krpq, k_aux, e_aux;
  static CVector kern_tab(i_kern);
  const double c=0.3*1.e+12/1;//nu here
  // local coordinate system
  xq = 0.;    yq = 0.;    zq = 0.;

  double  shapf[6]={shapf0[ig],shapf1[ig],shapf2[ig],shapf3[ig],shapf4[ig],shapf5[ig]};
  for(int ic=0; ic<dimension;++ic){
    
    iic=node(ielem,ic);
    xq += shapf[ic]*x[iic];
    yq += shapf[ic]*y[iic];
    zq += shapf[ic]*z[iic];
  }
  xpxq = xq-xp; ypyq = yq-yp; zpzq = zq-zp;
  rpq = sqrt(xpxq*xpxq+ypyq*ypyq+zpzq*zpzq); 
  
  krpq= -k[0]*rpq;
  e_aux=exp(krpq);


  k_aux= ( (-iwc[0]+ma[0])* e_aux *3. )*(1/(8. *c* rpq * pi));
  
  // kernel1 ->  kern_tab[0]; kern_tab[1]; kern_tab[2]; 
  kern_tab[0] = k_aux*xpxq;
  kern_tab[1] = k_aux*ypyq;
  kern_tab[2] = k_aux*zpzq;
  
  // the kernel2 -> kern_tab[3] 
  kern_tab[3] =  -( k[0] * e_aux)/(8*pi*(ma[0]+ms[0]));
  
  return kern_tab;
}

//####################################################################################

CVector bemdiff::Dkernel_dms(int i_region, int ielem, double ksi1, double ksi2, double xp, double yp, double zp)
{
  
  //==================================================================
  //  purpose: Calculate the kernals (Green's Function)
  //           Diffusion equation 
  //           only for the n_gauss : 7 gaussian points precalculated
  //==================================================================
  int iic;  
  double xq, yq, zq, xpxq, ypyq, zpzq, rpq;
  std::complex<double> krpq, k_aux, e_aux;
  static CVector kern_tab(i_kern);
  static RVector shapf(dimension);

  const double c=0.3*1.e+12/1;//nu here
  // local coordinate system
  xq = 0.;    yq = 0.;    zq = 0.;
  shapf=shapef(ksi1, ksi2);
  // double  shapf[6]={shapf0[ig],shapf1[ig],shapf2[ig],shapf3[ig],shapf4[ig],shapf5[ig]};
  for(int ic=0; ic<dimension;++ic){
    
    iic=node(ielem,ic);
    xq += shapf[ic]*x[iic];
    yq += shapf[ic]*y[iic];
    zq += shapf[ic]*z[iic];
  }
  xpxq = xq-xp; ypyq = yq-yp; zpzq = zq-zp;
  rpq = sqrt(xpxq*xpxq+ypyq*ypyq+zpzq*zpzq); 
  
  krpq= -k[0]*rpq;
  e_aux=exp(krpq);


  k_aux= ( (-iwc[0]+ma[0])* e_aux *3. )*(1/(8. *c* rpq * pi));
  
  // kernel1 ->  kern_tab[0]; kern_tab[1]; kern_tab[2]; 
  kern_tab[0] = k_aux*xpxq;
  kern_tab[1] = k_aux*ypyq;
  kern_tab[2] = k_aux*zpzq;
  
  // the kernel2 -> kern_tab[3] 
  kern_tab[3] =  -( k[0] * e_aux)/(8*pi*(ma[0]+ms[0]));
  
  return kern_tab;
}

//####################################################################################



CVector bemdiff::Dkernel_dma_mat( int ielem, int ig , double xp, double yp, double zp)
{
  
  //==================================================================
  //  purpose: Calculate the kernals (Green's Function)
  //           Diffusion equation 
  //           only for the n_gauss : 7 gaussian points precalculated
  //==================================================================
  int iic;  
  double xq, yq, zq, xpxq, ypyq, zpzq, rpq;
  std::complex<double> krpq, k_aux, e_aux, MaMsIwc;
  static CVector kern_tab(i_kern);
  const double c=0.3*1.e+12/1;//nu here
  // local coordinate system
  xq = 0.;    yq = 0.;    zq = 0.;

  double  shapf[6]={shapf0[ig],shapf1[ig],shapf2[ig],shapf3[ig],shapf4[ig],shapf5[ig]};
  for(int ic=0; ic<dimension;++ic){
    
    iic=node(ielem,ic);
    xq += shapf[ic]*x[iic];
    yq += shapf[ic]*y[iic];
    zq += shapf[ic]*z[iic];
  }
  xpxq = xq-xp; ypyq = yq-yp; zpzq = zq-zp;
  rpq = sqrt(xpxq*xpxq+ypyq*ypyq+zpzq*zpzq);
  MaMsIwc=(-iwc[0]+2*ma[0]+ms[0]);
  
  krpq= -k[0]*rpq;
  e_aux=exp(krpq);

  k_aux= ( MaMsIwc* e_aux * 3.0)/(8*c* rpq * pi);
  
  // kernel1 ->  kern_tab[0]; kern_tab[1]; kern_tab[2]; 
  kern_tab[0] = k_aux*xpxq;
  kern_tab[1] = k_aux*ypyq;
  kern_tab[2] = k_aux*zpzq;
  
  // the kernel2 -> kern_tab[3] 
  kern_tab[3] =  -( MaMsIwc * e_aux * 3.0)/(k[0]*8.0*pi);
  
  return kern_tab;
}

//####################################################################################


//####################################################################################



CVector bemdiff::Dkernel_dma(int i_region, int ielem, double ksi1, double ksi2, double xp, double yp, double zp)
{
  
  //==================================================================
  //  purpose: Calculate the kernals (Green's Function)
  //           Diffusion equation 
  //           only for the n_gauss : 7 gaussian points precalculated
  //==================================================================
  int iic;  
  double xq, yq, zq, xpxq, ypyq, zpzq, rpq;
  std::complex<double> krpq, k_aux, e_aux, MaMsIwc;
  static CVector kern_tab(i_kern);
  static RVector shapf(dimension);

  const double c=0.3*1.e+12/1;//nu here
  // local coordinate system
  xq = 0.;    yq = 0.;    zq = 0.;
  shapf=shapef(ksi1, ksi2);
  // double  shapf[6]={shapf0[ig],shapf1[ig],shapf2[ig],shapf3[ig],shapf4[ig],shapf5[ig]};
  for(int ic=0; ic<dimension;++ic){
    
    iic=node(ielem,ic);
    xq += shapf[ic]*x[iic];
    yq += shapf[ic]*y[iic];
    zq += shapf[ic]*z[iic];
  }
  xpxq = xq-xp; ypyq = yq-yp; zpzq = zq-zp;
  rpq = sqrt(xpxq*xpxq+ypyq*ypyq+zpzq*zpzq);
  MaMsIwc=(-iwc[0]+2*ma[0]+ms[0]);
  
  krpq= -k[0]*rpq;
  e_aux=exp(krpq);

  k_aux= ( MaMsIwc* e_aux * 3.0)/(8*c* rpq * pi);
  
  // kernel1 ->  kern_tab[0]; kern_tab[1]; kern_tab[2]; 
  kern_tab[0] = k_aux*xpxq;
  kern_tab[1] = k_aux*ypyq;
  kern_tab[2] = k_aux*zpzq;
  
  // the kernel2 -> kern_tab[3] 
  kern_tab[3] =  -( MaMsIwc * e_aux * 3.0)/(k[0]*8.0*pi);
  
  return kern_tab;
}

//####################################################################################




//####################################################################################





//####################################################################################


//####################################################################################






//####################################################################################

CVector bemdiff::singu2(int i_region, int ic, int nodep, int ielemq, 
			  double xp, double yp, double zp, int g0_ma1_ms2)
{  
  int icd, *norder = new int[dimension]; 
  double ksi1, ksi2, det_J_r;
  std::complex<double> kern1, kern2;

  static RVector jac_tab(dimension+1), shapf(dimension);
  static CVector kern_tab(i_kern), a_aux(dim2);
//==================================================================
//  purpose: Singular integrals in triangular element
//  element is divided by 2 sub-triangles
//==================================================================

//matrix a_aux initialization
  for(int i=0;i<dim2;++i){
    a_aux[i]=std::complex<double>(0.,0.);
  }
  icd= ic + dimension;

  //  node 0------------------------------------------------
  if(nodep==node(ielemq,0)){ 
    for(int ig=0; ig<n_gauss_o; ig++){
     
      det_J_r=(1.+xg_k_1[ig])/8.;
          
      for(int jg=0; jg<n_gauss_o; jg++){

	// Integral t1

	ksi1=(1.+xg_k_1[ig])*(1.-xg_k_2[jg])/4.; 
	ksi2=(1.+xg_k_1[ig])*(1.+xg_k_2[jg])/4.; 

	shapf=shapef(ksi1, ksi2);    

	jac_tab=Jacobi(i_region,ielemq, ksi1, ksi2);//, node);

	//	lomef<< jac_tab<<" : "<<ielemq<<" "<< ksi1<<" "<< ksi2<<endl;

	if (g0_ma1_ms2==0)
	{
	kern_tab=kernel(i_region, ielemq, ksi1, ksi2, xp, yp, zp);//, node);
	}
	if (g0_ma1_ms2==1)
	{
	kern_tab=Dkernel_dma(i_region, ielemq, ksi1, ksi2, xp, yp, zp);//, node);
	}
	if (g0_ma1_ms2==2)
	{
	kern_tab=Dkernel_dms(i_region, ielemq, ksi1, ksi2, xp, yp, zp);//, node);
	}

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
  if(nodep==node(ielemq,1)){ 
    for(int ig=0; ig<n_gauss_o; ig++){

      det_J_r=(1.+xg_k_1[ig])/16.;
     

      for(int jg=0; jg<n_gauss_o; jg++){

	// Integral t1

	ksi1=(1.-xg_k_1[ig])/4.;
	ksi2=(1.+xg_k_1[ig])*(1.-xg_k_2[jg])/4.;

	shapf=shapef(ksi1, ksi2);    

	jac_tab=Jacobi(i_region, ielemq, ksi1, ksi2);//, node);

	// differenct kernels for greens / ma and ms 
	if (g0_ma1_ms2==0)
	{
	kern_tab=kernel(i_region, ielemq, ksi1, ksi2, xp, yp, zp);//, node);
	}
	if (g0_ma1_ms2==1)
	{
	kern_tab=Dkernel_dma(i_region, ielemq, ksi1, ksi2, xp, yp, zp);//, node);
	}
	if (g0_ma1_ms2==2)
	{
	kern_tab=Dkernel_dms(i_region, ielemq, ksi1, ksi2, xp, yp, zp);//, node);
	}

	kern1 = kern_tab[0]*jac_tab[1]+kern_tab[1]*jac_tab[2]+
	  kern_tab[2]*jac_tab[3];//Simon!

	a_aux[ic] += kern1*wg_k[ig]*wg_k[jg]*shapf[ic]*
	  jac_tab[0]*det_J_r;//Simon!
	  
	a_aux[icd]+= kern_tab[3]*wg_k[ig]*wg_k[jg]*shapf[ic]*
	  jac_tab[0]*det_J_r;
	//a_aux[icd]+= wg_k[ig]*wg_k[jg]*jac_tab[0]*det_J_r;


	// Integral t2    

	ksi1=(2.-(1.+xg_k_1[ig])*xg_k_2[jg])/4.;
	ksi2=(1.+xg_k_1[ig])*(1.+xg_k_2[jg])/4.;

	shapf=shapef(ksi1, ksi2);    

	jac_tab=Jacobi(i_region, ielemq, ksi1, ksi2);//, node);
	// differenct kernels for greens / ma and ms 
	if (g0_ma1_ms2==0)
	{
	kern_tab=kernel(i_region, ielemq, ksi1, ksi2, xp, yp, zp);//, node);
	}
	if (g0_ma1_ms2==1)
	{
	kern_tab=Dkernel_dma(i_region, ielemq, ksi1, ksi2, xp, yp, zp);//, node);
	}
	if (g0_ma1_ms2==2)
	{
	kern_tab=Dkernel_dms(i_region, ielemq, ksi1, ksi2, xp, yp, zp);//, node);
	}

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
  if(nodep==node(ielemq,2)){  
    for(int ig=0; ig<n_gauss_o; ig++){

      det_J_r=(1.+xg_k_1[ig])/8.;
         
      for(int jg=0; jg<n_gauss_o; jg++){

	// Integral t1

	ksi1=(1.-xg_k_1[ig])/2.; 
	ksi2=(1.+xg_k_1[ig])*(1.-xg_k_2[jg])/4.; 

	shapf=shapef(ksi1, ksi2);    

	jac_tab=Jacobi(i_region, ielemq, ksi1, ksi2);//, node);
	// differenct kernels for greens / ma and ms 
	if (g0_ma1_ms2==0)
	{
	kern_tab=kernel(i_region, ielemq, ksi1, ksi2, xp, yp, zp);//, node);
	}
	if (g0_ma1_ms2==1)
	{
	kern_tab=Dkernel_dma(i_region, ielemq, ksi1, ksi2, xp, yp, zp);//, node);
	}
	if (g0_ma1_ms2==2)
	{
	kern_tab=Dkernel_dms(i_region, ielemq, ksi1, ksi2, xp, yp, zp);//, node);
	}

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
  if(nodep==node(ielemq,3)){ 
    for(int ig=0; ig<n_gauss_o; ig++){

      det_J_r=(1.+xg_k_1[ig])/16.;
     
      for(int jg=0; jg<n_gauss_o; jg++){
	
	// Integral t1

	ksi1=(2.+(1.+xg_k_1[ig])*xg_k_2[jg])/4.; 
	ksi2=(1.-xg_k_1[ig])/4.; 

	shapf=shapef(ksi1, ksi2);    

	jac_tab=Jacobi(i_region, ielemq, ksi1, ksi2);//, node);
	// differenct kernels for greens / ma and ms 
	if (g0_ma1_ms2==0)
	{
	kern_tab=kernel(i_region, ielemq, ksi1, ksi2, xp, yp, zp);//, node);
	}
	if (g0_ma1_ms2==1)
	{
	kern_tab=Dkernel_dma(i_region, ielemq, ksi1, ksi2, xp, yp, zp);//, node);
	}
	if (g0_ma1_ms2==2)
	{
	kern_tab=Dkernel_dms(i_region, ielemq, ksi1, ksi2, xp, yp, zp);//, node);
	}

	kern1 = kern_tab[0]*jac_tab[1]+kern_tab[1]*jac_tab[2]+
	  kern_tab[2]*jac_tab[3];//Simon!

	a_aux[ic] += kern1*wg_k[ig]*wg_k[jg]*shapf[ic]*
	  jac_tab[0]*det_J_r;//Simon!
	  
	a_aux[icd]+= kern_tab[3]*wg_k[ig]*wg_k[jg]*shapf[ic]*
	  jac_tab[0]*det_J_r;
	//a_aux[icd]+= wg_k[ig]*wg_k[jg]*jac_tab[0]*det_J_r;

	// Integral t2    

	ksi1=(1.-xg_k_1[ig])/4.; 
	ksi2=(2.-(1.+xg_k_1[ig])*xg_k_2[jg])/4.; 

	shapf=shapef(ksi1, ksi2);    

	jac_tab=Jacobi(i_region, ielemq, ksi1, ksi2);//, node);
	// differenct kernels for greens / ma and ms 
	if (g0_ma1_ms2==0)
	{
	kern_tab=kernel(i_region, ielemq, ksi1, ksi2, xp, yp, zp);//, node);
	}
	if (g0_ma1_ms2==1)
	{
	kern_tab=Dkernel_dma(i_region, ielemq, ksi1, ksi2, xp, yp, zp);//, node);
	}
	if (g0_ma1_ms2==2)
	{
	kern_tab=Dkernel_dms(i_region, ielemq, ksi1, ksi2, xp, yp, zp);//, node);
	}

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
  if(nodep==node(ielemq,4)){  
    for(int ig=0; ig<n_gauss_o; ig++){
 
      det_J_r=(1.+xg_k_1[ig])/8.;
     
      for(int jg=0; jg<n_gauss_o; jg++){

	// Integral t1

	ksi1=(1.+xg_k_1[ig])*(1.+xg_k_2[jg])/4.;
	ksi2=(1.-xg_k_1[ig])/2.; 

	shapf=shapef(ksi1, ksi2);    

	jac_tab=Jacobi(i_region, ielemq, ksi1, ksi2);//, node);
	// differenct kernels for greens / ma and ms 
	if (g0_ma1_ms2==0)
	{
	kern_tab=kernel(i_region, ielemq, ksi1, ksi2, xp, yp, zp);//, node);
	}
	if (g0_ma1_ms2==1)
	{
	kern_tab=Dkernel_dma(i_region, ielemq, ksi1, ksi2, xp, yp, zp);//, node);
	}
	if (g0_ma1_ms2==2)
	{
	kern_tab=Dkernel_dms(i_region, ielemq, ksi1, ksi2, xp, yp, zp);//, node);
	}
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
  if(nodep==node(ielemq,5)){  
    for(int ig=0; ig<n_gauss_o; ig++){

      det_J_r=(1.+xg_k_1[ig])/16.;
     
      for(int jg=0; jg<n_gauss_o; jg++){

	// Integral t1

	ksi1=(1.+xg_k_1[ig])*(1.+xg_k_2[jg])/4.;
	ksi2=(1.-xg_k_1[ig])/4.; 

	shapf=shapef(ksi1, ksi2);    

	jac_tab=Jacobi(i_region, ielemq, ksi1, ksi2);//, node);
	// differenct kernels for greens / ma and ms 
	if (g0_ma1_ms2==0)
	{
	kern_tab=kernel(i_region, ielemq, ksi1, ksi2, xp, yp, zp);//, node);
	}
	if (g0_ma1_ms2==1)
	{
	kern_tab=Dkernel_dma(i_region, ielemq, ksi1, ksi2, xp, yp, zp);//, node);
	}
	if (g0_ma1_ms2==2)
	{
	kern_tab=Dkernel_dms(i_region, ielemq, ksi1, ksi2, xp, yp, zp);//, node);
	}

	kern1 = kern_tab[0]*jac_tab[1]+kern_tab[1]*jac_tab[2]+
	  kern_tab[2]*jac_tab[3];//Simon!

	a_aux[ic] += kern1*wg_k[ig]*wg_k[jg]*shapf[ic]*
	  jac_tab[0]*det_J_r;//Simon!
	  
	a_aux[icd]+= kern_tab[3]*wg_k[ig]*wg_k[jg]*shapf[ic]*
	  jac_tab[0]*det_J_r;
	//a_aux[icd]+= wg_k[ig]*wg_k[jg]*jac_tab[0]*det_J_r;

	// Integral t2 
  
	ksi1=(1.+xg_k_1[ig])*(1.-xg_k_2[jg])/4.;
	ksi2=(2.+(1.+xg_k_1[ig])*xg_k_2[jg])/4.; 

	shapf=shapef(ksi1, ksi2);    

	jac_tab=Jacobi(i_region, ielemq, ksi1, ksi2);//, node);
	// differenct kernels for greens / ma and ms 
	if (g0_ma1_ms2==0)
	{
	kern_tab=kernel(i_region, ielemq, ksi1, ksi2, xp, yp, zp);//, node);
	}
	if (g0_ma1_ms2==1)
	{
	kern_tab=Dkernel_dma(i_region, ielemq, ksi1, ksi2, xp, yp, zp);//, node);
	}
	if (g0_ma1_ms2==2)
	{
	kern_tab=Dkernel_dms(i_region, ielemq, ksi1, ksi2, xp, yp, zp);//, node);
	}

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

  //-------------------------------------------------
  delete []norder;
  return a_aux;
}
//####################################################################################





CVector bemdiff::singu4(int i_region, int ic, int nodep, int ielemq, double xp, double yp, double zp)
{
  int icd, *norder = new int[dimension]; 
  
  double ksi1, ksi2, det_J_r;
  std::complex<double> kern1, kern2;
  
  static RVector jac_tab(dimension+1), shapf(dimension);
  static CVector kern_tab(i_kern), a_aux(dim2);
  //==================================================================
  //  purpose: Singular integrals in triangular element
  //  element is divided by 4 sub-triangles
  //==================================================================
  
  //matrix a_aux initialization
  for(int i=0;i<dim2;++i){
    a_aux[i]=std::complex<double>(0.,0.);
  }
  icd= ic + dimension;
  
  //  node 0------------------------------------------------
  if(nodep==node(ielemq,0)){ 
    for(int ig=0; ig<n_gauss_o; ig++){
      
      det_J_r=(1.+xg_k_1[ig])*0.03125;
      
      for(int jg=0; jg<n_gauss_o; jg++){
	
	// Integral t0 - singular point
	
	ksi1=det_J_r*(1.-xg_k_2[jg])*4.; 
	ksi2=det_J_r*(1.+xg_k_2[jg])*4.; 
	
	shapf=shapef(ksi1, ksi2);    
	
	jac_tab=Jacobi(i_region, ielemq, ksi1, ksi2);//, node);
	
	kern_tab=kernel(i_region, ielemq, ksi1, ksi2, xp, yp, zp);//, node);
	
	a_aux[icd]+= kern_tab[3]*wg_k[ig]*wg_k[jg]*shapf[ic]*
	  jac_tab[0]*det_J_r;
	
	//a_aux[icd]+= wg_k[ig]*wg_k[jg]*jac_tab[0]*det_J_r;
	
	//fileout<<"ielemq= "<<ielemq<<"\tic= "<<ic<<"\ta_aux[icd]= "<<
	//  a_aux[icd]<<endl;
      }//jg loop
    }//ig loop
    
    for(int ig=0; ig<n_gauss; ig++){
      
      // Integral t1
      
      ksi1=(1.+xg_1[ig])*.5; ksi2=xg_2[ig]*.5; 
      
      shapf=shapef(ksi1, ksi2);    
      
      jac_tab=Jacobi(i_region, ielemq, ksi1, ksi2);//, node);
      
      kern_tab=kernel(i_region, ielemq, ksi1, ksi2, xp, yp, zp);//, node);
      
      a_aux[icd]+= kern_tab[3]*wg[ig]*shapf[ic]*jac_tab[0]*.25;//.25=det_J_r;
      
      //a_aux[icd]+= wg[ig]*jac_tab[0]*.25;
      
      //fileout<<"ielemq= "<<ielemq<<"\tic= "<<ic<<"\ta_aux[icd]= "<<
      //  a_aux[icd]<<endl;
      
      // Integral t2
      
      ksi1=(1.-xg_2[ig])*.5; ksi2=(xg_1[ig]+xg_2[ig])*.5; 
      
      shapf=shapef(ksi1, ksi2);    
      
      jac_tab=Jacobi(i_region, ielemq, ksi1, ksi2);//, node);
      
      kern_tab=kernel(i_region, ielemq, ksi1, ksi2, xp, yp, zp);//, node);
      
      a_aux[icd]+= kern_tab[3]*wg[ig]*shapf[ic]*jac_tab[0]*.25;//.25=det_J_r;
      
	//a_aux[icd]+= wg[ig]*jac_tab[0]*.25;
      
	//fileout<<"ielemq= "<<ielemq<<"\tic= "<<ic<<"\ta_aux[icd]= "<<
	//  a_aux[icd]<<endl;
      
	// Integral t3
      
      ksi1=xg_1[ig]*.5; ksi2=(1.+xg_2[ig])*.5; 
      
      shapf=shapef(ksi1, ksi2);    
      
      jac_tab=Jacobi(i_region, ielemq, ksi1, ksi2);//, node);
      
      kern_tab=kernel(i_region, ielemq, ksi1, ksi2, xp, yp, zp);//, node);
      
      a_aux[icd]+= kern_tab[3]*wg[ig]*shapf[ic]*jac_tab[0]*.25;//.25=det_J_r;
      
      //a_aux[icd]+= wg[ig]*jac_tab[0]*.25;
      
      //fileout<<"ielemq= "<<ielemq<<"\tic= "<<ic<<"\ta_aux[icd]= "<<
      //  a_aux[icd]<<endl;
    }
  }
  
  //node 1-------------------------------------------------
  if(nodep==node(ielemq,1)){ 
    for(int ig=0; ig<n_gauss_o; ig++){
      
      det_J_r=(1.+xg_k_1[ig])*0.03125;
     
      
      for(int jg=0; jg<n_gauss_o; jg++){
	
	// Integral t0
	
	ksi1=(1.-xg_k_1[ig])*.25;
	ksi2=det_J_r*(1.-xg_k_2[jg])*4.;

	shapf=shapef(ksi1, ksi2);    
	
	jac_tab=Jacobi(i_region, ielemq, ksi1, ksi2);//, node);
	kern_tab=kernel(i_region, ielemq, ksi1, ksi2, xp, yp, zp);//, node);

	a_aux[icd]+= kern_tab[3]*wg_k[ig]*wg_k[jg]*shapf[ic]*
	  jac_tab[0]*det_J_r;
	//a_aux[icd]+= wg_k[ig]*wg_k[jg]*jac_tab[0]*det_J_r;
	
	
	// Integral t1    
	
	ksi1=det_J_r*(3.-xg_k_2[jg])*4.+(1.-xg_k_1[ig])*.25;
	ksi2=det_J_r*(1.+xg_k_2[jg])*4.;
	
	shapf=shapef(ksi1, ksi2);    
	
	jac_tab=Jacobi(i_region, ielemq, ksi1, ksi2);//, node);
	kern_tab=kernel(i_region, ielemq, ksi1, ksi2, xp, yp, zp);//, node);
	
	a_aux[icd]+= kern_tab[3]*wg_k[ig]*wg_k[jg]*shapf[ic]*
	  jac_tab[0]*det_J_r;
	//a_aux[icd]+= wg_k[ig]*wg_k[jg]*jac_tab[0]*det_J_r;

	
	//fileout<<"ielemq= "<<ielemq<<"\tic= "<<ic<<"\ta_aux[icd]= "<<
	//  a_aux[icd]<<endl;
	
	// Integral t2    
	
	ksi1=det_J_r*(1.-xg_k_2[jg])*4.+(1.-xg_k_1[ig])*.25;
	ksi2=det_J_r*8.;
	
	shapf=shapef(ksi1, ksi2);    
	
	jac_tab=Jacobi(i_region, ielemq, ksi1, ksi2);//, node);
	kern_tab=kernel(i_region, ielemq, ksi1, ksi2, xp, yp, zp);//, node);
	
	a_aux[icd]+= kern_tab[3]*wg_k[ig]*wg_k[jg]*shapf[ic]*
	  jac_tab[0]*det_J_r;
	//a_aux[icd]+= wg_k[ig]*wg_k[jg]*jac_tab[0]*det_J_r;
	
	
	//fileout<<"ielemq= "<<ielemq<<"\tic= "<<ic<<"\ta_aux[icd]= "<<
	//  a_aux[icd]<<endl;
	
      }
    }
    
    for(int ig=0; ig<n_gauss; ig++){
      
      // Integral t3
      
    	ksi1=xg_1[ig]*.5; ksi2=(1.+xg_2[ig])*.5; 
	
	shapf=shapef(ksi1, ksi2);    
	
	jac_tab=Jacobi( i_region, ielemq, ksi1, ksi2);//, node);
	
	kern_tab=kernel(i_region, ielemq, ksi1, ksi2, xp, yp, zp);//, node);
	
	a_aux[icd]+= kern_tab[3]*wg[ig]*shapf[ic]*jac_tab[0]*.25;//.25=det_J_r;
	
	//a_aux[icd]+= wg[ig]*jac_tab[0]*.25;
	
	//fileout<<"ielemq= "<<ielemq<<"\tic= "<<ic<<"\ta_aux[icd]= "<<
	//  a_aux[icd]<<endl;
    }
  }
  //node 2-------------------------------------------------
  if(nodep==node(ielemq,2)){  
    for(int ig=0; ig<n_gauss_o; ig++){
      
      det_J_r=(1.+xg_k_1[ig])*0.03125;
      
      for(int jg=0; jg<n_gauss_o; jg++){
	
	// Integral t1
	
	ksi1=(3.-xg_k_1[ig])/4.; 
	ksi2=det_J_r*(1.-xg_k_2[jg])*4.; 
	
	shapf=shapef(ksi1, ksi2);    
	
	jac_tab=Jacobi(i_region, ielemq, ksi1, ksi2);//, node);
	kern_tab=kernel(i_region, ielemq, ksi1, ksi2, xp, yp, zp);//, node);
	
	a_aux[icd]+= kern_tab[3]*wg_k[ig]*wg_k[jg]*shapf[ic]*
	  jac_tab[0]*det_J_r;
	//a_aux[icd]+= wg_k[ig]*wg_k[jg]*jac_tab[0]*det_J_r;
	
	//fileout<<"ielemq= "<<ielemq<<"\tic= "<<ic<<"\ta_aux[icd]= "<<
	// a_aux[icd]<<endl;
      }
    }
    
    for(int ig=0; ig<n_gauss; ig++){
      
      // Integral t0
      
      ksi1=xg_1[ig]*.5; ksi2=xg_2[ig]*.5; 
      
      shapf=shapef(ksi1, ksi2);    
      
      jac_tab=Jacobi(i_region, ielemq, ksi1, ksi2);//, node);
      
      kern_tab=kernel(i_region, ielemq, ksi1, ksi2, xp, yp, zp);//, node);

      a_aux[icd]+= kern_tab[3]*wg[ig]*shapf[ic]*jac_tab[0]*.25;//.25=det_J_r;
      
      // Integral t2
      
      ksi1=(1.-xg_2[ig])*.5; ksi2=(xg_1[ig]+xg_2[ig])*.5; 

      shapf=shapef(ksi1, ksi2);    
      
      jac_tab=Jacobi(i_region, ielemq, ksi1, ksi2);//, node);
      
      kern_tab=kernel(i_region, ielemq, ksi1, ksi2, xp, yp, zp);//, node);
      
      a_aux[icd]+= kern_tab[3]*wg[ig]*shapf[ic]*jac_tab[0]*.25;//.25=det_J_r;
      
      // Integral t3
      
      ksi1=xg_1[ig]*.5; ksi2=(1.+xg_2[ig])*.5; 
	
      shapf=shapef(ksi1, ksi2);    
      
      jac_tab=Jacobi(i_region, ielemq, ksi1, ksi2);//, node);
      
      kern_tab=kernel(i_region, ielemq, ksi1, ksi2, xp, yp, zp);//, node);
      
      a_aux[icd]+= kern_tab[3]*wg[ig]*shapf[ic]*jac_tab[0]*.25;//.25=det_J_r;
      
      //a_aux[icd]+= wg[ig]*jac_tab[0]*.25;
      
      //fileout<<"ielemq= "<<ielemq<<"\tic= "<<ic<<"\ta_aux[icd]= "<<
      //  a_aux[icd]<<endl;
    }
  }
  
  //node 3-------------------------------------------------
  if(nodep==node(ielemq,3)){ 
    for(int ig=0; ig<n_gauss_o; ig++){
      
      det_J_r=(1.+xg_k_1[ig])*0.03125;
      
      for(int jg=0; jg<n_gauss_o; jg++){
	
	// Integral t1
	
	ksi2=(1.-xg_k_1[ig])/4.; 
	ksi1=det_J_r*(3.+xg_k_2[jg])*4.+ksi2; 
	
	shapf=shapef(ksi1, ksi2);    
	
	jac_tab=Jacobi(i_region, ielemq, ksi1, ksi2);//, node);
	kern_tab=kernel(i_region, ielemq, ksi1, ksi2, xp, yp, zp);//, node);
	
	a_aux[icd]+= kern_tab[3]*wg_k[ig]*wg_k[jg]*shapf[ic]*
	  jac_tab[0]*det_J_r;
	//a_aux[icd]+= wg_k[ig]*wg_k[jg]*jac_tab[0]*det_J_r;

	// Integral t2    
	
	ksi1=det_J_r*(1.+xg_k_2[jg])*4.+(1.-xg_k_1[ig])*.25;
	ksi2=det_J_r*(1.-xg_k_2[jg])*4.+(1.-xg_k_1[ig])*.25;
	
	shapf=shapef(ksi1, ksi2);    
	
	jac_tab=Jacobi(i_region, ielemq, ksi1, ksi2);//, node);
	kern_tab=kernel(i_region, ielemq, ksi1, ksi2, xp, yp, zp);//, node);
	
	a_aux[icd]+= kern_tab[3]*wg_k[ig]*wg_k[jg]*shapf[ic]*
	  jac_tab[0]*det_J_r;
	//a_aux[icd]+= wg_k[ig]*wg_k[jg]*jac_tab[0]*det_J_r;
	

	//fileout<<"ielemq= "<<ielemq<<"\tic= "<<ic<<"\ta_aux[icd]= "<<
	//  a_aux[icd]<<endl;
	
	// Integral t3    
	
	ksi1=(1.-xg_k_1[ig])/4.; 
	ksi2=det_J_r*(3.-xg_k_2[jg])*4.+ksi1; 
	
	shapf=shapef(ksi1, ksi2);    
	
	jac_tab=Jacobi(i_region, ielemq, ksi1, ksi2);//, node);
	kern_tab=kernel(i_region, ielemq, ksi1, ksi2, xp, yp, zp);//, node);
	
	a_aux[icd]+= kern_tab[3]*wg_k[ig]*wg_k[jg]*shapf[ic]*
	  jac_tab[0]*det_J_r;
	//a_aux[icd]+= wg_k[ig]*wg_k[jg]*jac_tab[0]*det_J_r;
	
	
	//fileout<<"ielemq= "<<ielemq<<"\tic= "<<ic<<"\ta_aux[icd]= "<<
	//  a_aux[icd]<<endl;
	
      }
    }
    
    for(int ig=0; ig<n_gauss; ig++){
      
      // Integral t0
      
      ksi1=xg_1[ig]*.5; ksi2=xg_2[ig]*.5; 

      shapf=shapef(ksi1, ksi2);    
      
      jac_tab=Jacobi(i_region, ielemq, ksi1, ksi2);//, node);
      
      kern_tab=kernel(i_region, ielemq, ksi1, ksi2, xp, yp, zp);//, node);
      
      a_aux[icd]+= kern_tab[3]*wg[ig]*shapf[ic]*jac_tab[0]*.25;//.25=det_J_r;
    }
  }
  
  //node 4-------------------------------------------------
  if(nodep==node(ielemq,4)){  
    for(int ig=0; ig<n_gauss_o; ig++){
      
      det_J_r=(1.+xg_k_1[ig])*0.03125;
      
      for(int jg=0; jg<n_gauss_o; jg++){
	
	// Integral t3
	
	ksi1=det_J_r*(1.+xg_k_2[jg])*4.;
	ksi2=(3.-xg_k_1[ig])*.25; 
	
	shapf=shapef(ksi1, ksi2);    
	
	jac_tab=Jacobi(i_region, ielemq, ksi1, ksi2);//, node);
	kern_tab=kernel(i_region, ielemq, ksi1, ksi2, xp, yp, zp);//, node);
	
	a_aux[icd]+= kern_tab[3]*wg_k[ig]*wg_k[jg]*shapf[ic]*
	  jac_tab[0]*det_J_r;
	//a_aux[icd]+= wg_k[ig]*wg_k[jg]*jac_tab[0]*det_J_r;
	
	//fileout<<"ielemq= "<<ielemq<<"\tic= "<<ic<<"\ta_aux[icd]= "<<
        // a_aux[icd]<<endl;
      }
    }
    
    for(int ig=0; ig<n_gauss; ig++){
      
      // Integral t0
      
      ksi1=xg_1[ig]*.5; ksi2=xg_2[ig]*.5; 

      shapf=shapef(ksi1, ksi2);    
      
      jac_tab=Jacobi(i_region, ielemq, ksi1, ksi2);//, node);
      
      kern_tab=kernel(i_region, ielemq, ksi1, ksi2, xp, yp, zp);//, node);
      
      a_aux[icd]+= kern_tab[3]*wg[ig]*shapf[ic]*jac_tab[0]*.25;//.25=det_J_r;
      
      //a_aux[icd]+= wg[ig]*jac_tab[0]*.25;

      //fileout<<"ielemq= "<<ielemq<<"\tic= "<<ic<<"\ta_aux[icd]= "<<
	//  a_aux[icd]<<endl;
      
	// Integral t1
      
      ksi1=(1.+xg_1[ig])*.5; ksi2=xg_2[ig]*.5; 
      
      shapf=shapef(ksi1, ksi2);   
      
      jac_tab=Jacobi(i_region, ielemq, ksi1, ksi2);//, node);
      
      kern_tab=kernel(i_region, ielemq, ksi1, ksi2, xp, yp, zp);//, node);
      
      a_aux[icd]+= kern_tab[3]*wg[ig]*shapf[ic]*jac_tab[0]*.25;//.25=det_J_r;
      
      //a_aux[icd]+= wg[ig]*jac_tab[0]*.25;
      
      //fileout<<"ielemq= "<<ielemq<<"\tic= "<<ic<<"\ta_aux[icd]= "<<
      //  a_aux[icd]<<endl;
      
      // Integral t2
      
      ksi1=(1.-xg_2[ig])*.5; ksi2=(xg_1[ig]+xg_2[ig])*.5; 
      
      shapf=shapef(ksi1, ksi2);    
      
      jac_tab=Jacobi(i_region, ielemq, ksi1, ksi2);//, node);
      
      kern_tab=kernel(i_region, ielemq, ksi1, ksi2, xp, yp, zp);//, node);
      
      a_aux[icd]+= kern_tab[3]*wg[ig]*shapf[ic]*jac_tab[0]*.25;//.25=det_J_r;
      
      //a_aux[icd]+= wg[ig]*jac_tab[0]*.25;
	
      //fileout<<"ielemq= "<<ielemq<<"\tic= "<<ic<<"\ta_aux[icd]= "<<
      //  a_aux[icd]<<endl;
    }
  }

  //node 5-------------------------------------------------
  if(nodep==node(ielemq,5)){  
    for(int ig=0; ig<n_gauss_o; ig++){
      
      det_J_r=(1.+xg_k_1[ig])*0.03125;
      
      for(int jg=0; jg<n_gauss_o; jg++){
	
	// Integral t0
	
	ksi1=det_J_r*(1.+xg_k_2[jg])*4.;
	ksi2=(1.-xg_k_1[ig])*.25; 
	
	shapf=shapef(ksi1, ksi2);    
	
	jac_tab=Jacobi(i_region, ielemq, ksi1, ksi2);//, node);
	kern_tab=kernel(i_region, ielemq, ksi1, ksi2, xp, yp, zp);//, node);

	a_aux[icd]+= kern_tab[3]*wg_k[ig]*wg_k[jg]*shapf[ic]*
	  jac_tab[0]*det_J_r;
	//a_aux[icd]+= wg_k[ig]*wg_k[jg]*jac_tab[0]*det_J_r;
	
	// Integral t2 
	
	ksi1=det_J_r*8.;
	ksi2=det_J_r*(1.+xg_k_2[jg])*4.+(1.-xg_k_1[ig])/4.; 
	
	shapf=shapef(ksi1, ksi2);    
	
	jac_tab=Jacobi(i_region, ielemq, ksi1, ksi2);//, node);
	kern_tab=kernel(i_region, ielemq, ksi1, ksi2, xp, yp, zp);//, node);
	
	a_aux[icd]+= kern_tab[3]*wg_k[ig]*wg_k[jg]*shapf[ic]*
	  jac_tab[0]*det_J_r;
	//a_aux[icd]+= wg_k[ig]*wg_k[jg]*jac_tab[0]*det_J_r;
	
	//fileout<<"ielemq= "<<ielemq<<"\tic= "<<ic<<"\ta_aux[icd]= "<<
	//  a_aux[icd]<<endl;
	
	
	// Integral t3 
	
	ksi1=det_J_r*(1.-xg_k_2[jg])*4.;
	ksi2=det_J_r*(3.+xg_k_2[jg])*4.+(1.-xg_k_1[ig])*.25;
	
	shapf=shapef(ksi1, ksi2);    
	
	jac_tab=Jacobi(i_region,  ielemq, ksi1, ksi2);//, node);
	kern_tab=kernel(i_region, ielemq, ksi1, ksi2, xp, yp, zp);//, node);
	
	a_aux[icd]+= kern_tab[3]*wg_k[ig]*wg_k[jg]*shapf[ic]*
	  jac_tab[0]*det_J_r;
	//a_aux[icd]+= wg_k[ig]*wg_k[jg]*jac_tab[0]*det_J_r;
	
	//fileout<<"ielemq= "<<ielemq<<"\tic= "<<ic<<"\ta_aux[icd]= "<<
	//  a_aux[icd]<<endl;
      }
    }

    for(int ig=0; ig<n_gauss; ig++){
      
      // Integral t1
      
      ksi1=(1.+xg_1[ig])*.5; ksi2=xg_2[ig]*.5; 

      shapf=shapef(ksi1, ksi2);    
      
      jac_tab=Jacobi(i_region, ielemq, ksi1, ksi2);//, node);
      
      kern_tab=kernel(i_region, ielemq, ksi1, ksi2, xp, yp, zp);//, node);
      
      a_aux[icd]+= kern_tab[3]*wg[ig]*shapf[ic]*jac_tab[0]*.25;//.25=det_J_r;

      //a_aux[icd]+= wg[ig]*jac_tab[0]*.25;
      
      //fileout<<"ielemq= "<<ielemq<<"\tic= "<<ic<<"\ta_aux[icd]= "<<
      //  a_aux[icd]<<endl;
    }
  }
  
  //-------------------------------------------------
  delete []norder;
  return a_aux;
}

//####################################################################################

CVector bemdiff::nonsin_mat(int ic, int ielemq, double xp, double yp, double zp ,int g0_ma1_ms2)
{
   //==================================================================
  // purpose: to calculate nonsingular integrals 
  // 
  //==================================================================
  int icd;
  double  a_wsj;
  std::complex<double> kern1;
  static CVector kern_tab(i_kern),a_aux(dim2);
  
  //matrix a_aux initialization
  for(int ix=0;ix<dim2;++ix)
    {
      a_aux[ix]=std::complex<double>(0.,0.);
    }
  
  icd=ic+dimension;
 
  for(int ig=0; ig<n_gauss; ++ig)
    {
      // Use the precalculated shapef
      double  shapf[6]={shapf0[ig],shapf1[ig],shapf2[ig],shapf3[ig],shapf4[ig],shapf5[ig]};
      
      // Calculate the kernel with the precalculated shapef

      if (g0_ma1_ms2==0)
	{
	kern_tab= kernel_mat( ielemq, ig , xp, yp, zp);
	}
      if (g0_ma1_ms2==1)
	{
	kern_tab= Dkernel_dma_mat( ielemq, ig , xp, yp, zp);
	}
      if (g0_ma1_ms2==2)
	{
	kern_tab= Dkernel_dms_mat( ielemq, ig , xp, yp, zp);
	}


      
    //  kern_tab = kernel_mat( ielemq, ig , xp, yp, zp);
      
      kern1 = kern_tab[0]*x_nm1(ielemq,ig) + kern_tab[1]*y_nm1(ielemq,ig) + kern_tab[2]*z_nm1(ielemq,ig);
    
      //for the left hand side
      a_wsj=shapf[ic]*det_jac1(ielemq,ig);//already multiplied by wg[ig]
      a_aux[ic] += kern1*a_wsj;

      //for the right hand side
      a_aux[icd]+= kern_tab[3]*a_wsj;
    
    }
  return a_aux;
}


//####################################################################################



//####################################################################################

CVector bemdiff::nonsin(int ic, int ielemq, double xp, double yp, double zp)
{
  
  int icd;
  double  a_wsj;
  std::complex<double> kern1;
  static RVector jac_tab(dimension+1);//, shapf(dimension);
  static CVector kern_tab(i_kern),a_aux(dim2);
   
  //==================================================================
  // purpose: to calculate nonsingular integrals 
  //  node is in local numerating system
  //==================================================================
  
  //matrix a_aux initialization
  for(int ix=0;ix<dim2;++ix)
    {
      a_aux[ix]=std::complex<double>(0.,0.);
    }
  icd=ic+dimension;
  //made the variable ksi1, ksi2 equal to the ordinary Gaussian quadrature
  for(int ig=0; ig<n_gauss; ++ig)
    {
      // Lomef : here I got rid of the calculation of shape function in every iter
      //  ksi1=xg_1[ig]; ksi2=xg_2[ig];
      // shapf=shapef(ksi1, ksi2);
      //  kern_tab=kernel(i_region, ielemq, ksi1, ksi2, xp, yp, zp);//, node);
    

      // Use the precalculated shapef
      double  shapf[6]={shapf0[ig],shapf1[ig],shapf2[ig],shapf3[ig],shapf4[ig],shapf5[ig]};
   
      // Calculate the kernel with the precalculated shapef
      kern_tab = kernel_mat( ielemq, ig , xp, yp, zp);
      
      kern1 = kern_tab[0]*x_nm1(ielemq,ig) + kern_tab[1]*y_nm1(ielemq,ig) + kern_tab[2]*z_nm1(ielemq,ig);
    
      //for the left hand side
      a_wsj=shapf[ic]*det_jac1(ielemq,ig);//already multiplied by wg[ig]
      a_aux[ic] += kern1*a_wsj;

      //for the right hand side
      a_aux[icd]+= kern_tab[3]*a_wsj;
    
    }
  return a_aux;
}


//####################################################################################



//####################################################################################




//####################################################################################



//####################################################################################








//####################################################################################

void bemdiff:: int_cal(RDenseMatrix &S_nodes, IDenseMatrix &S_elements, 
			double mua, double mus, double freq, CVector fi, CVector dfi, RDenseMatrix &int_pos, CVector &fi_in)
{
  // Optical stuff----------------------------------------------
    
  //Initialise optical parameters vectors
  D.New(1);  //diffusion coefficient
  k.New(1);  //Wave number 

  const double nu=1.;                            // Refractive index
  const double c=0.3*1.e+12/nu;                  // Speed of light in the medium in mm/ps
  double omega=2.*pi*freq;                       // Angular speed
  D[0] = 1./(3.*(mua+mus));                      // Diffusion coefficient
  k[0] = sqrt(std::complex<double>(mua*c, -omega)/(c*D[0]));  // Wave number
 
  
  // The surfaces definition------------------------------------
 
  //number of elements of  region
  n_elements = S_elements.nRows();
  
  //number of nodes of the  region
  n_nodes  = S_nodes.nRows();
  
  // node - integer matrix stored as an n_elements by dimension
   
  node = S_elements;
   
  x=S_nodes.Col(0);
  y=S_nodes.Col(1);
  z=S_nodes.Col(2); 

  // The Jan definitions--------------------------------------  
  dimension = 6;
  dim2 = 12;
  i_kern = 4;
  
  int nodeq;
  double ksi1, ksi2, x_in, y_in, z_in;
  std::complex<double>  fi_temp, dfi_temp;// fi_in,
  std::complex<double> kern1;
  




  static RVector jac_tab(dimension+1), shapf(dimension);
  static CVector kern_tab(i_kern);
 
 //==================================================================
  // purpose: to calculate the internal values of the state function
  // real case
  //==================================================================
 
  //initialization of fi_in

  //number of internal points fro calculation
  int  n_int_points  = int_pos.nRows();
  fi_in.New(n_int_points);

  //for each internal position
  for(int ipos=0;ipos<n_int_points;ipos++)
    {
      
      x_in=int_pos(ipos,0);
      y_in=int_pos(ipos,1);
      z_in=int_pos(ipos,2);
      // take each element on the boundary in turn as the field element
      for(int ielem=0;ielem<n_elements;ielem++)
	{
    
	  for(int ic=0;ic<dimension;ic++)
	    {
	      nodeq=node(ielem,ic);
	      //---------------------------------------------------------------
	      fi_temp = fi[nodeq];
	      dfi_temp=dfi[nodeq];
	      //---------------------------------------------------------------
	      //made the variables ksi1 and ksi2 equal to the ordinary Gaussian quadrature
	      for(int ig=0; ig<n_gauss_i; ++ig)
		{
		  ksi1=xg_i_1[ig]; ksi2=xg_i_2[ig];
		  shapf=shapef(ksi1, ksi2);
		  jac_tab =Jacobi(0, ielem, ksi1, ksi2);
		  
		  kern_tab=kernel(0,ielem, ksi1, ksi2, x_in, y_in, z_in);
		  
		  kern1 = kern_tab[0]*jac_tab[1]+kern_tab[1]*jac_tab[2]+ kern_tab[2]*jac_tab[3];
		  
		  fi_in[ipos] += (-kern1*fi_temp+kern_tab[3]*dfi_temp)*wg_i[ig]*jac_tab[0]*shapf[ic];
		}
	    }
	}
    }
}
 

 //####################################################################################

void bemdiff:: int_cal_mat(RDenseMatrix &S_nodes, IDenseMatrix &S_elements,
			    double mua, double mus, double freq, RDenseMatrix &int_pos, CDenseMatrix &fi_in)
{
  // Optical stuff----------------------------------------------
    
  //Initialise optical parameters vectors
  D.New(1);  //diffusion coefficient
  k.New(1);  //Wave number 

  const double nu=1.;                            // Refractive index
  const double c=0.3*1.e+12/nu;                  // Speed of light in the medium in mm/ps
  double omega=2.*pi*freq;                       // Angular speed
  D[0] = 1./(3.*(mua+mus));                      // Diffusion coefficient
  k[0] = sqrt(std::complex<double>(mua*c, -omega)/(c*D[0]));  // Wave number
 
  
  // The surfaces definition------------------------------------
 
  //number of elements of  region
  n_elements = S_elements.nRows();
  
  //number of nodes of the  region
  n_nodes  = S_nodes.nRows();
  
  // node - integer matrix stored as an n_elements by dimension
   
  node = S_elements;
   
  x=S_nodes.Col(0);
  y=S_nodes.Col(1);
  z=S_nodes.Col(2); 

  // The Jan definitions--------------------------------------  
  dimension = 6;
  dim2 = 12;
  i_kern = 4;
  
  int nodeq;
  double ksi1, ksi2, x_in, y_in, z_in;
  std::complex<double> kern1;
  
  static RVector jac_tab(dimension+1), shapf(dimension);
  static CVector kern_tab(i_kern);
 
 //==================================================================
  // purpose: to calculate the internal values of the state function
  // real case
  //==================================================================
 
  //initialization of fi_in

  //number of internal points fro calculation
  int  n_int_points  = int_pos.nRows();
  fi_in.New(n_int_points,2*n_nodes);

 
  //for each internal position
  for(int ipos=0;ipos<n_int_points;ipos++)
    {
      
      x_in=int_pos(ipos,0);
      y_in=int_pos(ipos,1);
      z_in=int_pos(ipos,2);
      
      // take each element on the boundary in turn as the field element
      for(int ielem=0;ielem<n_elements;ielem++)
	{
    	  for(int ic=0;ic<dimension;ic++)
	    {
	      nodeq=node(ielem,ic);
	   
	      //made the variables ksi1 and ksi2 equal to the ordinary Gaussian quadrature
	      for(int ig=0; ig<n_gauss_i; ++ig)
		{
		  ksi1=xg_i_1[ig]; ksi2=xg_i_2[ig];
		  shapf=shapef(ksi1, ksi2);
		  jac_tab =Jacobi(0, ielem, ksi1, ksi2);
		  
		  kern_tab=kernel(0,ielem, ksi1, ksi2, x_in, y_in, z_in);
		  
		  kern1 = kern_tab[0]*jac_tab[1]+kern_tab[1]*jac_tab[2]+ kern_tab[2]*jac_tab[3];
		  
		  fi_in(ipos,nodeq) += ( -kern1   ) * wg_i[ig] * jac_tab[0] * shapf[ic];
		  fi_in(ipos,nodeq+n_nodes) += (  kern_tab[3] ) * wg_i[ig] * jac_tab[0] * shapf[ic];
		}
	    }
	}
    }
}
 
