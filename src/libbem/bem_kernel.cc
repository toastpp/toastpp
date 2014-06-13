#include "bem_kernel.h"

CVector BEM_Kernel_Helmholtz::Calculate (BEM_Element *el, Point2D &loc, const Point3D &load)
{
  
  //==================================================================
  //  purpose: Calculate the kernals (Green's Function)
  //  Diffusion equation
  //==================================================================
  int iic; 
  double xq, yq, zq, rpq1, rpq2, xpxq, ypyq, zpzq;
  std::complex<double> krpq1, k_aux, e_aux;
  static RVector shapf;
  static CVector kern_tab(4);

  // local coordinate system
  xq = 0.;    yq = 0.;    zq = 0.;
  shapf = el->ShapeF (loc);
  Point *nd = el->Surface()->NodeList();

  for(int ic=0; ic<el->nNode();++ic){
	  iic = el->NodeIndex(ic);
      xq += shapf[ic]*nd[iic][0];
      yq += shapf[ic]*nd[iic][1];
      zq += shapf[ic]*nd[iic][2];
  }

  double xp = load[0], yp = load[1], zp = load[2];
  xpxq=xq-xp; ypyq=yq-yp; zpzq=zq-zp;
  rpq2=1./(xpxq*xpxq+ypyq*ypyq+zpzq*zpzq);	// 1/|r-r'|^2
  rpq1=sqrt(rpq2);							// 1/|r-r'|
  
  krpq1= -k/rpq1;							// - k *|r-r'|
  e_aux=exp(krpq1);							// exp( - k *|r-r'|)
  k_aux=-(k+rpq1)*e_aux*pi4*rpq2;			

  // kernel1 ->  kern_tab[0]; kern_tab[1]; kern_tab[2]; 
  kern_tab[0] = k_aux*xpxq;
  kern_tab[1] = k_aux*ypyq;
  kern_tab[2] = k_aux*zpzq;
  
  // the kernel2 -> kern_tab[3] 
  kern_tab[3] = e_aux*rpq1*pi4;
  
  return kern_tab;
}



CVector BEM_Kernel_Domega::Calculate (BEM_Element *el, Point2D &loc, const Point3D &load)
{
  
  //==================================================================
  //  The kernel is the derivative of the Green's function with respect
  //  to frequency.
  //  Diffusion equation
  //==================================================================
  int iic; 
  double xq, yq, zq, xpxq, ypyq, zpzq;
  double dst, dst2;
  std::complex<double> krpq1, k_aux, e_aux, arg;
  static RVector shapf;
  static CVector kern_tab(4);
  static std::complex<double> iunit = std::complex<double>(0,1);

  // local coordinate system
  xq = 0.;    yq = 0.;    zq = 0.;
  shapf = el->ShapeF (loc);
  Point *nd = el->Surface()->NodeList();

  for(int ic=0; ic<el->nNode();++ic){
	  iic = el->NodeIndex(ic);
      xq += shapf[ic]*nd[iic][0];
      yq += shapf[ic]*nd[iic][1];
      zq += shapf[ic]*nd[iic][2];
  }

  double xp = load[0], yp = load[1], zp = load[2];
  xpxq=xq-xp; ypyq=yq-yp; zpzq=zq-zp;

  std::complex<double> dkdomega = 1; // TODO iunit / (2*c*kappa*k);

  dst2 = (xpxq*xpxq+ypyq*ypyq+zpzq*zpzq);
  dst  = sqrt(dst2);
  arg = -k*dst;

  k_aux = k/dst*pi4*dkdomega*exp(arg);

  // kernel1 ->  kern_tab[0]; kern_tab[1]; kern_tab[2]; 
  kern_tab[0] = k_aux*xpxq;
  kern_tab[1] = k_aux*ypyq;
  kern_tab[2] = k_aux*zpzq;
  
  // the kernel2 -> kern_tab[3] 
  // MS 24.4.09: is this required?
  //kern_tab[3] = e_aux*rpq1*pi4;
  kern_tab[3] = 0.0;

  return kern_tab;
}



CVector BEM_Kernel_Dk::Calculate (BEM_Element *el, Point2D &loc, const Point3D &load)
{
  
  //==================================================================
  //  The kernel is the derivative of the Green's function with respect
  //  to the wavenumber.
  //  Diffusion equation
  //==================================================================
  int iic; 
  double xq, yq, zq, rpq1, rpq2, xpxq, ypyq, zpzq;
  std::complex<double> krpq1, k_aux, e_aux,  dk_aux;
  static RVector shapf;
  static CVector kern_tab(4);

  // local coordinate system
  xq = 0.;    yq = 0.;    zq = 0.;
  shapf = el->ShapeF (loc);
  Point *nd = el->Surface()->NodeList();

  for(int ic=0; ic<el->nNode();++ic){
	  iic = el->NodeIndex(ic);
      xq += shapf[ic]*nd[iic][0];
      yq += shapf[ic]*nd[iic][1];
      zq += shapf[ic]*nd[iic][2];
  }

  double xp = load[0], yp = load[1], zp = load[2];
  xpxq=xq-xp; ypyq=yq-yp; zpzq=zq-zp;
  rpq2=1./(xpxq*xpxq+ypyq*ypyq+zpzq*zpzq);	// 1/|r-r'|^2
  rpq1=sqrt(rpq2);							// 1/|r-r'|
  
  krpq1= -k/rpq1;							// - k *|r-r'|
  e_aux=exp(krpq1);							// exp( - k *|r-r'|)
 // k_aux=-(k+rpq1)*e_aux*pi4*rpq2;			// exp( - k *|r-r'|)* (1/4Pi)*{ }
  dk_aux= k*e_aux*rpq1*pi4;   // check that
  // kernel1 ->  kern_tab[0]; kern_tab[1]; kern_tab[2]; 
  kern_tab[0] = dk_aux*xpxq;
  kern_tab[1] = dk_aux*ypyq;
  kern_tab[2] = dk_aux*zpzq;
  
  // the kernel2 -> kern_tab[3] 
  kern_tab[3] = - e_aux*pi4;
  
  return kern_tab;
}
