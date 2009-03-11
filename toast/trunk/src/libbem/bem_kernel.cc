#include "bem_kernel.h"

using namespace toast;

CVector BEM_Kernel_Helmholtz::Calculate (BEM_Element *el, Point2D &loc, const Point3D &load)
{
  
  //==================================================================
  //  purpose: Calculate the kernals (Green's Function)
  //  Diffusion equation
  //==================================================================
  int iic; 
  double xq, yq, zq, rpq1, rpq2, xpxq, ypyq, zpzq;
  complex krpq1, k_aux, e_aux;
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
  rpq2=1./(xpxq*xpxq+ypyq*ypyq+zpzq*zpzq);
  rpq1=sqrt(rpq2);
  
  krpq1= -k/rpq1;
  e_aux=exp(krpq1);
  k_aux=-(k+rpq1)*e_aux*pi4*rpq2;

  // kernel1 ->  kern_tab[0]; kern_tab[1]; kern_tab[2]; 
  kern_tab[0] = k_aux*xpxq;
  kern_tab[1] = k_aux*ypyq;
  kern_tab[2] = k_aux*zpzq;
  
  // the kernel2 -> kern_tab[3] 
  kern_tab[3] = e_aux*rpq1*pi4;
  
  return kern_tab;
}
