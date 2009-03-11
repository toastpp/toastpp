#ifndef __BEMDIFF_H
#define __BEMDIFF_H

// Symbol import/export direction
#ifdef BEMDIFFLIB_IMPLEMENTATION
#define BEMDIFFLIB DLLEXPORT
#else
#define BEMDIFFLIB DLLIMPORT
#endif

class BEMDIFFLIB bemdiff {
  
 public:

  bemdiff();
  ~bemdiff
();


 void single_region_3D(RDenseMatrix &S_nodes, IDenseMatrix &S_elements,double mua, double mus, double freq, CDenseMatrix &a, CDenseMatrix &b );
 void single_region_ma(RDenseMatrix &S_nodes, IDenseMatrix &S_elements,double mua, double mus, double freq, CDenseMatrix &a, CDenseMatrix &b );
 void single_region_ms(RDenseMatrix &S_nodes, IDenseMatrix &S_elements,double mua, double mus, double freq, CDenseMatrix &a, CDenseMatrix &b );
 void single_region_refidx_3D(RDenseMatrix &S_nodes, IDenseMatrix &S_elements,double mua, double mus, double freq, double nu, CDenseMatrix &a, CDenseMatrix &b );
 void int_cal(RDenseMatrix &S_nodes, IDenseMatrix &S_elements,double mua, double mus, double freq, CVector fi, CVector dfi, RDenseMatrix &int_pos, CVector &fi_in);
 void int_cal_mat(RDenseMatrix &S_nodes, IDenseMatrix &S_elements,double mua, double mus, double freq, RDenseMatrix &int_pos, CDenseMatrix  &fi_in);
 

 private:
  
  //mesh definitions
  IDenseMatrix node;//the mesh net. 
  RDenseMatrix shapd_mat;  
  int n_nodes, n_elements, dimension, dim2, i_kern;
  RVector x, y, z;// the coordinates of the nodes. 
  RDenseMatrix det_jac1, x_nm1, y_nm1, z_nm1;

  //OPTICAL PROPERTIES
  RVector D, ma, ms;
  CVector k, iwc;

  //Usefull functions 
 
  inline RVector shapef(double ksi1, double ksi2); 
  inline RVector shaped(double ksi1, double ksi2); 
  inline RVector shape0f(double ksi1, double ksi2); 
  inline RVector shape0d(); 
  RVector Jacobi(int i_region,int ielem, double ksi1, double ksi2);
  RVector Jacobi_mat(int ielem, int ig);
  CVector kernel(int i_region, int ielem, double ksi1, double ksi2,double xp, double yp, double zp);
  CVector Dkernel_dma(int i_region, int ielem, double ksi1, double ksi2,double xp, double yp, double zp);
  CVector Dkernel_dms(int i_region, int ielem, double ksi1, double ksi2,double xp, double yp, double zp);
  CVector kernel_mat( int ielem,int ig,double xp, double yp, double zp);
  CVector Dkernel_dma_mat( int ielem,int ig,double xp, double yp, double zp);
  CVector Dkernel_dms_mat( int ielem,int ig,double xp, double yp, double zp);


  CVector singu2(int i_region, int ic, int nodep, int ielemq, double xp, double yp, double zp, int g0_ma1_ms2);
  CVector singu4(int i_region, int ic, int nodep, int ielemq, double xp, double yp, double zp);
  CVector nonsin( int ic, int ielemq, double xp, double yp, double zp);
  CVector nonsin_mat( int ic, int ielemq, double xp, double yp, double zp, int g0_ma1_ms2 );
  //  CVector dfi_arrang(CVector x_sol);
  // CVector fi_arrang(CVector x_sol);


};
#endif
