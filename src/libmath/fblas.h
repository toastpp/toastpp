// FORTRAN BLAS interface

#ifndef __FBLAS_H
#define __FBLAS_H

//typedef int integer;
typedef struct { double r, i; } dcomplex;

extern "C" {

    // BLAS level 1 interfaces

    void    dscal_  (int&, const double&, double*, int&);
    void    sscal_  (int&, const float&, float*, int&);
    void    zscal_  (int&, const toast::complex&, toast::complex*, int&);

    void    daxpy_  (int&, double&, double*, int&, double*, int&);
    void    saxpy_  (int&, float&, float*, int&, float*, int&);

    void    dcopy_  (int&, double*, int&, double*, int&);
    void    scopy_  (int&, float*, int&, float*, int&);

    double  ddot_   (int&, double*, int&, double*, int&);
    float   sdot_   (int&, float*, int&, float*, int&);
    dcomplex zdotu_  (int*, dcomplex*, int*, dcomplex*, int*);
    dcomplex zdotc_  (int&, toast::complex*, int&, toast::complex*,
			    int&);

    double  dnrm2_  (int&, double*, int&);
    float   snrm2_  (int&, float*, int&);
    double  dznrm2_ (int&, toast::complex*, int&);

    int     dsyrk_  (char&, char&, int&, int&, double&, const double*, int&,
		     double&, double*, int&);
    int     ssyrk_  (char&, char&, int&, int&, float&, const float*, int&,
		     float&, float*, int&);
    int     zsyrk_  (char&, char&, int&, int&, toast::complex&,
		     const toast::complex*, int&, toast::complex&,
		     toast::complex*, int&);

    // BLAS level 2 interfaces

    void    dgemv_  (char&, int&, int&, double&, double*, int&,
		     const double*, int&, double&, double*, int&);
    void    sgemv_  (char&, int&, int&, float&, float*, int&,
		     const float*, int&, float&, float*, int&);
    void    zgemv_  (char&, int&, int&, toast::complex&, toast::complex*, int&,
		     const toast::complex*, int&, toast::complex&,
		     toast::complex*, int&);

    // BLAS level 3 interfaces

    void    dgemm_  (char&, char&, int&, int&, int&, double&, double*, int&,
		     double*, int&, double&, double*, int&);
    void    sgemm_  (char&, char&, int&, int&, int&, float&, float*, int&,
		     float*, int&, float&, float*, int&);
    void    zgemm_  (char&, char&, int&, int&, int&, toast::complex&,
		     toast::complex*, int&,
		     toast::complex*, int&, toast::complex&, toast::complex*,
		     int&);

    // sparse BLAS toolkit interface

    void    dcsrmmz_(int&, int&, int&, int&, double&, int*, double*, int*,
		     int*, int*, const double*, int&, double&, double*, int&,
		     double*, int&);
}

#endif // !__FBLAS_H
