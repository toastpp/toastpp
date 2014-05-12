#ifndef __TOASTSPMV_H
#define __TOASTSPMV_H

template<class T>
void Ax_spmv(int nr, int nc, T *Av, int *Ar, int *Ac,
	     const T *x, T *b);

template<class T>
void Ax_spmv(int nr, int nc, T *Av, int *Ar, int *Ac,
	     const T **x, T **b, int nrhs);

#endif // !__TOASTSPMV_H
