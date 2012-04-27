// ==========================================================================
// Module mathlib
// File vector.cc
// Definition of template class TVector ('template vector')
// ==========================================================================

#define __VECTOR_CC
#define MATHLIB_IMPLEMENTATION

//#if (defined(WIN32)||defined(WIN64))
//#include "stdafx.h"
//#endif

#include <iostream>
#include <string>
#include <sstream>
#include <iomanip>
#include <math.h>
#include <cstdlib>
#include "mathlib.h"

#ifdef USE_CUDA_FLOAT
#include "toastcuda.h"
#endif

using namespace std;
using namespace toast;

// ==========================================================================
// member definitions

// link to the data block of `vec'. This function must be called by a
// constructor, or after a call to Unlink()
template<class VT>
void TVector<VT>::Link (const TVector<VT> &vec)
{
    dASSERT (!base_nref, "Data block present. Use Unlink first.");
    base_nref = vec.base_nref;      // link to vec's data block
    if (base_nref) {
	base_size = vec.base_size;  // use vec's base size pointer
	base_data = vec.base_data;  // use vec's base data pointer
	size = *base_size;          // use full data block
	data = base_data;           // set local data pointer to start of base
	(*base_nref)++;             // increment data block's reference counter
    } else {
	data = 0;
	size = 0;
    }
}

// link to part of the data block of `vec'. This function must be called by
// a constructor, or after a call to Unlink()
template<class VT>
void TVector<VT>::Link (const TVector<VT> &vec, int ofs, int dim)
{
    dASSERT(ofs >= 0 && dim >= 0, "Invalid arguments");
    dASSERT(ofs+dim <= vec.size,
	    "Reference exceeds index range of base vector.");
    base_nref = vec.base_nref;	 // link to vec's data block
    if (base_nref) {
	base_size = vec.base_size; // use vec's base size pointer
	base_data = vec.base_data; // use vec's base data pointer
	size = dim;		   // local vector dimension
	data = base_data + ofs;	   // set local data pointer to ofs from start
	(*base_nref)++;		   // increment data block's reference counter
    } else {
	data = 0;
	size = 0;
    }
}

template<>
inline bool TVector<toast::complex>::Clip (toast::complex vmin, toast::complex vmax)
{
    // toast::complex version: clip real and imaginary parts separately

    bool clip = false;

    for (int i = 0; i < size; i++) {
        if      (data[i].re < vmin.re) data[i].re = vmin.re, clip = true;
	else if (data[i].re > vmax.re) data[i].re = vmax.re, clip = true;
        if      (data[i].im < vmin.im) data[i].im = vmin.im, clip = true;
	else if (data[i].im > vmax.im) data[i].im = vmax.im, clip = true;
    }
    return clip;
}

template<class VT>
inline bool TVector<VT>::Clip (VT vmin, VT vmax)
{
    bool clip = false;

    for (int i = 0; i < size; i++) {
        if      (data[i] < vmin) data[i] = vmin, clip = true;
	else if (data[i] > vmax) data[i] = vmax, clip = true;
    }
    return clip;
}

template<class VT>
void TVector<VT>::ShiftLeft (int n)
{
    int i;
    if (n > size) n = size;
    for (i = 0; i < size-n; i++) data[i] = data[i+n];
    for (; i < size; i++) data[i] = (VT)0;
}

template<class VT>
void TVector<VT>::ShiftRight (int n)
{
    int i;
    if (n > size) n = size;
    for (i = size-1; i >= n; i--) data[i] = data[i-n];
    for (; i >= 0; i--) data[i] = (VT)0;
}

// ==========================================================================

// vector-vector addition, general template

template<class VT>
TVector<VT> TVector<VT>::operator+ (const TVector<VT> &v) const
{
    dASSERT(size == v.size, "Vectors have different size.");
    TVector<VT> tmp(size);
    for (int i = 0; i < size; i++) tmp.data[i] = data[i] + v.data[i];
    return tmp;
}
#ifdef USE_BLAS_LEVEL1 // interface to BLAS level 1 xAXPY functions
template<>
MATHLIB TVector<double> TVector<double>::operator+ (const TVector<double> &v) const
{
    dASSERT(size == v.size, "Vectors have different size.");
    static int incr = 1;
    static double scale = 1.0;
    TVector<double> tmp(*this);
    daxpy_((int&)size, scale, v.data, incr, tmp.data, incr);
    return tmp;
}
#endif // BLAS

// vector-vector addition, single precision specialisation

#ifdef USE_BLAS_LEVEL1 // interface to BLAS level 1 SAXPY functions
template<>
MATHLIB TVector<float> TVector<float>::operator+(const TVector<float> &v) const
{
    dASSERT(size == v.size, "Vectors have different size.");
    static int incr = 1;
    static float scale = 1.0f;
    TVector<float> tmp(*this);
    saxpy_((int&)size, scale, v.data, incr, tmp.data, incr);
    return tmp;
}
#endif // BLAS

// ==========================================================================

template<class VT>
TVector<VT> TVector<VT>::operator+ (const VT &s) const
{
    TVector<VT> tmp(size);
    for (int i = 0; i < size; i++) tmp.data[i] = data[i] + s;
    return tmp;
}

#ifndef USE_CUDA_FLOAT
template<class VT>
TVector<VT> operator+ (const VT &s, const TVector<VT> &v)
{
    // pre-addition of scalar (implemented as friend)
    TVector<VT> tmp(v);
    for (int i = 0; i < v.size; i++) tmp.data[i] += s;
    return tmp;  
}
#endif

template<class VT>
TVector<VT> TVector<VT>::operator- (const TVector<VT> &v) const
{
    dASSERT(size == v.size, "Vectors have different size.");
    TVector<VT> tmp(size);
    for (int i = 0; i < size; i++) tmp.data[i] = data[i] - v.data[i];
    return tmp;
}
#ifdef USE_BLAS_LEVEL1 // interface to BLAS level 1 xAXPY functions
template<>
MATHLIB TVector<double> TVector<double>::operator- (const TVector<double> &v) const
{
    dASSERT(size == v.size, "Vectors have different size.");
    static int incr = 1;
    static double scale = -1.0;
    TVector<double> tmp(*this);
    daxpy_((int&)size, scale, v.data, incr, tmp.data, incr);
    return tmp;
}
template<>
MATHLIB TVector<float> TVector<float>::operator- (const TVector<float> &v) const
{
    dASSERT(size == v.size, "Vectors have different size.");
    static int incr = 1;
    static float scale = -1.0f;
    TVector<float> tmp(*this);
    saxpy_((int&)size, scale, v.data, incr, tmp.data, incr);
    return tmp;
}
#endif // USE_BLAS_LEVEL1

template<class VT>
TVector<VT> TVector<VT>::operator- (const VT &s) const
{
    TVector<VT> tmp(size);
    for (int i = 0; i < size; i++) tmp.data[i] = data[i] - s;
    return tmp;
}

template<class VT>
TVector<VT> operator- (const VT &s, const TVector<VT> &v)
{
    // subtraction from scalar (implemented as friend)
    TVector<VT> tmp(v.size);
    for (int i = 0; i < v.size; i++) tmp.data[i] = s - v.data[i];
    return tmp;  
}

template<class VT>
TVector<VT> TVector<VT>::operator- () const
{
    TVector<VT> tmp(size);
    for (int i = 0; i < size; i++) tmp.data[i] = -data[i];
    return tmp;
}

template<class VT>
TVector<VT> TVector<VT>::operator* (const TVector<VT> &v) const
{
    dASSERT(size == v.size, "Vectors have different size.");
    TVector<VT> tmp(size);
    for (int i = 0; i < size; i++) tmp.data[i] = data[i] * v.data[i];
    return tmp;
}

template<class VT>
TVector<VT> TVector<VT>::operator* (const VT &s) const
{
    TVector<VT> tmp(size);
    for (int i = 0; i < size; i++) tmp.data[i] = data[i] * s;
    return tmp;
}

template<class VT>
TVector<VT> operator* (const VT &s, const TVector<VT> &v)
{
    // pre-multiplication with scalar (implemented as friend)
    TVector<VT> tmp(v.size);
    for (int i = 0; i < v.size; i++) tmp.data[i] = s * v[i];
    return tmp;  
}

template<class VT>
TVector<VT> TVector<VT>::operator/ (const TVector<VT> &v) const
{
    dASSERT(size == v.size, "Vectors have different size.");
    TVector<VT> tmp(size);
    for (int i = 0; i < size; i++) {
        dASSERT (v.data[i], "Attempt to divide by zero");
        tmp.data[i] = data[i] / v.data[i];
    }
    return tmp;
}

template<class VT>
TVector<VT> TVector<VT>::operator/ (const VT &s) const
{
    dASSERT(s, "Attempt to divide by zero");
    TVector<VT> tmp(size);
    for (int i = 0; i < size; i++) tmp.data[i] = data[i] / s;
    return tmp;
}

template<class VT>
TVector<VT> operator/ (const VT &s, const TVector<VT> &v)
{
    TVector<VT> tmp(v.size);
    for (int i = 0; i < v.size; i++) {
        dASSERT (v.data[i], "Attempt to divide by zero");
        tmp.data[i] = s / v.data[i];
    }
    return tmp;  
}

#ifdef MATH_DEBUG // otherwise inline
template<class VT>
TVector<VT> &TVector<VT>::operator+= (const TVector<VT> &v)
{
    dASSERT(size == v.size, "Vectors have different size.");
    for (int i = 0; i < size; i++) data[i] += v.data[i];
    return *this;
}
#endif // MATH_DEBUG

#ifdef MATH_DEBUG // otherwise inline
template<class VT>
TVector<VT> &TVector<VT>::operator-= (const TVector<VT> &v)
{
    dASSERT(size == v.size, "Vectors have different size.");
    for (int i = 0; i < size; i++) data[i] -= v.data[i];
    return *this;
}
#endif // MATH_DEBUG

#ifdef MATH_DEBUG // otherwise inline
template<class VT>
TVector<VT> &TVector<VT>::operator*= (const TVector<VT> &v)
{
    dASSERT(size == v.size, "Vectors have different size.");
    for (int i = 0; i < size; i++) data[i] *= v.data[i];
    return *this;
}
#endif // MATH_DEBUG

#ifdef MATH_DEBUG // otherwise inline
template<class VT>
TVector<VT> &TVector<VT>::operator/= (const TVector<VT> &v)
{
    dASSERT(size == v.size, "Vectors have different size.");
    for (int i = 0; i < size; i++) {
        dASSERT (v.data[i], "Attempt to divide by zero");
        data[i] /= v.data[i];
    }
    return *this;
}
#endif // MATH_DEBUG


#ifdef NODEF

// I had to throw this out because a bug in gcc-2.8.1 does not allow
// structs inside template class methods

template<class VT>
void TVector<VT>::Read (istream &is, int start, int n)
{
    char c;
    int i;
    VT dummy;

    typedef struct vtbuf_tag {
	VT buf[128];
	struct vtbuf_tag *next;
    } vtbuf;
    //typedef struct vtbuf_tag vtbuf;

    do { is.get (c); } while (c != '[');	// beginning of vector
    for (i = 0; i < start; i++) is >> dummy;	// skip first elements
    if (n >= 0) {
	New (n);
	for (i = 0; i < n; i++)	is >> data[i];	// read elements
	do { is.get (c); } while (c != ']');	// end of vector	
    } else {					// vector size unknown;
	vtbuf *buf0 = 0, *buf;
	n = 0;
	while (is >> dummy) {			// read elements to buffers
	    if (!(n%128)) {
		if (buf0) {
		    buf->next = new vtbuf;
		    buf = buf->next;
		} else {
		    buf0 = new vtbuf;
		    buf = buf0;
		}
	    }
	    buf->buf[(n++)%128] = dummy;
	}
	is.clear(0);				// reset error state
	New (n);
	buf = buf0;
	for (i = 0; i < n; i++) {
	    data[i] = buf->buf[i%128];
	    if (!((i+1)%128) || i == n-1) {
		buf0 = buf; buf = buf->next;
		delete buf0;
	    }
	}
    }
}

#else

// this is a less efficient version which copies the buffer each time
// it has filled up. This works around the gcc bug.

template<class VT>
void TVector<VT>::Read (istream &is, int start, int n)
{
    const int blocksize = 256;
    char c;
    int i, bufsize, nbufsize;
    VT dummy, *buf, *nbuf;

    do {  // search for beginning of vector
	is.get (c);
	if (!is.good()) return;  // no vector found in stream
    } while (c != '[');
    for (i = 0; i < start; i++) is >> dummy;	// skip first elements
    if (n >= 0) {
	New (n);
	for (i = 0; i < n; i++)	is >> data[i];	// read elements
	do { is.get (c); } while (c != ']');	// end of vector	
    } else {					// vector size unknown;
	n = bufsize = 0;
	while (is >> dummy) {
	    if (n >= bufsize) {
		nbufsize = bufsize + blocksize;
		nbuf = new VT[nbufsize];
		if (bufsize)
    		    memcpy (nbuf, buf, bufsize * sizeof(VT));
		//for (i = 0; i < bufsize; i++) nbuf[i] = buf[i];
		if (bufsize) delete []buf;
		buf = nbuf;
		bufsize = nbufsize;
	    }
	    buf[n++] = dummy;
	}
	is.clear();  // clear fail flag
	New (n);
	for (i = 0; i < n; i++)
	    data[i] = buf[i];
	if (n) delete []buf;
    }
}

#endif

template<class VT>
void TVector<VT>::ReadIndexed (istream &is, int n)
{
    int i, nvals, index;
    is >> nvals;
    if (n >= 0) New (n);
    for (i = 0; i < nvals; i++) {
        is >> index;
	dASSERT (index >=0 && index < size, "Index out of range");
	is >> data[index];
    }
}

// ==========================================================================
// friend definitions

template<class VT>
bool operator== (const TVector<VT> &v1, const TVector<VT> &v2)
{
    int dim = v1.Dim();
    if (dim != v2.Dim()) return false;
    for (int i = 0; i < dim; i++)
        if (v1[i] != v2[i]) return false;
    return true;
}


template<class VT>
bool operator!= (const TVector<VT> &v1, const TVector<VT> &v2)
{
    int dim = v1.Dim();
    if (dim != v2.Dim()) return true;
    for (int i = 0; i < dim; i++)
        if (v1[i] != v2[i]) return true;
    return false;
}


template<class VT>
TVector<VT> inv (const TVector<VT> &v)
{
    const static VT one = (VT)1;
    TVector<VT> tmp(v.size);
    for (int i = 0; i < v.size; i++) {
        dASSERT (v[i], "Attempt to divide by zero");
        tmp[i] = one/v[i];
    }
    return tmp;
}

template<class VT>
TVector<VT> sqr (const TVector<VT> &v)
{
    TVector<VT> tmp(v.size);
    for (int i = 0; i < v.size; i++) tmp[i] = v[i]*v[i];
    return tmp;
}

template<class VT>
TVector<VT> sqrt (const TVector<VT> &v)
{
    TVector<VT> tmp(v.size);
    for (int i = 0; i < v.size; i++) tmp[i] = sqrt (v[i]);
    return tmp;
}

// ==========================================================================

template<class VT>
TVector<VT> log (const TVector<VT> &v)
{
    TVector<VT> tmp(v.size);
    for (int i = 0; i < v.size; i++) tmp[i] = log (v[i]);
    return tmp;
}

// ==========================================================================

template<class VT>
TVector<VT> exp (const TVector<VT> &v)
{
    TVector<VT> tmp(v.size);
    for (int i = 0; i < v.size; i++) tmp[i] = exp (v[i]);
    return tmp;
}

// ==========================================================================

template<class VT>
TVector<VT> pow (const TVector<VT> &v, const VT &s)
{
    TVector<VT> tmp(v.size);
    for (int i = 0; i < v.size; i++) tmp[i] = pow (v[i], s);
    return tmp;
}

// ==========================================================================

template<class VT>
TVector<VT> conj (const TVector<VT> &v)
{
    TVector<VT> tmp(v.size);
    for (int i = 0; i < v.size; i++) tmp[i] = conj (v[i]);
    return tmp;
}

// ==========================================================================
template<class VT>
TVector<VT> sin (const TVector<VT> &v)
{
    TVector<VT> tmp(v.size);
    for (int i = 0; i < v.size; i++) tmp[i] = sin (v[i]);
    return tmp;
}
// ==========================================================================
template<class VT>
TVector<VT> cos (const TVector<VT> &v)
{
    TVector<VT> tmp(v.size);
    for (int i = 0; i < v.size; i++) tmp[i] = cos (v[i]);
    return tmp;
}
// ==========================================================================
template<class VT>
TVector<VT> tan (const TVector<VT> &v)
{
    TVector<VT> tmp(v.size);
    for (int i = 0; i < v.size; i++) tmp[i] = tan (v[i]);
    return tmp;
}
// ==========================================================================
template<class VT>
TVector<VT> asin (const TVector<VT> &v)
{
    TVector<VT> tmp(v.size);
    for (int i = 0; i < v.size; i++) tmp[i] = asin (v[i]);
    return tmp;
}
// ==========================================================================
template<class VT>
TVector<VT> acos (const TVector<VT> &v)
{
    TVector<VT> tmp(v.size);
    for (int i = 0; i < v.size; i++) tmp[i] = acos (v[i]);
    return tmp;
}
// ==========================================================================
template<class VT>
TVector<VT> atan (const TVector<VT> &v)
{
    TVector<VT> tmp(v.size);
    for (int i = 0; i < v.size; i++) tmp[i] = atan (v[i]);
    return tmp;
}
// ==========================================================================
template<class VT>
TVector<VT> sinh (const TVector<VT> &v)
{
    TVector<VT> tmp(v.size);
    for (int i = 0; i < v.size; i++) tmp[i] = sinh (v[i]);
    return tmp;
}
// ==========================================================================
template<class VT>
TVector<VT> cosh (const TVector<VT> &v)
{
    TVector<VT> tmp(v.size);
    for (int i = 0; i < v.size; i++) tmp[i] = cosh (v[i]);
    return tmp;
}
// ==========================================================================
template<class VT>
TVector<VT> tanh (const TVector<VT> &v)
{
    TVector<VT> tmp(v.size);
    for (int i = 0; i < v.size; i++) tmp[i] = tanh (v[i]);
    return tmp;
}
// ==========================================================================
template<class VT>
TVector<VT> asinh (const TVector<VT> &v)
{
    TVector<VT> tmp(v.size);
    for (int i = 0; i < v.size; i++) tmp[i] = (VT)log ( v[i] + sqrt(sqr(v[i])+1.) );
    return tmp;
}
// ==========================================================================
template<class VT>
TVector<VT> acosh (const TVector<VT> &v)
{
    TVector<VT> tmp(v.size);
    for (int i = 0; i < v.size; i++) tmp[i] = (VT)log ( v[i] + sqrt(sqr(v[i])-1.) );
    return tmp;
}
// ==========================================================================
template<class VT>
TVector<VT> atanh (const TVector<VT> &v)
{
    TVector<VT> tmp(v.size);
    for (int i = 0; i < v.size; i++) tmp[i] =(VT)0.5 * (VT)log ( (1.+v[i]) / (1.-v[i]) );
    return tmp;
}

// ==========================================================================

#ifdef MATH_DEBUG // otherwise inline
template<class VT>
VT dot (const TVector<VT> &v1, const TVector<VT> &v2)
{
    dASSERT (v1.size == v2.size, "Vector dimensions incompatible");
    VT d = (VT)0;
    for (int i = 0; i < v1.size; i++) d += v1[i] * v2[i];
    return d;
}
#ifdef USE_BLAS_LEVEL1 // interface to BLAS level1 xDOT functions
template<>
double dot (const TVector<double> &v1, const TVector<double> &v2)
{
    dASSERT (v1.size == v2.size, "Vector dimensions incompatible");
    static int incr = 1;
    int size = v1.size;
    return ddot_(size, v1.data, incr, v2.data, incr);
}
template<>
float dot (const TVector<float> &v1, const TVector<float> &v2)
{
    dASSERT (v1.size == v2.size, "Vector dimensions incompatible");
    static int incr = 1;
    int size = v1.size;
    return sdot_(size, v1.data, incr, v2.data, incr);
}
template<>
toast::complex dot (const TVector<toast::complex> &v1, const TVector<toast::complex> &v2)
{
    dASSERT (v1.size == v2.size, "Vector dimensions incompatible");
    static int incr = 1;
    int size = v1.size;
    dcomplex z = zdotu_(&size, (dcomplex*)v1.data, &incr,
			       (dcomplex*)v2.data, &incr);
    return toast::complex(z.r, z.i);
}
#endif // USE_BLAS_LEVEL1
#endif // MATH_DEBUG

// ==========================================================================

#ifdef MATH_DEBUG // otherwise inline
template<class VT>
VT doth (const TVector<VT> &v1, const TVector<VT> &v2)
{
    return dot (v1, v2);
}
template<>
toast::complex doth (const TVector<toast::complex> &v1, const TVector<toast::complex> &v2)
{
    dASSERT (v1.size == v2.size, "Vector dimensions incompatible");
#ifdef USE_BLAS_LEVEL1
    static int incr = 1;
    int size = v1.size;
    dcomplex z = zdotc_(size, v1.data, incr, v2.data, incr);
    return toast::complex(z.r, z.i);
#else
    toast::complex d = (toast::complex)0;
    for (int i = 0; i < v1.size; i++) d += conj(v1[i]) * v2[i];
    return d;
#endif // USE_BLAS_LEVEL1
}
#endif // MATH_DEBUG

// ==========================================================================

template<class VT>
VT vmin (const TVector<VT> &v)
{
    if (!v.size) return (VT)0; // trap pathological case
    VT mn = v[0];
    for (int i = 1; i < v.size; i++)
        if (v[i] < mn) mn = v[i];
    return mn;
}

// ==========================================================================

template<class VT>
VT vmax (const TVector<VT> &v)
{
    if (!v.size) return (VT)0;
    VT mx = v[0];
    for (int i = 1; i < v.size; i++)
        if (v[i] > mx) mx = v[i];
    return mx;
}
// ==========================================================================

template<class VT>
TVector<VT> vsort (const TVector<VT> &v)
{
	TVector<VT> tmp(v.size);
	VT mx= (VT) 0;
    for (int i = 0; i < v.size; i++) tmp[i] = v[i] ;
	if (!tmp.size) return tmp;
	for (int i=0; i < tmp.size; ++i)
	{
        for (int j=i-1; j>=0 && (tmp[j+1] < tmp[j]); --j)
		{
            mx = tmp[j];         /* swap a[j] and a[j+1] */
            tmp[j] = tmp[j+1];
            tmp[j+1] = mx;
        }
    }
    return tmp;
}

template<class VT>
void vsort (const TVector<VT> &v, TVector<VT> &sorted_v, IVector &sort_order)
{
    sorted_v.New(v.Dim());
    sort_order.New(v.Dim());
    
    for(int i=0; i < v.Dim(); i++)
	sort_order[i] = i;

    int mi = 0;
    VT mx= (VT) 0;
    for (int i = 0; i < v.Dim(); i++) sorted_v[i] = v[i] ;
    if (!sorted_v.Dim()) return;
    for (int i=0; i < sorted_v.Dim(); ++i)
    {
      for (int j=i-1; j>=0 && (sorted_v[j+1] < sorted_v[j]); --j)
      {
            mx = sorted_v[j];         /* swap a[j] and a[j+1] */
	    mi = sort_order[j];

            sorted_v[j] = sorted_v[j+1];
	    sort_order[j] = sort_order[j+1];

            sorted_v[j+1] = mx;
	    sort_order[j+1] = mi;
        }
    }
    return;
}

// ==========================================================================

template<class VT>
VT sum (const TVector<VT> &v)
{
    if (v.size == 0) return 0;
    VT sum = 0;
    for (int i = 0; i < v.size; i++) sum += v[i];
    return sum;
}

// ==========================================================================

template<class VT>
VT mean (const TVector<VT> &v)
{
    if (v.size == 0) return 0;
    VT sum = 0;
    for (int i = 0; i < v.size; i++) sum += v[i];
    return sum / v.size;
}

// ==========================================================================

template<class VT>
VT median (const TVector<VT> &v)
{
	TVector<VT> tmp(v.size);
	int  oddtest = (int) ceil(fmod((double) v.size ,2.0));
	VT med=(VT) 0 ;

	tmp = vsort (v);
    if (v.size == 0)  med = (VT)0;
	if (oddtest == 0) med = (v[(v.size / 2) ]+v[(v.size / 2)+1]) / ((VT) 2.0);
	if (oddtest == 1) med = v[(v.size + 1) / 2];
    return med;
}

// ==========================================================================

template<class VT>
VT variance (const TVector<VT> &v)
{
    if (v.size == 0) return 0;
    VT sum2 = 0;
    VT mn = mean (v);
    for (int i = 0; i < v.size; i++) sum2 += v[i] * v[i];
    return sum2 / v.size - mn*mn;
}

// ==========================================================================

template<class VT>
VT stdv (const TVector<VT> &v)
{
    return sqrt (variance (v));
}

// ==========================================================================

template<class VT>
double l1norm (const TVector<VT> &v)
{
    double sum = 0.0;
    for (int i = 0; i < v.size; i++) sum += norm (v[i]);
    return sum;
}
#ifdef USE_CBLAS
double l1norm (const TVector<double> &v)
{
    return cblas_dasum (v.size, v.data, 1);
}
double l1norm (const TVector<float> &v)
{
    return cblas_sasum (v.size, v.data, 1);
}
#endif // USE_CBLAS

// ==========================================================================

#if THREAD_LEVEL!=1 // threaded version implemented in vector_MT.cc

template<class VT>
double l2normsq (const TVector<VT> &v)
{
    double term, sum = 0.0;
    for (int i = 0; i < v.size; i++) {
	term = norm (v[i]);
	sum += term*term;
    }
    return sum;
}

#endif // THREAD_LEVEL!=1

template<class VT>
double linfnorm (const TVector<VT> &v)
{
    double nm, mx = 0.0;
    for (int i = 0; i < v.size; i++)
	if ((nm = norm (v[i])) > mx) mx = nm;
    return mx;
}

// Return concatenation of v1 and v2
template<class VT>
inline TVector<VT> cat (const TVector<VT> &v1, const TVector<VT> &v2)
{
    TVector<VT> tmp (v1.size + v2.size);
    tmp.Copy (v1, 0, 0, v1.size);
    tmp.Copy (v2, v1.size, 0, v2.size);
    return tmp;
}

// Append v2 to v1 and return resulting v1
template<class VT>
inline TVector<VT> &append (TVector<VT> &v1, const TVector<VT> &v2)
{
    if (v2.size) {
        TVector<VT> tmp(v1);
	v1.New (v1.size+v2.size);
	v1.Copy (tmp, 0, 0, tmp.size);
	v1.Copy (v2, tmp.size, 0, v2.size);
    }
    return v1;
}

template<>
inline bool visnan (const CVector &v)
{
    for (int i = 0; i < v.size; i++)
	if (cisnan(v[i])) return true;
    return false;
}

template<>
inline bool visnan (const SCVector &v)
{
    for (int i = 0; i < v.size; i++)
	if (cisnan(v[i])) return true;
    return false;
}

template<class VT>
bool visnan (const TVector<VT> &v)
{
    for (int i = 0; i < v.size; i++)
	if (isnan(v[i])) return true;
    return false;
}

template<class VT>
ostream &operator<< (ostream &os, const TVector<VT> &v)
{
    os << '[';
    for (int i = 0; i < v.size; i++) {
	if (i) os << ' ';
	os << setprecision(10) << v[i];
    }
    os << ']';
    return os;
}

template<class VT>
istream &operator>> (istream &is, TVector<VT> &v)
{
    const int BUFSIZE = 1024;

    typedef struct _vb {
        VT buf[1024];
        struct _vb *next;
    } vecbuf;

    char c;
    int i, j, ilen, vlen = 0, blen = 0;
    VT s;
    vecbuf *vb0 = 0, *vbc;

    do {
        if (!is.get(c)) return is;
    } while (c != '[');  // find beginning of vector
    for (;;) {
        is >> s;
	if (!is.good()) {
	    is.clear();
	    break;
	}
	if (vlen == blen) {
	    vecbuf *vbn = new vecbuf;
	    if (vlen) vbc->next = vbn;
	    else      vb0 = vbn;
	    vbc = vbn;
	    vbc->next = 0;
	    blen += BUFSIZE;
	    ilen = 0;
	}
	vbc->buf[ilen++] = s;
	vlen++;
    }

    do {
        if (!is.get(c)) return is;
    } while (c != ']');  // find end of vector

    if (vlen != v.size) v.New (vlen);
    vbc = vb0;
    for (i = j = 0; i < vlen; i++) {
        if (j == BUFSIZE) {
	    vbc = vbc->next;
	    j = 0;
	}
	v[i] = vbc->buf[j++];
    }
    while (vb0) {
        vbc = vb0;
	vb0 = vb0->next;
	delete vbc;
    }
    return is;
}

template<class VT>
VT SparseDotp (const TVector<VT> &v1, idxtype *idx1, int nidx1,
	       const TVector<VT> &v2, idxtype *idx2, int nidx2,
	       int from, int to)
{
    VT sum = 0;
    if (!nidx1 || !nidx2) return sum; // sanity check
    int imin = (idx1[0] < idx2[0] ? idx2[0] : idx1[0]);
    int imax = (idx1[nidx1-1] > idx2[nidx2-1] ? idx2[nidx2-1] : idx1[nidx1-1]);
    if (imin > imax) return sum;      // no overlap
    int i1 = 0, i2 = 0, c1, c2;
    if ((c1 = idx1[i1]) > to) return sum;
    if ((c2 = idx2[i2]) > to) return sum;
    for (;;) {
	if (c1 < c2) {
	    do {
	        if (++i1 == nidx1) return sum;
	    } while ((c1 = idx1[i1]) < c2);
	    if (c1 > to) return sum;
	} else if (c1 > c2) {
	    do {
	        if (++i2 == nidx2) return sum;
	    } while ((c2 = idx2[i2]) < c1);
	    if (c2 > to) return sum;
	}
	if (c1 == c2) {
	    if (c1 >= from) sum += v1[c1]*v2[c1];
	    if (++i1 == nidx1 || ++i2 == nidx2) return sum;
	    c1 = idx1[i1];
	    c2 = idx2[i2];
	}
    }
    // should not get here
    return sum;
}

template<>
inline TVector<double> UnfoldComplex (const RVector &v)
{
    // nothing to do for real case
    return v;
}

template<>
inline TVector<double> UnfoldComplex (const CVector &v)
{
    int n = v.Dim();
    RVector rvec(n*2);
    RVector rvec_r(rvec,0,n); rvec_r = Re(v);
    RVector rvec_i(rvec,n,n); rvec_i = Im(v);
    return rvec;
}

template<class VT>
TVector<double> UnfoldComplex (const TVector<VT> &v)
{
	xERROR("Not implemented");
	return TVector<double>();
}

template<class VT>
void TVector<VT>::Scan (const char *cbuf, int nmax)
{
    if (!nmax || nmax > size) nmax = size;
    std::istringstream iss (cbuf);
    for (int i = 0; i < size; i++) iss >> data[i];
}

/* Explicit complex conversions */
/* These ought to be friends really, except that I can't see how to do that
when using template */
inline FVector Re (const SCVector &vec)
{
    FVector tmp(vec.Dim());
    for (int i = 0; i < vec.Dim(); i++) tmp[i] = vec[i].re;
    return tmp;
}

inline RVector Re (const CVector &vec)
{
    RVector tmp(vec.Dim());
    for (int i = 0; i < vec.Dim(); i++) tmp[i] = vec[i].re;
    return tmp;
}

inline FVector Im (const SCVector &vec)
{
    FVector tmp(vec.Dim());
    for (int i = 0; i < vec.Dim(); i++) tmp[i] = vec[i].im;
    return tmp;
}

inline RVector Im (const CVector &vec)
{
    RVector tmp(vec.Dim());
    for (int i = 0; i < vec.Dim(); i++) tmp[i] = vec[i].im;
    return tmp;
}

inline RVector Mod (const CVector &vec)
{
    RVector tmp(vec.Dim());
    for (int i = 0; i < vec.Dim(); i++) tmp[i] = mod(vec[i]);
    return tmp;
}

inline RVector LogMod (const CVector &vec)
{
    RVector tmp(vec.Dim());
    for (int i = 0; i < vec.Dim(); i++) tmp[i] = log(mod(vec[i]));
    return tmp;
}


inline RVector Arg (const CVector &vec)
{
    RVector tmp(vec.Dim());
    for (int i = 0; i < vec.Dim(); i++) tmp[i] = arg(vec[i]);
    return tmp;
}


inline CVector Conj (const CVector &vec)
{
    CVector tmp(vec.Dim());
    for (int i = 0; i < vec.Dim(); i++) tmp[i] = conj(vec[i]);
    return tmp;
}
inline void SelfConj (const CVector &vec)
{
  /* version converts this vector to conjugate */
    for (int i = 0; i < vec.Dim(); i++) vec[i] = conj(vec[i]);
}
inline CVector Hadamard (const CVector &a, const CVector &b)
{
    dASSERT(a.Dim() == b.Dim(), "Dimension mismatch");
    CVector tmp(a.Dim());
    for (int i = 0; i < a.Dim(); i++) tmp[i] = hadamard(a[i],b[i]);
    return tmp;
}

inline void SetReal (CVector &z, const RVector &zre)
{
    dASSERT(z.Dim() == zre.Dim(), "Dimension mismatch");
    for (int i = 0; i < z.Dim(); i++) z[i].re = zre[i];
}

inline void SetImag (CVector &z, const RVector &zim)
{
    dASSERT(z.Dim() == zim.Dim(), "Dimension mismatch");
    for (int i = 0; i < z.Dim(); i++) z[i].im = zim[i];
}

inline CVector MakeCVector (const SCVector &v)
{
    CVector c(v.Dim());
    for (int i = 0; i < v.Dim(); i++)
	c[i] = (toast::complex)v[i];
    return c;
}
