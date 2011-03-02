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

using namespace std;
using namespace toast;

// ==========================================================================
// member definitions

// link to the data block of `vec'. This function must be called by a
// constructor, or after a call to Unlink()
template<class VT>
void TVector<VT>::Link (const TVector<VT> &vec)
{
    dASSERT (!base_nref, Data block present. Use Unlink first.);
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
    dASSERT(ofs >= 0 && dim >= 0, Invalid arguments);
    dASSERT(ofs+dim <= vec.size, Reference exceeds index range of base vector.);
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

template<class VT>
bool TVector<VT>::Clip (VT vmin, VT vmax)
{
    bool clip = false;

    for (int i = 0; i < size; i++) {
        if      (data[i] < vmin) data[i] = vmin, clip = true;
	else if (data[i] > vmax) data[i] = vmax, clip = true;
    }
    return clip;
}

template<>
bool TVector<complex>::Clip (complex vmin, complex vmax)
{
    // complex version: clip real and imaginary parts separately

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

template<class VT>
TVector<VT> TVector<VT>::operator+ (const TVector<VT> &v) const
{
    dASSERT(size == v.size, Vectors have different size.);
    TVector<VT> tmp(size);
    for (int i = 0; i < size; i++) tmp.data[i] = data[i] + v.data[i];
    return tmp;
}
#ifdef USE_BLAS_LEVEL1 // interface to BLAS level 1 xAXPY functions
template<>
MATHLIB TVector<double> TVector<double>::operator+ (const TVector<double> &v) const
{
    dASSERT(size == v.size, Vectors have different size.);
    static int incr = 1;
    static double scale = 1.0;
    TVector<double> tmp(*this);
    daxpy_((int&)size, scale, v.data, incr, tmp.data, incr);
    return tmp;
}
template<>
MATHLIB TVector<float> TVector<float>::operator+ (const TVector<float> &v) const
{
    dASSERT(size == v.size, Vectors have different size.);
    static int incr = 1;
    static float scale = 1.0f;
    TVector<float> tmp(*this);
    saxpy_((int&)size, scale, v.data, incr, tmp.data, incr);
    return tmp;
}
#endif // USE_BLAS_LEVEL1

template<class VT>
TVector<VT> TVector<VT>::operator+ (const VT &s) const
{
    TVector<VT> tmp(size);
    for (int i = 0; i < size; i++) tmp.data[i] = data[i] + s;
    return tmp;
}

template<class VT>
MATHLIB TVector<VT> operator+ (const VT &s, const TVector<VT> &v)
{
    // pre-addition of scalar (implemented as friend)
    TVector<VT> tmp(v);
    for (int i = 0; i < v.size; i++) tmp.data[i] += s;
    return tmp;  
}

template<class VT>
TVector<VT> TVector<VT>::operator- (const TVector<VT> &v) const
{
    dASSERT(size == v.size, Vectors have different size.);
    TVector<VT> tmp(size);
    for (int i = 0; i < size; i++) tmp.data[i] = data[i] - v.data[i];
    return tmp;
}
#ifdef USE_BLAS_LEVEL1 // interface to BLAS level 1 xAXPY functions
template<>
MATHLIB TVector<double> TVector<double>::operator- (const TVector<double> &v) const
{
    dASSERT(size == v.size, Vectors have different size.);
    static int incr = 1;
    static double scale = -1.0;
    TVector<double> tmp(*this);
    daxpy_((int&)size, scale, v.data, incr, tmp.data, incr);
    return tmp;
}
template<>
MATHLIB TVector<float> TVector<float>::operator- (const TVector<float> &v) const
{
    dASSERT(size == v.size, Vectors have different size.);
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
    dASSERT(size == v.size, Vectors have different size.);
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
MATHLIB TVector<VT> operator* (const VT &s, const TVector<VT> &v)
{
    // pre-multiplication with scalar (implemented as friend)
    TVector<VT> tmp(v.size);
    for (int i = 0; i < v.size; i++) tmp.data[i] = s * v[i];
    return tmp;  
}

template<class VT>
TVector<VT> TVector<VT>::operator/ (const TVector<VT> &v) const
{
    dASSERT(size == v.size, Vectors have different size.);
    TVector<VT> tmp(size);
    for (int i = 0; i < size; i++) {
        dASSERT (v.data[i], Attempt to divide by zero);
        tmp.data[i] = data[i] / v.data[i];
    }
    return tmp;
}

template<class VT>
TVector<VT> TVector<VT>::operator/ (const VT &s) const
{
    dASSERT(s, Attempt to divide by zero);
    TVector<VT> tmp(size);
    for (int i = 0; i < size; i++) tmp.data[i] = data[i] / s;
    return tmp;
}

template<class VT>
TVector<VT> operator/ (const VT &s, const TVector<VT> &v)
{
    TVector<VT> tmp(v.size);
    for (int i = 0; i < v.size; i++) {
        dASSERT (v.data[i], Attempt to divide by zero);
        tmp.data[i] = s / v.data[i];
    }
    return tmp;  
}

#ifdef MATH_DEBUG // otherwise inline
template<class VT>
TVector<VT> &TVector<VT>::operator+= (const TVector<VT> &v)
{
    dASSERT(size == v.size, Vectors have different size.);
    for (int i = 0; i < size; i++) data[i] += v.data[i];
    return *this;
}
#endif // MATH_DEBUG

#ifdef MATH_DEBUG // otherwise inline
template<class VT>
TVector<VT> &TVector<VT>::operator-= (const TVector<VT> &v)
{
    dASSERT(size == v.size, Vectors have different size.);
    for (int i = 0; i < size; i++) data[i] -= v.data[i];
    return *this;
}
#endif // MATH_DEBUG

#ifdef MATH_DEBUG // otherwise inline
template<class VT>
TVector<VT> &TVector<VT>::operator*= (const TVector<VT> &v)
{
    dASSERT(size == v.size, Vectors have different size.);
    for (int i = 0; i < size; i++) data[i] *= v.data[i];
    return *this;
}
#endif // MATH_DEBUG

#ifdef MATH_DEBUG // otherwise inline
template<class VT>
TVector<VT> &TVector<VT>::operator/= (const TVector<VT> &v)
{
    dASSERT(size == v.size, Vectors have different size.);
    for (int i = 0; i < size; i++) {
        dASSERT (v.data[i], Attempt to divide by zero);
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
	dASSERT (index >=0 && index < size, Index out of range);
	is >> data[index];
    }
}

// ==========================================================================
// friend definitions

template<class VT>
MATHLIB bool operator== (const TVector<VT> &v1, const TVector<VT> &v2)
{
    int dim = v1.Dim();
    if (dim != v2.Dim()) return false;
    for (int i = 0; i < dim; i++)
        if (v1[i] != v2[i]) return false;
    return true;
}


template<class VT>
MATHLIB bool operator!= (const TVector<VT> &v1, const TVector<VT> &v2)
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
        dASSERT (v[i], Attempt to divide by zero);
        tmp[i] = one/v[i];
    }
    return tmp;
}

template<class VT>
MATHLIB TVector<VT> sqr (const TVector<VT> &v)
{
    TVector<VT> tmp(v.size);
    for (int i = 0; i < v.size; i++) tmp[i] = v[i]*v[i];
    return tmp;
}

template<class VT>
MATHLIB TVector<VT> sqrt (const TVector<VT> &v)
{
    TVector<VT> tmp(v.size);
    for (int i = 0; i < v.size; i++) tmp[i] = sqrt (v[i]);
    return tmp;
}

// ==========================================================================

template<class VT>
MATHLIB TVector<VT> log (const TVector<VT> &v)
{
    TVector<VT> tmp(v.size);
    for (int i = 0; i < v.size; i++) tmp[i] = log (v[i]);
    return tmp;
}

// ==========================================================================

template<class VT>
MATHLIB TVector<VT> exp (const TVector<VT> &v)
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
    dASSERT (v1.size == v2.size, Vector dimensions incompatible);
    VT d = (VT)0;
    for (int i = 0; i < v1.size; i++) d += v1[i] * v2[i];
    return d;
}
#ifdef USE_BLAS_LEVEL1 // interface to BLAS level1 xDOT functions
template<>
double dot (const TVector<double> &v1, const TVector<double> &v2)
{
    dASSERT (v1.size == v2.size, Vector dimensions incompatible);
    static int incr = 1;
    int size = v1.size;
    return ddot_(size, v1.data, incr, v2.data, incr);
}
template<>
float dot (const TVector<float> &v1, const TVector<float> &v2)
{
    dASSERT (v1.size == v2.size, Vector dimensions incompatible);
    static int incr = 1;
    int size = v1.size;
    return sdot_(size, v1.data, incr, v2.data, incr);
}
template<>
complex dot (const TVector<complex> &v1, const TVector<complex> &v2)
{
    dASSERT (v1.size == v2.size, Vector dimensions incompatible);
    static int incr = 1;
    int size = v1.size;
    dcomplex z = zdotu_(&size, (dcomplex*)v1.data, &incr,
			       (dcomplex*)v2.data, &incr);
    return complex(z.r, z.i);
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
complex doth (const TVector<complex> &v1, const TVector<complex> &v2)
{
    dASSERT (v1.size == v2.size, Vector dimensions incompatible);
#ifdef USE_BLAS_LEVEL1
    static int incr = 1;
    int size = v1.size;
    dcomplex z = zdotc_(size, v1.data, incr, v2.data, incr);
    return complex(z.r, z.i);
#else
    complex d = (complex)0;
    for (int i = 0; i < v1.size; i++) d += conj(v1[i]) * v2[i];
    return d;
#endif // USE_BLAS_LEVEL1
}
#endif // MATH_DEBUG

// ==========================================================================

template<class VT>
MATHLIB VT vmin (const TVector<VT> &v)
{
    if (!v.size) return (VT)0; // trap pathological case
    VT mn = v[0];
    for (int i = 1; i < v.size; i++)
        if (v[i] < mn) mn = v[i];
    return mn;
}

// ==========================================================================

template<class VT>
MATHLIB VT vmax (const TVector<VT> &v)
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
MATHLIB double l1norm (const TVector<VT> &v)
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

#ifdef TOAST_THREAD

template<class VT>
void *l2normsq_engine (void *context)
{
    typedef struct {
	const TVector<VT> *v;
	pthread_t ht;
	int i0, i1;
	double sum;
    } THDATA;
    THDATA *thdata = (THDATA*)context;

    double term;
    thdata->sum = 0.0;
    const TVector<VT> &v = *thdata->v;
    for (int i = thdata->i0; i < thdata->i1; i++) {
	term = norm(v[i]);
	thdata->sum += term*term;
    }
    return NULL;
}

#ifdef NEED_EXPLICIT_INSTANTIATION
template void *l2normsq_engine<double> (void *context);
template void *l2normsq_engine<float> (void *context);
template void *l2normsq_engine<complex> (void *context);
template void *l2normsq_engine<scomplex> (void *context);
template void *l2normsq_engine<int> (void *context);
#endif

template<class VT>
MATHLIB double l2normsq (const TVector<VT> &v)
{
    //static const int nthread = 2;
    static struct {
	const TVector<VT> *v;
	pthread_t ht;
	int i0, i1;
	double sum;
    } thdata[NUMTHREAD];

    int i;

    for (i = 0; i < NUMTHREAD; i++) {
	thdata[i].v = &v;
	thdata[i].i0 = (i*v.size)/NUMTHREAD;
	thdata[i].i1 = ((i+1)*v.size)/NUMTHREAD;
	pthread_create (&thdata[i].ht, NULL,
	    l2normsq_engine<VT>, (void*)(thdata+i));
    }
    double sum = 0.0;
    for (i = 0; i < NUMTHREAD; i++) {
	pthread_join (thdata[i].ht, NULL);
	sum += thdata[i].sum;
    }
    return sum;
}

#else

template<class VT>
MATHLIB double l2normsq (const TVector<VT> &v)
{
    double term, sum = 0.0;
    for (int i = 0; i < v.size; i++) {
	term = norm (v[i]);
	sum += term*term;
    }
    return sum;
}

#endif

template<class VT>
MATHLIB double linfnorm (const TVector<VT> &v)
{
    double nm, mx = 0.0;
    for (int i = 0; i < v.size; i++)
	if ((nm = norm (v[i])) > mx) mx = nm;
    return mx;
}

// Return concatenation of v1 and v2
template<class VT>
MATHLIB TVector<VT> cat (const TVector<VT> &v1, const TVector<VT> &v2)
{
    TVector<VT> tmp (v1.size + v2.size);
    tmp.Copy (v1, 0, 0, v1.size);
    tmp.Copy (v2, v1.size, 0, v2.size);
    return tmp;
}

// Append v2 to v1 and return resulting v1
template<class VT>
MATHLIB TVector<VT> &append (TVector<VT> &v1, const TVector<VT> &v2)
{
    if (v2.size) {
        TVector<VT> tmp(v1);
	v1.New (v1.size+v2.size);
	v1.Copy (tmp, 0, 0, tmp.size);
	v1.Copy (v2, tmp.size, 0, v2.size);
    }
    return v1;
}

template<class VT>
MATHLIB bool visnan (const TVector<VT> &v)
{
    for (int i = 0; i < v.size; i++)
	if (isnan(v[i])) return true;
    return false;
}

template<>
MATHLIB bool visnan (const CVector &v)
{
    for (int i = 0; i < v.size; i++)
	if (cisnan(v[i])) return true;
    return false;
}

template<>
MATHLIB bool visnan (const SCVector &v)
{
    for (int i = 0; i < v.size; i++)
	if (cisnan(v[i])) return true;
    return false;
}

template<class VT>
MATHLIB ostream &operator<< (ostream &os, const TVector<VT> &v)
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
MATHLIB istream &operator>> (istream &is, TVector<VT> &v)
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
MATHLIB TVector<double> UnfoldComplex (const RVector &v)
{
    // nothing to do for real case
    return v;
}

template<>
MATHLIB RVector UnfoldComplex (const CVector &v)
{
    int n = v.Dim();
    RVector rvec(n*2);
    RVector rvec_r(rvec,0,n); rvec_r = Re(v);
    RVector rvec_i(rvec,n,n); rvec_i = Im(v);
    return rvec;
}

template<class VT>
void TVector<VT>::Scan (const char *cbuf, int nmax)
{
    if (!nmax || nmax > size) nmax = size;
    std::istringstream iss (cbuf);
    for (int i = 0; i < size; i++) iss >> data[i];
}

// ==========================================================================
// class and friend instantiations

#ifdef NEED_EXPLICIT_INSTANTIATION

template class MATHLIB TVector<double>;
template class MATHLIB TVector<float>;
template class MATHLIB TVector<toast::complex>;
template class MATHLIB TVector<scomplex>;
template class MATHLIB TVector<int>;

template MATHLIB RVector  inv (const RVector &v);
template MATHLIB FVector  inv (const FVector &v);
template MATHLIB CVector  inv (const CVector &v);
template MATHLIB SCVector inv (const SCVector &v);
template MATHLIB IVector  inv (const IVector &v);

template MATHLIB RVector  sqr (const RVector &v);
template MATHLIB FVector  sqr (const FVector &v);
template MATHLIB CVector  sqr (const CVector &v);
template MATHLIB SCVector sqr (const SCVector &v);
template MATHLIB IVector  sqr (const IVector &v);

template MATHLIB RVector  sqrt (const RVector &v);
template MATHLIB FVector  sqrt (const FVector &v);
template MATHLIB CVector  sqrt (const CVector &v);
template MATHLIB SCVector sqrt (const SCVector &v);

template MATHLIB RVector  log (const RVector &v);
template MATHLIB FVector  log (const FVector &v);
template MATHLIB CVector  log (const CVector &v);
template MATHLIB SCVector log (const SCVector &v);

template MATHLIB RVector  exp (const RVector &v);
template MATHLIB FVector  exp (const FVector &v);
template MATHLIB CVector  exp (const CVector &v);
template MATHLIB SCVector exp (const SCVector &v);

template RVector  pow (const RVector &vec, const double &vt);
template FVector  pow (const FVector &vec, const float &vt);
template CVector  pow (const CVector &vec, const complex &vt);
template SCVector pow (const SCVector &vec, const scomplex &vt);
template IVector  pow (const IVector &vec, const int &vt);

template RVector  conj (const RVector &v);
template CVector  conj (const CVector &v);
template SCVector conj (const SCVector &v);

template RVector sin (const RVector &v);
template FVector sin (const FVector &v);

template RVector cos (const RVector &v);
template FVector cos (const FVector &v);

template RVector tan (const RVector &v);
template FVector tan (const FVector &v);

template RVector asin (const RVector &v);
template FVector asin (const FVector &v);

template RVector acos (const RVector &v);
template FVector acos (const FVector &v);

template RVector atan (const RVector &v);
template FVector atan (const FVector &v);

template RVector sinh (const RVector &v);
template FVector sinh (const FVector &v);

template RVector cosh (const RVector &v);
template FVector cosh (const FVector &v);

template RVector tanh (const RVector &v);
template FVector tanh (const FVector &v);

template RVector asinh (const RVector &v);
template FVector asinh (const FVector &v);

template RVector acosh (const RVector &v);
template FVector acosh (const FVector &v);

template RVector atanh (const RVector &v);
template FVector atanh (const FVector &v);

#ifndef USE_BLAS_LEVEL1
template double   dot (const RVector &v1, const RVector &v2);
template float    dot (const FVector &v1, const FVector &v2);
template complex  dot (const CVector &v1, const CVector &v2);
#endif // !USE_BLAS_LEVEL1
template int      dot (const IVector &v1, const IVector &v2);
template scomplex dot (const SCVector &v1, const SCVector &v2);

template double   doth (const RVector &v1, const RVector &v2);
template float    doth (const FVector &v1, const FVector &v2);
#ifndef USE_BLAS_LEVEL1
template complex  doth (const CVector &v1, const CVector &v2);
#endif // !USE_BLAS_LEVEL1
template scomplex doth (const SCVector &v1, const SCVector &v2);
template int      doth (const IVector &v1, const IVector &v2);

template MATHLIB double   vmin (const RVector &v);
template MATHLIB float    vmin (const FVector &v);
template MATHLIB complex  vmin (const CVector &v);
template MATHLIB scomplex vmin (const SCVector &v);
template MATHLIB int      vmin (const IVector &v);

template MATHLIB double   vmax (const RVector &v);
template MATHLIB float    vmax (const FVector &v);
template MATHLIB complex  vmax (const CVector &v);
template MATHLIB scomplex vmax (const SCVector &v);
template MATHLIB int      vmax (const IVector &v);

template RVector vsort (const RVector &v);
template FVector vsort (const FVector &v);
template CVector vsort (const CVector &v);
template IVector vsort (const IVector &v);

template void vsort (const RVector &v, RVector &sorted_v, IVector &sort_order);
template void vsort (const FVector &v, FVector &sorted_v, IVector &sort_order);
template void vsort (const CVector &v, CVector &sorted_v, IVector &sort_order);
template void vsort (const IVector &v, IVector &sorted_v, IVector &sort_order);


template double   sum  (const RVector &v);
template float    sum  (const FVector &v);
template complex  sum  (const CVector &v);
template scomplex sum  (const SCVector &v);
template int      sum  (const IVector &v);

template MATHLIB double   mean (const RVector &v);
template MATHLIB float    mean (const FVector &v);
template MATHLIB complex  mean (const CVector &v);
template MATHLIB scomplex mean (const SCVector &v);
template MATHLIB int      mean (const IVector &v);

template double  median (const RVector &v);
template float   median (const FVector &v);
template complex median (const CVector &v);
template int     median (const IVector &v);

template double   variance (const RVector &v);
template float    variance (const FVector &v);
template complex  variance (const CVector &v);
template scomplex variance (const SCVector &v);
template int      variance (const IVector &v);

template double   stdv (const RVector &v);
template float    stdv (const FVector &v);
template complex  stdv (const CVector &v);
template scomplex stdv (const SCVector &v);

#ifndef USE_CBLAS
template MATHLIB double   l1norm (const RVector &v);
template MATHLIB double   l1norm (const FVector &v);
#endif // !USE_CBLAS
template MATHLIB double   l1norm (const CVector &v);
template MATHLIB double   l1norm (const SCVector &v);
template MATHLIB double   l1norm (const IVector &v);

#ifndef USE_BLAS_LEVEL1
template MATHLIB double   l2norm (const RVector &v);
template MATHLIB double   l2norm (const FVector &v);
template MATHLIB double   l2norm (const CVector &v);
#endif // !USE_BLAS_LEVEL1
template MATHLIB double   l2norm (const SCVector &v);
template MATHLIB double   l2norm (const IVector &v);

template MATHLIB double   linfnorm (const RVector &v);
template MATHLIB double   linfnorm (const FVector &v);
template MATHLIB double   linfnorm (const CVector &v);
template MATHLIB double   linfnorm (const SCVector &v);
template MATHLIB double   linfnorm (const IVector &v);

double length (const RVector &vec);
double length (const FVector &vec);
double length (const CVector &vec);
double length (const SCVector &vec);
double length (const IVector &vec);

template MATHLIB double   l2normsq<double> (const RVector &vec);
template MATHLIB double   l2normsq (const FVector &vec);
template MATHLIB double   l2normsq (const CVector &vec);
template MATHLIB double   l2normsq (const SCVector &vec);
template MATHLIB double   l2normsq (const IVector &vec);

template MATHLIB RVector  &append (RVector &v1, const RVector &v2);
template MATHLIB FVector  &append (FVector &v1, const FVector &v2);
template MATHLIB CVector  &append (CVector &v1, const CVector &v2);
template MATHLIB SCVector &append (SCVector &v1, const SCVector &v2);
template MATHLIB IVector  &append (IVector &v1, const IVector &v2);

template MATHLIB RVector  cat (const RVector &v1, const RVector &v2);
template MATHLIB FVector  cat (const FVector &v1, const FVector &v2);
template MATHLIB CVector  cat (const CVector &v1, const CVector &v2);
template MATHLIB SCVector cat (const SCVector &v1, const SCVector &v2);
template MATHLIB IVector  cat (const IVector &v1, const IVector &v2);

template MATHLIB bool visnan (const RVector &v);
template MATHLIB bool visnan (const FVector &v);
template MATHLIB bool visnan (const IVector &v);

template MATHLIB RVector  operator+ (const double  &s, const RVector &v);
template MATHLIB FVector  operator+ (const float   &s, const FVector &v);
template MATHLIB CVector  operator+ (const complex &s, const CVector &v);
template MATHLIB SCVector operator+ (const scomplex &s, const SCVector &v);
template MATHLIB IVector  operator+ (const int     &s, const IVector &v);

template RVector  operator- (const double  &s, const RVector &v);
template FVector  operator- (const float   &s, const FVector &v);
template CVector  operator- (const complex &s, const CVector &v);
template SCVector operator- (const scomplex &s, const SCVector &v);
template IVector  operator- (const int     &s, const IVector &v);

template MATHLIB RVector  operator* (const double  &s, const RVector &v);
template MATHLIB FVector  operator* (const float   &s, const FVector &v);
template MATHLIB CVector  operator* (const complex &s, const CVector &v);
template MATHLIB SCVector operator* (const scomplex &s, const SCVector &v);
template MATHLIB IVector  operator* (const int     &s, const IVector &v);

template MATHLIB RVector  operator/ (const double  &s, const RVector &v);
template MATHLIB FVector  operator/ (const float   &s, const FVector &v);
template MATHLIB CVector  operator/ (const complex &s, const CVector &v);
template MATHLIB SCVector operator/ (const scomplex &s, const SCVector &v);
template MATHLIB IVector  operator/ (const int     &s, const IVector &v);

template MATHLIB bool operator== (const RVector &v1, const RVector &v2);
template MATHLIB bool operator== (const FVector &v1, const FVector &v2);
template MATHLIB bool operator== (const CVector &v1, const CVector &v2);
template MATHLIB bool operator== (const SCVector &v1, const SCVector &v2);
template MATHLIB bool operator== (const IVector &v1, const IVector &v2);

template MATHLIB bool operator!= (const RVector &v1, const RVector &v2);
template MATHLIB bool operator!= (const FVector &v1, const FVector &v2);
template MATHLIB bool operator!= (const CVector &v1, const CVector &v2);
template MATHLIB bool operator!= (const SCVector &v1, const SCVector &v2);
template MATHLIB bool operator!= (const IVector &v1, const IVector &v2);

template MATHLIB ostream  &operator<< (ostream &os, const RVector &v);
template MATHLIB ostream  &operator<< (ostream &os, const FVector &v);
template MATHLIB ostream  &operator<< (ostream &os, const CVector &v);
template MATHLIB ostream  &operator<< (ostream &os, const SCVector &v);
template MATHLIB ostream  &operator<< (ostream &os, const IVector &v);

template MATHLIB istream  &operator>> (istream &os, RVector &v);
template MATHLIB istream  &operator>> (istream &os, FVector &v);
template MATHLIB istream  &operator>> (istream &os, CVector &v);
template MATHLIB istream  &operator>> (istream &os, SCVector &v);
template MATHLIB istream  &operator>> (istream &os, IVector &v);

template double SparseDotp (const RVector &v1, idxtype *idx1, int nidx1,
			    const RVector &v2, idxtype *idx2, int nidx2,
			    int from, int to);

#endif // NEED_EXPLICIT_INSTANTIATION

/* Explicit complex conversions */
/* These ought to be friends really, except that I can't see how to do that
when using template */
MATHLIB FVector Re (const SCVector &vec)
{
    FVector tmp(vec.Dim());
    for (int i = 0; i < vec.Dim(); i++) tmp[i] = vec[i].re;
    return tmp;
}

MATHLIB RVector Re (const CVector &vec)
{
    RVector tmp(vec.Dim());
    for (int i = 0; i < vec.Dim(); i++) tmp[i] = vec[i].re;
    return tmp;
}

MATHLIB FVector Im (const SCVector &vec)
{
    FVector tmp(vec.Dim());
    for (int i = 0; i < vec.Dim(); i++) tmp[i] = vec[i].im;
    return tmp;
}

MATHLIB RVector Im (const CVector &vec)
{
    RVector tmp(vec.Dim());
    for (int i = 0; i < vec.Dim(); i++) tmp[i] = vec[i].im;
    return tmp;
}

RVector Mod (const CVector &vec)
{
    RVector tmp(vec.Dim());
    for (int i = 0; i < vec.Dim(); i++) tmp[i] = mod(vec[i]);
    return tmp;
}

MATHLIB RVector LogMod (const CVector &vec)
{
    RVector tmp(vec.Dim());
    for (int i = 0; i < vec.Dim(); i++) tmp[i] = log(mod(vec[i]));
    return tmp;
}


MATHLIB RVector Arg (const CVector &vec)
{
    RVector tmp(vec.Dim());
    for (int i = 0; i < vec.Dim(); i++) tmp[i] = arg(vec[i]);
    return tmp;
}


CVector Conj (const CVector &vec)
{
    CVector tmp(vec.Dim());
    for (int i = 0; i < vec.Dim(); i++) tmp[i] = conj(vec[i]);
    return tmp;
}
void SelfConj (const CVector &vec)
{
  /* version converts this vector to conjugate */
    for (int i = 0; i < vec.Dim(); i++) vec[i] = conj(vec[i]);
}
CVector Hadamard (const CVector &a, const CVector &b)
{
    dASSERT(a.Dim() == b.Dim(), Dimension mismatch);
    CVector tmp(a.Dim());
    for (int i = 0; i < a.Dim(); i++) tmp[i] = hadamard(a[i],b[i]);
    return tmp;
}

MATHLIB void SetReal (CVector &z, const RVector &zre)
{
    dASSERT(z.Dim() == zre.Dim(), Dimension mismatch);
    for (int i = 0; i < z.Dim(); i++) z[i].re = zre[i];
}

MATHLIB void SetImag (CVector &z, const RVector &zim)
{
    dASSERT(z.Dim() == zim.Dim(), Dimension mismatch);
    for (int i = 0; i < z.Dim(); i++) z[i].im = zim[i];
}

CVector MakeCVector (const SCVector &v)
{
    CVector c(v.Dim());
    for (int i = 0; i < v.Dim(); i++)
	c[i] = (complex)v[i];
    return c;
}
