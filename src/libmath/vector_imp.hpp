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

// ==========================================================================
// member definitions

// --------------------------------------------------------------------------
// constructor (vector of size 0)

template<class VT>
TVector<VT>::TVector ()
{
    base_nref = 0;
    data = 0;
    size = 0;
    ext_data = false;
}

// --------------------------------------------------------------------------
// constructor (vector of size 'dim' with zero elements)

template<class VT>
TVector<VT>::TVector (int dim)
{
    dASSERT(dim >= 0, "Parameter 1 must be >= 0");
    base_nref = 0;
    ext_data = false;
    Allocate (dim);
}

// --------------------------------------------------------------------------
// constructor (uniform element assignment from scalar)

template<class VT>
TVector<VT>::TVector (int dim, const VT s)
{
    dASSERT(dim >= 0, "Parameter 1 must be >= 0");
    base_nref = 0;
    ext_data = false;
    Allocate (dim);
    *this = s;
}

// --------------------------------------------------------------------------
// constructor (assign elements from array or use array directly)

template<class VT>
TVector<VT>::TVector (int dim, VT *values, CopyMode cmode)
{
    dASSERT(dim >= 0, "Parameter 1 must be >= 0");
    base_nref = 0;
    if ((ext_data = (cmode==SHALLOW_COPY))) {
	Link (values, dim);
    } else {
	Allocate (dim);
	memcpy (data, values, dim*sizeof(VT));
    }
}

// --------------------------------------------------------------------------
// constructor (element assignment from string)

template<class VT>
TVector<VT>::TVector (int dim, const char *init)
{
    dASSERT(dim >= 0, "Parameter 1 must be >= 0");
    base_nref = 0;
    ext_data = false;
    Allocate (dim);
    std::istringstream iss(init);
    for (int i = 0; i < dim; i++)
        iss >> data[i];
}

// --------------------------------------------------------------------------
// copy constructor

template<class VT>
TVector<VT>::TVector (const TVector<VT> &v)
{
    base_nref = 0;
    ext_data = false;
    Allocate (v.size);
    Copy (v);				// copy elements from vec
}

// --------------------------------------------------------------------------
// reference constructor

template<class VT>
TVector<VT>::TVector (const TVector<VT> &v, int ofs, int dim)
{
    dASSERT(ofs >= 0 && ofs < v.Dim(), "Parameter 1 index out of range");
    dASSERT(dim >= 0, "Parameter 3 must be >= 0");
    dASSERT(ofs+dim <= v.Dim(),
	"Data block of reference vector must be contained in original vector");
    base_nref = 0;
    if ((ext_data = v.ext_data)) {
	Link (v.data + ofs, dim);
	// raw link to v's data buffer without reference counting
	// i.e. v must remain valid throughout the lifetime of this
    } else {
	Link (v, ofs, dim);
	// this increments v's reference counter
    }
}

// --------------------------------------------------------------------------
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

// --------------------------------------------------------------------------
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

// --------------------------------------------------------------------------
// Link to an external data buffer

template<class VT>
void TVector<VT>::Link (VT *values, int dim)
{
    dASSERT (values, "Invalid data buffer");
    dASSERT (dim >= 0, "Invalid buffer size");
    data = base_data = values;
    size = dim;
    base_size = 0;
    base_nref = 0;
    ext_data = true;
}

// --------------------------------------------------------------------------
// Remove the vector's link to the data block

template<class VT>
void TVector<VT>::Unlink ()
{
    if (!ext_data && base_nref) {
	if (--(*base_nref) == 0)   // was last link
	    free (base_nref);      // -> deallocate base block
	base_nref = 0;             // invalidate base pointer
    }
    data = 0;		           // invalidate data pointer
    size = 0;			   // set dimension zero
    ext_data = false;
}

// --------------------------------------------------------------------------
// Link vector to a different data block

template<class VT>
void TVector<VT>::Relink (const TVector<VT> &v)
{
    Unlink ();				// remove existing link
    Link (v);                           // link to v's data block
}

// --------------------------------------------------------------------------
// Relink vector to part of a different data block

template<class VT>
void TVector<VT>::Relink (const TVector<VT> &v, int ofs, int dim)
{
    Unlink ();				// remove existing link
    Link (v, ofs, dim);                 // link into v's data block
}

// --------------------------------------------------------------------------
// Relink vector to an external buffer

template<class VT>
void TVector<VT>::Relink (VT *values, int dim)
{
    Unlink();
    Link (values, dim);
}

// --------------------------------------------------------------------------
// allocate a data block for the vector. This function must be called by
// a constructor, or after a call to Unlink()

template<class VT>
void TVector<VT>::Allocate (int dim)
{
    dASSERT(!base_nref, "Data block present. Use Unlink first.");
    if (dim) {
	char *base = (char*)malloc (2*sizeof(int) + dim*sizeof(VT));
	dASSERT(base, "Memory allocation failed.");
	memset (base, 0, 2*sizeof(int) + dim*sizeof(VT));
	base_nref = (int*)base;			 // 1st int: refcount
	base_size = base_nref + 1;		 // 2nd int: array size
	data = base_data = (VT*)(base_nref + 2); // data are offset by 2
	size = *base_size = dim;		 // init array size
	*base_nref = 1;				 // init ref count
    } else {
	data = 0;
	size = 0;
    }
}

// --------------------------------------------------------------------------

template<class VT>
void TVector<VT>::Clear ()
{
    memset ((void*)data, 0, size*sizeof(VT));
    // warning: this assumes that "(VT)0" is represented by a sequence of
    // "0" bytes - not necessarily correct
}

// --------------------------------------------------------------------------
// Vector->Vector copy: copy from v to *this

#ifdef USE_BLAS_LEVEL1
template<>
inline void TVector<double>::Copy (const TVector<double> &v)
{
    const int incr = 1;
    if (size != v.size) New (v.size); // reallocate
    dcopy_(size, v.data, incr, data, incr);
}
template<>
inline void TVector<float>::Copy (const TVector<float> &v)
{
    const int incr = 1;
    if (size != v.size) New (v.size); // reallocate
    scopy_(size, v.data, incr, data, incr);
}
#endif // USE_BLAS_LEVEL1

template<class VT>
void TVector<VT>::Copy (const TVector<VT> &v)
{
    if (size != v.size) New (v.size); // reallocate
    memcpy (data, v.data, size * sizeof (VT));
}

// --------------------------------------------------------------------------
// partial Vector->Vector copy: copy part of v into *this

#ifdef USE_BLAS_LEVEL1
template<>
inline void TVector<double>::Copy (const TVector<double> &v,
    int tofs, int sofs, int n)
{
    const int incr = 1;
    if (n < 0) n = v.size - sofs;
    if (n > size - tofs) n = size - tofs;
    dcopy_(n, v.data+sofs, incr, data+tofs, incr);
}
template<>
inline void TVector<float>::Copy (const TVector<float> &v,
    int tofs, int sofs, int n)
{
    const int incr = 1;
    if (n < 0) n = v.size - sofs;
    if (n > size - tofs) n = size - tofs;
    scopy_(n, v.data+sofs, incr, data+tofs, incr);
}
#endif // USE_BLAS_LEVEL1

template<class VT>
void TVector<VT>::Copy (const TVector<VT> &v,
    int tofs, int sofs, int n)
{
    if (n < 0) n = v.size - sofs;
    if (n > size - tofs) n = size - tofs;
    memcpy (data+tofs, v.data+sofs, n * sizeof (VT));
}

// --------------------------------------------------------------------------

template<>
inline bool TVector<std::complex<double> >::Clip (std::complex<double> vmin, std::complex<double> vmax)
{
    // toast::complex version: clip real and imaginary parts separately

    bool clip = false;
    double vmin_re = vmin.real();
    double vmin_im = vmin.imag();
    double vmax_re = vmax.real();
    double vmax_im = vmax.imag();

    for (int i = 0; i < size; i++) {
      double d_re = data[i].real();
      double d_im = data[i].imag();
      bool d_clip = false;
      if (d_re < vmin_re) d_re = vmin_re, d_clip = true;
      if (d_re > vmax_re) d_re = vmax_re, d_clip = true;
      if (d_im < vmin_im) d_im = vmin_im, d_clip = true;
      if (d_im > vmax_im) d_im = vmax_im, d_clip = true;
      if (d_clip) {
	  data[i] = std::complex<double>(d_re,d_im);
	  clip = true;
      }
    }
    return clip;
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

// --------------------------------------------------------------------------

template<class VT>
VT &TVector<VT>::operator[] (int i) const
{
    dASSERT(i >= 0 && i < size, "Index out of range");
    return data[i];
}

// --------------------------------------------------------------------------

template<class VT>
TVector<VT> &TVector<VT>::operator= (VT s)
{
    for (int i = 0; i < size; i++) data[i] = s;
    return *this;
}

// --------------------------------------------------------------------------

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
    const int incr = 1;
    const double scale = 1.0;
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
    const int incr = 1;
    const float scale = 1.0f;
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
    const int incr = 1;
    const double scale = -1.0;
    TVector<double> tmp(*this);
    daxpy_((int&)size, scale, v.data, incr, tmp.data, incr);
    return tmp;
}
template<>
MATHLIB TVector<float> TVector<float>::operator- (const TVector<float> &v) const
{
    dASSERT(size == v.size, "Vectors have different size.");
    const int incr = 1;
    const float scale = -1.0f;
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
    //dASSERT(s, "Attempt to divide by zero");
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

template<class VT>
TVector<VT> &TVector<VT>::operator+= (const TVector<VT> &v)
{
    dASSERT(size == v.size, "Vectors have different size.");
    for (int i = 0; i < size; i++) data[i] += v.data[i];
    return *this;
}

template<class VT>
TVector<VT> &TVector<VT>::operator+= (const VT &s)
{
    for (int i = 0; i < size; i++) data[i] += s;
    return *this;
}

template<class VT>
TVector<VT> &TVector<VT>::operator-= (const TVector<VT> &v)
{
    dASSERT(size == v.size, "Vectors have different size.");
    for (int i = 0; i < size; i++) data[i] -= v.data[i];
    return *this;
}

template<class VT>
TVector<VT> &TVector<VT>::operator-= (const VT &s)
{
    for (int i = 0; i < size; i++) data[i] -= s;
    return *this;
}

template<class VT>
TVector<VT> &TVector<VT>::operator*= (const TVector<VT> &v)
{
    dASSERT(size == v.size, "Vectors have different size.");
    for (int i = 0; i < size; i++) data[i] *= v.data[i];
    return *this;
}

// --------------------------------------------------------------------------

#ifdef USE_BLAS_LEVEL1
template<>
inline TVector<double> &TVector<double>::operator*= (const double &s)
{
    const int incr = 1;
    dscal_(size, s, data, incr);
    return *this;
}
template<>
inline TVector<float> &TVector<float>::operator*= (const float &s)
{
    const int incr = 1;
    sscal_(size, s, data, incr);
    return *this;
}
template<>
inline TVector<std::complex<double> > &TVector<std::complex<double> >::operator*=
(const std::complex<double> &s)
{
    const int incr = 1;
    zscal_(size, s, data, incr);
    return *this;
}
#endif // USE_BLAS_LEVEL1

template<class VT>
TVector<VT> &TVector<VT>::operator*= (const VT &s)
{
    for (int i = 0; i < size; i++) data[i] *= s;
    return *this;
}

// --------------------------------------------------------------------------

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

// --------------------------------------------------------------------------

template<class VT>
TVector<VT> &TVector<VT>::operator/= (const VT &s)
{
    dASSERT (s, "Attempt to divide by zero");
    *this *= ((VT)1/s);
    return *this;
}

// --------------------------------------------------------------------------

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

// --------------------------------------------------------------------------

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

// --------------------------------------------------------------------------

template<class VT>
bool operator!= (const TVector<VT> &v1, const TVector<VT> &v2)
{
    int dim = v1.Dim();
    if (dim != v2.Dim()) return true;
    for (int i = 0; i < dim; i++)
        if (v1[i] != v2[i]) return true;
    return false;
}

// --------------------------------------------------------------------------

template<class VT>
TVector<VT> inv (const TVector<VT> &v)
{
    const VT one = (VT)1;
    TVector<VT> tmp(v.size);
    for (int i = 0; i < v.size; i++) {
        //dASSERT (v[i], "Attempt to divide by zero");
        tmp[i] = one/v[i];
    }
    return tmp;
}

// --------------------------------------------------------------------------

template<class VT>
TVector<VT> sqr (const TVector<VT> &v)
{
    TVector<VT> tmp(v.size);
    for (int i = 0; i < v.size; i++) tmp[i] = v[i]*v[i];
    return tmp;
}

// --------------------------------------------------------------------------

template<class VT>
TVector<VT> sqrt (const TVector<VT> &v)
{
    TVector<VT> tmp(v.size);
    for (int i = 0; i < v.size; i++) tmp[i] = sqrt (v[i]);
    return tmp;
}

// --------------------------------------------------------------------------

template<class VT>
TVector<VT> log (const TVector<VT> &v)
{
    TVector<VT> tmp(v.size);
    for (int i = 0; i < v.size; i++) tmp[i] = log (v[i]);
    return tmp;
}

// --------------------------------------------------------------------------

template<class VT>
TVector<VT> exp (const TVector<VT> &v)
{
    TVector<VT> tmp(v.size);
    for (int i = 0; i < v.size; i++) tmp[i] = exp (v[i]);
    return tmp;
}

// --------------------------------------------------------------------------

template<class VT>
TVector<VT> pow (const TVector<VT> &v, const VT &s)
{
    TVector<VT> tmp(v.size);
    for (int i = 0; i < v.size; i++) tmp[i] = pow (v[i], s);
    return tmp;
}

// --------------------------------------------------------------------------

template<class VT>
TVector<VT> conj (const TVector<VT> &v)
{
    TVector<VT> tmp(v.size);
    for (int i = 0; i < v.size; i++) tmp[i] = toast::conj (v[i]);
    return tmp;
}

// --------------------------------------------------------------------------

template<class VT>
TVector<VT> sin (const TVector<VT> &v)
{
    TVector<VT> tmp(v.size);
    for (int i = 0; i < v.size; i++) tmp[i] = sin (v[i]);
    return tmp;
}

// --------------------------------------------------------------------------

template<class VT>
TVector<VT> cos (const TVector<VT> &v)
{
    TVector<VT> tmp(v.size);
    for (int i = 0; i < v.size; i++) tmp[i] = cos (v[i]);
    return tmp;
}

// --------------------------------------------------------------------------

template<class VT>
TVector<VT> tan (const TVector<VT> &v)
{
    TVector<VT> tmp(v.size);
    for (int i = 0; i < v.size; i++) tmp[i] = tan (v[i]);
    return tmp;
}

// --------------------------------------------------------------------------

template<class VT>
TVector<VT> asin (const TVector<VT> &v)
{
    TVector<VT> tmp(v.size);
    for (int i = 0; i < v.size; i++) tmp[i] = asin (v[i]);
    return tmp;
}

// --------------------------------------------------------------------------

template<class VT>
TVector<VT> acos (const TVector<VT> &v)
{
    TVector<VT> tmp(v.size);
    for (int i = 0; i < v.size; i++) tmp[i] = acos (v[i]);
    return tmp;
}

// --------------------------------------------------------------------------

template<class VT>
TVector<VT> atan (const TVector<VT> &v)
{
    TVector<VT> tmp(v.size);
    for (int i = 0; i < v.size; i++) tmp[i] = atan (v[i]);
    return tmp;
}

// --------------------------------------------------------------------------

template<class VT>
TVector<VT> sinh (const TVector<VT> &v)
{
    TVector<VT> tmp(v.size);
    for (int i = 0; i < v.size; i++) tmp[i] = sinh (v[i]);
    return tmp;
}

// --------------------------------------------------------------------------

template<class VT>
TVector<VT> cosh (const TVector<VT> &v)
{
    TVector<VT> tmp(v.size);
    for (int i = 0; i < v.size; i++) tmp[i] = cosh (v[i]);
    return tmp;
}

// --------------------------------------------------------------------------

template<class VT>
TVector<VT> tanh (const TVector<VT> &v)
{
    TVector<VT> tmp(v.size);
    for (int i = 0; i < v.size; i++) tmp[i] = tanh (v[i]);
    return tmp;
}

// --------------------------------------------------------------------------

template<class VT>
TVector<VT> asinh (const TVector<VT> &v)
{
    TVector<VT> tmp(v.size);
    for (int i = 0; i < v.size; i++) tmp[i] = (VT)log ( v[i] + sqrt(sqr(v[i])+1.) );
    return tmp;
}

// --------------------------------------------------------------------------

template<class VT>
TVector<VT> acosh (const TVector<VT> &v)
{
    TVector<VT> tmp(v.size);
    for (int i = 0; i < v.size; i++) tmp[i] = (VT)log ( v[i] + sqrt(sqr(v[i])-1.) );
    return tmp;
}

// --------------------------------------------------------------------------

template<class VT>
TVector<VT> atanh (const TVector<VT> &v)
{
    TVector<VT> tmp(v.size);
    for (int i = 0; i < v.size; i++) tmp[i] =(VT)0.5 * (VT)log ( (1.+v[i]) / (1.-v[i]) );
    return tmp;
}

// --------------------------------------------------------------------------

#ifdef USE_BLAS_LEVEL1 // interface to BLAS level1 xDOT functions
template<>
inline double dot (const TVector<double> &v1, const TVector<double> &v2)
{
    dASSERT (v1.size == v2.size, "Vector dimensions incompatible");
    const int incr = 1;
    int size = v1.size;
    return ddot_(size, v1.data, incr, v2.data, incr);
}
template<>
inline float dot (const TVector<float> &v1, const TVector<float> &v2)
{
    dASSERT (v1.size == v2.size, "Vector dimensions incompatible");
    const int incr = 1;
    int size = v1.size;
    return sdot_(size, v1.data, incr, v2.data, incr);
}
template<>
inline std::complex<double> dot (const TVector<std::complex<double> > &v1, const TVector<std::complex<double> > &v2)
{
    dASSERT (v1.size == v2.size, "Vector dimensions incompatible");
    const int incr = 1;
    int size = v1.size;
    dcomplex z = zdotu_(&size, (dcomplex*)v1.data, &incr,
			       (dcomplex*)v2.data, &incr);
    return std::complex<double>(z.r, z.i);
}
#endif // USE_BLAS_LEVEL1

template<class VT>
VT dot (const TVector<VT> &v1, const TVector<VT> &v2)
{
    dASSERT (v1.size == v2.size, "Vector dimensions incompatible");
    VT d = (VT)0;
    for (int i = 0; i < v1.size; i++) d += v1[i] * v2[i];
    return d;
}

// --------------------------------------------------------------------------

template<>
inline std::complex<double> doth (const TVector<std::complex<double> > &v1, const TVector<std::complex<double> > &v2)
{
    dASSERT (v1.size == v2.size, "Vector dimensions incompatible");
#ifdef USE_BLAS_LEVEL1
    const int incr = 1;
    int size = v1.size;
    dcomplex z = zdotc_(size, v1.data, incr, v2.data, incr);
    return std::complex<double>(z.r, z.i);
#else
    std::complex<double> d = (std::complex<double>)0;
    for (int i = 0; i < v1.size; i++) d += toast::conj(v1[i]) * v2[i];
    return d;
#endif // USE_BLAS_LEVEL1
}

template<class VT>
VT doth (const TVector<VT> &v1, const TVector<VT> &v2)
{
    return dot (v1, v2);
}

// --------------------------------------------------------------------------

template<class VT>
TVector<VT> cross (const TVector<VT> &v1, const TVector<VT> &v2)
{
    xASSERT(v1.Dim() == 3 && v2.Dim() == 3,
        "TVector::cross: invalid argument size");

    TVector<VT> res(3);
    res[0] = v1[1]*v2[2] - v1[2]*v2[1];
    res[1] = v1[2]*v2[0] - v1[0]*v2[2];
    res[2] = v1[0]*v2[1] - v1[1]*v2[0];

    return res;
}

// --------------------------------------------------------------------------
#ifdef USE_BLAS_LEVEL1
template<>
inline double l2norm (const TVector<double> &v)
{
    const int incr = 1;
    return dnrm2_((int&)v.size, v.data, incr);
}
template<>
inline double l2norm (const TVector<float> &v)
{
    const int incr = 1;
    return (double)snrm2_((int&)v.size, v.data, incr);
}
template<>
inline double l2norm (const TVector<std::complex<double> > &v)
{
    const int incr = 1;
    return dznrm2_((int&)v.size, v.data, incr);
}
#endif // USE_BLAS_LEVEL1

template<class VT>
double l2norm (const TVector<VT> &v)
{
	return sqrt (l2normsq (v));
}

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
VT prod (const TVector<VT> &v)
{
    if (v.size == 0) return 0;
    VT sum = (VT)1;
    for (int i = 0; i < v.size; i++) sum *= v[i];
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
    for (int i = 0; i < v.size; i++) sum += std::abs(v[i]);
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
        term = std::abs(v[i]);
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
TVector<VT> cat (const TVector<VT> &v1, const TVector<VT> &v2)
{
    TVector<VT> tmp (v1.size + v2.size);
    tmp.Copy (v1, 0, 0, v1.size);
    tmp.Copy (v2, v1.size, 0, v2.size);
    return tmp;
}

// Append v2 to v1 and return resulting v1
template<class VT>
TVector<VT> &append (TVector<VT> &v1, const TVector<VT> &v2)
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
        if (std::isnan(v[i].real()) || std::isnan(v[i].imag())) return true;
    return false;
}

template<>
inline bool visnan (const SCVector &v)
{
    for (int i = 0; i < v.size; i++)
        if (std::isnan(v[i].real()) || std::isnan(v[i].imag())) return true;
    return false;
}

template<class VT>
bool visnan (const TVector<VT> &v)
{
    for (int i = 0; i < v.size; i++)
	if (std::isnan(v[i])) return true;
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
