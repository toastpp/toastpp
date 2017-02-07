// -*-C++-*-
// ==========================================================================
// Module mathlib
// File vector.h
// Declaration of template class TVector ('template vector')
//
// The following typedefs for specific template types have been defined
// for convenience:
//	RVector = TVector<double>	('real')
//	FVector = TVector<float>	('float')
//	CVector = TVector<complex>	('complex')
//	IVector = TVector<int>		('integer')
//	Vector  = TVector<double>	for backward compatibility
// ==========================================================================

#ifndef __VECTOR_H
#define __VECTOR_H

#include <string.h>
#ifndef __APPLE__
#include <malloc.h>
#endif
#include <iostream>
#include <sstream>
#include "mathdef.h"
#include "error.h"
//#include "complex.h"
//#include "scomplex.h"
#include <complex>
#include "fblas.h"

#ifdef USE_CBLAS
#include <cblas++.h>
#endif // USE_CBLAS

const int END = -1; // "end" index flag

// ==========================================================================
// Nonmember declarations

template<class VT> class TVector;

template<class VT>
TVector<VT> inv (const TVector<VT> &v);  // 1/v

template<class VT>
TVector<VT> sqr (const TVector<VT> &v);

template<class VT>
TVector<VT> sqrt (const TVector<VT> &v);

template<class VT>
TVector<VT> log (const TVector<VT> &v);

template<class VT>
TVector<VT> exp (const TVector<VT> &v);

template<class VT>
TVector<VT> pow (const TVector<VT> &v, const VT &s);

template<class VT>
TVector<VT> conj (const TVector<VT> &v);

template<class VT>
TVector<VT> sin (const TVector<VT> &v);

template<class VT>
TVector<VT> cos (const TVector<VT> &v);

template<class VT>
TVector<VT> tan (const TVector<VT> &v);

template<class VT>
TVector<VT> asin (const TVector<VT> &v);

template<class VT>
TVector<VT> acos (const TVector<VT> &v);

template<class VT>
TVector<VT> atan (const TVector<VT> &v);

template<class VT>
TVector<VT> sinh (const TVector<VT> &v);

template<class VT>
TVector<VT> cosh (const TVector<VT> &v);

template<class VT>
TVector<VT> tanh (const TVector<VT> &v);

template<class VT>
TVector<VT> asinh (const TVector<VT> &v);

template<class VT>
TVector<VT> acosh (const TVector<VT> &v);

template<class VT>
TVector<VT> atanh (const TVector<VT> &v);

template<class VT>
VT dot (const TVector<VT> &v1, const TVector<VT> &v2);

template<class VT>
VT doth (const TVector<VT> &v1, const TVector<VT> &v2);

template<class VT>
TVector<VT> cross (const TVector<VT> &v1, const TVector<VT> &v2);

template<class VT>
VT vmin (const TVector<VT> &v);

template<class VT>
VT vmax (const TVector<VT> &v);

template<class VT>
TVector<VT> vsort (const TVector<VT> &v);

template<class VT>
void vsort (const TVector<VT> &v, TVector<VT> &sorted_v, TVector<int> &sort_order);

template<class VT>
VT sum (const TVector<VT> &v);

template<class VT>
VT prod (const TVector<VT> &v);

template<class VT>
VT mean (const TVector<VT> &v);

template<class VT>
VT median (const TVector<VT> &v);

template<class VT>
VT variance (const TVector<VT> &v);

template<class VT>
VT stdv (const TVector<VT> &v);

template<class VT>
double l1norm (const TVector<VT> &v);

template<class VT>
double l2norm (const TVector<VT> &v);

template<class VT>
double linfnorm (const TVector<VT> &v);

template<class VT>
double l2normsq (const TVector<VT> &v);

template<class VT>
TVector<VT> &append (TVector<VT> &v1, const TVector<VT> &v2);

template<class VT>
TVector<VT> cat (const TVector<VT> &v1, const TVector<VT> &v2);

template<class VT>
TVector<VT> operator+ (const VT &s, const TVector<VT> &v);

template<class VT>
TVector<VT> operator- (const VT &s, const TVector<VT> &v);

template<class VT>
TVector<VT> operator* (const VT &s, const TVector<VT> &v);

template<class VT>
TVector<VT> operator/ (const VT &s, const TVector<VT> &v);

template<class VT>
bool operator== (const TVector<VT> &v1, const TVector<VT> &v2);

template<class VT>
bool operator!= (const TVector<VT> &v1, const TVector<VT> &v2);

template<class VT>
bool visnan (const TVector<VT> &v);

template<class VT>
std::ostream &operator<< (std::ostream &os, const TVector<VT> &v);

template<class VT>
std::istream &operator>> (std::istream &is, TVector<VT> &v);

template<class VT>
VT SparseDotp (const TVector<VT> &v1, idxtype *idx1, int nidx1,
	       const TVector<VT> &v2, idxtype *idx2, int nidx2,
	       int from, int to);

template<class VT>
TVector<double> UnfoldComplex (const TVector<VT> &v);

/* add specific function for CVector */
MATHLIB TVector<float> Re(const TVector<std::complex<float> > &vec);
MATHLIB TVector<double> Re(const TVector<std::complex<double> > &vec);
MATHLIB TVector<float> Im(const TVector<std::complex<float> > &vec);
MATHLIB TVector<double> Im(const TVector<std::complex<double> > &vec);
MATHLIB TVector<double> Mod(const TVector<std::complex<double> > &vec);
MATHLIB TVector<double> LogMod(const TVector<std::complex<double> > &vec);
MATHLIB TVector<double> Arg (const TVector<std::complex<double> > &vec);
MATHLIB TVector<std::complex<double> > Conj(const TVector<std::complex<double> > &vec); // obsolete
MATHLIB void SelfConj(const TVector<std::complex<double> > &vec);
MATHLIB TVector<std::complex<double> > Hadamard(const TVector<std::complex<double> > &a, const TVector<std::complex<double> > &b);
MATHLIB void SetReal (TVector<std::complex<double> > &z, const TVector<double> &zre);
MATHLIB void SetImag (TVector<std::complex<double> > &z, const TVector<double> &zim);

MATHLIB TVector<std::complex<double> > MakeCVector (const TVector<std::complex<float> > &v);

// ==========================================================================
// class TVector

/**
 * \brief Templated vector class
 *
 * Vectors consist of a header block and a dynamically allocated data block.
 * Multiple vectors can share a common data block (or parts thereof). Any
 * changes to a common data block affect all associated vectors simultaneously.
 * 
 * The following template types are instantiated by default:
 * - TVector<double>  (RVector)
 * - TVector<float>   (FVector)
 * - TVector<complex> (CVector)
 * - TVector<int>     (IVector)
 */

template<class VT> class TVector {
public:
    /**
     * \brief Constructor. Creates a vector of length 0
     */
    TVector ();

    /**
     * \brief Constructor. Create a vector of length 'dim' and zero all elements.
     * \param dim vector size (>= 0)
     */
    TVector (int dim);

    /**
     * \brief Constructor. Creates a vector of length 'dim' and set all values to 's'
     * \param dim vector size (>= 0)
     * \param s element value
     */
    TVector (int dim, const VT s);

    /**
     * \brief Constructor. Create a vector of length 'dim' and initialise
     *   from 'values' array.
     * \param dim vector size (>= 0)
     * \param values array of element values
     * \param cmode if set to SHALLOW_COPY, the value array is not copied
     *   into a local buffer but used directly.
     * \note The 'values' array must contain at least 'dim' elements.
     * \note If cmode==SHALLOW_COPY,  the buffer must exist for the lifetime
     *   of the vector, and the caller is responsible for cleaning it
     *   up after use.
     */
    TVector (int dim, VT *values, CopyMode cmode=DEEP_COPY);

    /**
     * \brief Constructor. Create a vector of length 'dim' and initialise
     *   values from a text string.
     * \param dim vector size (>= 0)
     * \param init text string containing list of vector elements
     * \note \e init must contain a list of at least \e dim elements.
     * \note Elements must be in a format that is readable for the C++
     *   stream-in operator ('>>') for the appropriate template type.
     * \note Elements must be separated by whitespace.
     */
    TVector (int dim, const char *init);

    /**
     * \brief Copy constructor. Create vector as copy of another vector.
     * \param v original vector
     */
    TVector (const TVector<VT> &v);

    /**
     * \brief Reference constructor. Create a vector which shares its data
     *   block with vector v.
     * \param v original vector
     * \param ofs offset into the data block of the original vector that
     *   is to be used as the first element of the reference vector
     * \param dim size of the reference vector
     * \note The new vector has length 'dim' and shares all elements from
     *   ofs <= i < ofs+dim with vector v.
     * \note Any changes to the original vector in that range also affect the
     *   new vector (and vice versa).
     * \note ofs+dim <= v.Dim() is required.
     * \note This constructor is not threadsafe if used by multiple threads
     *   referencing the same vector v.
     */
    TVector (const TVector<VT> &v, int ofs, int dim);

    /**
     * \brief Destructor. Delete vector and deallocate data block.
     */
    ~TVector () { Unlink (); }

    /**
     * \brief Returns the size of the vector.
     * \return Size (number of elements) of the vector.
     */
    int Dim () const { return size; }

    /**
     * \brief Vector element access operator (read and write)
     * \param i element index (0 <= i < TVector::Dim())
     * \note The operator syntax is: x = v[i]
     */
    VT &operator[] (int i) const;

    /// \brief Assignment operator.
    ///
    /// - Sets the left-hand operator to 'v' and resizes, if required.
    /// - The return value of the operation is the copied vector, so that
    ///   constructs like v1 = v2 = v0 are possible.
    TVector &operator= (const TVector &v) { Copy (v); return *this; }

    /// Scalar assignment. All vector elements are set to 's'.
    TVector &operator= (VT s);

    /// \brief Vector copy. Replaces the vector with a copy of 'v'
    ///
    /// - This method calls the xCOPY BLAS functions for float and double
    ///   types, if USE_BLAS_LEVEL1 is defined in toastdef.h
    void Copy (const TVector &v);

    /// \brief Partial vector copy. Replaces the vector with a part of 'v'
    ///
    /// - A block of data of length n is copied from vector 'v' (starting
    ///   at index 'sofs', into the data block of *this, starting at index
    ///   'tofs': (*this)[i+tofs]=v[i+sofs], i=0...Dim()-1.
    /// - The size of the target vector (*this) is not changed.
    /// - If the data chunk exceeds the limits of the source or target vector,
    ///   then n is truncated as required.
    /// - If n < 0 then the chunk size is maximised.
    /// - This method calls the xCOPY BLAS functions for float and double
    ///   types, if USE_BLAS_LEVEL1 is defined.
    void Copy (const TVector &v, int tofs, int sofs, int n);

    /// \brief Resize the vector.
    /// \param dim New vector length
    /// \remarks Re-allocates the data block to size 'dim' and resets all
    ///   values to zero.
    /// \sa Allocate(), Unlink()
    void New (int dim) { Unlink (); Allocate (dim); }

    /**
     * \brief Zeroes all elements.
     * \remarks Resets all elements to (VT)0, but does not resize the vector.
     * \remarks Equivalent to, but more efficient than *this = 0
     */
    void Clear ();

    /// \brief Link the vector to the data block of another vector
    ///
    /// - Shares data with another vector.
    /// - This method unlinks the vector from its current data block and links
    ///   it to the data block of 'v'. The original data block is removed from
    ///   memory, unless it is referenced by another vector.
    void Relink (const TVector &v);

    /// \brief Link the vector to part of the data block of another vector
    ///
    /// - Shares data with another vector.
    /// - This method unlinks the vector from its current data block and links
    ///   it to a subset of the data block of 'v' of length 'n', starting at
    ///   index ofs.
    /// - Requires
    ///   - ofs >= 0
    ///   - n > 0
    ///   - ofs+n <= v.Dim()
    void Relink (const TVector &v, int ofs, int n);

    /**
     * \brief Link vector to an external data array.
     * \param values external array of length >= dim
     * \param dim vector size
     * \note The external buffer is not copied into a local buffer but used
     *   directly. It must exist during the lifetime of the vector, or until
     *   it is linked to a different buffer.
     */
    void Relink (VT *values, int dim);

    /// \brief Direct access to vector data array.
    ///
    /// - Returns a pointer to the vector's data block to provide
    ///   low-level access (read and write).
    /// - The pointer may become invalid if the vector re-allocates its
    ///   buffer, e.g. after a call to New().
    VT *data_buffer () { return data; }

    /// \brief Direct access to vector data array.
    ///
    /// - Returns a pointer to the vector's data block to provide
    ///   low-level access (read only).
    /// - The pointer may become invalid if the vector re-allocates its
    ///   buffer, e.g. after a call to New().
    const VT *data_buffer () const { return data; }

    /// \brief Truncate vector elements to specified range.
    ///
    /// - Truncate all elements of the vector to the range vmin ... vmax.
    /// - Return value is true if any elements have been truncated.
    /// - This method is only defined for data types which support boolean
    ///   comparisons.
    bool Clip (VT vmin, VT vmax);

    /// \brief Move vector elements left.
    ///
    /// - Shifts all elements of the vector n places left (toward lower
    ///   indices)
    /// - The n leftmost elements are dropped.
    /// - The n rightmost elements are reset to zero.
    void ShiftLeft (int n);

    /// \brief Move vector elements right.
    ///
    /// - Shifts all elements of the vector n places right (toward higher
    ///   indices)
    /// - The n rightmost elements are dropped.
    /// - The n leftmost elements are reset to zero.
    void ShiftRight (int n);

    // arithmetic operators

    /// \brief Element-wise addition of vectors
    ///
    /// - Return value is *this[i] + v[i]
    /// - The operands are unchanged.
    /// - Operands must have the same length
    /// - This method calls the xAXPY BLAS functions for float and double
    ///   types, if USE_BLAS_LEVEL1 is defined.
    TVector operator+ (const TVector &v) const;

    /// \brief Element-wise subtraction of vectors
    ///
    /// - Return value is *this[i] - v[i]
    /// - The operands are unchanged.
    /// - Operands must have the same length
    /// - This method calls the xAXPY BLAS functions for float and double
    ///   types, if USE_BLAS_LEVEL1 is defined.
    TVector operator- (const TVector &v) const;

    /// \brief Element-wise multiplication of vectors
    ///
    /// - Return value is *this[i] * v[i]
    /// - The operands are unchanged.
    /// - Operands must have the same length
    TVector operator* (const TVector &v) const;

    /// \brief Element-wise division of vectors
    ///
    /// - Return value is *this[i] / v[i]
    /// - The operands are unchanged.
    /// - Operands must have the same length
    TVector operator/ (const TVector &v) const;

    /// \brief Element-wise addition with a scalar
    ///
    /// - Return value is *this[i] + s
    /// - The operands are unchanged.
    TVector operator+ (const VT &s) const;

    /// \brief Element-wise subtraction of a scalar
    ///
    /// - Return value is *this[i] - s
    /// - The operands are unchanged.
    TVector operator- (const VT &s) const;

    /// \brief Element-wise multiplication with a scalar
    ///
    /// - Return value is *this[i] * s
    /// - The operands are unchanged.
    TVector operator* (const VT &s) const;

    /// \brief Element-wise division by a scalar
    ///
    /// - Return value is *this[i] / s
    /// - The operands are unchanged.
    TVector operator/ (const VT &s) const;

    /// \brief Unary minus.
    ///
    /// - Return value is -*this[i]
    /// - The operand is unchanged.
    TVector operator- () const;

    /**
     * \brief Vector comparison (relational operator)
     * \param v1 left-hand vector operand
     * \param v2 right-hand vector operand
     * \return Returns \e true if both vectors are identical, \e false
     *   otherwise.
     */
    friend bool operator== <> (const TVector<VT> &v1,
        const TVector<VT> &v2);

    /**
     * \brief Vector comparison (relational operator)
     * \param v1 left-hand vector operand
     * \param v2 right-hand vector operand
     * \return Returns \e false if both vectors are identical, \e true
     *   otherwise.
     */
    friend bool operator!=<> (const TVector<VT> &v1,
        const TVector<VT> &v2);

    /// \brief Dot product of two vectors
    /// \param v Right-hand operand
    /// \return dot product
    ///   \f$\vec{*this}\cdot\vec{v} = \sum_i *this_i \cdot v_i\f$
    /// \par Conditions:
    ///   this->Dim() == v.Dim()
    /// \remarks The operands are unchanged.
    /// \remarks If USE_BLAS_LEVEL1 is set in toastdef.h, this calls BLAS
    ///   functions ddot (for TVector<double>), sdot (for TVector<float>),
    ///   or zdotu (for TVector<complex>)
    /// \sa dot()
    VT operator& (const TVector& v) const { return dot (*this, v); }

    // **** arithmetic/assignment operators ****

    /// \brief Element-wise addition/assignment of vectors
    ///
    /// (*this)[i] += v[i]
    /// \param v Right-hand operand (vector)
    /// \return The result of the addition (the new value of *this)
    TVector &operator+= (const TVector &v);

    /// \brief Element-wise subtraction/assignment of vectors
    ///
    /// (*this)[i] -= v[i]
    /// \param v Right-hand operand (vector)
    /// \return The result of the subtraction (the new value of *this)
    TVector &operator-= (const TVector &v);

    /// \brief Element-wise multiplication/assignment of vectors
    ///
    /// (*this)[i] *= v[i]
    /// \param v Right-hand operand (vector)
    /// \return The result of the multiplication (the new value of *this)
    TVector &operator*= (const TVector &v);

    /// \brief Element-wise division/assignment of vectors
    ///
    /// (*this)[i] /= v[i]
    /// \param v Right-hand operand (vector)
    /// \return The result of the division (the new value of *this)
    TVector &operator/= (const TVector &v);

    /// \brief Element-wise addition/assignment with a scalar
    ///
    /// (*this)[i] += s
    /// \param s Right-hand operand (scalar)
    /// \return The result of the addition (the new value of *this)
    TVector &operator+= (const VT &s);

    /// \brief Element-wise subtraction/assignment of a scalar
    ///
    /// (*this)[i] -= s
    /// \param s Right-hand operand (scalar)
    /// \return The result of the subtraction (the new value of *this)
    TVector &operator-= (const VT &s);

    /// \brief Element-wise multiplication/assignment with a scalar
    ///
    /// (*this)[i] *= s
    /// \param s Right-hand operand (scalar)
    /// \return The result of the multiplication (the new value of *this)
    TVector &operator*= (const VT &s);

    /// \brief Element-wise division/assignment of a scalar
    ///
    /// (*this)[i] /= s
    /// \param s Right-hand operand (scalar)
    /// \return The result of the division (the new value of *this)
    TVector &operator/= (const VT &s);

    // IO

    /// \brief Read vector from stream
    /// \param is input stream
    /// \param start index of first element to read
    /// \param n max. number of elements to read (-1 for all)
    /// \remarks This method allows to read a partial vector from a stream,
    ///   by skipping the first elements (start > 0) and dropping the last
    ///   elements (n < stored vector length)
    /// \sa operator>>()
    void Read (std::istream &is, int start = 0, int n = -1);

    /// \brief Read indexed vector elements from stream
    /// \param is input stream
    /// \param n vector length
    /// \remarks Resizes *this to length n and reads index-value
    ///   pairs (white-space separated) from stream to fill it.
    /// \remarks Indices are zero-based.
    /// \remarks Unreferenced entries are set to 0.
    /// \remarks The first value read from the stream is interpreted as the
    ///   number of index-value pairs to follow.
    /// \remarks If n < 0, the vector is assumed to be already of correct
    ///   length (no resizing)
    void ReadIndexed (std::istream &is, int n = -1);

    // friends
    /// \brief Returns element-wise inverse vector 1/(*this)[i]
    friend TVector<VT> inv<> (const TVector<VT> &v);

    /// \brief Returns element-wise square vector (*this)[i]^2
    friend TVector<VT> sqr<> (const TVector<VT> &v);

    /// \brief Returns element-wise square-root vector (*this)[i]^(1/2)
    friend TVector<VT> sqrt<> (const TVector<VT> &v);

    /// \brief Returns element-wise power vector (*this)[i]^s
    friend TVector<VT> pow<> (const TVector<VT> &v, const VT &s);

    /// \brief Returns element-wise natural logarithm vector ln((*this)[i])
    friend TVector<VT> log<> (const TVector<VT> &v);

    /// \brief Returns element-wise base-e exponential vector exp((*this)[i])
    friend TVector<VT> exp<> (const TVector<VT> &v);

    /// \brief Returns element-wise complex conjugate vector ((*this([i].re,
    ///   -(*this)[i].im)
    friend TVector<VT> conj<> (const TVector<VT> &v);

    /// \brief Returns element-wise sine vector sin((*this)[i])
    friend TVector<VT> sin<> (const TVector<VT> &v);

    /// \brief Returns element-wise cosine vector cos((*this)[i])
    friend TVector<VT> cos<> (const TVector<VT> &v);

    /// \brief Returns element-wise tangent vector tan((*this)[i])
    friend TVector<VT> tan<> (const TVector<VT> &v);

    /// \brief Returns element-wise arc sine vector asin((*this)[i])
    friend TVector<VT> asin<> (const TVector<VT> &v);

    /// \brief Returns element-wise arc cosine vector acos((*this)[i])
    friend TVector<VT> acos<> (const TVector<VT> &v);

    /// \brief Returns element-wise arc tangent vector atan((*this)[i])
    friend TVector<VT> atan<> (const TVector<VT> &v);

    /// \brief Returns element-wise hyperbolic sine vector sinh((*this)[i])
    friend TVector<VT> sinh<> (const TVector<VT> &v);

    /// \brief Returns element-wise hyperbolic cosine vector cosh((*this)[i])
    friend TVector<VT> cosh<> (const TVector<VT> &v);

    /// \brief Returns element-wise hyperbolic tangent vector tanh((*this)[i])
    friend TVector<VT> tanh<> (const TVector<VT> &v);

    /// \brief Returns element-wise inverse hyperbolic sine vector
    ///   asinh((*this)[i])
    friend TVector<VT> asinh<> (const TVector<VT> &v);

    /// \brief Returns element-wise inverse hyperbolic cosine vector
    ///   acosh((*this)[i])
    friend TVector<VT> acosh<> (const TVector<VT> &v);

    /**
     * \brief Returns element-wise inverse hyperbolic tangent vector
     *   atanh(v[i])
     * \param v vector argument
     * \return result vector
     * \par Conditions:
     *   |v[i]| <= 1 is required for all elements of a real input vector
     * \sa sinh, cosh, tanh, asinh, acosh
     */
    friend TVector<VT> atanh<> (const TVector<VT> &v);    

    /**
     * \brief Dot product of two vectors
     * \param v1 first vector argument
     * \param v2 second vector argument
     * \return dot product \f$ \vec{v1}\cdot\vec{v2}=\sum_i v1_i \cdot v2_i\f$
     * \par Conditions:
     *   v1.Dim() == v2.Dim() (checked in DEBUG only)
     * \remarks If USE_BLAS_LEVEL1 is set in toastdef.h, this calls BLAS
     *   functions ddot (for TVector<double>), sdot (for TVector<float>),
     *   or zdotu (for TVector<complex>)
     * \sa operator&().
     */
    friend VT dot<> (const TVector<VT> &v1, const TVector<VT> &v2);

    /**
     * \brief Hermitian product of two vectors
     * \return Hermitian product \f$ \sum_i v1_i^* \cdot v2_i\f$
     * \par Conditions:
     *   v1.Dim() == v2.Dim() (checked in DEBUG only)
     * \remarks For all template types other than complex, this function is
     *   identical to dot(). For complex, this returns the inner product
     *   of the complex conjugate of the first argument, and the second
     *   argument.
     * \remarks If USE_BLAS_LEVEL1 is set in toastdef.h, the function calls
     *   BLAS function zdotc.
     */
    friend VT doth<> (const TVector<VT> &v1, const TVector<VT> &v2);

    /**
     * \brief Cross product of two vectors
     * \param v1 first vector argument
     * \param v2 second vector argument
     * \return Result of v1 x v2
     * \note v1 and v2 must be of size 3.
     */
    friend TVector<VT> cross<> (const TVector<VT> &v1,
        const TVector<VT> &v2);
    
    /**
     * \brief Extract smallest element in vector.
     * \return \f$ \min_i v_i \f$
     * \remarks Only works with template types that define relational
     *   operator <.
     */
    friend VT vmin<> (const TVector<VT> &v);

    /**
     * \brief Extract largest element in vector.
     * \return \f$ \max_i v_i \f$
     * \remarks Only works with template types that define relational
     *   operator >.
     */
    friend VT vmax<> (const TVector<VT> &v);

    /**
     * \brief Sort vector
     * \return Vector with elements of v sorted in ascending order.
     * \remarks Uses a simple bubble sort algorithm
     * \remarks Input argument unchanged.
     * \remarks Only works with template types that define relational
     *   operator <.
     */
    friend TVector<VT> vsort<> (const TVector<VT> &v);

    /**
     * \brief Sort vector
     * \return Vector with elements of v sorted in ascending order.
     * \return Sort order vector.
     * \remarks Uses a simple bubble sort algorithm
     * \remarks Input argument unchanged.
     * \remarks Only works with template types that define relational
     *   operator <.
     */
	
    friend void vsort<> (const TVector<VT> &v, TVector<VT> &sorted_v, TVector<int> &sort_order);

    /**
     * \brief Sum of elements
     * \return \f$ \sum_i^{v.Dim()} v_i \f$
     */
    friend VT sum<> (const TVector<VT> &v);

    /**
     * \brief Product of elements
     * \return \f$ \Prod_i^{v.Dim()} v_i \$f
     */
    friend VT prod<> (const TVector<VT> &v);

    /**
     * \brief Mean value of elements
     * \return \f$ E(v) = (1/v.Dim()) \sum_i^{v.Dim()} v_i \f$
     * \par Condition
     *   v.Dim >= 1
     * \sa median(), variance(), stdv()
     */
    friend VT mean<> (const TVector<VT> &v);

    /**
     * \brief Median value of elements
     *
     * Returns a value from v such that at most half the values of v have
     * values less than the median and at most half the have values greater
     * than the median.
     * \remarks If v.Dim() is even, returns the mean of the two middle values.
     * \remarks If v.Dim() == 0, returns 0.
     * \sa mean()
     */
    friend VT median<> (const TVector<VT> &v);

    /**
     * \brief Variance of element values
     * \return \f$ var(v) = E((v-E(v))^2) = E(v^2) - (E(v))^2 \f$
     * \sa mean(), stdv()
     */
    friend VT variance<> (const TVector<VT> &v);

    /**
     * \brief Standard deviation of element values
     * \return \f$ \sigma(v) = \sqrt{var(v)} \f$
     * \sa mean(), variance()
     */
    friend VT stdv<> (const TVector<VT> &v);

    /**
     * \brief L1 norm of vector elements
     * \return \f$ \sum_i^{v.Dim()} | v_i | \f$
     * \sa l2norm(), linfnorm()
     */
    friend double l1norm<> (const TVector<VT> &v);

    /**
     * \brief L2 norm of vector elements
     * \return \f$ \sqrt{\sum_i^{v.Dim()} v_i^2} \f$
     * \sa length(), l2normsq(), l1norm(), linfnorm()
     */
    friend double l2norm<> (const TVector<VT> &v);

    /**
     * \brief L-infinity norm of vector elements
     * \return \f$ \max_i | v_i | \f$
     * \sa l1norm(), l2norm()
     */
    friend double linfnorm<> (const TVector<VT> &v);

    /**
     * \brief Square of L2 norm of vector elements
     * \return \f$ \sum_i^{v.Dim()} v_i^2 \f$
     * \sa l2norm()
     */
    friend double l2normsq<> (const TVector<VT> &v);

    /**
     * \brief Length of a vector (L2 norm of elements)
     * \return \f$ \sqrt{\sum_i^{v.Dim()} v_i^2} \f$
     * \remarks Synonymous to l2norm().
     */
    friend double length (const TVector<VT> &vec)
	{ return l2norm (vec); }

    /**
     * \brief Concatenate two vectors
     *
     * Appends the elements of the second argument (v2) to the end of the
     * first argument (v1).
     * \return Result of the concatenation (the modified v1 argument)
     * \remarks Argument v2 unchanged.
     */
    friend TVector<VT> &append<> (TVector<VT> &v1, const TVector<VT> &v2);

    /**
     * \brief Concatenate two vectors
     *
     * Creates a new vector by concatenating the two arguments.
     * \return The result of the concatenation, {v1,v2}
     * \note Both arguments unchanged
     */
    friend TVector<VT> cat<> (const TVector<VT> &v1, const TVector<VT> &v2);

    /**
     * \brief Element-wise addition of a scalar from the left
     *
     * Allows adding scalars from the left, in constructs like
     * \code
     * w = s + v;
     * \endcode
     * \param s scalar value
     * \param v vector value
     * \return Result of the addition, {s+v[i]}
     * \sa TVector::operator+()
     */
#ifndef USE_CUDA_FLOAT
// none of the operator friend definitions seem to work with nvcc
#if defined(_WIN32)
    template <class VT> friend TVector<VT> operator+ (const VT &s,
        const TVector<VT> &v);
#elif (GCC_VERSION < 30404) // old-style definitions
    friend TVector<VT> operator+ <> (const VT &s, const TVector<VT> &v);
#else
    friend TVector<VT> (::operator+ <>) (const VT &s, const TVector<VT> &v);
#endif
#endif // !USE_CUDA_FLOAT

    /**
     * \brief Element-wise subtraction from a scalar from the left
     *
     * Allows subtraction from a scalar on the left, in constructs like
     * \code
     * w = s - v;
     * \endcode
     * \param s scalar value
     * \param v vector value
     * \return Result of the subtraction, {s-v[i]}
     * \sa TVector::operator-()
     */
#if defined(_WIN32)
    template<class VT> friend TVector<VT> operator- (const VT &s,
        const TVector<VT> &v);
#elif (GCC_VERSION < 30404) // old-style definitions
    friend TVector<VT> operator- <> (const VT &s, const TVector<VT> &v);
#else
    friend TVector<VT> (::operator- <>) (const VT &s, const TVector<VT> &v);
#endif

    /**
     * \brief Element-wise multiplication with a scalar from the left
     *
     * Allows multiplication with a scalar from the left, in constructs like
     * \code
     *   w = s * v;
     * \endcode
     * \param s scalar value
     * \param v vector value
     * \return Result of the multiplication, {s*v[i]}
     * \sa TVector::operator*()
     */
#if defined(_WIN32)
    template<class VT> friend TVector<VT> operator* (const VT &s,
        const TVector<VT> &v);
#elif (GCC_VERSION < 30404) // old-style definitions
    friend TVector<VT> operator* <> (const VT &s, const TVector<VT> &v);
#else
    friend TVector<VT> (::operator* <>) (const VT &s, const TVector<VT> &v);
#endif

    /**
     * \brief Element-wise division of a scalar by a vector
     *
     * Division of a scalar by a vector, in constructs like
     * \code
     *   w = s / v;
     * \endcode
     * \param s scalar value
     * \param v vector value
     * \return Result of the division, {s/v[i]}
     * \sa TVector::operator/()
     */
#if defined(_WIN32)
    template<class VT> friend TVector<VT> operator/ (const VT &s,
        const TVector<VT> &v);
    // JK also works (TVector::operator/ )
#elif (GCC_VERSION < 30404) // old-style definitions
    friend TVector<VT> operator/ <> (const VT &s, const TVector<VT> &v);
#else
    friend TVector<VT> (::operator/ <>) (const VT &s, const TVector<VT> &v);
#endif

    friend bool visnan<> (const TVector<VT> &v);

    /**
     * \brief Write vector to output stream
     * \param os output stream
     * \param v vector to be written
     * \return output stream (to allow chaining of output objects)
     * \note The output is whitespace-separated ASCII text, enclosed by
     *   square brackets. Example:
     * \code
     * [0 7.12 -3.1415]
     * \endcode
     * \note For complex vectors, each complex entry follows the standard
     *   toast notation of separating the real and imaginary part with
     *   whitespace and enclosing in a <> bracket pair. Example:
     * \code
     * [<0 0> <-3.1 44.123>]
     * \endcode
     */
    friend std::ostream &operator<< <> (std::ostream &os,
        const TVector<VT> &v);

    /**
     * \brief Read vector from input stream
     * \param is input stream
     * \param v vector to be read
     * \return input stream (to allow chaining of input objects)
     * \note The vector is resized to the appropriate length as required.
     * \note The expected format of the input stream is whitespace-separated
     *   ASCII text, enclosed by square brackets. Example:
     * \code
     * [0 7.12 -3.1415]
     * \endcode
     * \note For complex vectors, each complex entry follows the standard
     *   toast notation of separating the real and imaginary part with
     *   whitespace and enclosing in a <> bracket pair. Example:
     * \code
     * [<0 0> <-3.1 44.123>]
     * \endcode
     */
    friend std::istream &operator>> <> (std::istream &is, TVector<VT> &v);

    /**
     * \brief Sparse dot product of two vectors
     *
     * Given vectors v1 and v2 and index lists idx1, idx2 of length
     * nidx1 and nidx2, respectively, calculate the dot product using
     * only the vector indices from the lists.
     * \param v1 first vector argument
     * \param idx1 first index array
     * \param nidx1 length of first index array
     * \param v2 second vector argument
     * \param idx2 second index array
     * \param nidx2 length of second index array
     * \param from lower index limit
     * \param to high index limit
     * \return sparse dot product
     * \note idx1 and idx2 are assumed sorted (ascending)
     */
    friend VT SparseDotp<> (const TVector<VT> &v1, idxtype *idx1, int nidx1,
			    const TVector<VT> &v2, idxtype *idx2, int nidx2,
			    int from, int to);

    /**
     * \brief Unfold real and imaginary parts of a complex vector.
     * \param v input vector argument
     * \return Unfolded real vector
     * \note For complex template types, this method creates a real output
     *   vector of twice the size of the input vector, containing the real
     *   and imaginary parts of the vector.
     * \note For real template types, this simply returns the input vector.
     */
    friend TVector<double> UnfoldComplex<> (const TVector<VT> &v);

    /**
     * \brief Initialise vector elements from a string.
     * \param cbuf ASCII string of values (white-space separated)
     * \param nmax Max number of values to read (or 0 for all)
     * \note The string should NOT contain vector delimiters [].
     * \note Vector must be allocated beforehand (no resizing)
     */
    void Scan (const char *cbuf, int nmax = 0);

protected:
    /// \brief Allocate data buffer
    /// \param dim length of new buffer
    /// \remarks Resets the length of the vector to dim and allocates a data
    ///   buffer to store the element values.
    /// \remarks Values are initialised to zero.
    /// \remarks If the vector already had an associated data block, this
    ///   must be detached with Unlink first.
    /// \sa Unlink(), New()
    void Allocate (int dim);

    /// \brief Link the vector to the data block of another vector.
    /// \param vec referenced vector
    /// \remarks The two vectors share the same data block. Any changes to
    ///   element values apply to both simultaneously.
    /// \remarks A re-allocation of the data block of one of the vectors
    ///   (e.g. as a result of resizing) disconnects the two vectors.
    /// \remarks If the vector already had an associated data block, this
    ///   must be detached with Unlink first.
    /// \sa Unlink(), Relink()
    void Link (const TVector<VT> &vec);

    /// \brief Link the vector to part of the data block of another vector.
    /// \param vec referenced vector
    /// \param ofs starting index of the referenced block of elements
    ///   (zero-based)
    /// \param dim block length
    /// \remarks The two vectors share part of their data block. Any changes
    ///   to element values apply to boths simultaneously.
    /// \remarks A re-allocation of the data block of one of the vectors
    ///   (e.g. as a result of resizing) disconnects the two vectors.
    /// \remarks If the vector already had an associated data block, this
    ///   must be detached with Unlink first.
    /// \sa Unlink(), Relink()
    void Link (const TVector<VT> &vec, int ofs, int dim);

    /**
     * \brief Link vector to an external data buffer
     * \param values external data array of length >= dim
     * \param dim vector size
     * \note The external buffer must be guaranteed to exist while linked
     *   by the vector.
     */
    void Link (VT *values, int dim);

    /// \brief Unlink vector from its data block.
    /// \remarks This results in a (valid) vector of length zero.
    /// \sa New(), Clear()
    void Unlink ();

    int size;		 ///< vector length
    int *base_size;	 ///< pointer to length of linked data block
    int *base_nref;	 ///< pointer to data block reference counter
    VT *data;		 ///< pointer to first vector element
    VT *base_data;	 ///< pointer to linked data block

private:
    bool ext_data;       ///< data buffer is external
};

// ==========================================================================
// template typedefs

typedef TVector<double>          RVector;
typedef TVector<float>           FVector;
typedef TVector<std::complex<double> >  CVector;
typedef TVector<std::complex<float> >   SCVector;
typedef TVector<int>             IVector;

#endif // !__VECTOR_H
