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
#include <malloc.h>
#include <iostream>
#include <sstream>
#include "mathdef.h"
#include "error.h"
#include "complex.h"
#include "scomplex.h"
#include "fblas.h"

#ifdef USE_CBLAS
#include <cblas++.h>
#endif // USE_CBLAS

const int END = -1; // "end" index flag

// ==========================================================================
// Nonmember declarations

template<class VT> class TVector;

template<class VT>
MATHLIB TVector<VT> inv (const TVector<VT> &v);  // 1/v

template<class VT>
MATHLIB TVector<VT> sqr (const TVector<VT> &v);

template<class VT>
MATHLIB TVector<VT> sqrt (const TVector<VT> &v);

template<class VT>
MATHLIB TVector<VT> log (const TVector<VT> &v);

template<class VT>
MATHLIB TVector<VT> exp (const TVector<VT> &v);

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
MATHLIB VT vmin (const TVector<VT> &v);

template<class VT>
MATHLIB VT vmax (const TVector<VT> &v);

template<class VT>
TVector<VT> vsort (const TVector<VT> &v);

template<class VT>
void vsort (const TVector<VT> &v, TVector<VT> &sorted_v, TVector<int> &sort_order);

template<class VT>
VT sum (const TVector<VT> &v);

template<class VT>
MATHLIB VT mean (const TVector<VT> &v);

template<class VT>
VT median (const TVector<VT> &v);

template<class VT>
VT variance (const TVector<VT> &v);

template<class VT>
VT stdv (const TVector<VT> &v);

template<class VT>
MATHLIB double l1norm (const TVector<VT> &v);

template<class VT>
MATHLIB double l2norm (const TVector<VT> &v);

template<class VT>
MATHLIB double linfnorm (const TVector<VT> &v);

template<class VT>
MATHLIB double l2normsq (const TVector<VT> &v);

template<class VT>
MATHLIB TVector<VT> &append (TVector<VT> &v1, const TVector<VT> &v2);

template<class VT>
MATHLIB TVector<VT> cat (const TVector<VT> &v1, const TVector<VT> &v2);

template<class VT>
MATHLIB TVector<VT> operator+ (const VT &s, const TVector<VT> &v);

template<class VT>
MATHLIB TVector<VT> operator- (const VT &s, const TVector<VT> &v);

template<class VT>
MATHLIB TVector<VT> operator* (const VT &s, const TVector<VT> &v);

template<class VT>
MATHLIB TVector<VT> operator/ (const VT &s, const TVector<VT> &v);

template<class VT>
MATHLIB bool operator== (const TVector<VT> &v1, const TVector<VT> &v2);

template<class VT>
MATHLIB bool operator!= (const TVector<VT> &v1, const TVector<VT> &v2);

template<class VT>
MATHLIB bool visnan (const TVector<VT> &v);

template<class VT>
MATHLIB std::ostream &operator<< (std::ostream &os, const TVector<VT> &v);

template<class VT>
MATHLIB std::istream &operator>> (std::istream &is, TVector<VT> &v);

template<class VT>
VT SparseDotp (const TVector<VT> &v1, idxtype *idx1, int nidx1,
	       const TVector<VT> &v2, idxtype *idx2, int nidx2,
	       int from, int to);

template<class VT>
MATHLIB TVector<double> UnfoldComplex (const TVector<VT> &v);

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
    inline TVector ();

    /**
     * \brief Constructor. Create a vector of length 'dim' and zero all elements.
     * \param dim vector size (>= 0)
     */
    inline TVector (int dim);

    /**
     * \brief Constructor. Creates a vector of length 'dim' and set all values to 's'
     * \param dim vector size (>= 0)
     * \param s element value
     */
    inline TVector (int dim, const VT s);

    /**
     * \brief Constructor. Create a vector of length 'dim' and initialise
     *   from 'values' array.
     * \param dim vector size (>= 0)
     * \param values array of element values
     * \note The 'values' array must contain at least 'dim' elements.
     */
    inline TVector (int dim, VT *values);

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
    inline TVector (int dim, const char *init);

    /**
     * \brief Copy constructor. Create vector as copy of another vector.
     * \param v original vector
     */
    inline TVector (const TVector<VT> &v);

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
     */
    inline TVector (const TVector<VT> &v, int ofs, int dim);

    /**
     * \brief Destructor. Delete vector and deallocate data block.
     */
    inline ~TVector () { Unlink (); }

    /**
     * \brief Returns the size of the vector.
     * \return Size (number of elements) of the vector.
     */
    inline int Dim () const { return size; }

    /**
     * \brief Vector element access operator (read and write)
     * \param i element index (0 <= i < TVector::Dim())
     * \note The operator syntax is: x = v[i]
     */
    inline VT &operator[] (int i) const;

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

    //inline VT Vmin () const { return vmin (*this); }
    // returns smallest element

    //inline VT Vmax () const { return vmax (*this); }
    // returns largest element

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
    friend MATHLIB bool operator== <> (const TVector<VT> &v1,
        const TVector<VT> &v2);

    /**
     * \brief Vector comparison (relational operator)
     * \param v1 left-hand vector operand
     * \param v2 right-hand vector operand
     * \return Returns \e false if both vectors are identical, \e true
     *   otherwise.
     */
    friend MATHLIB bool operator!=<> (const TVector<VT> &v1,
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
    friend MATHLIB TVector<VT> inv<> (const TVector<VT> &v);

    /// \brief Returns element-wise square vector (*this)[i]^2
    friend MATHLIB TVector<VT> sqr<> (const TVector<VT> &v);

    /// \brief Returns element-wise square-root vector (*this)[i]^(1/2)
    friend MATHLIB TVector<VT> sqrt<> (const TVector<VT> &v);

    /// \brief Returns element-wise power vector (*this)[i]^s
    friend TVector<VT> pow<> (const TVector<VT> &v, const VT &s);

    /// \brief Returns element-wise natural logarithm vector ln((*this)[i])
    friend MATHLIB TVector<VT> log<> (const TVector<VT> &v);

    /// \brief Returns element-wise base-e exponential vector exp((*this)[i])
    friend MATHLIB TVector<VT> exp<> (const TVector<VT> &v);

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
     * \brief Extract smallest element in vector.
     * \return \f$ \min_i v_i \f$
     * \remarks Only works with template types that define relational
     *   operator <.
     */
    friend MATHLIB VT vmin<> (const TVector<VT> &v);

    /**
     * \brief Extract largest element in vector.
     * \return \f$ \max_i v_i \f$
     * \remarks Only works with template types that define relational
     *   operator >.
     */
    friend MATHLIB VT vmax<> (const TVector<VT> &v);

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
     * \brief Mean value of elements
     * \return \f$ E(v) = (1/v.Dim()) \sum_i^{v.Dim()} v_i \f$
     * \par Condition
     *   v.Dim >= 1
     * \sa median(), variance(), stdv()
     */
    friend MATHLIB VT mean<> (const TVector<VT> &v);

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
    friend MATHLIB double l1norm<> (const TVector<VT> &v);

    /**
     * \brief L2 norm of vector elements
     * \return \f$ \sqrt{\sum_i^{v.Dim()} v_i^2} \f$
     * \sa length(), l2normsq(), l1norm(), linfnorm()
     */
    friend MATHLIB double l2norm<> (const TVector<VT> &v);

    /**
     * \brief L-infinity norm of vector elements
     * \return \f$ \max_i | v_i | \f$
     * \sa l1norm(), l2norm()
     */
    friend MATHLIB double linfnorm<> (const TVector<VT> &v);

    /**
     * \brief Square of L2 norm of vector elements
     * \return \f$ \sum_i^{v.Dim()} v_i^2 \f$
     * \sa l2norm()
     */
    friend MATHLIB double l2normsq<> (const TVector<VT> &v);

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
    friend MATHLIB TVector<VT> &append<> (TVector<VT> &v1, const TVector<VT> &v2);

    /**
     * \brief Concatenate two vectors
     *
     * Creates a new vector by concatenating the two arguments.
     * \return The result of the concatenation, {v1,v2}
     * \note Both arguments unchanged
     */
    friend MATHLIB TVector<VT> cat<> (const TVector<VT> &v1, const TVector<VT> &v2);

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
#if (defined(WIN32)||defined(WIN64))
    template <class VT> friend MATHLIB TVector<VT> operator+ (const VT &s,
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
#if (defined(WIN32)||defined(WIN64))
    template<class VT> friend MATHLIB TVector<VT> operator- (const VT &s,
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
#if (defined(WIN32)||defined(WIN64))
    template<class VT> friend MATHLIB TVector<VT> operator* (const VT &s,
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
#if (defined(WIN32)||defined(WIN64))
    template<class VT> friend MATHLIB TVector<VT> operator/ (const VT &s,
        const TVector<VT> &v);
    // JK also works (TVector::operator/ )
#elif (GCC_VERSION < 30404) // old-style definitions
    friend TVector<VT> operator/ <> (const VT &s, const TVector<VT> &v);
#else
    friend TVector<VT> (::operator/ <>) (const VT &s, const TVector<VT> &v);
#endif

    friend MATHLIB bool visnan<> (const TVector<VT> &v);

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
    friend MATHLIB std::ostream &operator<< <> (std::ostream &os,
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
    friend MATHLIB std::istream &operator>> <> (std::istream &is, TVector<VT> &v);

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
    friend MATHLIB TVector<double> UnfoldComplex<> (const TVector<VT> &v);

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
    inline void Allocate (int dim);

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

    /// \brief Unlink vector from its data block.
    /// \remarks This results in a (valid) vector of length zero.
    /// \sa New(), Clear()
    inline void Unlink ();

    int size;		 ///< vector length
    int *base_size;	 ///< pointer to length of linked data block
    int *base_nref;	 ///< pointer to data block reference counter
    VT *data;		 ///< pointer to first vector element
    VT *base_data;	 ///< pointer to linked data block
};

// ==========================================================================
// ==========================================================================
// Member definitions

// --------------------------------------------------------------------------
// constructor (vector of size 0)

template<class VT>
TVector<VT>::TVector ()
{
    base_nref = 0;
    data = 0;
    size = 0;
}

// --------------------------------------------------------------------------
// constructor (vector of size 'dim' with zero elements)

template<class VT>
TVector<VT>::TVector (int dim)
{
    dASSERT(dim >= 0, Parameter 1 must be >= 0);
    base_nref = 0;
    Allocate (dim);
}

// --------------------------------------------------------------------------
// constructor (uniform element assignment from scalar)

template<class VT>
TVector<VT>::TVector (int dim, const VT s)
{
    dASSERT(dim >= 0, Parameter 1 must be >= 0);
    base_nref = 0;
    Allocate (dim);
    *this = s;
}

// --------------------------------------------------------------------------
// constructor (assign elements from array)

template<class VT>
TVector<VT>::TVector (int dim, VT *values)
{
    dASSERT(dim >= 0, Parameter 1 must be >= 0);
    base_nref = 0;
    Allocate (dim);
    memcpy (data, values, dim*sizeof(VT));
}

// --------------------------------------------------------------------------
// constructor (element assignment from string)

template<class VT>
TVector<VT>::TVector (int dim, const char *init)
{
    dASSERT(dim >= 0, Parameter 1 must be >= 0);
    base_nref = 0;
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
    Allocate (v.size);
    Copy (v);				// copy elements from vec
}

// --------------------------------------------------------------------------
// reference constructor

template<class VT>
TVector<VT>::TVector (const TVector<VT> &v, int ofs, int dim)
{
    dASSERT(ofs >= 0 && ofs < v.Dim(), Parameter 1 index out of range);
    dASSERT(dim >= 0, Parameter 3 must be >= 0);
    dASSERT(ofs+dim <= v.Dim(), Data block of reference vector must be contained in original vector);
    base_nref = 0;
    Link (v, ofs, dim);
}

// --------------------------------------------------------------------------
// allocate a data block for the vector. This function must be called by
// a constructor, or after a call to Unlink()
template<class VT>
inline void TVector<VT>::Allocate (int dim)
{
    dASSERT(!base_nref, Data block present. Use Unlink first.);
    if (dim) {
	char *base = (char*)malloc (2*sizeof(int) + dim*sizeof(VT));
	dASSERT(base, Memory allocation failed.);
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
// Remove the vector's link to the data block
template<class VT>
inline void TVector<VT>::Unlink ()
{
    if (base_nref) {
	if (--(*base_nref) == 0)   // was last link
	    free (base_nref);      // -> deallocate base block
	base_nref = 0;             // invalidate base pointer
    }
    data = 0;		           // invalidate data pointer
    size = 0;			   // set dimension zero
}

// --------------------------------------------------------------------------
// Vector->Vector copy: copy from v to *this

template<class VT>
inline void TVector<VT>::Copy (const TVector<VT> &v)
{
    if (size != v.size) New (v.size); // reallocate
    memcpy (data, v.data, size * sizeof (VT));
}
#ifdef USE_BLAS_LEVEL1
template<>
inline void TVector<double>::Copy (const TVector<double> &v)
{
    static int incr = 1;
    if (size != v.size) New (v.size); // reallocate
    dcopy_(size, v.data, incr, data, incr);
}
template<>
inline void TVector<float>::Copy (const TVector<float> &v)
{
    static int incr = 1;
    if (size != v.size) New (v.size); // reallocate
    scopy_(size, v.data, incr, data, incr);
}
#endif // USE_BLAS_LEVEL1

// --------------------------------------------------------------------------
// partial Vector->Vector copy: copy part of v into *this

template<class VT>
inline void TVector<VT>::Copy (const TVector<VT> &v,
    int tofs, int sofs, int n)
{
    if (n < 0) n = v.size - sofs;
    if (n > size - tofs) n = size - tofs;
    memcpy (data+tofs, v.data+sofs, n * sizeof (VT));
}
#ifdef USE_BLAS_LEVEL1
template<>
inline void TVector<double>::Copy (const TVector<double> &v,
    int tofs, int sofs, int n)
{
    static int incr = 1;
    if (n < 0) n = v.size - sofs;
    if (n > size - tofs) n = size - tofs;
    dcopy_(n, v.data+sofs, incr, data+tofs, incr);
}
template<>
inline void TVector<float>::Copy (const TVector<float> &v,
    int tofs, int sofs, int n)
{
    static int incr = 1;
    if (n < 0) n = v.size - sofs;
    if (n > size - tofs) n = size - tofs;
    scopy_(n, v.data+sofs, incr, data+tofs, incr);
}
#endif // USE_BLAS_LEVEL1

// --------------------------------------------------------------------------

template<class VT>
inline TVector<VT> &TVector<VT>::operator= (VT s)
{
    for (int i = 0; i < size; i++) data[i] = s;
    return *this;
}

// --------------------------------------------------------------------------

template<class VT>
inline VT &TVector<VT>::operator[] (int i) const
{
    dASSERT(i >= 0 && i < size, Index out of range);
    return data[i];
}

// --------------------------------------------------------------------------

#ifndef MATH_DEBUG
template<class VT>
inline TVector<VT> &TVector<VT>::operator+= (const TVector<VT> &v)
{
    for (int i = 0; i < size; i++) data[i] += v.data[i];
    return *this;
}
#endif // !MATH_DEBUG

template<class VT>
inline TVector<VT> &TVector<VT>::operator+= (const VT &s)
{
    for (int i = 0; i < size; i++) data[i] += s;
    return *this;
}

#ifndef MATH_DEBUG
template<class VT>
inline TVector<VT> &TVector<VT>::operator-= (const TVector<VT> &v)
{
    for (int i = 0; i < size; i++) data[i] -= v.data[i];
    return *this;
}
#endif // !MATH_DEBUG

template<class VT>
inline TVector<VT> &TVector<VT>::operator-= (const VT &s)
{
    for (int i = 0; i < size; i++) data[i] -= s;
    return *this;
}

#ifndef MATH_DEBUG
template<class VT>
inline TVector<VT> &TVector<VT>::operator*= (const TVector<VT> &v)
{
    for (int i = 0; i < size; i++) data[i] *= v.data[i];
    return *this;
}
#endif // !MATH_DEBUG

// ==========================================================================

template<class VT>
inline TVector<VT> &TVector<VT>::operator*= (const VT &s)
{
    for (int i = 0; i < size; i++) data[i] *= s;
    return *this;
}
#ifdef USE_BLAS_LEVEL1
template<>
inline TVector<double> &TVector<double>::operator*= (const double &s)
{
    static int incr = 1;
    dscal_(size, s, data, incr);
    return *this;
}
template<>
inline TVector<float> &TVector<float>::operator*= (const float &s)
{
    static int incr = 1;
    sscal_(size, s, data, incr);
    return *this;
}
template<>
inline TVector<toast::complex> &TVector<toast::complex>::operator*=
    (const toast::complex &s)
{
    static int incr = 1;
    zscal_(size, s, data, incr);
    return *this;
}
#endif // USE_BLAS_LEVEL1

// ==========================================================================

#ifndef MATH_DEBUG
template<class VT>
inline TVector<VT> &TVector<VT>::operator/= (const TVector<VT> &v)
{
    for (int i = 0; i < size; i++) {
        data[i] /= v.data[i];
    }
    return *this;
}
#endif // !MATH_DEBUG

template<class VT>
inline TVector<VT> &TVector<VT>::operator/= (const VT &s)
{
    dASSERT (s, Attempt to divide by zero);
    *this *= ((VT)1/s);
    return *this;
}

// Link vector to a different data block
template<class VT>
inline void TVector<VT>::Relink (const TVector<VT> &v)
{
    Unlink ();				// remove existing link
    Link (v);                           // link to v's data block
}

// Relink vector to part of a different data block
template<class VT>
inline void TVector<VT>::Relink (const TVector<VT> &v, int ofs, int dim)
{
    Unlink ();				// remove existing link
    Link (v, ofs, dim);                 // link into v's data block
}

template<class VT>
inline void TVector<VT>::Clear ()
{
    memset ((void*)data, 0, size*sizeof(VT));
    // warning: this assumes that "(VT)0" is represented by a sequence of
    // "0" bytes - not necessarily correct
}

#ifndef MATH_DEBUG
template<class VT>
inline VT dot (const TVector<VT> &v1, const TVector<VT> &v2)
{
    VT d = (VT)0;
    for (int i = 0; i < v1.size; i++) d += v1[i] * v2[i];
    return d;
}
#ifdef USE_BLAS_LEVEL1 // interface to BLAS level1 xDOT/ZDOTU functions
template<>
inline double dot (const TVector<double> &v1, const TVector<double> &v2)
{
    static int incr = 1;
    return ddot_((int&)v1.size, v1.data, incr, v2.data, incr);
}
template<>
inline float dot (const TVector<float> &v1, const TVector<float> &v2)
{
    static int incr = 1;
    return sdot_((int&)v1.size, v1.data, incr, v2.data, incr);
}
template<>
inline toast::complex dot (const TVector<toast::complex> &v1,
    const TVector<toast::complex> &v2)
{
    static int incr = 1;
    int size = v1.size;
    dcomplex z = zdotu_(&size, (dcomplex*)v1.data, &incr,
			       (dcomplex*)v2.data, &incr);
    return toast::complex(z.r, z.i);
}
#endif // USE_BLAS_LEVEL1
#endif // !MATH_DEBUG

#ifndef MATH_DEBUG
template<class VT>
inline VT doth (const TVector<VT> &v1, const TVector<VT> &v2)
{
    return dot (v1, v2);
}
template<class VT>
inline toast::complex doth (const TVector<toast::complex> &v1,
    const TVector<toast::complex> &v2)
{
#ifdef USE_BLAS_LEVEL1
    static int incr = 1;
    return zdotc_((int&)v1.size, v1.data, incr, v2.data, incr);
#else
    toast::complex d = (toast::complex)0;
    for (int i = 0; i < v1.size; i++) d += conj(v1[i]) * v2[i];
    return d;
#endif // USE_BLAS_LEVEL1
}
#endif // !MATH_DEBUG

template<class VT>
inline MATHLIB double l2norm (const TVector<VT> &v)
{
    return sqrt (l2normsq (v));
}
#ifdef USE_BLAS_LEVEL1
template<>
inline double l2norm (const TVector<double> &v)
{
    static int incr = 1;
    return dnrm2_((int&)v.size, v.data, incr);
}
template<>
inline double l2norm (const TVector<float> &v)
{
    static int incr = 1;
    return (double)snrm2_((int&)v.size, v.data, incr);
}
template<>
inline double l2norm (const TVector<toast::complex> &v)
{
    static int incr = 1;
    return dznrm2_((int&)v.size, v.data, incr);
}
#endif // USE_BLAS_LEVEL1

// ==========================================================================
// template typedefs

typedef TVector<double>   RVector;
typedef TVector<float>    FVector;
typedef TVector<toast::complex>  CVector;
typedef TVector<scomplex> SCVector;
typedef TVector<int>      IVector;

/* add specific function for CVector */
MATHLIB FVector Re(const SCVector &);
MATHLIB RVector Re(const CVector &);
MATHLIB FVector Im(const SCVector &);
MATHLIB RVector Im(const CVector &);
RVector Mod(const CVector &);
MATHLIB RVector LogMod(const CVector &);
MATHLIB RVector Arg (const CVector &);
CVector Conj(const CVector &); // obsolete
void SelfConj(const CVector &);
CVector Hadamard(const CVector &, const CVector &);
MATHLIB void SetReal (CVector&, const RVector&);
MATHLIB void SetImag (CVector&, const RVector&);

CVector MakeCVector (const SCVector &v);

// ==========================================================================
// extern declarations of TVector (only required for VS)

#ifndef __VECTOR_CC
extern template class MATHLIB TVector<double>;
extern template class MATHLIB TVector<float>;
extern template class MATHLIB TVector<toast::complex>;
extern template class MATHLIB TVector<scomplex>;
extern template class MATHLIB TVector<int>;
#endif // !__VECTOR_CC

#endif // !__VECTOR_H
