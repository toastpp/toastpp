// -*-C++-*-
// ==========================================================================
// Module mathlib
// File cdmatrix.h
// Declaration of template class TCoordMatrix ('template coordinate storage
// matrix')
//
// The following typedefs for specific template types have been defined
// for convenience:
//	RCoordMatrix = TCoordMatrix<double>   ('real')
//	FCoordMatrix = TCoordMatrix<float>    ('float')
//	CCoordMatrix = TCoordMatrix<complex>  ('complex')
//	ICoordMatrix = TCoordMatrix<int>      ('integer')
//
// Inheritance:
// ------------
// TGenericSparseMatrix ----> TCoordMatrix
// ==========================================================================

#ifndef __CDMATRIX_H
#define __CDMATRIX_H

// ==========================================================================
// Nonmember declarations

template<class MT> class TCoordMatrix;
template<class MT> class TCompRowMatrix;

template<class MT>
TCoordMatrix<MT> cath (const TCoordMatrix<MT> &A,
    const TCoordMatrix<MT> &B);

template<class MT>
TCoordMatrix<MT> catv (const TCoordMatrix<MT> &A,
    const TCoordMatrix<MT> &B);

template<class MT>
std::istream &operator>> (std::istream &is, TCoordMatrix<MT> &m);

template<class MT>
std::ostream &operator<< (std::ostream &os, const TCoordMatrix<MT> &m);

// ==========================================================================
// class TCoordMatrix
// ==========================================================================
/**
 * \brief Coordinate-storage sparse matrix class
 *
 * This matrix class represents its data by a data array 'val' of length nz, a
 * column index array 'colidx' of length nz, and a row index array 'rowidx'
 * of length nz, where nz is the number of allocated elements.
 * Elements can be stored unsorted, which makes it easier to dynamically
 * re-allocate storage space, but reduces efficiency of various methods.
 * Sorting can be applied once the matrix structure is fixed to improve
 * performance.
 */
template<class MT> class TCoordMatrix: public TGenericSparseMatrix<MT> {

    friend class TCompRowMatrix<MT>;

public:

    /**
     * \brief Creates a coordinate matrix of dimension 0 x 0.
     */
    TCoordMatrix ();

    /**
     * \brief Creates a coordinate matrix of dimension rows x cols.
     * \param rows number of matrix rows
     * \param cols number of matrix columns
     * \note No data block is allocated, i.e. all matrix elements have a
     *   value of zero.
     */
    TCoordMatrix (int rows, int cols);

    TCoordMatrix (int rows, int cols, int vals, idxtype *_rowidx, idxtype *_colidx,
		  MT *_data = 0);
    TCoordMatrix (const TCoordMatrix<MT> &m);
    ~TCoordMatrix ();

    MatrixStorage StorageType() const { return MATRIX_COORD; }

    void New (int rows, int cols);
    // reset logical dimension to rows x cols
    // deallocate data and index lists, i.e. create an empty matrix

    void Unlink ();

    void Initialise (int nv, idxtype *_rowidx, idxtype *_colidx, MT *data = 0);
    // reallocate data block to size nv, and assign data and index lists

    TCoordMatrix &operator= (const TCoordMatrix<MT> &m);
    // assignment operator

    MT &operator() (int r, int c);
    MT Get (int r, int c) const;
    // access to element in row r, col c.
    // First version for write acces, second for read access

    bool Exists (int r, int c) const;
    // true if an entry for the element is allocated

    TVector<MT> Row (int r) const;
    // Returns row r as a vector

    TVector<MT> Col (int c) const;
    // Returns column c as a vector

    int SparseRow (int r, idxtype *ci, MT *rv) const;
    // See TRootMatrix
    // Note: 1. returned column indices may be unsorted
    //       2. This is expensive since complete matrix must be traversed

    void ColScale (const TVector<MT> &scale);
    // scales the columns with 'scale'

    void RowScale (const TVector<MT> &scale);
    // scales the rows with 'scale'

    MT GetNext (int &r, int &c) const;

    void Ax (const TVector<MT> &x, TVector<MT> &b) const;
    void Ax (const TVector<MT> &x, TVector<MT> &b, int r1, int r2) const;
    void ATx (const TVector<MT> &x, TVector<MT> &b) const;

    void Transpone ();
    // Replace matrix with its transpose

    inline MT &Val (int i) { return this->val[i]; }
    inline const MT &Val (int i) const { return this->val[i]; }
    // return i-th value of data vector

    void Sort (bool roworder=true) const;

    inline const idxtype *RowIndex () const { return rowidx; }
    inline const idxtype *ColIndex () const { return colidx; }

    inline idxtype RowIndex (int i) const { return rowidx[i]; }
    inline idxtype ColIndex (int i) const { return colidx[i]; }
    // Index list entries. These should probably not be public, but
    // can't be avoided now

    /**
     * \brief Concatenate two matrices horizontally
     * \param A first matrix argument
     * \param B second matrix argument
     * \return Result of concatenation [A B]
     * \note A and B must have the same number of rows.
     * \note This method has the functionality of the MATLAB construct
     *   C = (A,B)
     */
    friend TCoordMatrix<MT> cath<> (const TCoordMatrix<MT> &A,
				    const TCoordMatrix<MT> &B);
    
    /**
     * \brief Concatenate two matrices vertically
     * \param A first matrix argument
     * \param B second matrix argument
     * \return Result of concatenation
     *    \f[ C = \left[ \begin{array}{c} A\\B \end{array} \right] \f]
     * \note A and B must have the same number of columns
     * \note This method has the functionality of the MATLAB construct
     *   C = (A;B)
     */
    friend TCoordMatrix<MT> catv<> (const TCoordMatrix<MT> &A,
				    const TCoordMatrix<MT> &B);

    friend std::istream &operator>> <> (std::istream &is, TCoordMatrix<MT> &m);
    friend std::ostream &operator<< <> (std::ostream &os,
        const TCoordMatrix<MT> &m);

private:
    int Get_index (int r, int c) const;

    int Insert (int r, int c, MT v = 0);
    // create an entry for row r, col c, set value to v,
    // and return the index into the data array
    // Assumes that an entry for (r,c) does not exist yet!

    idxtype *rowidx;
    idxtype *colidx;
    mutable int iterator_pos;
};

// ==========================================================================
// typedefs for specific instances of `TCoordMatrix'

typedef TCoordMatrix<double>	RCoordMatrix;	// 'real'
typedef TCoordMatrix<float>	FCoordMatrix;	// 'float'
typedef TCoordMatrix<std::complex<double> >	CCoordMatrix;	// 'complex'
typedef TCoordMatrix<std::complex<float> >	SCCoordMatrix;	// 'single complex'
typedef TCoordMatrix<int>	ICoordMatrix;	// 'integer'

// ==========================================================================
// extern declarations of TCoordMatrix (only required for VS)

#ifndef __CDMATRIX_CC
//extern template class MATHLIB TCoordMatrix<double>;
//extern template class MATHLIB TCoordMatrix<float>;
//extern template class MATHLIB TCoordMatrix<toast::complex>;
//extern template class MATHLIB TCoordMatrix<scomplex>;
//extern template class MATHLIB TCoordMatrix<int>;
#endif // !__CDMATRIX_CC

#endif // !__CDMATRIX_H
