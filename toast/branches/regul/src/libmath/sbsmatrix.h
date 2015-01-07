// -*-C++-*-

#ifndef __SBSMATRIX_H
#define __SBSMATRIX_H

template<class MT> class TSparseBandSymMatrix: public TBandSymMatrix<MT> {
public:
    TSparseBandSymMatrix ();
    TSparseBandSymMatrix (int rc, int hb);
    TSparseBandSymMatrix (const TSparseBandSymMatrix<MT> &A);

    virtual ~TSparseBandSymMatrix ();

    void Register ();
    virtual void Ax (const TVector<MT> &x, TVector<MT> &b) const;
    friend void CHdecomp (TSparseBandSymMatrix<MT> &a);

private:
    int *rowfill;
    int *colfill;
    int **rowindex;
    int **colindex;
};

// ==========================================================================
// inline functions

// ==========================================================================
// typedefs for specific instances of `TSparseBandSymMatrix'

typedef TSparseBandSymMatrix<double>  RSparseBandSymMatrix;	// 'real'
typedef TSparseBandSymMatrix<float>   FSparseBandSymMatrix;	// 'float'
typedef TSparseBandSymMatrix<complex> CSparseBandSymMatrix;	// 'complex'
typedef TSparseBandSymMatrix<int>     ISparseBandSymMatrix;	// 'integer'
typedef TSparseBandSymMatrix<double>  SparseBandSymMatrix;	// for backward
                                                               // compatibility

#endif // !__SBSMATRIX_H
