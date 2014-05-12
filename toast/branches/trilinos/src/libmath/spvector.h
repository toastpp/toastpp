// -*-C++-*-
// ==========================================================================
// Module mathlib
// File spvector.h
// Declaration of template class TSparseVector ('template sparse vector')
//
// The following typedefs for specific template types have been defined
// for convenience:
//	RSparseVector = TSparseVector<double>	('real')
//	FSparseVector = TSparseVector<float>	('float')
//	CSparseVector = TSparseVector<complex>	('complex')
//	ISparseVector = TSparseVector<int>	('integer')
//	SparseVector  = TSparseVector<double>	for backward compatibility
// ==========================================================================

#ifndef __SPVECTOR_H
#define __SPVECTOR_H

template<class VT> class TSparseMatrix;

// ==========================================================================
// struct for defining a sparse vector data block

template<class VT>
struct spvecbase {
    VT *data;  		// data array
    int *index;		// index array
    int lsize;		// logical vector size
    int psize;		// physical vector size (= number of used entries)
    int bufsize;        // buffer size (= total number of entries)
    int nref;		// reference counter
};

// ==========================================================================
// class TSparseVector

template<class VT> class TSparseVector {
public:
    TSparseVector ();

    TSparseVector (int ldim);
    // creates a sparse vector with logical dimension ldim and physical
    // dimension 0

    TSparseVector (int ldim, int pdim);
    // creates a sparse vector with logical dimension ldim and physical
    // dimension pdim (all entries unused)

    TSparseVector (const TSparseVector<VT> &vec);
    // copy constructor; *this is created as a copy of, not a reference to
    // `vec'

    ~TSparseVector ();

    operator TVector<VT> ();
    // cast to full vector

    void Copy (const TSparseVector<VT> &vec);
    void Copy (const TVector<VT> &vec);
    // copy from sparse or full vector

    void New (int ldim, int nbuf = 0);
    // reallocates vector; sets logical size to `ldim' and allocates memory
    // for `nbuf' entries

    void Allocate (int nbuf);  // previously `Initialise'
    // reallocates memory for `nbuf' entries, but keeps logical size

    void Relink (const TSparseVector<VT> &vec);

    void Shrink ();
    // removes all unused entries from the vector

    void Initialise (int nz, int *idx);
    // initialise vector to `nz' entries and index list `idx'
    // values are set to zero

    void Clear ();
    // sets all entries to zero but without deallocating space

    TSparseVector &operator*= (const VT &vt);
    // multiplication with scalar

    VT operator& (const TVector<VT> &vec) const;
    // scalar product with a full vector

    friend VT iprod (const TSparseVector<VT> &v1, const TSparseVector<VT> &v2,
        int from = 0, int to = -1);
    // scalar product of two sparse vectors with optional range limit
    // The range specification [from,to] denotes logical indices

    int LDim() const { return base->lsize; }
    // returns logical dimension of vector

    int PDim() const { return base->psize; }
    // returns physical dimension of data array

    VT Get (int i) const;		// i: logical index

    void Put (int i, VT val);

    void Add (int i, VT val);

    void Erase (int i);
    // deletes entry for logical index i

    inline int PIndex (int i) const;
    // returns physical index to logical index `i', or -1 if `i' is not
    // in the index list

    int Insert (int i);
    // inserts a new entry for logical index `i' and returns the
    // corresponding physical index
    // assumes that an entry for `i' is not already present

    void SparseOutput (ostream &os) const;

// conditional inline functions
#ifdef MATH_DEBUG
    VT &operator[] (int i) const;	// i: physical index

    int LIndex (int i) const;
    // returns the logical index of physical index `i'.

    void SetIndex (int pi, int li);
    // sets index at position pi to logical index li.
    // Warning: this is a low-level routine. Use it only when you know
    // what you are doing.
#else
    VT &operator[] (int i) const { return data[i]; }

    int LIndex (int i) const { return index[i]; }

    void SetIndex (int pi, int li) { index[pi] = li; }
    // sets index at position pi to logical index li.
    // Warning: this is a low-level routine. Use it only when you know
    // what you are doing.
#endif

    // friends
    friend TVector<VT> ToVector (const TSparseVector<VT> &vec);
    friend ostream &operator<< (ostream &os, const TSparseVector<VT> &vec);

private:
    spvecbase<VT> *base;
    VT *data;
    int *index;

    friend VT iprod_sorted (const TSparseVector<VT> &v1,
        const TSparseVector<VT> &v2, int from = 0, int to = -1);
    // as iprod, but assumes index lists of v1 and v2 are sorted in
    // ascending order

    friend bool overlap_sorted (const TSparseVector<VT> &v1,
        const TSparseVector<VT> &v2, int from = 0, int to = -1);
    // returns TRUE if sorted vectors v1 and v2 share common nonzero entries
    // in the (logical) range between 'from' and 'to'

    friend class TSparseMatrix<VT>;
};


// ==========================================================================
// inline functions

template<class VT>
inline TSparseVector<VT>::TSparseVector ()
{
    base = new spvecbase<VT>;
    data = base->data;
    index = base->index;
    base->nref = 1;
    base->bufsize = base->psize = base->lsize = 0;
}

template<class VT>
inline TSparseVector<VT>::TSparseVector (int ldim)
{
    base = new spvecbase<VT>;
    data = base->data;
    index = base->index;
    base->nref = 1;
    base->bufsize = base->psize = 0;
    base->lsize = ldim;
}

template<class VT>
inline TSparseVector<VT>::TSparseVector (int ldim, int pdim)
{
    base = new spvecbase<VT>;
    data = base->data;
    index = base->index;
    base->nref = 1;
    base->bufsize = 0;
    New (ldim, pdim);
}

template<class VT>
inline TSparseVector<VT>::~TSparseVector ()
{
    if (--base->nref == 0) {  // destroy base only if this was the last link
	if (base->bufsize) {
	    delete []base->data;
	    delete []base->index;
	}
	delete base;
    }
}

template<class VT>
inline int TSparseVector<VT>::PIndex (int i) const
{
    for (int j = 0; j < base->psize; j++)
	if (index[j] == i) return j;
    return -1;
}

#ifndef MATH_DEBUG
template<class VT>
inline VT TSparseVector<VT>::operator& (const TVector<VT> &vec) const
{
    VT sum = (VT)0;
    int i, ps = base->psize;

    for (i = 0; i < ps; i++)
        sum += data[i] * vec[index[i]];
    return sum;
}
#endif // !MATH_DEBUG

// ==========================================================================
// friend prototypes

#ifdef NEED_FRIEND_PT

template<class VT>
bool overlap_sorted (const TSparseVector<VT> &v1, const TSparseVector<VT> &v2,
        int from = 0, int to = -1);

template<class VT>
VT iprod (const TSparseVector<VT> &v1, const TSparseVector<VT> &v2,
        int from = 0, int to = -1);

template<class VT>
VT iprod_sorted (const TSparseVector<VT> &v1, const TSparseVector<VT> &v2,
        int from = 0, int to = -1);

template<class VT>
TVector<VT> ToVector (const TSparseVector<VT> &vec);

template<class VT>
ostream &operator<< (ostream &os, const TSparseVector<VT> &vec);

#endif

// ==========================================================================
// typedefs for specific instances of `TSparseVector'

typedef TSparseVector<double>	RSparseVector;	// 'real'
typedef TSparseVector<float>	FSparseVector;	// 'float'
typedef TSparseVector<complex>	CSparseVector;	// 'complex'
typedef TSparseVector<int>	ISparseVector;	// 'integer'
typedef TSparseVector<double>	SparseVector;	// for backward compatibility

#endif // !__SPVECTOR_H
