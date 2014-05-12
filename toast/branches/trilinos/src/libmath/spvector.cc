// ==========================================================================
// Module mathlib
// File spvector.cc
// Definition of template class TSparseVector ('template sparse vector')
// ==========================================================================

#include <iostream.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "error.h"
#include "mathdef.h"
#include "complex.h"
#include "vector.h"
#include "spvector.h"

#define SPBLOCKSIZE 256  // chunk size for data block allocation

// ==========================================================================
// member definitions

template<class VT>
TSparseVector<VT>::TSparseVector (const TSparseVector<VT> &vec)
{
    base = new spvecbase<VT>;
    data = base->data;
    index = base->index;
    base->nref = 1;
    base->bufsize = 0;
    New (vec.base->lsize, vec.base->bufsize);
    base->psize = vec.base->psize;
    memcpy (data, vec.data, base->psize * sizeof (VT));
    memcpy (index, vec.index, base->psize * sizeof (int));
}

template<class VT>
TSparseVector<VT>::operator TVector<VT> ()
{
    TVector<VT> tmp (base->lsize);
    for (int i = 0; i < base->psize; i++)
	tmp[index[i]] = data[i];
    return tmp;
}

template<class VT>
void TSparseVector<VT>::Copy (const TSparseVector<VT> &vec)
{
    dASSERT(base->lsize == vec.base->lsize, Vectors have different size);
    Allocate (vec.base->bufsize);
    base->psize = vec.base->psize;
    memcpy (data, vec.data, base->psize * sizeof (VT));
    memcpy (index, vec.index, base->psize * sizeof (int));
}

template<class VT>
void TSparseVector<VT>::Copy (const TVector<VT> &vec)
{
    int i, j, n = vec.Dim();

    // pass 1: find sparseness of vec
    for (i = j = 0; i < n; i++)
        if (vec[i] != 0) j++;
    // re-allocate sparse vector if necessary
    base->lsize = n;
    if (base->psize != j) {
        Allocate (j);
	base->psize = j;
    }
    // pass 2: copy nonzero entries
    for (i = j = 0; i < vec.Dim(); i++) {
        if (vec[i] != 0) {
	    data[j] = vec[i];
	    index[j++] = i;
	}
    }
}

template<class VT>
void TSparseVector<VT>::New (int ldim, int nbuf)
{
    dASSERT(nbuf <= ldim, Buffer size greater than logical vector size.);
    base->lsize = ldim;
    Allocate (nbuf);
}

template<class VT>
void TSparseVector<VT>::Allocate (int nbuf)
{
    dASSERT(nbuf >= 0 && nbuf <= base->lsize, Buffer size out of range);
    if (base->bufsize) {
	delete []base->data;
	delete []base->index;
    }
    if ((base->bufsize = nbuf) != 0) {
	data  = base->data  = new VT[nbuf];
	index = base->index = new int[nbuf];
    }
    base->psize = 0;       // no entries used
}

template<class VT>
void TSparseVector<VT>::Relink (const TSparseVector<VT> &vec)
{
    // unlink current base
    if (--base->nref == 0) {
	if (base->bufsize) {
	    delete []base->data;
	    delete []base->index;
	}
	delete base;
    }

    // link new base
    base = vec.base;
    data = base->data;
    index = base->index;
    base->nref++;
}

template<class VT>
void TSparseVector<VT>::Shrink ()
{
    if (base->psize == base->bufsize) return;

    VT *tmpdata = new VT[base->psize];
    memcpy (tmpdata, data, base->psize * sizeof(VT));
    delete []base->data;
    data = base->data = tmpdata;

    int *tmpindex = new int[base->psize];
    memcpy (tmpindex, index, base->psize * sizeof(int));
    delete []base->index;
    index = base->index = tmpindex;

    base->bufsize = base->psize;
}

template<class VT>
void TSparseVector<VT>::Initialise (int nz, int *idx)
{
    dASSERT(nz <= base->lsize, Allocation request exceeds logical vector size);
    Allocate (nz);
    for (int i = 0; i < nz; i++) {
        index[i] = idx[i];
	data[i] = (VT)0;
    }
    base->psize = nz;
}

template<class VT>
void TSparseVector<VT>::Clear ()
{
    for (int i = 0; i < base->psize; i++)
        data[i] = (VT)0;
}

template<class VT>
TSparseVector<VT> &TSparseVector<VT>::operator*= (const VT &vt)
{
    for (int i = 0; i < base->psize; i++)
	data[i] *= vt;
    return *this;
}

#ifdef MATH_DEBUG // otherwise inline
template<class VT>
VT TSparseVector<VT>::operator& (const TVector<VT> &vec) const
{
    dASSERT (base->lsize == vec.Dim(), Vector sizes differ.);
    VT sum = (VT)0;
    int i, ps = base->psize;
    for (i = 0; i < ps; i++)
	sum += data[i] * vec[index[i]];
    return sum;
}
#endif // MATH_DEBUG

template<class VT>
VT TSparseVector<VT>::Get (int i) const
{
    dASSERT(i >= 0 && i < base->lsize, Index out of range.);
    int j = PIndex (i);
    if (j >= 0) return data[j];
    else return (VT)0;
}

template<class VT>
void TSparseVector<VT>::Put (int i, VT val)
{
    int j;

    if ((j = PIndex(i)) < 0) j = Insert (i);
    data[j] = val;
}

template<class VT>
void TSparseVector<VT>::Add (int i, VT val)
{
    int j;

    if ((j = PIndex(i)) < 0) {
	j = Insert (i);
	data[j] = val;
    } else
	data[j] += val;
}

template<class VT>
void TSparseVector<VT>::Erase (int i)
{
    int j;
    if ((j = PIndex(i)) >= 0) {
	if (j < base->psize-1) {
	    data[j]  = data[base->psize-1];
	    index[j] = index[base->psize-1];
	}
	base->psize--;
    }
}

template<class VT>
int TSparseVector<VT>::Insert (int i)
{
    dASSERT(i >= 0 && i < base->lsize, Index out of range);
    int ps = base->psize;

    // unused buffer area available?
    if (ps < base->bufsize) {
	data[ps] = (VT)0;
	index[ps] = i;
	base->psize++;
	return ps;
    }

    // allocate new chunk
    if (base->bufsize) {
        VT *tmpdata = new VT[base->bufsize*2];
	memcpy (tmpdata, base->data, ps * sizeof (VT));
	delete []base->data;
	data = base->data = tmpdata;
	int *tmpindex = new int[base->bufsize*2];
	memcpy (tmpindex, base->index, ps * sizeof (VT));
	delete []base->index;
	index = base->index = tmpindex;
	base->bufsize *= 2;
    } else {
        data = base->data = new VT[SPBLOCKSIZE];
	index = base->index = new int[SPBLOCKSIZE];
        base->bufsize = SPBLOCKSIZE;
    }

    data[ps] = (VT)0;
    index[ps] = i;
    base->psize++;
    return ps;
}

template<class VT>
void TSparseVector<VT>::SparseOutput (ostream &os) const
{
    // the entries in the index list are unsorted, so we sort them first
    // (bubble sort - replace!)

    int i, j, tmp;
    int *iindex = new int[base->psize];
    for (i = 0; i < base->psize; i++) iindex[i] = i;
    for (i = 0; i < base->psize-1; i++) {
	for (j = i+1; j < base->psize; j++) {
	    if (index[iindex[i]] > index[iindex[j]]) {
		tmp = iindex[i];
		iindex[i] = iindex[j];
		iindex[j] = tmp;
	    }
	}
    }
    for (i = 0; i < base->psize; i++) {
	os << index[iindex[i]] << " "
	   << data[iindex[i]] << " ";
    }
    delete []iindex;
}

#ifdef MATH_DEBUG
// the following functions are inline in the non-debugging version

template<class VT>
VT &TSparseVector<VT>::operator[] (int i) const
{
    dASSERT(i >= 0 && i < base->psize, Index out of range);
    return data[i];
}

template<class VT>
int TSparseVector<VT>::LIndex (int i) const
{
    dASSERT(i >= 0 && i < base->lsize, Index out of range.);
    return index[i];
}

template<class VT>
void TSparseVector<VT>::SetIndex (int pi, int li)
{
    dASSERT(pi >= 0 && pi < base->psize, Physical index out of range.);
    dASSERT(li >= 0 && li < base->lsize, Logical index out of range.);
    index[pi] = li;
}

#endif // MATH_DEBUG

// ==========================================================================
// friend definitions

template<class VT>
VT iprod (const TSparseVector<VT> &v1, const TSparseVector<VT> &v2,
        int from, int to)
{
    VT sum = (VT)0;
    int pdim1, i1, i2, ind;

    if (to < 0) to = (v1.LDim() < v2.LDim() ? v1.LDim() : v2.LDim());
    pdim1 = v1.PDim();
    for (i1 = 0; i1 < pdim1; i1++) {
        ind = v1.LIndex(i1);
	if (ind < from || ind > to) continue; // outside range
	if ((i2 = v2.PIndex (ind)) >= 0) sum += v1[i1]*v2[i2];
    }
    return sum;
}

// this version assumes index lists of v1 and v2 are sorted in
// ascending order
template<class VT>
VT iprod_sorted (const TSparseVector<VT> &v1, const TSparseVector<VT> &v2,
    int from, int to)
{
    VT sum = (VT)0;
    int pdim1, pdim2, i1, i2, ind1, ind2;

    if (to >= 0) {
        dASSERT(to < v1.LDim() && to < v2.LDim(), Index out of range);
    } else {
        to = min (v1.LDim(), v2.LDim());
    }

    pdim1 = v1.PDim();
    pdim2 = v2.PDim();
    if (!pdim1 || !pdim2) return sum; // nothing to do

    if ((ind1 = v1.index[i1=0]) > to) return sum;
    if ((ind2 = v2.index[i2=0]) > to) return sum;
    for (;;) {
        if (ind1 < ind2) {
	    do {
	        if (++i1 == pdim1) return sum; // beyond phys. end of v1
	    } while ((ind1 = v1.index[i1]) < ind2);
	    if (ind1 > to) return sum;         // beyond upper bound
	} else if (ind1 > ind2) {
	    do {
	        if (++i2 == pdim2) return sum; // beyond phys. end of v2 
	    } while ((ind2 = v2.index[i2]) < ind1);
	    if (ind2 > to) return sum;         // beyond upper bound
	}
	if (ind1 == ind2 && ind1 >= from) {
	    sum += v1[i1++]  * v2[i2++];
	    if (i1 == pdim1 || i2 == pdim2) return sum;
	    ind1 = v1.index[i1];
	    ind2 = v2.index[i2];
	}
    }
}

template<class VT>
bool overlap_sorted (const TSparseVector<VT> &v1,
		     const TSparseVector<VT> &v2, int from, int to)
{
    int pdim1, pdim2, i1, i2, ind1, ind2;

    if (to >= 0) {
        dASSERT(to < v1.LDim() && to < v2.LDim(), Index out of range);
    } else {
        to = min (v1.LDim(), v2.LDim());
    }

    pdim1 = v1.PDim();
    pdim2 = v2.PDim();
    if (!pdim1 || !pdim2) return false; // zero vector

    if ((ind1 = v1.index[i1=0]) > to) return false;
    if ((ind2 = v2.index[i2=0]) > to) return false;
    for (;;) {
        if (ind1 < ind2) {
	    do {
	        if (++i1 == pdim1) return false; // beyond phys. end of v1
	    } while ((ind1 = v1.index[i1]) < ind2);
	    if (ind1 > to) return false;         // beyond upper bound
	} else if (ind1 > ind2) {
	    do {
	        if (++i2 == pdim2) return false; // beyond phys. end of v2 
	    } while ((ind2 = v2.index[i2]) < ind1);
	    if (ind2 > to) return false;         // beyond upper bound
	}
	if (ind1 == ind2 && ind1 >= from) return true;
    }
}

template<class VT>
TVector<VT> ToVector (const TSparseVector<VT> &vec)
{
    TVector<VT> fvec(vec.base->lsize);
    for (int i = 0; i < vec.base->psize; i++)
	fvec[vec.index[i]] = vec.data[i];
    return fvec;
}

template<class VT>
ostream &operator<< (ostream &os, const TSparseVector<VT> &vec)
{
    int i, ind;
    os << '[';
    for (i = ind = 0; i < vec.base->lsize; i++) {
	if (i) os << ' ';
	if (ind < vec.base->psize && vec.index[ind] == i)
	    os << vec.data[ind++];
	else
	    os << 0;
    }
    os << ']';
    return os;
}

// ==========================================================================
// class and friend instantiations

#ifdef NEED_EXPLICIT_INSTANTIATION

template class TSparseVector<double>;
template class TSparseVector<float>;
template class TSparseVector<complex>;
template class TSparseVector<int>;

template RVector ToVector (const RSparseVector &vec);
template FVector ToVector (const FSparseVector &vec);
template CVector ToVector (const CSparseVector &vec);
template IVector ToVector (const ISparseVector &vec);

template double iprod (const RSparseVector &v1, const RSparseVector &v2,
        int from, int to);

template double iprod_sorted (const RSparseVector &v1, const RSparseVector &v2,
        int from, int to);

template bool overlap_sorted (const ISparseVector &v1, const ISparseVector &v2,
        int from, int to);

template ostream &operator<< (ostream &os, const RSparseVector &vec);

#endif // NEED_EXPLICIT_INSTANTIATION

