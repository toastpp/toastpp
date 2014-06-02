#include "crmatrix_trilinos.h"

class zoltanMesh;


// --------------------------------------------------------------------------

template<class MT>
TCompRowMatrixTrilinos<MT>::TCompRowMatrixTrilinos ()
  : TGenericSparseMatrix<MT> ()
{
    rowptr = 0;
    colidx = 0;
    numMyElements = 0;
    map = Teuchos::null;
}

// --------------------------------------------------------------------------

template<class MT>
TCompRowMatrixTrilinos<MT>::TCompRowMatrixTrilinos (
    int rows, int cols, int *_rowptr, int *_colidx,
    Teuchos::RCP<const map_type>_map)
  : TGenericSparseMatrix<MT> (rows, cols, _rowptr[rows])
{
    rowptr = new idxtype[rows+1];
    colidx = new idxtype[this->nval];
    memcpy (rowptr, _rowptr, (rows+1)*sizeof(idxtype));
    memcpy (colidx, _colidx, this->nval*sizeof(idxtype));
    SetPartition (_map);
}

// --------------------------------------------------------------------------

template<class MT>
TCompRowMatrixTrilinos<MT>::~TCompRowMatrixTrilinos ()
{
    if (rowptr) delete []rowptr;
    if (colidx) delete []colidx;
    NumNz = Teuchos::null; // is this necessary for a managed pointer?
}

// --------------------------------------------------------------------------

template<class MT>
void TCompRowMatrixTrilinos<MT>::SetPartition (
    Teuchos::RCP<const map_type> _map)
{
    map = _map;
    numMyElements = map->getNodeNumElements();
    numGlobElements = map->getGlobalNumElements();

    // global indices of my nodes
    myGlobalElements = map->getNodeElementList();

    // number of nonzeros on each of my matrix rows
    numNzero = 0;
    NumNz = Teuchos::arcp<size_t> (numMyElements);
    for (int i = 0; i < numMyElements; ++i) {
	int r = myGlobalElements[i];
	NumNz[i] = (size_t)(rowptr[r+1]-rowptr[r]);
	numNzero += NumNz[i];
    }

    // Create the toast interface for the local map
    int *myRowptr = new int[numMyElements+1];
    int *myColidx = new int[numNzero];
    myRowptr[0] = 0;
    for (int i = 0; i < numMyElements; i++) {
	myRowptr[i+1] = myRowptr[i] + NumNz[i];
	for (int j = 0; j < NumNz[i]; j++)
	    myColidx[myRowptr[i]+j] = colidx[rowptr[myGlobalElements[i]]+j];
    }
    toastMatrix.New (numMyElements, numGlobElements);
    toastMatrix.Initialise (myRowptr, myColidx);

    // create the distributed Teuchos matrix
    teuchosMatrix = rcp (new matrix_type (map, NumNz, Tpetra::StaticProfile));
}

// --------------------------------------------------------------------------

template<class MT>
void TCompRowMatrixTrilinos<MT>::Assemble (zoltanMesh *mesh, void *data,
    void(*clbkAssemble)(TCompRowMatrix<MT>&,zoltanMesh*,void*))
{
    using Teuchos::arcp;
    using Teuchos::ArrayRCP;
    using Teuchos::ArrayView;

    size_t i, j;
    const idxtype *myRowptr, *myColidx;

    // Assemble my part of the matrix
    clbkAssemble (toastMatrix, mesh, data);

    // Copy to the distributed Teuchos matrix
    MT *val = toastMatrix.ValPtr();
    toastMatrix.GetSparseStructure (&myRowptr, &myColidx);

    for (i = 0; i < numMyElements; i++) {
	int rp = myRowptr[i];
	ArrayView<const idxtype> cols (myColidx+rp, NumNz[i]);
	ArrayView<MT> vals (val+rp, NumNz[i]);
	teuchosMatrix->insertGlobalValues (myGlobalElements[i], cols, vals);
    }

    // Finish up the matrix.
    teuchosMatrix->fillComplete ();
}

// ==========================================================================
// class and friend instantiations

#ifdef NEED_EXPLICIT_INSTANTIATION

template class MATHLIB TCompRowMatrixTrilinos<double>;
//template class MATHLIB TCompRowMatrixTrilinos<float>;
//template class MATHLIB TCompRowMatrixTrilinos<toast::complex>;
//template class MATHLIB TCompRowMatrixTrilinos<scomplex>;
//template class MATHLIB TCompRowMatrixTrilinos<int>;

#endif // NEED_EXPLICIT_INSTANTIATION
