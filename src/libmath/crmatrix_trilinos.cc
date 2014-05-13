#include "crmatrix_trilinos.h"

class zoltanMesh;


template<class MT>
TCompRowMatrixTrilinos<MT>::TCompRowMatrixTrilinos ()
  : TGenericSparseMatrix<MT> ()
{
    rowptr = 0;
    colidx = 0;
    map = Teuchos::null;
}

// --------------------------------------------------------------------------

template<class MT>
TCompRowMatrixTrilinos<MT>::TCompRowMatrixTrilinos (int rows, int cols,
    int *_rowptr, int *_colidx,  Teuchos::RCP<const map_type> _map)
  : TGenericSparseMatrix<MT> (rows, cols, _rowptr[rows])
{
    rowptr = new idxtype[rows+1];
    colidx = new idxtype[this->nval];

    memcpy (rowptr, _rowptr, (rows+1)*sizeof(idxtype));
    memcpy (colidx, _colidx, (this->nval)*sizeof(idxtype));
    map = _map;
}

// --------------------------------------------------------------------------

template<class MT>
TCompRowMatrixTrilinos<MT>::~TCompRowMatrixTrilinos ()
{
    if (rowptr) delete []rowptr;
    if (colidx) delete []colidx;
}

// --------------------------------------------------------------------------

template<class MT>
void TCompRowMatrixTrilinos<MT>::Assemble (
    const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
    const Teuchos::RCP<typename matrix_type::node_type>& node,
    zoltanMesh *mesh, void *data,
    void (*assemble_clbk)(zoltanMesh *mesh, void *data, 
                          TCompRowMatrix<MT> &localMatrix))
{
    using Teuchos::ArrayView;
    using Teuchos::ArrayRCP;
    using Teuchos::RCP;
    using Teuchos::arcp;

    int i, j;

    const int numGlobElements = map->getGlobalNumElements();
    const int numMyElements = map->getNodeNumElements();
    const int indexBase = 0;

    ArrayView<const int>myGlobalElements =
        map->getNodeElementList();
    ArrayRCP<size_t> NumNz = arcp<size_t>(numMyElements);
    int numNzero = 0;

    // Assign number of elements for each row we are owning
    for (i = 0; i < numMyElements; i++) {
        int r = myGlobalElements[i];
	NumNz[i] = rowptr[r+1]-rowptr[r];
	numNzero += NumNz[i];
    }

    // Create a local matrix for the assembly
    int *myRowptr = new int[numMyElements+1];
    int *myColidx = new int[numNzero];
    myRowptr[0] = 0;
    for (i = 0; i < numMyElements; i++) {
	myRowptr[i+1] = myRowptr[i] + NumNz[i];
	for (j = 0; j < NumNz[i]; j++)
	    myColidx[myRowptr[i]+j] = colidx[rowptr[myGlobalElements[i]]+j];
    }
    TCompRowMatrix<MT> Alocal(numMyElements, numGlobElements,
        myRowptr, myColidx);
    assemble_clbk(mesh, data, Alocal);
    
    teuchosMatrix = rcp (new matrix_type (map, NumNz,
					       Tpetra::StaticProfile));

    // assemble local matrix into Tpetra matrix
    double *Alocal_val = Alocal.ValPtr();
    for (i = 0; i < numMyElements; i++) {
	int rp = myRowptr[i];
	ArrayView<int> cols (myColidx+rp, NumNz[i]);
	ArrayView<double> vals (Alocal_val+rp, NumNz[i]);
	teuchosMatrix->insertGlobalValues (myGlobalElements[i], cols, vals);
    }

    // We are done with NumNZ; free it.
    NumNz = Teuchos::null;
    delete []rowptr;
    delete []colidx;
    delete []myRowptr;
    delete []myColidx;

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
