#include "toast.h"
using namespace toast;
/** Thresholds and shrinks a real dense matrix to give a real sparse matrix
* NOTE!! The threshold is set to 1e-15
**/
RCompRowMatrix shrink(const RDenseMatrix &dnsmat)
{
    int i, j;
    int *rowptr, *colidx;
    double *val;
    int m =  dnsmat.nRows(), n= dnsmat.nCols();

    rowptr = new int[m+1];
    rowptr[0] = 0;
    for (i = 0; i < m; i++) {
	int nz=0;
	for (j = 0; j < n; j++)
		if(!(fabs(dnsmat.Get(i, j)) < 1e-15))
			nz++;
	rowptr[i+1] = rowptr[i] + nz;
    }
    
    colidx = new int[rowptr[m]];
    val = new double[rowptr[m]];
    int k=0;
    for(i=0; i < m; i++){
	for(j=0; j<n; j++){
		if(!(fabs(dnsmat.Get(i, j)) < 1e-15)){
			colidx[k] = j; 
			val[k] = dnsmat.Get(i, j);
			k++;
		}
	}
    }
    RCompRowMatrix C (m, n, rowptr, colidx, val);

    delete []rowptr;
    delete []colidx;
    delete []val;	
    return C;
}


