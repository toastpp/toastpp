/*
 *  Copyright 2008-2009 NVIDIA Corporation
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 */



#pragma once

////////////////////////////////////////////////////////////////////////////////
//! Basic sparse matrix arithmetic and transformations
////////////////////////////////////////////////////////////////////////////////

#include <algorithm>
#include "sparse_formats.h"
#include "mem.h"


////////////////////////////////////////////////////////////////////////////////
//! Sum together the duplicate nonzeros in a CSR format
//! CSR format will be modified *in place*
//! @param num_rows       number of rows
//! @param num_cols       number of columns
//! @param Ap             CSR pointer array
//! @param Ai             CSR index array
//! @param Ax             CSR data array
////////////////////////////////////////////////////////////////////////////////
template <class IndexType, class ValueType>
void sum_csr_duplicates(const IndexType num_rows,
                        const IndexType num_cols, 
                              IndexType * Ap, 
                              IndexType * Aj, 
                              ValueType * Ax)
{
    IndexType * next = new_host_array<IndexType>(num_cols);
    ValueType * sums = new_host_array<ValueType>(num_cols);

    for(IndexType i = 0; i < num_cols; i++){
        next[i] = (IndexType) -1;
        sums[i] = (ValueType)   0;
    }

    IndexType NNZ = 0;

    IndexType row_start = 0;
    IndexType row_end   = 0;


    for(IndexType i = 0; i < num_rows; i++){
        IndexType head = (IndexType)-2;

        row_start = row_end; //Ap[i] may have been changed
        row_end   = Ap[i+1]; //Ap[i+1] is safe

        for(IndexType jj = row_start; jj < row_end; jj++){
            IndexType j = Aj[jj];

            sums[j] += Ax[jj];
            if(next[j] == (IndexType)-1){
                next[j] = head;                        
                head    = j;
            }
        }


        while(head != (IndexType)-2){
            IndexType curr = head; //current column
            head   = next[curr];

            if(sums[curr] != 0){
                Aj[NNZ] = curr;
                Ax[NNZ] = sums[curr];
                NNZ++;
            }

            next[curr] = (IndexType)-1;
            sums[curr] =  0;
        }
        Ap[i+1] = NNZ;
    }

    delete_host_array(next);
    delete_host_array(sums);
}
template <class IndexType, class ValueType>
void sum_csr_duplicates(csr_matrix<IndexType,ValueType>& A){
    sum_csr_duplicates(A.num_rows, A.num_cols, A.Ap, A.Aj, A.Ax);
    A.num_nonzeros = A.Ap[A.num_rows];
}


////////////////////////////////////////////////////////////////////////////////
//! Transpose a matrix in CSR format
//! Storage for B is assumed to have been allocated
//! @param Ap         CSR pointer array
//! @param Aj         CSR column index array
//! @param Ax         CSR data array
//! @param num_rows   number of rows in A
//! @param num_cols   number of columns in A
//! @param Bp         CSR pointer array
//! @param Bi         CSR column index array
//! @param Bx         CSR data array
////////////////////////////////////////////////////////////////////////////////
template <class IndexType, class ValueType>
void csr_transpose(const IndexType * Ap, 
                   const IndexType * Aj, 
                   const ValueType * Ax,
				   const IndexType num_rows, 
                   const IndexType num_cols, 
                         IndexType * Bp, 
                         IndexType * Bj, 
                         ValueType * Bx)
{	
    //TODO use temp-free method 
	IndexType * temp = new_host_array<IndexType>(num_cols);

    for(IndexType i = 0; i < num_cols; i++){
        temp[i] = 0;
    }

	IndexType num_nonzeros = Ap[num_rows];

	for(IndexType i = 0; i < num_nonzeros; i++) //count number of entries in each column
		temp[Aj[i]]++;

	Bp[0] = 0;
	for(IndexType i = 0; i < num_cols; i++){		
		Bp[i+1] = Bp[i] + temp[i];    //cumsum number column entries to form Bp
		temp[i] = 0;                  //count number of entries in each column
	}
	
	for(IndexType i = 0; i < num_rows; i++){
		IndexType row_start = Ap[i];
		IndexType row_end   = Ap[i+1];
		for(IndexType jj = row_start; jj < row_end; jj++){
			IndexType col    = Aj[jj];
			IndexType offset = temp[col] + Bp[col];

			Bj[offset] = i;
			Bx[offset] = Ax[jj];

			temp[col]++;
		}
	}		

    delete_host_array(temp);
}

template <typename IndexType, typename ValueType>
csr_matrix<IndexType, ValueType>
 csr_transpose(const csr_matrix<IndexType,ValueType>& csr){	
	csr_matrix<IndexType, ValueType> csr_t;

	csr_t.num_rows = csr.num_cols;
	csr_t.num_cols = csr.num_rows;
	csr_t.num_nonzeros = csr.num_nonzeros;

	csr_t.Ap = new_host_array<IndexType>(csr.num_cols + 1);
	csr_t.Aj = new_host_array<IndexType>(csr.num_nonzeros);
	csr_t.Ax = new_host_array<ValueType>(csr.num_nonzeros);

	csr_transpose(csr.Ap, csr.Aj, csr.Ax, 
                  csr.num_rows, csr.num_cols, 
                  csr_t.Ap, csr_t.Aj, csr_t.Ax);

	return csr_t;
}




////////////////////////////////////////////////////////////////////////////////
//! Compute Optimal Number of Columns per Row in the ELL part of the HYB format
//! Examines the distribution of nonzeros per row of the input CSR matrix to find
//! the optimal tradeoff between the ELL and COO portions of the hybrid (HYB)
//! sparse matrix format under the assumption that ELL performance is a fixed
//! multiple of COO performance.  Furthermore, since ELL performance is also
//! sensitive to the absolute number of rows (and COO is not), a threshold is
//! used to ensure that the ELL portion contains enough rows to be worthwhile.
//! The default values were chosen empirically for a GTX280.
//!
//! @param csr                  CSR matrix
//! @param relative_speed       Speed of ELL relative to COO (e.g. 2.0 -> ELL is twice as fast)
//! @param breakeven_threshold  Minimum threshold at which ELL is faster than COO
////////////////////////////////////////////////////////////////////////////////

// relative speed of full ELL vs. COO (full = no padding)
template <typename IndexType, typename ValueType>
IndexType compute_hyb_cols_per_row(const csr_matrix<IndexType,ValueType>& csr,
                                   float relative_speed = 3.0, IndexType breakeven_threshold = 4096)
{
    // compute maximum row length
    IndexType max_cols_per_row = 0;
    for(IndexType i = 0; i < csr.num_rows; i++)
        max_cols_per_row = std::max(max_cols_per_row, csr.Ap[i+1] - csr.Ap[i]); 

    // compute distribution of nnz per row
    IndexType * histogram = new_host_array<IndexType>(max_cols_per_row + 1);
    std::fill(histogram, histogram + max_cols_per_row + 1, 0);
    for(IndexType i = 0; i < csr.num_rows; i++)
        histogram[csr.Ap[i+1] - csr.Ap[i]]++;

    // compute optimal ELL column size 
    IndexType num_cols_per_row = max_cols_per_row;
    for(IndexType i = 0, rows = csr.num_rows; i < max_cols_per_row; i++)
    {
        rows -= histogram[i];  //number of rows of length > i
        if(relative_speed * rows < csr.num_rows || rows < breakeven_threshold)
        {
            num_cols_per_row = i;
            break;
        }
    }

    delete_host_array(histogram);

    return num_cols_per_row;
}


