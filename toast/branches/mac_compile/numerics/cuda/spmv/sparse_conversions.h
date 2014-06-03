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

#include <algorithm>
#include "sparse_operations.h"

////////////////////////////////////////////////////////////////////////////////
//! Convert CSR format to DIA format
// If the matrix has more than 'max_diags' occupied diagonals, then a dia_matrix
// with dimensions (0,0) and 0 nonzeros is returned.
////////////////////////////////////////////////////////////////////////////////
template <class IndexType, class ValueType>
dia_matrix<IndexType, ValueType>
 csr_to_dia(const csr_matrix<IndexType,ValueType>& csr, const IndexType max_diags, const IndexType alignment = 16)
{
    dia_matrix<IndexType, ValueType> dia;
    dia.num_rows     = csr.num_rows;
    dia.num_cols     = csr.num_cols;
    dia.num_nonzeros = csr.num_nonzeros;
    dia.diag_offsets = NULL;
    dia.diag_data    = NULL;

    // compute number of occupied diagonals and enumerate them
    IndexType num_diags = 0;
    const IndexType unmarked = (IndexType) -1; // works for both signed and unsigned
    IndexType * diag_map = new_host_array<IndexType>(csr.num_rows + csr.num_cols);
    std::fill(diag_map, diag_map + csr.num_rows + csr.num_cols, unmarked);

    for(IndexType i = 0; i < csr.num_rows; i++){
        for(IndexType jj = csr.Ap[i]; jj < csr.Ap[i+1]; jj++){
            IndexType j = csr.Aj[jj];
            IndexType map_index = (csr.num_rows - i) + j; //offset shifted by + num_rows
            if(diag_map[map_index] == unmarked)
                diag_map[map_index] = num_diags++;
        }
    }
    
    dia.num_diags = num_diags;

    // Too many occupied diagonals
    if(num_diags > max_diags){
        dia.num_rows     = 0;
        dia.num_cols     = 0;
        dia.num_nonzeros = 0;
        dia.stride       = 0;
        delete_host_array(diag_map);
        return dia;
    }

    // length of each diagonal in memory
    dia.stride = alignment * ((dia.num_rows + alignment - 1)/ alignment);

    dia.diag_offsets = new_host_array<int>(dia.num_diags);
    dia.diag_data    = new_host_array<ValueType>(dia.num_diags * dia.stride);

    for(IndexType n = 0; n < csr.num_rows + csr.num_cols; n++)
        if(diag_map[n] != unmarked)
            dia.diag_offsets[diag_map[n]] = (int) n - (int) csr.num_rows;

    std::fill(dia.diag_data, dia.diag_data + dia.num_diags * dia.stride, ValueType(0));

    for(IndexType i = 0; i < csr.num_rows; i++){
        for(IndexType jj = csr.Ap[i]; jj < csr.Ap[i+1]; jj++){
            IndexType j = csr.Aj[jj];
            IndexType map_index = (csr.num_rows - i) + j; //offset shifted by + num_rows
            IndexType diag = diag_map[map_index];
            dia.diag_data[diag * dia.stride + i] = csr.Ax[jj];
        }
    }

    delete_host_array(diag_map);

    return dia;
}



////////////////////////////////////////////////////////////////////////////////
//! Convert CSR format to HYB (hybrid ELL/COO) format
// If the ELL portion of the HYB matrix will have 'num_cols_per_row' columns.
// Nonzero values that do not fit within the ELL structure are placed in the 
// COO format portion of the HYB matrix.
////////////////////////////////////////////////////////////////////////////////
template <class IndexType, class ValueType>
hyb_matrix<IndexType, ValueType>
 csr_to_hyb(const csr_matrix<IndexType,ValueType>& csr, const IndexType num_cols_per_row, const IndexType alignment = 16)
{
    hyb_matrix<IndexType, ValueType> hyb;

    ell_matrix<IndexType, ValueType> & ell = hyb.ell;
    coo_matrix<IndexType, ValueType> & coo = hyb.coo;

    hyb.num_rows = csr.num_rows;
    hyb.num_cols = csr.num_cols;
    hyb.num_nonzeros = csr.num_nonzeros;

    //initialize shapes
    ell.num_rows = csr.num_rows;
    ell.num_cols = csr.num_cols;
    coo.num_rows = csr.num_rows;
    coo.num_cols = csr.num_cols;
   
    ell.stride = alignment * ((ell.num_rows + alignment - 1)/ alignment);
    ell.num_cols_per_row = num_cols_per_row;

    // compute number of nonzeros in the ELL and COO portions
    ell.num_nonzeros = 0;
    for(IndexType i = 0; i < csr.num_rows; i++)
        ell.num_nonzeros += std::min(ell.num_cols_per_row, csr.Ap[i+1] - csr.Ap[i]); 

    coo.num_nonzeros = csr.num_nonzeros - ell.num_nonzeros;

    // allocate storage for ELL and COO matrices
    ell.Aj = new_host_array<IndexType>(ell.num_cols_per_row * ell.stride);
    ell.Ax = new_host_array<ValueType>(ell.num_cols_per_row * ell.stride);

    if(coo.num_nonzeros > 0){
        coo.I = new_host_array<IndexType>(coo.num_nonzeros);
        coo.J = new_host_array<IndexType>(coo.num_nonzeros);
        coo.V = new_host_array<ValueType>(coo.num_nonzeros);
    } else {
        coo.I = NULL;
        coo.J = NULL;
        coo.V = NULL;
    }

    // pad out ELL format with zeros
    std::fill(ell.Aj, ell.Aj + ell.num_cols_per_row * ell.stride, 0);
    std::fill(ell.Ax, ell.Ax + ell.num_cols_per_row * ell.stride, 0);

    for(IndexType i = 0, coo_nnz = 0; i < csr.num_rows; i++){
        IndexType n = 0;
        IndexType jj = csr.Ap[i];

        // copy up to num_cols_per_row values of row i into the ELL
        while(jj < csr.Ap[i+1] && n < ell.num_cols_per_row){
            ell.Aj[ell.stride * n + i] = csr.Aj[jj];
            ell.Ax[ell.stride * n + i] = csr.Ax[jj];
            jj++, n++;
        }

        // copy any remaining values in row i into the COO
        while(jj < csr.Ap[i+1]){
            coo.I[coo_nnz] = i;
            coo.J[coo_nnz] = csr.Aj[jj];
            coo.V[coo_nnz] = csr.Ax[jj];
            jj++; coo_nnz++;
        }
    }

    return hyb;
}


////////////////////////////////////////////////////////////////////////////////
//! Convert CSR format to ELL format
// If the matrix has more than 'max_cols_per_row' columns in any row, then 
// an ell_matrix with dimensions (0,0) and 0 nonzeros is returned. Rows with 
// fewer than 'num_cols_per_row' columns are padded with zeros.
////////////////////////////////////////////////////////////////////////////////
template <class IndexType, class ValueType>
ell_matrix<IndexType, ValueType>
 csr_to_ell(const csr_matrix<IndexType,ValueType>& csr, const IndexType max_cols_per_row, const IndexType alignment = 16)
{
    // compute maximum number of columns in any row
    IndexType num_cols_per_row = 0;
    for(IndexType i = 0; i < csr.num_rows; i++)
        num_cols_per_row = std::max(num_cols_per_row, csr.Ap[i+1] - csr.Ap[i]); 
    
    if(num_cols_per_row >= max_cols_per_row){
        //too many columns
        ell_matrix<IndexType, ValueType> ell;
        ell.Aj = NULL;
        ell.Ax = NULL;
        ell.num_rows = 0;
        ell.num_cols = 0;
        ell.num_nonzeros = 0;
        ell.stride = 0;
        ell.num_cols_per_row = num_cols_per_row;
        return ell;
    } else {
        // use CSR->HYB and grab the ELL portion
        return csr_to_hyb(csr, num_cols_per_row, alignment).ell;
    }
}


////////////////////////////////////////////////////////////////////////////////
//! Convert CSR format to COO format
// Storage for output is assumed to have been allocated
//! @param Ap             CSR pointer array
//! @param Aj             CSR index array
//! @param Ax             CSR data array
//! @param num_rows       number of rows
//! @param num_cols       number of columns
//! @param num_nonzeros   number of nonzeros
//! @param rows           COO row array
//! @param cols           COO column array
//! @param data           COO data array
////////////////////////////////////////////////////////////////////////////////
template <class IndexType, class ValueType>
void csr_to_coo(const IndexType * Ap,
                const IndexType * Aj,
                const ValueType  * Ax,
                const IndexType num_rows, 
                const IndexType num_cols, 
                const IndexType num_nonzeros,
                      IndexType * rows,
                      IndexType * cols,
                      ValueType * data)
{
    for(IndexType i = 0; i < num_rows; i++){        
        IndexType row_start = Ap[i];
        IndexType row_end   = Ap[i+1];
        for(IndexType jj = row_start; jj < row_end; jj++){
            rows[jj] = i;
        }
    }

    for(IndexType i = 0; i < num_nonzeros; i++){
        cols[i] = Aj[i];
        data[i] = Ax[i];
    }
}

template <class IndexType, class ValueType>
coo_matrix<IndexType, ValueType>
 csr_to_coo(const csr_matrix<IndexType,ValueType>& csr)
{   
    coo_matrix<IndexType, ValueType> coo;

    coo.num_rows     = csr.num_rows;
    coo.num_cols     = csr.num_cols;
    coo.num_nonzeros = csr.num_nonzeros;

    coo.I = new_host_array<IndexType>(csr.num_nonzeros);
    coo.J = new_host_array<IndexType>(csr.num_nonzeros);
    coo.V = new_host_array<ValueType>(csr.num_nonzeros);

    csr_to_coo(csr.Ap,csr.Aj,csr.Ax,
               coo.num_rows,coo.num_cols,coo.num_nonzeros,
               coo.I,coo.J,coo.V);


    return coo;
}


////////////////////////////////////////////////////////////////////////////////
//! Convert COO format to CSR format
// Storage for output is assumed to have been allocated
//! @param rows           COO row array
//! @param cols           COO column array
//! @param data           COO data array
//! @param num_rows       number of rows
//! @param num_cols       number of columns
//! @param num_nonzeros   number of nonzeros
//! @param Ap             CSR pointer array
//! @param Ai             CSR index array
//! @param Ax             CSR data array
////////////////////////////////////////////////////////////////////////////////
template <class IndexType, class ValueType>
void coo_to_csr(const IndexType * rows,
                const IndexType * cols,
                const ValueType * data,
                const IndexType num_rows, 
                const IndexType num_cols, 
                const IndexType num_nonzeros,
                      IndexType * Ap,
                      IndexType * Aj,
                      ValueType * Ax)
{
    for (IndexType i = 0; i < num_rows; i++)
        Ap[i] = 0;

    for (IndexType i = 0; i < num_nonzeros; i++)
        Ap[rows[i]]++;


    //cumsum the nnz per row to get Bp[]
    for(IndexType i = 0, cumsum = 0; i < num_rows; i++){     
        IndexType temp = Ap[i];
        Ap[i] = cumsum;
        cumsum += temp;
    }
    Ap[num_rows] = num_nonzeros;

    //write Aj,Ax into Bj,Bx
    for(IndexType i = 0; i < num_nonzeros; i++){
        IndexType row  = rows[i];
        IndexType dest = Ap[row];

        Aj[dest] = cols[i];
        Ax[dest] = data[i];

        Ap[row]++;
    }

    for(IndexType i = 0, last = 0; i <= num_rows; i++){
        IndexType temp = Ap[i];
        Ap[i]  = last;
        last   = temp;
    }
    
}


////////////////////////////////////////////////////////////////////////////////
//! Convert COOrdinate format (triplet) to CSR format
//! @param coo        coo_matrix
////////////////////////////////////////////////////////////////////////////////
template <class IndexType, class ValueType>
csr_matrix<IndexType, ValueType>
 coo_to_csr(const coo_matrix<IndexType,ValueType>& coo, bool compact = false){  

    csr_matrix<IndexType, ValueType> csr;

    csr.num_rows     = coo.num_rows;
    csr.num_cols     = coo.num_cols;
    csr.num_nonzeros = coo.num_nonzeros;

    csr.Ap = new_host_array<IndexType>(csr.num_rows + 1);
    csr.Aj = new_host_array<IndexType>(csr.num_nonzeros);
    csr.Ax = new_host_array<ValueType>(csr.num_nonzeros);

    coo_to_csr(coo.I, coo.J, coo.V,
               coo.num_rows, coo.num_cols, coo.num_nonzeros,
               csr.Ap, csr.Aj, csr.Ax);
    
    if (compact) {
        //sum duplicates together
        sum_csr_duplicates(csr.num_rows, csr.num_cols, csr.Ap, csr.Aj, csr.Ax);
        csr.num_nonzeros = csr.Ap[csr.num_rows];
    }

    return csr;
}

////////////////////////////////////////////////////////////////////////////////
//! Convert a csr_matrix to pkt_matrix format
////////////////////////////////////////////////////////////////////////////////

#include "csr_to_pkt.h"

