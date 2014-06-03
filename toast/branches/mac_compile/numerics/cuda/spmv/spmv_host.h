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
//! CPU SpMV kernels
////////////////////////////////////////////////////////////////////////////////

#include "sparse_formats.h"

////////////////////////////////////////////////////////////////////////////////
//! Compute y += A*x for a sparse DIA matrix A and dense vectors x and y
//! @param num_rows   number of rows in A
//! @param num_cols   number of columns in A
//! @param num_diags  number of diagonals 
//! @param offsets    DIA offset array
//! @param data       DIA data array
//! @param x          column vector
//! @param y          column vector
////////////////////////////////////////////////////////////////////////////////
template <typename IndexType, typename ValueType>
void __spmv_dia_serial_host(const IndexType num_rows,
                            const IndexType num_cols,
                            const IndexType num_diags,
                            const IndexType stride,
                            const int       * offsets, 
                            const ValueType * data, 
                            const ValueType * x,
                                  ValueType * y)
{
    for(IndexType i = 0; i < num_diags; i++){
        const int k = offsets[i];  //diagonal offset

        const IndexType i_start = std::max((int)0,-k);
        const IndexType j_start = std::max((int)0, k);

        //number of elements to process
        const IndexType N = std::min(num_rows - i_start, num_cols - j_start);

        const ValueType * d_ = data + i*stride + i_start;
        const ValueType * x_ = x + j_start;
              ValueType * y_ = y + i_start;

        for(IndexType n = 0; n < N; n++){
            y_[n] += d_[n] * x_[n]; 
        }
    }
}

template <typename IndexType, typename ValueType>
void spmv_dia_serial_host(const dia_matrix<IndexType, ValueType>& dia, 
                          const ValueType * x,  
                                ValueType * y)
{
    __spmv_dia_serial_host(dia.num_rows, dia.num_cols, dia.num_diags, dia.stride,
                           dia.diag_offsets, dia.diag_data, 
                           x, y);
}

////////////////////////////////////////////////////////////////////////////////
//! Compute y += A*x for a sparse ELL matrix A and column vectors x and y
//! @param num_rows          number of rows in A
//! @param num_cols          number of columns in A
//! @param num_cols_per_row  number columns in each row (smaller rows are zero padded)
//! @param stride            seperation between row entries (stride >= num_rows, for alignment)
//! @param Aj                ELL column indices
//! @param Ax                ELL nonzero values
//! @param x                 column vector
//! @param y                 column vector
////////////////////////////////////////////////////////////////////////////////
template <typename IndexType, typename ValueType>
void __spmv_ell_serial_host(const IndexType num_rows,
                            const IndexType num_cols,
                            const IndexType num_cols_per_row,
                            const IndexType stride,
                            const IndexType * Aj, 
                            const ValueType * Ax, 
                            const ValueType * x,
                                  ValueType * y)
{
    for(IndexType n = 0; n < num_cols_per_row; n++){
        const IndexType * Aj_n = Aj + n * stride;
        const ValueType * Ax_n = Ax + n * stride;
        for(IndexType i = 0; i < num_rows; i++){
            y[i] += Ax_n[i] * x[Aj_n[i]];
        }
    }
}

template <typename IndexType, typename ValueType>
void spmv_ell_serial_host(const ell_matrix<IndexType, ValueType>& ell, 
                          const ValueType * x,  
                                ValueType * y)
{
    __spmv_ell_serial_host(ell.num_rows, ell.num_cols, ell.num_cols_per_row, ell.stride,
                           ell.Aj, ell.Ax,
                           x, y);
}



////////////////////////////////////////////////////////////////////////////////
//! Compute y += A*x for a sparse CSR matrix A and column vectors x and y
//! @param num_rows   number of rows in A
//! @param Ap         CSR pointer array
//! @param Aj         CSR index array
//! @param Ax         CSR data array
//! @param x          column vector
//! @param y          column vector
////////////////////////////////////////////////////////////////////////////////
template <typename IndexType, typename ValueType>
void __spmv_csr_serial_host(const IndexType num_rows, 
                            const IndexType * Ap, 
                            const IndexType * Aj, 
                            const ValueType * Ax, 
                            const ValueType * x,    
                                  ValueType * y)    
{
    for (IndexType i = 0; i < num_rows; i++){
        const IndexType row_start = Ap[i];
        const IndexType row_end   = Ap[i+1];
        ValueType sum = y[i];
        for (IndexType jj = row_start; jj < row_end; jj++) {            
            const IndexType j = Aj[jj];  //column index
            sum += x[j] * Ax[jj];
        }
        y[i] = sum; 
    }
}

template <typename IndexType, typename ValueType>
void spmv_csr_serial_host(const csr_matrix<IndexType, ValueType>& csr, 
                          const ValueType * x,  
                                ValueType * y)
{
    __spmv_csr_serial_host(csr.num_rows, csr.Ap, csr.Aj, csr.Ax, x, y);
}


////////////////////////////////////////////////////////////////////////////////
//! Compute y += A*x for a sparse COO matrix A and column vectors x and y
//! @param num_nonzeros   number of nonzeros in A
//! @param rows           COO row indices array
//! @param cols           COO column index array
//! @param data           COO data array
//! @param x              column vector
//! @param y              column vector
////////////////////////////////////////////////////////////////////////////////
template <typename IndexType, typename ValueType>
void __spmv_coo_serial_host(const IndexType num_nonzeros,
                            const IndexType * rows, 
                            const IndexType * cols, 
                            const ValueType * data, 
                            const ValueType * x,  
                                  ValueType * y)
{
    for (IndexType i = 0; i < num_nonzeros; i++){   
        y[rows[i]] += data[i] * x[cols[i]];
    }
}

template <typename IndexType, typename ValueType>
void spmv_coo_serial_host(const coo_matrix<IndexType, ValueType>& coo, 
                          const ValueType * x,  
                                ValueType * y)
{
    __spmv_coo_serial_host(coo.num_nonzeros, coo.I, coo.J, coo.V, x, y);
}


////////////////////////////////////////////////////////////////////////////////
//! Compute y += A*x for a hybrid ELL/COO matrix A and column vectors x and y
//! @param hyb        hyb_matrix
//! @param x          column vector
//! @param y          column vector
////////////////////////////////////////////////////////////////////////////////
template <typename IndexType, typename ValueType>
void spmv_hyb_serial_host(const hyb_matrix<IndexType,ValueType>& hyb, 
                          const ValueType * x, 
                                ValueType * y)
{
    spmv_ell_serial_host(hyb.ell, x, y);
    spmv_coo_serial_host(hyb.coo, x, y);
}


////////////////////////////////////////////////////////////////////////////////
//! Compute y += A*x for a Packet Matrix matrix A and column vectors x and y
//! @param pkt        pkt_matrix
//! @param x          column vector
//! @param y          column vector
////////////////////////////////////////////////////////////////////////////////
template <typename IndexType, typename ValueType>
void spmv_pkt_serial_host(const pkt_matrix<IndexType,ValueType>& pkt, 
                          const ValueType * x, 
                                ValueType * y)
{
    //initialize temp arrays for permuted x and y
    ValueType * x_perm = new_host_array<ValueType>(pkt.num_cols); 
    ValueType * y_perm = new_host_array<ValueType>(pkt.num_rows); 
    
    for(IndexType i = 0; i < pkt.num_rows; i++){
        //permute x and y    
        x_perm[i] = x[pkt.permute_new_to_old[i]];
        y_perm[i] = y[pkt.permute_new_to_old[i]];
    }

    for(IndexType N = 0; N < pkt.num_packets; N++){
        IndexType base_row = pkt.row_ptr[N];

        for(IndexType H = 0; H < pkt.threads_per_packet; H++){
            IndexType thread_id = H + N*pkt.threads_per_packet;

            IndexType pos = pkt.packets.pos_start[thread_id];
            IndexType end = pkt.packets.pos_end[thread_id];
            
            while(pos != end){
                const IndexType packed_index = pkt.packets.index_array[pos];

                const IndexType row = base_row + pkt_unpack_row_index(packed_index);
                const IndexType col = base_row + pkt_unpack_col_index(packed_index);
                const ValueType val = pkt.packets.data_array[pos];          

                y_perm[row] += val*x_perm[col];         

                pos    += pkt.threads_per_packet;
            }

        }
    }

    // remaining entries are stored in COO format
    spmv_coo_serial_host<IndexType,ValueType>(pkt.coo, x_perm, y_perm);
   
    for(IndexType i = 0; i < pkt.num_rows; i++){
        //unpermute y    
        y[i] = y_perm[pkt.permute_old_to_new[i]];
    }
    

    //delete temp arrays
    delete_host_array(x_perm);
    delete_host_array(y_perm);
}

