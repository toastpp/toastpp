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

#include <stdio.h>

#include "sparse_formats.h"
#include "utils.h"
#include "texture.h"

////////////////////////////////////////////////////////////////////////
// CSR SpMV kernels based on a scalar model (one thread per row)
///////////////////////////////////////////////////////////////////////
//
// spmv_csr_scalar_device
//   Straightforward translation of standard CSR SpMV to CUDA
//   where each thread computes y[i] += A[i,:] * x 
//   (the dot product of the i-th row of A with the x vector)
//
// spmv_csr_scalar_tex_device
//   Same as spmv_csr_scalar_device, except x is accessed via texture cache.
//

template <typename IndexType, typename ValueType, bool UseCache>
__global__ void
spmv_csr_scalar_kernel(const IndexType num_rows,
                       const IndexType * Ap, 
                       const IndexType * Aj, 
                       const ValueType * Ax, 
                       const ValueType * x, 
                             ValueType * y)
{

    // row index
    const IndexType row = large_grid_thread_id();
    
    if(row < num_rows){     
        ValueType sum = y[row];

        const IndexType row_start = Ap[row];
        const IndexType row_end   = Ap[row+1];
    
        for (IndexType jj = row_start; jj < row_end; jj++){             
            sum += Ax[jj] * fetch_x<UseCache>(Aj[jj], x);       
        }

        y[row] = sum;
    }
}

template <typename IndexType, typename ValueType>
void spmv_csr_scalar_device(const csr_matrix<IndexType,ValueType>& d_csr, 
                            const ValueType * d_x, 
                                  ValueType * d_y)
{
    const unsigned int BLOCK_SIZE = 256; 
    
    const dim3 grid = make_large_grid(d_csr.num_rows, BLOCK_SIZE);

    spmv_csr_scalar_kernel<IndexType, ValueType, false> <<<grid, BLOCK_SIZE>>> 
        (d_csr.num_rows, d_csr.Ap, d_csr.Aj, d_csr.Ax, d_x, d_y);   
}

template <typename IndexType, typename ValueType>
void spmv_csr_scalar_tex_device(const csr_matrix<IndexType,ValueType>& d_csr, 
                                const ValueType * d_x, 
                                      ValueType * d_y)
{
    const unsigned int BLOCK_SIZE = 256;
    
    const dim3 grid = make_large_grid(d_csr.num_rows, BLOCK_SIZE);
    
    bind_x(d_x);

    spmv_csr_scalar_kernel<IndexType, ValueType, true> <<<grid, BLOCK_SIZE>>> 
        (d_csr.num_rows, d_csr.Ap, d_csr.Aj, d_csr.Ax, d_x, d_y);   

    unbind_x(d_x);
}

