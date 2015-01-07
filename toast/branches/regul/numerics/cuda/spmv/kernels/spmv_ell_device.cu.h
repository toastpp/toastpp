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

#include "sparse_formats.h"
#include "utils.h"
#include "texture.h"

// SpMV kernel for the ELLPACK/ITPACK matrix format.

template <typename IndexType, typename ValueType, bool UseCache>
__global__ void
spmv_ell_kernel(const IndexType num_rows, 
                const IndexType num_cols, 
                const IndexType num_cols_per_row,
                const IndexType stride,
                const IndexType * Aj,
                const ValueType * Ax, 
                const ValueType * x, 
                      ValueType * y)
{
    const IndexType row = large_grid_thread_id();

    if(row >= num_rows){ return; }

    ValueType sum = y[row];

    Aj += row;
    Ax += row;

    for(IndexType n = 0; n < num_cols_per_row; n++){
        const ValueType A_ij = *Ax;

        if (A_ij != 0){
            const IndexType col = *Aj;
            sum += A_ij * fetch_x<UseCache>(col, x);
        }

        Aj += stride;
        Ax += stride;
    }

    y[row] = sum;
}

template <typename IndexType, typename ValueType>
void spmv_ell_device(const ell_matrix<IndexType,ValueType>& d_ell, 
                     const ValueType * d_x, 
                           ValueType * d_y)
{
    const unsigned int BLOCK_SIZE = 256;
    const dim3 grid = make_large_grid(d_ell.num_rows, BLOCK_SIZE);
   
    spmv_ell_kernel<IndexType,ValueType,false> <<<grid, BLOCK_SIZE>>>
        (d_ell.num_rows, d_ell.num_cols, d_ell.num_cols_per_row, d_ell.stride,
         d_ell.Aj, d_ell.Ax,
         d_x, d_y);
}

template <typename IndexType, typename ValueType>
void spmv_ell_tex_device(const ell_matrix<IndexType,ValueType>& d_ell, 
                         const ValueType * d_x, 
                               ValueType * d_y)
{
    const unsigned int BLOCK_SIZE = 256;
    const dim3 grid = make_large_grid(d_ell.num_rows, BLOCK_SIZE);
  
    bind_x(d_x);

    spmv_ell_kernel<IndexType,ValueType,true> <<<grid, BLOCK_SIZE>>>
        (d_ell.num_rows, d_ell.num_cols, d_ell.num_cols_per_row, d_ell.stride,
         d_ell.Aj, d_ell.Ax,
         d_x, d_y);

    unbind_x(d_x);
}

