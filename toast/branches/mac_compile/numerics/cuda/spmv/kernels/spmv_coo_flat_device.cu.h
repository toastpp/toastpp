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
#include "kernels/spmv_coo_serial_device.cu.h"
#include "kernels/spmv_common_device.cu.h"


// segmented reduction in shared memory
template <typename IndexType, typename ValueType>
__device__ ValueType segreduce_warp(const IndexType thread_lane, IndexType row, ValueType val, IndexType * rows, ValueType * vals)
{
    rows[threadIdx.x] = row;
    vals[threadIdx.x] = val;

    if( thread_lane >=  1 && row == rows[threadIdx.x -  1] ) { vals[threadIdx.x] = val = val + vals[threadIdx.x -  1]; } 
    if( thread_lane >=  2 && row == rows[threadIdx.x -  2] ) { vals[threadIdx.x] = val = val + vals[threadIdx.x -  2]; }
    if( thread_lane >=  4 && row == rows[threadIdx.x -  4] ) { vals[threadIdx.x] = val = val + vals[threadIdx.x -  4]; }
    if( thread_lane >=  8 && row == rows[threadIdx.x -  8] ) { vals[threadIdx.x] = val = val + vals[threadIdx.x -  8]; }
    if( thread_lane >= 16 && row == rows[threadIdx.x - 16] ) { vals[threadIdx.x] = val = val + vals[threadIdx.x - 16]; }

    return val;
}

template <typename IndexType, typename ValueType>
__device__ void segreduce_block(const IndexType * idx, ValueType * val)
{
    ValueType left = 0;
    if( threadIdx.x >=   1 && idx[threadIdx.x] == idx[threadIdx.x -   1] ) { left = val[threadIdx.x -   1]; } __syncthreads(); val[threadIdx.x] += left; left = 0; __syncthreads();  
    if( threadIdx.x >=   2 && idx[threadIdx.x] == idx[threadIdx.x -   2] ) { left = val[threadIdx.x -   2]; } __syncthreads(); val[threadIdx.x] += left; left = 0; __syncthreads();
    if( threadIdx.x >=   4 && idx[threadIdx.x] == idx[threadIdx.x -   4] ) { left = val[threadIdx.x -   4]; } __syncthreads(); val[threadIdx.x] += left; left = 0; __syncthreads();
    if( threadIdx.x >=   8 && idx[threadIdx.x] == idx[threadIdx.x -   8] ) { left = val[threadIdx.x -   8]; } __syncthreads(); val[threadIdx.x] += left; left = 0; __syncthreads();
    if( threadIdx.x >=  16 && idx[threadIdx.x] == idx[threadIdx.x -  16] ) { left = val[threadIdx.x -  16]; } __syncthreads(); val[threadIdx.x] += left; left = 0; __syncthreads();
    if( threadIdx.x >=  32 && idx[threadIdx.x] == idx[threadIdx.x -  32] ) { left = val[threadIdx.x -  32]; } __syncthreads(); val[threadIdx.x] += left; left = 0; __syncthreads();  
    if( threadIdx.x >=  64 && idx[threadIdx.x] == idx[threadIdx.x -  64] ) { left = val[threadIdx.x -  64]; } __syncthreads(); val[threadIdx.x] += left; left = 0; __syncthreads();
    if( threadIdx.x >= 128 && idx[threadIdx.x] == idx[threadIdx.x - 128] ) { left = val[threadIdx.x - 128]; } __syncthreads(); val[threadIdx.x] += left; left = 0; __syncthreads();
    if( threadIdx.x >= 256 && idx[threadIdx.x] == idx[threadIdx.x - 256] ) { left = val[threadIdx.x - 256]; } __syncthreads(); val[threadIdx.x] += left; left = 0; __syncthreads();
}

//////////////////////////////////////////////////////////////////////////////
// COO SpMV kernel which flattens data irregularity (segmented reduction)
//////////////////////////////////////////////////////////////////////////////
//
// spmv_coo_flat_device
//   The input coo_matrix must be sorted by row.  Columns within each row
//   may appear in any order and duplicate entries are also acceptable.
//   This sorted COO format is easily obtained by expanding the row pointer
//   of a CSR matrix (csr.Ap) into proper row indices and then copying 
//   the arrays containing the CSR column indices (csr.Aj) and nonzero values
//   (csr.Ax) verbatim.  A segmented reduction is used to compute the per-row
//   sums.
//
// spmv_coo_flat_tex_device
//   Same as spmv_coo_flat_device, except that the texture cache is 
//   used for accessing the x vector.
//

template <typename IndexType, typename ValueType, unsigned int BLOCK_SIZE, bool UseCache>
__global__ void
spmv_coo_flat_kernel(const IndexType num_nonzeros,
                     const IndexType interval_size,
                     const IndexType * I, 
                     const IndexType * J, 
                     const ValueType * V, 
                     const ValueType * x, 
                           ValueType * y,
                           IndexType * temp_rows,
                           ValueType * temp_vals)
{
    __shared__ IndexType rows[BLOCK_SIZE];
    __shared__ ValueType vals[BLOCK_SIZE];

    const IndexType thread_id   = BLOCK_SIZE * blockIdx.x + threadIdx.x;                 // global thread index
    const IndexType thread_lane = threadIdx.x & (WARP_SIZE-1);                           // thread index within the warp
    const IndexType warp_id     = thread_id   / WARP_SIZE;                               // global warp index

    const IndexType interval_begin = warp_id * interval_size;                            // warp's offset into I,J,V
    const IndexType interval_end   = min(interval_begin + interval_size, num_nonzeros);  // end of warps's work

    if(interval_begin >= interval_end)                                                   // warp has no work to do 
        return;

    if (thread_lane == 31){
        // initialize the carry in values
        rows[threadIdx.x] = I[interval_begin]; 
        vals[threadIdx.x] = 0;
    }
  
    for(IndexType n = interval_begin + thread_lane; n < interval_end; n += WARP_SIZE)
    {
        IndexType row = I[n];                                         // row index (i)
        ValueType val = V[n] * fetch_x<UseCache>(J[n], x);            // A(i,j) * x(j)
        
        if (thread_lane == 0)
        {
            if(row == rows[threadIdx.x + 31])
                val += vals[threadIdx.x + 31];                        // row continues
            else
                y[rows[threadIdx.x + 31]] += vals[threadIdx.x + 31];  // row terminated
        }
        
        val = segreduce_warp(thread_lane, row, val, rows, vals);      // segmented reduction in shared memory

        if(thread_lane < 31 && row != rows[threadIdx.x + 1])
            y[row] += val;                                            // row terminated
    }

    if(thread_lane == 31)
    {
        // write the carry out values
        temp_rows[warp_id] = rows[threadIdx.x];
        temp_vals[warp_id] = vals[threadIdx.x];
    }
}

// The second level of the segmented reduction operation
template <typename IndexType, typename ValueType, unsigned int BLOCK_SIZE>
__global__ void
spmv_coo_reduce_update_kernel(const IndexType num_warps,
                              const IndexType * temp_rows,
                              const ValueType * temp_vals,
                                    ValueType * y)
{
    __shared__ IndexType rows[BLOCK_SIZE + 1];    
    __shared__ ValueType vals[BLOCK_SIZE + 1];    

    const IndexType end = num_warps - (num_warps & (BLOCK_SIZE - 1));

    if (threadIdx.x == 0){
        rows[BLOCK_SIZE] = (IndexType) -1;
        vals[BLOCK_SIZE] = (ValueType)  0;
    }
    
    __syncthreads();

    IndexType i = threadIdx.x;

    while (i < end){
        // do full blocks
        rows[threadIdx.x] = temp_rows[i];
        vals[threadIdx.x] = temp_vals[i];

        __syncthreads();

        segreduce_block(rows, vals);

        if (rows[threadIdx.x] != rows[threadIdx.x + 1])
            y[rows[threadIdx.x]] += vals[threadIdx.x];

        __syncthreads();

        i += BLOCK_SIZE; 
    }

    if (end < num_warps){
        if (i < num_warps){
            rows[threadIdx.x] = temp_rows[i];
            vals[threadIdx.x] = temp_vals[i];
        } else {
            rows[threadIdx.x] = (IndexType) -1;
            vals[threadIdx.x] = (ValueType)  0;
        }

        __syncthreads();
   
        segreduce_block(rows, vals);

        if (i < num_warps)
            if (rows[threadIdx.x] != rows[threadIdx.x + 1])
                y[rows[threadIdx.x]] += vals[threadIdx.x];
    }
}


template <typename IndexType, typename ValueType, bool UseCache>
void __spmv_coo_flat_device(const coo_matrix<IndexType,ValueType>& d_coo, 
                            const ValueType * d_x, 
                                  ValueType * d_y)
{
    if(d_coo.num_nonzeros == 0)
    {
        // empty matrix
        return;
    }
    else if (d_coo.num_nonzeros < WARP_SIZE)
    {
        // small matrix
        spmv_coo_serial_kernel<IndexType,ValueType> <<<1,1>>>
            (d_coo.num_nonzeros, d_coo.I, d_coo.J, d_coo.V, d_x, d_y);
        return;
    }

    //TODO Determine optimal BLOCK_SIZE and MAX_BLOCKS
    const unsigned int BLOCK_SIZE      = 256;
    const unsigned int MAX_BLOCKS      = MAX_THREADS / (2 * BLOCK_SIZE);
    const unsigned int WARPS_PER_BLOCK = BLOCK_SIZE / WARP_SIZE;

    const unsigned int num_units  = d_coo.num_nonzeros / WARP_SIZE; 
    const unsigned int num_warps  = std::min(num_units, WARPS_PER_BLOCK * MAX_BLOCKS);
    const unsigned int num_blocks = DIVIDE_INTO(num_warps, WARPS_PER_BLOCK);
    const unsigned int num_iters  = DIVIDE_INTO(num_units, num_warps);
    
    const unsigned int interval_size = WARP_SIZE * num_iters;

    const IndexType tail = num_units * WARP_SIZE; // do the last few nonzeros separately (fewer than WARP_SIZE elements)

    const unsigned int active_warps = (interval_size == 0) ? 0 : DIVIDE_INTO(tail, interval_size);

    if (UseCache)
        bind_x(d_x);

    IndexType * temp_rows = new_device_array<IndexType>(active_warps);
    ValueType * temp_vals = new_device_array<ValueType>(active_warps);

    spmv_coo_flat_kernel<IndexType, ValueType, BLOCK_SIZE, UseCache> <<<num_blocks, BLOCK_SIZE>>>
        (tail, interval_size, d_coo.I, d_coo.J, d_coo.V, d_x, d_y, temp_rows, temp_vals);

    spmv_coo_serial_kernel<IndexType,ValueType> <<<1,1>>>
        (d_coo.num_nonzeros - tail, d_coo.I + tail, d_coo.J + tail, d_coo.V + tail, d_x, d_y);

    spmv_coo_reduce_update_kernel<IndexType, ValueType, 512> <<<1, 512>>> (active_warps, temp_rows, temp_vals, d_y);

    delete_device_array(temp_rows);
    delete_device_array(temp_vals);

    if (UseCache)
        unbind_x(d_x);
}

template <typename IndexType, typename ValueType>
void spmv_coo_flat_device(const coo_matrix<IndexType,ValueType>& d_coo, 
                          const ValueType * d_x, 
                                ValueType * d_y)
{ 
    __spmv_coo_flat_device<IndexType, ValueType, false>(d_coo, d_x, d_y);
}


template <typename IndexType, typename ValueType>
void spmv_coo_flat_tex_device(const coo_matrix<IndexType,ValueType>& d_coo, 
                              const ValueType * d_x, 
                                    ValueType * d_y)
{ 
    __spmv_coo_flat_device<IndexType, ValueType, true>(d_coo, d_x, d_y);
}

