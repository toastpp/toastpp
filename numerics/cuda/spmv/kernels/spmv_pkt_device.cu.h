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
#include "kernels/spmv_coo_flat_device.cu.h"

//about half of smem
#define MAX_BYTES_PER_PACKET 8128

template <typename IndexType, typename ValueType>
__global__ void
spmv_pkt_kernel(const IndexType * row_ptr,
                const IndexType * pos_start, 
                const IndexType * pos_end, 
                const PackedIndexType * index_array, 
                const ValueType       * data_array,
                const ValueType * x, 
                      ValueType * y)
{
    __shared__ ValueType s_x[MAX_BYTES_PER_PACKET/sizeof(ValueType)];  // input x-values
    __shared__ ValueType s_y[MAX_BYTES_PER_PACKET/sizeof(ValueType)];  // output y-values

    const IndexType thread_id = small_grid_thread_id();

    // base index of the submatrix corresponding to this packet
    const IndexType packet_base_row = row_ptr[blockIdx.x]; 
    const IndexType packet_num_rows = row_ptr[blockIdx.x+1] - packet_base_row;
    
    // copy local x and y values from global memory into shared memory
    memcpy_device(s_x, x + packet_base_row, packet_num_rows);
    memcpy_device(s_y, y + packet_base_row, packet_num_rows);
    
    __syncthreads();

    ///////////////////////
    // Process Packet
    const IndexType packet_start = pos_start[thread_id];
    const IndexType packet_end   = pos_end[thread_id];

    for(IndexType pos = packet_start; pos != packet_end; pos += blockDim.x){
        //row and column indices are stored in the same 32-bit word        
        const IndexType packed_index = index_array[pos];  

        const IndexType row = pkt_unpack_row_index(packed_index);
        const IndexType col = pkt_unpack_col_index(packed_index);
        const ValueType val = data_array[pos];  

        s_y[row] += val * s_x[col]; 
    }

    __syncthreads();

    //copy y-values from shared memory to global array
    memcpy_device(y + packet_base_row, s_y, packet_num_rows);   
}



template <typename IndexType, typename ValueType>
void spmv_pkt_device(const pkt_matrix<IndexType,ValueType>& d_pkt, 
                     const ValueType * d_x, 
                           ValueType * d_y)
{
    const unsigned int BLOCK_SIZE = d_pkt.threads_per_packet;    
    const dim3 grid = make_small_grid(d_pkt.num_packets);

    assert(d_pkt.max_rows_per_packet <= MAX_BYTES_PER_PACKET/sizeof(ValueType));

    spmv_pkt_kernel <IndexType,ValueType> <<<grid,BLOCK_SIZE>>> 
        (d_pkt.row_ptr, 
         d_pkt.packets.pos_start,   d_pkt.packets.pos_end, 
         d_pkt.packets.index_array, d_pkt.packets.data_array,
         d_x, d_y);

    spmv_coo_flat_tex_device<IndexType,ValueType>(d_pkt.coo, d_x, d_y);
}






/*
 * Implements a complete SpMV, including permutation of input and output vectors
 *
 * This is only used for testing correctness, since in practice, we would
 * permute the variables after paritioning the matrix.  In an iterative solver, 
 * we would permute only once at the beginning and end of the solve, rather than 
 * at each iteration of the method.
 *
 */
template <typename IndexType, typename ValueType>
void full_spmv_pkt_device(const pkt_matrix<IndexType,ValueType>& d_pkt, 
                          const ValueType * x, 
                                ValueType * y)
{
    //initialize temp arrays for permuted x and y
    ValueType * x_perm = new_device_array<ValueType>(d_pkt.num_cols); 
    ValueType * y_perm = new_device_array<ValueType>(d_pkt.num_rows); 
        
    //permute x and y    
    gather_device<IndexType,ValueType>(x_perm, x, d_pkt.permute_new_to_old, d_pkt.num_cols);
    gather_device<IndexType,ValueType>(y_perm, y, d_pkt.permute_new_to_old, d_pkt.num_rows); // necessary if y += A*x rather than y = A*x
    
    spmv_pkt_device(d_pkt, x_perm, y_perm); 

    //unpermute result
    gather_device<IndexType,ValueType>(y, y_perm, d_pkt.permute_old_to_new, d_pkt.num_rows);

    //delete temp arrays
    delete_device_array(x_perm);
    delete_device_array(y_perm);
}

