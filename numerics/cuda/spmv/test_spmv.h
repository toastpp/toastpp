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

// Functions to test SpMV kernels

#include <algorithm>
#include <limits>
#include <cmath>
#include "mem.h"

template <typename T>
T maximum_relative_error(const T * A, const T * B, const size_t N)
{
    T max_error = 0;

    T eps = std::sqrt( std::numeric_limits<T>::epsilon() );

    for(size_t i = 0; i < N; i++)
    {
        const T a = A[i];
        const T b = B[i];
        const T error = std::abs(a - b);
        if (error != 0){
            max_error = std::max(max_error, error/(std::abs(a) + std::abs(b) + eps) );
        }
    }
 
    return max_error;
}


// Compare two SpMV kernels
template <typename SparseMatrix1, typename SpMV1,
          typename SparseMatrix2, typename SpMV2>
void compare_spmv_kernels(const SparseMatrix1 & sm1_host, SpMV1 spmv1, const memory_location loc1,
                          const SparseMatrix2 & sm2_host, SpMV2 spmv2, const memory_location loc2)
{
    // sanity checking
    assert(sm1_host.num_rows == sm2_host.num_rows);
    assert(sm1_host.num_cols == sm2_host.num_cols);
    assert(sm1_host.num_nonzeros == sm2_host.num_nonzeros);

    typedef typename SparseMatrix1::index_type IndexType;
    typedef typename SparseMatrix2::value_type ValueType;

    const IndexType num_rows = sm1_host.num_rows;
    const IndexType num_cols = sm1_host.num_cols;

    // transfer matrices from host to destination location
    SparseMatrix1 sm1_loc1 = (loc1 == DEVICE_MEMORY) ? copy_matrix_to_device(sm1_host) : sm1_host;
    SparseMatrix2 sm2_loc2 = (loc2 == DEVICE_MEMORY) ? copy_matrix_to_device(sm2_host) : sm2_host;
    
    // initialize host vectors
    ValueType * x_host = new_host_array<ValueType>(num_cols);
    ValueType * y_host = new_host_array<ValueType>(num_rows);
    
    for(IndexType i = 0; i < num_cols; i++)
        x_host[i] = rand() / (RAND_MAX + 1.0); 
    for(IndexType i = 0; i < num_rows; i++)
        y_host[i] = rand() / (RAND_MAX + 1.0);

    // create vectors in appropriate locations
    ValueType * x_loc1 = copy_array(x_host, num_cols, HOST_MEMORY, loc1);
    ValueType * y_loc1 = copy_array(y_host, num_rows, HOST_MEMORY, loc1);
    ValueType * x_loc2 = copy_array(x_host, num_cols, HOST_MEMORY, loc2);
    ValueType * y_loc2 = copy_array(y_host, num_rows, HOST_MEMORY, loc2);
    
    // compute y = A*x
    spmv1(sm1_loc1, x_loc1, y_loc1);
    spmv2(sm2_loc2, x_loc2, y_loc2);
   
    // transfer results to host
    ValueType * y_sm1_result = copy_array(y_loc1, num_rows, loc1, HOST_MEMORY);
    ValueType * y_sm2_result = copy_array(y_loc2, num_rows, loc2, HOST_MEMORY);

    ValueType max_error = maximum_relative_error(y_sm1_result, y_sm2_result, num_rows);
    printf(" [max error %9f]", max_error);
    
    if ( max_error > 5 * std::sqrt( std::numeric_limits<ValueType>::epsilon() ) )
        printf(" POSSIBLE FAILURE");
               
    // cleanup
    if (loc1 == DEVICE_MEMORY) delete_device_matrix(sm1_loc1);
    if (loc2 == DEVICE_MEMORY) delete_device_matrix(sm2_loc2);
    delete_host_array(x_host);
    delete_host_array(y_host);
    delete_array(x_loc1, loc1);
    delete_array(y_loc1, loc1);
    delete_array(x_loc2, loc2);
    delete_array(y_loc2, loc2);
    delete_host_array(y_sm1_result);
    delete_host_array(y_sm2_result);
}


// SpMV1 is the reference
template <typename SparseMatrix1, typename SpMV1,
          typename SparseMatrix2, typename SpMV2>
void test_spmv_kernel(const SparseMatrix1 & sm1_host, SpMV1 spmv1, const memory_location loc1,
                      const SparseMatrix2 & sm2_host, SpMV2 spmv2, const memory_location loc2,
                      const char * method_name)
{
    printf("\ttesting %-26s", method_name);
    if(loc2 == HOST_MEMORY)
        printf("[cpu]:");
    else
        printf("[gpu]:");
    
    compare_spmv_kernels( sm1_host, spmv1, loc1, sm2_host, spmv2, loc2);

    printf("\n");
}

