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
#include "timer.h"

char * BENCHMARK_OUTPUT_FILE_NAME = "benchmark_output.log";
   
template <typename IndexType, typename ValueType>
size_t bytes_per_spmv(const dia_matrix<IndexType,ValueType>& mtx)
{
    // note: this neglects diag_offsets, which is < 1% of other parts
    size_t bytes = 0;
    bytes += 2*sizeof(ValueType) * mtx.num_nonzeros; // A[i,j] and x[j]
    bytes += 2*sizeof(ValueType) * mtx.num_rows;     // y[i] = y[i] + ...
    return bytes;
}

template <typename IndexType, typename ValueType>
size_t bytes_per_spmv(const ell_matrix<IndexType,ValueType>& mtx)
{
    size_t bytes = 0;
    bytes += 1*sizeof(IndexType) * mtx.num_nonzeros; // column index
    bytes += 1*sizeof(ValueType) * mtx.stride * mtx.num_cols_per_row; // A[i,j] and padding
    bytes += 1*sizeof(ValueType) * mtx.num_nonzeros; // x[j]
    bytes += 2*sizeof(ValueType) * mtx.num_rows;     // y[i] = y[i] + ...
    return bytes;
}

template <typename IndexType, typename ValueType>
size_t bytes_per_spmv(const csr_matrix<IndexType,ValueType>& mtx)
{
    size_t bytes = 0;
    bytes += 2*sizeof(IndexType) * mtx.num_rows;     // row pointer
    bytes += 1*sizeof(IndexType) * mtx.num_nonzeros; // column index
    bytes += 2*sizeof(ValueType) * mtx.num_nonzeros; // A[i,j] and x[j]
    bytes += 2*sizeof(ValueType) * mtx.num_rows;     // y[i] = y[i] + ...
    return bytes;
}

template <typename IndexType, typename ValueType>
size_t bytes_per_spmv(const coo_matrix<IndexType,ValueType>& mtx)
{
    size_t bytes = 0;
    bytes += 2*sizeof(IndexType) * mtx.num_nonzeros; // row and column indices
    bytes += 2*sizeof(ValueType) * mtx.num_nonzeros; // A[i,j] and x[j]

    std::vector<size_t> occupied_rows(mtx.num_rows, 0);
    for(size_t n = 0; n < mtx.num_nonzeros; n++)
        occupied_rows[mtx.I[n]] = 1;
    for(size_t n = 0; n < mtx.num_rows; n++)
        if(occupied_rows[n] == 1)
            bytes += 2*sizeof(ValueType);            // y[i] = y[i] + ...
    return bytes;
}

template <typename IndexType, typename ValueType>
size_t bytes_per_spmv(const hyb_matrix<IndexType,ValueType>& mtx)
{
    return bytes_per_spmv(mtx.ell) + bytes_per_spmv(mtx.coo);
}

template <typename IndexType, typename ValueType>
size_t bytes_per_spmv(const pkt_matrix<IndexType,ValueType>& mtx)
{
    size_t bytes = 0;
    bytes +=   sizeof(PackedIndexType) * mtx.num_nonzeros; // row and column indices
    bytes += 1*sizeof(ValueType)       * mtx.num_nonzeros; // A[i,j] 
    bytes += 1*sizeof(ValueType)       * mtx.num_cols;     // x[j]
    bytes += 2*sizeof(ValueType)       * mtx.num_rows;     // y[i] = y[i] + ...
    bytes += 2*sizeof(IndexType)       * mtx.threads_per_packet * mtx.num_packets;  // pos_start, pos_end
    bytes += 2*sizeof(IndexType)       * mtx.threads_per_packet * mtx.num_packets;  // row_ptr
    bytes += bytes_per_spmv(mtx.coo);
    return bytes;
}


template <typename SparseMatrix, typename SpMV>
double benchmark_spmv(const SparseMatrix & sp_host, SpMV spmv, const memory_location loc, const char * method_name, const size_t min_iterations = 1, const size_t max_iterations = 1000, const double seconds = 3.0)
{
    SparseMatrix sp_loc = (loc == DEVICE_MEMORY) ? copy_matrix_to_device(sp_host) : sp_host;

    typedef typename SparseMatrix::value_type ValueType;
    typedef typename SparseMatrix::index_type IndexType;

    //initialize host arrays
    ValueType * x_host = new_host_array<ValueType>(sp_host.num_cols);
    ValueType * y_host = new_host_array<ValueType>(sp_host.num_rows);
    for(IndexType i = 0; i < sp_host.num_cols; i++)
        x_host[i] = rand() / (RAND_MAX + 1.0); 
    std::fill(y_host, y_host + sp_host.num_rows, 0);
   
    //initialize device arrays
    ValueType * x_loc = copy_array(x_host, sp_host.num_cols, HOST_MEMORY, loc);
    ValueType * y_loc = copy_array(y_host, sp_host.num_rows, HOST_MEMORY, loc);
        
    // warmup    
    timer time_one_iteration;
    spmv(sp_loc, x_loc, y_loc);
    cudaThreadSynchronize();
    double estimated_time = time_one_iteration.seconds_elapsed();

    //printf("estimated time %f\n", (float) estimated_time);

    // determine # of seconds dynamically
    size_t num_iterations;
    if (estimated_time == 0)
        num_iterations = max_iterations;
    else
        num_iterations = std::min(max_iterations, std::max(min_iterations, (size_t) (seconds / estimated_time)) ); 
   
    //printf("performing %d iterations\n", num_iterations);

    // time several SpMV iterations
    timer t;
    for(size_t i = 0; i < num_iterations; i++)
        spmv(sp_loc, x_loc, y_loc);
    cudaThreadSynchronize();
    double msec_per_iteration = t.milliseconds_elapsed() / (double) num_iterations;
    double sec_per_iteration = msec_per_iteration / 1000.0;
    double GFLOPs = (sec_per_iteration == 0) ? 0 : (2.0 * (double) sp_host.num_nonzeros / sec_per_iteration) / 1e9;
    double GBYTEs = (sec_per_iteration == 0) ? 0 : ((double) bytes_per_spmv(sp_host) / sec_per_iteration) / 1e9;

    const char * location = (loc == HOST_MEMORY) ? "cpu" : "gpu";
    printf("\tbenchmarking %-20s [%s]: %8.4f ms ( %5.2f GFLOP/s %5.1f GB/s)\n", \
            method_name, location, msec_per_iteration, GFLOPs, GBYTEs); 

    //record results to file
    FILE * fid = fopen(BENCHMARK_OUTPUT_FILE_NAME, "a");
    fprintf(fid, "kernel=%s gflops=%f gbytes=%f msec=%f\n", method_name, GFLOPs, GBYTEs, msec_per_iteration); 
    fclose(fid);

    //deallocate buffers
    delete_host_array(x_host);
    delete_host_array(y_host);
    delete_array(x_loc, loc);
    delete_array(y_loc, loc);

    if (loc == DEVICE_MEMORY) delete_device_matrix(sp_loc);

    return msec_per_iteration;
}

template <typename SparseMatrix, typename SpMV>
double benchmark_spmv_on_device(const SparseMatrix & sp_host, SpMV spmv, const char * method_name = NULL)
{
    return benchmark_spmv<SparseMatrix, SpMV>(sp_host, spmv, DEVICE_MEMORY, method_name);
}

template <typename SparseMatrix, typename SpMV>
double benchmark_spmv_on_host(const SparseMatrix & sp_host, SpMV spmv, const char * method_name = NULL, const size_t num_iterations = 500)
{
    return benchmark_spmv<SparseMatrix, SpMV>(sp_host, spmv, HOST_MEMORY, method_name);
}

