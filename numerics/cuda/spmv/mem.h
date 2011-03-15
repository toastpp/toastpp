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

#include <stdlib.h>

#ifdef __CUDACC__
#include <cuda.h>
#  define CUDA_SAFE_CALL_NO_SYNC( call) do {                                 \
    cudaError err = call;                                                    \
    if( cudaSuccess != err) {                                                \
        fprintf(stderr, "Cuda error in file '%s' in line %i : %s.\n",        \
                __FILE__, __LINE__, cudaGetErrorString( err) );              \
        exit(EXIT_FAILURE);                                                  \
    } } while (0)

#  define CUDA_SAFE_CALL( call) do {                                         \
    CUDA_SAFE_CALL_NO_SYNC(call);                                            \
    cudaError err = cudaThreadSynchronize();                                 \
    if( cudaSuccess != err) {                                                \
        fprintf(stderr, "Cuda error in file '%s' in line %i : %s.\n",        \
                __FILE__, __LINE__, cudaGetErrorString( err) );              \
        exit(EXIT_FAILURE);                                                  \
    } } while (0)

#else
// drop all CUDA calls
#define CUDA_SAFE_CALL_NO_SYNC(x)
#define CUDA_SAFE_CALL(x)
#endif



/////////////////////////////////////////////////////////////////////
// allocate memory on host and device
/////////////////////////////////////////////////////////////////////

enum memory_location { HOST_MEMORY, DEVICE_MEMORY };

template <typename T>
T * new_array(const size_t N, const memory_location loc) 
{ 
    //dispatch on location
    if (loc == HOST_MEMORY){
        return (T*) malloc(N * sizeof(T));
    }
    else{
        void * ptr = 0;
        CUDA_SAFE_CALL(cudaMalloc(&ptr, sizeof(T)*N));
        return static_cast<T*>(ptr);
    }
}

template <typename T>
void delete_array(T * p, const memory_location loc) 
{ 
    //dispatch on location
    if (loc == HOST_MEMORY){
        free(p);
    }
    else {
        CUDA_SAFE_CALL(cudaFree(p));
    };
}



template<typename T>  T * new_host_array(const size_t N)   { return new_array<T>(N, HOST_MEMORY); }
template<typename T>  T * new_device_array(const size_t N) { return new_array<T>(N, DEVICE_MEMORY); }
template<typename T>  void delete_host_array(T *p)   { delete_array(p, HOST_MEMORY); }
template<typename T>  void delete_device_array(T *p) { delete_array(p, DEVICE_MEMORY); }



/////////////////////////////////////////////////////////////////////
// transfer data between host and device
/////////////////////////////////////////////////////////////////////

template<typename T>
void memcpy_array(T * dst, const T * src, const size_t N, 
                  const memory_location src_loc,
                  const memory_location dst_loc)
{
    if(src_loc == HOST_MEMORY   && dst_loc == HOST_MEMORY  )
        CUDA_SAFE_CALL(cudaMemcpy(dst, src, sizeof(T)*N, cudaMemcpyHostToHost));
    if(src_loc == HOST_MEMORY   && dst_loc == DEVICE_MEMORY)
        CUDA_SAFE_CALL(cudaMemcpy(dst, src, sizeof(T)*N, cudaMemcpyHostToDevice));
    if(src_loc == DEVICE_MEMORY && dst_loc == HOST_MEMORY  )
        CUDA_SAFE_CALL(cudaMemcpy(dst, src, sizeof(T)*N, cudaMemcpyDeviceToHost));
    if(src_loc == DEVICE_MEMORY && dst_loc == DEVICE_MEMORY)
        CUDA_SAFE_CALL(cudaMemcpy(dst, src, sizeof(T)*N, cudaMemcpyDeviceToDevice));
}

template<typename T>
 void memcpy_to_device(T *dp, const T *hp, const size_t N)
{
    memcpy_array(dp, hp, N, HOST_MEMORY, DEVICE_MEMORY);
}

template<typename T>
 void memcpy_to_host(T *hp, const T *dp, const size_t N)
{
    memcpy_array(hp, dp, N, DEVICE_MEMORY, HOST_MEMORY);
}

template<typename T>
 void memcpy_on_device(T *dp2, const T *dp1, const size_t N)
{
    memcpy_array(dp2, dp1, N, DEVICE_MEMORY, DEVICE_MEMORY);
}

template<typename T>
 void memcpy_on_host(T *hp2, const T *hp1, const size_t N)
{
    memcpy_array(hp2, hp1, N, HOST_MEMORY, HOST_MEMORY);
}

template<typename T, typename S>
 void memcpy_to_symbol(const S& s, const T *hp, const size_t N)
{
    cudaMemcpyToSymbol(s, hp, sizeof(T)*N);
}

template<typename T, typename S>
 void memcpy_from_symbol(T *hp, const S& s, const size_t N)
{
    cudaMemcpyFromSymbol(hp, s, sizeof(T)*N);
}

/////////////////////////////////////////////////////////////////////
// allocate and transfer data
/////////////////////////////////////////////////////////////////////
template <typename T>
T * copy_array_to_host(const T * d_ptr, const size_t N)
{
    T * h_ptr = new_host_array<T>(N);
    memcpy_to_host(h_ptr, d_ptr, N);
    return h_ptr;
}

template <typename T>
T * copy_array_to_device(const T * h_ptr, const size_t N)
{
    T * d_ptr = new_device_array<T>(N);
    memcpy_to_device(d_ptr, h_ptr, N);
    return d_ptr;
}

template <typename T>
T * copy_array_on_host(const T * h_src, const size_t N)
{
    T * h_dst = new_host_array<T>(N);
    memcpy_on_host(h_dst, h_src, N);
    return h_dst;
}

template <typename T>
T * copy_array_on_device(const T * d_src, const size_t N)
{
    T * d_dst = new_device_array<T>(N);
    memcpy_on_device(d_dst, d_src, N);
    return d_dst;
}

template <typename T>
T * copy_array(const T * src, const size_t N, 
               const memory_location src_loc, 
               const memory_location dst_loc)
{
    T * dst = new_array<T>(N, dst_loc);
    memcpy_array(dst, src, N, src_loc, dst_loc);
    return dst;
}

