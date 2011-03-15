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


#include <iostream>
#include <stdio.h>
#include "cmdline.h"
#include "gallery.h"
#include "tests.h"

void usage(int argc, char** argv)
{
    std::cout << "Usage:\n";
    std::cout << "\t" << argv[0] << "\n";
    std::cout << "\t" << argv[0] << " my_matrix.mtx\n";
    std::cout << "\t" << argv[0] << " my_matrix.mtx --device=1\n";
    std::cout << "\t" << argv[0] << " my_matrix.mtx --precision=64\n";
    std::cout << "\t" << argv[0] << " my_matrix.mtx --partition=somefile.txt\n\n";
    std::cout << "Note: my_matrix.mtx must be real-valued sparse matrix in the MatrixMarket file format.\n"; 
    std::cout << "      If no matrix file is provided then a simple example is created.\n";  
}


void list_devices(void)
{
    int deviceCount;
    CUDA_SAFE_CALL(cudaGetDeviceCount(&deviceCount));
    if (deviceCount == 0)
        std::cout << "There is no device supporting CUDA" << std::endl;

    for (int dev = 0; dev < deviceCount; ++dev) {
        cudaDeviceProp deviceProp;
        CUDA_SAFE_CALL(cudaGetDeviceProperties(&deviceProp, dev));

        if (dev == 0) {
            if (deviceProp.major == 9999 && deviceProp.minor == 9999)
                std::cout << "There is no device supporting CUDA." << std::endl;
            else if (deviceCount == 1)
                std::cout << "There is 1 device supporting CUDA" << std:: endl;
            else
                std::cout << "There are " << deviceCount <<  " devices supporting CUDA" << std:: endl;
        }

        std::cout << "\nDevice " << dev << ": \"" << deviceProp.name << "\"" << std::endl;
        std::cout << "  Major revision number:                         " << deviceProp.major << std::endl;
        std::cout << "  Minor revision number:                         " << deviceProp.minor << std::endl;
        std::cout << "  Total amount of global memory:                 " << deviceProp.totalGlobalMem << " bytes" << std::endl;
    }
    std::cout << std::endl;
}


void set_device(int argc, char ** argv)
{
    int dev = 0;

    char * dev_str = get_argval(argc, argv, "device");
    if(dev_str != NULL)
        dev = atoi(dev_str);

    cudaSetDevice(dev);

    cudaDeviceProp deviceProp;
    CUDA_SAFE_CALL(cudaGetDeviceProperties(&deviceProp, dev));

    std::cout << "\nRunning on Device " << dev << ": \"" << deviceProp.name << "\"\n";;
}


template <typename IndexType, typename ValueType>
void run_all_kernels(int argc, char **argv)
{
    char * mm_filename = NULL;
    for(int i = 1; i < argc; i++){
        if(argv[i][0] != '-'){
            mm_filename = argv[i];
            break;
        }
    }
    

    csr_matrix<IndexType,ValueType> csr;
    
    if(mm_filename == NULL)
        csr = laplacian_5pt<IndexType,ValueType>(512);
    else
        csr = read_csr_matrix<IndexType,ValueType>(mm_filename);
        
    printf("Using %d-by-%d matrix with %d nonzero values\n", csr.num_rows, csr.num_cols, csr.num_nonzeros); 

    // fill matrix with random values: some matrices have extreme values, 
    // which makes correctness testing difficult, especially in single precision
    srand(13);
    for(IndexType i = 0; i < csr.num_nonzeros; i++){
        csr.Ax[i] = 1.0 - 2.0 * (rand() / (RAND_MAX + 1.0)); 
    }
    
    FILE * fid = fopen(BENCHMARK_OUTPUT_FILE_NAME, "a");
    fprintf(fid, "file=%s rows=%d cols=%d nonzeros=%d\n", mm_filename, csr.num_rows, csr.num_cols, csr.num_nonzeros);
    fclose(fid);

    // specify a device (e.g. --device=1 ) : default 0
    set_device(argc, argv);
    list_devices();

    test_dia_matrix_kernels(csr);
    test_ell_matrix_kernels(csr);
    test_csr_matrix_kernels(csr);
    test_coo_matrix_kernels(csr);
    test_hyb_matrix_kernels(csr);

    
    // specify a paritioning : default None
    char * partition_file = get_argval(argc, argv, "partition");
    if(partition_file != NULL){
        // a partition file was supplied
        IndexType * partition = load_partition<IndexType>(partition_file, csr.num_rows);
        if(partition != NULL){
            test_pkt_matrix_kernels(csr, partition);
            delete_host_array(partition);
        }
    }

    delete_host_matrix(csr);
}

int main(int argc, char** argv)
{
    if (get_arg(argc, argv, "help") != NULL){
        usage(argc, argv);
        return EXIT_SUCCESS;
    }
        

    int precision = 32;
    char * precision_str = get_argval(argc, argv, "precision");
    if(precision_str != NULL)
        precision = atoi(precision_str);
    printf("\nUsing %d-bit floating point precision\n\n", precision);

    if(precision ==  32){
        run_all_kernels<unsigned int, float>(argc,argv);
    }
    else if(precision == 64)
    {
        int current_device = -1;
        cudaDeviceProp properties;
        cudaGetDevice(&current_device);
        cudaGetDeviceProperties(&properties, current_device);
        if (properties.major == 1 && properties.minor < 3)
            std::cerr << "ERROR: Support for \'double\' requires Compute Capability 1.3 or greater\n\n";
        else
        run_all_kernels<unsigned int, double>(argc,argv);
    }
    else{
        usage(argc, argv);
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}

