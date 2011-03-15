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

#include "sparse_io.h"
#include "sparse_operations.h"
#include "spmv_host.h"
#include "spmv_device.h"
#include "test_spmv.h"
#include "benchmark_spmv.h"


template <typename IndexType, typename ValueType>
void test_dia_matrix_kernels(const csr_matrix<IndexType,ValueType>& csr, float max_fill = 3.0)
{
    printf("\n####  Testing DIA Kernels  ####\n");

    IndexType max_diags = static_cast<IndexType>( (max_fill * csr.num_nonzeros) / csr.num_rows + 1 );

    // CREATE DIA MATRIX
    printf("\tcreating dia_matrix:");
    dia_matrix<IndexType,ValueType> dia = csr_to_dia<IndexType,ValueType>(csr, max_diags);
    if (dia.num_nonzeros == 0 && csr.num_nonzeros != 0){
        printf(" number of diagonals (%d) excedes limit (%d)\n", dia.num_diags, max_diags);
        return;
    }
    printf(" found %d diagonals\n", dia.num_diags);

    // TEST FORMAT
    test_spmv_kernel(csr, spmv_csr_serial_host<IndexType,ValueType>, HOST_MEMORY,
                     dia, spmv_dia_serial_host<IndexType,ValueType>, HOST_MEMORY,
                     "dia_serial");

    // TEST KERNELS
    test_spmv_kernel(csr, spmv_csr_serial_host<IndexType,ValueType>, HOST_MEMORY,
                     dia, spmv_dia_device<IndexType,ValueType>,      DEVICE_MEMORY,
                     "dia");
    
    test_spmv_kernel(csr, spmv_csr_serial_host<IndexType,ValueType>, HOST_MEMORY,
                     dia, spmv_dia_tex_device<IndexType,ValueType>,  DEVICE_MEMORY,
                     "dia_tex");
    
    // BENCHMARK KERNELS
    benchmark_spmv_on_device(dia, spmv_dia_device<IndexType, ValueType>,     "dia" );
    benchmark_spmv_on_device(dia, spmv_dia_tex_device<IndexType, ValueType>, "dia_tex" );
    
    delete_host_matrix(dia);
}

template <typename IndexType, typename ValueType>
void test_ell_matrix_kernels(const csr_matrix<IndexType,ValueType>& csr, float max_fill = 3.0)
{
    printf("\n####  Testing ELL Kernels  ####\n");
   
    IndexType max_cols_per_row = static_cast<IndexType>( (max_fill * csr.num_nonzeros) / csr.num_rows + 1 );
  
    // CREATE ELL MATRIX
    printf("\tcreating ell_matrix:");
    ell_matrix<IndexType,ValueType> ell = csr_to_ell<IndexType,ValueType>(csr, max_cols_per_row);
    if (ell.num_nonzeros == 0 && csr.num_nonzeros != 0){
        printf(" num_cols_per_row (%d) excedes limit (%d)\n", ell.num_cols_per_row, max_cols_per_row);
        return;
    }
    printf(" ELL has %d columns per row\n", ell.num_cols_per_row);

    // TEST FORMAT
    test_spmv_kernel(csr, spmv_csr_serial_host<IndexType,ValueType>, HOST_MEMORY,
                     ell, spmv_ell_serial_host<IndexType,ValueType>, HOST_MEMORY,
                     "ell_serial");

    // TEST KERNELS
    test_spmv_kernel(csr, spmv_csr_serial_host<IndexType,ValueType>, HOST_MEMORY,
                     ell, spmv_ell_device<IndexType,ValueType>,      DEVICE_MEMORY,
                     "ell");
    
    test_spmv_kernel(csr, spmv_csr_serial_host<IndexType,ValueType>, HOST_MEMORY,
                     ell, spmv_ell_tex_device<IndexType,ValueType>,  DEVICE_MEMORY,
                     "ell_tex");
    
    // BENCHMARK KERNELS
    benchmark_spmv_on_device(ell, spmv_ell_device<IndexType, ValueType>,     "ell" );
    benchmark_spmv_on_device(ell, spmv_ell_tex_device<IndexType, ValueType>, "ell_tex" );
    
    delete_host_matrix(ell);
}

template <typename IndexType, typename ValueType>
void test_csr_matrix_kernels(const csr_matrix<IndexType,ValueType>& csr)
{
    printf("\n####  Testing CSR Kernels  ####\n");

    // TEST KERNELS
    test_spmv_kernel(csr, spmv_csr_serial_host<IndexType,ValueType>,   HOST_MEMORY,
                     csr, spmv_csr_scalar_device<IndexType,ValueType>, DEVICE_MEMORY,
                     "csr_scalar");
    
    test_spmv_kernel(csr, spmv_csr_serial_host<IndexType,ValueType>,       HOST_MEMORY,
                     csr, spmv_csr_scalar_tex_device<IndexType,ValueType>, DEVICE_MEMORY,
                     "csr_scalar_tex");
    
    test_spmv_kernel(csr, spmv_csr_serial_host<IndexType,ValueType>,   HOST_MEMORY,
                     csr, spmv_csr_vector_device<IndexType,ValueType>, DEVICE_MEMORY,
                     "csr_vector");
        
    test_spmv_kernel(csr, spmv_csr_serial_host<IndexType,ValueType>,   HOST_MEMORY,
                     csr, spmv_csr_vector_tex_device<IndexType,ValueType>, DEVICE_MEMORY,
                     "csr_vector_tex");

    // BENCHMARK KERNELS
    benchmark_spmv_on_host(csr,   spmv_csr_serial_host<IndexType, ValueType>,       "csr_serial" );
    benchmark_spmv_on_device(csr, spmv_csr_scalar_device<IndexType, ValueType>,     "csr_scalar" );
    benchmark_spmv_on_device(csr, spmv_csr_scalar_tex_device<IndexType, ValueType>, "csr_scalar_tex" );
    benchmark_spmv_on_device(csr, spmv_csr_vector_device<IndexType, ValueType>,     "csr_vector" );
    benchmark_spmv_on_device(csr, spmv_csr_vector_tex_device<IndexType, ValueType>, "csr_vector_tex" );
}

template <typename IndexType, typename ValueType>
void test_coo_matrix_kernels(const csr_matrix<IndexType,ValueType>& csr)
{
    printf("\n####  Testing COO Kernels  ####\n");

    // CREATE COO MATRIX
    printf("\tcreating coo_matrix:");
    coo_matrix<IndexType,ValueType> coo = csr_to_coo<IndexType,ValueType>(csr);
    printf("\n");
    
    // TEST FORMAT
    test_spmv_kernel(csr, spmv_csr_serial_host<IndexType,ValueType>, HOST_MEMORY,
                     coo, spmv_coo_serial_host<IndexType,ValueType>, HOST_MEMORY,
                     "coo_serial");

    // TEST KERNELS
    test_spmv_kernel(csr, spmv_csr_serial_host<IndexType,ValueType>, HOST_MEMORY,
                     coo, spmv_coo_flat_device<IndexType,ValueType>, DEVICE_MEMORY,
                     "coo_flat");

    test_spmv_kernel(csr, spmv_csr_serial_host<IndexType,ValueType>,     HOST_MEMORY,
                     coo, spmv_coo_flat_tex_device<IndexType,ValueType>, DEVICE_MEMORY,
                     "coo_flat_tex");
    
    // BENCHMARK KERNELS
    benchmark_spmv_on_device(coo, spmv_coo_flat_device<IndexType, ValueType>,     "coo_flat");
    benchmark_spmv_on_device(coo, spmv_coo_flat_tex_device<IndexType, ValueType>, "coo_flat_tex");
    
    delete_host_matrix(coo);
}


template <typename IndexType, typename ValueType>
void test_hyb_matrix_kernels(const csr_matrix<IndexType,ValueType>& csr)
{
    printf("\n####  Testing HYB Kernels  ####\n");
    
    const IndexType num_cols_per_row = compute_hyb_cols_per_row(csr);

    // CREATE HYB MATRIX
    printf("\tcreating hyb_matrix:");
    hyb_matrix<IndexType,ValueType> hyb = csr_to_hyb<IndexType,ValueType>(csr, num_cols_per_row);
    printf(" ELL has %d columns. NNZ: [ELL %4.1f%%]  [COO %4.1f%%]\n", 
                    num_cols_per_row,
                    100 * (float) hyb.ell.num_nonzeros / (float) csr.num_nonzeros, 
                    100 * (float) hyb.coo.num_nonzeros / (float) csr.num_nonzeros);

    // TEST FORMAT
    test_spmv_kernel(csr, spmv_csr_serial_host<IndexType,ValueType>, HOST_MEMORY,
                     hyb, spmv_hyb_serial_host<IndexType,ValueType>, HOST_MEMORY,
                     "hyb_serial");

    // TEST KERNELS
    test_spmv_kernel(csr, spmv_csr_serial_host<IndexType,ValueType>, HOST_MEMORY,
                     hyb, spmv_hyb_device<IndexType,ValueType>,      DEVICE_MEMORY,
                     "hyb");
    
    test_spmv_kernel(csr, spmv_csr_serial_host<IndexType,ValueType>, HOST_MEMORY,
                     hyb, spmv_hyb_tex_device<IndexType,ValueType>,  DEVICE_MEMORY,
                     "hyb_tex");
    
    // BENCHMARK KERNELS
    benchmark_spmv_on_device(hyb, spmv_hyb_device<IndexType, ValueType>, "hyb" );
    benchmark_spmv_on_device(hyb, spmv_hyb_tex_device<IndexType, ValueType>, "hyb_tex" );
    
    delete_host_matrix(hyb);
}


template <typename IndexType, typename ValueType>
void test_pkt_matrix_kernels(const csr_matrix<IndexType,ValueType>& csr, 
                                const IndexType * partition)
{
    printf("\n####  Testing PKT Kernels  ####\n");

    if(csr.num_rows != csr.num_cols){
        printf("\tPacket matrix does not support non-square matrices\n");
        return;
    }
  
    // CREATE PACKET MATRIX
    printf("\tcreating pkt_matrix:");
    pkt_matrix<IndexType,ValueType> pkt = csr_to_pkt<IndexType,ValueType>(csr, partition);
    printf("\n");
    
    // TEST FORMAT
    test_spmv_kernel(csr, spmv_csr_serial_host<IndexType,ValueType>, HOST_MEMORY,
                     pkt, spmv_pkt_serial_host<IndexType,ValueType>,        HOST_MEMORY,
                     "pkt_serial");

    // TEST KERNELS
    test_spmv_kernel(csr, spmv_csr_serial_host<IndexType,ValueType>, HOST_MEMORY,
                     pkt, full_spmv_pkt_device<IndexType,ValueType>, DEVICE_MEMORY,
                     "pkt");
    
    // BENCHMARK KERNELS
    benchmark_spmv_on_device(pkt, spmv_pkt_device<IndexType, ValueType>, "pkt" );
    
    delete_host_matrix(pkt);
}



template <typename IndexType>
IndexType * load_partition(const char * partition_filename, IndexType N)
{

    FILE * partition_file = fopen(partition_filename, "r");

    if (partition_file == NULL){
        fprintf(stderr,"Error opening parition file (%s)\n", partition_filename);
        return NULL;
    }

    IndexType * partition = new_host_array<IndexType>(N);
   
    IndexType i = 0;
    while(feof(partition_file)== 0){
        unsigned int partnum;

        if( fscanf(partition_file," %u \n",&partnum) != 1) { 
            fclose(partition_file);
            delete_host_array(partition); 
            return NULL; 
         }

        partition[i++] = (IndexType) partnum;
    }

    fclose(partition_file);

    if(i != N){
        fprintf(stderr, "partiton file contained %u values, expected %u\n", i, N);
        delete_host_array(partition); 
        return NULL; 
    } else {
        return partition;
    }
}



