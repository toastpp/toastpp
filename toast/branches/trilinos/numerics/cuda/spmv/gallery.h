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
#include "mem.h"

/*
 * The standard 5-point finite difference approximation
 * to the Laplacian operator on a regular N-by-N grid.
 */
template <typename IndexType, typename ValueType>
csr_matrix<IndexType,ValueType> laplacian_5pt(const IndexType N)
{
    csr_matrix<IndexType,ValueType> csr;
    csr.num_rows = N*N;
    csr.num_cols = N*N;
    csr.num_nonzeros = 5*N*N - 4*N; 

    csr.Ap = new_host_array<IndexType>(csr.num_rows + 1);
    csr.Aj = new_host_array<IndexType>(csr.num_nonzeros);
    csr.Ax = new_host_array<ValueType>(csr.num_nonzeros);

    IndexType nz = 0;

    for(IndexType i = 0; i < N; i++){
        for(IndexType j = 0; j < N; j++){
            IndexType indx = N*i + j;

            if (i > 0){
                csr.Aj[nz] = indx - N;
                csr.Ax[nz] = -1;
                nz++;
            }

            if (j > 0){
                csr.Aj[nz] = indx - 1;
                csr.Ax[nz] = -1;
                nz++;
            }

            csr.Aj[nz] = indx;
            csr.Ax[nz] = 4;
            nz++;

            if (j < N - 1){
                csr.Aj[nz] = indx + 1;
                csr.Ax[nz] = -1;
                nz++;
            }

            if (i < N - 1){
                csr.Aj[nz] = indx + N;
                csr.Ax[nz] = -1;
                nz++;
            }
            
            csr.Ap[indx + 1] = nz;
        }
    }

    return csr;
}


