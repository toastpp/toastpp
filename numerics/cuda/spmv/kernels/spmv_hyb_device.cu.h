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
#include "kernels/spmv_ell_device.cu.h"
#include "kernels/spmv_coo_flat_device.cu.h"

// SpMV kernels for the hybrid ELL/COO matrix format.

template <typename IndexType, typename ValueType>
void spmv_hyb_device(const hyb_matrix<IndexType,ValueType>& d_hyb, 
                     const ValueType * d_x, 
                           ValueType * d_y)
{
    spmv_ell_device(d_hyb.ell, d_x, d_y);
    spmv_coo_flat_device(d_hyb.coo, d_x, d_y);
}

template <typename IndexType, typename ValueType>
void spmv_hyb_tex_device(const hyb_matrix<IndexType,ValueType>& d_hyb, 
                         const ValueType * d_x, 
                               ValueType * d_y)
{
    spmv_ell_tex_device(d_hyb.ell, d_x, d_y);
    spmv_coo_flat_tex_device(d_hyb.coo, d_x, d_y);
}

