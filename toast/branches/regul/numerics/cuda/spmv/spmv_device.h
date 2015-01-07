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

////////////////////////////////////////////////////////////////////////////////
//! GPU SpMV kernels
//////////////////////////////////////////////////////////////////////////////////

#include "kernels/spmv_dia_device.cu.h"
#include "kernels/spmv_ell_device.cu.h"
#include "kernels/spmv_coo_flat_device.cu.h"
#include "kernels/spmv_csr_scalar_device.cu.h"
#include "kernels/spmv_csr_vector_device.cu.h"
#include "kernels/spmv_hyb_device.cu.h"
#include "kernels/spmv_pkt_device.cu.h"

