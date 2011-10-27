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

#include <thrust/version.h>

#if THRUST_VERSION >= 100500
#include <thrust/detail/backend/cuda/arch.h>
#else
#include <thrust/detail/device/cuda/arch.h>
namespace thrust{
  namespace detail{
    namespace backend = device;
  }
}
#endif

namespace cusp
{
namespace detail
{
namespace device
{
namespace arch
{

template <typename KernelFunction>
size_t max_active_blocks(KernelFunction kernel, const size_t CTA_SIZE, const size_t dynamic_smem_bytes)
{
    return thrust::detail::backend::cuda::arch::max_active_blocks(kernel, CTA_SIZE, dynamic_smem_bytes);
}

} // end namespace arch
} // end namespace device
} // end namespace detail
} // end namespace cusp

