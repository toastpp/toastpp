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

#include <vector>
#include <algorithm>
#include <assert.h>

#include "mem.h"

template<class T>
T * vector_to_array(const std::vector<T>& V){
    T * arr = new_host_array<T>(V.size());
    std::copy(V.begin(),V.end(),arr);
    return arr;
}


template<class IndexType>
std::vector<IndexType> bincount(const IndexType * V, const IndexType N){
    IndexType maximum = *std::max_element(V, V + N);

    std::vector<IndexType> bc(maximum+1,(IndexType) 0);

    for(IndexType i = 0; i < N; i++)
    {
        bc[V[i]]++;
    }

    return bc;
}

template<class T>
std::vector<T> cumsum(std::vector<T> V){
    std::vector<T> cs(V.size()+1,0);

    for(size_t i = 0; i < V.size(); i++){
        cs[i+1] = cs[i] + V[i];
    }

    return cs;
}


