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

#include <algorithm>
#include <vector>
#include "sparse_formats.h"
#include "sparse_operations.h"
#include "mem.h"

////////////////////////////////////////////////////////////////////////////////
//! Partition a graph given in csr_matrix format into a given number of parts
//!
//!  Note: IndexType must be 'int' or 'unsigned int'
//!
////////////////////////////////////////////////////////////////////////////////

extern "C" {
#include <metis.h>
extern void METIS_PartGraphRecursive(int*, idxtype*, idxtype*, idxtype*, idxtype*, int*, int*, int*, int*, int*, idxtype*);
extern void METIS_PartGraphKway(int*, idxtype*, idxtype*, idxtype*, idxtype*, int*, int*, int*, int*, int*, idxtype*);
};


/*
 *  Return true if CSR matrix A has a symmetric nonzero storage pattern and false otherwise.
 *
 *  Note: the matrix values are irrelevant, i.e. [0 1] is structurally symmetric
 *                                               [2 0]  
 */
template <class IndexType, class ValueType>
bool csr_is_structurally_symmetric(const IndexType * Ap, 
                                   const IndexType * Aj, 
                                   const ValueType * Ax, 
                                   const IndexType num_rows, 
                                   const IndexType num_cols)
{
    if(num_rows != num_cols)
        return false;
    
    IndexType num_nonzeros = Ap[num_rows];

    //Compute B = A^T  (B will have sorted indices)
    IndexType * Bp = new_host_array<IndexType>(num_rows + 1);
    IndexType * Bj = new_host_array<IndexType>(num_nonzeros);
    ValueType * Bx = new_host_array<ValueType>(num_nonzeros);

    csr_transpose<IndexType,ValueType>(Ap,Aj,Ax,num_rows,num_cols,Bp,Bj,Bx);

    //see if Ap == Bp
    for(IndexType i = 0; i < num_rows; i++){ //Ap[num_rows] == Bp[num_rows] by construction
        if(Ap[i] != Bp[i]){
            delete_host_array(Bp); delete_host_array(Bj); delete_host_array(Bx);
            return false;
        }
    }

    //Compute C = A^T^T = A (C will have sorted indices)
    IndexType * Cp = new_host_array<IndexType>(num_rows + 1);
    IndexType * Cj = new_host_array<IndexType>(num_nonzeros);
    ValueType * Cx = new_host_array<ValueType>(num_nonzeros);

    csr_transpose(Bp,Bj,Bx,num_rows,num_cols,Cp,Cj,Cx);

    for(IndexType i = 0; i < num_nonzeros; i++){
        if(Bj[i] != Cj[i]){         
            delete_host_array(Bp); delete_host_array(Bj); delete_host_array(Bx);
            delete_host_array(Cp); delete_host_array(Cj); delete_host_array(Cx);
            return false;
        }
    }

    delete_host_array(Bp); delete_host_array(Bj); delete_host_array(Bx);
    delete_host_array(Cp); delete_host_array(Cj); delete_host_array(Cx);
    return true;    
}

template <class IndexType, class ValueType>
bool csr_is_structurally_symmetric(const csr_matrix<IndexType,ValueType>& csr){
    return csr_is_structurally_symmetric(csr.Ap,csr.Aj,csr.Ax,csr.num_rows,csr.num_cols);
}




template <class IndexType, class ValueType>
IndexType partition_csr(const csr_matrix<IndexType,ValueType>& graph, 
                        int num_parts,  
                        std::vector<IndexType>& partition,                 
                        bool Kway = false)
{

    if(graph.num_rows != graph.num_cols){
        printf("matrix is nonsquare\n");
        exit(EXIT_FAILURE);
    }

    if(num_parts < 2){  
        // only one partition
        std::fill(partition.begin(), partition.end(), 0);
        return 0;
    }
    
    int * int_partition = new_host_array<int>(partition.size()); //because METIS expects an int array
        
    int options[10];
    int numflag = 0, wgtflag = 0, edgecut;
    options[0] = 0; 
        
    printf(" [A=A^T:");
    bool is_symmetric = csr_is_structurally_symmetric(graph);
    if (is_symmetric)
        printf(" yes]");
    else
        printf(" no]");


    if(is_symmetric){
        if(Kway)
            METIS_PartGraphKway((int *)&graph.num_rows, (int *)graph.Ap, (int *)graph.Aj, NULL, NULL, 
                                &wgtflag, &numflag, &num_parts, options, &edgecut, int_partition);      
        else
            METIS_PartGraphRecursive((int *)&graph.num_rows, (int *)graph.Ap, (int *)graph.Aj, NULL, NULL, 
                                     &wgtflag, &numflag, &num_parts, options, &edgecut, int_partition);     
    } else {
        //printf("Matrix graph is nonsymmetric, using A^T + A instead\n");
        
        //coo_matrix<unsigned int,float> coo_A  = csr_to_coo(graph);

        coo_matrix<IndexType,ValueType> coo_AtA;
        coo_AtA.num_rows = graph.num_rows;
        coo_AtA.num_cols = graph.num_cols;
        coo_AtA.num_nonzeros = 2*graph.num_nonzeros;
        coo_AtA.I = new_host_array<IndexType>(2*graph.num_nonzeros);
        coo_AtA.J = new_host_array<IndexType>(2*graph.num_nonzeros);
        coo_AtA.V = new_host_array<ValueType>(2*graph.num_nonzeros);
        
        csr_to_coo(graph.Ap,graph.Aj,graph.Ax,
                   graph.num_rows,graph.num_cols,graph.num_nonzeros,
                   coo_AtA.I,coo_AtA.J,coo_AtA.V);
        csr_to_coo(graph.Ap,graph.Aj,graph.Ax,
                   graph.num_rows,graph.num_cols,graph.num_nonzeros,
                   coo_AtA.J+graph.num_nonzeros,coo_AtA.I+graph.num_nonzeros,coo_AtA.V+graph.num_nonzeros);
            

        std::fill(coo_AtA.V,coo_AtA.V + 2*graph.num_nonzeros,(float)1);
        
        csr_matrix<IndexType,ValueType> csr_AtA = coo_to_csr(coo_AtA);

        assert(csr_is_structurally_symmetric(csr_AtA));
        if (Kway)
            METIS_PartGraphKway((int *)&csr_AtA.num_rows, (int *)csr_AtA.Ap, (int *)csr_AtA.Aj, NULL, NULL, 
                                &wgtflag, &numflag, &num_parts, options, &edgecut, int_partition);          
        else
            METIS_PartGraphRecursive((int *)&csr_AtA.num_rows, (int *)csr_AtA.Ap, (int *)csr_AtA.Aj, NULL, NULL, 
                                     &wgtflag, &numflag, &num_parts, options, &edgecut, int_partition);         

        delete_host_array(coo_AtA.I);  delete_host_array(coo_AtA.J);  delete_host_array(coo_AtA.V);
        delete_host_array(csr_AtA.Ap); delete_host_array(csr_AtA.Aj); delete_host_array(csr_AtA.Ax);
    }

    printf(" [edgecut: %4.1f%%]", 100* (double) edgecut/ (double) graph.num_nonzeros);

    //copy partition to vector
    std::copy(int_partition, int_partition + partition.size(), partition.begin());
    delete_host_array(int_partition);

    return (IndexType) edgecut;
}


