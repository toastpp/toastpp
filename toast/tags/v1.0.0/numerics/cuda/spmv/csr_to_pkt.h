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

#include <assert.h>
#include <vector>
#include <queue>
#include <algorithm>
#include <stdio.h>

#include "sparse_operations.h"
#include "array_utils.h"

////////////////////////////////////////////////////////////////////////////////
//! Convert a csr_matrix to pkt_matrix format
////////////////////////////////////////////////////////////////////////////////

/*
 *  Assign rows to threads using Heuristic #1 (least occupied bin first)
 *
 *  Returns the maximum bin size (total work assigned to any thread).
 */
template <class IndexType>
IndexType
compute_thread_assignments_H1(IndexType part_num, 
                              IndexType threads_per_packet,
                              const std::vector<IndexType>& partition_ptr, 
                              const std::vector<IndexType>& permute_new_to_old, 
                              const std::vector<IndexType>& row_costs, 
                                    std::vector<IndexType>& thread_assignments)
{
    //priority queue of tuples <bin size,bin number>
    std::priority_queue< std::pair<IndexType, IndexType>, 
        std::vector<std::pair<IndexType, IndexType> >,
          std::greater< std::pair<IndexType, IndexType> > > bin_queue;


    for(IndexType i = 0; i < threads_per_packet; i++){
        bin_queue.push(std::make_pair(0,i + part_num*threads_per_packet));
    }   
    
    IndexType max_bin_size = 0;
    for(IndexType i = partition_ptr[part_num]; i < partition_ptr[part_num+1]; i++){
        IndexType row = permute_new_to_old[i];      
        if( row_costs[row] > 0 ){           
            thread_assignments[row] = bin_queue.top().second;

            max_bin_size = std::max(max_bin_size,bin_queue.top().first + row_costs[row]);

            bin_queue.push(std::make_pair(bin_queue.top().first + row_costs[row],bin_queue.top().second));
            bin_queue.pop();
        }
    }   

    return max_bin_size;
}



/*
 *  Assign rows to threads using Heuristic #2 (sort then zigzag placement)
 *
 *  Returns the maximum bin size (total work assigned to any thread).
 */
template <class IndexType>
IndexType
compute_thread_assignments_H2(IndexType part_num, 
                              IndexType threads_per_packet,
                              const std::vector<IndexType>& partition_ptr, 
                              const std::vector<IndexType>& permute_new_to_old, 
                              const std::vector<IndexType>& row_costs, 
                                    std::vector<IndexType>& thread_assignments)
{
    // pairs are of the form <cost,row#>
    std::vector< std::pair<IndexType, IndexType> > cost_pairs;
    std::vector< IndexType > costs(threads_per_packet,0);

    //IndexType max_bin_size = 0;
    for(IndexType i = partition_ptr[part_num]; i < partition_ptr[part_num+1]; i++){
        IndexType row = permute_new_to_old[i];  
        IndexType row_cost = row_costs[row];

        cost_pairs.push_back(std::make_pair(row_cost,row));
    }

    std::sort(cost_pairs.begin(),cost_pairs.end(),std::greater<std::pair<IndexType, IndexType> >());


    IndexType base_thread_index =  part_num*threads_per_packet;

    IndexType pos  = 0;
    int step = 1;
    for(IndexType i = 0; i < (IndexType) cost_pairs.size(); i++){
        thread_assignments[cost_pairs[i].second] = pos + base_thread_index;
        costs[pos] += cost_pairs[i].first;

        if (pos == threads_per_packet-1 && step == 1){          
            step = -1;
        } else if(pos == 0 && step == -1){
            step = 1;
        } else {
            pos += step;
        }
        
    }

    return *(std::max_element(costs.begin(),costs.end()));  
}


template <class IndexType, class ValueType>
void schedule_work(pkt_matrix<IndexType,ValueType>& pm, 
                   packet_array<IndexType,ValueType>&  pa, 
                   const csr_matrix<IndexType,ValueType>& csr, 
                   const std::vector<IndexType>& thread_assignments,
                   const std::vector<IndexType>& packet_ptr,
                   const IndexType * partition)
{   
    //initialize thread offset pointers into row,col,data arrays
    for(IndexType i = 0, thread_id = 0; i < pm.num_packets; i++){
        IndexType offset = pm.threads_per_packet * packet_ptr[i];
        for(IndexType j = 0; j < pm.threads_per_packet; j++){
            pa.pos_start[thread_id] = offset;
            pa.pos_end[thread_id]   = offset;
            thread_id++; offset++;
        }
    }

    for(IndexType packet_num = 0; packet_num < pm.num_packets; packet_num++){
        IndexType base_row = pm.row_ptr[packet_num];        

        for(IndexType new_i = pm.row_ptr[packet_num]; new_i < pm.row_ptr[packet_num+1]; new_i++){
            IndexType old_i = pm.permute_new_to_old[new_i]; //old row index 
            IndexType thread_id  = thread_assignments[old_i];       

            IndexType& offset = pa.pos_end[thread_id];          
            
            for(IndexType jj = csr.Ap[old_i]; jj < csr.Ap[old_i+1]; jj++){
                IndexType old_j = csr.Aj[jj]; //old column index
                
                if(packet_num == partition[old_j]){
                    const IndexType row = new_i - base_row;                        // relative row index
                    const IndexType col = pm.permute_old_to_new[old_j] - base_row; // relative col index

                    pa.index_array[offset] = pkt_pack_indices(row,col);
                    pa.data_array[offset]  = csr.Ax[jj];

                    offset += pm.threads_per_packet;
                }
            }
        }
    }


}


template <class IndexType, class ValueType>
pkt_matrix<IndexType,ValueType> 
csr_to_pkt(const csr_matrix<IndexType,ValueType>& csr, 
           const IndexType * partition,
           const IndexType threads_per_packet = 512,
           int packing_heuristic = 2)
{
    if(csr.num_rows != csr.num_cols){
        printf("matrix must be square\n");
        exit(EXIT_FAILURE);
    }

    std::vector<IndexType> partition_sizes = bincount(partition, csr.num_rows); // number of rows in each partition

    const IndexType N = csr.num_rows;
    const IndexType num_parts = partition_sizes.size();
    
    std::vector<IndexType> permute_old_to_new(csr.num_rows);
    std::vector<IndexType> permute_new_to_old(csr.num_rows);
    std::vector<IndexType> partition_ptr  = cumsum(partition_sizes);    

    IndexType max_rows_per_packet = 0;
    {   
        //Given a partition, find a permutation of the rows such that
        //all the rows in partition P proceed those in partition P+1
        std::vector<IndexType> part_count(num_parts,0);

        for(IndexType i = 0; i < N; i++){
            IndexType part = partition[i];
            IndexType new_pos = part_count[part] + partition_ptr[part];
            permute_old_to_new[i]       = new_pos;
            permute_new_to_old[new_pos] = i;
            part_count[part]++;
        }       

        for(IndexType i = 0; i < (IndexType) num_parts; i++)
            max_rows_per_packet = std::max(max_rows_per_packet, part_count[i]);
    }


    //compute the number of local and global members in each row
    std::vector<IndexType> local_row_costs(csr.num_rows);
        
    IndexType local_nonzeros = 0;
    IndexType global_nonzeros = 0;
    for(IndexType i = 0; i < csr.num_rows; i++){
        IndexType local_cost = 0;

        for(IndexType jj = csr.Ap[i]; jj < csr.Ap[i+1]; jj++){  
            IndexType j = csr.Aj[jj];
            if( partition[i] == partition[j] ){
                local_cost++;
            } else {
                global_nonzeros++;
            }
        }

        local_nonzeros     += local_cost;
        local_row_costs[i]  = local_cost;
    }


    // assign each local/global row of the matrix to a thread
    std::vector<IndexType> local_thread_assignments(csr.num_rows);
    std::vector<IndexType> local_packet_sizes(num_parts);

    assert(packing_heuristic == 1 || packing_heuristic == 2);

    for(IndexType i = 0; i < (IndexType) num_parts; i++){
        if(packing_heuristic == 1){
            local_packet_sizes[i]  = compute_thread_assignments_H1(i,threads_per_packet,
                                                                   partition_ptr, permute_new_to_old,
                                                                   local_row_costs, local_thread_assignments);
        } else {
            local_packet_sizes[i]  = compute_thread_assignments_H2(i,threads_per_packet,
                                                                   partition_ptr, permute_new_to_old,
                                                                   local_row_costs, local_thread_assignments);
        }
    }
    std::vector<IndexType> local_packet_ptr  = cumsum(local_packet_sizes);
    
    //initialize packet matrix
    pkt_matrix<IndexType,ValueType> p_mtx;
    p_mtx.num_rows = csr.num_rows;
    p_mtx.num_cols = csr.num_cols;
    p_mtx.num_nonzeros = csr.num_nonzeros;
    p_mtx.num_packets  = num_parts;
    p_mtx.threads_per_packet = threads_per_packet;
    p_mtx.max_rows_per_packet = max_rows_per_packet;

    p_mtx.row_ptr = vector_to_array(partition_ptr);

    p_mtx.permute_new_to_old = vector_to_array(permute_new_to_old);
    p_mtx.permute_old_to_new = vector_to_array(permute_old_to_new);
    
    //initialize packet array
    //IndexType total_local_work  = std::accumulate(local_packet_sizes.begin(),local_packet_sizes.end(),0);
    IndexType total_local_work  = 0;
    for(IndexType i = 0; i < (IndexType) local_packet_sizes.size(); i++)
        total_local_work += local_packet_sizes[i];

    p_mtx.packets.total_cycles = total_local_work;
    p_mtx.packets.pos_start    = new_host_array<IndexType>(threads_per_packet * p_mtx.num_packets);
    p_mtx.packets.pos_end      = new_host_array<IndexType>(threads_per_packet * p_mtx.num_packets);
    p_mtx.packets.index_array  = new_host_array<PackedIndexType>(threads_per_packet * total_local_work);
    p_mtx.packets.data_array   = new_host_array<ValueType>(threads_per_packet * total_local_work);

    schedule_work(p_mtx, p_mtx.packets, csr, local_thread_assignments, local_packet_ptr, partition);

    //initialize coo_matrix
    p_mtx.coo.num_rows = csr.num_rows;
    p_mtx.coo.num_cols = csr.num_cols;
    p_mtx.coo.num_nonzeros = global_nonzeros;
    
    p_mtx.coo.I = new_host_array<IndexType>(global_nonzeros); 
    p_mtx.coo.J = new_host_array<IndexType>(global_nonzeros); 
    p_mtx.coo.V = new_host_array<ValueType>(global_nonzeros); 
            
    for(IndexType i_new = 0, nnz = 0; i_new < csr.num_rows; i_new++){
        IndexType i_old = permute_new_to_old[i_new];
        for(IndexType jj = csr.Ap[i_old]; jj < csr.Ap[i_old + 1]; jj++){
            IndexType j_old = csr.Aj[jj];

            if( partition[i_old] == partition[j_old] ){ continue; } // nz stored in packet

            IndexType j_new = permute_old_to_new[j_old];

            p_mtx.coo.I[nnz] = i_new;
            p_mtx.coo.J[nnz] = j_new;
            p_mtx.coo.V[nnz] = csr.Ax[jj];

            nnz++;
        }
    }


    IndexType total_values = total_local_work * threads_per_packet + global_nonzeros;
    printf(" [waste: %4.1f%% edgecut %4.1f%%]", 
            100 * (double) (total_values - csr.num_nonzeros) / (double) total_values,
            100 * (double) global_nonzeros / (double) csr.num_nonzeros);

    return p_mtx;
}

