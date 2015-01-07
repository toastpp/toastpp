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

#include "mem.h"

////////////////////////////////////////////////////////////////////////////////
//! Defines the following sparse matrix formats
//
// DIA - Diagonal
// ELL - ELLPACK/ITPACK
// CSR - Compressed Sparse Row
// CSC - Compressed Sparse Column
// COO - Coordinate
// PKT - Packet
////////////////////////////////////////////////////////////////////////////////

template<typename IndexType>
struct matrix_shape
{
    typedef IndexType index_type;
    IndexType num_rows, num_cols, num_nonzeros;
};

// DIAgonal matrix
template <typename IndexType, typename ValueType>
struct dia_matrix : public matrix_shape<IndexType> 
{
    typedef IndexType index_type;
    typedef ValueType value_type;

    IndexType stride;
    IndexType num_diags;

    int       * diag_offsets;  //diagonal offsets (must be a signed type)
    ValueType * diag_data;     //nonzero values stored in a (num_diags x num_cols) matrix 
};

// ELLPACK/ITPACK matrix format
template <typename IndexType, typename ValueType>
struct ell_matrix : public matrix_shape<IndexType> 
{
    typedef IndexType index_type;
    typedef ValueType value_type;

    IndexType stride;
    IndexType num_cols_per_row;

    IndexType * Aj;           //column indices stored in a (cols_per_row x stride) matrix
    ValueType * Ax;           //nonzero values stored in a (cols_per_row x stride) matrix
};

// COOrdinate matrix (aka IJV or Triplet format)
template <typename IndexType, typename ValueType>
struct coo_matrix : public matrix_shape<IndexType> 
{
    typedef IndexType index_type;
    typedef ValueType value_type;

    IndexType * I;  //row indices
    IndexType * J;  //column indices
    ValueType * V;  //nonzero values
};


/*
 *  Compressed Sparse Row matrix (aka CRS)
 */
template <typename IndexType, typename ValueType>
struct csr_matrix : public matrix_shape<IndexType>
{
    typedef IndexType index_type;
    typedef ValueType value_type;

    IndexType * Ap;  //row pointer
    IndexType * Aj;  //column indices
    ValueType * Ax;  //nonzeros
};



/*
 *  Hybrid ELL/COO format
 */
template <typename IndexType, typename ValueType>
struct hyb_matrix : public matrix_shape<IndexType>
{
    typedef IndexType index_type;
    typedef ValueType value_type;

    ell_matrix<IndexType,ValueType> ell; //ELL portion
    coo_matrix<IndexType,ValueType> coo; //COO portion
};



/*
 *  Packet matrix
 */
typedef unsigned int PackedIndexType;
template <typename IndexType, typename ValueType>
struct packet_array : public matrix_shape<IndexType>
{
    typedef IndexType index_type;
    typedef ValueType value_type;

    PackedIndexType * index_array;  // compressed row/col indices
    ValueType * data_array;         // nonzero values

    IndexType * pos_start;          // start ptr into index and data arrays for each thread
    IndexType * pos_end;      

    IndexType total_cycles;         // total amount of work in each thread lane
};

template <typename IndexType, typename ValueType>
struct pkt_matrix : public matrix_shape<IndexType>
{ 
    typedef IndexType index_type;
    typedef ValueType value_type;

    IndexType threads_per_packet;    // # of threads in a block, e.g. 256
    IndexType max_rows_per_packet;   // maximum over all packets

    IndexType num_packets;
    IndexType * row_ptr;             // packet i corresponds to rows row_ptr[i] through row_ptr[i+1] - 1
    IndexType * permute_old_to_new;  
    IndexType * permute_new_to_old;

    packet_array<IndexType,ValueType> packets;
    coo_matrix<IndexType,ValueType> coo;
};

// store row index in upper 16 bits and col index in lower 16 bits
#define pkt_pack_indices(row,col)          (  (row << 16) + col  )
#define pkt_unpack_row_index(packed_index) ( packed_index >> 16  )  
#define pkt_unpack_col_index(packed_index) (packed_index & 0xFFFF)


////////////////////////////////////////////////////////////////////////////////
//! sparse matrix memory management 
////////////////////////////////////////////////////////////////////////////////

template <typename IndexType, typename ValueType>
void delete_dia_matrix(dia_matrix<IndexType,ValueType>& dia, const memory_location loc){
    delete_array(dia.diag_offsets, loc);  delete_array(dia.diag_data, loc);
}

template <typename IndexType, typename ValueType>
void delete_ell_matrix(ell_matrix<IndexType,ValueType>& ell, const memory_location loc){
    delete_array(ell.Aj, loc);  delete_array(ell.Ax, loc);
}

template <typename IndexType, typename ValueType>
void delete_csr_matrix(csr_matrix<IndexType,ValueType>& csr, const memory_location loc){
    delete_array(csr.Ap, loc);  delete_array(csr.Aj, loc);   delete_array(csr.Ax, loc);
}

template <typename IndexType, typename ValueType>
void delete_coo_matrix(coo_matrix<IndexType,ValueType>& coo, const memory_location loc){
    delete_array(coo.I, loc);   delete_array(coo.J, loc);   delete_array(coo.V, loc);
}

template <typename IndexType, typename ValueType>
void delete_packet_array(packet_array<IndexType,ValueType>& pa, const memory_location loc){
    delete_array(pa.index_array, loc); delete_array(pa.data_array, loc);
    delete_array(pa.pos_start, loc);   delete_array(pa.pos_end, loc);
}

template <typename IndexType, typename ValueType>
void delete_hyb_matrix(hyb_matrix<IndexType,ValueType>& hyb, const memory_location loc){
    delete_ell_matrix(hyb.ell, loc);
    delete_coo_matrix(hyb.coo, loc);
}

template <typename IndexType, typename ValueType>
void delete_pkt_matrix(pkt_matrix<IndexType, ValueType>& pm, const memory_location loc)
{
    delete_array(pm.row_ptr, loc);   
    delete_array(pm.permute_new_to_old, loc);   
    delete_array(pm.permute_old_to_new, loc);

    delete_packet_array(pm.packets,  loc);  
    delete_coo_matrix(pm.coo, loc);
}


////////////////////////////////////////////////////////////////////////////////
//! host functions
////////////////////////////////////////////////////////////////////////////////

template <typename IndexType, typename ValueType>
void delete_host_matrix(dia_matrix<IndexType,ValueType>& dia){ delete_dia_matrix(dia, HOST_MEMORY); }

template <typename IndexType, typename ValueType>
void delete_host_matrix(ell_matrix<IndexType,ValueType>& ell){ delete_ell_matrix(ell, HOST_MEMORY); }

template <typename IndexType, typename ValueType>
void delete_host_matrix(coo_matrix<IndexType,ValueType>& coo){ delete_coo_matrix(coo, HOST_MEMORY); }

template <typename IndexType, typename ValueType>
void delete_host_matrix(csr_matrix<IndexType,ValueType>& csr){ delete_csr_matrix(csr, HOST_MEMORY); }

template <class IndexType, class ValueType>
void delete_host_matrix(hyb_matrix<IndexType,ValueType>& hyb){  delete_hyb_matrix(hyb, HOST_MEMORY); }

template <typename IndexType, typename ValueType>
void delete_host_matrix(pkt_matrix<IndexType, ValueType>& pm){ delete_pkt_matrix(pm, HOST_MEMORY); }

////////////////////////////////////////////////////////////////////////////////
//! device functions
////////////////////////////////////////////////////////////////////////////////

template <typename IndexType, typename ValueType>
void delete_device_matrix(dia_matrix<IndexType,ValueType>& dia){ delete_dia_matrix(dia, DEVICE_MEMORY); }

template <typename IndexType, typename ValueType>
void delete_device_matrix(ell_matrix<IndexType,ValueType>& ell){ delete_ell_matrix(ell, DEVICE_MEMORY); }

template <typename IndexType, typename ValueType>
void delete_device_matrix(coo_matrix<IndexType,ValueType>& coo){ delete_coo_matrix(coo, DEVICE_MEMORY); }

template <typename IndexType, typename ValueType>
void delete_device_matrix(csr_matrix<IndexType,ValueType>& csr){ delete_csr_matrix(csr, DEVICE_MEMORY); }

template <class IndexType, class ValueType>
void delete_device_matrix(hyb_matrix<IndexType,ValueType>& hyb){  delete_hyb_matrix(hyb, DEVICE_MEMORY); }

template <typename IndexType, typename ValueType>
void delete_device_matrix(pkt_matrix<IndexType, ValueType>& pm){ delete_pkt_matrix(pm, DEVICE_MEMORY); }

////////////////////////////////////////////////////////////////////////////////
//! copy to device
////////////////////////////////////////////////////////////////////////////////

// generalize with:
// template <typename T>
// T * copy_to(const T *, const size_t n, const memory_location loc)

template <typename IndexType, typename ValueType>
dia_matrix<IndexType, ValueType> copy_matrix_to_device(const dia_matrix<IndexType, ValueType>& h_dia)
{
    dia_matrix<IndexType, ValueType> d_dia = h_dia; //copy fields
    d_dia.diag_offsets = copy_array_to_device(h_dia.diag_offsets, h_dia.num_diags);
    d_dia.diag_data    = copy_array_to_device(h_dia.diag_data,    h_dia.num_diags * h_dia.stride);
    return d_dia;
}

template <typename IndexType, typename ValueType>
ell_matrix<IndexType, ValueType> copy_matrix_to_device(const ell_matrix<IndexType, ValueType>& h_ell)
{
    ell_matrix<IndexType, ValueType> d_ell = h_ell; //copy fields
    d_ell.Aj = copy_array_to_device(h_ell.Aj, h_ell.stride * h_ell.num_cols_per_row);
    d_ell.Ax = copy_array_to_device(h_ell.Ax, h_ell.stride * h_ell.num_cols_per_row);
    return d_ell;
}

template <typename IndexType, typename ValueType>
csr_matrix<IndexType, ValueType> copy_matrix_to_device(const csr_matrix<IndexType, ValueType>& h_csr)
{
    csr_matrix<IndexType, ValueType> d_csr = h_csr; //copy fields
    d_csr.Ap = copy_array_to_device(h_csr.Ap, h_csr.num_rows + 1);
    d_csr.Aj = copy_array_to_device(h_csr.Aj, h_csr.num_nonzeros);
    d_csr.Ax = copy_array_to_device(h_csr.Ax, h_csr.num_nonzeros);
    return d_csr;
}

template <typename IndexType, typename ValueType>
coo_matrix<IndexType, ValueType> copy_matrix_to_device(const coo_matrix<IndexType, ValueType>& h_coo)
{
    coo_matrix<IndexType, ValueType> d_coo = h_coo; //copy fields
    d_coo.I = copy_array_to_device(h_coo.I, h_coo.num_nonzeros);
    d_coo.J = copy_array_to_device(h_coo.J, h_coo.num_nonzeros);
    d_coo.V = copy_array_to_device(h_coo.V, h_coo.num_nonzeros);
    return d_coo;
}

template <typename IndexType, typename ValueType>
hyb_matrix<IndexType, ValueType> copy_matrix_to_device(const hyb_matrix<IndexType, ValueType>& h_hyb)
{

    hyb_matrix<IndexType, ValueType> d_hyb = h_hyb; //copy fields
    d_hyb.ell = copy_matrix_to_device(h_hyb.ell);
    d_hyb.coo = copy_matrix_to_device(h_hyb.coo);
    return d_hyb;
}

template <class IndexType, class ValueType>
pkt_matrix<IndexType,ValueType> copy_matrix_to_device(const pkt_matrix<IndexType,ValueType>& h_pkt)
{
    pkt_matrix<IndexType,ValueType> d_pkt = h_pkt; //copy fields
    
    d_pkt.row_ptr  = copy_array_to_device(h_pkt.row_ptr, h_pkt.num_packets + 1);

    d_pkt.permute_old_to_new = copy_array_to_device(h_pkt.permute_old_to_new, h_pkt.num_rows);
    d_pkt.permute_new_to_old = copy_array_to_device(h_pkt.permute_new_to_old, h_pkt.num_rows);

    const IndexType num_threads = h_pkt.num_packets * h_pkt.threads_per_packet;

    const IndexType num_local_elements = h_pkt.packets.total_cycles * h_pkt.threads_per_packet; 
    d_pkt.packets.pos_start  = copy_array_to_device(h_pkt.packets.pos_start,  num_threads);
    d_pkt.packets.pos_end    = copy_array_to_device(h_pkt.packets.pos_end,    num_threads);
    
    d_pkt.packets.index_array = copy_array_to_device(h_pkt.packets.index_array, num_local_elements);
    d_pkt.packets.data_array  = copy_array_to_device(h_pkt.packets.data_array,  num_local_elements);
    
    d_pkt.coo = copy_matrix_to_device(h_pkt.coo);
    
    return d_pkt;   
}
