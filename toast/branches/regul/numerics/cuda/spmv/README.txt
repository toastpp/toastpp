Sparse matrix-vector multiplication (SpMV) CUDA kernels.

Please see LICENSE.txt for full source license.

Description
-----------

The SpMV operation is 
    y += A*x, 
where A is a sparse matrix and x and y are column vectors.  

Kernels and Formats
-------------------

The following SpMV kernels are provided:
 - csr_scalar : simple assignment of 1 thread per row for matrix in CSR format
 - csr_vector : 1 warp per row for matrix in CSR format
 - coo_flat : a segmented reduction of matrix in a sorted COO format
 - ell : 1 thread per row for a matrix in ELLPACK/ITPACK format
 - dia : 1 thread per row for a matrix in DIA format
 - pkt : kernel for the "packet" format

The sparse matrix storage formats are:
 - CSR : Compressed Sparse Row 
 - COO : Coordinate (aka Triplet or IJV format)
 - ELL : ELLPACK/ITPACK (fixed number of nz per row, zero padded)
 - DIA : Sparse diagonal format (e.g. spdiags() in MATLAB)
 - HYB : hybrid ELL/COO combination
 - PKT : Packet format

Other than HYB and PKT, these are all fairly standard formats.  For 
references consider [2] and Section 3.4 of [3].  Note that DIA and ELL can't 
efficiently represent all matrices, so they will print an error if 
the storage required exceeds a predefined threshold. HYB attempts 
to store the majority of the matrix in ELL format (faster) and only 
uses COO format (slower) for the parts of the matrix that aren't 
efficiently represented in ELL.

All kernels have a _tex variant that uses texture memory for gathering
values from the x vector.  On reasonably ordered matrices this can produce 
a noticeable improvement.  We implement texture<double> using texture<int2> 
as shown in texture.h.

Refer to the following papers for additional details (see [6,7] below):
 	"Efficient Sparse Matrix-Vector Multiplication on CUDA"
    Nathan Bell and Michael Garland, "NVIDIA Technical Report NVR-2008-004", December 2008

  	"Implementing Sparse Matrix-Vector Multiplication on Throughput-Oriented Processors"
    Nathan Bell and Michael Garland, in "Proc. Supercomputing '09", November 2009


Build Instructions
------------------

The following commands should compile the driver program:

  For GT200 processors (Tesla 1060 and Geforce 200 series)
    $ nvcc  -arch=sm_13 -O3 -o spmv driver.cu mmio.c
  On G8x and G9x processors (e.g. Geforce 8800 GT)
    $ nvcc  -arch=sm_11 -O3 -o spmv driver.cu mmio.c
  On G80 processors
    $ nvcc  -arch=sm_10 -O3 -o spmv driver.cu mmio.c

This code has been tested extensively under CUDA 2.3 and CUDA 3.0b on 64-bit Linux.

Hardware Assumptions
--------------------

These kernels have been tested on primarly GT200.  Note that unlike a previous
release of this code, the COO kernel *does not* require the use of atomic 
operations (atomicCAS). All kernels support double precision and the 
use of single vs. double precision is determined by template parameter.  Other 
data types such as integers could be supported by adding the appropriate support
in texture.h.

The MAX_THREADS macro in kernels/spmv_common_device.cu.h is hard coded for the 
maximum number of threads on a full GT200 chip.  For best performance on
on devices with fewer SMs or fewer registers per SM this number would require 
some adjustment. Presently, only the coo_flat and csr_vector kernels make use 
MAX_THREADS.

The csr_vector kernel requires GT200 coalescing rules for reasonable performance. 
Performance on G8x and G9x will be considerably lower.  The DIA, ELL, COO, and PKT
formats agree with G80 coalescing rules, so their performance will be comparably
better.

Compiling with -arch=sm_10 or -arch=sm_11 should work, but double precision 
support will be disabled.

Program Usage
-------------

The program is executed using,

  $ ./spmv example.mtx

where example.mtx is the file name of a sparse matrix in MatrixMarket format [1].
For reference, there are several small examples of sparse matrices in the examples/ 
sub-directory.  The kernels are not written for best performance for small matrices 
(e.g. NNZ < 500k), but it is possible to do so.  The UF Sparse Matrix Collection
is a good source for larger examples [5].  The collection of unstructured matrices
used in the report "Efficient Sparse Matrix-Vector Multiplication on CUDA" is
available online [6].  If no matrix file is provided in the arguments to the 
program, a simple example is generated on the fly.

The program also supports the following optional arguments  
  $ ./spmv example.mtx --device=1            # any valid device number
  $ ./spmv example.mtx --precision=64        # use 32-bit (single) or 64-bit (double) FP
  $ ./spmv example.mtx --partition=part.txt  # where part.txt is a partition file

The partition file is required for use of the PKT matrix format.  A partition file
must contain a partition number for each row of the (square) matrix.  For instance,
sample.txt below defines a partitioning of a 5-by-5 matrix into 3 parts.  In
this example, the first partition (index 0) contains the first and last row 
of the matrix.

------- part.txt -------
0
2
1
1
0
------------------------

The included partitioning code (partition.cxx and partition.h) require the METIS
partitioning library.  See metis/README.txt for additional information.  With
a compiled METIS library residing in spmv/metis/, the partitioner is compiled
as follows:

    $ gcc -o mmio.o -c -O3 mmio.c
    $ g++ -o partition.o -c -O3 -Wall -I. -Imetis/include partition.cxx
    $ g++ -o partition partition.o mmio.o -Lmetis/build/Linux-x86_64 -lmetis


Source Code
-----------

All source code was written by Nathan Bell and Michael Garland at NVIDIA, 
except for the MatrixMarket IO files (mmio.c and mmio.h) which are provided
by NIST [1].


Examples
--------

There are several example matrices in the examples/ sub-directory.  Note that these
are trivially small, so they aren't useful for performance testing.  Nathan
Bell (nbell@nvidia.com) can provide you with a more comprehensive set of 
examples with 500K-10M nonzero values.  For instance, we observe very competitive
performance on the matrices used in [4] which are available online [6].


Performance
-----------

We measure performance in GFLOP/s where:
   GFLOP/s = (2 * nonzeros(A)) / (time of one SpMV)
The factor of two in the numerator is due the the fact that each nonzero costs
exactly one multiply and one add (i.e. y[i] += A[i,j] * x[j]).

The _effective_ memory bandwidth used by each kernel is also outputted.  This
performance metric is explained in references [6,7].

Further Information
-------------------

A tech report with a comprehensive discussion of this code is available 
online [6], as well as a subsequent Supercomputing '09 paper [7]. Feel
free to contact Nathan Bell (nbell@nvidia.com) with any questions or comments.

Note that kernels are also available in the open-source Cusp library [8].

[1] MatrixMarket file format: http://math.nist.gov/MatrixMarket/formats.html#MMformat
[2] SPARSKIT documentation: http://www-users.cs.umn.edu/~saad/software/SPARSKIT/sparskit.html
[3] Saad's Book : http://www-users.cs.umn.edu/~saad/books.html
[4] http://www.cs.berkeley.edu/~samw/research/papers/sc07.pdf
[5] http://www.cis.ufl.edu/research/sparse/matrices/Boeing/pwtk.html 
[6] http://www.nvidia.com/object/nvidia_research_pub_001.html
[7] http://www.nvidia.com/object/nvidia_research_pub_013.html
[8] http://code.google.com/p/cusp-library/

