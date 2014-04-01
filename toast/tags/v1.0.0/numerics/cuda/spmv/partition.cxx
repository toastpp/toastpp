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



#include "sparse_io.h"
#include "partition.h"

#include <string.h>

#define SMEM_SIZE 16384 //size of shared memory in bytes

template <typename IndexType, typename ValueType>
              
void partition_matrix(const char * mm_filename, 
                      const char * output_filename, 
                      bool Kway)
{
    csr_matrix<IndexType,ValueType> csr = read_csr_matrix<IndexType,ValueType>(mm_filename);
	printf("Read matrix (%d,%d) with %d nonzeros\n", csr.num_rows, csr.num_cols, csr.num_nonzeros);

	if(csr.num_rows != csr.num_cols){
		printf("matrix must be square\n");
		exit(EXIT_FAILURE);
	}

    const IndexType rows_per_packet = (sizeof(ValueType) == 4) ? 1940 : 970; // METIS k-way guarantees +/- 3%
    int num_parts = (csr.num_rows + rows_per_packet - 1) / rows_per_packet;  //number of partitions

	std::vector<IndexType> partition(csr.num_rows);

	printf(" [partitioning matrix:");
	partition_csr<IndexType,ValueType>(csr, num_parts, partition, Kway);
	printf(" done]\n");
    
    printf("writing partition to file: %s\n", output_filename);

    FILE * fid = fopen(output_filename, "w");
    for(IndexType i = 0; i < partition.size(); i++)
        fprintf(fid, "%d\n", partition[i]);
    fclose(fid);

    delete_host_matrix(csr);
}



void usage(int argc, char** argv)
{
	printf("Usage:\n");
	printf("\t%s  my_matrix.mtx\n\n",argv[0]);
	printf("Note: my_matrix.mtx must be real-valued sparse matrix in the MatrixMarket file format.\n");	
}


int main(int argc, char** argv)
{
	if (argc != 2)
    {
		usage(argc,argv);
        return EXIT_FAILURE;
	}

    // Use standard recursive METIS partitioner
    //{
    //    const char * single_suffix = ".partition.2way.fp32";
    //    char * single_filename = (char *) malloc(sizeof(char) * (strlen(argv[1]) + strlen(single_suffix) + 1));
    //    strcpy(single_filename, argv[1]);
    //    strcat(single_filename, single_suffix);

    //    partition_matrix<unsigned int, float>(argv[1], single_filename, false);
    //    
    //    const char * double_suffix = ".partition.2way.fp64";
    //    char * double_filename = (char *) malloc(sizeof(char) * (strlen(argv[1]) + strlen(double_suffix) + 1));
    //    strcpy(double_filename, argv[1]);
    //    strcat(double_filename, double_suffix);
    //    
    //    partition_matrix<unsigned int, double>(argv[1], double_filename, false);
    //}

    // Use k-Way METIS partitioner
    {
        const char * single_suffix = ".partition.kway.fp32";
        char * single_filename = (char *) malloc(sizeof(char) * (strlen(argv[1]) + strlen(single_suffix) + 1));
        strcpy(single_filename, argv[1]);
        strcat(single_filename, single_suffix);

        partition_matrix<unsigned int, float>(argv[1], single_filename, true);

        const char * double_suffix = ".partition.kway.fp64";
        char * double_filename = (char *) malloc(sizeof(char) * (strlen(argv[1]) + strlen(double_suffix) + 1));
        strcpy(double_filename, argv[1]);
        strcat(double_filename, double_suffix);

        partition_matrix<unsigned int, double>(argv[1], double_filename, true);
    }

	return EXIT_SUCCESS;
}
