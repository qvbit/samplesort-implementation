#include <mpi.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdio.h>
#include <string.h>

#include "terarec.h"


int binary_search(terarec_t *local_data, terarec_t *local_splitters, int l, int r, int j) {
    int m;

    while (l < r) {
        m = l + (r-l)/2;
        if (teraCompare(&local_data[m], &local_splitters[j]) < 0) {
            l = m+1;
        }
        else {
            r = m;
        }
    }
    return l;
}


void terasort(terarec_t *local_data, int  local_len, 
			  terarec_t **sorted_data, int* sorted_counts, long* sorted_displs){		  
	int rank, P;
	int root = 0;
	MPI_Comm_size (MPI_COMM_WORLD, &P);
	MPI_Comm_rank (MPI_COMM_WORLD, &rank);

	qsort(local_data, local_len, sizeof(terarec_t), teraCompare);

	int nsplitters = P-1;
    int i;
	terarec_t* local_splitters;

	local_splitters = (terarec_t*)malloc(sizeof(terarec_t)*nsplitters);
	for (i=0; i < nsplitters; i++) {
		local_splitters[i] = local_data[(local_len/P)*(i+1)];
	}

	int ntotalsplitters = P*nsplitters;
	terarec_t *all_splitters;

	all_splitters = (terarec_t*)malloc(sizeof(terarec_t)*ntotalsplitters);
	MPI_Gather(local_splitters, nsplitters, mpi_tera_type,
			   all_splitters, nsplitters, mpi_tera_type,
			   root, MPI_COMM_WORLD);

	if (rank == root) {
		qsort(all_splitters, ntotalsplitters, sizeof(terarec_t), teraCompare);
		
		for (int i = 0; i < nsplitters; i++) {
			local_splitters[i] = all_splitters[nsplitters*(i+1)];
		}
	}

	MPI_Bcast(local_splitters, nsplitters, mpi_tera_type, root, MPI_COMM_WORLD);

    int j, curr, prev;
    int* scounts;

    scounts = (int*) malloc(sizeof(int)*P);
    for (i = 0; i < P; i++) {
        scounts[i] = 0;
    }

    // for (i=j=0; i < local_len; i++) {
    //     if (j == nsplitters || teraCompare(&local_data[i], &local_splitters[j]) < 0) {
    //         scounts[j]++;
    //     }
    //     else {
    //         j++;
    //         i--;
    //     }
    // }

    prev = 0;
    for (j = 0; j < P-1; j++) {
        curr = binary_search(local_data, local_splitters, prev, local_len-1, j);
        scounts[j] = curr-prev;
        prev = curr;
    }
    scounts[P-1] = local_len-prev;

    int *sdispls;

    sdispls = (int*)malloc(P*sizeof(int));
    sdispls[0] = 0;
	for (i = 1; i < P; i++) {
		sdispls[i] = sdispls[i-1] + scounts[i-1];
	}

    int nsorted;
    int *rcounts, *rdispls;

    rcounts = (int*) malloc(sizeof(int)*P);
    MPI_Alltoall(scounts, 1, MPI_INT,
                 rcounts, 1, MPI_INT,
                 MPI_COMM_WORLD);

    rdispls = (int*)malloc(P*sizeof(int));
    rdispls[0] = 0;
    for (i = 1; i < P; i++) {
        rdispls[i] = rdispls[i-1] + rcounts[i-1];
    }

    nsorted = rdispls[P-1] + rcounts[i-1];

    *sorted_data = (terarec_t*)malloc(sizeof(terarec_t)*nsorted);

    MPI_Alltoallv(local_data, scounts, sdispls, mpi_tera_type,
                  *sorted_data, rcounts, rdispls, mpi_tera_type,
                  MPI_COMM_WORLD);
    
    qsort((*sorted_data), nsorted, sizeof(terarec_t), teraCompare);

	MPI_Allgather(&nsorted, 1, MPI_INT,
				  sorted_counts, 1, MPI_INT,
				  MPI_COMM_WORLD);

	sorted_displs[0] = 0;
	for (int i = 1; i < P; i++) {
		sorted_displs[i] = sorted_displs[i-1] + sorted_counts[i-1];
	}

    free(local_splitters);
    free(all_splitters);
    free(scounts);
    free(sdispls);
    free(rcounts);
    free(rdispls);
}
