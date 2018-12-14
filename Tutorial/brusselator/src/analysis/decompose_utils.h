#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#ifndef __DECOMP_UTILS__
#define __DECOMP_UTILS__

void get_starts_counts_3d_decomp (size_t nx, size_t ny, size_t nz, size_t *starts, size_t *counts, int comm_size, int rank);

#endif

