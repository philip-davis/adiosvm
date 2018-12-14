#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

/*
 * Function that calculates the starting offsets and counts for a 3D decomposition.
 * The 'starts' and 'counts' can be used in the adios2_define_variable call.
 * The decomposition is done such that the first (non-contiguous) dimension is decomposed first,
 * followed by the second dimension, and so on.
 *
 * Inputs:
 * Global array dimensions nx,ny,nz
 * Communicator size comm_size
 * Process rank 'rank'
 *
 * Output:
 * starting offsets for the rank in 'starts'
 * no. of elements in each direction for the rank in 'counts'
 *
 * This is a collective call. It must be called by all processes of the participating communicator.
 *
 * Assumptions and caveats:
 * The global array size is perfectly divisible by the communicator size.
 * Additionally, because of the way the dimensions are split,
 *  the dimension x must be perfectly divisible by the communicator size, and so on.
 */
void get_starts_counts_3d_decomp (size_t nx, size_t ny, size_t nz, size_t *starts, size_t *counts, int comm_size, int rank) {
    size_t yz_plane_idx, x_idx, y_idx, z_idx;
    size_t x_count, y_count, z_count;

    // Calculate the starting offsets
    // Strategy:
    // linearize the 3D array into a 1D array and divide it equally amongst all processes to get the local starting offset.
    // Then get the 3D index of this offset as follows:
    // First calculate which x-plane the element will be in. This will be its x index.
    // Then get the yz plane that the element will fall in, and then calculate its y and z indices.
 
    size_t global_arr_size = nx*ny*nz;
    size_t my_global_index = global_arr_size/comm_size * rank;
    x_idx = my_global_index/(nz*ny);
    yz_plane_idx = my_global_index - (nz*ny*x_idx);
    y_idx = yz_plane_idx/ny;
    z_idx = yz_plane_idx%ny;

    starts[0] = x_idx;
    starts[1] = y_idx;
    starts[2] = z_idx;

    // Calculate the counts
    // Split the first dimension (non-contigous dimension), then the second one, and so on.
    x_count = nx;
    y_count = ny;
    z_count = nz;

    x_count = nx/comm_size;
    if (x_count == 0) {
        x_count = 1;

        y_count = ny/(comm_size/nx);
        if (y_count == 0) {
            y_count = 1;
            z_count = nz/(comm_size /(nx*ny));

            if(z_count == 0) {
                fprintf(stderr, "ERROR: Cannot decompose this array. "
                        "Ensure no. of processes is <= total no. of elements in the array\n");
                exit(-1);
            }
        }
    }

    counts[0] = x_count;
    counts[1] = y_count;
    counts[2] = z_count;
}

