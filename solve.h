#ifndef _SOLVER_H_
#define _SOLVER_H_

#include "cholmod.h"
#include <stdbool.h>

#ifndef __cplusplus
extern "C" {
#endif

// Every DECOMPOSITION_FREQUENCY iterations, we recalculate the LDL decomposition
// from the A matrix.
#define DECOMPOSITION_FREQUENCY 1000

// triplet used to supply connection data
typedef cholmod_triplet solver_triplet_t;

typedef struct {
    cholmod_common *common;

    // lookup table: is_input[i] is true if group i is an input else it is false
    bool *is_input;
    bool *is_output;  // same as for input, but these are the output nodes
    size_t n_input_groups;
    size_t n_output_groups;
    size_t n_groups;  // use int over size_t for consistancy with cholmod
    size_t nnz;  // number of non-zero values in A

    // since the matrix doesn't contain input groups (knowns), we need mappings
    size_t *group_to_mat_map;  // mapping from group index to A matrix index
    size_t *mat_to_group_map;  // mapping from A matrix index to group index
    size_t *group_to_g_map;  // mapping from group index to G matrix index
    size_t n_mat;  // matrix size (should be n_groups - n_input_groups) size(A)

    // stuff needed to solve Ax = b, were x [volts] and b [amps] are vectors
    cholmod_sparse *A;  // the full conductance matrix
    cholmod_sparse *G;  // matrix of conductances on the right-hand-side
    cholmod_factor *LD;  // the LD decomposition
    cholmod_dense *b;  // vector of injected currents b = G v [known voltages]
    cholmod_sparse *b_set;  // pattern of groups which have injected currents
    cholmod_dense *x;  // vector of solved voltages
    cholmod_sparse *x_set;  // pattern of groups voltages with defined solutions

    // vectors of length N, if A is NxN, which are allocated for intermediate
    // steps in the solve process and reused by the cholmod api
    cholmod_dense *y;
    cholmod_dense *e;
} solver_state_t;

// warning: don't print large matrices
void print_sparse(solver_state_t *state, cholmod_sparse *A, const char* name);
#define TRIPLET_DEBUG_LIMIT ((size_t) 10)
void print_triplet(solver_state_t* state, solver_triplet_t* triplet, const char* name);

/**
 * Three-way insertion sort for sorting triplets by column then row.
 * 
 * Since I assume that the indicies will be some what in the correct order, or
 * prehaps even already ordered, using insertion sort isn't so bad. In most cases
 * we are O(n) time with O(1) extra memory. O(n^2) worst case time though.
 * 
 * If this function is becomming a drag, consider quicksort or mergesort.
 * */
void sort_triplet(solver_triplet_t* triplet);

/**
 * Allocate memory for a state object needed for cholmod operations, configure,
 * and allow chomod to initalise itself (clears pointers, stats, default params).
 * Must be called before anything else to obtain our state object.
 * 
 * Must only every be one state, because of cholmod_start/finish.
 * */
solver_state_t* solver_create_state();

/**
 * Compute the A matrix and b vector from a triplet and input vector.
 * 
 * This is prep for solving Ax = Gv = b, where x is a vector of group voltages
 * and b is a vector of injected currents obtained from the known voltages.
 * 
 * The input triplet is a list where each entry is a coordinate between two
 * groups and the conductance of that connection. The input triplet is expected
 * to be symmetric, only specify entries where row > col, and to specify NO
 * self conductance (ie no conduntance between the same group a->a). The network
 * must also be fully connected - ie have only a single component which consists
 * of more than one group. Some of these groups will be input groups.
 * 
 * The A matrix is a conductance matrix, stored as a sparse matrix in
 * compressed column form. All rows/columns corrisponding to input voltages
 * are moved to the right-hand-side (the b matrix), since the corrisponding
 * voltages are known at each timestep, and removed from A.
 * 
 * The b vector is a vector of injected currents such that the input voltages
 * are maintained at the specified groups. If no input voltages are specified,
 * the system has the trivial solution (x = 0). If only one input voltage is
 * specified, then the system has the solution x = a, where a is some voltage.
 * Hence, the number of input groups should be at-least 2.
 * 
 * Memory requirements:
 * 2*nnz for triplet
 * 1.5*(nnz + n_groups) for adj_c
 * 1.5*nnz for adj_r
 * triplet dealloced => -2*nnz
 * 1.5*nnz + n_mat for A
 * < 1.5*n_input_groups*n_mat + n_input_groups for G
 * 
 * 0.5 for int type (assuming 4 bytes) and 1 for double type (8 bytes)
 * 
 * Assuming nnz >> n_groups, n_mat, n_input_groups and the number of inputs is
 * small then the most memory needed at any one time is roughly:
 * => 8 * 5 * nnz (bytes)
 * but we probably don't need matrices big enough to care about this. Other
 * memory will be allocated into state for input/output lookup tables and
 * group-to-matrix index mappings (size proportonal to n_groups).
 * 
 * Time complexity:
 * O(nnz)
 * 
 * nnz: the number of non-zero elements.
 * 
 * @param state: where all our data/results are kept
 * @param triplet: a cholmod_triplet describing the network (this will get freed!)
 * @param input_groups: a vector of groups which have a known voltage
 * @param n_input_groups: number of input groups
 * @param output_groups: a vector of groups we want to know the voltage of
 * @param n_output_groups: number of output groups
 * 
 * @returns: 0 on completion, error otherwise
 * */
int solver_initalise_network
(
    solver_state_t* state,
    // a list of connections between groups and their conductances, there are
    // helper function to help construct this object. This triplet will be freed
    solver_triplet_t *triplet,  // Note: will sort the triplets
    const int *input_groups,  // array of indices which are input groups
    const size_t n_input_groups,  // number of input groups
    const int *output_groups,  // array of indices which are output groups
    const size_t n_output_groups   // number of output groups
);

/**
 * Deallocate memory and tell cholmod we are done.
 * Must be the very last call to this module.
 * */
void solver_destroy_state(solver_state_t* state);

// pointer to update function
typedef int (*update_func_t)(solver_state_t*, solver_triplet_t*, double**, int i, void* data);

/**
 * Solve a system where the input voltages can vary with time.
 * */
int solver_iterate_ac
(
    solver_state_t *state,
    const size_t n_iters,
    solver_triplet_t* connections,
    const int* input_groups,
    const size_t n_input_groups,
    const int* output_groups,
    const size_t n_outputs,
    update_func_t update_func,
    void* data
);

/**
 * Read in a .mtx file. The file must be in symmetric lower triangular real
 * coordinate form.
 * */
cholmod_triplet* triplet_from_file
(
    solver_state_t* state,
    // the filename to read from, as a string
    const char* filename
);

#ifndef __cplusplus
}
#endif
#endif
