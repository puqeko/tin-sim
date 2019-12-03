/* make_symmetric.cpp
 * 
 * Checks to see if a matrix, given in .mtx format to stdin, is stored as a
 * symmetric matrix in lower triangular form. If it is not, check if we can
 * safely convert to lower triangular symmetric form, do the conversion, then
 * save out to stdout.
 * 
 * Info is printed to stderr.
 * 
 * example use: convert the simple.mtx file to symmetric format
 * ./make_symmetric.o < matrices/simple.mtx > out.mtx
 */

#include <stdint.h>
#include <stdbool.h>
#include <stdlib.h>
#include <cholmod.h>
#include "sort_r.h"  // cause qsort_r works different between operating systems

#define LOWER_TRIANGULAR -1
#define UPPER_TRIANGULAR 1

#define min(a, b) (((a) < (b)) ? (a) : (b))

typedef struct {
    cholmod_common *common;

    // lookup table: is_input[i] is true if group i is an input else it is false
    bool *is_input;
    bool *is_output;  // same as for input, but these are the output nodes

    // TODO: test if we can use x to store input voltages + solutions togther
    // or is it not safe?
    // double *voltages;

    // a map between a group index and the index in A
    // ie j = group_to_mat_map[i] means that group i is at index j in A
    // we must do this since groups of known voltage are removed from A
    int *group_to_mat_map, *mat_to_group_map;

    // stuff needed to solve Ax = b, were x [volts] and b [amps] are vectors
    cholmod_sparse *A;  // the full conductance matrix
    cholmod_factor *LD;  // the LD decomposition
    cholmod_dense *g;  // vector of conductances on right-hand-side
    cholmod_dense *b;  // vector of injected currents b = g x [known voltages]
    cholmod_sparse *b_set;  // pattern of groups which have injected currents
    cholmod_dense *x;  // vector of solved voltages
    cholmod_sparse *x_set;  // pattern of groups voltages with defined solutions

    // vectors of length N, if A is NxN, which are allocated for intermediate
    // steps in the solve process and reused by the cholmod api
    cholmod_dense *y_handle;
    cholmod_dense *e_handle;
} solver_state_t;



/**
 * Allocate memory for a state object needed for cholmod operations, configure,
 * and allow chomod to initalise itself (clears pointers, stats, default params).
 * Must be called before anything else to obtain our state object.
 * */
solver_state_t* solver_create_state(void)
{
    solver_state_t* state;
    state = calloc(1, sizeof(solver_state_t));
    state->common = calloc(1, sizeof(cholmod_common));

    // use lower triangular in the default case when storing sparse matrices.
    state->common->prefer_upper = false;

    cholmod_start(state->common);
    return state;
}

/**
 * Deallocate memory and tell cholmod we are done.
 * Must be the very last call to this module.
 * */
void solver_destroy_state(solver_state_t* state)
{
    // deallocate all pointed heap memory
    if (state->is_input) free(state->is_input);
    if (state->is_output) free(state->is_output);
    if (state->group_to_mat_map) free(state->group_to_mat_map);
    if (state->mat_to_group_map) free(state->mat_to_group_map);
    if (state->A) cholmod_free_sparse(&state->A, state->common);
    if (state->LD) cholmod_free_factor(&state->LD, state->common);
    if (state->g) cholmod_free_dense(&state->b, state->common);
    if (state->b) cholmod_free_dense(&state->b, state->common);
    if (state->b_set) cholmod_free_sparse(&state->b_set, state->common);
    if (state->x) cholmod_free_dense(&state->x, state->common);
    if (state->x_set) cholmod_free_sparse(&state->x_set, state->common);
    if (state->y_handle) cholmod_free_dense(&state->y_handle, state->common);
    if (state->e_handle) cholmod_free_dense(&state->e_handle, state->common);

    // finish up cholmod and free state
    cholmod_finish(state->common);
    free(state->common);
    free(state);
}


// Used in insertion sort for triplets to return row followed by column ordering
// Returns -1 if left < right, 1 if left > right, else 0
int _triplet_comparitor(cholmod_triplet *triplet, const int l, const int r)
{
    // sort out the types of our inputs
    // the reason this code looks like cancer is because we must cast all the
    // void pointers to types so that c knows how to deal with them
    // this happens for triplet properties and for this function's arguments

    if (((int*)triplet->j)[l] < ((int*)triplet->j)[r]) return -1;
    else if (((int*)triplet->j)[l] > ((int*)triplet->j)[r]) return 1;
    else {
        // the column number is the same, so compare row
        if (((int*)triplet->i)[l] < ((int*)triplet->i)[r]) return -1;
        else if (((int*)triplet->i)[l] > ((int*)triplet->i)[r]) return 1;
    }
    return 0;  // default
}


/**
 * Three way insertion sort for sorting triplets by row then column.
 * 
 * Since I assume that the indicies will be some what in the correct order, or
 * prehaps even already ordered, using insertion sort isn't so bad. In most cases
 * we are O(n) time with O(1) extra memory. O(n^2) worst case time.
 * */
void sort_triplet(cholmod_triplet* triplet)
{
    // start with the first item 'sorted' then find a place to put the item at
    // position next.
    for (int next = 1; next < triplet->nnz; next++) {
        printf("next: %d\n", next);
        for (int i = next - 1; i >= 0; i--) {
            printf("%d\n", i);
            if (_triplet_comparitor(triplet, i, i+1) > 0) {
                int tempi = ((int*)triplet->i)[i];
                int tempj = ((int*)triplet->j)[i];
                double tempx = ((double*)triplet->x)[i];

                ((int*)triplet->i)[i] = ((int*)triplet->i)[i+1];
                ((int*)triplet->j)[i] = ((int*)triplet->j)[i+1];
                ((double*)triplet->x)[i] = ((double*)triplet->x)[i+1];

                ((int*)triplet->i)[i+1] = tempi;
                ((int*)triplet->j)[i+1] = tempj;
                ((double*)triplet->x)[i+1] = tempx;
            }
            else break;  // otherwise, these are sorted so move to next one
        }
    }
}

/**
 * Compute the A matrix and b vector from a triplet and input vector.
 * 
 * This is prep for solving Ax = b, where x is a vector of group voltages and b
 * is a vector of injected currents.
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
 * About 2*nnz + 3*nrows memory is needed, where nnz is the number of non-zero
 * entries (excluding globals).
 * Other memory will be allocated into state if needed.
 * 
 * @param state: where all our data/results are kept
 * @param triplet: a cholmod_triplet describing the network (this will get freed!)
 * @param input_groups: a vector of groups which have a known voltage
 * 
 * @returns: 0 on success, error otherwise
 * */
int solver_initalise_network
(
    solver_state_t* state,
    // a list of connections between groups and their conductances, there are
    // helper function to help construct this object. This triplet will be freed
    cholmod_triplet *triplet,
    const int *input_groups,  // array of indices which are input groups
    const int n_input_groups,  // number of input groups
    const int *output_groups,  // array of indices which are output groups
    const int n_output_groups   // number of output groups
)
{
    if (n_output_groups < 1) {
        fprintf(stderr, "If you don't specify any output voltages, then I have "
                "nothing to calculate. Try adding some output groups.\n");
        return -1;
    }

    // see main comment above
    if (n_input_groups < 2) {
        fprintf(stderr, "Warning: You probably intended to specify more than one "
                "input voltage. All the voltages will be the same otherwise. "
                "Try adding some input groups.\n");
        // still try it though, just a warning
    }
    
    // make some assertions about symmetry, type, etc of triplet
    if (triplet->ncol != triplet->nrow || triplet->stype != LOWER_TRIANGULAR ||
        triplet->xtype != CHOLMOD_REAL) {
        fprintf(stderr, "Triplet must be symmetric, lower triangular, and of type"
                " float.\nDEBUG INFO: %ldx%ld mat, stype:%d, xtype:%d\n",
                triplet->ncol, triplet->nrow, triplet->stype, triplet->xtype);
        return -1;
    }

    size_t n_groups = triplet->nrow;  // ie, the total number of groups

    // make set of input and output groups for faster O(1) lookup
    if (!state->is_input)
        state->is_input = calloc(n_groups, sizeof(*state->is_input));
    for (int i = 0; i < n_input_groups; i++) {
        // assume that the group indices are not out of bounds 
        state->is_input[input_groups[i]] = true;
    }
    if (!state->is_output)
        state->is_output = calloc(n_groups, sizeof(*state->is_output));
    for (int i =0; i < n_output_groups; i++) {
        state->is_output[output_groups[i]] = true;
    }

    // map from group index to matrix index
    state->group_to_mat_map = calloc(n_groups, sizeof(int));
    state->mat_to_group_map = calloc(n_groups, sizeof(int));

    // lets form the sparse matrix A in compressed column form...
    // cholmod_sparse* A = cholmod_allocate_sparse(
    //     n_groups, n_groups, triplet->nnz + n_groups, true, true,
    //     LOWER_TRIANGULAR, CHOLMOD_REAL, state->common);

    // make sure triplet is sorted by row then column in accending order
    sort_triplet(triplet);
    return 0;
}

/**
 * Read in a .mtx file. The file must be in symmetric lower triangular real
 * coordinate form.
 * */
cholmod_triplet* triplet_from_file
(
    solver_state_t* state,
    // the filename to read from, as a string
    const char* filename
)
{
    // read file into triplet form
    FILE *file = fopen(filename, "r");
    cholmod_triplet* conductance_coords = cholmod_read_triplet(file, state->common);

    // TODO: check if real && coords && lower triangular

    return conductance_coords;
}

// /**
//  * Compute the A matrix and b vector from a triplet and input vector.
//  * 
//  * This is prep for solving Ax = b, where x is a vector of group voltages and b
//  * is a vector of injected currents.
//  * 
//  * The input triplet is a list where each entry is a coordinate between two
//  * groups and the conductance of that connection. The input triplet is expected
//  * to be symmetric, only specify entries where row > col, and to specify NO
//  * self conductance (ie no conduntance between the same group a->a). The network
//  * must also be fully connected - ie have only a single component which consists
//  * of more than one group. Some of these groups will be input groups.
//  * 
//  * The A matrix is a conductance matrix, stored as a sparse matrix in
//  * compressed column form. All rows/columns corrisponding to input voltages
//  * are moved to the right-hand-side (the b matrix), since the corrisponding
//  * voltages are known at each timestep, and removed from A.
//  * 
//  * The b vector is a vector of injected currents such that the input voltages
//  * are maintained at the specified groups. If no input voltages are specified,
//  * the system has the trivial solution (x = 0). If only one input voltage is
//  * specified, then the system has the solution x = a, where a is some voltage.
//  * Hence, the number of input groups should be at-least 2.
//  * 
//  * About 2*nnz + 3*nrows memory is needed, where nnz is the number of non-zero
//  * entries (excluding globals)
//  * 
//  * @param state: where all our data/results are kept
//  * @param triplet: a cholmod_triplet describing the network (this will get freed!)
//  * @param input_groups: a vector of groups which have a known voltage
//  * 
//  * @returns: 0 on success, error otherwise
//  * */
// int solver_initalise_network2
// (
//     solver_state_t* state,
//     // a list of connections between groups and their conductances, there are
//     // helper function to help construct this object. This triplet will be freed
//     const cholmod_triplet *triplet,
//     // the indices of groups that we will know the voltage of (ie group number)
//     const int *input_groups,
//     const int n_input_groups,
//     const int *output_groups
// )
// {
//     // TODO: potential optimisation. instead of creating a new_triplet structure and
//     // converting to a sparse matrix, one could compute the compressed column
//     // form of the sparse matrix directly by first sorting the entries in
//     // triplet then adding only the approporate entries to the structure. This
//     // reduces the amount of copying needed.

//     if (n_input_groups == 0)
//     {
//         fprintf(stderr, "Must call set_input_groups(...) first.\n");
//         return -1;
//     }

//     // see main comment above
//     if (n_input_groups < 2) {
//         fprintf(stderr, "Warning: You should specify more than two input voltages. "
//                      "The system isn't going to be very interesting otherwise.\n");
//         // still try it though, just a warning
//     }

//     // make some assertions about symmetry, type, etc of triplet
//     if (triplet->ncol != triplet->nrow || triplet->xtype != CHOLMOD_DOUBLE ||
//         triplet->stype != LOWER_TRIANGULAR) {
//         fprintf(stderr, "Triplet must be symmetric, lower triangular, and of type"
//                      " float.\n");
//         return -1;
//     }

//     // make set of input groups for faster O(1) lookup
//     int *input_node_lookup = calloc(triplet->nrow, sizeof(int));
//     for (int i = 0; i < n_input_groups; i++) {
//         input_node_lookup[input_groups[i]] = true;
//     }

//     // map from group index to matrix index
//     state->group_to_mat_map = calloc(triplet->nrow, sizeof(int));
//     state->mat_to_group_map = calloc(triplet->nrow, sizeof(int));

//     // build new triplet to construct a conductance matrix A
//     // assuming lower triangular form (row > col) with no diagonal entries, then
//     // we must add the self conductance terms to the existing matrix (the diagonls)
//     // and remove entries corrisponding to know input votlages.

//     // allocate space to store our now coordinates
//     // use the nnz + extra diagonal entries as a starting estimate for the memory
//     // we need. It will be less than this since we remove some entries.
//     int len_est = triplet->nnz + triplet->nrow;

//     cholmod_triplet* new_triplet;
//     new_triplet = cholmod_allocate_triplet(triplet->nrow, triplet->nrow, len_est,
//                                            LOWER_TRIANGULAR, CHOLMOD_DOUBLE, state->common);

//     // total self conductance (diagonal entries)
//     double* g_total = calloc(triplet->nrow, sizeof(double));
//     double* b_temp = calloc(triplet->nrow, sizeof(double));  // b vector
//     // By kirchoffs current law, the sum of currents in/out of a group is zero
//     // assuming no injected current. g_total is initalised to zero also.

//     int t = 0;  // the current index in the new_triplet
//     for (int i = 0; i < triplet->nnz; i++) {
//         // get row, col coords and conductance value of current triple
//         int ri = ((int*)triplet->i)[i];
//         int ci = ((int*)triplet->j)[i];
//         double gi = ((double*)triplet->x)[i];

//         // get total conductance into/out-of a group for the diagonals of
//         // the conductance matrix we are building.
//         g_total[ri] += gi;
//         g_total[ci] += gi;

//         // if the row or column belongs to our known input groups then
//         // adjust b matrix, ie our injected currents, accordingly.
//         // Remember, the only the lower triangular tripples are given.
//         if (state->is_input[ri]) {
//             b_temp[ci] += gi;  // move -gi (upper triangular term) to other side
//         } else if (state->is_input[ci]) {
//             b_temp[ri] += gi;  // move -gi (lower triangular term) to other side
//         } else {
//             // this isn't a known input row/col
//             ((int*)new_triplet->i)[t] = ri;
//             ((int*)new_triplet->j)[t] = ci;
//             ((double*)new_triplet->x)[t] = -gi;  // -ive off diagonals
//             t++;  // this will lag behind i
//         }
//     }

//     // add diagonal self conductance entries
//     for (int j = 0; j < new_triplet->nrow; j++) {
//         if (!state->is_input[j]) {
//             // ignoring input nodes
//             ((int*)new_triplet->i)[t] = j;
//             ((int*)new_triplet->j)[t] = j;
//             ((double*)new_triplet->x)[t] = g_total[j];

//             // map between matrix coords and actual group indices
//             state->group_to_mat_map[j] = t;
//             state->mat_to_group_map[t] = j;
//             t++;
//         }
//     }
//     cholmod_free_triplet(&triplet, state->common);
//     free(g_total);
//     // resize now that we know exactly how big new_triplet is
//     cholmod_reallocate_triplet(t, new_triplet, state->common);

//     // convert A to sparse
//     cholmod_sparse *A;
//     A = cholmod_triplet_to_sparse(new_triplet, new_triplet->nnz, state->common);
//     cholmod_free_triplet(&new_triplet, state->common);
//     state->A = A;  // save pointer to state

//     // convert b to sparse
//     cholmod_triplet *b_triplet;
//     b_triplet = cholmod_allocate_triplet(new_triplet->nrow, new_triplet->nrow,
//                                          new_triplet->nrow, 0, CHOLMOD_DOUBLE, state->common);

//     // populate with non-zero b entries
//     int u = 0;  // b_triplet index
//     for (int k = 0; k < b_triplet->nrow; k++) {
//         if (b_temp[k] != 0.0) {
//             ((int*)b_triplet->i)[u] = k;
//             ((int*)b_triplet->j)[u] = 1;  // only one column
//             ((double*)b_triplet->x)[u] = b_temp[k];
//             u++;
//         }
//     }
//     cholmod_reallocate_triplet(u, b_triplet, state->common);  // shrink to right size
//     free(b_temp);

//     // convert b to sparse
//     // we do this not because it saves memory, but because it reduces the
//     // number of computations required every time we multiply an element of b.
//     cholmod_sparse *b;
//     b = cholmod_triplet_to_sparse(b_triplet, u, state->common);
//     cholmod_free_triplet(&b_triplet, state->common);
//     state->b = b;  // save to state

//     // TODO: return object containing relevant stuff.
//     return 0;
// }

// https://stackoverflow.com/questions/2744181/how-to-call-c-function-from-c


#define TRIPLET_DEBUG_LIMIT 10
void print_triplet(solver_state_t* state, cholmod_triplet* triplet)
{
    cholmod_print_triplet(triplet, NULL, state->common);

    // If TRIPLET_DEBUG_LIMIT == -1, print all triplet entries. Otherwise, just
    // print up to TRIPLET_DEBUG_LIMIT number of entries
    int n_print;
    if (TRIPLET_DEBUG_LIMIT > 0) n_print = min(triplet->nnz, TRIPLET_DEBUG_LIMIT);
    else n_print = triplet->nnz;

    // print each entry
    for (int i =0; i < n_print; i++) {
        printf("%d %d %f\n", ((int*)triplet->i)[i], ((int*)triplet->j)[i], ((double*)triplet->x)[i]);
    }

    // Show ... if we skipped some entires
    if (n_print < triplet->nnz) printf("...\n\n");
    else printf("\n");
}


// Load a real simple example matrix and solve.
// We're using the folloiwng a 4x4 conductance matrix. The interconductances
// are -ve and the mutual conductances are +ve as per convention.
//
//       column
//         0:    1:   2:   3:
//  row   ----------------------
//   0:  |  0.7  -0.4 -0.1 -0.2
//   1:  | -0.4   0.7 -0.3  0
//   2:  | -0.1  -0.3  0.9 -0.5
//   3:  | -0.2   0   -0.5  0.7
//
void testcase_0()
{
    // initalise module
    solver_state_t* state = solver_create_state();

    // Load the above matrix in triplet form, as below
    // 
    // from | to  | conductance of connection
    //-------------------------
    // 0    | 1   | 0.4
    // 0    | 2   | 0.1
    // 1    | 2   | 0.3
    // 0    | 3   | 0.2
    // 2    | 3   | 0.5
    //
    // Notice that only the entries with to > from are specified, since the
    // matrix is symmetric (connections go both ways).
    cholmod_triplet *triplet = triplet_from_file(state, "matrices/testcase_0.mtx");
    print_triplet(state, triplet);

    // Set which indices we want to be inputs and outputs. This allows us to
    // deterwmine which voltages are known and unknown and which solutions we
    // care about
    int input_groups[] = {1, 2}; int n_inputs = 2;  // set v_1 and v_2 as known
    int output_groups[] = {0, 1, 2, 3}; int n_ouputs = 4;  // get all voltages

    // // Initalise the solver. This will construct the A matrix and b vector used
    // // to represent the system.
    solver_initalise_network(state, triplet,
                             &input_groups[0], n_inputs,
                             &output_groups[0], n_ouputs);

    print_triplet(state, triplet);

    solver_destroy_state(state);
    // state can no longer be used
}


int main ()
{
    testcase_0();
}
