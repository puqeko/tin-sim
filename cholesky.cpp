/* cholesky.cpp */

#include <stdint.h>
#include <stdbool.h>
#include <stdlib.h>
#include <cholmod.h>
#include "sort_r.h"  // cause qsort_r works different between operating systems

#include <iostream>
#include <vector>
#include <tuple>

#define LOWER_TRIANGULAR -1
#define UPPER_TRIANGULAR 1

#define min(a, b) (((a) < (b)) ? (a) : (b))

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
    size_t *group_to_mat_map;  // mapping from group index to matrix index
    size_t *mat_to_group_map;  // mapping from matrix index to group index
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

// a vector or vectors each containing pairs (the row index and conductance)
typedef std::pair<size_t, double> link_t;
typedef std::vector<std::vector<link_t> > adj_list_t;

/**
 * Allocate memory for a state object needed for cholmod operations, configure,
 * and allow chomod to initalise itself (clears pointers, stats, default params).
 * Must be called before anything else to obtain our state object.
 * 
 * Must only every be one state, because of cholmod_start/finish.
 * */
solver_state_t* solver_create_state()
{
    solver_state_t* state = new solver_state_t();
    state->common = new cholmod_common();

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
    if (state->is_input) delete state->is_input;
    if (state->is_output) delete state->is_output;
    if (state->group_to_mat_map) delete state->group_to_mat_map;
    if (state->mat_to_group_map) delete state->mat_to_group_map;
    if (state->A) cholmod_free_sparse(&state->A, state->common);
    if (state->LD) cholmod_free_factor(&state->LD, state->common);
    if (state->G) cholmod_free_sparse(&state->G, state->common);
    if (state->b) cholmod_free_dense(&state->b, state->common);
    if (state->b_set) cholmod_free_sparse(&state->b_set, state->common);
    if (state->x) cholmod_free_dense(&state->x, state->common);
    if (state->x_set) cholmod_free_sparse(&state->x_set, state->common);
    if (state->y) cholmod_free_dense(&state->y, state->common);
    if (state->e) cholmod_free_dense(&state->e, state->common);

    cholmod_finish(state->common);

    // finish up cholmod and free state
    delete state->common;
    delete state;
}


// Used in insertion sort for triplets to return column followed by row ordering
// Returns -1 if left < right, 1 if left > right, else 0
int _triplet_comparitor(cholmod_triplet *triplet, const int &l, const int &r)
{
    // get pointers to triplet storage elements
    int* ti = (int*)triplet->i;  // row index
    int* tj = (int*)triplet->j;  // col index

    // sort out the types of our inputs
    // the reason this code looks like cancer is because we must cast all the
    // void pointers to types so that c knows how to deal with them
    // this happens for triplet properties and for this function's arguments
    if (tj[l] < tj[r]) return -1;
    else if (tj[l] > tj[r]) return 1;
    else {
        // the column number is the same, so compare row
        if (ti[l] < ti[r]) return -1;
        else if (ti[l] > ti[r]) return 1;
    }
    return 0;  // default
}


/**
 * Three-way insertion sort for sorting triplets by column then row.
 * 
 * Since I assume that the indicies will be some what in the correct order, or
 * prehaps even already ordered, using insertion sort isn't so bad. In most cases
 * we are O(n) time with O(1) extra memory. O(n^2) worst case time though.
 * 
 * If this function is becomming a drag, consider quicksort or mergesort.
 * */
void sort_triplet(cholmod_triplet* triplet)
{
    // get pointers to triplet storage elements
    int* ti = (int*)triplet->i;  // row index
    int* tj = (int*)triplet->j;  // col index
    double* tx = (double*)triplet->x;  // value

    // start with the first item 'sorted' then find a place to put the item at
    // position next.
    for (size_t next = 1; next < triplet->nnz; next++) {
        for (size_t i = next - 1; i >= 0; i--) {
            if (_triplet_comparitor(triplet, i, i+1) > 0) {

                // swap the ith and (i+1)th items, but we have to do this for
                // all three lists (the row/col index and conductance value)
                int tempi = ti[i];
                int tempj = tj[i];
                double tempx = tx[i];

                ti[i] = ti[i+1];
                tj[i] = tj[i+1];
                tx[i] = tx[i+1];

                ti[i+1] = tempi;
                tj[i+1] = tempj;
                tx[i+1] = tempx;
            }
            else break;  // otherwise, these are sorted so move to next one
        }
    }
}

// warning: don't print large matrices
void print_adj(const adj_list_t &adj_c, const char *name)
{
    std::cout << name << std::endl;
    int i = 0;
    for (auto a : adj_c) {
         std::cout << i++ << ": ";
         for (auto b : a) {
             int c; double d;
             std::tie(c, d) = b;
             std::cout << "(" << c << ", " << d << ") ";
         }
         std::cout << "\n";
     }
     std::cout << "\n";
}

#define TRIPLET_DEBUG_LIMIT 10
void print_triplet(solver_state_t* state, cholmod_triplet* triplet, const char* name)
{
    printf("%s:\n", name);
    cholmod_print_triplet(triplet, NULL, state->common);

    // If TRIPLET_DEBUG_LIMIT == -1, print all triplet entries. Otherwise, just
    // print up to TRIPLET_DEBUG_LIMIT number of entries
    size_t n_print;
    if (TRIPLET_DEBUG_LIMIT > 0) n_print = min(triplet->nnz, TRIPLET_DEBUG_LIMIT);
    else n_print = triplet->nnz;

    // print each entry
    for (size_t i =0; i < n_print; i++) {
        printf("%d %d %f\n", ((int*)triplet->i)[i], ((int*)triplet->j)[i], ((double*)triplet->x)[i]);
    }

    // Show ... if we skipped some entires
    if (n_print < triplet->nnz) printf("...\n\n");
    else printf("\n");
}

// warning: don't print large matrices
void print_sparse(solver_state_t *state, cholmod_sparse *A, const char* name)
{
    printf("%s:\n", name);
    // print info about matrix
    cholmod_print_sparse(A, NULL, state->common);
    size_t nnz = cholmod_nnz(A, state->common);
    // print values
    for (size_t i = 0; i < A->ncol+1; i++) {
        printf("%d ", ((int*)A->p)[i]);
    }
    printf("\n");
    for (size_t i = 0; i < nnz; i++) {
        printf("%d ", ((int*)A->i)[i]);
    }
    printf("\n");
    if (A->xtype == CHOLMOD_REAL) {
        for (size_t i = 0; i < nnz; i++) {
            printf("%f ", ((double*)A->x)[i]);
        }
        printf("\n");
    }
    printf("\n");
}

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
    cholmod_triplet *triplet,  // Note: will sort the triplets
    const int *input_groups,  // array of indices which are input groups
    const size_t n_input_groups,  // number of input groups
    const int *output_groups,  // array of indices which are output groups
    const size_t n_output_groups   // number of output groups
)
{
    if (n_output_groups < 1) {
        fprintf(stderr, "If you don't specify any output voltages, then I have "
                "nothing to calculate. Try adding some output groups.\n");
        return -1;
    }

    size_t n_groups = triplet->nrow;  // ie, the total number of groups

    if (n_input_groups == n_groups) {
        fprintf(stderr, "Nothing to do. All the voltages are known.\n");
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

    // remember these
    state->n_groups = n_groups;
    state->n_input_groups = n_input_groups;
    state->n_output_groups = n_output_groups;

    // make sure triplet is sorted by column then row in increasing order
    // the adjacancy list will also be sorted then
    sort_triplet(triplet);

    // make set of input and output groups for faster O(1) lookup
    if (!state->is_input)
        state->is_input = new bool[n_groups]();
    for (size_t i = 0; i < n_input_groups; i++) {
        // assume that the group indices are not out of bounds 
        state->is_input[input_groups[i]] = true;
    }
    if (!state->is_output)
        state->is_output = new bool[n_groups]();
    for (size_t i =0; i < n_output_groups; i++) {
        state->is_output[output_groups[i]] = true;
    }

    // column first adjacancy list
    size_t nnz = triplet->nnz;
    // column wise adjacancy list (+ diagonals)
    // keeps all column entries down from (and including) the diagonals
    adj_list_t adj_c(n_groups);
    // row wise adjacancy list
    // keps all the row entries to the left of the diagonal
    adj_list_t adj_r(n_groups);
    // we keep two lists so that the values remain in order and sorted
    for (size_t ig = 0, it = 0; ig < n_groups; ig++) {
        // placeholder for diagonals only in column wise one
        adj_c[ig].push_back(std::make_pair(ig, 0.0));
        while (it < nnz && (size_t)(((int*)triplet->j)[it]) == ig) {
            // append row index and conductance value to column vector
            size_t irow = ((int*)triplet->i)[it];
            double cond = ((double*)triplet->x)[it];
            // -ve off diagonals by convension
            // do both directions, stored in seporate adj lists
            adj_c[ig].push_back(std::make_pair(irow, -cond));
            adj_r[irow].push_back(std::make_pair(ig, -cond));
            it++;
        }
    }

    // we don't need the triplet anymore
    cholmod_free_triplet(&triplet, state->common);

    // sum conductances for diagonal entries
    // use -= since we are summing -ve entries but want a +ve diagonal
    for (size_t i = 0; i < n_groups; i++) {
        // column wise
        for (auto &pair : adj_c[i]) {
            double &cond = pair.second;
            adj_c[i][0].second -= cond;
        }
        // row wise
        for (auto &pair : adj_r[i]) {
            double &cond = pair.second;

            // add only to column wise entries since the row wise ones don't
            // store diagonals
            adj_c[i][0].second -= cond;
        }
    }

    // print_adj(adj_c, "Column wise adjacancy:");
    // print_adj(adj_r, "Row wise adjacancy:");

    // A matrix size
    size_t n_mat = n_groups - n_input_groups;
    state->n_mat = n_mat;

    // map from group index to matrix index
    state->group_to_mat_map = new size_t[n_groups];
    state->mat_to_group_map = new size_t[n_mat];

    size_t ig = 0;  // current group index
    size_t i = 0;  // current matrix index
    for (; ig < n_groups; ig++) {
        if (state->is_input[ig]) {
            // mark as input group
            state->group_to_mat_map[ig] = -1;
        } else {
            // create group mapping from group index to matrix index and visa versa
            state->group_to_mat_map[ig] = i;
            state->mat_to_group_map[i] = ig;
            i++;  // next matrix index
        }
    }

    // the empty matrix we are going to fill
    cholmod_sparse* A = cholmod_allocate_sparse(
        n_mat, n_mat, nnz + n_groups, true, true,
        LOWER_TRIANGULAR, CHOLMOD_REAL, state->common);

    state->A = A;  // save a reference in state

    // the empty matrix we are going to fill
    cholmod_sparse* G = cholmod_allocate_sparse(
        n_mat, n_input_groups, n_input_groups * n_mat, true, true,
        LOWER_TRIANGULAR, CHOLMOD_REAL, state->common);

    state->G = G;  // save a reference in state

    // fill in A matrix as reduced symmetric column form
    size_t ix = 0;  // current A matrix value index
    size_t jx = 0;  // current G matrix value index 
    ig = 0;  // current group
    i = 0;  // current matrix index A (unknown)
    size_t j = 0;  // current matrix index G (input)

    // get pointers to sparse storage elements of A and G for compressed col form
    int* Gp = (int*)(G->p);
    int* Gi = (int*)(G->i);
    double* Gx = (double*)(G->x);

    int* Ap = (int*)(A->p);
    int* Ai = (int*)(A->i);
    double* Ax = (double*)(A->x);

    // now let's actually construct the A and G matrices
    for (ig = 0; ig < n_groups; ig++) {
        if (state->is_input[ig]) {
            Gp[j] = jx;

            // Add all non-input column entries above the diagonal.
            // Since the matrix is symmetric, this is the same as
            // adding all row entries to the left, which we known.
            for (auto &pair : adj_r[ig]) {
                size_t &jr = pair.first;  // get the row index
                double &g = pair.second;  // conductance

                if (!state->is_input[jr]) {
                    // save matrix value for half of the matrix
                    Gx[jx] = -g;  // -ve since on rhs
                    Gi[jx] = state->group_to_mat_map[jr];  // row index in column
                    jx++;
                }
            }
            // add all non-input col entries below (and including) the diagonal
            for (auto &pair : adj_c[ig]) {
                size_t &jr = pair.first;  // get the row index
                double &g = pair.second;  // conductance

                if (!state->is_input[jr]) {
                    // save matrix value for half of the matrix
                    Gx[jx] = -g;  // -ve since on rhs
                    Gi[jx] = state->group_to_mat_map[jr];  // row index in column
                    jx++;
                }
            }
            j++;
        } else {
            Ap[i] = ix;  // pointer to column start

            for (auto &pair : adj_c[ig]) {
                size_t &ir = pair.first;  // get the row index
                double &g = pair.second;  // conductance

                if (!state->is_input[ir]) {
                    // save matrix value for half of the matrix
                    Ax[ix] = g;
                    Ai[ix] = state->group_to_mat_map[ir];  // row index in column
                    ix++;
                }
            }
            // move to next matrix index
            i++;
        }
    }
    // the number of nonzeros
    // the array A->p is of length [n columns] + 1 where the last is nnz
    size_t nz = ix;
    state->nnz = nz;
    Ap[i] = nz;  // for A matrix
    Gp[j] = jx;  // for G matrix

    // now we know the actual nnz, so shead unused memory
    cholmod_reallocate_sparse(ix, A, state->common);
    cholmod_reallocate_sparse(jx, G, state->common);

    // create b_set pattern
    // This is the pattern of output voltages we care about, which means we can
    // reduce the number of multiplications. It doesn't store values, only indices
    // in the b vector
    cholmod_sparse *b_set = cholmod_allocate_sparse(n_mat, 1, n_mat, true, true, 0,
                                                    CHOLMOD_PATTERN, state->common);

    state->b_set = b_set;
    int *b_set_p = (int*)(b_set->p);
    int *b_set_i = (int*)(b_set->i);

    // std::sort(output_groups, output_groups + n_output_groups);
    b_set_p[0] = 0;  // point to first column
    ig = 0;
    i = 0;
    for(; ig < n_output_groups; ig++) {
        if (!state->is_input[output_groups[ig]]) {
            b_set_i[i] = state->group_to_mat_map[output_groups[ig]];
            i++;
        }
    }
    b_set_p[1] = i;  // nnz

    print_sparse(state, state->b_set, "b_set");

    return 0;
}

/**
 * Create the inital decomposition and configure cholmod ready for some
 * speedy calculations. Must call solver_initalise_network on a state first.
 * 
 * @param state: the state object of the system we are solving
 * */
int solver_begin(solver_state_t *state)
{
    if (!state->A || !state->G) {
        std::cout << "Must call solver_initalise_network() before starting.\n";
        return -1;
    }

    // get symbolic factor and compute permutation
    state->LD = cholmod_analyze(state->A, state->common);
    // get actual factorisation for LDLT
    cholmod_factorize(state->A, state->LD, state->common);

    // vectors allocated for the injected current vector b and solution vector x
    // each vector is stored in dense form with a sparse pattern which tells
    // cholmod which rows to compute (hense reducing the number of computations)
    const size_t &n_mat = state->n_mat;
    state->x = cholmod_allocate_dense(n_mat, 1, n_mat, CHOLMOD_REAL, state->common);
    state->x_set = cholmod_allocate_sparse(n_mat, 1, n_mat, true, true, 0,
                                           CHOLMOD_PATTERN, state->common);

    state->b = cholmod_allocate_dense(n_mat, 1, n_mat, CHOLMOD_REAL, state->common);

    // vectors allocated for intermediate steps and used by cholmod solve2
    // eg the y = L\b step etc
    state->y = cholmod_allocate_dense(n_mat, 1, n_mat, CHOLMOD_REAL, state->common);
    state->e = cholmod_allocate_dense(n_mat, 1, n_mat, CHOLMOD_REAL, state->common);

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
    print_triplet(state, triplet, "Input");

    // Set which indices we want to be inputs and outputs. This allows us to
    // deterwmine which voltages are known and unknown and which solutions we
    // care about
    int input_groups[] = {1, 2}; int n_inputs = 2;  // set v_1 and v_2 as known
    int output_groups[] = {0, 1, 2, 3}; int n_ouputs = 4;  // get all voltages

    // Initalise the solver. This will construct the A matrix and b vector used
    // to represent the system.
    solver_initalise_network(state, triplet,
                             &input_groups[0], n_inputs,
                             &output_groups[0], n_ouputs);

    print_sparse(state, state->A, "A");
    print_sparse(state, state->G, "G");

    solver_destroy_state(state);
    // state can no longer be used
}

void testcase_1(size_t iter)
{
    // initalise module
    solver_state_t* state = solver_create_state();

    cholmod_triplet *triplet = triplet_from_file(state, "matrices/testcase_0.mtx");
    print_triplet(state, triplet, "Input");

    // Set which indices we want to be inputs and outputs. This allows us to
    // deterwmine which voltages are known and unknown and which solutions we
    // care about
    int input_groups[] = {1, 2}; int n_inputs = 2;  // set v_1 and v_2 as known
    int output_groups[] = {0, 1, 2, 3}; int n_ouputs = 4;  // get all voltages

    // Initalise the solver. This will construct the A matrix and b vector used
    // to represent the system.
    solver_initalise_network(state, triplet,
                             &input_groups[0], n_inputs,
                             &output_groups[0], n_ouputs);

    print_sparse(state, state->A, "A");
    print_sparse(state, state->G, "G");

    solver_begin(state);

    double known_voltages[] = {1, 0};
    cholmod_dense *v = cholmod_zeros(2, 1, CHOLMOD_REAL, state->common);
    double *vx = (double*)(v->x);
    for (size_t i = 0; i < 2; i++) vx[i] = known_voltages[i];

    for (size_t i = 0; i < iter; i++) {
        // compute b matrix
        // b = alpha*(A*x) + beta*b
        double alpha[] = {1, 1};
        double beta[] = {0, 0};
        cholmod_sdmult(state->G, false, alpha, beta, v, state->b, state->common);
        cholmod_print_dense(state->b, NULL, state->common);
        double *bx = (double*)(state->b->x);
        for (size_t j = 0; j < state->n_mat; j++) {
            std::cout << bx[j] << ' ';
        }
        std::cout << std::endl;

        cholmod_solve2(CHOLMOD_A, state->LD, state->b, state->b_set, &state->x,
                       &state->x_set, &state->y, &state->e, state->common);
        
        cholmod_print_dense(state->x, NULL, state->common);
        double *xx = (double*)(state->x->x);
        for (size_t j = 0; j < state->n_mat; j++) {
            std::cout << xx[j] << ' ';
        } 
        std::cout << std::endl;
    }

    solver_destroy_state(state);
}

int main ()
{
    testcase_1(2);
}
