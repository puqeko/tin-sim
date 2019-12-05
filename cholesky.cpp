/* make_symmetric.cpp
 * 
 * Checks to see if a matrix, given in .mtx format to stdin, is stored as a
 * symmetric matrix in lower triangular form. If it is not, check if we can
 * safely convert to lower triangular symmetric form, do the conversion, then
 * save out to stdout.
 * 
 * Info is printed to stderr.
 * 
 * Althought cpp is used, the interface is kept to c.
 * 
 * example use: convert the simple.mtx file to symmetric format
 * ./make_symmetric.o < matrices/simple.mtx > out.mtx
 */

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
    size_t n_groups;  // use int over size_t for consistancy with cholmod
    size_t nnz;  // number of non-zero values

    // TODO: test if we can use x to store input voltages + solutions togther
    // or is it not safe?
    // double *voltages;

    // stuff needed to solve Ax = b, were x [volts] and b [amps] are vectors
    cholmod_sparse *A;  // the full conductance matrix
    cholmod_factor *LD;  // the LD decomposition
    cholmod_sparse *G;  // matrix of conductances on the right-hand-side
    cholmod_dense *b;  // vector of injected currents b = G v [known voltages]
    cholmod_sparse *b_set;  // pattern of groups which have injected currents
    cholmod_dense *x;  // vector of solved voltages
    cholmod_sparse *x_set;  // pattern of groups voltages with defined solutions

    // vectors of length N, if A is NxN, which are allocated for intermediate
    // steps in the solve process and reused by the cholmod api
    cholmod_dense *y_handle;
    cholmod_dense *e_handle;
} solver_state_t;

// a vector or vectors each containing pairs (the row index and conductance)
typedef std::pair<size_t, double> link_t;
typedef std::vector<std::vector<link_t> > adj_list_t;

/**
 * Allocate memory for a state object needed for cholmod operations, configure,
 * and allow chomod to initalise itself (clears pointers, stats, default params).
 * Must be called before anything else to obtain our state object.
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
    if (state->A) cholmod_free_sparse(&state->A, state->common);
    if (state->LD) cholmod_free_factor(&state->LD, state->common);
    if (state->G) cholmod_free_sparse(&state->G, state->common);
    if (state->b) cholmod_free_dense(&state->b, state->common);
    if (state->b_set) cholmod_free_sparse(&state->b_set, state->common);
    if (state->x) cholmod_free_dense(&state->x, state->common);
    if (state->x_set) cholmod_free_sparse(&state->x_set, state->common);
    if (state->y_handle) cholmod_free_dense(&state->y_handle, state->common);
    if (state->e_handle) cholmod_free_dense(&state->e_handle, state->common);

    // finish up cholmod and free state
    cholmod_finish(state->common);
    delete state->common;
    delete state;
}


// Used in insertion sort for triplets to return column followed by row ordering
// Returns -1 if left < right, 1 if left > right, else 0
int _triplet_comparitor(cholmod_triplet *triplet, const int &l, const int &r)
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
    // start with the first item 'sorted' then find a place to put the item at
    // position next.
    for (size_t next = 1; next < triplet->nnz; next++) {
        for (size_t i = next - 1; i >= 0; i--) {
            if (_triplet_comparitor(triplet, i, i+1) > 0) {

                // swap the ith and (i+1)th items, but we have to do this for
                // all three lists (the row/col index and conductance value)
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

// print an adjacancy list
void print_adj(adj_list_t &adj_c)
{
    std::cout << "Adjacancy list:\n";
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

void print_sparse(solver_state_t *state, cholmod_sparse *A)
{
    printf("Sparse:\n");
    for (size_t i = 0; i < state->n_groups; i++) {
        printf("%d ", ((int*)A->p)[i]);
    }
    printf("\n");
    for (size_t i = 0; i < state->nnz; i++) {
        printf("%d ", ((int*)A->i)[i]);
    }
    printf("\n");
    for (size_t i = 0; i < state->nnz; i++) {
        printf("%f ", ((double*)A->x)[i]);
    }
    printf("\n");
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
 * @param n_input_groups: number of input groups
 * @param output_groups: a vector of groups we want to know the voltage of
 * @param n_output_groups: number of output groups
 * 
 * @returns: 0 on success, error otherwise
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

    // make sure triplet is sorted by column then row in increasing order
    // the adjacancy list will also be sorted then
    sort_triplet(triplet);

    size_t n_groups = triplet->nrow;  // ie, the total number of groups
    state->n_groups = n_groups;

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

    // the empty matrix we are going to fill
    cholmod_sparse* A = cholmod_allocate_sparse(
        n_groups, n_groups, nnz + n_groups, true, true,
        LOWER_TRIANGULAR, CHOLMOD_REAL, state->common);

    state->A = A;  // save a reference in state

    print_adj(adj_c);
    print_adj(adj_r);

    // fill in A matrix as reduced symmetric column form
    size_t ix = 0;  // current A matrix value index
    size_t ig = 0;  // current group
    size_t i = 0; // current matrix index

    for (ig = 0; ig < n_groups; ig++) {
        // use inputs for sum but don't add to sparse matrix A
        // otherwise, this group is an unknown and should be added to A
        ((int*)A->p)[ig] = ix;  // pointer to column start

        for (auto &pair : adj_c[ig]) {
            size_t &ir = pair.first;  // get the row index
            double &g = pair.second;  // conductance

            // save matrix value for half of the matrix
            ((double*)A->x)[ix] = g;
            ((int*)A->i)[ix] = ir;  // row index in column
            ix++;
        }
    }
    // the number of nonzeros
    // the array A->p is of length [n columns] + 1 where the last is nnz
    size_t nz = ix;
    state->nnz = nz;
    ((int*)A->p)[ig] = nz;

    cholmod_print_sparse(A, NULL, state->common);
    print_sparse(state, A);

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


#define TRIPLET_DEBUG_LIMIT 10
void print_triplet(solver_state_t* state, cholmod_triplet* triplet)
{
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

    solver_destroy_state(state);
    // state can no longer be used
}


int main ()
{
    testcase_0();
}
