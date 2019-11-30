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

#include <iostream>

#include <cholmod.h>
#include <stdint.h>

#include <string>
#include <vector>
#include <unordered_set>

#define LOWER_TRIANGULAR -1
#define UPPER_TRIANGULAR 1


static cholmod_common g_common;

/**
 * Create g_common object needed for cholmod settings which are used by pretty
 * much all cholmod functions. 'g_common' is a global of this module and must
 * be initalised before any other call to cholmod.
 * */
void initalise_cholmod()
{
    cholmod_start(&g_common);

    // use lower triangular in the default case when storing sparse matrices.
    g_common.prefer_upper = false;
}


cholmod_triplet* triplet_from_file
(
    // the filename to read from, as a string
    const std::string &filename
)
{
    // read file into triplet form
    FILE *file = fopen(filename.c_str, "r");
    cholmod_triplet* conductance_coords = cholmod_read_triplet(file, &g_common);

    if (conductance_coords->ncol != conductance_coords->nrow) {
        std::cerr << "Must be a sqaure matrix.\n";
        cholmod_free_triplet(&conductance_coords, &g_common);
        return NULL;
    }

    return conductance_coords;
}

// when we construct the matrix, we remove some rows which means that the
// ith group may not corrispond to the ith row in the matrix anymore. Hence,
// j = g_input_groups[i] tells us that group i is row j in the matrix.
static std::vector<int> g_group_to_mat_map;
static std::vector<int> g_mat_to_group_map;

// outputs
// the sparse conductance matrix, or the A matrix
static cholmod_sparse* g_A;
// the injected currents, or the b vector
static cholmod_sparse* g_b;

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
 * entries (excluding globals)
 * 
 * @param triplet: a cholmod_triplet describing the network (this will get freed!)
 * @param input_groups: a vector of groups which have a known voltage
 * 
 * @returns: 0 on success, error otherwise
 * */
int build_system_from_network
(
    // a list of connections between groups and their conductances, there are
    // helper function to help construct this object. This triplet will be freed
    cholmod_triplet *triplet,
    // the indices of groups that we will know the voltage of (ie group number)
    const std::vector<int> &input_groups
)
{
    if (input_groups.size() == 0)
    {
        std::cerr << "Must call set_input_groups(...) first.\n";
        return -1;
    }

    // see main comment above
    if (input_groups.size() < 2) {
        std::cerr << "Warning: You should specify more than two input voltages. "
                     "The system isn't going to be very interesting otherwise.\n";
        // still do it though
    }

    // make some assertions about symmetry, type, etc of triplet
    if (triplet->ncol != triplet->nrow || triplet->xtype != CHOLMOD_DOUBLE ||
        triplet->stype != LOWER_TRIANGULAR) {
        std::cerr << "Triplet must be symmetric, lower triangular, and of type"
                     " float.\n";
        return -1;
    }

    // make set of input groups for fast O(1) lookup
    std::unordered_set<int> input_node_lookup (input_groups.begin(), input_groups.end());

    // map from group index to matrix index
    g_group_to_mat_map = std::vector<int>(triplet->nrow);
    g_mat_to_group_map = std::vector<int>(triplet->nrow);

    // build new triplet to construct a conductance matrix A
    // assuming lower triangular form (row > col) with no diagonal entries, then
    // we must add the self conductance terms to the existing matrix (the diagonls)
    // and remove entries corrisponding to know input votlages.

    // allocate space to store our now coordinates
    // use the nnz + extra diagonal entries as a starting estimate for the memory
    // we need. It will be less than this since we remove some entries.
    int len_est = triplet->nnz + triplet->nrow;

    cholmod_triplet* new_triplet;
    new_triplet = cholmod_allocate_triplet(triplet->nrow, triplet->nrow, len_est,
                                           LOWER_TRIANGULAR, CHOLMOD_DOUBLE, &g_common);

    // total self conductance (diagonal entries)
    double *g_total = new double[triplet->nrow];
    double *b_temp = new double[triplet->nrow];  // b vector
    // By kirchoffs current law, the sum of currents in/out of a group is zero
    // assuming no injected current.
    for (int i = 0; i < triplet->nnz; i++) b_temp[i] = 0;

    int t = 0;  // the current index in the new_triplet
    for (int i = 0; i < triplet->nnz; i++) {
        // get row, col coords and conductance value of current triple
        const int &ri = ((int*)triplet->i)[i];
        const int &ci = ((int*)triplet->j)[i];
        const double &gi = ((double*)triplet->x)[i];

        // get total conductance into/out-of a group for the diagonals of
        // the conductance matrix we are building.
        g_total[ri] += gi;
        g_total[ci] += gi;

        // if the row or column belongs to our known input groups then
        // adjust b matrix, ie our injected currents, accordingly.
        // Remember, the only the lower triangular tripples are given.
        if (input_node_lookup.count(ri)) {
            b_temp[ci] += gi;  // move -gi (upper triangular term) to other side
        } else if (input_node_lookup.count(ci)) {
            b_temp[ri] += gi;  // move -gi (lower triangular term) to other side
        } else {
            // this isn't a known input row/col
            ((int*)new_triplet->i)[t] = ri;
            ((int*)new_triplet->j)[t] = ci;
            ((double*)new_triplet->x)[t] = -gi;  // -ive off diagonals
            t++;  // this will lag behind i
        }
    }

    // add diagonal self conductance entries
    for (int j = 0; j < new_triplet->nrow; j++) {
        if (!input_node_lookup.count(j)) {
            // ignoring input nodes
            ((int*)new_triplet->i)[t] = j;
            ((int*)new_triplet->j)[t] = j;
            ((double*)new_triplet->x)[t] = g_total[j];

            // map between matrix coords and actual group indices
            g_group_to_mat_map[j] = t;
            g_mat_to_group_map[t] = j;
            t++;
        }
    }
    cholmod_free_triplet(&triplet, &g_common);
    delete g_total;
    // resize now that we know exactly how big new_triplet is
    cholmod_reallocate_triplet(t, new_triplet, &g_common);

    // convert A to sparse
    cholmod_sparse *A;
    A = cholmod_triplet_to_sparse(new_triplet, new_triplet->nnz, &g_common);
    cholmod_free_triplet(&new_triplet, &g_common);
    g_A = A;  // save pointer to global

    // convert b to sparse
    cholmod_triplet *b_triplet;
    b_triplet = cholmod_allocate_triplet(new_triplet->nrow, new_triplet->nrow,
                                         new_triplet->nrow, 0, CHOLMOD_DOUBLE, &g_common);

    // populate with non-zero b entries
    int u = 0;  // b_triplet index
    for (int k = 0; k < b_triplet->nrow; k++) {
        if (b_temp[k] != 0.0) {
            ((int*)b_triplet->i)[u] = k;
            ((int*)b_triplet->j)[u] = 1;  // only one column
            ((double*)b_triplet->x)[u] = b_temp[k];
            u++;
        }
    }
    cholmod_reallocate_triplet(u, b_triplet, &g_common);  // shrink to right size
    delete b_temp;

    // convert b to sparse
    // we do this not because it saves memory, but because it reduces the
    // number of computations required every time we multiply an element of b.
    cholmod_sparse *b;
    b = cholmod_triplet_to_sparse(b_triplet, u, &g_common);
    cholmod_free_triplet(&b_triplet, &g_common);
    g_b = b;  // save to global

    return 0;
}

// https://stackoverflow.com/questions/2744181/how-to-call-c-function-from-c

int main ()
{
    cholmod_sparse *Ain, *A;

    initalise_cholmod();
    Ain = cholmod_read_sparse(stdin, &g_common);

    // is this matrix already valid?
    if (Ain->stype != LOWER_TRIANGULAR) {
        std::cerr << "This matrix is not in lower triangular form. Use make_symmetric.\n";
        cholmod_free_sparse(&Ain, &g_common);
        cholmod_finish(&g_common);
        return 0;
    }


    // cholmod_sort(Ain, &c);  // needed for cholmod_symmetry
    // int symm = cholmod_symmetry(Ain, 0, NULL, NULL, NULL, NULL, &c);

    // // requires positive diagonals
    // // requires symmetry
    // // if symm == CHOLMOD_MM_SYMMETRIC_POSDIAG then matrix suitible for LDLT
    // if (symm != CHOLMOD_MM_SYMMETRIC_POSDIAG) {
    //     std::cerr << "Oh no, this matrix is not positive-definate-symmetric.\n";
    //     cholmod_free_sparse(&Ain, &c);
    //     cholmod_finish(&c);
    //     return 0;
    // }

    // std::cerr << "Great, this matrix is positive-definate-symmetric.\n";
    // std::cerr << "Storing as triangular matrix.\n";

    // // create new matrix A that stores lower triangular values only
    // // value 1 means use values mode
    // A = cholmod_copy(Ain, LOWER_TRIANGULAR, 1, &c);

    // cholmod_write_sparse(stdout, A, NULL, NULL, &c);
    // std::cerr << "Sent to stdout.\n";

    // cholmod_free_sparse(&Ain, &c);
    // cholmod_free_sparse(&A, &c);
    // cholmod_finish(&c);
}
