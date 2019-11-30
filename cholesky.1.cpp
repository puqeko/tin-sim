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
#include <fstream>
#include <string>
#include <vector>
#include <unordered_set>
#include <cholmod.h>
#include <stdint.h>

#define LOWER_TRIANGULAR -1
#define UPPER_TRIANGULAR 1

/**
 * Assume all groups are connected.
 * */
int read_mtx_file(std::string filename, size_t n_nodes)
{
    // allocate our adjacency list
    // [
    //     group n_0: [(n_1, g_1), (n_2, g_2), ...]
    // ]
    // such that group n_0 connects to n_1 with conductance g_1 etc
    adj_list adj(n_nodes);
    std::ifstream file;
    file.open(filename);

    if (!file.is_open()) {
        std::cerr << "Could not open '" << filename << "'.\n";
        return -1; 
    }

    // look for header of form '%%MatrixMarket
    std::string line; 
    while (std::getline(file, line))
    {

    }

    file.close();
}

static cholmod_common common;

/**
 * Create common object needed for cholmod settings which are used by pretty
 * much all cholmod functions. 'common' is a global of this module and must
 * be initalised before any other call to cholmod.
 * */
void initalise_cholmod()
{
    cholmod_start(&common);

    // use lower triangular in the default case when storing sparse matrices.
    common.prefer_upper = false;
}

/**
 * Construct and preorder a sparse conductance matrix.
 * Uses CHOLMOD's preording techniques based on graph structure and the
 * input groups.
 * */
cholmod_triplet* triplets_from_file
(
    // the filename to read from, as a string
    const std::string &filename
)
{
    // read file into triplet form
    FILE *file = fopen(filename.c_str, "r");
    cholmod_triplet* conductance_coords = cholmod_read_triplet(file, &common);

    if (conductance_coords->ncol != conductance_coords->nrow) {
        std::cerr << "Must be a sqaure matrix.\n";
        cholmod_free_triplet(&conductance_coords, &common);
        return NULL;
    }

    return conductance_coords;
}

/**
 * Get the A matrix and b vector from a triplet.
 * 
 * The input triplet is a list where each entry is a coordinate between two
 * groups and the conductance of that connection. The input triplet is expected
 * to be symmetric (ie only specify entries where row > col) and specify NO
 * self conductance (ie no conduntance between the same group a->a). The network
 * must also be fully connected - ie have only a single component which consists
 * of more than one group.
 * 
 * The A matrix is a conductance matrix, stored as a sparse matrix in
 * compressed column form. All rows/columns corrisponding to input voltages
 * are moved to the right-hand-side (the b matrix), since the corrisponding
 * voltages are known at each timestep, and removed from A.
 * 
 * The b vector is a vector of injected currents such that the input voltages
 * are maintained at the specified groups. If no input voltages are specified,
 * the system has the trivial solution (b = 0). If only one input voltage is
 * specified, then the system has the solution b = a, where a is the specified
 * input voltage at some group. Hence, input voltages should be at-least 2.
 * 
 * @param triplet: a cholmod_triplet describing the network (this will get freed)
 * @param input_groups: a vector of groups which have a known voltage
 * 
 * outputs:
 * @param A_out: the sparse A matrix
 * @param b_out: the sparse b matrix
 * @returns: 0 on success, error otherwise
 * */
int system_from_triplet
(
    // a list of connections between groups and their conductances
    // this triplet will be freed
    const cholmod_triplet *triplet,
    // groups that we know the voltage of
    const std::vector<int> &input_groups,

    // outputs
    // the sparse conductance matrix, or the A matrix
    cholmod_sparse* A_out,
    // the injected currents, or the b vector
    cholmod_sparse* b_out
)
{
    // see main comment above
    if (input_groups.size() < 2) {
        std::cerr << "Warning: You should specify more than two input voltages."
                     "The system isn't going to be very interesting otherwise.\n";
    }

    // make set of input groups for fast O(1) lookup
    std::unordered_set<int> input_node_lookup (input_groups.begin(), input_groups.end());

    // build new triplet to construct a conductance matrix A
    // assuming lower triangular form (row > col) with no diagonal entries, then
    // we must add the self conductance terms to the existing matrix (the diagonls)
    // and remove entries corrisponding to know input votlages.

    // allocate space to store our now coordinates
    // use the nnz + extra diagonal entries as a starting estimate for the memory
    // we need. It will be less than this since we remove some entries, but
    // we can get std vector to do the resizing.
    int len_est = triplet->nnz + triplet->nrow;
    int *i_row = new int[len_est]; // row coordinates
    int *i_col = new int[len_est]; // column coordinates
    double *g = new double[len_est];  // conductance value

    // total self conductance (diagonal entries)
    double *g_total = new double[triplet->nrow];

    double *b = new double[triplet->nnz];  // b vector

    // TODO: make some assertions about symmetry, type, etc

    int i = 0;
    for (; i < triplet->nnz; i++) {
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
            b[ci] += gi;  // move -gi (upper triangular term) to other side
        } else if (input_node_lookup.count(ci)) {
            b[ri] += gi;  // move -gi (lower triangular term) to other side
        } else {
            // this isn't a known input row/col
            i_row[i] = ri;
            i_col[i] = ci;
            g[i] = -gi;  // -ive off diagonals
        }
    }

    // add diagonal self conductance entries
    int j = 0;
    for (; j < triplet->nrow; j++, i++) {
        if (!input_node_lookup.count(j)) {
            // ignoring input nodes
            i_row[i] = j;
            i_col[i] = j;
            g[i] = g_total[j];
        }
    }

    cholmod

    cholmod_sparse *A = cholmod_triplet_to_sparse(new_triplet, new_triplet->nnz, &common);
    

    return 0;
}

// https://stackoverflow.com/questions/2744181/how-to-call-c-function-from-c

int main ()
{
    cholmod_sparse *Ain, *A;

    cholmod_start(&c);
    Ain = cholmod_read_sparse(stdin, &c);

    // is this matrix already valid?
    if (Ain->stype != LOWER_TRIANGULAR) {
        std::cerr << "This matrix is not in lower triangular form. Use make_symmetric.\n";
        cholmod_free_sparse(&Ain, &c);
        cholmod_finish(&c);
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
