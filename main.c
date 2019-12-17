#include "solve.h"

static double *known_voltages;

int network_update(solver_state_t* state, solver_triplet_t* triplet,
                   double** voltages, int iter, void* data)
{
    // change voltages pointer to point to our known voltages (fixed voltages)
    // alternativly, we could write new voltages each iter to (*voltages)[i] = vi
    // these voltages are in index order, corrisponding to the input_groups passed
    // in at startup. Ie, the ith value in *voltages is the ith group in input_groups.
    *voltages = known_voltages;
    printf("Iter: %d\n", iter);

    // Build a triplet, which will update the conductance of some connections.
    // These conductances must be deltas.
    // We are only going to change one link in this example.

    // Get pointers to update triplet
    int *Trow = (int*)(triplet->i);
    int *Tcol = (int*)(triplet->j);
    double *Tcond = (double*)(triplet->x);

    // When specifing indices, it is important to make sure each connection is
    // added only once. ie don't add 3->0 and 0->3. We always specify an update
    // such that row > col (lower triangular).
    Trow[0] = 3;
    Tcol[0] = 0;
    // alternate between adding and removing 0.1 to conductance of connection from 0 to 3
    Tcond[0] = (iter%2==0) ? +0.1 : -0.1;
    triplet->nnz = 1;  // specify that we have 1 update in the triplet

    // add another update as such
    // Trow[1] = 1;
    // Tcol[1] = 0;
    // Tcond[1] = -0.3;
    // triplet->nnz = 2;

    // data is a pointer to whatever is passed in to solver_iterate_ac if you
    // want to remember stuff between calls
    return 0;
}

int main (void)
{
    solver_state_t* state = solver_create_state();
    solver_triplet_t* triplet = triplet_from_file(state, "matrices/testcase_0.mtx");
    print_triplet(state, triplet, "Input");

    // allocate know voltages on the heap
    known_voltages = calloc(2, sizeof(double));
    known_voltages[0] = 1;
    known_voltages[1] = 0;

    // Set which indices we want to be inputs and outputs. This allows us to
    // deterwmine which voltages are known and unknown and which solutions we
    // care about
    int input_groups[] = {1, 2}; int n_inputs = 2;  // set v_1 and v_2 as known
    int output_groups[] = {0, 1, 2, 3}; int n_ouputs = 4;  // get all voltages

    const int n_iters = 1;  // change this to the number of iterations
    
    solver_iterate_ac(
        state, n_iters, triplet, input_groups, n_inputs, output_groups, n_ouputs,
        network_update, NULL
    );
    solver_destroy_state(state);
    return 0;
}
