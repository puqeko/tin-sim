/*
 * 
 * 
 */

#include <iostream>
#include <cholmod.h>

int main ()
{
    cholmod_sparse *A;
    cholmod_common c;

    cholmod_start(&c);
    A = cholmod_read_sparse(stdin, &c);

    if (A->stype == 0) std::cout << "stype=0\n";
    else std::cout << ":(\n";

    // B = cholmod_copy(A, )

    // use simple.mtx
    // convert to symettric and check for symmetry
    // sort A
    // use cholmod_symmetry

    cholmod_free_sparse(&A, &c);
    cholmod_finish(&c);
}
