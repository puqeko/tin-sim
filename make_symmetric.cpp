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

#define LOWER_TRIANGULAR -1
#define UPPER_TRIANGULAR 1

int main ()
{
    cholmod_sparse *Ain, *A;
    cholmod_common c;
    c.prefer_upper = false;  // use lower triangular in the default case

    cholmod_start(&c);
    Ain = cholmod_read_sparse(stdin, &c);

    // is this matrix already valid?
    if (Ain->stype == LOWER_TRIANGULAR) {
        std::cerr << "This matrix is already in lower triangular form.\n";
        cholmod_free_sparse(&Ain, &c);
        cholmod_finish(&c);
        return 0;
    }

    cholmod_sort(Ain, &c);  // needed for cholmod_symmetry
    int symm = cholmod_symmetry(Ain, 0, NULL, NULL, NULL, NULL, &c);

    // requires positive diagonals
    // requires symmetry
    // if symm == CHOLMOD_MM_SYMMETRIC_POSDIAG then matrix suitible for LDLT
    if (symm != CHOLMOD_MM_SYMMETRIC_POSDIAG) {
        std::cerr << "Oh no, this matrix is not positive-definate-symmetric.\n";
        cholmod_free_sparse(&Ain, &c);
        cholmod_finish(&c);
        return 0;
    }

    std::cerr << "Great, this matrix is positive-definate-symmetric.\n";
    std::cerr << "Storing as triangular matrix.\n";

    // create new matrix A that stores lower triangular values only
    // value 1 means use values mode
    A = cholmod_copy(Ain, LOWER_TRIANGULAR, 1, &c);

    cholmod_write_sparse(stdout, A, NULL, NULL, &c);
    std::cerr << "Sent to stdout.\n";

    cholmod_free_sparse(&Ain, &c);
    cholmod_free_sparse(&A, &c);
    cholmod_finish(&c);
}
