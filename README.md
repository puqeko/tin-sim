Testing matrix factorisation methods as part of a deterministic simulation of electrical interactions between tin particles.

The methods examined are:

* Cholosky Factorisation
* Cholosky Factorisation with L matrix updates
* + Other stuff to be done...

## Installation

* Download, compile, and install [SuiteSparse](http://faculty.cse.tamu.edu/davis/suitesparse.html). Only the CHOLMOD module (and its dependances) is used.
* Set the environment variable SSPATH="/path/to/SuiteSparse-X.X.X". You may also need to change the compiler command you are using at the `CC=..` line in the MakeFile. I used GNU.
* Run `make` to compile.
