default: example.o

# Must set SSPATH=path/to/SuiteSparse-X.X.X/
include $(SSPATH)/SuiteSparse_config/SuiteSparse_config.mk

# compiler commands (change this for your compiler)
Cpp = g++-9 -std=c++17 -Wall -O2  # needs C++17
CC = gcc-9 -std=c11 -Wall -O2  # use C11

# compiler flags
FLAGS = -I$(SSPATH)/include -L$(SSPATH)/lib -lcholmod -lsuitesparseconfig -lm

# make object file
solve.o: solve.cpp
	$(Cpp) -o solve.o -c solve.cpp $(FLAGS)

example.o: solve.o
	$(CC) -o example.o solve.o example.c $(FLAGS) -lstdc++

clean:
	rm *.o *.so *.out
