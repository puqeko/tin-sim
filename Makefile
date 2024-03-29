default: example.o make_symmetric.o

# Must set SSPATH=path/to/SuiteSparse-X.X.X/
include $(SSPATH)/SuiteSparse_config/SuiteSparse_config.mk

# compiler commands (change this for your compiler)
Cpp = g++-9 -std=c++17 -Wall -O2  # needs C++17
CC = gcc-9 -std=c11 -Wall -O2  # use C11

# compiler flags
FLAGS = -I$(SSPATH)/include -L$(SSPATH)/lib -lcholmod -lsuitesparseconfig -lm

# make object file
solve.o: solve.cpp solve.h
	$(Cpp) -o solve.o -c solve.cpp $(FLAGS)

example.o: solve.o example.c solve.h
	$(CC) -o example.o solve.o example.c $(FLAGS) -lstdc++

make_symmetric.o: make_symmetric.cpp
	$(Cpp) -o make_symmetric.o make_symmetric.cpp $(FLAGS)

clean:
	rm *.o *.so *.out
