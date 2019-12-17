default: main

# Must set SSPATH=path/to/SuiteSparse-X.X.X/
include $(SSPATH)/SuiteSparse_config/SuiteSparse_config.mk

# compiler commands (change this for your compiler)
Cpp = g++-9 -std=c++17 -Wall -O2
CC = gcc-9 -std=c11 -Wall -O2

# compiler flags
FLAGS = -I$(SSPATH)/include -L$(SSPATH)/lib -lumfpack -lamd -lcholmod -lsuitesparseconfig -lm

solve.o: solve.cpp solve.h
	$(Cpp) -c -o solve.o solve.cpp $(FLAGS)

main.o: solve.o solve.h
	$(CC) -o main.o solve.o $(FLAGS)

clean:
	rm main.o solve.o
