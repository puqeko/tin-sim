default: all

# Must set SSPATH=path/to/SuiteSparse-X.X.X/
include $(SSPATH)/SuiteSparse_config/SuiteSparse_config.mk

# compiler commands (change this for your compiler)
CC = g++-9 -std=c++17 -Wall
CCC = gcc-9 -std=c11 -Wall

# compiler flags
FLAGS = -I$(SSPATH)/include -L$(SSPATH)/lib -lumfpack -lamd -lcholmod -lsuitesparseconfig -lm

symmetric: make_symmetric.cpp
	$(CC) -o make_symmetric.o make_symmetric.cpp $(FLAGS)

all: cholesky
cholesky: cholesky.c
	$(CCC) -o cholesky.o cholesky.c $(FLAGS)
