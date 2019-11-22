default: all

# Must set SSPATH=path/to/SuiteSparse-X.X.X/
include $(SSPATH)/SuiteSparse_config/SuiteSparse_config.mk

# compiler commands (change this for your compiler)
CC = g++-9 -std=c++17 -Wall

# compiler flags
FLAGS = -I$(SSPATH)/include -L$(SSPATH)/lib -lumfpack -lamd -lcholmod -lsuitesparseconfig -lm

all: make_symmetric
make_symmetric: make_symmetric.cpp
	$(CC) -o make_symmetric.o make_symmetric.cpp $(FLAGS)
