
CFLAGS= -Wall -Wextra -g
CXXFLAGS=-std=gnu++11 ${CFLAGS}

all : ConstructMatricies SliceMatrix BasisTest

.PHONY : clean

clean :
	rm -f ConstructMatricies SliceMatrix BasisTest Basis.o ConstructMatrix.o SliceMatrix.o BasisTest.o Combinatorics.o

ConstructMatricies : ConstructMatricies.o Basis.o Combinatorics.o
	g++ ${CXXFLAGS} -o $@ $^

SliceMatrix : SliceMatrix.o Basis.o Combinatorics.o
	g++ ${CXXFLAGS} -o $@ $^

BasisTest : BasisTest.o Basis.o Combinatorics.o
	g++ ${CXXFLAGS} -o $@ $^
