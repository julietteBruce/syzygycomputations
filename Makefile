
CFLAGS= -Wall -Wextra -g -O2
CXXFLAGS=-std=gnu++11 $(CFLAGS)

srcdir=src/CreateMatricies
builddir=build

vpath %.h src/CreateMatricies/
vpath %.c src/CreateMatricies/
vpath %.cpp src/CreateMatricies/

TARGETS=$(addprefix $(builddir)/, ConstructMatricies SliceMatrix BasisTest DirectConstructMatrices)

all : $(TARGETS)

.PHONY : clean

$(builddir)/%.o : %.cpp | $(builddir)
	$(CXX) $(CXXFLAGS) -c -o $@ $<

$(builddir)/ConstructMatricies : $(addprefix $(builddir)/, ConstructMatricies.o Basis.o Combinatorics.o)
	$(CXX) $(CXXFLAGS) -o $@ $^

$(builddir)/SliceMatrix : $(addprefix $(builddir)/, SliceMatrix.o Basis.o Combinatorics.o)
	$(CXX) $(CXXFLAGS) -o $@ $^

$(builddir)/DirectConstructMatrices : $(addprefix $(builddir)/, DirectConstructMatrices.o Basis.o Combinatorics.o)
	$(CXX) $(CXXFLAGS) -o $@ $^

$(builddir)/BasisTest : $(addprefix $(builddir)/, BasisTest.o Basis.o Combinatorics.o)
	$(CXX) $(CXXFLAGS) -o $@ $^

clean :
	rm -rf $(builddir)

$(builddir):
	mkdir $(builddir)
