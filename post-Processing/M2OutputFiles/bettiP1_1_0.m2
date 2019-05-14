--This file computes Betti tables for P^1 for d = 1 and b = 0
--
--Need to initialize ring A for multigraded Hilbert series
A = QQ[t_0,t_1];
--
--tb stands for Total Betti numbers
tb10 = new HashTable from {(0,0) => 1, (0,1) => 0};
--mb stands for Multigraded Betti numbers
mb10 = new HashTable from {(0,0) => 1, (0,1) => 0};
--sb represents the betti numbers as sums of Schur functors
sb10 = new HashTable from {(0,0) => {({0,0},1)}, (0,1) => {}};
--dw encodes the dominant weights in each entry
dw10 = new HashTable from {(0,0) => {{0,0}}, (0,1) => {}};
--lw encodes the lex leading weight in each entry
lw10 = new HashTable from {(0,0) => {0,0}, (0,1) => {}};
--nr encodes the number of disctinct reprsentations in each entry
nr10 = new HashTable from {(0,0) => 1, (0,1) => 0};
--nrm encodes the number of representations with multiplicity in each entry
nrm10 = new HashTable from {(0,0) => 1, (0,1) => 0};
--er encodes the errors in the computed multigraded Hilbert series via our Schur method in each entry
er10 = new HashTable from {(0,0) => 0, (0,1) => 0};
end;