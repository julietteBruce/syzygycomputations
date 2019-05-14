--This file computes Betti tables for P^2 for d = 3 and b = 0
A := QQ[t_0,t_1,t_2];
--tb stands for Total Betti numbers
tb30 = new HashTable from {(0,0) => 1, (0,1) => 0, (1,0) => 0, (2,0) => 0, (1,1) => 3, (2,1) => 2};
--mb stands for Multigraded Betti numbers
mb30 = new HashTable from {(0,0) => 1, (0,1) => 0, (1,0) => 0, (1,1) => t_0^4*t_1^2+t_0^3*t_1^3+t_0^2*t_1^4, (2,0) => 0, (2,1) => t_0^5*t_1^4+t_0^4*t_1^5};
--sb represents the betti numbers as sums of Schur functors
sb30 = new HashTable from {(0,0) => {({0,0},1)}, (0,1) => {}, (1,0) => {}, (2,0) => {}, (1,1) => {({4,2},1)}, (2,1) => {({5,4},1)}};
--dw encodes the dominant weights in each entry
dw30 = new HashTable from {(0,0) => {{0,0}}, (0,1) => {}, (1,0) => {}, (2,0) => {}, (1,1) => {{4,2}}, (2,1) => {{5,4}}};
--lw encodes the lex leading weight in each entry
lw30 = new HashTable from {(0,0) => {0,0}, (0,1) => {}, (1,0) => {}, (2,0) => {}, (1,1) => {4,2}, (2,1) => {5,4}};
--nr encodes the number of disctinct reprsentations in each entry
nr30 = new HashTable from {(0,0) => 1, (0,1) => 0, (1,0) => 0, (2,0) => 0, (1,1) => 1, (2,1) => 1};
--nrm encodes the number of representations with multiplicity in each entry
nrm30 = new HashTable from {(0,0) => 1, (0,1) => 0, (1,0) => 0, (2,0) => 0, (1,1) => 1, (2,1) => 1};
--er encodes the errors in the computed multigraded Hilbert series via our Schur method in each entry
er30 = new HashTable from {(0,0) => 0, (0,1) => 0, (1,0) => 0, (2,0) => 0, (1,1) => 0, (2,1) => 0};
--bs encodes the Boij-Soederberg coefficients each entry
bs30 = {6/1};
end;