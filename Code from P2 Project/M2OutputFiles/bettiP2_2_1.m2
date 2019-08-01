--This file computes Betti tables for P^2 for d = 2 and b = 1
A := QQ[t_0,t_1,t_2];
--tb stands for Total Betti numbers
tb21 = new HashTable from {(0,0) => 3, (0,1) => 0, (1,0) => 8, (2,0) => 6, (0,2) => 0, (1,1) => 0, (3,0) => 0, (2,1) => 0, (1,2) => 0, (3,1) => 1, (2,2) => 0, (3,2) => 0};
--mb stands for Multigraded Betti numbers
mb21 = new HashTable from {(0,0) => t_0+t_1+t_2, (1,0) => t_0^2*t_1+t_0^2*t_2+t_0*t_1^2+2*t_0*t_1*t_2+t_0*t_2^2+t_1^2*t_2+t_1*t_2^2, (0,1) => 0, (2,0) => t_0^3*t_1*t_2+t_0^2*t_1^2*t_2+t_0^2*t_1*t_2^2+t_0*t_1^3*t_2+t_0*t_1^2*t_2^2+t_0*t_1*t_2^3, (0,2) => 0, (1,1) => 0, (1,2) => 0, (2,1) => 0, (3,0) => 0, (2,2) => 0, (3,1) => t_0^3*t_1^3*t_2^3, (3,2) => 0};
--sb represents the betti numbers as sums of Schur functors
sb21 = new HashTable from {(0,0) => {({1,0,0},1/1)}, (1,0) => {({2,1,0},1/1)}, (0,1) => {}, (2,0) => {({3,1,1},1/1)}, (0,2) => {}, (1,1) => {}, (1,2) => {}, (2,1) => {}, (3,0) => {}, (2,2) => {}, (3,1) => {({3,3,3},1/1)}, (3,2) => {}};
--dw encodes the dominant weights in each entry
dw21 = new HashTable from {(0,0) => {{1,0,0}}, (0,1) => {}, (1,0) => {{2,1,0}}, (2,0) => {{3,1,1}}, (0,2) => {}, (1,1) => {}, (3,0) => {}, (2,1) => {}, (1,2) => {}, (3,1) => {{3,3,3}}, (2,2) => {}, (3,2) => {}};
--lw encodes the lex leading weight in each entry
lw21 = new HashTable from {(0,0) => {1,0,0}, (0,1) => {}, (1,0) => {2,1,0}, (2,0) => {3,1,1}, (0,2) => {}, (1,1) => {}, (3,0) => {}, (2,1) => {}, (1,2) => {}, (3,1) => {3,3,3}, (2,2) => {}, (3,2) => {}};
--nr encodes the number of disctinct reprsentations in each entry
nr21 = new HashTable from {(0,0) => 1, (0,1) => 0, (1,0) => 1, (2,0) => 1, (0,2) => 0, (1,1) => 0, (3,0) => 0, (2,1) => 0, (1,2) => 0, (3,1) => 1, (2,2) => 0, (3,2) => 0};
--nrm encodes the number of representations with multiplicity in each entry
nrm21 = new HashTable from {(0,0) => 1/1, (0,1) => 0, (1,0) => 1/1, (2,0) => 1/1, (0,2) => 0, (1,1) => 0, (3,0) => 0, (2,1) => 0, (1,2) => 0, (3,1) => 1/1, (2,2) => 0, (3,2) => 0};
--er encodes the errors in the computed multigraded Hilbert series via our Schur method in each entry
er21 = new HashTable from {(0,0) => 0, (0,1) => 0, (1,0) => 0, (2,0) => 0, (0,2) => 0, (1,1) => 0, (3,0) => 0, (2,1) => 0, (1,2) => 0, (3,1) => 0, (2,2) => 0, (3,2) => 0};
--bs encodes the Boij-Soederberg coefficients each entry
bs21 = {24/1};
end;