--This file computes Betti tables for P^2 for d = 4 and b = 2
A := QQ[t_0,t_1,t_2];
--tb stands for Total Betti numbers
tb42 = new HashTable from {(0,0) => 3, (1,0) => 8, (0,1) => 0, (1,1) => 0, (2,0) => 6, (3,0) => 0, (2,1) => 0, (3,1) => 1};
--mb stands for Multigraded Betti numbers
mb42 = new HashTable from {(0,0) => t_0^2+t_0*t_1+t_1^2, (0,1) => 0, (1,0) => t_0^5*t_1+2*t_0^4*t_1^2+2*t_0^3*t_1^3+2*t_0^2*t_1^4+t_0*t_1^5, (1,1) => 0, (2,0) => t_0^7*t_1^3+t_0^6*t_1^4+2*t_0^5*t_1^5+t_0^4*t_1^6+t_0^3*t_1^7, (2,1) => 0, (3,0) => 0, (3,1) => t_0^9*t_1^9};
--sb represents the betti numbers as sums of Schur functors
sb42 = new HashTable from {(0,0) => {({2,0},1)}, (1,0) => {({5,1},1)}, (0,1) => {}, (1,1) => {}, (2,0) => {({7,3},1)}, (3,0) => {}, (2,1) => {}, (3,1) => {({9,9},1)}};
--dw encodes the dominant weights in each entry
dw42 = new HashTable from {(0,0) => {{2,0}}, (1,0) => {{5,1}}, (0,1) => {}, (1,1) => {}, (2,0) => {{7,3}}, (3,0) => {}, (2,1) => {}, (3,1) => {{9,9}}};
--lw encodes the lex leading weight in each entry
lw42 = new HashTable from {(0,0) => {2,0}, (1,0) => {5,1}, (0,1) => {}, (1,1) => {}, (2,0) => {7,3}, (3,0) => {}, (2,1) => {}, (3,1) => {9,9}};
--nr encodes the number of disctinct reprsentations in each entry
nr42 = new HashTable from {(0,0) => 1, (1,0) => 1, (0,1) => 0, (1,1) => 0, (2,0) => 1, (3,0) => 0, (2,1) => 0, (3,1) => 1};
--nrm encodes the number of representations with multiplicity in each entry
nrm42 = new HashTable from {(0,0) => 1, (1,0) => 1, (0,1) => 0, (1,1) => 0, (2,0) => 1, (3,0) => 0, (2,1) => 0, (3,1) => 1};
--er encodes the errors in the computed multigraded Hilbert series via our Schur method in each entry
er42 = new HashTable from {(0,0) => 0, (1,0) => 0, (0,1) => 0, (1,1) => 0, (2,0) => 0, (3,0) => 0, (2,1) => 0, (3,1) => 0};
--bs encodes the Boij-Soederberg coefficients each entry
bs42 = {24/1};
end;