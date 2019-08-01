--This file computes Betti tables for P^2 for d = 4 and b = 3
A := QQ[t_0,t_1,t_2];
--tb stands for Total Betti numbers
tb43 = new HashTable from {(0,0) => 4, (1,0) => 12, (0,1) => 0, (1,1) => 0, (2,0) => 12, (3,0) => 4, (2,1) => 0, (3,1) => 0};
--mb stands for Multigraded Betti numbers
mb43 = new HashTable from {(0,0) => t_0^3+t_0^2*t_1+t_0*t_1^2+t_1^3, (0,1) => 0, (1,0) => t_0^6*t_1+2*t_0^5*t_1^2+3*t_0^4*t_1^3+3*t_0^3*t_1^4+2*t_0^2*t_1^5+t_0*t_1^6, (1,1) => 0, (2,0) => t_0^8*t_1^3+2*t_0^7*t_1^4+3*t_0^6*t_1^5+3*t_0^5*t_1^6+2*t_0^4*t_1^7+t_0^3*t_1^8, (2,1) => 0, (3,0) => t_0^9*t_1^6+t_0^8*t_1^7+t_0^7*t_1^8+t_0^6*t_1^9, (3,1) => 0};
--sb represents the betti numbers as sums of Schur functors
sb43 = new HashTable from {(0,0) => {({3,0},1)}, (1,0) => {({6,1},1)}, (0,1) => {}, (1,1) => {}, (2,0) => {({8,3},1)}, (3,0) => {({9,6},1)}, (2,1) => {}, (3,1) => {}};
--dw encodes the dominant weights in each entry
dw43 = new HashTable from {(0,0) => {{3,0}}, (1,0) => {{6,1}}, (0,1) => {}, (1,1) => {}, (2,0) => {{8,3}}, (3,0) => {{9,6}}, (2,1) => {}, (3,1) => {}};
--lw encodes the lex leading weight in each entry
lw43 = new HashTable from {(0,0) => {3,0}, (1,0) => {6,1}, (0,1) => {}, (1,1) => {}, (2,0) => {8,3}, (3,0) => {9,6}, (2,1) => {}, (3,1) => {}};
--nr encodes the number of disctinct reprsentations in each entry
nr43 = new HashTable from {(0,0) => 1, (1,0) => 1, (0,1) => 0, (1,1) => 0, (2,0) => 1, (3,0) => 1, (2,1) => 0, (3,1) => 0};
--nrm encodes the number of representations with multiplicity in each entry
nrm43 = new HashTable from {(0,0) => 1, (1,0) => 1, (0,1) => 0, (1,1) => 0, (2,0) => 1, (3,0) => 1, (2,1) => 0, (3,1) => 0};
--er encodes the errors in the computed multigraded Hilbert series via our Schur method in each entry
er43 = new HashTable from {(0,0) => 0, (1,0) => 0, (0,1) => 0, (1,1) => 0, (2,0) => 0, (3,0) => 0, (2,1) => 0, (3,1) => 0};
--bs encodes the Boij-Soederberg coefficients each entry
bs43 = {24/1};
end;