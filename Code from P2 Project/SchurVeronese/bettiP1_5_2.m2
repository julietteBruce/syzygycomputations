--This file computes Betti tables for P^2 for d = 5 and b = 2
A := QQ[t_0,t_1,t_2];
--tb stands for Total Betti numbers
tb52 = new HashTable from {(0,0) => 3, (1,0) => 10, (0,1) => 0, (1,1) => 0, (2,0) => 10, (3,0) => 0, (2,1) => 0, (4,0) => 0, (3,1) => 5, (4,1) => 2};
--mb stands for Multigraded Betti numbers
mb52 = new HashTable from {(0,0) => t_0^2+t_0*t_1+t_1^2, (0,1) => 0, (1,0) => t_0^6*t_1+2*t_0^5*t_1^2+2*t_0^4*t_1^3+2*t_0^3*t_1^4+2*t_0^2*t_1^5+t_0*t_1^6, (1,1) => 0, (2,0) => t_0^9*t_1^3+t_0^8*t_1^4+2*t_0^7*t_1^5+2*t_0^6*t_1^6+2*t_0^5*t_1^7+t_0^4*t_1^8+t_0^3*t_1^9, (2,1) => 0, (3,0) => 0, (3,1) => t_0^13*t_1^9+t_0^12*t_1^10+t_0^11*t_1^11+t_0^10*t_1^12+t_0^9*t_1^13, (4,0) => 0, (4,1) => t_0^14*t_1^13+t_0^13*t_1^14};
--sb represents the betti numbers as sums of Schur functors
sb52 = new HashTable from {(0,0) => {({2,0},1)}, (1,0) => {({6,1},1)}, (0,1) => {}, (1,1) => {}, (2,0) => {({9,3},1)}, (3,0) => {}, (2,1) => {}, (4,0) => {}, (3,1) => {({13,9},1)}, (4,1) => {({14,13},1)}};
--dw encodes the dominant weights in each entry
dw52 = new HashTable from {(0,0) => {{2,0}}, (1,0) => {{6,1}}, (0,1) => {}, (1,1) => {}, (2,0) => {{9,3}}, (3,0) => {}, (2,1) => {}, (4,0) => {}, (3,1) => {{13,9}}, (4,1) => {{14,13}}};
--lw encodes the lex leading weight in each entry
lw52 = new HashTable from {(0,0) => {2,0}, (1,0) => {6,1}, (0,1) => {}, (1,1) => {}, (2,0) => {9,3}, (3,0) => {}, (2,1) => {}, (4,0) => {}, (3,1) => {13,9}, (4,1) => {14,13}};
--nr encodes the number of disctinct reprsentations in each entry
nr52 = new HashTable from {(0,0) => 1, (1,0) => 1, (0,1) => 0, (1,1) => 0, (2,0) => 1, (3,0) => 0, (2,1) => 0, (4,0) => 0, (3,1) => 1, (4,1) => 1};
--nrm encodes the number of representations with multiplicity in each entry
nrm52 = new HashTable from {(0,0) => 1, (1,0) => 1, (0,1) => 0, (1,1) => 0, (2,0) => 1, (3,0) => 0, (2,1) => 0, (4,0) => 0, (3,1) => 1, (4,1) => 1};
--er encodes the errors in the computed multigraded Hilbert series via our Schur method in each entry
er52 = new HashTable from {(0,0) => 0, (1,0) => 0, (0,1) => 0, (1,1) => 0, (2,0) => 0, (3,0) => 0, (2,1) => 0, (4,0) => 0, (3,1) => 0, (4,1) => 0};
--bs encodes the Boij-Soederberg coefficients each entry
bs52 = {120/1};
end;