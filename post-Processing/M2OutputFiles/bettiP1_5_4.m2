--This file computes Betti tables for P^2 for d = 5 and b = 4
A := QQ[t_0,t_1,t_2];
--tb stands for Total Betti numbers
tb54 = new HashTable from {(0,0) => 5, (1,0) => 20, (0,1) => 0, (1,1) => 0, (2,0) => 30, (3,0) => 20, (2,1) => 0, (4,0) => 5, (3,1) => 0, (4,1) => 0};
--mb stands for Multigraded Betti numbers
mb54 = new HashTable from {(0,0) => t_0^4+t_0^3*t_1+t_0^2*t_1^2+t_0*t_1^3+t_1^4, (0,1) => 0, (1,0) => t_0^8*t_1+2*t_0^7*t_1^2+3*t_0^6*t_1^3+4*t_0^5*t_1^4+4*t_0^4*t_1^5+3*t_0^3*t_1^6+2*t_0^2*t_1^7+t_0*t_1^8, (1,1) => 0, (2,0) => t_0^11*t_1^3+2*t_0^10*t_1^4+4*t_0^9*t_1^5+5*t_0^8*t_1^6+6*t_0^7*t_1^7+5*t_0^6*t_1^8+4*t_0^5*t_1^9+2*t_0^4*t_1^10+t_0^3*t_1^11, (2,1) => 0, (3,0) => t_0^13*t_1^6+2*t_0^12*t_1^7+3*t_0^11*t_1^8+4*t_0^10*t_1^9+4*t_0^9*t_1^10+3*t_0^8*t_1^11+2*t_0^7*t_1^12+t_0^6*t_1^13, (3,1) => 0, (4,0) => t_0^14*t_1^10+t_0^13*t_1^11+t_0^12*t_1^12+t_0^11*t_1^13+t_0^10*t_1^14, (4,1) => 0};
--sb represents the betti numbers as sums of Schur functors
sb54 = new HashTable from {(0,0) => {({4,0},1)}, (1,0) => {({8,1},1)}, (0,1) => {}, (1,1) => {}, (2,0) => {({11,3},1)}, (3,0) => {({13,6},1)}, (2,1) => {}, (4,0) => {({14,10},1)}, (3,1) => {}, (4,1) => {}};
--dw encodes the dominant weights in each entry
dw54 = new HashTable from {(0,0) => {{4,0}}, (1,0) => {{8,1}}, (0,1) => {}, (1,1) => {}, (2,0) => {{11,3}}, (3,0) => {{13,6}}, (2,1) => {}, (4,0) => {{14,10}}, (3,1) => {}, (4,1) => {}};
--lw encodes the lex leading weight in each entry
lw54 = new HashTable from {(0,0) => {4,0}, (1,0) => {8,1}, (0,1) => {}, (1,1) => {}, (2,0) => {11,3}, (3,0) => {13,6}, (2,1) => {}, (4,0) => {14,10}, (3,1) => {}, (4,1) => {}};
--nr encodes the number of disctinct reprsentations in each entry
nr54 = new HashTable from {(0,0) => 1, (1,0) => 1, (0,1) => 0, (1,1) => 0, (2,0) => 1, (3,0) => 1, (2,1) => 0, (4,0) => 1, (3,1) => 0, (4,1) => 0};
--nrm encodes the number of representations with multiplicity in each entry
nrm54 = new HashTable from {(0,0) => 1, (1,0) => 1, (0,1) => 0, (1,1) => 0, (2,0) => 1, (3,0) => 1, (2,1) => 0, (4,0) => 1, (3,1) => 0, (4,1) => 0};
--er encodes the errors in the computed multigraded Hilbert series via our Schur method in each entry
er54 = new HashTable from {(0,0) => 0, (1,0) => 0, (0,1) => 0, (1,1) => 0, (2,0) => 0, (3,0) => 0, (2,1) => 0, (4,0) => 0, (3,1) => 0, (4,1) => 0};
--bs encodes the Boij-Soederberg coefficients each entry
bs54 = {120/1};
end;