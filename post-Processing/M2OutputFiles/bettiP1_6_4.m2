--This file computes Betti tables for P^2 for d = 6 and b = 4
A := QQ[t_0,t_1,t_2];
--tb stands for Total Betti numbers
tb64 = new HashTable from {(0,0) => 5, (1,0) => 24, (0,1) => 0, (1,1) => 0, (2,0) => 45, (3,0) => 40, (2,1) => 0, (4,0) => 15, (3,1) => 0, (5,0) => 0, (4,1) => 0, (5,1) => 1};
--mb stands for Multigraded Betti numbers
mb64 = new HashTable from {(0,0) => t_0^4+t_0^3*t_1+t_0^2*t_1^2+t_0*t_1^3+t_1^4, (0,1) => 0, (1,0) => t_0^9*t_1+2*t_0^8*t_1^2+3*t_0^7*t_1^3+4*t_0^6*t_1^4+4*t_0^5*t_1^5+4*t_0^4*t_1^6+3*t_0^3*t_1^7+2*t_0^2*t_1^8+t_0*t_1^9, (1,1) => 0, (2,0) => t_0^13*t_1^3+2*t_0^12*t_1^4+4*t_0^11*t_1^5+5*t_0^10*t_1^6+7*t_0^9*t_1^7+7*t_0^8*t_1^8+7*t_0^7*t_1^9+5*t_0^6*t_1^10+4*t_0^5*t_1^11+2*t_0^4*t_1^12+t_0^3*t_1^13, (2,1) => 0, (3,0) => t_0^16*t_1^6+2*t_0^15*t_1^7+3*t_0^14*t_1^8+5*t_0^13*t_1^9+6*t_0^12*t_1^10+6*t_0^11*t_1^11+6*t_0^10*t_1^12+5*t_0^9*t_1^13+3*t_0^8*t_1^14+2*t_0^7*t_1^15+t_0^6*t_1^16, (3,1) => 0, (4,0) => t_0^18*t_1^10+t_0^17*t_1^11+2*t_0^16*t_1^12+2*t_0^15*t_1^13+3*t_0^14*t_1^14+2*t_0^13*t_1^15+2*t_0^12*t_1^16+t_0^11*t_1^17+t_0^10*t_1^18, (4,1) => 0, (5,0) => 0, (5,1) => t_0^20*t_1^20};
--sb represents the betti numbers as sums of Schur functors
sb64 = new HashTable from {(0,0) => {({4,0},1)}, (1,0) => {({9,1},1)}, (0,1) => {}, (1,1) => {}, (2,0) => {({13,3},1)}, (3,0) => {({16,6},1)}, (2,1) => {}, (4,0) => {({18,10},1)}, (3,1) => {}, (5,0) => {}, (4,1) => {}, (5,1) => {({20,20},1)}};
--dw encodes the dominant weights in each entry
dw64 = new HashTable from {(0,0) => {{4,0}}, (1,0) => {{9,1}}, (0,1) => {}, (1,1) => {}, (2,0) => {{13,3}}, (3,0) => {{16,6}}, (2,1) => {}, (4,0) => {{18,10}}, (3,1) => {}, (5,0) => {}, (4,1) => {}, (5,1) => {{20,20}}};
--lw encodes the lex leading weight in each entry
lw64 = new HashTable from {(0,0) => {4,0}, (1,0) => {9,1}, (0,1) => {}, (1,1) => {}, (2,0) => {13,3}, (3,0) => {16,6}, (2,1) => {}, (4,0) => {18,10}, (3,1) => {}, (5,0) => {}, (4,1) => {}, (5,1) => {20,20}};
--nr encodes the number of disctinct reprsentations in each entry
nr64 = new HashTable from {(0,0) => 1, (1,0) => 1, (0,1) => 0, (1,1) => 0, (2,0) => 1, (3,0) => 1, (2,1) => 0, (4,0) => 1, (3,1) => 0, (5,0) => 0, (4,1) => 0, (5,1) => 1};
--nrm encodes the number of representations with multiplicity in each entry
nrm64 = new HashTable from {(0,0) => 1, (1,0) => 1, (0,1) => 0, (1,1) => 0, (2,0) => 1, (3,0) => 1, (2,1) => 0, (4,0) => 1, (3,1) => 0, (5,0) => 0, (4,1) => 0, (5,1) => 1};
--er encodes the errors in the computed multigraded Hilbert series via our Schur method in each entry
er64 = new HashTable from {(0,0) => 0, (1,0) => 0, (0,1) => 0, (1,1) => 0, (2,0) => 0, (3,0) => 0, (2,1) => 0, (4,0) => 0, (3,1) => 0, (5,0) => 0, (4,1) => 0, (5,1) => 0};
--bs encodes the Boij-Soederberg coefficients each entry
bs64 = {720/1};
end;