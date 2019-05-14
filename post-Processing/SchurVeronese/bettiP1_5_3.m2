--This file computes Betti tables for P^2 for d = 5 and b = 3
A := QQ[t_0,t_1,t_2];
--tb stands for Total Betti numbers
tb53 = new HashTable from {(0,0) => 4, (1,0) => 15, (0,1) => 0, (1,1) => 0, (2,0) => 20, (3,0) => 10, (2,1) => 0, (4,0) => 0, (3,1) => 0, (4,1) => 1};
--mb stands for Multigraded Betti numbers
mb53 = new HashTable from {(0,0) => t_0^3+t_0^2*t_1+t_0*t_1^2+t_1^3, (0,1) => 0, (1,0) => t_0^7*t_1+2*t_0^6*t_1^2+3*t_0^5*t_1^3+3*t_0^4*t_1^4+3*t_0^3*t_1^5+2*t_0^2*t_1^6+t_0*t_1^7, (1,1) => 0, (2,0) => t_0^10*t_1^3+2*t_0^9*t_1^4+3*t_0^8*t_1^5+4*t_0^7*t_1^6+4*t_0^6*t_1^7+3*t_0^5*t_1^8+2*t_0^4*t_1^9+t_0^3*t_1^10, (2,1) => 0, (3,0) => t_0^12*t_1^6+t_0^11*t_1^7+2*t_0^10*t_1^8+2*t_0^9*t_1^9+2*t_0^8*t_1^10+t_0^7*t_1^11+t_0^6*t_1^12, (3,1) => 0, (4,0) => 0, (4,1) => t_0^14*t_1^14};
--sb represents the betti numbers as sums of Schur functors
sb53 = new HashTable from {(0,0) => {({3,0},1)}, (1,0) => {({7,1},1)}, (0,1) => {}, (1,1) => {}, (2,0) => {({10,3},1)}, (3,0) => {({12,6},1)}, (2,1) => {}, (4,0) => {}, (3,1) => {}, (4,1) => {({14,14},1)}};
--dw encodes the dominant weights in each entry
dw53 = new HashTable from {(0,0) => {{3,0}}, (1,0) => {{7,1}}, (0,1) => {}, (1,1) => {}, (2,0) => {{10,3}}, (3,0) => {{12,6}}, (2,1) => {}, (4,0) => {}, (3,1) => {}, (4,1) => {{14,14}}};
--lw encodes the lex leading weight in each entry
lw53 = new HashTable from {(0,0) => {3,0}, (1,0) => {7,1}, (0,1) => {}, (1,1) => {}, (2,0) => {10,3}, (3,0) => {12,6}, (2,1) => {}, (4,0) => {}, (3,1) => {}, (4,1) => {14,14}};
--nr encodes the number of disctinct reprsentations in each entry
nr53 = new HashTable from {(0,0) => 1, (1,0) => 1, (0,1) => 0, (1,1) => 0, (2,0) => 1, (3,0) => 1, (2,1) => 0, (4,0) => 0, (3,1) => 0, (4,1) => 1};
--nrm encodes the number of representations with multiplicity in each entry
nrm53 = new HashTable from {(0,0) => 1, (1,0) => 1, (0,1) => 0, (1,1) => 0, (2,0) => 1, (3,0) => 1, (2,1) => 0, (4,0) => 0, (3,1) => 0, (4,1) => 1};
--er encodes the errors in the computed multigraded Hilbert series via our Schur method in each entry
er53 = new HashTable from {(0,0) => 0, (1,0) => 0, (0,1) => 0, (1,1) => 0, (2,0) => 0, (3,0) => 0, (2,1) => 0, (4,0) => 0, (3,1) => 0, (4,1) => 0};
--bs encodes the Boij-Soederberg coefficients each entry
bs53 = {120/1};
end;