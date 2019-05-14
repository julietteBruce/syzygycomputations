--This file computes Betti tables for P^2 for d = 3 and b = 2
A := QQ[t_0,t_1,t_2];
--tb stands for Total Betti numbers
tb32 = new HashTable from {(0,0) => 3, (0,1) => 0, (1,0) => 6, (2,0) => 3, (1,1) => 0, (2,1) => 0};
--mb stands for Multigraded Betti numbers
mb32 = new HashTable from {(0,0) => t_0^2+t_0*t_1+t_1^2, (0,1) => 0, (1,0) => t_0^4*t_1+2*t_0^3*t_1^2+2*t_0^2*t_1^3+t_0*t_1^4, (1,1) => 0, (2,0) => t_0^5*t_1^3+t_0^4*t_1^4+t_0^3*t_1^5, (2,1) => 0};
--sb represents the betti numbers as sums of Schur functors
sb32 = new HashTable from {(0,0) => {({2,0},1)}, (0,1) => {}, (1,0) => {({4,1},1)}, (2,0) => {({5,3},1)}, (1,1) => {}, (2,1) => {}};
--dw encodes the dominant weights in each entry
dw32 = new HashTable from {(0,0) => {{2,0}}, (0,1) => {}, (1,0) => {{4,1}}, (2,0) => {{5,3}}, (1,1) => {}, (2,1) => {}};
--lw encodes the lex leading weight in each entry
lw32 = new HashTable from {(0,0) => {2,0}, (0,1) => {}, (1,0) => {4,1}, (2,0) => {5,3}, (1,1) => {}, (2,1) => {}};
--nr encodes the number of disctinct reprsentations in each entry
nr32 = new HashTable from {(0,0) => 1, (0,1) => 0, (1,0) => 1, (2,0) => 1, (1,1) => 0, (2,1) => 0};
--nrm encodes the number of representations with multiplicity in each entry
nrm32 = new HashTable from {(0,0) => 1, (0,1) => 0, (1,0) => 1, (2,0) => 1, (1,1) => 0, (2,1) => 0};
--er encodes the errors in the computed multigraded Hilbert series via our Schur method in each entry
er32 = new HashTable from {(0,0) => 0, (0,1) => 0, (1,0) => 0, (2,0) => 0, (1,1) => 0, (2,1) => 0};
--bs encodes the Boij-Soederberg coefficients each entry
bs32 = {6/1};
end;