--This file computes Betti tables for P^2 for d = 6 and b = 2
A := QQ[t_0,t_1,t_2];
--tb stands for Total Betti numbers
tb62 = new HashTable from {(0,0) => 3, (1,0) => 12, (0,1) => 0, (1,1) => 0, (2,0) => 15, (3,0) => 0, (2,1) => 0, (4,0) => 0, (3,1) => 15, (5,0) => 0, (4,1) => 12, (5,1) => 3};
--mb stands for Multigraded Betti numbers
mb62 = new HashTable from {(0,0) => t_0^2+t_0*t_1+t_1^2, (0,1) => 0, (1,0) => t_0^7*t_1+2*t_0^6*t_1^2+2*t_0^5*t_1^3+2*t_0^4*t_1^4+2*t_0^3*t_1^5+2*t_0^2*t_1^6+t_0*t_1^7, (1,1) => 0, (2,0) => t_0^11*t_1^3+t_0^10*t_1^4+2*t_0^9*t_1^5+2*t_0^8*t_1^6+3*t_0^7*t_1^7+2*t_0^6*t_1^8+2*t_0^5*t_1^9+t_0^4*t_1^10+t_0^3*t_1^11, (2,1) => 0, (3,0) => 0, (3,1) => t_0^17*t_1^9+t_0^16*t_1^10+2*t_0^15*t_1^11+2*t_0^14*t_1^12+3*t_0^13*t_1^13+2*t_0^12*t_1^14+2*t_0^11*t_1^15+t_0^10*t_1^16+t_0^9*t_1^17, (4,0) => 0, (4,1) => t_0^19*t_1^13+2*t_0^18*t_1^14+2*t_0^17*t_1^15+2*t_0^16*t_1^16+2*t_0^15*t_1^17+2*t_0^14*t_1^18+t_0^13*t_1^19, (5,0) => 0, (5,1) => t_0^20*t_1^18+t_0^19*t_1^19+t_0^18*t_1^20};
--sb represents the betti numbers as sums of Schur functors
sb62 = new HashTable from {(0,0) => {({2,0},1)}, (1,0) => {({7,1},1)}, (0,1) => {}, (1,1) => {}, (2,0) => {({11,3},1)}, (3,0) => {}, (2,1) => {}, (4,0) => {}, (3,1) => {({17,9},1)}, (5,0) => {}, (4,1) => {({19,13},1)}, (5,1) => {({20,18},1)}};
--dw encodes the dominant weights in each entry
dw62 = new HashTable from {(0,0) => {{2,0}}, (1,0) => {{7,1}}, (0,1) => {}, (1,1) => {}, (2,0) => {{11,3}}, (3,0) => {}, (2,1) => {}, (4,0) => {}, (3,1) => {{17,9}}, (5,0) => {}, (4,1) => {{19,13}}, (5,1) => {{20,18}}};
--lw encodes the lex leading weight in each entry
lw62 = new HashTable from {(0,0) => {2,0}, (1,0) => {7,1}, (0,1) => {}, (1,1) => {}, (2,0) => {11,3}, (3,0) => {}, (2,1) => {}, (4,0) => {}, (3,1) => {17,9}, (5,0) => {}, (4,1) => {19,13}, (5,1) => {20,18}};
--nr encodes the number of disctinct reprsentations in each entry
nr62 = new HashTable from {(0,0) => 1, (1,0) => 1, (0,1) => 0, (1,1) => 0, (2,0) => 1, (3,0) => 0, (2,1) => 0, (4,0) => 0, (3,1) => 1, (5,0) => 0, (4,1) => 1, (5,1) => 1};
--nrm encodes the number of representations with multiplicity in each entry
nrm62 = new HashTable from {(0,0) => 1, (1,0) => 1, (0,1) => 0, (1,1) => 0, (2,0) => 1, (3,0) => 0, (2,1) => 0, (4,0) => 0, (3,1) => 1, (5,0) => 0, (4,1) => 1, (5,1) => 1};
--er encodes the errors in the computed multigraded Hilbert series via our Schur method in each entry
er62 = new HashTable from {(0,0) => 0, (1,0) => 0, (0,1) => 0, (1,1) => 0, (2,0) => 0, (3,0) => 0, (2,1) => 0, (4,0) => 0, (3,1) => 0, (5,0) => 0, (4,1) => 0, (5,1) => 0};
--bs encodes the Boij-Soederberg coefficients each entry
bs62 = {720/1};
end;