--This file computes Betti tables for P^2 for d = 6 and b = 0
A := QQ[t_0,t_1,t_2];
--tb stands for Total Betti numbers
tb60 = new HashTable from {(0,0) => 1, (1,0) => 0, (0,1) => 0, (1,1) => 15, (2,0) => 0, (3,0) => 0, (2,1) => 40, (4,0) => 0, (3,1) => 45, (5,0) => 0, (4,1) => 24, (5,1) => 5};
--mb stands for Multigraded Betti numbers
mb60 = new HashTable from {(0,0) => 1, (0,1) => 0, (1,0) => 0, (1,1) => t_0^10*t_1^2+t_0^9*t_1^3+2*t_0^8*t_1^4+2*t_0^7*t_1^5+3*t_0^6*t_1^6+2*t_0^5*t_1^7+2*t_0^4*t_1^8+t_0^3*t_1^9+t_0^2*t_1^10, (2,0) => 0, (2,1) => t_0^14*t_1^4+2*t_0^13*t_1^5+3*t_0^12*t_1^6+5*t_0^11*t_1^7+6*t_0^10*t_1^8+6*t_0^9*t_1^9+6*t_0^8*t_1^10+5*t_0^7*t_1^11+3*t_0^6*t_1^12+2*t_0^5*t_1^13+t_0^4*t_1^14, (3,0) => 0, (3,1) => t_0^17*t_1^7+2*t_0^16*t_1^8+4*t_0^15*t_1^9+5*t_0^14*t_1^10+7*t_0^13*t_1^11+7*t_0^12*t_1^12+7*t_0^11*t_1^13+5*t_0^10*t_1^14+4*t_0^9*t_1^15+2*t_0^8*t_1^16+t_0^7*t_1^17, (4,0) => 0, (4,1) => t_0^19*t_1^11+2*t_0^18*t_1^12+3*t_0^17*t_1^13+4*t_0^16*t_1^14+4*t_0^15*t_1^15+4*t_0^14*t_1^16+3*t_0^13*t_1^17+2*t_0^12*t_1^18+t_0^11*t_1^19, (5,0) => 0, (5,1) => t_0^20*t_1^16+t_0^19*t_1^17+t_0^18*t_1^18+t_0^17*t_1^19+t_0^16*t_1^20};
--sb represents the betti numbers as sums of Schur functors
sb60 = new HashTable from {(0,0) => {({0,0},1)}, (1,0) => {}, (0,1) => {}, (1,1) => {({10,2},1)}, (2,0) => {}, (3,0) => {}, (2,1) => {({14,4},1)}, (4,0) => {}, (3,1) => {({17,7},1)}, (5,0) => {}, (4,1) => {({19,11},1)}, (5,1) => {({20,16},1)}};
--dw encodes the dominant weights in each entry
dw60 = new HashTable from {(0,0) => {{0,0}}, (1,0) => {}, (0,1) => {}, (1,1) => {{10,2}}, (2,0) => {}, (3,0) => {}, (2,1) => {{14,4}}, (4,0) => {}, (3,1) => {{17,7}}, (5,0) => {}, (4,1) => {{19,11}}, (5,1) => {{20,16}}};
--lw encodes the lex leading weight in each entry
lw60 = new HashTable from {(0,0) => {0,0}, (1,0) => {}, (0,1) => {}, (1,1) => {10,2}, (2,0) => {}, (3,0) => {}, (2,1) => {14,4}, (4,0) => {}, (3,1) => {17,7}, (5,0) => {}, (4,1) => {19,11}, (5,1) => {20,16}};
--nr encodes the number of disctinct reprsentations in each entry
nr60 = new HashTable from {(0,0) => 1, (1,0) => 0, (0,1) => 0, (1,1) => 1, (2,0) => 0, (3,0) => 0, (2,1) => 1, (4,0) => 0, (3,1) => 1, (5,0) => 0, (4,1) => 1, (5,1) => 1};
--nrm encodes the number of representations with multiplicity in each entry
nrm60 = new HashTable from {(0,0) => 1, (1,0) => 0, (0,1) => 0, (1,1) => 1, (2,0) => 0, (3,0) => 0, (2,1) => 1, (4,0) => 0, (3,1) => 1, (5,0) => 0, (4,1) => 1, (5,1) => 1};
--er encodes the errors in the computed multigraded Hilbert series via our Schur method in each entry
er60 = new HashTable from {(0,0) => 0, (1,0) => 0, (0,1) => 0, (1,1) => 0, (2,0) => 0, (3,0) => 0, (2,1) => 0, (4,0) => 0, (3,1) => 0, (5,0) => 0, (4,1) => 0, (5,1) => 0};
--bs encodes the Boij-Soederberg coefficients each entry
bs60 = {720/1};
end;