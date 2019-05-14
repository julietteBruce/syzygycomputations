--This file computes Betti tables for P^2 for d = 6 and b = 5
A := QQ[t_0,t_1,t_2];
--tb stands for Total Betti numbers
tb65 = new HashTable from {(0,0) => 6, (1,0) => 30, (0,1) => 0, (1,1) => 0, (2,0) => 60, (3,0) => 60, (2,1) => 0, (4,0) => 30, (3,1) => 0, (5,0) => 6, (4,1) => 0, (5,1) => 0};
--mb stands for Multigraded Betti numbers
mb65 = new HashTable from {(0,0) => t_0^5+t_0^4*t_1+t_0^3*t_1^2+t_0^2*t_1^3+t_0*t_1^4+t_1^5, (0,1) => 0, (1,0) => t_0^10*t_1+2*t_0^9*t_1^2+3*t_0^8*t_1^3+4*t_0^7*t_1^4+5*t_0^6*t_1^5+5*t_0^5*t_1^6+4*t_0^4*t_1^7+3*t_0^3*t_1^8+2*t_0^2*t_1^9+t_0*t_1^10, (1,1) => 0, (2,0) => t_0^14*t_1^3+2*t_0^13*t_1^4+4*t_0^12*t_1^5+6*t_0^11*t_1^6+8*t_0^10*t_1^7+9*t_0^9*t_1^8+9*t_0^8*t_1^9+8*t_0^7*t_1^10+6*t_0^6*t_1^11+4*t_0^5*t_1^12+2*t_0^4*t_1^13+t_0^3*t_1^14, (2,1) => 0, (3,0) => t_0^17*t_1^6+2*t_0^16*t_1^7+4*t_0^15*t_1^8+6*t_0^14*t_1^9+8*t_0^13*t_1^10+9*t_0^12*t_1^11+9*t_0^11*t_1^12+8*t_0^10*t_1^13+6*t_0^9*t_1^14+4*t_0^8*t_1^15+2*t_0^7*t_1^16+t_0^6*t_1^17, (3,1) => 0, (4,0) => t_0^19*t_1^10+2*t_0^18*t_1^11+3*t_0^17*t_1^12+4*t_0^16*t_1^13+5*t_0^15*t_1^14+5*t_0^14*t_1^15+4*t_0^13*t_1^16+3*t_0^12*t_1^17+2*t_0^11*t_1^18+t_0^10*t_1^19, (4,1) => 0, (5,0) => t_0^20*t_1^15+t_0^19*t_1^16+t_0^18*t_1^17+t_0^17*t_1^18+t_0^16*t_1^19+t_0^15*t_1^20, (5,1) => 0};
--sb represents the betti numbers as sums of Schur functors
sb65 = new HashTable from {(0,0) => {({5,0},1)}, (1,0) => {({10,1},1)}, (0,1) => {}, (1,1) => {}, (2,0) => {({14,3},1)}, (3,0) => {({17,6},1)}, (2,1) => {}, (4,0) => {({19,10},1)}, (3,1) => {}, (5,0) => {({20,15},1)}, (4,1) => {}, (5,1) => {}};
--dw encodes the dominant weights in each entry
dw65 = new HashTable from {(0,0) => {{5,0}}, (1,0) => {{10,1}}, (0,1) => {}, (1,1) => {}, (2,0) => {{14,3}}, (3,0) => {{17,6}}, (2,1) => {}, (4,0) => {{19,10}}, (3,1) => {}, (5,0) => {{20,15}}, (4,1) => {}, (5,1) => {}};
--lw encodes the lex leading weight in each entry
lw65 = new HashTable from {(0,0) => {5,0}, (1,0) => {10,1}, (0,1) => {}, (1,1) => {}, (2,0) => {14,3}, (3,0) => {17,6}, (2,1) => {}, (4,0) => {19,10}, (3,1) => {}, (5,0) => {20,15}, (4,1) => {}, (5,1) => {}};
--nr encodes the number of disctinct reprsentations in each entry
nr65 = new HashTable from {(0,0) => 1, (1,0) => 1, (0,1) => 0, (1,1) => 0, (2,0) => 1, (3,0) => 1, (2,1) => 0, (4,0) => 1, (3,1) => 0, (5,0) => 1, (4,1) => 0, (5,1) => 0};
--nrm encodes the number of representations with multiplicity in each entry
nrm65 = new HashTable from {(0,0) => 1, (1,0) => 1, (0,1) => 0, (1,1) => 0, (2,0) => 1, (3,0) => 1, (2,1) => 0, (4,0) => 1, (3,1) => 0, (5,0) => 1, (4,1) => 0, (5,1) => 0};
--er encodes the errors in the computed multigraded Hilbert series via our Schur method in each entry
er65 = new HashTable from {(0,0) => 0, (1,0) => 0, (0,1) => 0, (1,1) => 0, (2,0) => 0, (3,0) => 0, (2,1) => 0, (4,0) => 0, (3,1) => 0, (5,0) => 0, (4,1) => 0, (5,1) => 0};
--bs encodes the Boij-Soederberg coefficients each entry
bs65 = {720/1};
end;