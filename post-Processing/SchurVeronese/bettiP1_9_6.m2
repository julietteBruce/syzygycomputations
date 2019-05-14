--This file computes Betti tables for P^2 for d = 9 and b = 6
A := QQ[t_0,t_1,t_2];
--tb stands for Total Betti numbers
tb96 = new HashTable from {(7,0) => 0, (6,1) => 0, (8,0) => 0, (7,1) => 9, (8,1) => 2, (0,0) => 7, (0,1) => 0, (1,0) => 54, (2,0) => 180, (1,1) => 0, (2,1) => 0, (3,0) => 336, (3,1) => 0, (4,0) => 378, (4,1) => 0, (5,0) => 252, (6,0) => 84, (5,1) => 0};
--mb stands for Multigraded Betti numbers
mb96 = new HashTable from {(7,0) => 0, (6,1) => 0, (7,1) => t_0^43*t_1^35+t_0^42*t_1^36+t_0^41*t_1^37+t_0^40*t_1^38+t_0^39*t_1^39+t_0^38*t_1^40+t_0^37*t_1^41+t_0^36*t_1^42+t_0^35*t_1^43, (8,0) => 0, (8,1) => t_0^44*t_1^43+t_0^43*t_1^44, (0,0) => t_0^6+t_0^5*t_1+t_0^4*t_1^2+t_0^3*t_1^3+t_0^2*t_1^4+t_0*t_1^5+t_1^6, (0,1) => 0, (1,0) => t_0^14*t_1+2*t_0^13*t_1^2+3*t_0^12*t_1^3+4*t_0^11*t_1^4+5*t_0^10*t_1^5+6*t_0^9*t_1^6+6*t_0^8*t_1^7+6*t_0^7*t_1^8+6*t_0^6*t_1^9+5*t_0^5*t_1^10+4*t_0^4*t_1^11+3*t_0^3*t_1^12+2*t_0^2*t_1^13+t_0*t_1^14, (2,0) => t_0^21*t_1^3+2*t_0^20*t_1^4+4*t_0^19*t_1^5+6*t_0^18*t_1^6+9*t_0^17*t_1^7+11*t_0^16*t_1^8+14*t_0^15*t_1^9+16*t_0^14*t_1^10+18*t_0^13*t_1^11+18*t_0^12*t_1^12+18*t_0^11*t_1^13+16*t_0^10*t_1^14+14*t_0^9*t_1^15+11*t_0^8*t_1^16+9*t_0^7*t_1^17+6*t_0^6*t_1^18+4*t_0^5*t_1^19+2*t_0^4*t_1^20+t_0^3*t_1^21, (1,1) => 0, (2,1) => 0, (3,0) => t_0^27*t_1^6+2*t_0^26*t_1^7+4*t_0^25*t_1^8+7*t_0^24*t_1^9+10*t_0^23*t_1^10+14*t_0^22*t_1^11+19*t_0^21*t_1^12+23*t_0^20*t_1^13+27*t_0^19*t_1^14+30*t_0^18*t_1^15+31*t_0^17*t_1^16+31*t_0^16*t_1^17+30*t_0^15*t_1^18+27*t_0^14*t_1^19+23*t_0^13*t_1^20+19*t_0^12*t_1^21+14*t_0^11*t_1^22+10*t_0^10*t_1^23+7*t_0^9*t_1^24+4*t_0^8*t_1^25+2*t_0^7*t_1^26+t_0^6*t_1^27, (4,0) => t_0^32*t_1^10+2*t_0^31*t_1^11+4*t_0^30*t_1^12+6*t_0^29*t_1^13+10*t_0^28*t_1^14+14*t_0^27*t_1^15+19*t_0^26*t_1^16+23*t_0^25*t_1^17+28*t_0^24*t_1^18+31*t_0^23*t_1^19+34*t_0^22*t_1^20+34*t_0^21*t_1^21+34*t_0^20*t_1^22+31*t_0^19*t_1^23+28*t_0^18*t_1^24+23*t_0^17*t_1^25+19*t_0^16*t_1^26+14*t_0^15*t_1^27+10*t_0^14*t_1^28+6*t_0^13*t_1^29+4*t_0^12*t_1^30+2*t_0^11*t_1^31+t_0^10*t_1^32, (3,1) => 0, (5,0) => t_0^36*t_1^15+2*t_0^35*t_1^16+3*t_0^34*t_1^17+5*t_0^33*t_1^18+8*t_0^32*t_1^19+11*t_0^31*t_1^20+14*t_0^30*t_1^21+17*t_0^29*t_1^22+20*t_0^28*t_1^23+22*t_0^27*t_1^24+23*t_0^26*t_1^25+23*t_0^25*t_1^26+22*t_0^24*t_1^27+20*t_0^23*t_1^28+17*t_0^22*t_1^29+14*t_0^21*t_1^30+11*t_0^20*t_1^31+8*t_0^19*t_1^32+5*t_0^18*t_1^33+3*t_0^17*t_1^34+2*t_0^16*t_1^35+t_0^15*t_1^36, (4,1) => 0, (5,1) => 0, (6,0) => t_0^39*t_1^21+t_0^38*t_1^22+2*t_0^37*t_1^23+3*t_0^36*t_1^24+4*t_0^35*t_1^25+5*t_0^34*t_1^26+7*t_0^33*t_1^27+7*t_0^32*t_1^28+8*t_0^31*t_1^29+8*t_0^30*t_1^30+8*t_0^29*t_1^31+7*t_0^28*t_1^32+7*t_0^27*t_1^33+5*t_0^26*t_1^34+4*t_0^25*t_1^35+3*t_0^24*t_1^36+2*t_0^23*t_1^37+t_0^22*t_1^38+t_0^21*t_1^39};
--sb represents the betti numbers as sums of Schur functors
sb96 = new HashTable from {(7,0) => {}, (6,1) => {}, (8,0) => {}, (7,1) => {({43,35},1)}, (8,1) => {({44,43},1)}, (0,0) => {({6,0},1)}, (0,1) => {}, (1,0) => {({14,1},1)}, (2,0) => {({21,3},1)}, (1,1) => {}, (2,1) => {}, (3,0) => {({27,6},1)}, (3,1) => {}, (4,0) => {({32,10},1)}, (4,1) => {}, (5,0) => {({36,15},1)}, (6,0) => {({39,21},1)}, (5,1) => {}};
--dw encodes the dominant weights in each entry
dw96 = new HashTable from {(7,0) => {}, (6,1) => {}, (8,0) => {}, (7,1) => {{43,35}}, (8,1) => {{44,43}}, (0,0) => {{6,0}}, (0,1) => {}, (1,0) => {{14,1}}, (2,0) => {{21,3}}, (1,1) => {}, (2,1) => {}, (3,0) => {{27,6}}, (3,1) => {}, (4,0) => {{32,10}}, (4,1) => {}, (5,0) => {{36,15}}, (6,0) => {{39,21}}, (5,1) => {}};
--lw encodes the lex leading weight in each entry
lw96 = new HashTable from {(7,0) => {}, (6,1) => {}, (8,0) => {}, (7,1) => {43,35}, (8,1) => {44,43}, (0,0) => {6,0}, (0,1) => {}, (1,0) => {14,1}, (2,0) => {21,3}, (1,1) => {}, (2,1) => {}, (3,0) => {27,6}, (3,1) => {}, (4,0) => {32,10}, (4,1) => {}, (5,0) => {36,15}, (6,0) => {39,21}, (5,1) => {}};
--nr encodes the number of disctinct reprsentations in each entry
nr96 = new HashTable from {(7,0) => 0, (6,1) => 0, (8,0) => 0, (7,1) => 1, (8,1) => 1, (0,0) => 1, (0,1) => 0, (1,0) => 1, (2,0) => 1, (1,1) => 0, (2,1) => 0, (3,0) => 1, (3,1) => 0, (4,0) => 1, (4,1) => 0, (5,0) => 1, (6,0) => 1, (5,1) => 0};
--nrm encodes the number of representations with multiplicity in each entry
nrm96 = new HashTable from {(7,0) => 0, (6,1) => 0, (8,0) => 0, (7,1) => 1, (8,1) => 1, (0,0) => 1, (0,1) => 0, (1,0) => 1, (2,0) => 1, (1,1) => 0, (2,1) => 0, (3,0) => 1, (3,1) => 0, (4,0) => 1, (4,1) => 0, (5,0) => 1, (6,0) => 1, (5,1) => 0};
--er encodes the errors in the computed multigraded Hilbert series via our Schur method in each entry
er96 = new HashTable from {(7,0) => 0, (6,1) => 0, (8,0) => 0, (7,1) => 0, (8,1) => 0, (0,0) => 0, (0,1) => 0, (1,0) => 0, (2,0) => 0, (1,1) => 0, (2,1) => 0, (3,0) => 0, (3,1) => 0, (4,0) => 0, (4,1) => 0, (5,0) => 0, (6,0) => 0, (5,1) => 0};
--bs encodes the Boij-Soederberg coefficients each entry
bs96 = {362880/1};
end;