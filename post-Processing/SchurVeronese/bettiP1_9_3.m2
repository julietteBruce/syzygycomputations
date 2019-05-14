--This file computes Betti tables for P^2 for d = 9 and b = 3
A := QQ[t_0,t_1,t_2];
--tb stands for Total Betti numbers
tb93 = new HashTable from {(7,0) => 0, (6,1) => 108, (8,0) => 0, (7,1) => 36, (8,1) => 5, (0,0) => 4, (0,1) => 0, (1,0) => 27, (2,0) => 72, (1,1) => 0, (2,1) => 0, (3,0) => 84, (3,1) => 0, (4,0) => 0, (4,1) => 126, (5,0) => 0, (6,0) => 0, (5,1) => 168};
--mb stands for Multigraded Betti numbers
mb93 = new HashTable from {(7,0) => 0, (6,1) => t_0^41*t_1^25+2*t_0^40*t_1^26+4*t_0^39*t_1^27+5*t_0^38*t_1^28+7*t_0^37*t_1^29+8*t_0^36*t_1^30+10*t_0^35*t_1^31+11*t_0^34*t_1^32+12*t_0^33*t_1^33+11*t_0^32*t_1^34+10*t_0^31*t_1^35+8*t_0^30*t_1^36+7*t_0^29*t_1^37+5*t_0^28*t_1^38+4*t_0^27*t_1^39+2*t_0^26*t_1^40+t_0^25*t_1^41, (7,1) => t_0^43*t_1^32+2*t_0^42*t_1^33+3*t_0^41*t_1^34+4*t_0^40*t_1^35+4*t_0^39*t_1^36+4*t_0^38*t_1^37+4*t_0^37*t_1^38+4*t_0^36*t_1^39+4*t_0^35*t_1^40+3*t_0^34*t_1^41+2*t_0^33*t_1^42+t_0^32*t_1^43, (8,0) => 0, (8,1) => t_0^44*t_1^40+t_0^43*t_1^41+t_0^42*t_1^42+t_0^41*t_1^43+t_0^40*t_1^44, (0,0) => t_0^3+t_0^2*t_1+t_0*t_1^2+t_1^3, (0,1) => 0, (1,0) => t_0^11*t_1+2*t_0^10*t_1^2+3*t_0^9*t_1^3+3*t_0^8*t_1^4+3*t_0^7*t_1^5+3*t_0^6*t_1^6+3*t_0^5*t_1^7+3*t_0^4*t_1^8+3*t_0^3*t_1^9+2*t_0^2*t_1^10+t_0*t_1^11, (2,0) => t_0^18*t_1^3+2*t_0^17*t_1^4+3*t_0^16*t_1^5+4*t_0^15*t_1^6+5*t_0^14*t_1^7+6*t_0^13*t_1^8+7*t_0^12*t_1^9+8*t_0^11*t_1^10+8*t_0^10*t_1^11+7*t_0^9*t_1^12+6*t_0^8*t_1^13+5*t_0^7*t_1^14+4*t_0^6*t_1^15+3*t_0^5*t_1^16+2*t_0^4*t_1^17+t_0^3*t_1^18, (1,1) => 0, (2,1) => 0, (3,0) => t_0^24*t_1^6+t_0^23*t_1^7+2*t_0^22*t_1^8+3*t_0^21*t_1^9+4*t_0^20*t_1^10+5*t_0^19*t_1^11+7*t_0^18*t_1^12+7*t_0^17*t_1^13+8*t_0^16*t_1^14+8*t_0^15*t_1^15+8*t_0^14*t_1^16+7*t_0^13*t_1^17+7*t_0^12*t_1^18+5*t_0^11*t_1^19+4*t_0^10*t_1^20+3*t_0^9*t_1^21+2*t_0^8*t_1^22+t_0^7*t_1^23+t_0^6*t_1^24, (4,0) => 0, (3,1) => 0, (5,0) => 0, (4,1) => t_0^34*t_1^14+t_0^33*t_1^15+2*t_0^32*t_1^16+3*t_0^31*t_1^17+5*t_0^30*t_1^18+6*t_0^29*t_1^19+8*t_0^28*t_1^20+9*t_0^27*t_1^21+11*t_0^26*t_1^22+11*t_0^25*t_1^23+12*t_0^24*t_1^24+11*t_0^23*t_1^25+11*t_0^22*t_1^26+9*t_0^21*t_1^27+8*t_0^20*t_1^28+6*t_0^19*t_1^29+5*t_0^18*t_1^30+3*t_0^17*t_1^31+2*t_0^16*t_1^32+t_0^15*t_1^33+t_0^14*t_1^34, (5,1) => t_0^38*t_1^19+2*t_0^37*t_1^20+3*t_0^36*t_1^21+5*t_0^35*t_1^22+7*t_0^34*t_1^23+9*t_0^33*t_1^24+12*t_0^32*t_1^25+14*t_0^31*t_1^26+15*t_0^30*t_1^27+16*t_0^29*t_1^28+16*t_0^28*t_1^29+15*t_0^27*t_1^30+14*t_0^26*t_1^31+12*t_0^25*t_1^32+9*t_0^24*t_1^33+7*t_0^23*t_1^34+5*t_0^22*t_1^35+3*t_0^21*t_1^36+2*t_0^20*t_1^37+t_0^19*t_1^38, (6,0) => 0};
--sb represents the betti numbers as sums of Schur functors
sb93 = new HashTable from {(7,0) => {}, (6,1) => {({41,25},1)}, (8,0) => {}, (7,1) => {({43,32},1)}, (8,1) => {({44,40},1)}, (0,0) => {({3,0},1)}, (0,1) => {}, (1,0) => {({11,1},1)}, (2,0) => {({18,3},1)}, (1,1) => {}, (2,1) => {}, (3,0) => {({24,6},1)}, (3,1) => {}, (4,0) => {}, (4,1) => {({34,14},1)}, (5,0) => {}, (6,0) => {}, (5,1) => {({38,19},1)}};
--dw encodes the dominant weights in each entry
dw93 = new HashTable from {(7,0) => {}, (6,1) => {{41,25}}, (8,0) => {}, (7,1) => {{43,32}}, (8,1) => {{44,40}}, (0,0) => {{3,0}}, (0,1) => {}, (1,0) => {{11,1}}, (2,0) => {{18,3}}, (1,1) => {}, (2,1) => {}, (3,0) => {{24,6}}, (3,1) => {}, (4,0) => {}, (4,1) => {{34,14}}, (5,0) => {}, (6,0) => {}, (5,1) => {{38,19}}};
--lw encodes the lex leading weight in each entry
lw93 = new HashTable from {(7,0) => {}, (6,1) => {41,25}, (8,0) => {}, (7,1) => {43,32}, (8,1) => {44,40}, (0,0) => {3,0}, (0,1) => {}, (1,0) => {11,1}, (2,0) => {18,3}, (1,1) => {}, (2,1) => {}, (3,0) => {24,6}, (3,1) => {}, (4,0) => {}, (4,1) => {34,14}, (5,0) => {}, (6,0) => {}, (5,1) => {38,19}};
--nr encodes the number of disctinct reprsentations in each entry
nr93 = new HashTable from {(7,0) => 0, (6,1) => 1, (8,0) => 0, (7,1) => 1, (8,1) => 1, (0,0) => 1, (0,1) => 0, (1,0) => 1, (2,0) => 1, (1,1) => 0, (2,1) => 0, (3,0) => 1, (3,1) => 0, (4,0) => 0, (4,1) => 1, (5,0) => 0, (6,0) => 0, (5,1) => 1};
--nrm encodes the number of representations with multiplicity in each entry
nrm93 = new HashTable from {(7,0) => 0, (6,1) => 1, (8,0) => 0, (7,1) => 1, (8,1) => 1, (0,0) => 1, (0,1) => 0, (1,0) => 1, (2,0) => 1, (1,1) => 0, (2,1) => 0, (3,0) => 1, (3,1) => 0, (4,0) => 0, (4,1) => 1, (5,0) => 0, (6,0) => 0, (5,1) => 1};
--er encodes the errors in the computed multigraded Hilbert series via our Schur method in each entry
er93 = new HashTable from {(7,0) => 0, (6,1) => 0, (8,0) => 0, (7,1) => 0, (8,1) => 0, (0,0) => 0, (0,1) => 0, (1,0) => 0, (2,0) => 0, (1,1) => 0, (2,1) => 0, (3,0) => 0, (3,1) => 0, (4,0) => 0, (4,1) => 0, (5,0) => 0, (6,0) => 0, (5,1) => 0};
--bs encodes the Boij-Soederberg coefficients each entry
bs93 = {362880/1};
end;