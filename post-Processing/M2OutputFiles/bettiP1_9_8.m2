--This file computes Betti tables for P^2 for d = 9 and b = 8
A := QQ[t_0,t_1,t_2];
--tb stands for Total Betti numbers
tb98 = new HashTable from {(7,0) => 72, (6,1) => 0, (8,0) => 9, (7,1) => 0, (8,1) => 0, (0,0) => 9, (0,1) => 0, (1,0) => 72, (2,0) => 252, (1,1) => 0, (2,1) => 0, (3,0) => 504, (3,1) => 0, (4,0) => 630, (4,1) => 0, (5,0) => 504, (6,0) => 252, (5,1) => 0};
--mb stands for Multigraded Betti numbers
mb98 = new HashTable from {(7,0) => t_0^43*t_1^28+2*t_0^42*t_1^29+3*t_0^41*t_1^30+4*t_0^40*t_1^31+5*t_0^39*t_1^32+6*t_0^38*t_1^33+7*t_0^37*t_1^34+8*t_0^36*t_1^35+8*t_0^35*t_1^36+7*t_0^34*t_1^37+6*t_0^33*t_1^38+5*t_0^32*t_1^39+4*t_0^31*t_1^40+3*t_0^30*t_1^41+2*t_0^29*t_1^42+t_0^28*t_1^43, (6,1) => 0, (7,1) => 0, (8,0) => t_0^44*t_1^36+t_0^43*t_1^37+t_0^42*t_1^38+t_0^41*t_1^39+t_0^40*t_1^40+t_0^39*t_1^41+t_0^38*t_1^42+t_0^37*t_1^43+t_0^36*t_1^44, (8,1) => 0, (0,0) => t_0^8+t_0^7*t_1+t_0^6*t_1^2+t_0^5*t_1^3+t_0^4*t_1^4+t_0^3*t_1^5+t_0^2*t_1^6+t_0*t_1^7+t_1^8, (0,1) => 0, (1,0) => t_0^16*t_1+2*t_0^15*t_1^2+3*t_0^14*t_1^3+4*t_0^13*t_1^4+5*t_0^12*t_1^5+6*t_0^11*t_1^6+7*t_0^10*t_1^7+8*t_0^9*t_1^8+8*t_0^8*t_1^9+7*t_0^7*t_1^10+6*t_0^6*t_1^11+5*t_0^5*t_1^12+4*t_0^4*t_1^13+3*t_0^3*t_1^14+2*t_0^2*t_1^15+t_0*t_1^16, (2,0) => t_0^23*t_1^3+2*t_0^22*t_1^4+4*t_0^21*t_1^5+6*t_0^20*t_1^6+9*t_0^19*t_1^7+12*t_0^18*t_1^8+16*t_0^17*t_1^9+19*t_0^16*t_1^10+22*t_0^15*t_1^11+23*t_0^14*t_1^12+24*t_0^13*t_1^13+23*t_0^12*t_1^14+22*t_0^11*t_1^15+19*t_0^10*t_1^16+16*t_0^9*t_1^17+12*t_0^8*t_1^18+9*t_0^7*t_1^19+6*t_0^6*t_1^20+4*t_0^5*t_1^21+2*t_0^4*t_1^22+t_0^3*t_1^23, (1,1) => 0, (2,1) => 0, (3,0) => t_0^29*t_1^6+2*t_0^28*t_1^7+4*t_0^27*t_1^8+7*t_0^26*t_1^9+11*t_0^25*t_1^10+16*t_0^24*t_1^11+22*t_0^23*t_1^12+28*t_0^22*t_1^13+34*t_0^21*t_1^14+39*t_0^20*t_1^15+43*t_0^19*t_1^16+45*t_0^18*t_1^17+45*t_0^17*t_1^18+43*t_0^16*t_1^19+39*t_0^15*t_1^20+34*t_0^14*t_1^21+28*t_0^13*t_1^22+22*t_0^12*t_1^23+16*t_0^11*t_1^24+11*t_0^10*t_1^25+7*t_0^9*t_1^26+4*t_0^8*t_1^27+2*t_0^7*t_1^28+t_0^6*t_1^29, (4,0) => t_0^34*t_1^10+2*t_0^33*t_1^11+4*t_0^32*t_1^12+7*t_0^31*t_1^13+12*t_0^30*t_1^14+17*t_0^29*t_1^15+24*t_0^28*t_1^16+31*t_0^27*t_1^17+39*t_0^26*t_1^18+45*t_0^25*t_1^19+51*t_0^24*t_1^20+54*t_0^23*t_1^21+56*t_0^22*t_1^22+54*t_0^21*t_1^23+51*t_0^20*t_1^24+45*t_0^19*t_1^25+39*t_0^18*t_1^26+31*t_0^17*t_1^27+24*t_0^16*t_1^28+17*t_0^15*t_1^29+12*t_0^14*t_1^30+7*t_0^13*t_1^31+4*t_0^12*t_1^32+2*t_0^11*t_1^33+t_0^10*t_1^34, (3,1) => 0, (5,0) => t_0^38*t_1^15+2*t_0^37*t_1^16+4*t_0^36*t_1^17+7*t_0^35*t_1^18+11*t_0^34*t_1^19+16*t_0^33*t_1^20+22*t_0^32*t_1^21+28*t_0^31*t_1^22+34*t_0^30*t_1^23+39*t_0^29*t_1^24+43*t_0^28*t_1^25+45*t_0^27*t_1^26+45*t_0^26*t_1^27+43*t_0^25*t_1^28+39*t_0^24*t_1^29+34*t_0^23*t_1^30+28*t_0^22*t_1^31+22*t_0^21*t_1^32+16*t_0^20*t_1^33+11*t_0^19*t_1^34+7*t_0^18*t_1^35+4*t_0^17*t_1^36+2*t_0^16*t_1^37+t_0^15*t_1^38, (4,1) => 0, (5,1) => 0, (6,0) => t_0^41*t_1^21+2*t_0^40*t_1^22+4*t_0^39*t_1^23+6*t_0^38*t_1^24+9*t_0^37*t_1^25+12*t_0^36*t_1^26+16*t_0^35*t_1^27+19*t_0^34*t_1^28+22*t_0^33*t_1^29+23*t_0^32*t_1^30+24*t_0^31*t_1^31+23*t_0^30*t_1^32+22*t_0^29*t_1^33+19*t_0^28*t_1^34+16*t_0^27*t_1^35+12*t_0^26*t_1^36+9*t_0^25*t_1^37+6*t_0^24*t_1^38+4*t_0^23*t_1^39+2*t_0^22*t_1^40+t_0^21*t_1^41};
--sb represents the betti numbers as sums of Schur functors
sb98 = new HashTable from {(7,0) => {({43,28},1)}, (6,1) => {}, (8,0) => {({44,36},1)}, (7,1) => {}, (8,1) => {}, (0,0) => {({8,0},1)}, (0,1) => {}, (1,0) => {({16,1},1)}, (2,0) => {({23,3},1)}, (1,1) => {}, (2,1) => {}, (3,0) => {({29,6},1)}, (3,1) => {}, (4,0) => {({34,10},1)}, (4,1) => {}, (5,0) => {({38,15},1)}, (6,0) => {({41,21},1)}, (5,1) => {}};
--dw encodes the dominant weights in each entry
dw98 = new HashTable from {(7,0) => {{43,28}}, (6,1) => {}, (8,0) => {{44,36}}, (7,1) => {}, (8,1) => {}, (0,0) => {{8,0}}, (0,1) => {}, (1,0) => {{16,1}}, (2,0) => {{23,3}}, (1,1) => {}, (2,1) => {}, (3,0) => {{29,6}}, (3,1) => {}, (4,0) => {{34,10}}, (4,1) => {}, (5,0) => {{38,15}}, (6,0) => {{41,21}}, (5,1) => {}};
--lw encodes the lex leading weight in each entry
lw98 = new HashTable from {(7,0) => {43,28}, (6,1) => {}, (8,0) => {44,36}, (7,1) => {}, (8,1) => {}, (0,0) => {8,0}, (0,1) => {}, (1,0) => {16,1}, (2,0) => {23,3}, (1,1) => {}, (2,1) => {}, (3,0) => {29,6}, (3,1) => {}, (4,0) => {34,10}, (4,1) => {}, (5,0) => {38,15}, (6,0) => {41,21}, (5,1) => {}};
--nr encodes the number of disctinct reprsentations in each entry
nr98 = new HashTable from {(7,0) => 1, (6,1) => 0, (8,0) => 1, (7,1) => 0, (8,1) => 0, (0,0) => 1, (0,1) => 0, (1,0) => 1, (2,0) => 1, (1,1) => 0, (2,1) => 0, (3,0) => 1, (3,1) => 0, (4,0) => 1, (4,1) => 0, (5,0) => 1, (6,0) => 1, (5,1) => 0};
--nrm encodes the number of representations with multiplicity in each entry
nrm98 = new HashTable from {(7,0) => 1, (6,1) => 0, (8,0) => 1, (7,1) => 0, (8,1) => 0, (0,0) => 1, (0,1) => 0, (1,0) => 1, (2,0) => 1, (1,1) => 0, (2,1) => 0, (3,0) => 1, (3,1) => 0, (4,0) => 1, (4,1) => 0, (5,0) => 1, (6,0) => 1, (5,1) => 0};
--er encodes the errors in the computed multigraded Hilbert series via our Schur method in each entry
er98 = new HashTable from {(7,0) => 0, (6,1) => 0, (8,0) => 0, (7,1) => 0, (8,1) => 0, (0,0) => 0, (0,1) => 0, (1,0) => 0, (2,0) => 0, (1,1) => 0, (2,1) => 0, (3,0) => 0, (3,1) => 0, (4,0) => 0, (4,1) => 0, (5,0) => 0, (6,0) => 0, (5,1) => 0};
--bs encodes the Boij-Soederberg coefficients each entry
bs98 = {362880/1};
end;