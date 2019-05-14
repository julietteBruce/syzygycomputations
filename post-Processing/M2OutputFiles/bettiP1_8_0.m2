--This file computes Betti tables for P^2 for d = 8 and b = 0
A := QQ[t_0,t_1,t_2];
--tb stands for Total Betti numbers
tb80 = new HashTable from {(7,0) => 0, (6,1) => 48, (7,1) => 7, (0,0) => 1, (0,1) => 0, (1,0) => 0, (2,0) => 0, (1,1) => 28, (2,1) => 112, (3,0) => 0, (3,1) => 210, (4,0) => 0, (4,1) => 224, (5,0) => 0, (6,0) => 0, (5,1) => 140};
--mb stands for Multigraded Betti numbers
mb80 = new HashTable from {(7,0) => 0, (6,1) => t_0^34*t_1^22+2*t_0^33*t_1^23+3*t_0^32*t_1^24+4*t_0^31*t_1^25+5*t_0^30*t_1^26+6*t_0^29*t_1^27+6*t_0^28*t_1^28+6*t_0^27*t_1^29+5*t_0^26*t_1^30+4*t_0^25*t_1^31+3*t_0^24*t_1^32+2*t_0^23*t_1^33+t_0^22*t_1^34, (7,1) => t_0^35*t_1^29+t_0^34*t_1^30+t_0^33*t_1^31+t_0^32*t_1^32+t_0^31*t_1^33+t_0^30*t_1^34+t_0^29*t_1^35, (0,0) => 1, (1,0) => 0, (0,1) => 0, (2,0) => 0, (1,1) => t_0^14*t_1^2+t_0^13*t_1^3+2*t_0^12*t_1^4+2*t_0^11*t_1^5+3*t_0^10*t_1^6+3*t_0^9*t_1^7+4*t_0^8*t_1^8+3*t_0^7*t_1^9+3*t_0^6*t_1^10+2*t_0^5*t_1^11+2*t_0^4*t_1^12+t_0^3*t_1^13+t_0^2*t_1^14, (2,1) => t_0^20*t_1^4+2*t_0^19*t_1^5+3*t_0^18*t_1^6+5*t_0^17*t_1^7+7*t_0^16*t_1^8+9*t_0^15*t_1^9+11*t_0^14*t_1^10+12*t_0^13*t_1^11+12*t_0^12*t_1^12+12*t_0^11*t_1^13+11*t_0^10*t_1^14+9*t_0^9*t_1^15+7*t_0^8*t_1^16+5*t_0^7*t_1^17+3*t_0^6*t_1^18+2*t_0^5*t_1^19+t_0^4*t_1^20, (3,0) => 0, (3,1) => t_0^25*t_1^7+2*t_0^24*t_1^8+4*t_0^23*t_1^9+6*t_0^22*t_1^10+10*t_0^21*t_1^11+13*t_0^20*t_1^12+17*t_0^19*t_1^13+19*t_0^18*t_1^14+22*t_0^17*t_1^15+22*t_0^16*t_1^16+22*t_0^15*t_1^17+19*t_0^14*t_1^18+17*t_0^13*t_1^19+13*t_0^12*t_1^20+10*t_0^11*t_1^21+6*t_0^10*t_1^22+4*t_0^9*t_1^23+2*t_0^8*t_1^24+t_0^7*t_1^25, (4,0) => 0, (5,0) => 0, (4,1) => t_0^29*t_1^11+2*t_0^28*t_1^12+4*t_0^27*t_1^13+7*t_0^26*t_1^14+10*t_0^25*t_1^15+14*t_0^24*t_1^16+18*t_0^23*t_1^17+21*t_0^22*t_1^18+23*t_0^21*t_1^19+24*t_0^20*t_1^20+23*t_0^19*t_1^21+21*t_0^18*t_1^22+18*t_0^17*t_1^23+14*t_0^16*t_1^24+10*t_0^15*t_1^25+7*t_0^14*t_1^26+4*t_0^13*t_1^27+2*t_0^12*t_1^28+t_0^11*t_1^29, (5,1) => t_0^32*t_1^16+2*t_0^31*t_1^17+4*t_0^30*t_1^18+6*t_0^29*t_1^19+9*t_0^28*t_1^20+11*t_0^27*t_1^21+14*t_0^26*t_1^22+15*t_0^25*t_1^23+16*t_0^24*t_1^24+15*t_0^23*t_1^25+14*t_0^22*t_1^26+11*t_0^21*t_1^27+9*t_0^20*t_1^28+6*t_0^19*t_1^29+4*t_0^18*t_1^30+2*t_0^17*t_1^31+t_0^16*t_1^32, (6,0) => 0};
--sb represents the betti numbers as sums of Schur functors
sb80 = new HashTable from {(7,0) => {}, (6,1) => {({34,22},1)}, (7,1) => {({35,29},1)}, (0,0) => {({0,0},1)}, (0,1) => {}, (1,0) => {}, (2,0) => {}, (1,1) => {({14,2},1)}, (2,1) => {({20,4},1)}, (3,0) => {}, (3,1) => {({25,7},1)}, (4,0) => {}, (4,1) => {({29,11},1)}, (5,0) => {}, (6,0) => {}, (5,1) => {({32,16},1)}};
--dw encodes the dominant weights in each entry
dw80 = new HashTable from {(7,0) => {}, (6,1) => {{34,22}}, (7,1) => {{35,29}}, (0,0) => {{0,0}}, (0,1) => {}, (1,0) => {}, (2,0) => {}, (1,1) => {{14,2}}, (2,1) => {{20,4}}, (3,0) => {}, (3,1) => {{25,7}}, (4,0) => {}, (4,1) => {{29,11}}, (5,0) => {}, (6,0) => {}, (5,1) => {{32,16}}};
--lw encodes the lex leading weight in each entry
lw80 = new HashTable from {(7,0) => {}, (6,1) => {34,22}, (7,1) => {35,29}, (0,0) => {0,0}, (0,1) => {}, (1,0) => {}, (2,0) => {}, (1,1) => {14,2}, (2,1) => {20,4}, (3,0) => {}, (3,1) => {25,7}, (4,0) => {}, (4,1) => {29,11}, (5,0) => {}, (6,0) => {}, (5,1) => {32,16}};
--nr encodes the number of disctinct reprsentations in each entry
nr80 = new HashTable from {(7,0) => 0, (6,1) => 1, (7,1) => 1, (0,0) => 1, (0,1) => 0, (1,0) => 0, (2,0) => 0, (1,1) => 1, (2,1) => 1, (3,0) => 0, (3,1) => 1, (4,0) => 0, (4,1) => 1, (5,0) => 0, (6,0) => 0, (5,1) => 1};
--nrm encodes the number of representations with multiplicity in each entry
nrm80 = new HashTable from {(7,0) => 0, (6,1) => 1, (7,1) => 1, (0,0) => 1, (0,1) => 0, (1,0) => 0, (2,0) => 0, (1,1) => 1, (2,1) => 1, (3,0) => 0, (3,1) => 1, (4,0) => 0, (4,1) => 1, (5,0) => 0, (6,0) => 0, (5,1) => 1};
--er encodes the errors in the computed multigraded Hilbert series via our Schur method in each entry
er80 = new HashTable from {(7,0) => 0, (6,1) => 0, (7,1) => 0, (0,0) => 0, (0,1) => 0, (1,0) => 0, (2,0) => 0, (1,1) => 0, (2,1) => 0, (3,0) => 0, (3,1) => 0, (4,0) => 0, (4,1) => 0, (5,0) => 0, (6,0) => 0, (5,1) => 0};
--bs encodes the Boij-Soederberg coefficients each entry
bs80 = {40320/1};
end;