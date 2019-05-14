--This file computes Betti tables for P^2 for d = 8 and b = 4
A := QQ[t_0,t_1,t_2];
--tb stands for Total Betti numbers
tb84 = new HashTable from {(7,0) => 0, (6,1) => 16, (7,1) => 3, (0,0) => 5, (0,1) => 0, (1,0) => 32, (2,0) => 84, (1,1) => 0, (2,1) => 0, (3,0) => 112, (3,1) => 0, (4,0) => 70, (4,1) => 0, (5,0) => 0, (6,0) => 0, (5,1) => 28};
--mb stands for Multigraded Betti numbers
mb84 = new HashTable from {(7,0) => 0, (6,1) => t_0^34*t_1^26+2*t_0^33*t_1^27+2*t_0^32*t_1^28+2*t_0^31*t_1^29+2*t_0^30*t_1^30+2*t_0^29*t_1^31+2*t_0^28*t_1^32+2*t_0^27*t_1^33+t_0^26*t_1^34, (7,1) => t_0^35*t_1^33+t_0^34*t_1^34+t_0^33*t_1^35, (0,0) => t_0^4+t_0^3*t_1+t_0^2*t_1^2+t_0*t_1^3+t_1^4, (1,0) => t_0^11*t_1+2*t_0^10*t_1^2+3*t_0^9*t_1^3+4*t_0^8*t_1^4+4*t_0^7*t_1^5+4*t_0^6*t_1^6+4*t_0^5*t_1^7+4*t_0^4*t_1^8+3*t_0^3*t_1^9+2*t_0^2*t_1^10+t_0*t_1^11, (0,1) => 0, (2,0) => t_0^17*t_1^3+2*t_0^16*t_1^4+4*t_0^15*t_1^5+5*t_0^14*t_1^6+7*t_0^13*t_1^7+8*t_0^12*t_1^8+10*t_0^11*t_1^9+10*t_0^10*t_1^10+10*t_0^9*t_1^11+8*t_0^8*t_1^12+7*t_0^7*t_1^13+5*t_0^6*t_1^14+4*t_0^5*t_1^15+2*t_0^4*t_1^16+t_0^3*t_1^17, (1,1) => 0, (2,1) => 0, (3,0) => t_0^22*t_1^6+2*t_0^21*t_1^7+3*t_0^20*t_1^8+5*t_0^19*t_1^9+7*t_0^18*t_1^10+9*t_0^17*t_1^11+11*t_0^16*t_1^12+12*t_0^15*t_1^13+12*t_0^14*t_1^14+12*t_0^13*t_1^15+11*t_0^12*t_1^16+9*t_0^11*t_1^17+7*t_0^10*t_1^18+5*t_0^9*t_1^19+3*t_0^8*t_1^20+2*t_0^7*t_1^21+t_0^6*t_1^22, (3,1) => 0, (4,0) => t_0^26*t_1^10+t_0^25*t_1^11+2*t_0^24*t_1^12+3*t_0^23*t_1^13+5*t_0^22*t_1^14+5*t_0^21*t_1^15+7*t_0^20*t_1^16+7*t_0^19*t_1^17+8*t_0^18*t_1^18+7*t_0^17*t_1^19+7*t_0^16*t_1^20+5*t_0^15*t_1^21+5*t_0^14*t_1^22+3*t_0^13*t_1^23+2*t_0^12*t_1^24+t_0^11*t_1^25+t_0^10*t_1^26, (5,0) => 0, (4,1) => 0, (5,1) => t_0^32*t_1^20+t_0^31*t_1^21+2*t_0^30*t_1^22+2*t_0^29*t_1^23+3*t_0^28*t_1^24+3*t_0^27*t_1^25+4*t_0^26*t_1^26+3*t_0^25*t_1^27+3*t_0^24*t_1^28+2*t_0^23*t_1^29+2*t_0^22*t_1^30+t_0^21*t_1^31+t_0^20*t_1^32, (6,0) => 0};
--sb represents the betti numbers as sums of Schur functors
sb84 = new HashTable from {(7,0) => {}, (6,1) => {({34,26},1)}, (7,1) => {({35,33},1)}, (0,0) => {({4,0},1)}, (0,1) => {}, (1,0) => {({11,1},1)}, (2,0) => {({17,3},1)}, (1,1) => {}, (2,1) => {}, (3,0) => {({22,6},1)}, (3,1) => {}, (4,0) => {({26,10},1)}, (4,1) => {}, (5,0) => {}, (6,0) => {}, (5,1) => {({32,20},1)}};
--dw encodes the dominant weights in each entry
dw84 = new HashTable from {(7,0) => {}, (6,1) => {{34,26}}, (7,1) => {{35,33}}, (0,0) => {{4,0}}, (0,1) => {}, (1,0) => {{11,1}}, (2,0) => {{17,3}}, (1,1) => {}, (2,1) => {}, (3,0) => {{22,6}}, (3,1) => {}, (4,0) => {{26,10}}, (4,1) => {}, (5,0) => {}, (6,0) => {}, (5,1) => {{32,20}}};
--lw encodes the lex leading weight in each entry
lw84 = new HashTable from {(7,0) => {}, (6,1) => {34,26}, (7,1) => {35,33}, (0,0) => {4,0}, (0,1) => {}, (1,0) => {11,1}, (2,0) => {17,3}, (1,1) => {}, (2,1) => {}, (3,0) => {22,6}, (3,1) => {}, (4,0) => {26,10}, (4,1) => {}, (5,0) => {}, (6,0) => {}, (5,1) => {32,20}};
--nr encodes the number of disctinct reprsentations in each entry
nr84 = new HashTable from {(7,0) => 0, (6,1) => 1, (7,1) => 1, (0,0) => 1, (0,1) => 0, (1,0) => 1, (2,0) => 1, (1,1) => 0, (2,1) => 0, (3,0) => 1, (3,1) => 0, (4,0) => 1, (4,1) => 0, (5,0) => 0, (6,0) => 0, (5,1) => 1};
--nrm encodes the number of representations with multiplicity in each entry
nrm84 = new HashTable from {(7,0) => 0, (6,1) => 1, (7,1) => 1, (0,0) => 1, (0,1) => 0, (1,0) => 1, (2,0) => 1, (1,1) => 0, (2,1) => 0, (3,0) => 1, (3,1) => 0, (4,0) => 1, (4,1) => 0, (5,0) => 0, (6,0) => 0, (5,1) => 1};
--er encodes the errors in the computed multigraded Hilbert series via our Schur method in each entry
er84 = new HashTable from {(7,0) => 0, (6,1) => 0, (7,1) => 0, (0,0) => 0, (0,1) => 0, (1,0) => 0, (2,0) => 0, (1,1) => 0, (2,1) => 0, (3,0) => 0, (3,1) => 0, (4,0) => 0, (4,1) => 0, (5,0) => 0, (6,0) => 0, (5,1) => 0};
--bs encodes the Boij-Soederberg coefficients each entry
bs84 = {40320/1};
end;