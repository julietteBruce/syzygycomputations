--This file computes Betti tables for P^2 for d = 10 and b = 4
A := QQ[t_0,t_1,t_2];
--tb stands for Total Betti numbers
tb104 = new HashTable from {(7,0) => 0, (6,1) => 240, (8,0) => 0, (7,1) => 135, (9,0) => 0, (8,1) => 40, (9,1) => 5, (0,0) => 5, (0,1) => 0, (1,0) => 40, (2,0) => 135, (1,1) => 0, (2,1) => 0, (3,0) => 240, (3,1) => 0, (4,0) => 210, (4,1) => 0, (5,0) => 0, (6,0) => 0, (5,1) => 210};
--mb stands for Multigraded Betti numbers
mb104 = new HashTable from {(7,0) => 0, (6,1) => t_0^48*t_1^26+2*t_0^47*t_1^27+3*t_0^46*t_1^28+5*t_0^45*t_1^29+7*t_0^44*t_1^30+9*t_0^43*t_1^31+12*t_0^42*t_1^32+15*t_0^41*t_1^33+17*t_0^40*t_1^34+19*t_0^39*t_1^35+20*t_0^38*t_1^36+20*t_0^37*t_1^37+20*t_0^36*t_1^38+19*t_0^35*t_1^39+17*t_0^34*t_1^40+15*t_0^33*t_1^41+12*t_0^32*t_1^42+9*t_0^31*t_1^43+7*t_0^30*t_1^44+5*t_0^29*t_1^45+3*t_0^28*t_1^46+2*t_0^27*t_1^47+t_0^26*t_1^48, (7,1) => t_0^51*t_1^33+2*t_0^50*t_1^34+4*t_0^49*t_1^35+5*t_0^48*t_1^36+7*t_0^47*t_1^37+8*t_0^46*t_1^38+10*t_0^45*t_1^39+11*t_0^44*t_1^40+13*t_0^43*t_1^41+13*t_0^42*t_1^42+13*t_0^41*t_1^43+11*t_0^40*t_1^44+10*t_0^39*t_1^45+8*t_0^38*t_1^46+7*t_0^37*t_1^47+5*t_0^36*t_1^48+4*t_0^35*t_1^49+2*t_0^34*t_1^50+t_0^33*t_1^51, (8,0) => 0, (8,1) => t_0^53*t_1^41+2*t_0^52*t_1^42+3*t_0^51*t_1^43+4*t_0^50*t_1^44+4*t_0^49*t_1^45+4*t_0^48*t_1^46+4*t_0^47*t_1^47+4*t_0^46*t_1^48+4*t_0^45*t_1^49+4*t_0^44*t_1^50+3*t_0^43*t_1^51+2*t_0^42*t_1^52+t_0^41*t_1^53, (9,0) => 0, (9,1) => t_0^54*t_1^50+t_0^53*t_1^51+t_0^52*t_1^52+t_0^51*t_1^53+t_0^50*t_1^54, (0,0) => t_0^4+t_0^3*t_1+t_0^2*t_1^2+t_0*t_1^3+t_1^4, (0,1) => 0, (1,0) => t_0^13*t_1+2*t_0^12*t_1^2+3*t_0^11*t_1^3+4*t_0^10*t_1^4+4*t_0^9*t_1^5+4*t_0^8*t_1^6+4*t_0^7*t_1^7+4*t_0^6*t_1^8+4*t_0^5*t_1^9+4*t_0^4*t_1^10+3*t_0^3*t_1^11+2*t_0^2*t_1^12+t_0*t_1^13, (2,0) => t_0^21*t_1^3+2*t_0^20*t_1^4+4*t_0^19*t_1^5+5*t_0^18*t_1^6+7*t_0^17*t_1^7+8*t_0^16*t_1^8+10*t_0^15*t_1^9+11*t_0^14*t_1^10+13*t_0^13*t_1^11+13*t_0^12*t_1^12+13*t_0^11*t_1^13+11*t_0^10*t_1^14+10*t_0^9*t_1^15+8*t_0^8*t_1^16+7*t_0^7*t_1^17+5*t_0^6*t_1^18+4*t_0^5*t_1^19+2*t_0^4*t_1^20+t_0^3*t_1^21, (1,1) => 0, (3,0) => t_0^28*t_1^6+2*t_0^27*t_1^7+3*t_0^26*t_1^8+5*t_0^25*t_1^9+7*t_0^24*t_1^10+9*t_0^23*t_1^11+12*t_0^22*t_1^12+15*t_0^21*t_1^13+17*t_0^20*t_1^14+19*t_0^19*t_1^15+20*t_0^18*t_1^16+20*t_0^17*t_1^17+20*t_0^16*t_1^18+19*t_0^15*t_1^19+17*t_0^14*t_1^20+15*t_0^13*t_1^21+12*t_0^12*t_1^22+9*t_0^11*t_1^23+7*t_0^10*t_1^24+5*t_0^9*t_1^25+3*t_0^8*t_1^26+2*t_0^7*t_1^27+t_0^6*t_1^28, (2,1) => 0, (4,0) => t_0^34*t_1^10+t_0^33*t_1^11+2*t_0^32*t_1^12+3*t_0^31*t_1^13+5*t_0^30*t_1^14+6*t_0^29*t_1^15+9*t_0^28*t_1^16+10*t_0^27*t_1^17+13*t_0^26*t_1^18+14*t_0^25*t_1^19+16*t_0^24*t_1^20+16*t_0^23*t_1^21+18*t_0^22*t_1^22+16*t_0^21*t_1^23+16*t_0^20*t_1^24+14*t_0^19*t_1^25+13*t_0^18*t_1^26+10*t_0^17*t_1^27+9*t_0^16*t_1^28+6*t_0^15*t_1^29+5*t_0^14*t_1^30+3*t_0^13*t_1^31+2*t_0^12*t_1^32+t_0^11*t_1^33+t_0^10*t_1^34, (3,1) => 0, (5,0) => 0, (4,1) => 0, (5,1) => t_0^44*t_1^20+t_0^43*t_1^21+2*t_0^42*t_1^22+3*t_0^41*t_1^23+5*t_0^40*t_1^24+6*t_0^39*t_1^25+9*t_0^38*t_1^26+10*t_0^37*t_1^27+13*t_0^36*t_1^28+14*t_0^35*t_1^29+16*t_0^34*t_1^30+16*t_0^33*t_1^31+18*t_0^32*t_1^32+16*t_0^31*t_1^33+16*t_0^30*t_1^34+14*t_0^29*t_1^35+13*t_0^28*t_1^36+10*t_0^27*t_1^37+9*t_0^26*t_1^38+6*t_0^25*t_1^39+5*t_0^24*t_1^40+3*t_0^23*t_1^41+2*t_0^22*t_1^42+t_0^21*t_1^43+t_0^20*t_1^44, (6,0) => 0};
--sb represents the betti numbers as sums of Schur functors
sb104 = new HashTable from {(7,0) => {}, (6,1) => {({48,26},1)}, (8,0) => {}, (7,1) => {({51,33},1)}, (9,0) => {}, (8,1) => {({53,41},1)}, (9,1) => {({54,50},1)}, (0,0) => {({4,0},1)}, (0,1) => {}, (1,0) => {({13,1},1)}, (2,0) => {({21,3},1)}, (1,1) => {}, (2,1) => {}, (3,0) => {({28,6},1)}, (3,1) => {}, (4,0) => {({34,10},1)}, (4,1) => {}, (5,0) => {}, (6,0) => {}, (5,1) => {({44,20},1)}};
--dw encodes the dominant weights in each entry
dw104 = new HashTable from {(7,0) => {}, (6,1) => {{48,26}}, (8,0) => {}, (7,1) => {{51,33}}, (9,0) => {}, (8,1) => {{53,41}}, (9,1) => {{54,50}}, (0,0) => {{4,0}}, (0,1) => {}, (1,0) => {{13,1}}, (2,0) => {{21,3}}, (1,1) => {}, (2,1) => {}, (3,0) => {{28,6}}, (3,1) => {}, (4,0) => {{34,10}}, (4,1) => {}, (5,0) => {}, (6,0) => {}, (5,1) => {{44,20}}};
--lw encodes the lex leading weight in each entry
lw104 = new HashTable from {(7,0) => {}, (6,1) => {48,26}, (8,0) => {}, (7,1) => {51,33}, (9,0) => {}, (8,1) => {53,41}, (9,1) => {54,50}, (0,0) => {4,0}, (0,1) => {}, (1,0) => {13,1}, (2,0) => {21,3}, (1,1) => {}, (2,1) => {}, (3,0) => {28,6}, (3,1) => {}, (4,0) => {34,10}, (4,1) => {}, (5,0) => {}, (6,0) => {}, (5,1) => {44,20}};
--nr encodes the number of disctinct reprsentations in each entry
nr104 = new HashTable from {(7,0) => 0, (6,1) => 1, (8,0) => 0, (7,1) => 1, (9,0) => 0, (8,1) => 1, (9,1) => 1, (0,0) => 1, (0,1) => 0, (1,0) => 1, (2,0) => 1, (1,1) => 0, (2,1) => 0, (3,0) => 1, (3,1) => 0, (4,0) => 1, (4,1) => 0, (5,0) => 0, (6,0) => 0, (5,1) => 1};
--nrm encodes the number of representations with multiplicity in each entry
nrm104 = new HashTable from {(7,0) => 0, (6,1) => 1, (8,0) => 0, (7,1) => 1, (9,0) => 0, (8,1) => 1, (9,1) => 1, (0,0) => 1, (0,1) => 0, (1,0) => 1, (2,0) => 1, (1,1) => 0, (2,1) => 0, (3,0) => 1, (3,1) => 0, (4,0) => 1, (4,1) => 0, (5,0) => 0, (6,0) => 0, (5,1) => 1};
--er encodes the errors in the computed multigraded Hilbert series via our Schur method in each entry
er104 = new HashTable from {(7,0) => 0, (6,1) => 0, (8,0) => 0, (7,1) => 0, (9,0) => 0, (8,1) => 0, (9,1) => 0, (0,0) => 0, (0,1) => 0, (1,0) => 0, (2,0) => 0, (1,1) => 0, (2,1) => 0, (3,0) => 0, (3,1) => 0, (4,0) => 0, (4,1) => 0, (5,0) => 0, (6,0) => 0, (5,1) => 0};
--bs encodes the Boij-Soederberg coefficients each entry
bs104 = {3628800/1};
end;