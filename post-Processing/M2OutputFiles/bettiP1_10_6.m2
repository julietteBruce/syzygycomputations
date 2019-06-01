--This file computes Betti tables for P^2 for d = 10 and b = 6
A := QQ[t_0,t_1,t_2];
--tb stands for Total Betti numbers
tb106 = new HashTable from {(7,0) => 0, (6,1) => 0, (8,0) => 0, (7,1) => 45, (9,0) => 0, (8,1) => 20, (9,1) => 3, (0,0) => 7, (0,1) => 0, (1,0) => 60, (2,0) => 225, (1,1) => 0, (2,1) => 0, (3,0) => 480, (3,1) => 0, (4,0) => 630, (4,1) => 0, (5,0) => 504, (6,0) => 210, (5,1) => 0};
--mb stands for Multigraded Betti numbers
mb106 = new HashTable from {(7,0) => 0, (6,1) => 0, (7,1) => t_0^51*t_1^35+t_0^50*t_1^36+2*t_0^49*t_1^37+2*t_0^48*t_1^38+3*t_0^47*t_1^39+3*t_0^46*t_1^40+4*t_0^45*t_1^41+4*t_0^44*t_1^42+5*t_0^43*t_1^43+4*t_0^42*t_1^44+4*t_0^41*t_1^45+3*t_0^40*t_1^46+3*t_0^39*t_1^47+2*t_0^38*t_1^48+2*t_0^37*t_1^49+t_0^36*t_1^50+t_0^35*t_1^51, (8,0) => 0, (8,1) => t_0^53*t_1^43+2*t_0^52*t_1^44+2*t_0^51*t_1^45+2*t_0^50*t_1^46+2*t_0^49*t_1^47+2*t_0^48*t_1^48+2*t_0^47*t_1^49+2*t_0^46*t_1^50+2*t_0^45*t_1^51+2*t_0^44*t_1^52+t_0^43*t_1^53, (9,0) => 0, (9,1) => t_0^54*t_1^52+t_0^53*t_1^53+t_0^52*t_1^54, (0,0) => t_0^6+t_0^5*t_1+t_0^4*t_1^2+t_0^3*t_1^3+t_0^2*t_1^4+t_0*t_1^5+t_1^6, (0,1) => 0, (1,0) => t_0^15*t_1+2*t_0^14*t_1^2+3*t_0^13*t_1^3+4*t_0^12*t_1^4+5*t_0^11*t_1^5+6*t_0^10*t_1^6+6*t_0^9*t_1^7+6*t_0^8*t_1^8+6*t_0^7*t_1^9+6*t_0^6*t_1^10+5*t_0^5*t_1^11+4*t_0^4*t_1^12+3*t_0^3*t_1^13+2*t_0^2*t_1^14+t_0*t_1^15, (2,0) => t_0^23*t_1^3+2*t_0^22*t_1^4+4*t_0^21*t_1^5+6*t_0^20*t_1^6+9*t_0^19*t_1^7+11*t_0^18*t_1^8+14*t_0^17*t_1^9+16*t_0^16*t_1^10+19*t_0^15*t_1^11+20*t_0^14*t_1^12+21*t_0^13*t_1^13+20*t_0^12*t_1^14+19*t_0^11*t_1^15+16*t_0^10*t_1^16+14*t_0^9*t_1^17+11*t_0^8*t_1^18+9*t_0^7*t_1^19+6*t_0^6*t_1^20+4*t_0^5*t_1^21+2*t_0^4*t_1^22+t_0^3*t_1^23, (1,1) => 0, (3,0) => t_0^30*t_1^6+2*t_0^29*t_1^7+4*t_0^28*t_1^8+7*t_0^27*t_1^9+10*t_0^26*t_1^10+14*t_0^25*t_1^11+19*t_0^24*t_1^12+24*t_0^23*t_1^13+29*t_0^22*t_1^14+34*t_0^21*t_1^15+37*t_0^20*t_1^16+39*t_0^19*t_1^17+40*t_0^18*t_1^18+39*t_0^17*t_1^19+37*t_0^16*t_1^20+34*t_0^15*t_1^21+29*t_0^14*t_1^22+24*t_0^13*t_1^23+19*t_0^12*t_1^24+14*t_0^11*t_1^25+10*t_0^10*t_1^26+7*t_0^9*t_1^27+4*t_0^8*t_1^28+2*t_0^7*t_1^29+t_0^6*t_1^30, (2,1) => 0, (4,0) => t_0^36*t_1^10+2*t_0^35*t_1^11+4*t_0^34*t_1^12+6*t_0^33*t_1^13+10*t_0^32*t_1^14+14*t_0^31*t_1^15+20*t_0^30*t_1^16+25*t_0^29*t_1^17+32*t_0^28*t_1^18+37*t_0^27*t_1^19+43*t_0^26*t_1^20+46*t_0^25*t_1^21+50*t_0^24*t_1^22+50*t_0^23*t_1^23+50*t_0^22*t_1^24+46*t_0^21*t_1^25+43*t_0^20*t_1^26+37*t_0^19*t_1^27+32*t_0^18*t_1^28+25*t_0^17*t_1^29+20*t_0^16*t_1^30+14*t_0^15*t_1^31+10*t_0^14*t_1^32+6*t_0^13*t_1^33+4*t_0^12*t_1^34+2*t_0^11*t_1^35+t_0^10*t_1^36, (3,1) => 0, (5,0) => t_0^41*t_1^15+2*t_0^40*t_1^16+3*t_0^39*t_1^17+5*t_0^38*t_1^18+8*t_0^37*t_1^19+12*t_0^36*t_1^20+16*t_0^35*t_1^21+20*t_0^34*t_1^22+25*t_0^33*t_1^23+30*t_0^32*t_1^24+34*t_0^31*t_1^25+37*t_0^30*t_1^26+39*t_0^29*t_1^27+40*t_0^28*t_1^28+39*t_0^27*t_1^29+37*t_0^26*t_1^30+34*t_0^25*t_1^31+30*t_0^24*t_1^32+25*t_0^23*t_1^33+20*t_0^22*t_1^34+16*t_0^21*t_1^35+12*t_0^20*t_1^36+8*t_0^19*t_1^37+5*t_0^18*t_1^38+3*t_0^17*t_1^39+2*t_0^16*t_1^40+t_0^15*t_1^41, (4,1) => 0, (5,1) => 0, (6,0) => t_0^45*t_1^21+t_0^44*t_1^22+2*t_0^43*t_1^23+3*t_0^42*t_1^24+5*t_0^41*t_1^25+6*t_0^40*t_1^26+9*t_0^39*t_1^27+10*t_0^38*t_1^28+13*t_0^37*t_1^29+14*t_0^36*t_1^30+16*t_0^35*t_1^31+16*t_0^34*t_1^32+18*t_0^33*t_1^33+16*t_0^32*t_1^34+16*t_0^31*t_1^35+14*t_0^30*t_1^36+13*t_0^29*t_1^37+10*t_0^28*t_1^38+9*t_0^27*t_1^39+6*t_0^26*t_1^40+5*t_0^25*t_1^41+3*t_0^24*t_1^42+2*t_0^23*t_1^43+t_0^22*t_1^44+t_0^21*t_1^45};
--sb represents the betti numbers as sums of Schur functors
sb106 = new HashTable from {(7,0) => {}, (6,1) => {}, (8,0) => {}, (7,1) => {({51,35},1)}, (9,0) => {}, (8,1) => {({53,43},1)}, (9,1) => {({54,52},1)}, (0,0) => {({6,0},1)}, (0,1) => {}, (1,0) => {({15,1},1)}, (2,0) => {({23,3},1)}, (1,1) => {}, (2,1) => {}, (3,0) => {({30,6},1)}, (3,1) => {}, (4,0) => {({36,10},1)}, (4,1) => {}, (5,0) => {({41,15},1)}, (6,0) => {({45,21},1)}, (5,1) => {}};
--dw encodes the dominant weights in each entry
dw106 = new HashTable from {(7,0) => {}, (6,1) => {}, (8,0) => {}, (7,1) => {{51,35}}, (9,0) => {}, (8,1) => {{53,43}}, (9,1) => {{54,52}}, (0,0) => {{6,0}}, (0,1) => {}, (1,0) => {{15,1}}, (2,0) => {{23,3}}, (1,1) => {}, (2,1) => {}, (3,0) => {{30,6}}, (3,1) => {}, (4,0) => {{36,10}}, (4,1) => {}, (5,0) => {{41,15}}, (6,0) => {{45,21}}, (5,1) => {}};
--lw encodes the lex leading weight in each entry
lw106 = new HashTable from {(7,0) => {}, (6,1) => {}, (8,0) => {}, (7,1) => {51,35}, (9,0) => {}, (8,1) => {53,43}, (9,1) => {54,52}, (0,0) => {6,0}, (0,1) => {}, (1,0) => {15,1}, (2,0) => {23,3}, (1,1) => {}, (2,1) => {}, (3,0) => {30,6}, (3,1) => {}, (4,0) => {36,10}, (4,1) => {}, (5,0) => {41,15}, (6,0) => {45,21}, (5,1) => {}};
--nr encodes the number of disctinct reprsentations in each entry
nr106 = new HashTable from {(7,0) => 0, (6,1) => 0, (8,0) => 0, (7,1) => 1, (9,0) => 0, (8,1) => 1, (9,1) => 1, (0,0) => 1, (0,1) => 0, (1,0) => 1, (2,0) => 1, (1,1) => 0, (2,1) => 0, (3,0) => 1, (3,1) => 0, (4,0) => 1, (4,1) => 0, (5,0) => 1, (6,0) => 1, (5,1) => 0};
--nrm encodes the number of representations with multiplicity in each entry
nrm106 = new HashTable from {(7,0) => 0, (6,1) => 0, (8,0) => 0, (7,1) => 1, (9,0) => 0, (8,1) => 1, (9,1) => 1, (0,0) => 1, (0,1) => 0, (1,0) => 1, (2,0) => 1, (1,1) => 0, (2,1) => 0, (3,0) => 1, (3,1) => 0, (4,0) => 1, (4,1) => 0, (5,0) => 1, (6,0) => 1, (5,1) => 0};
--er encodes the errors in the computed multigraded Hilbert series via our Schur method in each entry
er106 = new HashTable from {(7,0) => 0, (6,1) => 0, (8,0) => 0, (7,1) => 0, (9,0) => 0, (8,1) => 0, (9,1) => 0, (0,0) => 0, (0,1) => 0, (1,0) => 0, (2,0) => 0, (1,1) => 0, (2,1) => 0, (3,0) => 0, (3,1) => 0, (4,0) => 0, (4,1) => 0, (5,0) => 0, (6,0) => 0, (5,1) => 0};
--bs encodes the Boij-Soederberg coefficients each entry
bs106 = {3628800/1};
end;