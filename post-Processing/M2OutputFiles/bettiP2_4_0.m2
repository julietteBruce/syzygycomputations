--This file computes Betti tables for P^2 for d = 4 and b = 0
A := QQ[t_0,t_1,t_2];
--tb stands for Total Betti numbers
tb40 = new HashTable from {(5,2) => 0, (6,1) => 7920, (7,0) => 0, (8,0) => 0, (6,2) => 0, (7,1) => 6237, (9,0) => 0, (8,1) => 3344, (7,2) => 0, (10,0) => 0, (9,1) => 1089, (8,2) => 0, (11,0) => 0, (10,1) => 120, (9,2) => 0, (12,0) => 0, (11,1) => 0, (10,2) => 55, (12,1) => 0, (11,2) => 24, (12,2) => 3, (0,0) => 1, (0,1) => 0, (1,0) => 0, (2,0) => 0, (0,2) => 0, (1,1) => 75, (3,0) => 0, (2,1) => 536, (1,2) => 0, (3,1) => 1947, (2,2) => 0, (4,0) => 0, (3,2) => 0, (4,1) => 4488, (5,0) => 0, (4,2) => 0, (5,1) => 7095, (6,0) => 0};
--mb stands for Multigraded Betti numbers
mb40 = new HashTable from {(7,0) => 0, (6,1) => t_0^16*t_1^8*t_2^4+2*t_0^16*t_1^7*t_2^5+2*t_0^16*t_1^6*t_2^6+2*t_0^16*t_1^5*t_2^7+t_0^16*t_1^4*t_2^8+t_0^15*t_1^10*t_2^3+4*t_0^15*t_1^9*t_2^4+9*t_0^15*t_1^8*t_2^5+12*t_0^15*t_1^7*t_2^6+12*t_0^15*t_1^6*t_2^7+9*t_0^15*t_1^5*t_2^8+4*t_0^15*t_1^4*t_2^9+t_0^15*t_1^3*t_2^10+2*t_0^14*t_1^11*t_2^3+9*t_0^14*t_1^10*t_2^4+21*t_0^14*t_1^9*t_2^5+33*t_0^14*t_1^8*t_2^6+38*t_0^14*t_1^7*t_2^7+33*t_0^14*t_1^6*t_2^8+21*t_0^14*t_1^5*t_2^9+9*t_0^14*t_1^4*t_2^10+2*t_0^14*t_1^3*t_2^11+3*t_0^13*t_1^12*t_2^3+14*t_0^13*t_1^11*t_2^4+36*t_0^13*t_1^10*t_2^5+62*t_0^13*t_1^9*t_2^6+82*t_0^13*t_1^8*t_2^7+82*t_0^13*t_1^7*t_2^8+62*t_0^13*t_1^6*t_2^9+36*t_0^13*t_1^5*t_2^10+14*t_0^13*t_1^4*t_2^11+3*t_0^13*t_1^3*t_2^12+3*t_0^12*t_1^13*t_2^3+16*t_0^12*t_1^12*t_2^4+46*t_0^12*t_1^11*t_2^5+89*t_0^12*t_1^10*t_2^6+130*t_0^12*t_1^9*t_2^7+148*t_0^12*t_1^8*t_2^8+130*t_0^12*t_1^7*t_2^9+89*t_0^12*t_1^6*t_2^10+46*t_0^12*t_1^5*t_2^11+16*t_0^12*t_1^4*t_2^12+3*t_0^12*t_1^3*t_2^13+2*t_0^11*t_1^14*t_2^3+14*t_0^11*t_1^13*t_2^4+46*t_0^11*t_1^12*t_2^5+100*t_0^11*t_1^11*t_2^6+163*t_0^11*t_1^10*t_2^7+206*t_0^11*t_1^9*t_2^8+206*t_0^11*t_1^8*t_2^9+163*t_0^11*t_1^7*t_2^10+100*t_0^11*t_1^6*t_2^11+46*t_0^11*t_1^5*t_2^12+14*t_0^11*t_1^4*t_2^13+2*t_0^11*t_1^3*t_2^14+t_0^10*t_1^15*t_2^3+9*t_0^10*t_1^14*t_2^4+36*t_0^10*t_1^13*t_2^5+89*t_0^10*t_1^12*t_2^6+163*t_0^10*t_1^11*t_2^7+230*t_0^10*t_1^10*t_2^8+256*t_0^10*t_1^9*t_2^9+230*t_0^10*t_1^8*t_2^10+163*t_0^10*t_1^7*t_2^11+89*t_0^10*t_1^6*t_2^12+36*t_0^10*t_1^5*t_2^13+9*t_0^10*t_1^4*t_2^14+t_0^10*t_1^3*t_2^15+4*t_0^9*t_1^15*t_2^4+21*t_0^9*t_1^14*t_2^5+62*t_0^9*t_1^13*t_2^6+130*t_0^9*t_1^12*t_2^7+206*t_0^9*t_1^11*t_2^8+256*t_0^9*t_1^10*t_2^9+256*t_0^9*t_1^9*t_2^10+206*t_0^9*t_1^8*t_2^11+130*t_0^9*t_1^7*t_2^12+62*t_0^9*t_1^6*t_2^13+21*t_0^9*t_1^5*t_2^14+4*t_0^9*t_1^4*t_2^15+t_0^8*t_1^16*t_2^4+9*t_0^8*t_1^15*t_2^5+33*t_0^8*t_1^14*t_2^6+82*t_0^8*t_1^13*t_2^7+148*t_0^8*t_1^12*t_2^8+206*t_0^8*t_1^11*t_2^9+230*t_0^8*t_1^10*t_2^10+206*t_0^8*t_1^9*t_2^11+148*t_0^8*t_1^8*t_2^12+82*t_0^8*t_1^7*t_2^13+33*t_0^8*t_1^6*t_2^14+9*t_0^8*t_1^5*t_2^15+t_0^8*t_1^4*t_2^16+2*t_0^7*t_1^16*t_2^5+12*t_0^7*t_1^15*t_2^6+38*t_0^7*t_1^14*t_2^7+82*t_0^7*t_1^13*t_2^8+130*t_0^7*t_1^12*t_2^9+163*t_0^7*t_1^11*t_2^10+163*t_0^7*t_1^10*t_2^11+130*t_0^7*t_1^9*t_2^12+82*t_0^7*t_1^8*t_2^13+38*t_0^7*t_1^7*t_2^14+12*t_0^7*t_1^6*t_2^15+2*t_0^7*t_1^5*t_2^16+2*t_0^6*t_1^16*t_2^6+12*t_0^6*t_1^15*t_2^7+33*t_0^6*t_1^14*t_2^8+62*t_0^6*t_1^13*t_2^9+89*t_0^6*t_1^12*t_2^10+100*t_0^6*t_1^11*t_2^11+89*t_0^6*t_1^10*t_2^12+62*t_0^6*t_1^9*t_2^13+33*t_0^6*t_1^8*t_2^14+12*t_0^6*t_1^7*t_2^15+2*t_0^6*t_1^6*t_2^16+2*t_0^5*t_1^16*t_2^7+9*t_0^5*t_1^15*t_2^8+21*t_0^5*t_1^14*t_2^9+36*t_0^5*t_1^13*t_2^10+46*t_0^5*t_1^12*t_2^11+46*t_0^5*t_1^11*t_2^12+36*t_0^5*t_1^10*t_2^13+21*t_0^5*t_1^9*t_2^14+9*t_0^5*t_1^8*t_2^15+2*t_0^5*t_1^7*t_2^16+t_0^4*t_1^16*t_2^8+4*t_0^4*t_1^15*t_2^9+9*t_0^4*t_1^14*t_2^10+14*t_0^4*t_1^13*t_2^11+16*t_0^4*t_1^12*t_2^12+14*t_0^4*t_1^11*t_2^13+9*t_0^4*t_1^10*t_2^14+4*t_0^4*t_1^9*t_2^15+t_0^4*t_1^8*t_2^16+t_0^3*t_1^15*t_2^10+2*t_0^3*t_1^14*t_2^11+3*t_0^3*t_1^13*t_2^12+3*t_0^3*t_1^12*t_2^13+2*t_0^3*t_1^11*t_2^14+t_0^3*t_1^10*t_2^15, (5,2) => 0, (8,0) => 0, (6,2) => 0, (7,1) => t_0^17*t_1^10*t_2^5+2*t_0^17*t_1^9*t_2^6+3*t_0^17*t_1^8*t_2^7+3*t_0^17*t_1^7*t_2^8+2*t_0^17*t_1^6*t_2^9+t_0^17*t_1^5*t_2^10+3*t_0^16*t_1^11*t_2^5+7*t_0^16*t_1^10*t_2^6+12*t_0^16*t_1^9*t_2^7+14*t_0^16*t_1^8*t_2^8+12*t_0^16*t_1^7*t_2^9+7*t_0^16*t_1^6*t_2^10+3*t_0^16*t_1^5*t_2^11+t_0^15*t_1^13*t_2^4+6*t_0^15*t_1^12*t_2^5+16*t_0^15*t_1^11*t_2^6+29*t_0^15*t_1^10*t_2^7+39*t_0^15*t_1^9*t_2^8+39*t_0^15*t_1^8*t_2^9+29*t_0^15*t_1^7*t_2^10+16*t_0^15*t_1^6*t_2^11+6*t_0^15*t_1^5*t_2^12+t_0^15*t_1^4*t_2^13+t_0^14*t_1^14*t_2^4+8*t_0^14*t_1^13*t_2^5+23*t_0^14*t_1^12*t_2^6+48*t_0^14*t_1^11*t_2^7+71*t_0^14*t_1^10*t_2^8+82*t_0^14*t_1^9*t_2^9+71*t_0^14*t_1^8*t_2^10+48*t_0^14*t_1^7*t_2^11+23*t_0^14*t_1^6*t_2^12+8*t_0^14*t_1^5*t_2^13+t_0^14*t_1^4*t_2^14+t_0^13*t_1^15*t_2^4+8*t_0^13*t_1^14*t_2^5+28*t_0^13*t_1^13*t_2^6+62*t_0^13*t_1^12*t_2^7+103*t_0^13*t_1^11*t_2^8+131*t_0^13*t_1^10*t_2^9+131*t_0^13*t_1^9*t_2^10+103*t_0^13*t_1^8*t_2^11+62*t_0^13*t_1^7*t_2^12+28*t_0^13*t_1^6*t_2^13+8*t_0^13*t_1^5*t_2^14+t_0^13*t_1^4*t_2^15+6*t_0^12*t_1^15*t_2^5+23*t_0^12*t_1^14*t_2^6+62*t_0^12*t_1^13*t_2^7+113*t_0^12*t_1^12*t_2^8+163*t_0^12*t_1^11*t_2^9+181*t_0^12*t_1^10*t_2^10+163*t_0^12*t_1^9*t_2^11+113*t_0^12*t_1^8*t_2^12+62*t_0^12*t_1^7*t_2^13+23*t_0^12*t_1^6*t_2^14+6*t_0^12*t_1^5*t_2^15+3*t_0^11*t_1^16*t_2^5+16*t_0^11*t_1^15*t_2^6+48*t_0^11*t_1^14*t_2^7+103*t_0^11*t_1^13*t_2^8+163*t_0^11*t_1^12*t_2^9+204*t_0^11*t_1^11*t_2^10+204*t_0^11*t_1^10*t_2^11+163*t_0^11*t_1^9*t_2^12+103*t_0^11*t_1^8*t_2^13+48*t_0^11*t_1^7*t_2^14+16*t_0^11*t_1^6*t_2^15+3*t_0^11*t_1^5*t_2^16+t_0^10*t_1^17*t_2^5+7*t_0^10*t_1^16*t_2^6+29*t_0^10*t_1^15*t_2^7+71*t_0^10*t_1^14*t_2^8+131*t_0^10*t_1^13*t_2^9+181*t_0^10*t_1^12*t_2^10+204*t_0^10*t_1^11*t_2^11+181*t_0^10*t_1^10*t_2^12+131*t_0^10*t_1^9*t_2^13+71*t_0^10*t_1^8*t_2^14+29*t_0^10*t_1^7*t_2^15+7*t_0^10*t_1^6*t_2^16+t_0^10*t_1^5*t_2^17+2*t_0^9*t_1^17*t_2^6+12*t_0^9*t_1^16*t_2^7+39*t_0^9*t_1^15*t_2^8+82*t_0^9*t_1^14*t_2^9+131*t_0^9*t_1^13*t_2^10+163*t_0^9*t_1^12*t_2^11+163*t_0^9*t_1^11*t_2^12+131*t_0^9*t_1^10*t_2^13+82*t_0^9*t_1^9*t_2^14+39*t_0^9*t_1^8*t_2^15+12*t_0^9*t_1^7*t_2^16+2*t_0^9*t_1^6*t_2^17+3*t_0^8*t_1^17*t_2^7+14*t_0^8*t_1^16*t_2^8+39*t_0^8*t_1^15*t_2^9+71*t_0^8*t_1^14*t_2^10+103*t_0^8*t_1^13*t_2^11+113*t_0^8*t_1^12*t_2^12+103*t_0^8*t_1^11*t_2^13+71*t_0^8*t_1^10*t_2^14+39*t_0^8*t_1^9*t_2^15+14*t_0^8*t_1^8*t_2^16+3*t_0^8*t_1^7*t_2^17+3*t_0^7*t_1^17*t_2^8+12*t_0^7*t_1^16*t_2^9+29*t_0^7*t_1^15*t_2^10+48*t_0^7*t_1^14*t_2^11+62*t_0^7*t_1^13*t_2^12+62*t_0^7*t_1^12*t_2^13+48*t_0^7*t_1^11*t_2^14+29*t_0^7*t_1^10*t_2^15+12*t_0^7*t_1^9*t_2^16+3*t_0^7*t_1^8*t_2^17+2*t_0^6*t_1^17*t_2^9+7*t_0^6*t_1^16*t_2^10+16*t_0^6*t_1^15*t_2^11+23*t_0^6*t_1^14*t_2^12+28*t_0^6*t_1^13*t_2^13+23*t_0^6*t_1^12*t_2^14+16*t_0^6*t_1^11*t_2^15+7*t_0^6*t_1^10*t_2^16+2*t_0^6*t_1^9*t_2^17+t_0^5*t_1^17*t_2^10+3*t_0^5*t_1^16*t_2^11+6*t_0^5*t_1^15*t_2^12+8*t_0^5*t_1^14*t_2^13+8*t_0^5*t_1^13*t_2^14+6*t_0^5*t_1^12*t_2^15+3*t_0^5*t_1^11*t_2^16+t_0^5*t_1^10*t_2^17+t_0^4*t_1^15*t_2^13+t_0^4*t_1^14*t_2^14+t_0^4*t_1^13*t_2^15, (7,2) => 0, (8,1) => t_0^18*t_1^11*t_2^7+2*t_0^18*t_1^10*t_2^8+2*t_0^18*t_1^9*t_2^9+2*t_0^18*t_1^8*t_2^10+t_0^18*t_1^7*t_2^11+t_0^17*t_1^13*t_2^6+3*t_0^17*t_1^12*t_2^7+7*t_0^17*t_1^11*t_2^8+10*t_0^17*t_1^10*t_2^9+10*t_0^17*t_1^9*t_2^10+7*t_0^17*t_1^8*t_2^11+3*t_0^17*t_1^7*t_2^12+t_0^17*t_1^6*t_2^13+2*t_0^16*t_1^14*t_2^6+7*t_0^16*t_1^13*t_2^7+15*t_0^16*t_1^12*t_2^8+24*t_0^16*t_1^11*t_2^9+28*t_0^16*t_1^10*t_2^10+24*t_0^16*t_1^9*t_2^11+15*t_0^16*t_1^8*t_2^12+7*t_0^16*t_1^7*t_2^13+2*t_0^16*t_1^6*t_2^14+2*t_0^15*t_1^15*t_2^6+9*t_0^15*t_1^14*t_2^7+23*t_0^15*t_1^13*t_2^8+39*t_0^15*t_1^12*t_2^9+52*t_0^15*t_1^11*t_2^10+52*t_0^15*t_1^10*t_2^11+39*t_0^15*t_1^9*t_2^12+23*t_0^15*t_1^8*t_2^13+9*t_0^15*t_1^7*t_2^14+2*t_0^15*t_1^6*t_2^15+2*t_0^14*t_1^16*t_2^6+9*t_0^14*t_1^15*t_2^7+26*t_0^14*t_1^14*t_2^8+51*t_0^14*t_1^13*t_2^9+74*t_0^14*t_1^12*t_2^10+84*t_0^14*t_1^11*t_2^11+74*t_0^14*t_1^10*t_2^12+51*t_0^14*t_1^9*t_2^13+26*t_0^14*t_1^8*t_2^14+9*t_0^14*t_1^7*t_2^15+2*t_0^14*t_1^6*t_2^16+t_0^13*t_1^17*t_2^6+7*t_0^13*t_1^16*t_2^7+23*t_0^13*t_1^15*t_2^8+51*t_0^13*t_1^14*t_2^9+84*t_0^13*t_1^13*t_2^10+105*t_0^13*t_1^12*t_2^11+105*t_0^13*t_1^11*t_2^12+84*t_0^13*t_1^10*t_2^13+51*t_0^13*t_1^9*t_2^14+23*t_0^13*t_1^8*t_2^15+7*t_0^13*t_1^7*t_2^16+t_0^13*t_1^6*t_2^17+3*t_0^12*t_1^17*t_2^7+15*t_0^12*t_1^16*t_2^8+39*t_0^12*t_1^15*t_2^9+74*t_0^12*t_1^14*t_2^10+105*t_0^12*t_1^13*t_2^11+116*t_0^12*t_1^12*t_2^12+105*t_0^12*t_1^11*t_2^13+74*t_0^12*t_1^10*t_2^14+39*t_0^12*t_1^9*t_2^15+15*t_0^12*t_1^8*t_2^16+3*t_0^12*t_1^7*t_2^17+t_0^11*t_1^18*t_2^7+7*t_0^11*t_1^17*t_2^8+24*t_0^11*t_1^16*t_2^9+52*t_0^11*t_1^15*t_2^10+84*t_0^11*t_1^14*t_2^11+105*t_0^11*t_1^13*t_2^12+105*t_0^11*t_1^12*t_2^13+84*t_0^11*t_1^11*t_2^14+52*t_0^11*t_1^10*t_2^15+24*t_0^11*t_1^9*t_2^16+7*t_0^11*t_1^8*t_2^17+t_0^11*t_1^7*t_2^18+2*t_0^10*t_1^18*t_2^8+10*t_0^10*t_1^17*t_2^9+28*t_0^10*t_1^16*t_2^10+52*t_0^10*t_1^15*t_2^11+74*t_0^10*t_1^14*t_2^12+84*t_0^10*t_1^13*t_2^13+74*t_0^10*t_1^12*t_2^14+52*t_0^10*t_1^11*t_2^15+28*t_0^10*t_1^10*t_2^16+10*t_0^10*t_1^9*t_2^17+2*t_0^10*t_1^8*t_2^18+2*t_0^9*t_1^18*t_2^9+10*t_0^9*t_1^17*t_2^10+24*t_0^9*t_1^16*t_2^11+39*t_0^9*t_1^15*t_2^12+51*t_0^9*t_1^14*t_2^13+51*t_0^9*t_1^13*t_2^14+39*t_0^9*t_1^12*t_2^15+24*t_0^9*t_1^11*t_2^16+10*t_0^9*t_1^10*t_2^17+2*t_0^9*t_1^9*t_2^18+2*t_0^8*t_1^18*t_2^10+7*t_0^8*t_1^17*t_2^11+15*t_0^8*t_1^16*t_2^12+23*t_0^8*t_1^15*t_2^13+26*t_0^8*t_1^14*t_2^14+23*t_0^8*t_1^13*t_2^15+15*t_0^8*t_1^12*t_2^16+7*t_0^8*t_1^11*t_2^17+2*t_0^8*t_1^10*t_2^18+t_0^7*t_1^18*t_2^11+3*t_0^7*t_1^17*t_2^12+7*t_0^7*t_1^16*t_2^13+9*t_0^7*t_1^15*t_2^14+9*t_0^7*t_1^14*t_2^15+7*t_0^7*t_1^13*t_2^16+3*t_0^7*t_1^12*t_2^17+t_0^7*t_1^11*t_2^18+t_0^6*t_1^17*t_2^13+2*t_0^6*t_1^16*t_2^14+2*t_0^6*t_1^15*t_2^15+2*t_0^6*t_1^14*t_2^16+t_0^6*t_1^13*t_2^17, (9,0) => 0, (8,2) => 0, (9,1) => t_0^19*t_1^11*t_2^10+t_0^19*t_1^10*t_2^11+t_0^18*t_1^14*t_2^8+2*t_0^18*t_1^13*t_2^9+4*t_0^18*t_1^12*t_2^10+5*t_0^18*t_1^11*t_2^11+4*t_0^18*t_1^10*t_2^12+2*t_0^18*t_1^9*t_2^13+t_0^18*t_1^8*t_2^14+t_0^17*t_1^15*t_2^8+4*t_0^17*t_1^14*t_2^9+8*t_0^17*t_1^13*t_2^10+11*t_0^17*t_1^12*t_2^11+11*t_0^17*t_1^11*t_2^12+8*t_0^17*t_1^10*t_2^13+4*t_0^17*t_1^9*t_2^14+t_0^17*t_1^8*t_2^15+2*t_0^16*t_1^16*t_2^8+5*t_0^16*t_1^15*t_2^9+13*t_0^16*t_1^14*t_2^10+19*t_0^16*t_1^13*t_2^11+22*t_0^16*t_1^12*t_2^12+19*t_0^16*t_1^11*t_2^13+13*t_0^16*t_1^10*t_2^14+5*t_0^16*t_1^9*t_2^15+2*t_0^16*t_1^8*t_2^16+t_0^15*t_1^17*t_2^8+5*t_0^15*t_1^16*t_2^9+13*t_0^15*t_1^15*t_2^10+24*t_0^15*t_1^14*t_2^11+30*t_0^15*t_1^13*t_2^12+30*t_0^15*t_1^12*t_2^13+24*t_0^15*t_1^11*t_2^14+13*t_0^15*t_1^10*t_2^15+5*t_0^15*t_1^9*t_2^16+t_0^15*t_1^8*t_2^17+t_0^14*t_1^18*t_2^8+4*t_0^14*t_1^17*t_2^9+13*t_0^14*t_1^16*t_2^10+24*t_0^14*t_1^15*t_2^11+36*t_0^14*t_1^14*t_2^12+39*t_0^14*t_1^13*t_2^13+36*t_0^14*t_1^12*t_2^14+24*t_0^14*t_1^11*t_2^15+13*t_0^14*t_1^10*t_2^16+4*t_0^14*t_1^9*t_2^17+t_0^14*t_1^8*t_2^18+2*t_0^13*t_1^18*t_2^9+8*t_0^13*t_1^17*t_2^10+19*t_0^13*t_1^16*t_2^11+30*t_0^13*t_1^15*t_2^12+39*t_0^13*t_1^14*t_2^13+39*t_0^13*t_1^13*t_2^14+30*t_0^13*t_1^12*t_2^15+19*t_0^13*t_1^11*t_2^16+8*t_0^13*t_1^10*t_2^17+2*t_0^13*t_1^9*t_2^18+4*t_0^12*t_1^18*t_2^10+11*t_0^12*t_1^17*t_2^11+22*t_0^12*t_1^16*t_2^12+30*t_0^12*t_1^15*t_2^13+36*t_0^12*t_1^14*t_2^14+30*t_0^12*t_1^13*t_2^15+22*t_0^12*t_1^12*t_2^16+11*t_0^12*t_1^11*t_2^17+4*t_0^12*t_1^10*t_2^18+t_0^11*t_1^19*t_2^10+5*t_0^11*t_1^18*t_2^11+11*t_0^11*t_1^17*t_2^12+19*t_0^11*t_1^16*t_2^13+24*t_0^11*t_1^15*t_2^14+24*t_0^11*t_1^14*t_2^15+19*t_0^11*t_1^13*t_2^16+11*t_0^11*t_1^12*t_2^17+5*t_0^11*t_1^11*t_2^18+t_0^11*t_1^10*t_2^19+t_0^10*t_1^19*t_2^11+4*t_0^10*t_1^18*t_2^12+8*t_0^10*t_1^17*t_2^13+13*t_0^10*t_1^16*t_2^14+13*t_0^10*t_1^15*t_2^15+13*t_0^10*t_1^14*t_2^16+8*t_0^10*t_1^13*t_2^17+4*t_0^10*t_1^12*t_2^18+t_0^10*t_1^11*t_2^19+2*t_0^9*t_1^18*t_2^13+4*t_0^9*t_1^17*t_2^14+5*t_0^9*t_1^16*t_2^15+5*t_0^9*t_1^15*t_2^16+4*t_0^9*t_1^14*t_2^17+2*t_0^9*t_1^13*t_2^18+t_0^8*t_1^18*t_2^14+t_0^8*t_1^17*t_2^15+2*t_0^8*t_1^16*t_2^16+t_0^8*t_1^15*t_2^17+t_0^8*t_1^14*t_2^18, (10,0) => 0, (11,0) => 0, (10,1) => t_0^19*t_1^14*t_2^11+t_0^19*t_1^13*t_2^12+t_0^19*t_1^12*t_2^13+t_0^19*t_1^11*t_2^14+t_0^18*t_1^15*t_2^11+2*t_0^18*t_1^14*t_2^12+2*t_0^18*t_1^13*t_2^13+2*t_0^18*t_1^12*t_2^14+t_0^18*t_1^11*t_2^15+t_0^17*t_1^16*t_2^11+2*t_0^17*t_1^15*t_2^12+3*t_0^17*t_1^14*t_2^13+3*t_0^17*t_1^13*t_2^14+2*t_0^17*t_1^12*t_2^15+t_0^17*t_1^11*t_2^16+t_0^16*t_1^17*t_2^11+2*t_0^16*t_1^16*t_2^12+3*t_0^16*t_1^15*t_2^13+4*t_0^16*t_1^14*t_2^14+3*t_0^16*t_1^13*t_2^15+2*t_0^16*t_1^12*t_2^16+t_0^16*t_1^11*t_2^17+t_0^15*t_1^18*t_2^11+2*t_0^15*t_1^17*t_2^12+3*t_0^15*t_1^16*t_2^13+4*t_0^15*t_1^15*t_2^14+4*t_0^15*t_1^14*t_2^15+3*t_0^15*t_1^13*t_2^16+2*t_0^15*t_1^12*t_2^17+t_0^15*t_1^11*t_2^18+t_0^14*t_1^19*t_2^11+2*t_0^14*t_1^18*t_2^12+3*t_0^14*t_1^17*t_2^13+4*t_0^14*t_1^16*t_2^14+4*t_0^14*t_1^15*t_2^15+4*t_0^14*t_1^14*t_2^16+3*t_0^14*t_1^13*t_2^17+2*t_0^14*t_1^12*t_2^18+t_0^14*t_1^11*t_2^19+t_0^13*t_1^19*t_2^12+2*t_0^13*t_1^18*t_2^13+3*t_0^13*t_1^17*t_2^14+3*t_0^13*t_1^16*t_2^15+3*t_0^13*t_1^15*t_2^16+3*t_0^13*t_1^14*t_2^17+2*t_0^13*t_1^13*t_2^18+t_0^13*t_1^12*t_2^19+t_0^12*t_1^19*t_2^13+2*t_0^12*t_1^18*t_2^14+2*t_0^12*t_1^17*t_2^15+2*t_0^12*t_1^16*t_2^16+2*t_0^12*t_1^15*t_2^17+2*t_0^12*t_1^14*t_2^18+t_0^12*t_1^13*t_2^19+t_0^11*t_1^19*t_2^14+t_0^11*t_1^18*t_2^15+t_0^11*t_1^17*t_2^16+t_0^11*t_1^16*t_2^17+t_0^11*t_1^15*t_2^18+t_0^11*t_1^14*t_2^19, (9,2) => 0, (12,0) => 0, (11,1) => 0, (10,2) => t_0^18*t_1^18*t_2^12+t_0^18*t_1^17*t_2^13+2*t_0^18*t_1^16*t_2^14+2*t_0^18*t_1^15*t_2^15+2*t_0^18*t_1^14*t_2^16+t_0^18*t_1^13*t_2^17+t_0^18*t_1^12*t_2^18+t_0^17*t_1^18*t_2^13+2*t_0^17*t_1^17*t_2^14+3*t_0^17*t_1^16*t_2^15+3*t_0^17*t_1^15*t_2^16+2*t_0^17*t_1^14*t_2^17+t_0^17*t_1^13*t_2^18+2*t_0^16*t_1^18*t_2^14+3*t_0^16*t_1^17*t_2^15+4*t_0^16*t_1^16*t_2^16+3*t_0^16*t_1^15*t_2^17+2*t_0^16*t_1^14*t_2^18+2*t_0^15*t_1^18*t_2^15+3*t_0^15*t_1^17*t_2^16+3*t_0^15*t_1^16*t_2^17+2*t_0^15*t_1^15*t_2^18+2*t_0^14*t_1^18*t_2^16+2*t_0^14*t_1^17*t_2^17+2*t_0^14*t_1^16*t_2^18+t_0^13*t_1^18*t_2^17+t_0^13*t_1^17*t_2^18+t_0^12*t_1^18*t_2^18, (12,1) => 0, (11,2) => t_0^19*t_1^18*t_2^15+t_0^19*t_1^17*t_2^16+t_0^19*t_1^16*t_2^17+t_0^19*t_1^15*t_2^18+t_0^18*t_1^19*t_2^15+2*t_0^18*t_1^18*t_2^16+2*t_0^18*t_1^17*t_2^17+2*t_0^18*t_1^16*t_2^18+t_0^18*t_1^15*t_2^19+t_0^17*t_1^19*t_2^16+2*t_0^17*t_1^18*t_2^17+2*t_0^17*t_1^17*t_2^18+t_0^17*t_1^16*t_2^19+t_0^16*t_1^19*t_2^17+2*t_0^16*t_1^18*t_2^18+t_0^16*t_1^17*t_2^19+t_0^15*t_1^19*t_2^18+t_0^15*t_1^18*t_2^19, (12,2) => t_0^19*t_1^19*t_2^18+t_0^19*t_1^18*t_2^19+t_0^18*t_1^19*t_2^19, (0,0) => 1, (0,1) => 0, (1,0) => 0, (1,1) => t_0^6*t_1^2+t_0^6*t_1*t_2+t_0^6*t_2^2+t_0^5*t_1^3+2*t_0^5*t_1^2*t_2+2*t_0^5*t_1*t_2^2+t_0^5*t_2^3+2*t_0^4*t_1^4+3*t_0^4*t_1^3*t_2+4*t_0^4*t_1^2*t_2^2+3*t_0^4*t_1*t_2^3+2*t_0^4*t_2^4+t_0^3*t_1^5+3*t_0^3*t_1^4*t_2+4*t_0^3*t_1^3*t_2^2+4*t_0^3*t_1^2*t_2^3+3*t_0^3*t_1*t_2^4+t_0^3*t_2^5+t_0^2*t_1^6+2*t_0^2*t_1^5*t_2+4*t_0^2*t_1^4*t_2^2+4*t_0^2*t_1^3*t_2^3+4*t_0^2*t_1^2*t_2^4+2*t_0^2*t_1*t_2^5+t_0^2*t_2^6+t_0*t_1^6*t_2+2*t_0*t_1^5*t_2^2+3*t_0*t_1^4*t_2^3+3*t_0*t_1^3*t_2^4+2*t_0*t_1^2*t_2^5+t_0*t_1*t_2^6+t_1^6*t_2^2+t_1^5*t_2^3+2*t_1^4*t_2^4+t_1^3*t_2^5+t_1^2*t_2^6, (0,2) => 0, (2,0) => 0, (1,2) => 0, (2,1) => t_0^9*t_1^2*t_2+t_0^9*t_1*t_2^2+t_0^8*t_1^4+3*t_0^8*t_1^3*t_2+4*t_0^8*t_1^2*t_2^2+3*t_0^8*t_1*t_2^3+t_0^8*t_2^4+2*t_0^7*t_1^5+6*t_0^7*t_1^4*t_2+9*t_0^7*t_1^3*t_2^2+9*t_0^7*t_1^2*t_2^3+6*t_0^7*t_1*t_2^4+2*t_0^7*t_2^5+2*t_0^6*t_1^6+8*t_0^6*t_1^5*t_2+14*t_0^6*t_1^4*t_2^2+16*t_0^6*t_1^3*t_2^3+14*t_0^6*t_1^2*t_2^4+8*t_0^6*t_1*t_2^5+2*t_0^6*t_2^6+2*t_0^5*t_1^7+8*t_0^5*t_1^6*t_2+16*t_0^5*t_1^5*t_2^2+22*t_0^5*t_1^4*t_2^3+22*t_0^5*t_1^3*t_2^4+16*t_0^5*t_1^2*t_2^5+8*t_0^5*t_1*t_2^6+2*t_0^5*t_2^7+t_0^4*t_1^8+6*t_0^4*t_1^7*t_2+14*t_0^4*t_1^6*t_2^2+22*t_0^4*t_1^5*t_2^3+26*t_0^4*t_1^4*t_2^4+22*t_0^4*t_1^3*t_2^5+14*t_0^4*t_1^2*t_2^6+6*t_0^4*t_1*t_2^7+t_0^4*t_2^8+3*t_0^3*t_1^8*t_2+9*t_0^3*t_1^7*t_2^2+16*t_0^3*t_1^6*t_2^3+22*t_0^3*t_1^5*t_2^4+22*t_0^3*t_1^4*t_2^5+16*t_0^3*t_1^3*t_2^6+9*t_0^3*t_1^2*t_2^7+3*t_0^3*t_1*t_2^8+t_0^2*t_1^9*t_2+4*t_0^2*t_1^8*t_2^2+9*t_0^2*t_1^7*t_2^3+14*t_0^2*t_1^6*t_2^4+16*t_0^2*t_1^5*t_2^5+14*t_0^2*t_1^4*t_2^6+9*t_0^2*t_1^3*t_2^7+4*t_0^2*t_1^2*t_2^8+t_0^2*t_1*t_2^9+t_0*t_1^9*t_2^2+3*t_0*t_1^8*t_2^3+6*t_0*t_1^7*t_2^4+8*t_0*t_1^6*t_2^5+8*t_0*t_1^5*t_2^6+6*t_0*t_1^4*t_2^7+3*t_0*t_1^3*t_2^8+t_0*t_1^2*t_2^9+t_1^8*t_2^4+2*t_1^7*t_2^5+2*t_1^6*t_2^6+2*t_1^5*t_2^7+t_1^4*t_2^8, (3,0) => 0, (4,0) => 0, (2,2) => 0, (3,1) => t_0^11*t_1^4*t_2+2*t_0^11*t_1^3*t_2^2+2*t_0^11*t_1^2*t_2^3+t_0^11*t_1*t_2^4+3*t_0^10*t_1^5*t_2+6*t_0^10*t_1^4*t_2^2+8*t_0^10*t_1^3*t_2^3+6*t_0^10*t_1^2*t_2^4+3*t_0^10*t_1*t_2^5+t_0^9*t_1^7+6*t_0^9*t_1^6*t_2+14*t_0^9*t_1^5*t_2^2+20*t_0^9*t_1^4*t_2^3+20*t_0^9*t_1^3*t_2^4+14*t_0^9*t_1^2*t_2^5+6*t_0^9*t_1*t_2^6+t_0^9*t_2^7+t_0^8*t_1^8+8*t_0^8*t_1^7*t_2+20*t_0^8*t_1^6*t_2^2+34*t_0^8*t_1^5*t_2^3+39*t_0^8*t_1^4*t_2^4+34*t_0^8*t_1^3*t_2^5+20*t_0^8*t_1^2*t_2^6+8*t_0^8*t_1*t_2^7+t_0^8*t_2^8+t_0^7*t_1^9+8*t_0^7*t_1^8*t_2+24*t_0^7*t_1^7*t_2^2+44*t_0^7*t_1^6*t_2^3+59*t_0^7*t_1^5*t_2^4+59*t_0^7*t_1^4*t_2^5+44*t_0^7*t_1^3*t_2^6+24*t_0^7*t_1^2*t_2^7+8*t_0^7*t_1*t_2^8+t_0^7*t_2^9+6*t_0^6*t_1^9*t_2+20*t_0^6*t_1^8*t_2^2+44*t_0^6*t_1^7*t_2^3+65*t_0^6*t_1^6*t_2^4+76*t_0^6*t_1^5*t_2^5+65*t_0^6*t_1^4*t_2^6+44*t_0^6*t_1^3*t_2^7+20*t_0^6*t_1^2*t_2^8+6*t_0^6*t_1*t_2^9+3*t_0^5*t_1^10*t_2+14*t_0^5*t_1^9*t_2^2+34*t_0^5*t_1^8*t_2^3+59*t_0^5*t_1^7*t_2^4+76*t_0^5*t_1^6*t_2^5+76*t_0^5*t_1^5*t_2^6+59*t_0^5*t_1^4*t_2^7+34*t_0^5*t_1^3*t_2^8+14*t_0^5*t_1^2*t_2^9+3*t_0^5*t_1*t_2^10+t_0^4*t_1^11*t_2+6*t_0^4*t_1^10*t_2^2+20*t_0^4*t_1^9*t_2^3+39*t_0^4*t_1^8*t_2^4+59*t_0^4*t_1^7*t_2^5+65*t_0^4*t_1^6*t_2^6+59*t_0^4*t_1^5*t_2^7+39*t_0^4*t_1^4*t_2^8+20*t_0^4*t_1^3*t_2^9+6*t_0^4*t_1^2*t_2^10+t_0^4*t_1*t_2^11+2*t_0^3*t_1^11*t_2^2+8*t_0^3*t_1^10*t_2^3+20*t_0^3*t_1^9*t_2^4+34*t_0^3*t_1^8*t_2^5+44*t_0^3*t_1^7*t_2^6+44*t_0^3*t_1^6*t_2^7+34*t_0^3*t_1^5*t_2^8+20*t_0^3*t_1^4*t_2^9+8*t_0^3*t_1^3*t_2^10+2*t_0^3*t_1^2*t_2^11+2*t_0^2*t_1^11*t_2^3+6*t_0^2*t_1^10*t_2^4+14*t_0^2*t_1^9*t_2^5+20*t_0^2*t_1^8*t_2^6+24*t_0^2*t_1^7*t_2^7+20*t_0^2*t_1^6*t_2^8+14*t_0^2*t_1^5*t_2^9+6*t_0^2*t_1^4*t_2^10+2*t_0^2*t_1^3*t_2^11+t_0*t_1^11*t_2^4+3*t_0*t_1^10*t_2^5+6*t_0*t_1^9*t_2^6+8*t_0*t_1^8*t_2^7+8*t_0*t_1^7*t_2^8+6*t_0*t_1^6*t_2^9+3*t_0*t_1^5*t_2^10+t_0*t_1^4*t_2^11+t_1^9*t_2^7+t_1^8*t_2^8+t_1^7*t_2^9, (5,0) => 0, (4,1) => t_0^13*t_1^5*t_2^2+2*t_0^13*t_1^4*t_2^3+2*t_0^13*t_1^3*t_2^4+t_0^13*t_1^2*t_2^5+t_0^12*t_1^7*t_2+4*t_0^12*t_1^6*t_2^2+8*t_0^12*t_1^5*t_2^3+10*t_0^12*t_1^4*t_2^4+8*t_0^12*t_1^3*t_2^5+4*t_0^12*t_1^2*t_2^6+t_0^12*t_1*t_2^7+2*t_0^11*t_1^8*t_2+9*t_0^11*t_1^7*t_2^2+19*t_0^11*t_1^6*t_2^3+27*t_0^11*t_1^5*t_2^4+27*t_0^11*t_1^4*t_2^5+19*t_0^11*t_1^3*t_2^6+9*t_0^11*t_1^2*t_2^7+2*t_0^11*t_1*t_2^8+3*t_0^10*t_1^9*t_2+14*t_0^10*t_1^8*t_2^2+33*t_0^10*t_1^7*t_2^3+52*t_0^10*t_1^6*t_2^4+60*t_0^10*t_1^5*t_2^5+52*t_0^10*t_1^4*t_2^6+33*t_0^10*t_1^3*t_2^7+14*t_0^10*t_1^2*t_2^8+3*t_0^10*t_1*t_2^9+3*t_0^9*t_1^10*t_2+16*t_0^9*t_1^9*t_2^2+42*t_0^9*t_1^8*t_2^3+75*t_0^9*t_1^7*t_2^4+98*t_0^9*t_1^6*t_2^5+98*t_0^9*t_1^5*t_2^6+75*t_0^9*t_1^4*t_2^7+42*t_0^9*t_1^3*t_2^8+16*t_0^9*t_1^2*t_2^9+3*t_0^9*t_1*t_2^10+2*t_0^8*t_1^11*t_2+14*t_0^8*t_1^10*t_2^2+42*t_0^8*t_1^9*t_2^3+84*t_0^8*t_1^8*t_2^4+124*t_0^8*t_1^7*t_2^5+140*t_0^8*t_1^6*t_2^6+124*t_0^8*t_1^5*t_2^7+84*t_0^8*t_1^4*t_2^8+42*t_0^8*t_1^3*t_2^9+14*t_0^8*t_1^2*t_2^10+2*t_0^8*t_1*t_2^11+t_0^7*t_1^12*t_2+9*t_0^7*t_1^11*t_2^2+33*t_0^7*t_1^10*t_2^3+75*t_0^7*t_1^9*t_2^4+124*t_0^7*t_1^8*t_2^5+158*t_0^7*t_1^7*t_2^6+158*t_0^7*t_1^6*t_2^7+124*t_0^7*t_1^5*t_2^8+75*t_0^7*t_1^4*t_2^9+33*t_0^7*t_1^3*t_2^10+9*t_0^7*t_1^2*t_2^11+t_0^7*t_1*t_2^12+4*t_0^6*t_1^12*t_2^2+19*t_0^6*t_1^11*t_2^3+52*t_0^6*t_1^10*t_2^4+98*t_0^6*t_1^9*t_2^5+140*t_0^6*t_1^8*t_2^6+158*t_0^6*t_1^7*t_2^7+140*t_0^6*t_1^6*t_2^8+98*t_0^6*t_1^5*t_2^9+52*t_0^6*t_1^4*t_2^10+19*t_0^6*t_1^3*t_2^11+4*t_0^6*t_1^2*t_2^12+t_0^5*t_1^13*t_2^2+8*t_0^5*t_1^12*t_2^3+27*t_0^5*t_1^11*t_2^4+60*t_0^5*t_1^10*t_2^5+98*t_0^5*t_1^9*t_2^6+124*t_0^5*t_1^8*t_2^7+124*t_0^5*t_1^7*t_2^8+98*t_0^5*t_1^6*t_2^9+60*t_0^5*t_1^5*t_2^10+27*t_0^5*t_1^4*t_2^11+8*t_0^5*t_1^3*t_2^12+t_0^5*t_1^2*t_2^13+2*t_0^4*t_1^13*t_2^3+10*t_0^4*t_1^12*t_2^4+27*t_0^4*t_1^11*t_2^5+52*t_0^4*t_1^10*t_2^6+75*t_0^4*t_1^9*t_2^7+84*t_0^4*t_1^8*t_2^8+75*t_0^4*t_1^7*t_2^9+52*t_0^4*t_1^6*t_2^10+27*t_0^4*t_1^5*t_2^11+10*t_0^4*t_1^4*t_2^12+2*t_0^4*t_1^3*t_2^13+2*t_0^3*t_1^13*t_2^4+8*t_0^3*t_1^12*t_2^5+19*t_0^3*t_1^11*t_2^6+33*t_0^3*t_1^10*t_2^7+42*t_0^3*t_1^9*t_2^8+42*t_0^3*t_1^8*t_2^9+33*t_0^3*t_1^7*t_2^10+19*t_0^3*t_1^6*t_2^11+8*t_0^3*t_1^5*t_2^12+2*t_0^3*t_1^4*t_2^13+t_0^2*t_1^13*t_2^5+4*t_0^2*t_1^12*t_2^6+9*t_0^2*t_1^11*t_2^7+14*t_0^2*t_1^10*t_2^8+16*t_0^2*t_1^9*t_2^9+14*t_0^2*t_1^8*t_2^10+9*t_0^2*t_1^7*t_2^11+4*t_0^2*t_1^6*t_2^12+t_0^2*t_1^5*t_2^13+t_0*t_1^12*t_2^7+2*t_0*t_1^11*t_2^8+3*t_0*t_1^10*t_2^9+3*t_0*t_1^9*t_2^10+2*t_0*t_1^8*t_2^11+t_0*t_1^7*t_2^12, (3,2) => 0, (6,0) => 0, (5,1) => t_0^15*t_1^5*t_2^4+t_0^15*t_1^4*t_2^5+t_0^14*t_1^8*t_2^2+3*t_0^14*t_1^7*t_2^3+6*t_0^14*t_1^6*t_2^4+7*t_0^14*t_1^5*t_2^5+6*t_0^14*t_1^4*t_2^6+3*t_0^14*t_1^3*t_2^7+t_0^14*t_1^2*t_2^8+2*t_0^13*t_1^9*t_2^2+8*t_0^13*t_1^8*t_2^3+17*t_0^13*t_1^7*t_2^4+23*t_0^13*t_1^6*t_2^5+23*t_0^13*t_1^5*t_2^6+17*t_0^13*t_1^4*t_2^7+8*t_0^13*t_1^3*t_2^8+2*t_0^13*t_1^2*t_2^9+4*t_0^12*t_1^10*t_2^2+14*t_0^12*t_1^9*t_2^3+34*t_0^12*t_1^8*t_2^4+52*t_0^12*t_1^7*t_2^5+60*t_0^12*t_1^6*t_2^6+52*t_0^12*t_1^5*t_2^7+34*t_0^12*t_1^4*t_2^8+14*t_0^12*t_1^3*t_2^9+4*t_0^12*t_1^2*t_2^10+4*t_0^11*t_1^11*t_2^2+19*t_0^11*t_1^10*t_2^3+48*t_0^11*t_1^9*t_2^4+84*t_0^11*t_1^8*t_2^5+109*t_0^11*t_1^7*t_2^6+109*t_0^11*t_1^6*t_2^7+84*t_0^11*t_1^5*t_2^8+48*t_0^11*t_1^4*t_2^9+19*t_0^11*t_1^3*t_2^10+4*t_0^11*t_1^2*t_2^11+4*t_0^10*t_1^12*t_2^2+19*t_0^10*t_1^11*t_2^3+56*t_0^10*t_1^10*t_2^4+106*t_0^10*t_1^9*t_2^5+156*t_0^10*t_1^8*t_2^6+175*t_0^10*t_1^7*t_2^7+156*t_0^10*t_1^6*t_2^8+106*t_0^10*t_1^5*t_2^9+56*t_0^10*t_1^4*t_2^10+19*t_0^10*t_1^3*t_2^11+4*t_0^10*t_1^2*t_2^12+2*t_0^9*t_1^13*t_2^2+14*t_0^9*t_1^12*t_2^3+48*t_0^9*t_1^11*t_2^4+106*t_0^9*t_1^10*t_2^5+171*t_0^9*t_1^9*t_2^6+218*t_0^9*t_1^8*t_2^7+218*t_0^9*t_1^7*t_2^8+171*t_0^9*t_1^6*t_2^9+106*t_0^9*t_1^5*t_2^10+48*t_0^9*t_1^4*t_2^11+14*t_0^9*t_1^3*t_2^12+2*t_0^9*t_1^2*t_2^13+t_0^8*t_1^14*t_2^2+8*t_0^8*t_1^13*t_2^3+34*t_0^8*t_1^12*t_2^4+84*t_0^8*t_1^11*t_2^5+156*t_0^8*t_1^10*t_2^6+218*t_0^8*t_1^9*t_2^7+246*t_0^8*t_1^8*t_2^8+218*t_0^8*t_1^7*t_2^9+156*t_0^8*t_1^6*t_2^10+84*t_0^8*t_1^5*t_2^11+34*t_0^8*t_1^4*t_2^12+8*t_0^8*t_1^3*t_2^13+t_0^8*t_1^2*t_2^14+3*t_0^7*t_1^14*t_2^3+17*t_0^7*t_1^13*t_2^4+52*t_0^7*t_1^12*t_2^5+109*t_0^7*t_1^11*t_2^6+175*t_0^7*t_1^10*t_2^7+218*t_0^7*t_1^9*t_2^8+218*t_0^7*t_1^8*t_2^9+175*t_0^7*t_1^7*t_2^10+109*t_0^7*t_1^6*t_2^11+52*t_0^7*t_1^5*t_2^12+17*t_0^7*t_1^4*t_2^13+3*t_0^7*t_1^3*t_2^14+6*t_0^6*t_1^14*t_2^4+23*t_0^6*t_1^13*t_2^5+60*t_0^6*t_1^12*t_2^6+109*t_0^6*t_1^11*t_2^7+156*t_0^6*t_1^10*t_2^8+171*t_0^6*t_1^9*t_2^9+156*t_0^6*t_1^8*t_2^10+109*t_0^6*t_1^7*t_2^11+60*t_0^6*t_1^6*t_2^12+23*t_0^6*t_1^5*t_2^13+6*t_0^6*t_1^4*t_2^14+t_0^5*t_1^15*t_2^4+7*t_0^5*t_1^14*t_2^5+23*t_0^5*t_1^13*t_2^6+52*t_0^5*t_1^12*t_2^7+84*t_0^5*t_1^11*t_2^8+106*t_0^5*t_1^10*t_2^9+106*t_0^5*t_1^9*t_2^10+84*t_0^5*t_1^8*t_2^11+52*t_0^5*t_1^7*t_2^12+23*t_0^5*t_1^6*t_2^13+7*t_0^5*t_1^5*t_2^14+t_0^5*t_1^4*t_2^15+t_0^4*t_1^15*t_2^5+6*t_0^4*t_1^14*t_2^6+17*t_0^4*t_1^13*t_2^7+34*t_0^4*t_1^12*t_2^8+48*t_0^4*t_1^11*t_2^9+56*t_0^4*t_1^10*t_2^10+48*t_0^4*t_1^9*t_2^11+34*t_0^4*t_1^8*t_2^12+17*t_0^4*t_1^7*t_2^13+6*t_0^4*t_1^6*t_2^14+t_0^4*t_1^5*t_2^15+3*t_0^3*t_1^14*t_2^7+8*t_0^3*t_1^13*t_2^8+14*t_0^3*t_1^12*t_2^9+19*t_0^3*t_1^11*t_2^10+19*t_0^3*t_1^10*t_2^11+14*t_0^3*t_1^9*t_2^12+8*t_0^3*t_1^8*t_2^13+3*t_0^3*t_1^7*t_2^14+t_0^2*t_1^14*t_2^8+2*t_0^2*t_1^13*t_2^9+4*t_0^2*t_1^12*t_2^10+4*t_0^2*t_1^11*t_2^11+4*t_0^2*t_1^10*t_2^12+2*t_0^2*t_1^9*t_2^13+t_0^2*t_1^8*t_2^14, (4,2) => 0};
--sb represents the betti numbers as sums of Schur functors
sb40 = new HashTable from {(7,0) => {}, (6,1) => {({16,8,4},1/1),({16,7,5},1/1),({15,10,3},1/1),({15,9,4},2/1),({15,8,5},3/1),({15,7,6},2/1),({14,11,3},1/1),({14,10,4},3/1),({14,9,5},5/1),({14,8,6},5/1),({14,7,7},2/1),({13,12,3},1/1),({13,11,4},3/1),({13,10,5},6/1),({13,9,6},7/1),({13,8,7},6/1),({12,12,4},1/1),({12,11,5},4/1),({12,10,6},7/1),({12,9,7},7/1),({12,8,8},3/1),({11,11,6},3/1),({11,10,7},5/1),({11,9,8},4/1),({10,10,8},2/1),({10,9,9},1/1)}, (5,2) => {}, (8,0) => {}, (6,2) => {}, (7,1) => {({17,10,5},1/1),({17,9,6},1/1),({17,8,7},1/1),({16,11,5},2/1),({16,10,6},2/1),({16,9,7},3/1),({16,8,8},1/1),({15,13,4},1/1),({15,12,5},2/1),({15,11,6},4/1),({15,10,7},5/1),({15,9,8},4/1),({14,13,5},1/1),({14,12,6},3/1),({14,11,7},6/1),({14,10,8},5/1),({14,9,9},3/1),({13,13,6},3/1),({13,12,7},4/1),({13,11,8},6/1),({13,10,9},4/1),({12,12,8},1/1),({12,11,9},4/1),({12,10,10},1/1),({11,11,10},1/1)}, (7,2) => {}, (8,1) => {({18,11,7},1/1),({18,10,8},1/1),({17,13,6},1/1),({17,12,7},1/1),({17,11,8},2/1),({17,10,9},2/1),({16,14,6},1/1),({16,13,7},2/1),({16,12,8},3/1),({16,11,9},3/1),({16,10,10},1/1),({15,14,7},1/1),({15,13,8},3/1),({15,12,9},3/1),({15,11,10},3/1),({14,14,8},1/1),({14,13,9},3/1),({14,12,10},3/1),({14,11,11},1/1),({13,13,10},1/1),({13,12,11},1/1)}, (9,0) => {}, (8,2) => {}, (9,1) => {({19,11,10},1/1),({18,14,8},1/1),({18,13,9},1/1),({18,12,10},1/1),({17,14,9},1/1),({17,13,10},1/1),({17,12,11},1/1),({16,16,8},1/1),({16,14,10},2/1),({16,13,11},1/1),({16,12,12},1/1),({15,14,11},1/1),({14,14,12},1/1)}, (10,0) => {}, (11,0) => {}, (10,1) => {({19,14,11},1/1)}, (9,2) => {}, (12,0) => {}, (11,1) => {}, (10,2) => {({18,18,12},1/1),({18,16,14},1/1)}, (12,1) => {}, (11,2) => {({19,18,15},1/1)}, (12,2) => {({19,19,18},1/1)}, (0,0) => {({0,0,0},1/1)}, (0,1) => {}, (1,0) => {}, (1,1) => {({6,2,0},1/1),({4,4,0},1/1)}, (0,2) => {}, (2,0) => {}, (1,2) => {}, (2,1) => {({9,2,1},1/1),({8,4,0},1/1),({8,3,1},1/1),({7,5,0},1/1),({7,4,1},1/1),({7,3,2},1/1),({6,5,1},1/1),({6,4,2},1/1),({5,4,3},1/1)}, (3,0) => {}, (4,0) => {}, (2,2) => {}, (3,1) => {({11,4,1},1/1),({11,3,2},1/1),({10,5,1},2/1),({10,4,2},1/1),({10,3,3},1/1),({9,7,0},1/1),({9,6,1},2/1),({9,5,2},3/1),({9,4,3},2/1),({8,7,1},1/1),({8,6,2},2/1),({8,5,3},3/1),({8,4,4},1/1),({7,7,2},2/1),({7,6,3},2/1),({7,5,4},2/1),({6,5,5},1/1)}, (5,0) => {}, (4,1) => {({13,5,2},1/1),({13,4,3},1/1),({12,7,1},1/1),({12,6,2},2/1),({12,5,3},2/1),({12,4,4},1/1),({11,8,1},1/1),({11,7,2},3/1),({11,6,3},4/1),({11,5,4},3/1),({10,9,1},1/1),({10,8,2},3/1),({10,7,3},5/1),({10,6,4},5/1),({10,5,5},2/1),({9,9,2},1/1),({9,8,3},3/1),({9,7,4},5/1),({9,6,5},4/1),({8,8,4},2/1),({8,7,5},3/1),({8,6,6},1/1),({7,7,6},1/1)}, (3,2) => {}, (6,0) => {}, (5,1) => {({15,5,4},1/1),({14,8,2},1/1),({14,7,3},2/1),({14,6,4},2/1),({13,9,2},1/1),({13,8,3},3/1),({13,7,4},4/1),({13,6,5},3/1),({12,10,2},2/1),({12,9,3},3/1),({12,8,4},7/1),({12,7,5},6/1),({12,6,6},3/1),({11,10,3},3/1),({11,9,4},5/1),({11,8,5},7/1),({11,7,6},5/1),({10,10,4},3/1),({10,9,5},5/1),({10,8,6},7/1),({10,7,7},2/1),({9,9,6},1/1),({9,8,7},3/1)}, (4,2) => {}};
--dw encodes the dominant weights in each entry
dw40 = new HashTable from {(5,2) => {}, (6,1) => {{16,8,4},{15,10,3}}, (7,0) => {}, (8,0) => {}, (6,2) => {}, (7,1) => {{17,10,5},{15,13,4}}, (9,0) => {}, (8,1) => {{18,11,7},{17,13,6}}, (7,2) => {}, (10,0) => {}, (9,1) => {{19,11,10},{18,14,8}}, (8,2) => {}, (11,0) => {}, (10,1) => {{19,14,11}}, (9,2) => {}, (12,0) => {}, (11,1) => {}, (10,2) => {{18,18,12}}, (12,1) => {}, (11,2) => {{19,18,15}}, (12,2) => {{19,19,18}}, (0,0) => {{0,0,0}}, (0,1) => {}, (1,0) => {}, (2,0) => {}, (0,2) => {}, (1,1) => {{6,2,0}}, (3,0) => {}, (2,1) => {{9,2,1},{8,4,0}}, (1,2) => {}, (3,1) => {{11,4,1},{9,7,0}}, (2,2) => {}, (4,0) => {}, (3,2) => {}, (4,1) => {{13,5,2},{12,7,1}}, (5,0) => {}, (4,2) => {}, (5,1) => {{15,5,4},{14,8,2}}, (6,0) => {}};
--lw encodes the lex leading weight in each entry
lw40 = new HashTable from {(5,2) => {}, (6,1) => {16,8,4}, (7,0) => {}, (8,0) => {}, (6,2) => {}, (7,1) => {17,10,5}, (9,0) => {}, (8,1) => {18,11,7}, (7,2) => {}, (10,0) => {}, (9,1) => {19,11,10}, (8,2) => {}, (11,0) => {}, (10,1) => {19,14,11}, (9,2) => {}, (12,0) => {}, (11,1) => {}, (10,2) => {18,18,12}, (12,1) => {}, (11,2) => {19,18,15}, (12,2) => {19,19,18}, (0,0) => {0,0,0}, (0,1) => {}, (1,0) => {}, (2,0) => {}, (0,2) => {}, (1,1) => {6,2,0}, (3,0) => {}, (2,1) => {9,2,1}, (1,2) => {}, (3,1) => {11,4,1}, (2,2) => {}, (4,0) => {}, (3,2) => {}, (4,1) => {13,5,2}, (5,0) => {}, (4,2) => {}, (5,1) => {15,5,4}, (6,0) => {}};
--nr encodes the number of disctinct reprsentations in each entry
nr40 = new HashTable from {(5,2) => 0, (6,1) => 26, (7,0) => 0, (8,0) => 0, (6,2) => 0, (7,1) => 25, (9,0) => 0, (8,1) => 21, (7,2) => 0, (10,0) => 0, (9,1) => 13, (8,2) => 0, (11,0) => 0, (10,1) => 1, (9,2) => 0, (12,0) => 0, (11,1) => 0, (10,2) => 2, (12,1) => 0, (11,2) => 1, (12,2) => 1, (0,0) => 1, (0,1) => 0, (1,0) => 0, (2,0) => 0, (0,2) => 0, (1,1) => 2, (3,0) => 0, (2,1) => 9, (1,2) => 0, (3,1) => 17, (2,2) => 0, (4,0) => 0, (3,2) => 0, (4,1) => 23, (5,0) => 0, (4,2) => 0, (5,1) => 23, (6,0) => 0};
--nrm encodes the number of representations with multiplicity in each entry
nrm40 = new HashTable from {(5,2) => 0, (6,1) => 86/1, (7,0) => 0, (8,0) => 0, (6,2) => 0, (7,1) => 69/1, (9,0) => 0, (8,1) => 38/1, (7,2) => 0, (10,0) => 0, (9,1) => 14/1, (8,2) => 0, (11,0) => 0, (10,1) => 1/1, (9,2) => 0, (12,0) => 0, (11,1) => 0, (10,2) => 2/1, (12,1) => 0, (11,2) => 1/1, (12,2) => 1/1, (0,0) => 1/1, (0,1) => 0, (1,0) => 0, (2,0) => 0, (0,2) => 0, (1,1) => 2/1, (3,0) => 0, (2,1) => 9/1, (1,2) => 0, (3,1) => 28/1, (2,2) => 0, (4,0) => 0, (3,2) => 0, (4,1) => 55/1, (5,0) => 0, (4,2) => 0, (5,1) => 79/1, (6,0) => 0};
--er encodes the errors in the computed multigraded Hilbert series via our Schur method in each entry
er40 = new HashTable from {(5,2) => 0, (6,1) => 0, (7,0) => 0, (8,0) => 0, (6,2) => 0, (7,1) => 0, (9,0) => 0, (8,1) => 0, (7,2) => 0, (10,0) => 0, (9,1) => 0, (8,2) => 0, (11,0) => 0, (10,1) => 0, (9,2) => 0, (12,0) => 0, (11,1) => 0, (10,2) => 0, (12,1) => 0, (11,2) => 0, (12,2) => 0, (0,0) => 0, (0,1) => 0, (1,0) => 0, (2,0) => 0, (0,2) => 0, (1,1) => 0, (3,0) => 0, (2,1) => 0, (1,2) => 0, (3,1) => 0, (2,2) => 0, (4,0) => 0, (3,2) => 0, (4,1) => 0, (5,0) => 0, (4,2) => 0, (5,1) => 0, (6,0) => 0};
--bs encodes the Boij-Soederberg coefficients each entry
bs40 = {2874009600/1,4790016000/1};
end;