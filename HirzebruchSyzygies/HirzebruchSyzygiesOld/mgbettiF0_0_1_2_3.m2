A := QQ[t_0,t_1,t_2,t_3];
--mb stands for Multigraded Betti numbers
mb0123 = new HashTable from {(7,0) => 0, (6,1) => t_0^10*t_1^4*t_2^13*t_3^9+2*t_0^10*t_1^4*t_2^12*t_3^10+3*t_0^10*t_1^4*t_2^11*t_3^11+2*t_0^10*t_1^4*t_2^10*t_3^12+t_0^10*t_1^4*t_2^9*t_3^13+2*t_0^9*t_1^5*t_2^14*t_3^8+6*t_0^9*t_1^5*t_2^13*t_3^9+11*t_0^9*t_1^5*t_2^12*t_3^10+13*t_0^9*t_1^5*t_2^11*t_3^11+11*t_0^9*t_1^5*t_2^10*t_3^12+6*t_0^9*t_1^5*t_2^9*t_3^13+2*t_0^9*t_1^5*t_2^8*t_3^14+2*t_0^8*t_1^6*t_2^15*t_3^7+7*t_0^8*t_1^6*t_2^14*t_3^8+17*t_0^8*t_1^6*t_2^13*t_3^9+26*t_0^8*t_1^6*t_2^12*t_3^10+31*t_0^8*t_1^6*t_2^11*t_3^11+26*t_0^8*t_1^6*t_2^10*t_3^12+17*t_0^8*t_1^6*t_2^9*t_3^13+7*t_0^8*t_1^6*t_2^8*t_3^14+2*t_0^8*t_1^6*t_2^7*t_3^15+2*t_0^7*t_1^7*t_2^15*t_3^7+9*t_0^7*t_1^7*t_2^14*t_3^8+21*t_0^7*t_1^7*t_2^13*t_3^9+33*t_0^7*t_1^7*t_2^12*t_3^10+38*t_0^7*t_1^7*t_2^11*t_3^11+33*t_0^7*t_1^7*t_2^10*t_3^12+21*t_0^7*t_1^7*t_2^9*t_3^13+9*t_0^7*t_1^7*t_2^8*t_3^14+2*t_0^7*t_1^7*t_2^7*t_3^15+2*t_0^6*t_1^8*t_2^15*t_3^7+7*t_0^6*t_1^8*t_2^14*t_3^8+17*t_0^6*t_1^8*t_2^13*t_3^9+26*t_0^6*t_1^8*t_2^12*t_3^10+31*t_0^6*t_1^8*t_2^11*t_3^11+26*t_0^6*t_1^8*t_2^10*t_3^12+17*t_0^6*t_1^8*t_2^9*t_3^13+7*t_0^6*t_1^8*t_2^8*t_3^14+2*t_0^6*t_1^8*t_2^7*t_3^15+2*t_0^5*t_1^9*t_2^14*t_3^8+6*t_0^5*t_1^9*t_2^13*t_3^9+11*t_0^5*t_1^9*t_2^12*t_3^10+13*t_0^5*t_1^9*t_2^11*t_3^11+11*t_0^5*t_1^9*t_2^10*t_3^12+6*t_0^5*t_1^9*t_2^9*t_3^13+2*t_0^5*t_1^9*t_2^8*t_3^14+t_0^4*t_1^10*t_2^13*t_3^9+2*t_0^4*t_1^10*t_2^12*t_3^10+3*t_0^4*t_1^10*t_2^11*t_3^11+2*t_0^4*t_1^10*t_2^10*t_3^12+t_0^4*t_1^10*t_2^9*t_3^13, (5,2) => 0, (8,0) => 0, (7,1) => t_0^11*t_1^5*t_2^13*t_3^12+t_0^11*t_1^5*t_2^12*t_3^13+t_0^10*t_1^6*t_2^15*t_3^10+3*t_0^10*t_1^6*t_2^14*t_3^11+5*t_0^10*t_1^6*t_2^13*t_3^12+5*t_0^10*t_1^6*t_2^12*t_3^13+3*t_0^10*t_1^6*t_2^11*t_3^14+t_0^10*t_1^6*t_2^10*t_3^15+t_0^9*t_1^7*t_2^16*t_3^9+4*t_0^9*t_1^7*t_2^15*t_3^10+9*t_0^9*t_1^7*t_2^14*t_3^11+13*t_0^9*t_1^7*t_2^13*t_3^12+13*t_0^9*t_1^7*t_2^12*t_3^13+9*t_0^9*t_1^7*t_2^11*t_3^14+4*t_0^9*t_1^7*t_2^10*t_3^15+t_0^9*t_1^7*t_2^9*t_3^16+2*t_0^8*t_1^8*t_2^16*t_3^9+6*t_0^8*t_1^8*t_2^15*t_3^10+12*t_0^8*t_1^8*t_2^14*t_3^11+17*t_0^8*t_1^8*t_2^13*t_3^12+17*t_0^8*t_1^8*t_2^12*t_3^13+12*t_0^8*t_1^8*t_2^11*t_3^14+6*t_0^8*t_1^8*t_2^10*t_3^15+2*t_0^8*t_1^8*t_2^9*t_3^16+t_0^7*t_1^9*t_2^16*t_3^9+4*t_0^7*t_1^9*t_2^15*t_3^10+9*t_0^7*t_1^9*t_2^14*t_3^11+13*t_0^7*t_1^9*t_2^13*t_3^12+13*t_0^7*t_1^9*t_2^12*t_3^13+9*t_0^7*t_1^9*t_2^11*t_3^14+4*t_0^7*t_1^9*t_2^10*t_3^15+t_0^7*t_1^9*t_2^9*t_3^16+t_0^6*t_1^10*t_2^15*t_3^10+3*t_0^6*t_1^10*t_2^14*t_3^11+5*t_0^6*t_1^10*t_2^13*t_3^12+5*t_0^6*t_1^10*t_2^12*t_3^13+3*t_0^6*t_1^10*t_2^11*t_3^14+t_0^6*t_1^10*t_2^10*t_3^15+t_0^5*t_1^11*t_2^13*t_3^12+t_0^5*t_1^11*t_2^12*t_3^13, (6,2) => 0, (7,2) => 0, (9,0) => 0, (8,1) => t_0^11*t_1^7*t_2^15*t_3^13+t_0^11*t_1^7*t_2^14*t_3^14+t_0^11*t_1^7*t_2^13*t_3^15+t_0^10*t_1^8*t_2^16*t_3^12+2*t_0^10*t_1^8*t_2^15*t_3^13+3*t_0^10*t_1^8*t_2^14*t_3^14+2*t_0^10*t_1^8*t_2^13*t_3^15+t_0^10*t_1^8*t_2^12*t_3^16+t_0^9*t_1^9*t_2^17*t_3^11+2*t_0^9*t_1^9*t_2^16*t_3^12+4*t_0^9*t_1^9*t_2^15*t_3^13+5*t_0^9*t_1^9*t_2^14*t_3^14+4*t_0^9*t_1^9*t_2^13*t_3^15+2*t_0^9*t_1^9*t_2^12*t_3^16+t_0^9*t_1^9*t_2^11*t_3^17+t_0^8*t_1^10*t_2^16*t_3^12+2*t_0^8*t_1^10*t_2^15*t_3^13+3*t_0^8*t_1^10*t_2^14*t_3^14+2*t_0^8*t_1^10*t_2^13*t_3^15+t_0^8*t_1^10*t_2^12*t_3^16+t_0^7*t_1^11*t_2^15*t_3^13+t_0^7*t_1^11*t_2^14*t_3^14+t_0^7*t_1^11*t_2^13*t_3^15, (8,2) => 0, (9,1) => 0, (9,2) => t_0^11*t_1^11*t_2^17*t_3^17, (10,2) => 0, (0,0) => t_2+t_3, (1,0) => t_0^2*t_2^3*t_3+t_0^2*t_2^2*t_3^2+t_0^2*t_2*t_3^3+t_0*t_1*t_2^3*t_3+t_0*t_1*t_2^2*t_3^2+t_0*t_1*t_2*t_3^3+t_1^2*t_2^3*t_3+t_1^2*t_2^2*t_3^2+t_1^2*t_2*t_3^3, (1,1) => t_0^2*t_1^2*t_2^7+t_0^2*t_1^2*t_2^6*t_3+t_0^2*t_1^2*t_2^5*t_3^2+t_0^2*t_1^2*t_2^4*t_3^3+t_0^2*t_1^2*t_2^3*t_3^4+t_0^2*t_1^2*t_2^2*t_3^5+t_0^2*t_1^2*t_2*t_3^6+t_0^2*t_1^2*t_3^7, (2,0) => 0, (2,1) => t_0^6*t_2^5*t_3^5+t_0^5*t_1*t_2^7*t_3^3+2*t_0^5*t_1*t_2^6*t_3^4+3*t_0^5*t_1*t_2^5*t_3^5+2*t_0^5*t_1*t_2^4*t_3^6+t_0^5*t_1*t_2^3*t_3^7+t_0^4*t_1^2*t_2^9*t_3+2*t_0^4*t_1^2*t_2^8*t_3^2+5*t_0^4*t_1^2*t_2^7*t_3^3+7*t_0^4*t_1^2*t_2^6*t_3^4+9*t_0^4*t_1^2*t_2^5*t_3^5+7*t_0^4*t_1^2*t_2^4*t_3^6+5*t_0^4*t_1^2*t_2^3*t_3^7+2*t_0^4*t_1^2*t_2^2*t_3^8+t_0^4*t_1^2*t_2*t_3^9+t_0^3*t_1^3*t_2^9*t_3+3*t_0^3*t_1^3*t_2^8*t_3^2+6*t_0^3*t_1^3*t_2^7*t_3^3+9*t_0^3*t_1^3*t_2^6*t_3^4+11*t_0^3*t_1^3*t_2^5*t_3^5+9*t_0^3*t_1^3*t_2^4*t_3^6+6*t_0^3*t_1^3*t_2^3*t_3^7+3*t_0^3*t_1^3*t_2^2*t_3^8+t_0^3*t_1^3*t_2*t_3^9+t_0^2*t_1^4*t_2^9*t_3+2*t_0^2*t_1^4*t_2^8*t_3^2+5*t_0^2*t_1^4*t_2^7*t_3^3+7*t_0^2*t_1^4*t_2^6*t_3^4+9*t_0^2*t_1^4*t_2^5*t_3^5+7*t_0^2*t_1^4*t_2^4*t_3^6+5*t_0^2*t_1^4*t_2^3*t_3^7+2*t_0^2*t_1^4*t_2^2*t_3^8+t_0^2*t_1^4*t_2*t_3^9+t_0*t_1^5*t_2^7*t_3^3+2*t_0*t_1^5*t_2^6*t_3^4+3*t_0*t_1^5*t_2^5*t_3^5+2*t_0*t_1^5*t_2^4*t_3^6+t_0*t_1^5*t_2^3*t_3^7+t_1^6*t_2^5*t_3^5, (3,0) => 0, (4,0) => 0, (3,1) => t_0^7*t_1*t_2^8*t_3^5+2*t_0^7*t_1*t_2^7*t_3^6+2*t_0^7*t_1*t_2^6*t_3^7+t_0^7*t_1*t_2^5*t_3^8+t_0^6*t_1^2*t_2^10*t_3^3+3*t_0^6*t_1^2*t_2^9*t_3^4+7*t_0^6*t_1^2*t_2^8*t_3^5+10*t_0^6*t_1^2*t_2^7*t_3^6+10*t_0^6*t_1^2*t_2^6*t_3^7+7*t_0^6*t_1^2*t_2^5*t_3^8+3*t_0^6*t_1^2*t_2^4*t_3^9+t_0^6*t_1^2*t_2^3*t_3^10+t_0^5*t_1^3*t_2^11*t_3^2+4*t_0^5*t_1^3*t_2^10*t_3^3+10*t_0^5*t_1^3*t_2^9*t_3^4+18*t_0^5*t_1^3*t_2^8*t_3^5+24*t_0^5*t_1^3*t_2^7*t_3^6+24*t_0^5*t_1^3*t_2^6*t_3^7+18*t_0^5*t_1^3*t_2^5*t_3^8+10*t_0^5*t_1^3*t_2^4*t_3^9+4*t_0^5*t_1^3*t_2^3*t_3^10+t_0^5*t_1^3*t_2^2*t_3^11+t_0^4*t_1^4*t_2^11*t_3^2+5*t_0^4*t_1^4*t_2^10*t_3^3+13*t_0^4*t_1^4*t_2^9*t_3^4+23*t_0^4*t_1^4*t_2^8*t_3^5+30*t_0^4*t_1^4*t_2^7*t_3^6+30*t_0^4*t_1^4*t_2^6*t_3^7+23*t_0^4*t_1^4*t_2^5*t_3^8+13*t_0^4*t_1^4*t_2^4*t_3^9+5*t_0^4*t_1^4*t_2^3*t_3^10+t_0^4*t_1^4*t_2^2*t_3^11+t_0^3*t_1^5*t_2^11*t_3^2+4*t_0^3*t_1^5*t_2^10*t_3^3+10*t_0^3*t_1^5*t_2^9*t_3^4+18*t_0^3*t_1^5*t_2^8*t_3^5+24*t_0^3*t_1^5*t_2^7*t_3^6+24*t_0^3*t_1^5*t_2^6*t_3^7+18*t_0^3*t_1^5*t_2^5*t_3^8+10*t_0^3*t_1^5*t_2^4*t_3^9+4*t_0^3*t_1^5*t_2^3*t_3^10+t_0^3*t_1^5*t_2^2*t_3^11+t_0^2*t_1^6*t_2^10*t_3^3+3*t_0^2*t_1^6*t_2^9*t_3^4+7*t_0^2*t_1^6*t_2^8*t_3^5+10*t_0^2*t_1^6*t_2^7*t_3^6+10*t_0^2*t_1^6*t_2^6*t_3^7+7*t_0^2*t_1^6*t_2^5*t_3^8+3*t_0^2*t_1^6*t_2^4*t_3^9+t_0^2*t_1^6*t_2^3*t_3^10+t_0*t_1^7*t_2^8*t_3^5+2*t_0*t_1^7*t_2^7*t_3^6+2*t_0*t_1^7*t_2^6*t_3^7+t_0*t_1^7*t_2^5*t_3^8, (2,2) => 0, (5,0) => 0, (4,1) => 2*t_0^8*t_1^2*t_2^10*t_3^6+3*t_0^8*t_1^2*t_2^9*t_3^7+4*t_0^8*t_1^2*t_2^8*t_3^8+3*t_0^8*t_1^2*t_2^7*t_3^9+2*t_0^8*t_1^2*t_2^6*t_3^10+t_0^7*t_1^3*t_2^12*t_3^4+4*t_0^7*t_1^3*t_2^11*t_3^5+10*t_0^7*t_1^3*t_2^10*t_3^6+16*t_0^7*t_1^3*t_2^9*t_3^7+19*t_0^7*t_1^3*t_2^8*t_3^8+16*t_0^7*t_1^3*t_2^7*t_3^9+10*t_0^7*t_1^3*t_2^6*t_3^10+4*t_0^7*t_1^3*t_2^5*t_3^11+t_0^7*t_1^3*t_2^4*t_3^12+3*t_0^6*t_1^4*t_2^12*t_3^4+10*t_0^6*t_1^4*t_2^11*t_3^5+23*t_0^6*t_1^4*t_2^10*t_3^6+34*t_0^6*t_1^4*t_2^9*t_3^7+40*t_0^6*t_1^4*t_2^8*t_3^8+34*t_0^6*t_1^4*t_2^7*t_3^9+23*t_0^6*t_1^4*t_2^6*t_3^10+10*t_0^6*t_1^4*t_2^5*t_3^11+3*t_0^6*t_1^4*t_2^4*t_3^12+t_0^5*t_1^5*t_2^13*t_3^3+5*t_0^5*t_1^5*t_2^12*t_3^4+15*t_0^5*t_1^5*t_2^11*t_3^5+31*t_0^5*t_1^5*t_2^10*t_3^6+46*t_0^5*t_1^5*t_2^9*t_3^7+52*t_0^5*t_1^5*t_2^8*t_3^8+46*t_0^5*t_1^5*t_2^7*t_3^9+31*t_0^5*t_1^5*t_2^6*t_3^10+15*t_0^5*t_1^5*t_2^5*t_3^11+5*t_0^5*t_1^5*t_2^4*t_3^12+t_0^5*t_1^5*t_2^3*t_3^13+3*t_0^4*t_1^6*t_2^12*t_3^4+10*t_0^4*t_1^6*t_2^11*t_3^5+23*t_0^4*t_1^6*t_2^10*t_3^6+34*t_0^4*t_1^6*t_2^9*t_3^7+40*t_0^4*t_1^6*t_2^8*t_3^8+34*t_0^4*t_1^6*t_2^7*t_3^9+23*t_0^4*t_1^6*t_2^6*t_3^10+10*t_0^4*t_1^6*t_2^5*t_3^11+3*t_0^4*t_1^6*t_2^4*t_3^12+t_0^3*t_1^7*t_2^12*t_3^4+4*t_0^3*t_1^7*t_2^11*t_3^5+10*t_0^3*t_1^7*t_2^10*t_3^6+16*t_0^3*t_1^7*t_2^9*t_3^7+19*t_0^3*t_1^7*t_2^8*t_3^8+16*t_0^3*t_1^7*t_2^7*t_3^9+10*t_0^3*t_1^7*t_2^6*t_3^10+4*t_0^3*t_1^7*t_2^5*t_3^11+t_0^3*t_1^7*t_2^4*t_3^12+2*t_0^2*t_1^8*t_2^10*t_3^6+3*t_0^2*t_1^8*t_2^9*t_3^7+4*t_0^2*t_1^8*t_2^8*t_3^8+3*t_0^2*t_1^8*t_2^7*t_3^9+2*t_0^2*t_1^8*t_2^6*t_3^10, (3,2) => 0, (6,0) => 0, (5,1) => t_0^9*t_1^3*t_2^12*t_3^7+3*t_0^9*t_1^3*t_2^11*t_3^8+4*t_0^9*t_1^3*t_2^10*t_3^9+4*t_0^9*t_1^3*t_2^9*t_3^10+3*t_0^9*t_1^3*t_2^8*t_3^11+t_0^9*t_1^3*t_2^7*t_3^12+2*t_0^8*t_1^4*t_2^13*t_3^6+7*t_0^8*t_1^4*t_2^12*t_3^7+14*t_0^8*t_1^4*t_2^11*t_3^8+19*t_0^8*t_1^4*t_2^10*t_3^9+19*t_0^8*t_1^4*t_2^9*t_3^10+14*t_0^8*t_1^4*t_2^8*t_3^11+7*t_0^8*t_1^4*t_2^7*t_3^12+2*t_0^8*t_1^4*t_2^6*t_3^13+t_0^7*t_1^5*t_2^14*t_3^5+6*t_0^7*t_1^5*t_2^13*t_3^6+17*t_0^7*t_1^5*t_2^12*t_3^7+31*t_0^7*t_1^5*t_2^11*t_3^8+41*t_0^7*t_1^5*t_2^10*t_3^9+41*t_0^7*t_1^5*t_2^9*t_3^10+31*t_0^7*t_1^5*t_2^8*t_3^11+17*t_0^7*t_1^5*t_2^7*t_3^12+6*t_0^7*t_1^5*t_2^6*t_3^13+t_0^7*t_1^5*t_2^5*t_3^14+2*t_0^6*t_1^6*t_2^14*t_3^5+9*t_0^6*t_1^6*t_2^13*t_3^6+23*t_0^6*t_1^6*t_2^12*t_3^7+41*t_0^6*t_1^6*t_2^11*t_3^8+53*t_0^6*t_1^6*t_2^10*t_3^9+53*t_0^6*t_1^6*t_2^9*t_3^10+41*t_0^6*t_1^6*t_2^8*t_3^11+23*t_0^6*t_1^6*t_2^7*t_3^12+9*t_0^6*t_1^6*t_2^6*t_3^13+2*t_0^6*t_1^6*t_2^5*t_3^14+t_0^5*t_1^7*t_2^14*t_3^5+6*t_0^5*t_1^7*t_2^13*t_3^6+17*t_0^5*t_1^7*t_2^12*t_3^7+31*t_0^5*t_1^7*t_2^11*t_3^8+41*t_0^5*t_1^7*t_2^10*t_3^9+41*t_0^5*t_1^7*t_2^9*t_3^10+31*t_0^5*t_1^7*t_2^8*t_3^11+17*t_0^5*t_1^7*t_2^7*t_3^12+6*t_0^5*t_1^7*t_2^6*t_3^13+t_0^5*t_1^7*t_2^5*t_3^14+2*t_0^4*t_1^8*t_2^13*t_3^6+7*t_0^4*t_1^8*t_2^12*t_3^7+14*t_0^4*t_1^8*t_2^11*t_3^8+19*t_0^4*t_1^8*t_2^10*t_3^9+19*t_0^4*t_1^8*t_2^9*t_3^10+14*t_0^4*t_1^8*t_2^8*t_3^11+7*t_0^4*t_1^8*t_2^7*t_3^12+2*t_0^4*t_1^8*t_2^6*t_3^13+t_0^3*t_1^9*t_2^12*t_3^7+3*t_0^3*t_1^9*t_2^11*t_3^8+4*t_0^3*t_1^9*t_2^10*t_3^9+4*t_0^3*t_1^9*t_2^9*t_3^10+3*t_0^3*t_1^9*t_2^8*t_3^11+t_0^3*t_1^9*t_2^7*t_3^12, (4,2) => 0};
end;