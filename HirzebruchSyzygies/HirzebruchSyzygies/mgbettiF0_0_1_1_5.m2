A := QQ[t_0,t_1,t_2,t_3];
--mb stands for Multigraded Betti numbers
mb0115 = new HashTable from {(7,0) => 0, (6,1) => t_0^5*t_1^2*t_2^23*t_3^13+2*t_0^5*t_1^2*t_2^22*t_3^14+4*t_0^5*t_1^2*t_2^21*t_3^15+6*t_0^5*t_1^2*t_2^20*t_3^16+8*t_0^5*t_1^2*t_2^19*t_3^17+8*t_0^5*t_1^2*t_2^18*t_3^18+8*t_0^5*t_1^2*t_2^17*t_3^19+6*t_0^5*t_1^2*t_2^16*t_3^20+4*t_0^5*t_1^2*t_2^15*t_3^21+2*t_0^5*t_1^2*t_2^14*t_3^22+t_0^5*t_1^2*t_2^13*t_3^23+t_0^4*t_1^3*t_2^25*t_3^11+3*t_0^4*t_1^3*t_2^24*t_3^12+7*t_0^4*t_1^3*t_2^23*t_3^13+13*t_0^4*t_1^3*t_2^22*t_3^14+21*t_0^4*t_1^3*t_2^21*t_3^15+28*t_0^4*t_1^3*t_2^20*t_3^16+34*t_0^4*t_1^3*t_2^19*t_3^17+36*t_0^4*t_1^3*t_2^18*t_3^18+34*t_0^4*t_1^3*t_2^17*t_3^19+28*t_0^4*t_1^3*t_2^16*t_3^20+21*t_0^4*t_1^3*t_2^15*t_3^21+13*t_0^4*t_1^3*t_2^14*t_3^22+7*t_0^4*t_1^3*t_2^13*t_3^23+3*t_0^4*t_1^3*t_2^12*t_3^24+t_0^4*t_1^3*t_2^11*t_3^25+t_0^3*t_1^4*t_2^25*t_3^11+3*t_0^3*t_1^4*t_2^24*t_3^12+7*t_0^3*t_1^4*t_2^23*t_3^13+13*t_0^3*t_1^4*t_2^22*t_3^14+21*t_0^3*t_1^4*t_2^21*t_3^15+28*t_0^3*t_1^4*t_2^20*t_3^16+34*t_0^3*t_1^4*t_2^19*t_3^17+36*t_0^3*t_1^4*t_2^18*t_3^18+34*t_0^3*t_1^4*t_2^17*t_3^19+28*t_0^3*t_1^4*t_2^16*t_3^20+21*t_0^3*t_1^4*t_2^15*t_3^21+13*t_0^3*t_1^4*t_2^14*t_3^22+7*t_0^3*t_1^4*t_2^13*t_3^23+3*t_0^3*t_1^4*t_2^12*t_3^24+t_0^3*t_1^4*t_2^11*t_3^25+t_0^2*t_1^5*t_2^23*t_3^13+2*t_0^2*t_1^5*t_2^22*t_3^14+4*t_0^2*t_1^5*t_2^21*t_3^15+6*t_0^2*t_1^5*t_2^20*t_3^16+8*t_0^2*t_1^5*t_2^19*t_3^17+8*t_0^2*t_1^5*t_2^18*t_3^18+8*t_0^2*t_1^5*t_2^17*t_3^19+6*t_0^2*t_1^5*t_2^16*t_3^20+4*t_0^2*t_1^5*t_2^15*t_3^21+2*t_0^2*t_1^5*t_2^14*t_3^22+t_0^2*t_1^5*t_2^13*t_3^23, (5,2) => 0, (8,0) => 0, (7,1) => t_0^5*t_1^3*t_2^26*t_3^15+2*t_0^5*t_1^3*t_2^25*t_3^16+4*t_0^5*t_1^3*t_2^24*t_3^17+6*t_0^5*t_1^3*t_2^23*t_3^18+8*t_0^5*t_1^3*t_2^22*t_3^19+9*t_0^5*t_1^3*t_2^21*t_3^20+9*t_0^5*t_1^3*t_2^20*t_3^21+8*t_0^5*t_1^3*t_2^19*t_3^22+6*t_0^5*t_1^3*t_2^18*t_3^23+4*t_0^5*t_1^3*t_2^17*t_3^24+2*t_0^5*t_1^3*t_2^16*t_3^25+t_0^5*t_1^3*t_2^15*t_3^26+t_0^4*t_1^4*t_2^27*t_3^14+3*t_0^4*t_1^4*t_2^26*t_3^15+6*t_0^4*t_1^4*t_2^25*t_3^16+10*t_0^4*t_1^4*t_2^24*t_3^17+15*t_0^4*t_1^4*t_2^23*t_3^18+19*t_0^4*t_1^4*t_2^22*t_3^19+21*t_0^4*t_1^4*t_2^21*t_3^20+21*t_0^4*t_1^4*t_2^20*t_3^21+19*t_0^4*t_1^4*t_2^19*t_3^22+15*t_0^4*t_1^4*t_2^18*t_3^23+10*t_0^4*t_1^4*t_2^17*t_3^24+6*t_0^4*t_1^4*t_2^16*t_3^25+3*t_0^4*t_1^4*t_2^15*t_3^26+t_0^4*t_1^4*t_2^14*t_3^27+t_0^3*t_1^5*t_2^26*t_3^15+2*t_0^3*t_1^5*t_2^25*t_3^16+4*t_0^3*t_1^5*t_2^24*t_3^17+6*t_0^3*t_1^5*t_2^23*t_3^18+8*t_0^3*t_1^5*t_2^22*t_3^19+9*t_0^3*t_1^5*t_2^21*t_3^20+9*t_0^3*t_1^5*t_2^20*t_3^21+8*t_0^3*t_1^5*t_2^19*t_3^22+6*t_0^3*t_1^5*t_2^18*t_3^23+4*t_0^3*t_1^5*t_2^17*t_3^24+2*t_0^3*t_1^5*t_2^16*t_3^25+t_0^3*t_1^5*t_2^15*t_3^26, (6,2) => 0, (7,2) => 0, (9,0) => 0, (8,1) => t_0^5*t_1^4*t_2^28*t_3^18+2*t_0^5*t_1^4*t_2^27*t_3^19+3*t_0^5*t_1^4*t_2^26*t_3^20+4*t_0^5*t_1^4*t_2^25*t_3^21+5*t_0^5*t_1^4*t_2^24*t_3^22+5*t_0^5*t_1^4*t_2^23*t_3^23+5*t_0^5*t_1^4*t_2^22*t_3^24+4*t_0^5*t_1^4*t_2^21*t_3^25+3*t_0^5*t_1^4*t_2^20*t_3^26+2*t_0^5*t_1^4*t_2^19*t_3^27+t_0^5*t_1^4*t_2^18*t_3^28+t_0^4*t_1^5*t_2^28*t_3^18+2*t_0^4*t_1^5*t_2^27*t_3^19+3*t_0^4*t_1^5*t_2^26*t_3^20+4*t_0^4*t_1^5*t_2^25*t_3^21+5*t_0^4*t_1^5*t_2^24*t_3^22+5*t_0^4*t_1^5*t_2^23*t_3^23+5*t_0^4*t_1^5*t_2^22*t_3^24+4*t_0^4*t_1^5*t_2^21*t_3^25+3*t_0^4*t_1^5*t_2^20*t_3^26+2*t_0^4*t_1^5*t_2^19*t_3^27+t_0^4*t_1^5*t_2^18*t_3^28, (8,2) => 0, (9,1) => t_0^5*t_1^5*t_2^29*t_3^22+t_0^5*t_1^5*t_2^28*t_3^23+t_0^5*t_1^5*t_2^27*t_3^24+t_0^5*t_1^5*t_2^26*t_3^25+t_0^5*t_1^5*t_2^25*t_3^26+t_0^5*t_1^5*t_2^24*t_3^27+t_0^5*t_1^5*t_2^23*t_3^28+t_0^5*t_1^5*t_2^22*t_3^29, (10,1) => 0, (9,2) => 0, (11,1) => 0, (0,0) => t_2+t_3, (1,0) => t_0*t_2^5*t_3+t_0*t_2^4*t_3^2+t_0*t_2^3*t_3^3+t_0*t_2^2*t_3^4+t_0*t_2*t_3^5+t_1*t_2^5*t_3+t_1*t_2^4*t_3^2+t_1*t_2^3*t_3^3+t_1*t_2^2*t_3^4+t_1*t_2*t_3^5, (1,1) => 0, (2,0) => 0, (2,1) => t_0^3*t_2^11*t_3^5+t_0^3*t_2^10*t_3^6+2*t_0^3*t_2^9*t_3^7+2*t_0^3*t_2^8*t_3^8+2*t_0^3*t_2^7*t_3^9+t_0^3*t_2^6*t_3^10+t_0^3*t_2^5*t_3^11+t_0^2*t_1*t_2^13*t_3^3+2*t_0^2*t_1*t_2^12*t_3^4+4*t_0^2*t_1*t_2^11*t_3^5+6*t_0^2*t_1*t_2^10*t_3^6+8*t_0^2*t_1*t_2^9*t_3^7+8*t_0^2*t_1*t_2^8*t_3^8+8*t_0^2*t_1*t_2^7*t_3^9+6*t_0^2*t_1*t_2^6*t_3^10+4*t_0^2*t_1*t_2^5*t_3^11+2*t_0^2*t_1*t_2^4*t_3^12+t_0^2*t_1*t_2^3*t_3^13+t_0*t_1^2*t_2^13*t_3^3+2*t_0*t_1^2*t_2^12*t_3^4+4*t_0*t_1^2*t_2^11*t_3^5+6*t_0*t_1^2*t_2^10*t_3^6+8*t_0*t_1^2*t_2^9*t_3^7+8*t_0*t_1^2*t_2^8*t_3^8+8*t_0*t_1^2*t_2^7*t_3^9+6*t_0*t_1^2*t_2^6*t_3^10+4*t_0*t_1^2*t_2^5*t_3^11+2*t_0*t_1^2*t_2^4*t_3^12+t_0*t_1^2*t_2^3*t_3^13+t_1^3*t_2^11*t_3^5+t_1^3*t_2^10*t_3^6+2*t_1^3*t_2^9*t_3^7+2*t_1^3*t_2^8*t_3^8+2*t_1^3*t_2^7*t_3^9+t_1^3*t_2^6*t_3^10+t_1^3*t_2^5*t_3^11, (3,0) => 0, (4,0) => 0, (3,1) => t_0^4*t_2^13*t_3^8+2*t_0^4*t_2^12*t_3^9+2*t_0^4*t_2^11*t_3^10+2*t_0^4*t_2^10*t_3^11+2*t_0^4*t_2^9*t_3^12+t_0^4*t_2^8*t_3^13+t_0^3*t_1*t_2^16*t_3^5+3*t_0^3*t_1*t_2^15*t_3^6+6*t_0^3*t_1*t_2^14*t_3^7+10*t_0^3*t_1*t_2^13*t_3^8+14*t_0^3*t_1*t_2^12*t_3^9+16*t_0^3*t_1*t_2^11*t_3^10+16*t_0^3*t_1*t_2^10*t_3^11+14*t_0^3*t_1*t_2^9*t_3^12+10*t_0^3*t_1*t_2^8*t_3^13+6*t_0^3*t_1*t_2^7*t_3^14+3*t_0^3*t_1*t_2^6*t_3^15+t_0^3*t_1*t_2^5*t_3^16+t_0^2*t_1^2*t_2^17*t_3^4+3*t_0^2*t_1^2*t_2^16*t_3^5+7*t_0^2*t_1^2*t_2^15*t_3^6+13*t_0^2*t_1^2*t_2^14*t_3^7+20*t_0^2*t_1^2*t_2^13*t_3^8+26*t_0^2*t_1^2*t_2^12*t_3^9+30*t_0^2*t_1^2*t_2^11*t_3^10+30*t_0^2*t_1^2*t_2^10*t_3^11+26*t_0^2*t_1^2*t_2^9*t_3^12+20*t_0^2*t_1^2*t_2^8*t_3^13+13*t_0^2*t_1^2*t_2^7*t_3^14+7*t_0^2*t_1^2*t_2^6*t_3^15+3*t_0^2*t_1^2*t_2^5*t_3^16+t_0^2*t_1^2*t_2^4*t_3^17+t_0*t_1^3*t_2^16*t_3^5+3*t_0*t_1^3*t_2^15*t_3^6+6*t_0*t_1^3*t_2^14*t_3^7+10*t_0*t_1^3*t_2^13*t_3^8+14*t_0*t_1^3*t_2^12*t_3^9+16*t_0*t_1^3*t_2^11*t_3^10+16*t_0*t_1^3*t_2^10*t_3^11+14*t_0*t_1^3*t_2^9*t_3^12+10*t_0*t_1^3*t_2^8*t_3^13+6*t_0*t_1^3*t_2^7*t_3^14+3*t_0*t_1^3*t_2^6*t_3^15+t_0*t_1^3*t_2^5*t_3^16+t_1^4*t_2^13*t_3^8+2*t_1^4*t_2^12*t_3^9+2*t_1^4*t_2^11*t_3^10+2*t_1^4*t_2^10*t_3^11+2*t_1^4*t_2^9*t_3^12+t_1^4*t_2^8*t_3^13, (2,2) => 0, (5,0) => 0, (4,1) => t_0^5*t_2^14*t_3^12+t_0^5*t_2^13*t_3^13+t_0^5*t_2^12*t_3^14+t_0^4*t_1*t_2^18*t_3^8+3*t_0^4*t_1*t_2^17*t_3^9+6*t_0^4*t_1*t_2^16*t_3^10+9*t_0^4*t_1*t_2^15*t_3^11+12*t_0^4*t_1*t_2^14*t_3^12+13*t_0^4*t_1*t_2^13*t_3^13+12*t_0^4*t_1*t_2^12*t_3^14+9*t_0^4*t_1*t_2^11*t_3^15+6*t_0^4*t_1*t_2^10*t_3^16+3*t_0^4*t_1*t_2^9*t_3^17+t_0^4*t_1*t_2^8*t_3^18+t_0^3*t_1^2*t_2^20*t_3^6+3*t_0^3*t_1^2*t_2^19*t_3^7+8*t_0^3*t_1^2*t_2^18*t_3^8+15*t_0^3*t_1^2*t_2^17*t_3^9+25*t_0^3*t_1^2*t_2^16*t_3^10+34*t_0^3*t_1^2*t_2^15*t_3^11+42*t_0^3*t_1^2*t_2^14*t_3^12+44*t_0^3*t_1^2*t_2^13*t_3^13+42*t_0^3*t_1^2*t_2^12*t_3^14+34*t_0^3*t_1^2*t_2^11*t_3^15+25*t_0^3*t_1^2*t_2^10*t_3^16+15*t_0^3*t_1^2*t_2^9*t_3^17+8*t_0^3*t_1^2*t_2^8*t_3^18+3*t_0^3*t_1^2*t_2^7*t_3^19+t_0^3*t_1^2*t_2^6*t_3^20+t_0^2*t_1^3*t_2^20*t_3^6+3*t_0^2*t_1^3*t_2^19*t_3^7+8*t_0^2*t_1^3*t_2^18*t_3^8+15*t_0^2*t_1^3*t_2^17*t_3^9+25*t_0^2*t_1^3*t_2^16*t_3^10+34*t_0^2*t_1^3*t_2^15*t_3^11+42*t_0^2*t_1^3*t_2^14*t_3^12+44*t_0^2*t_1^3*t_2^13*t_3^13+42*t_0^2*t_1^3*t_2^12*t_3^14+34*t_0^2*t_1^3*t_2^11*t_3^15+25*t_0^2*t_1^3*t_2^10*t_3^16+15*t_0^2*t_1^3*t_2^9*t_3^17+8*t_0^2*t_1^3*t_2^8*t_3^18+3*t_0^2*t_1^3*t_2^7*t_3^19+t_0^2*t_1^3*t_2^6*t_3^20+t_0*t_1^4*t_2^18*t_3^8+3*t_0*t_1^4*t_2^17*t_3^9+6*t_0*t_1^4*t_2^16*t_3^10+9*t_0*t_1^4*t_2^15*t_3^11+12*t_0*t_1^4*t_2^14*t_3^12+13*t_0*t_1^4*t_2^13*t_3^13+12*t_0*t_1^4*t_2^12*t_3^14+9*t_0*t_1^4*t_2^11*t_3^15+6*t_0*t_1^4*t_2^10*t_3^16+3*t_0*t_1^4*t_2^9*t_3^17+t_0*t_1^4*t_2^8*t_3^18+t_1^5*t_2^14*t_3^12+t_1^5*t_2^13*t_3^13+t_1^5*t_2^12*t_3^14, (3,2) => 0, (6,0) => 0, (5,1) => t_0^5*t_1*t_2^19*t_3^12+2*t_0^5*t_1*t_2^18*t_3^13+3*t_0^5*t_1*t_2^17*t_3^14+4*t_0^5*t_1*t_2^16*t_3^15+4*t_0^5*t_1*t_2^15*t_3^16+3*t_0^5*t_1*t_2^14*t_3^17+2*t_0^5*t_1*t_2^13*t_3^18+t_0^5*t_1*t_2^12*t_3^19+t_0^4*t_1^2*t_2^22*t_3^9+3*t_0^4*t_1^2*t_2^21*t_3^10+7*t_0^4*t_1^2*t_2^20*t_3^11+13*t_0^4*t_1^2*t_2^19*t_3^12+20*t_0^4*t_1^2*t_2^18*t_3^13+26*t_0^4*t_1^2*t_2^17*t_3^14+30*t_0^4*t_1^2*t_2^16*t_3^15+30*t_0^4*t_1^2*t_2^15*t_3^16+26*t_0^4*t_1^2*t_2^14*t_3^17+20*t_0^4*t_1^2*t_2^13*t_3^18+13*t_0^4*t_1^2*t_2^12*t_3^19+7*t_0^4*t_1^2*t_2^11*t_3^20+3*t_0^4*t_1^2*t_2^10*t_3^21+t_0^4*t_1^2*t_2^9*t_3^22+t_0^3*t_1^3*t_2^23*t_3^8+3*t_0^3*t_1^3*t_2^22*t_3^9+8*t_0^3*t_1^3*t_2^21*t_3^10+16*t_0^3*t_1^3*t_2^20*t_3^11+27*t_0^3*t_1^3*t_2^19*t_3^12+39*t_0^3*t_1^3*t_2^18*t_3^13+50*t_0^3*t_1^3*t_2^17*t_3^14+56*t_0^3*t_1^3*t_2^16*t_3^15+56*t_0^3*t_1^3*t_2^15*t_3^16+50*t_0^3*t_1^3*t_2^14*t_3^17+39*t_0^3*t_1^3*t_2^13*t_3^18+27*t_0^3*t_1^3*t_2^12*t_3^19+16*t_0^3*t_1^3*t_2^11*t_3^20+8*t_0^3*t_1^3*t_2^10*t_3^21+3*t_0^3*t_1^3*t_2^9*t_3^22+t_0^3*t_1^3*t_2^8*t_3^23+t_0^2*t_1^4*t_2^22*t_3^9+3*t_0^2*t_1^4*t_2^21*t_3^10+7*t_0^2*t_1^4*t_2^20*t_3^11+13*t_0^2*t_1^4*t_2^19*t_3^12+20*t_0^2*t_1^4*t_2^18*t_3^13+26*t_0^2*t_1^4*t_2^17*t_3^14+30*t_0^2*t_1^4*t_2^16*t_3^15+30*t_0^2*t_1^4*t_2^15*t_3^16+26*t_0^2*t_1^4*t_2^14*t_3^17+20*t_0^2*t_1^4*t_2^13*t_3^18+13*t_0^2*t_1^4*t_2^12*t_3^19+7*t_0^2*t_1^4*t_2^11*t_3^20+3*t_0^2*t_1^4*t_2^10*t_3^21+t_0^2*t_1^4*t_2^9*t_3^22+t_0*t_1^5*t_2^19*t_3^12+2*t_0*t_1^5*t_2^18*t_3^13+3*t_0*t_1^5*t_2^17*t_3^14+4*t_0*t_1^5*t_2^16*t_3^15+4*t_0*t_1^5*t_2^15*t_3^16+3*t_0*t_1^5*t_2^14*t_3^17+2*t_0*t_1^5*t_2^13*t_3^18+t_0*t_1^5*t_2^12*t_3^19, (4,2) => 0};
end;