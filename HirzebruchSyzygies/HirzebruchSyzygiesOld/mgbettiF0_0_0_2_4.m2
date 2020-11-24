A := QQ[t_0,t_1,t_2,t_3];
--mb stands for Multigraded Betti numbers
mb0024 = new HashTable from {(5,2) => 0, (7,0) => 0, (6,1) => t_0^11*t_1^3*t_2^18*t_3^10+3*t_0^11*t_1^3*t_2^17*t_3^11+5*t_0^11*t_1^3*t_2^16*t_3^12+7*t_0^11*t_1^3*t_2^15*t_3^13+8*t_0^11*t_1^3*t_2^14*t_3^14+7*t_0^11*t_1^3*t_2^13*t_3^15+5*t_0^11*t_1^3*t_2^12*t_3^16+3*t_0^11*t_1^3*t_2^11*t_3^17+t_0^11*t_1^3*t_2^10*t_3^18+t_0^10*t_1^4*t_2^20*t_3^8+4*t_0^10*t_1^4*t_2^19*t_3^9+11*t_0^10*t_1^4*t_2^18*t_3^10+22*t_0^10*t_1^4*t_2^17*t_3^11+34*t_0^10*t_1^4*t_2^16*t_3^12+44*t_0^10*t_1^4*t_2^15*t_3^13+48*t_0^10*t_1^4*t_2^14*t_3^14+44*t_0^10*t_1^4*t_2^13*t_3^15+34*t_0^10*t_1^4*t_2^12*t_3^16+22*t_0^10*t_1^4*t_2^11*t_3^17+11*t_0^10*t_1^4*t_2^10*t_3^18+4*t_0^10*t_1^4*t_2^9*t_3^19+t_0^10*t_1^4*t_2^8*t_3^20+t_0^9*t_1^5*t_2^21*t_3^7+6*t_0^9*t_1^5*t_2^20*t_3^8+18*t_0^9*t_1^5*t_2^19*t_3^9+40*t_0^9*t_1^5*t_2^18*t_3^10+71*t_0^9*t_1^5*t_2^17*t_3^11+104*t_0^9*t_1^5*t_2^16*t_3^12+130*t_0^9*t_1^5*t_2^15*t_3^13+140*t_0^9*t_1^5*t_2^14*t_3^14+130*t_0^9*t_1^5*t_2^13*t_3^15+104*t_0^9*t_1^5*t_2^12*t_3^16+71*t_0^9*t_1^5*t_2^11*t_3^17+40*t_0^9*t_1^5*t_2^10*t_3^18+18*t_0^9*t_1^5*t_2^9*t_3^19+6*t_0^9*t_1^5*t_2^8*t_3^20+t_0^9*t_1^5*t_2^7*t_3^21+t_0^8*t_1^6*t_2^22*t_3^6+5*t_0^8*t_1^6*t_2^21*t_3^7+17*t_0^8*t_1^6*t_2^20*t_3^8+42*t_0^8*t_1^6*t_2^19*t_3^9+84*t_0^8*t_1^6*t_2^18*t_3^10+140*t_0^8*t_1^6*t_2^17*t_3^11+198*t_0^8*t_1^6*t_2^16*t_3^12+243*t_0^8*t_1^6*t_2^15*t_3^13+260*t_0^8*t_1^6*t_2^14*t_3^14+243*t_0^8*t_1^6*t_2^13*t_3^15+198*t_0^8*t_1^6*t_2^12*t_3^16+140*t_0^8*t_1^6*t_2^11*t_3^17+84*t_0^8*t_1^6*t_2^10*t_3^18+42*t_0^8*t_1^6*t_2^9*t_3^19+17*t_0^8*t_1^6*t_2^8*t_3^20+5*t_0^8*t_1^6*t_2^7*t_3^21+t_0^8*t_1^6*t_2^6*t_3^22+t_0^7*t_1^7*t_2^22*t_3^6+6*t_0^7*t_1^7*t_2^21*t_3^7+21*t_0^7*t_1^7*t_2^20*t_3^8+52*t_0^7*t_1^7*t_2^19*t_3^9+103*t_0^7*t_1^7*t_2^18*t_3^10+170*t_0^7*t_1^7*t_2^17*t_3^11+239*t_0^7*t_1^7*t_2^16*t_3^12+292*t_0^7*t_1^7*t_2^15*t_3^13+312*t_0^7*t_1^7*t_2^14*t_3^14+292*t_0^7*t_1^7*t_2^13*t_3^15+239*t_0^7*t_1^7*t_2^12*t_3^16+170*t_0^7*t_1^7*t_2^11*t_3^17+103*t_0^7*t_1^7*t_2^10*t_3^18+52*t_0^7*t_1^7*t_2^9*t_3^19+21*t_0^7*t_1^7*t_2^8*t_3^20+6*t_0^7*t_1^7*t_2^7*t_3^21+t_0^7*t_1^7*t_2^6*t_3^22+t_0^6*t_1^8*t_2^22*t_3^6+5*t_0^6*t_1^8*t_2^21*t_3^7+17*t_0^6*t_1^8*t_2^20*t_3^8+42*t_0^6*t_1^8*t_2^19*t_3^9+84*t_0^6*t_1^8*t_2^18*t_3^10+140*t_0^6*t_1^8*t_2^17*t_3^11+198*t_0^6*t_1^8*t_2^16*t_3^12+243*t_0^6*t_1^8*t_2^15*t_3^13+260*t_0^6*t_1^8*t_2^14*t_3^14+243*t_0^6*t_1^8*t_2^13*t_3^15+198*t_0^6*t_1^8*t_2^12*t_3^16+140*t_0^6*t_1^8*t_2^11*t_3^17+84*t_0^6*t_1^8*t_2^10*t_3^18+42*t_0^6*t_1^8*t_2^9*t_3^19+17*t_0^6*t_1^8*t_2^8*t_3^20+5*t_0^6*t_1^8*t_2^7*t_3^21+t_0^6*t_1^8*t_2^6*t_3^22+t_0^5*t_1^9*t_2^21*t_3^7+6*t_0^5*t_1^9*t_2^20*t_3^8+18*t_0^5*t_1^9*t_2^19*t_3^9+40*t_0^5*t_1^9*t_2^18*t_3^10+71*t_0^5*t_1^9*t_2^17*t_3^11+104*t_0^5*t_1^9*t_2^16*t_3^12+130*t_0^5*t_1^9*t_2^15*t_3^13+140*t_0^5*t_1^9*t_2^14*t_3^14+130*t_0^5*t_1^9*t_2^13*t_3^15+104*t_0^5*t_1^9*t_2^12*t_3^16+71*t_0^5*t_1^9*t_2^11*t_3^17+40*t_0^5*t_1^9*t_2^10*t_3^18+18*t_0^5*t_1^9*t_2^9*t_3^19+6*t_0^5*t_1^9*t_2^8*t_3^20+t_0^5*t_1^9*t_2^7*t_3^21+t_0^4*t_1^10*t_2^20*t_3^8+4*t_0^4*t_1^10*t_2^19*t_3^9+11*t_0^4*t_1^10*t_2^18*t_3^10+22*t_0^4*t_1^10*t_2^17*t_3^11+34*t_0^4*t_1^10*t_2^16*t_3^12+44*t_0^4*t_1^10*t_2^15*t_3^13+48*t_0^4*t_1^10*t_2^14*t_3^14+44*t_0^4*t_1^10*t_2^13*t_3^15+34*t_0^4*t_1^10*t_2^12*t_3^16+22*t_0^4*t_1^10*t_2^11*t_3^17+11*t_0^4*t_1^10*t_2^10*t_3^18+4*t_0^4*t_1^10*t_2^9*t_3^19+t_0^4*t_1^10*t_2^8*t_3^20+t_0^3*t_1^11*t_2^18*t_3^10+3*t_0^3*t_1^11*t_2^17*t_3^11+5*t_0^3*t_1^11*t_2^16*t_3^12+7*t_0^3*t_1^11*t_2^15*t_3^13+8*t_0^3*t_1^11*t_2^14*t_3^14+7*t_0^3*t_1^11*t_2^13*t_3^15+5*t_0^3*t_1^11*t_2^12*t_3^16+3*t_0^3*t_1^11*t_2^11*t_3^17+t_0^3*t_1^11*t_2^10*t_3^18, (6,2) => 0, (8,0) => 0, (7,1) => 2*t_0^12*t_1^4*t_2^19*t_3^13+3*t_0^12*t_1^4*t_2^18*t_3^14+5*t_0^12*t_1^4*t_2^17*t_3^15+5*t_0^12*t_1^4*t_2^16*t_3^16+5*t_0^12*t_1^4*t_2^15*t_3^17+3*t_0^12*t_1^4*t_2^14*t_3^18+2*t_0^12*t_1^4*t_2^13*t_3^19+t_0^11*t_1^5*t_2^22*t_3^10+3*t_0^11*t_1^5*t_2^21*t_3^11+8*t_0^11*t_1^5*t_2^20*t_3^12+16*t_0^11*t_1^5*t_2^19*t_3^13+25*t_0^11*t_1^5*t_2^18*t_3^14+33*t_0^11*t_1^5*t_2^17*t_3^15+36*t_0^11*t_1^5*t_2^16*t_3^16+33*t_0^11*t_1^5*t_2^15*t_3^17+25*t_0^11*t_1^5*t_2^14*t_3^18+16*t_0^11*t_1^5*t_2^13*t_3^19+8*t_0^11*t_1^5*t_2^12*t_3^20+3*t_0^11*t_1^5*t_2^11*t_3^21+t_0^11*t_1^5*t_2^10*t_3^22+t_0^10*t_1^6*t_2^23*t_3^9+5*t_0^10*t_1^6*t_2^22*t_3^10+15*t_0^10*t_1^6*t_2^21*t_3^11+31*t_0^10*t_1^6*t_2^20*t_3^12+56*t_0^10*t_1^6*t_2^19*t_3^13+80*t_0^10*t_1^6*t_2^18*t_3^14+102*t_0^10*t_1^6*t_2^17*t_3^15+108*t_0^10*t_1^6*t_2^16*t_3^16+102*t_0^10*t_1^6*t_2^15*t_3^17+80*t_0^10*t_1^6*t_2^14*t_3^18+56*t_0^10*t_1^6*t_2^13*t_3^19+31*t_0^10*t_1^6*t_2^12*t_3^20+15*t_0^10*t_1^6*t_2^11*t_3^21+5*t_0^10*t_1^6*t_2^10*t_3^22+t_0^10*t_1^6*t_2^9*t_3^23+t_0^9*t_1^7*t_2^24*t_3^8+4*t_0^9*t_1^7*t_2^23*t_3^9+14*t_0^9*t_1^7*t_2^22*t_3^10+34*t_0^9*t_1^7*t_2^21*t_3^11+67*t_0^9*t_1^7*t_2^20*t_3^12+111*t_0^9*t_1^7*t_2^19*t_3^13+156*t_0^9*t_1^7*t_2^18*t_3^14+191*t_0^9*t_1^7*t_2^17*t_3^15+204*t_0^9*t_1^7*t_2^16*t_3^16+191*t_0^9*t_1^7*t_2^15*t_3^17+156*t_0^9*t_1^7*t_2^14*t_3^18+111*t_0^9*t_1^7*t_2^13*t_3^19+67*t_0^9*t_1^7*t_2^12*t_3^20+34*t_0^9*t_1^7*t_2^11*t_3^21+14*t_0^9*t_1^7*t_2^10*t_3^22+4*t_0^9*t_1^7*t_2^9*t_3^23+t_0^9*t_1^7*t_2^8*t_3^24+t_0^8*t_1^8*t_2^24*t_3^8+6*t_0^8*t_1^8*t_2^23*t_3^9+18*t_0^8*t_1^8*t_2^22*t_3^10+44*t_0^8*t_1^8*t_2^21*t_3^11+83*t_0^8*t_1^8*t_2^20*t_3^12+138*t_0^8*t_1^8*t_2^19*t_3^13+190*t_0^8*t_1^8*t_2^18*t_3^14+234*t_0^8*t_1^8*t_2^17*t_3^15+247*t_0^8*t_1^8*t_2^16*t_3^16+234*t_0^8*t_1^8*t_2^15*t_3^17+190*t_0^8*t_1^8*t_2^14*t_3^18+138*t_0^8*t_1^8*t_2^13*t_3^19+83*t_0^8*t_1^8*t_2^12*t_3^20+44*t_0^8*t_1^8*t_2^11*t_3^21+18*t_0^8*t_1^8*t_2^10*t_3^22+6*t_0^8*t_1^8*t_2^9*t_3^23+t_0^8*t_1^8*t_2^8*t_3^24+t_0^7*t_1^9*t_2^24*t_3^8+4*t_0^7*t_1^9*t_2^23*t_3^9+14*t_0^7*t_1^9*t_2^22*t_3^10+34*t_0^7*t_1^9*t_2^21*t_3^11+67*t_0^7*t_1^9*t_2^20*t_3^12+111*t_0^7*t_1^9*t_2^19*t_3^13+156*t_0^7*t_1^9*t_2^18*t_3^14+191*t_0^7*t_1^9*t_2^17*t_3^15+204*t_0^7*t_1^9*t_2^16*t_3^16+191*t_0^7*t_1^9*t_2^15*t_3^17+156*t_0^7*t_1^9*t_2^14*t_3^18+111*t_0^7*t_1^9*t_2^13*t_3^19+67*t_0^7*t_1^9*t_2^12*t_3^20+34*t_0^7*t_1^9*t_2^11*t_3^21+14*t_0^7*t_1^9*t_2^10*t_3^22+4*t_0^7*t_1^9*t_2^9*t_3^23+t_0^7*t_1^9*t_2^8*t_3^24+t_0^6*t_1^10*t_2^23*t_3^9+5*t_0^6*t_1^10*t_2^22*t_3^10+15*t_0^6*t_1^10*t_2^21*t_3^11+31*t_0^6*t_1^10*t_2^20*t_3^12+56*t_0^6*t_1^10*t_2^19*t_3^13+80*t_0^6*t_1^10*t_2^18*t_3^14+102*t_0^6*t_1^10*t_2^17*t_3^15+108*t_0^6*t_1^10*t_2^16*t_3^16+102*t_0^6*t_1^10*t_2^15*t_3^17+80*t_0^6*t_1^10*t_2^14*t_3^18+56*t_0^6*t_1^10*t_2^13*t_3^19+31*t_0^6*t_1^10*t_2^12*t_3^20+15*t_0^6*t_1^10*t_2^11*t_3^21+5*t_0^6*t_1^10*t_2^10*t_3^22+t_0^6*t_1^10*t_2^9*t_3^23+t_0^5*t_1^11*t_2^22*t_3^10+3*t_0^5*t_1^11*t_2^21*t_3^11+8*t_0^5*t_1^11*t_2^20*t_3^12+16*t_0^5*t_1^11*t_2^19*t_3^13+25*t_0^5*t_1^11*t_2^18*t_3^14+33*t_0^5*t_1^11*t_2^17*t_3^15+36*t_0^5*t_1^11*t_2^16*t_3^16+33*t_0^5*t_1^11*t_2^15*t_3^17+25*t_0^5*t_1^11*t_2^14*t_3^18+16*t_0^5*t_1^11*t_2^13*t_3^19+8*t_0^5*t_1^11*t_2^12*t_3^20+3*t_0^5*t_1^11*t_2^11*t_3^21+t_0^5*t_1^11*t_2^10*t_3^22+2*t_0^4*t_1^12*t_2^19*t_3^13+3*t_0^4*t_1^12*t_2^18*t_3^14+5*t_0^4*t_1^12*t_2^17*t_3^15+5*t_0^4*t_1^12*t_2^16*t_3^16+5*t_0^4*t_1^12*t_2^15*t_3^17+3*t_0^4*t_1^12*t_2^14*t_3^18+2*t_0^4*t_1^12*t_2^13*t_3^19, (7,2) => 0, (9,0) => 0, (8,1) => t_0^13*t_1^5*t_2^20*t_3^16+2*t_0^13*t_1^5*t_2^19*t_3^17+2*t_0^13*t_1^5*t_2^18*t_3^18+2*t_0^13*t_1^5*t_2^17*t_3^19+t_0^13*t_1^5*t_2^16*t_3^20+t_0^12*t_1^6*t_2^23*t_3^13+3*t_0^12*t_1^6*t_2^22*t_3^14+6*t_0^12*t_1^6*t_2^21*t_3^15+11*t_0^12*t_1^6*t_2^20*t_3^16+15*t_0^12*t_1^6*t_2^19*t_3^17+16*t_0^12*t_1^6*t_2^18*t_3^18+15*t_0^12*t_1^6*t_2^17*t_3^19+11*t_0^12*t_1^6*t_2^16*t_3^20+6*t_0^12*t_1^6*t_2^15*t_3^21+3*t_0^12*t_1^6*t_2^14*t_3^22+t_0^12*t_1^6*t_2^13*t_3^23+t_0^11*t_1^7*t_2^25*t_3^11+3*t_0^11*t_1^7*t_2^24*t_3^12+8*t_0^11*t_1^7*t_2^23*t_3^13+17*t_0^11*t_1^7*t_2^22*t_3^14+29*t_0^11*t_1^7*t_2^21*t_3^15+43*t_0^11*t_1^7*t_2^20*t_3^16+54*t_0^11*t_1^7*t_2^19*t_3^17+58*t_0^11*t_1^7*t_2^18*t_3^18+54*t_0^11*t_1^7*t_2^17*t_3^19+43*t_0^11*t_1^7*t_2^16*t_3^20+29*t_0^11*t_1^7*t_2^15*t_3^21+17*t_0^11*t_1^7*t_2^14*t_3^22+8*t_0^11*t_1^7*t_2^13*t_3^23+3*t_0^11*t_1^7*t_2^12*t_3^24+t_0^11*t_1^7*t_2^11*t_3^25+2*t_0^10*t_1^8*t_2^25*t_3^11+7*t_0^10*t_1^8*t_2^24*t_3^12+18*t_0^10*t_1^8*t_2^23*t_3^13+36*t_0^10*t_1^8*t_2^22*t_3^14+59*t_0^10*t_1^8*t_2^21*t_3^15+84*t_0^10*t_1^8*t_2^20*t_3^16+103*t_0^10*t_1^8*t_2^19*t_3^17+110*t_0^10*t_1^8*t_2^18*t_3^18+103*t_0^10*t_1^8*t_2^17*t_3^19+84*t_0^10*t_1^8*t_2^16*t_3^20+59*t_0^10*t_1^8*t_2^15*t_3^21+36*t_0^10*t_1^8*t_2^14*t_3^22+18*t_0^10*t_1^8*t_2^13*t_3^23+7*t_0^10*t_1^8*t_2^12*t_3^24+2*t_0^10*t_1^8*t_2^11*t_3^25+t_0^9*t_1^9*t_2^26*t_3^10+4*t_0^9*t_1^9*t_2^25*t_3^11+11*t_0^9*t_1^9*t_2^24*t_3^12+26*t_0^9*t_1^9*t_2^23*t_3^13+49*t_0^9*t_1^9*t_2^22*t_3^14+78*t_0^9*t_1^9*t_2^21*t_3^15+109*t_0^9*t_1^9*t_2^20*t_3^16+132*t_0^9*t_1^9*t_2^19*t_3^17+140*t_0^9*t_1^9*t_2^18*t_3^18+132*t_0^9*t_1^9*t_2^17*t_3^19+109*t_0^9*t_1^9*t_2^16*t_3^20+78*t_0^9*t_1^9*t_2^15*t_3^21+49*t_0^9*t_1^9*t_2^14*t_3^22+26*t_0^9*t_1^9*t_2^13*t_3^23+11*t_0^9*t_1^9*t_2^12*t_3^24+4*t_0^9*t_1^9*t_2^11*t_3^25+t_0^9*t_1^9*t_2^10*t_3^26+2*t_0^8*t_1^10*t_2^25*t_3^11+7*t_0^8*t_1^10*t_2^24*t_3^12+18*t_0^8*t_1^10*t_2^23*t_3^13+36*t_0^8*t_1^10*t_2^22*t_3^14+59*t_0^8*t_1^10*t_2^21*t_3^15+84*t_0^8*t_1^10*t_2^20*t_3^16+103*t_0^8*t_1^10*t_2^19*t_3^17+110*t_0^8*t_1^10*t_2^18*t_3^18+103*t_0^8*t_1^10*t_2^17*t_3^19+84*t_0^8*t_1^10*t_2^16*t_3^20+59*t_0^8*t_1^10*t_2^15*t_3^21+36*t_0^8*t_1^10*t_2^14*t_3^22+18*t_0^8*t_1^10*t_2^13*t_3^23+7*t_0^8*t_1^10*t_2^12*t_3^24+2*t_0^8*t_1^10*t_2^11*t_3^25+t_0^7*t_1^11*t_2^25*t_3^11+3*t_0^7*t_1^11*t_2^24*t_3^12+8*t_0^7*t_1^11*t_2^23*t_3^13+17*t_0^7*t_1^11*t_2^22*t_3^14+29*t_0^7*t_1^11*t_2^21*t_3^15+43*t_0^7*t_1^11*t_2^20*t_3^16+54*t_0^7*t_1^11*t_2^19*t_3^17+58*t_0^7*t_1^11*t_2^18*t_3^18+54*t_0^7*t_1^11*t_2^17*t_3^19+43*t_0^7*t_1^11*t_2^16*t_3^20+29*t_0^7*t_1^11*t_2^15*t_3^21+17*t_0^7*t_1^11*t_2^14*t_3^22+8*t_0^7*t_1^11*t_2^13*t_3^23+3*t_0^7*t_1^11*t_2^12*t_3^24+t_0^7*t_1^11*t_2^11*t_3^25+t_0^6*t_1^12*t_2^23*t_3^13+3*t_0^6*t_1^12*t_2^22*t_3^14+6*t_0^6*t_1^12*t_2^21*t_3^15+11*t_0^6*t_1^12*t_2^20*t_3^16+15*t_0^6*t_1^12*t_2^19*t_3^17+16*t_0^6*t_1^12*t_2^18*t_3^18+15*t_0^6*t_1^12*t_2^17*t_3^19+11*t_0^6*t_1^12*t_2^16*t_3^20+6*t_0^6*t_1^12*t_2^15*t_3^21+3*t_0^6*t_1^12*t_2^14*t_3^22+t_0^6*t_1^12*t_2^13*t_3^23+t_0^5*t_1^13*t_2^20*t_3^16+2*t_0^5*t_1^13*t_2^19*t_3^17+2*t_0^5*t_1^13*t_2^18*t_3^18+2*t_0^5*t_1^13*t_2^17*t_3^19+t_0^5*t_1^13*t_2^16*t_3^20, (8,2) => 0, (10,0) => 0, (9,1) => t_0^14*t_1^6*t_2^20*t_3^20+t_0^13*t_1^7*t_2^23*t_3^17+2*t_0^13*t_1^7*t_2^22*t_3^18+3*t_0^13*t_1^7*t_2^21*t_3^19+4*t_0^13*t_1^7*t_2^20*t_3^20+3*t_0^13*t_1^7*t_2^19*t_3^21+2*t_0^13*t_1^7*t_2^18*t_3^22+t_0^13*t_1^7*t_2^17*t_3^23+t_0^12*t_1^8*t_2^26*t_3^14+2*t_0^12*t_1^8*t_2^25*t_3^15+5*t_0^12*t_1^8*t_2^24*t_3^16+8*t_0^12*t_1^8*t_2^23*t_3^17+13*t_0^12*t_1^8*t_2^22*t_3^18+15*t_0^12*t_1^8*t_2^21*t_3^19+18*t_0^12*t_1^8*t_2^20*t_3^20+15*t_0^12*t_1^8*t_2^19*t_3^21+13*t_0^12*t_1^8*t_2^18*t_3^22+8*t_0^12*t_1^8*t_2^17*t_3^23+5*t_0^12*t_1^8*t_2^16*t_3^24+2*t_0^12*t_1^8*t_2^15*t_3^25+t_0^12*t_1^8*t_2^14*t_3^26+t_0^11*t_1^9*t_2^27*t_3^13+3*t_0^11*t_1^9*t_2^26*t_3^14+7*t_0^11*t_1^9*t_2^25*t_3^15+13*t_0^11*t_1^9*t_2^24*t_3^16+21*t_0^11*t_1^9*t_2^23*t_3^17+29*t_0^11*t_1^9*t_2^22*t_3^18+35*t_0^11*t_1^9*t_2^21*t_3^19+38*t_0^11*t_1^9*t_2^20*t_3^20+35*t_0^11*t_1^9*t_2^19*t_3^21+29*t_0^11*t_1^9*t_2^18*t_3^22+21*t_0^11*t_1^9*t_2^17*t_3^23+13*t_0^11*t_1^9*t_2^16*t_3^24+7*t_0^11*t_1^9*t_2^15*t_3^25+3*t_0^11*t_1^9*t_2^14*t_3^26+t_0^11*t_1^9*t_2^13*t_3^27+t_0^10*t_1^10*t_2^27*t_3^13+4*t_0^10*t_1^10*t_2^26*t_3^14+9*t_0^10*t_1^10*t_2^25*t_3^15+18*t_0^10*t_1^10*t_2^24*t_3^16+27*t_0^10*t_1^10*t_2^23*t_3^17+38*t_0^10*t_1^10*t_2^22*t_3^18+44*t_0^10*t_1^10*t_2^21*t_3^19+49*t_0^10*t_1^10*t_2^20*t_3^20+44*t_0^10*t_1^10*t_2^19*t_3^21+38*t_0^10*t_1^10*t_2^18*t_3^22+27*t_0^10*t_1^10*t_2^17*t_3^23+18*t_0^10*t_1^10*t_2^16*t_3^24+9*t_0^10*t_1^10*t_2^15*t_3^25+4*t_0^10*t_1^10*t_2^14*t_3^26+t_0^10*t_1^10*t_2^13*t_3^27+t_0^9*t_1^11*t_2^27*t_3^13+3*t_0^9*t_1^11*t_2^26*t_3^14+7*t_0^9*t_1^11*t_2^25*t_3^15+13*t_0^9*t_1^11*t_2^24*t_3^16+21*t_0^9*t_1^11*t_2^23*t_3^17+29*t_0^9*t_1^11*t_2^22*t_3^18+35*t_0^9*t_1^11*t_2^21*t_3^19+38*t_0^9*t_1^11*t_2^20*t_3^20+35*t_0^9*t_1^11*t_2^19*t_3^21+29*t_0^9*t_1^11*t_2^18*t_3^22+21*t_0^9*t_1^11*t_2^17*t_3^23+13*t_0^9*t_1^11*t_2^16*t_3^24+7*t_0^9*t_1^11*t_2^15*t_3^25+3*t_0^9*t_1^11*t_2^14*t_3^26+t_0^9*t_1^11*t_2^13*t_3^27+t_0^8*t_1^12*t_2^26*t_3^14+2*t_0^8*t_1^12*t_2^25*t_3^15+5*t_0^8*t_1^12*t_2^24*t_3^16+8*t_0^8*t_1^12*t_2^23*t_3^17+13*t_0^8*t_1^12*t_2^22*t_3^18+15*t_0^8*t_1^12*t_2^21*t_3^19+18*t_0^8*t_1^12*t_2^20*t_3^20+15*t_0^8*t_1^12*t_2^19*t_3^21+13*t_0^8*t_1^12*t_2^18*t_3^22+8*t_0^8*t_1^12*t_2^17*t_3^23+5*t_0^8*t_1^12*t_2^16*t_3^24+2*t_0^8*t_1^12*t_2^15*t_3^25+t_0^8*t_1^12*t_2^14*t_3^26+t_0^7*t_1^13*t_2^23*t_3^17+2*t_0^7*t_1^13*t_2^22*t_3^18+3*t_0^7*t_1^13*t_2^21*t_3^19+4*t_0^7*t_1^13*t_2^20*t_3^20+3*t_0^7*t_1^13*t_2^19*t_3^21+2*t_0^7*t_1^13*t_2^18*t_3^22+t_0^7*t_1^13*t_2^17*t_3^23+t_0^6*t_1^14*t_2^20*t_3^20, (9,2) => 0, (11,0) => 0, (10,1) => t_0^12*t_1^10*t_2^28*t_3^16+2*t_0^12*t_1^10*t_2^27*t_3^17+3*t_0^12*t_1^10*t_2^26*t_3^18+4*t_0^12*t_1^10*t_2^25*t_3^19+4*t_0^12*t_1^10*t_2^24*t_3^20+4*t_0^12*t_1^10*t_2^23*t_3^21+4*t_0^12*t_1^10*t_2^22*t_3^22+4*t_0^12*t_1^10*t_2^21*t_3^23+4*t_0^12*t_1^10*t_2^20*t_3^24+4*t_0^12*t_1^10*t_2^19*t_3^25+3*t_0^12*t_1^10*t_2^18*t_3^26+2*t_0^12*t_1^10*t_2^17*t_3^27+t_0^12*t_1^10*t_2^16*t_3^28+t_0^11*t_1^11*t_2^28*t_3^16+2*t_0^11*t_1^11*t_2^27*t_3^17+3*t_0^11*t_1^11*t_2^26*t_3^18+4*t_0^11*t_1^11*t_2^25*t_3^19+4*t_0^11*t_1^11*t_2^24*t_3^20+4*t_0^11*t_1^11*t_2^23*t_3^21+4*t_0^11*t_1^11*t_2^22*t_3^22+4*t_0^11*t_1^11*t_2^21*t_3^23+4*t_0^11*t_1^11*t_2^20*t_3^24+4*t_0^11*t_1^11*t_2^19*t_3^25+3*t_0^11*t_1^11*t_2^18*t_3^26+2*t_0^11*t_1^11*t_2^17*t_3^27+t_0^11*t_1^11*t_2^16*t_3^28+t_0^10*t_1^12*t_2^28*t_3^16+2*t_0^10*t_1^12*t_2^27*t_3^17+3*t_0^10*t_1^12*t_2^26*t_3^18+4*t_0^10*t_1^12*t_2^25*t_3^19+4*t_0^10*t_1^12*t_2^24*t_3^20+4*t_0^10*t_1^12*t_2^23*t_3^21+4*t_0^10*t_1^12*t_2^22*t_3^22+4*t_0^10*t_1^12*t_2^21*t_3^23+4*t_0^10*t_1^12*t_2^20*t_3^24+4*t_0^10*t_1^12*t_2^19*t_3^25+3*t_0^10*t_1^12*t_2^18*t_3^26+2*t_0^10*t_1^12*t_2^17*t_3^27+t_0^10*t_1^12*t_2^16*t_3^28, (10,2) => t_0^14*t_1^10*t_2^26*t_3^22+t_0^14*t_1^10*t_2^25*t_3^23+2*t_0^14*t_1^10*t_2^24*t_3^24+t_0^14*t_1^10*t_2^23*t_3^25+t_0^14*t_1^10*t_2^22*t_3^26+t_0^13*t_1^11*t_2^27*t_3^21+2*t_0^13*t_1^11*t_2^26*t_3^22+3*t_0^13*t_1^11*t_2^25*t_3^23+4*t_0^13*t_1^11*t_2^24*t_3^24+3*t_0^13*t_1^11*t_2^23*t_3^25+2*t_0^13*t_1^11*t_2^22*t_3^26+t_0^13*t_1^11*t_2^21*t_3^27+t_0^12*t_1^12*t_2^27*t_3^21+3*t_0^12*t_1^12*t_2^26*t_3^22+4*t_0^12*t_1^12*t_2^25*t_3^23+6*t_0^12*t_1^12*t_2^24*t_3^24+4*t_0^12*t_1^12*t_2^23*t_3^25+3*t_0^12*t_1^12*t_2^22*t_3^26+t_0^12*t_1^12*t_2^21*t_3^27+t_0^11*t_1^13*t_2^27*t_3^21+2*t_0^11*t_1^13*t_2^26*t_3^22+3*t_0^11*t_1^13*t_2^25*t_3^23+4*t_0^11*t_1^13*t_2^24*t_3^24+3*t_0^11*t_1^13*t_2^23*t_3^25+2*t_0^11*t_1^13*t_2^22*t_3^26+t_0^11*t_1^13*t_2^21*t_3^27+t_0^10*t_1^14*t_2^26*t_3^22+t_0^10*t_1^14*t_2^25*t_3^23+2*t_0^10*t_1^14*t_2^24*t_3^24+t_0^10*t_1^14*t_2^23*t_3^25+t_0^10*t_1^14*t_2^22*t_3^26, (12,0) => 0, (11,1) => t_0^12*t_1^12*t_2^29*t_3^19+t_0^12*t_1^12*t_2^28*t_3^20+t_0^12*t_1^12*t_2^27*t_3^21+t_0^12*t_1^12*t_2^26*t_3^22+t_0^12*t_1^12*t_2^25*t_3^23+t_0^12*t_1^12*t_2^24*t_3^24+t_0^12*t_1^12*t_2^23*t_3^25+t_0^12*t_1^12*t_2^22*t_3^26+t_0^12*t_1^12*t_2^21*t_3^27+t_0^12*t_1^12*t_2^20*t_3^28+t_0^12*t_1^12*t_2^19*t_3^29, (11,2) => t_0^14*t_1^12*t_2^28*t_3^24+2*t_0^14*t_1^12*t_2^27*t_3^25+2*t_0^14*t_1^12*t_2^26*t_3^26+2*t_0^14*t_1^12*t_2^25*t_3^27+t_0^14*t_1^12*t_2^24*t_3^28+t_0^13*t_1^13*t_2^28*t_3^24+2*t_0^13*t_1^13*t_2^27*t_3^25+2*t_0^13*t_1^13*t_2^26*t_3^26+2*t_0^13*t_1^13*t_2^25*t_3^27+t_0^13*t_1^13*t_2^24*t_3^28+t_0^12*t_1^14*t_2^28*t_3^24+2*t_0^12*t_1^14*t_2^27*t_3^25+2*t_0^12*t_1^14*t_2^26*t_3^26+2*t_0^12*t_1^14*t_2^25*t_3^27+t_0^12*t_1^14*t_2^24*t_3^28, (12,1) => 0, (12,2) => t_0^14*t_1^14*t_2^29*t_3^27+t_0^14*t_1^14*t_2^28*t_3^28+t_0^14*t_1^14*t_2^27*t_3^29, (13,2) => 0, (0,0) => 1, (1,0) => 0, (2,0) => 0, (1,1) => t_0^4*t_2^6*t_3^2+t_0^4*t_2^5*t_3^3+2*t_0^4*t_2^4*t_3^4+t_0^4*t_2^3*t_3^5+t_0^4*t_2^2*t_3^6+t_0^3*t_1*t_2^7*t_3+2*t_0^3*t_1*t_2^6*t_3^2+3*t_0^3*t_1*t_2^5*t_3^3+4*t_0^3*t_1*t_2^4*t_3^4+3*t_0^3*t_1*t_2^3*t_3^5+2*t_0^3*t_1*t_2^2*t_3^6+t_0^3*t_1*t_2*t_3^7+t_0^2*t_1^2*t_2^8+2*t_0^2*t_1^2*t_2^7*t_3+4*t_0^2*t_1^2*t_2^6*t_3^2+5*t_0^2*t_1^2*t_2^5*t_3^3+7*t_0^2*t_1^2*t_2^4*t_3^4+5*t_0^2*t_1^2*t_2^3*t_3^5+4*t_0^2*t_1^2*t_2^2*t_3^6+2*t_0^2*t_1^2*t_2*t_3^7+t_0^2*t_1^2*t_3^8+t_0*t_1^3*t_2^7*t_3+2*t_0*t_1^3*t_2^6*t_3^2+3*t_0*t_1^3*t_2^5*t_3^3+4*t_0*t_1^3*t_2^4*t_3^4+3*t_0*t_1^3*t_2^3*t_3^5+2*t_0*t_1^3*t_2^2*t_3^6+t_0*t_1^3*t_2*t_3^7+t_1^4*t_2^6*t_3^2+t_1^4*t_2^5*t_3^3+2*t_1^4*t_2^4*t_3^4+t_1^4*t_2^3*t_3^5+t_1^4*t_2^2*t_3^6, (3,0) => 0, (2,1) => t_0^6*t_2^8*t_3^4+2*t_0^6*t_2^7*t_3^5+2*t_0^6*t_2^6*t_3^6+2*t_0^6*t_2^5*t_3^7+t_0^6*t_2^4*t_3^8+t_0^5*t_1*t_2^10*t_3^2+3*t_0^5*t_1*t_2^9*t_3^3+6*t_0^5*t_1*t_2^8*t_3^4+9*t_0^5*t_1*t_2^7*t_3^5+10*t_0^5*t_1*t_2^6*t_3^6+9*t_0^5*t_1*t_2^5*t_3^7+6*t_0^5*t_1*t_2^4*t_3^8+3*t_0^5*t_1*t_2^3*t_3^9+t_0^5*t_1*t_2^2*t_3^10+t_0^4*t_1^2*t_2^11*t_3+4*t_0^4*t_1^2*t_2^10*t_3^2+9*t_0^4*t_1^2*t_2^9*t_3^3+16*t_0^4*t_1^2*t_2^8*t_3^4+22*t_0^4*t_1^2*t_2^7*t_3^5+24*t_0^4*t_1^2*t_2^6*t_3^6+22*t_0^4*t_1^2*t_2^5*t_3^7+16*t_0^4*t_1^2*t_2^4*t_3^8+9*t_0^4*t_1^2*t_2^3*t_3^9+4*t_0^4*t_1^2*t_2^2*t_3^10+t_0^4*t_1^2*t_2*t_3^11+2*t_0^3*t_1^3*t_2^11*t_3+6*t_0^3*t_1^3*t_2^10*t_3^2+12*t_0^3*t_1^3*t_2^9*t_3^3+21*t_0^3*t_1^3*t_2^8*t_3^4+28*t_0^3*t_1^3*t_2^7*t_3^5+30*t_0^3*t_1^3*t_2^6*t_3^6+28*t_0^3*t_1^3*t_2^5*t_3^7+21*t_0^3*t_1^3*t_2^4*t_3^8+12*t_0^3*t_1^3*t_2^3*t_3^9+6*t_0^3*t_1^3*t_2^2*t_3^10+2*t_0^3*t_1^3*t_2*t_3^11+t_0^2*t_1^4*t_2^11*t_3+4*t_0^2*t_1^4*t_2^10*t_3^2+9*t_0^2*t_1^4*t_2^9*t_3^3+16*t_0^2*t_1^4*t_2^8*t_3^4+22*t_0^2*t_1^4*t_2^7*t_3^5+24*t_0^2*t_1^4*t_2^6*t_3^6+22*t_0^2*t_1^4*t_2^5*t_3^7+16*t_0^2*t_1^4*t_2^4*t_3^8+9*t_0^2*t_1^4*t_2^3*t_3^9+4*t_0^2*t_1^4*t_2^2*t_3^10+t_0^2*t_1^4*t_2*t_3^11+t_0*t_1^5*t_2^10*t_3^2+3*t_0*t_1^5*t_2^9*t_3^3+6*t_0*t_1^5*t_2^8*t_3^4+9*t_0*t_1^5*t_2^7*t_3^5+10*t_0*t_1^5*t_2^6*t_3^6+9*t_0*t_1^5*t_2^5*t_3^7+6*t_0*t_1^5*t_2^4*t_3^8+3*t_0*t_1^5*t_2^3*t_3^9+t_0*t_1^5*t_2^2*t_3^10+t_1^6*t_2^8*t_3^4+2*t_1^6*t_2^7*t_3^5+2*t_1^6*t_2^6*t_3^6+2*t_1^6*t_2^5*t_3^7+t_1^6*t_2^4*t_3^8, (2,2) => 0, (4,0) => 0, (3,1) => t_0^8*t_2^9*t_3^7+t_0^8*t_2^8*t_3^8+t_0^8*t_2^7*t_3^9+t_0^7*t_1*t_2^12*t_3^4+3*t_0^7*t_1*t_2^11*t_3^5+6*t_0^7*t_1*t_2^10*t_3^6+9*t_0^7*t_1*t_2^9*t_3^7+10*t_0^7*t_1*t_2^8*t_3^8+9*t_0^7*t_1*t_2^7*t_3^9+6*t_0^7*t_1*t_2^6*t_3^10+3*t_0^7*t_1*t_2^5*t_3^11+t_0^7*t_1*t_2^4*t_3^12+2*t_0^6*t_1^2*t_2^13*t_3^3+6*t_0^6*t_1^2*t_2^12*t_3^4+15*t_0^6*t_1^2*t_2^11*t_3^5+24*t_0^6*t_1^2*t_2^10*t_3^6+34*t_0^6*t_1^2*t_2^9*t_3^7+36*t_0^6*t_1^2*t_2^8*t_3^8+34*t_0^6*t_1^2*t_2^7*t_3^9+24*t_0^6*t_1^2*t_2^6*t_3^10+15*t_0^6*t_1^2*t_2^5*t_3^11+6*t_0^6*t_1^2*t_2^4*t_3^12+2*t_0^6*t_1^2*t_2^3*t_3^13+2*t_0^5*t_1^3*t_2^14*t_3^2+7*t_0^5*t_1^3*t_2^13*t_3^3+18*t_0^5*t_1^3*t_2^12*t_3^4+35*t_0^5*t_1^3*t_2^11*t_3^5+54*t_0^5*t_1^3*t_2^10*t_3^6+70*t_0^5*t_1^3*t_2^9*t_3^7+76*t_0^5*t_1^3*t_2^8*t_3^8+70*t_0^5*t_1^3*t_2^7*t_3^9+54*t_0^5*t_1^3*t_2^6*t_3^10+35*t_0^5*t_1^3*t_2^5*t_3^11+18*t_0^5*t_1^3*t_2^4*t_3^12+7*t_0^5*t_1^3*t_2^3*t_3^13+2*t_0^5*t_1^3*t_2^2*t_3^14+2*t_0^4*t_1^4*t_2^14*t_3^2+9*t_0^4*t_1^4*t_2^13*t_3^3+22*t_0^4*t_1^4*t_2^12*t_3^4+44*t_0^4*t_1^4*t_2^11*t_3^5+66*t_0^4*t_1^4*t_2^10*t_3^6+87*t_0^4*t_1^4*t_2^9*t_3^7+93*t_0^4*t_1^4*t_2^8*t_3^8+87*t_0^4*t_1^4*t_2^7*t_3^9+66*t_0^4*t_1^4*t_2^6*t_3^10+44*t_0^4*t_1^4*t_2^5*t_3^11+22*t_0^4*t_1^4*t_2^4*t_3^12+9*t_0^4*t_1^4*t_2^3*t_3^13+2*t_0^4*t_1^4*t_2^2*t_3^14+2*t_0^3*t_1^5*t_2^14*t_3^2+7*t_0^3*t_1^5*t_2^13*t_3^3+18*t_0^3*t_1^5*t_2^12*t_3^4+35*t_0^3*t_1^5*t_2^11*t_3^5+54*t_0^3*t_1^5*t_2^10*t_3^6+70*t_0^3*t_1^5*t_2^9*t_3^7+76*t_0^3*t_1^5*t_2^8*t_3^8+70*t_0^3*t_1^5*t_2^7*t_3^9+54*t_0^3*t_1^5*t_2^6*t_3^10+35*t_0^3*t_1^5*t_2^5*t_3^11+18*t_0^3*t_1^5*t_2^4*t_3^12+7*t_0^3*t_1^5*t_2^3*t_3^13+2*t_0^3*t_1^5*t_2^2*t_3^14+2*t_0^2*t_1^6*t_2^13*t_3^3+6*t_0^2*t_1^6*t_2^12*t_3^4+15*t_0^2*t_1^6*t_2^11*t_3^5+24*t_0^2*t_1^6*t_2^10*t_3^6+34*t_0^2*t_1^6*t_2^9*t_3^7+36*t_0^2*t_1^6*t_2^8*t_3^8+34*t_0^2*t_1^6*t_2^7*t_3^9+24*t_0^2*t_1^6*t_2^6*t_3^10+15*t_0^2*t_1^6*t_2^5*t_3^11+6*t_0^2*t_1^6*t_2^4*t_3^12+2*t_0^2*t_1^6*t_2^3*t_3^13+t_0*t_1^7*t_2^12*t_3^4+3*t_0*t_1^7*t_2^11*t_3^5+6*t_0*t_1^7*t_2^10*t_3^6+9*t_0*t_1^7*t_2^9*t_3^7+10*t_0*t_1^7*t_2^8*t_3^8+9*t_0*t_1^7*t_2^7*t_3^9+6*t_0*t_1^7*t_2^6*t_3^10+3*t_0*t_1^7*t_2^5*t_3^11+t_0*t_1^7*t_2^4*t_3^12+t_1^8*t_2^9*t_3^7+t_1^8*t_2^8*t_3^8+t_1^8*t_2^7*t_3^9, (3,2) => 0, (5,0) => 0, (4,1) => t_0^9*t_1*t_2^13*t_3^7+2*t_0^9*t_1*t_2^12*t_3^8+3*t_0^9*t_1*t_2^11*t_3^9+4*t_0^9*t_1*t_2^10*t_3^10+3*t_0^9*t_1*t_2^9*t_3^11+2*t_0^9*t_1*t_2^8*t_3^12+t_0^9*t_1*t_2^7*t_3^13+t_0^8*t_1^2*t_2^15*t_3^5+4*t_0^8*t_1^2*t_2^14*t_3^6+10*t_0^8*t_1^2*t_2^13*t_3^7+17*t_0^8*t_1^2*t_2^12*t_3^8+23*t_0^8*t_1^2*t_2^11*t_3^9+26*t_0^8*t_1^2*t_2^10*t_3^10+23*t_0^8*t_1^2*t_2^9*t_3^11+17*t_0^8*t_1^2*t_2^8*t_3^12+10*t_0^8*t_1^2*t_2^7*t_3^13+4*t_0^8*t_1^2*t_2^6*t_3^14+t_0^8*t_1^2*t_2^5*t_3^15+2*t_0^7*t_1^3*t_2^16*t_3^4+8*t_0^7*t_1^3*t_2^15*t_3^5+20*t_0^7*t_1^3*t_2^14*t_3^6+39*t_0^7*t_1^3*t_2^13*t_3^7+60*t_0^7*t_1^3*t_2^12*t_3^8+77*t_0^7*t_1^3*t_2^11*t_3^9+84*t_0^7*t_1^3*t_2^10*t_3^10+77*t_0^7*t_1^3*t_2^9*t_3^11+60*t_0^7*t_1^3*t_2^8*t_3^12+39*t_0^7*t_1^3*t_2^7*t_3^13+20*t_0^7*t_1^3*t_2^6*t_3^14+8*t_0^7*t_1^3*t_2^5*t_3^15+2*t_0^7*t_1^3*t_2^4*t_3^16+t_0^6*t_1^4*t_2^17*t_3^3+6*t_0^6*t_1^4*t_2^16*t_3^4+19*t_0^6*t_1^4*t_2^15*t_3^5+43*t_0^6*t_1^4*t_2^14*t_3^6+78*t_0^6*t_1^4*t_2^13*t_3^7+116*t_0^6*t_1^4*t_2^12*t_3^8+146*t_0^6*t_1^4*t_2^11*t_3^9+158*t_0^6*t_1^4*t_2^10*t_3^10+146*t_0^6*t_1^4*t_2^9*t_3^11+116*t_0^6*t_1^4*t_2^8*t_3^12+78*t_0^6*t_1^4*t_2^7*t_3^13+43*t_0^6*t_1^4*t_2^6*t_3^14+19*t_0^6*t_1^4*t_2^5*t_3^15+6*t_0^6*t_1^4*t_2^4*t_3^16+t_0^6*t_1^4*t_2^3*t_3^17+2*t_0^5*t_1^5*t_2^17*t_3^3+9*t_0^5*t_1^5*t_2^16*t_3^4+26*t_0^5*t_1^5*t_2^15*t_3^5+57*t_0^5*t_1^5*t_2^14*t_3^6+100*t_0^5*t_1^5*t_2^13*t_3^7+146*t_0^5*t_1^5*t_2^12*t_3^8+182*t_0^5*t_1^5*t_2^11*t_3^9+196*t_0^5*t_1^5*t_2^10*t_3^10+182*t_0^5*t_1^5*t_2^9*t_3^11+146*t_0^5*t_1^5*t_2^8*t_3^12+100*t_0^5*t_1^5*t_2^7*t_3^13+57*t_0^5*t_1^5*t_2^6*t_3^14+26*t_0^5*t_1^5*t_2^5*t_3^15+9*t_0^5*t_1^5*t_2^4*t_3^16+2*t_0^5*t_1^5*t_2^3*t_3^17+t_0^4*t_1^6*t_2^17*t_3^3+6*t_0^4*t_1^6*t_2^16*t_3^4+19*t_0^4*t_1^6*t_2^15*t_3^5+43*t_0^4*t_1^6*t_2^14*t_3^6+78*t_0^4*t_1^6*t_2^13*t_3^7+116*t_0^4*t_1^6*t_2^12*t_3^8+146*t_0^4*t_1^6*t_2^11*t_3^9+158*t_0^4*t_1^6*t_2^10*t_3^10+146*t_0^4*t_1^6*t_2^9*t_3^11+116*t_0^4*t_1^6*t_2^8*t_3^12+78*t_0^4*t_1^6*t_2^7*t_3^13+43*t_0^4*t_1^6*t_2^6*t_3^14+19*t_0^4*t_1^6*t_2^5*t_3^15+6*t_0^4*t_1^6*t_2^4*t_3^16+t_0^4*t_1^6*t_2^3*t_3^17+2*t_0^3*t_1^7*t_2^16*t_3^4+8*t_0^3*t_1^7*t_2^15*t_3^5+20*t_0^3*t_1^7*t_2^14*t_3^6+39*t_0^3*t_1^7*t_2^13*t_3^7+60*t_0^3*t_1^7*t_2^12*t_3^8+77*t_0^3*t_1^7*t_2^11*t_3^9+84*t_0^3*t_1^7*t_2^10*t_3^10+77*t_0^3*t_1^7*t_2^9*t_3^11+60*t_0^3*t_1^7*t_2^8*t_3^12+39*t_0^3*t_1^7*t_2^7*t_3^13+20*t_0^3*t_1^7*t_2^6*t_3^14+8*t_0^3*t_1^7*t_2^5*t_3^15+2*t_0^3*t_1^7*t_2^4*t_3^16+t_0^2*t_1^8*t_2^15*t_3^5+4*t_0^2*t_1^8*t_2^14*t_3^6+10*t_0^2*t_1^8*t_2^13*t_3^7+17*t_0^2*t_1^8*t_2^12*t_3^8+23*t_0^2*t_1^8*t_2^11*t_3^9+26*t_0^2*t_1^8*t_2^10*t_3^10+23*t_0^2*t_1^8*t_2^9*t_3^11+17*t_0^2*t_1^8*t_2^8*t_3^12+10*t_0^2*t_1^8*t_2^7*t_3^13+4*t_0^2*t_1^8*t_2^6*t_3^14+t_0^2*t_1^8*t_2^5*t_3^15+t_0*t_1^9*t_2^13*t_3^7+2*t_0*t_1^9*t_2^12*t_3^8+3*t_0*t_1^9*t_2^11*t_3^9+4*t_0*t_1^9*t_2^10*t_3^10+3*t_0*t_1^9*t_2^9*t_3^11+2*t_0*t_1^9*t_2^8*t_3^12+t_0*t_1^9*t_2^7*t_3^13, (4,2) => 0, (6,0) => 0, (5,1) => t_0^10*t_1^2*t_2^16*t_3^8+2*t_0^10*t_1^2*t_2^15*t_3^9+5*t_0^10*t_1^2*t_2^14*t_3^10+6*t_0^10*t_1^2*t_2^13*t_3^11+7*t_0^10*t_1^2*t_2^12*t_3^12+6*t_0^10*t_1^2*t_2^11*t_3^13+5*t_0^10*t_1^2*t_2^10*t_3^14+2*t_0^10*t_1^2*t_2^9*t_3^15+t_0^10*t_1^2*t_2^8*t_3^16+3*t_0^9*t_1^3*t_2^17*t_3^7+9*t_0^9*t_1^3*t_2^16*t_3^8+18*t_0^9*t_1^3*t_2^15*t_3^9+30*t_0^9*t_1^3*t_2^14*t_3^10+39*t_0^9*t_1^3*t_2^13*t_3^11+42*t_0^9*t_1^3*t_2^12*t_3^12+39*t_0^9*t_1^3*t_2^11*t_3^13+30*t_0^9*t_1^3*t_2^10*t_3^14+18*t_0^9*t_1^3*t_2^9*t_3^15+9*t_0^9*t_1^3*t_2^8*t_3^16+3*t_0^9*t_1^3*t_2^7*t_3^17+t_0^8*t_1^4*t_2^19*t_3^5+5*t_0^8*t_1^4*t_2^18*t_3^6+15*t_0^8*t_1^4*t_2^17*t_3^7+36*t_0^8*t_1^4*t_2^16*t_3^8+63*t_0^8*t_1^4*t_2^15*t_3^9+96*t_0^8*t_1^4*t_2^14*t_3^10+119*t_0^8*t_1^4*t_2^13*t_3^11+130*t_0^8*t_1^4*t_2^12*t_3^12+119*t_0^8*t_1^4*t_2^11*t_3^13+96*t_0^8*t_1^4*t_2^10*t_3^14+63*t_0^8*t_1^4*t_2^9*t_3^15+36*t_0^8*t_1^4*t_2^8*t_3^16+15*t_0^8*t_1^4*t_2^7*t_3^17+5*t_0^8*t_1^4*t_2^6*t_3^18+t_0^8*t_1^4*t_2^5*t_3^19+3*t_0^7*t_1^5*t_2^19*t_3^5+12*t_0^7*t_1^5*t_2^18*t_3^6+34*t_0^7*t_1^5*t_2^17*t_3^7+72*t_0^7*t_1^5*t_2^16*t_3^8+123*t_0^7*t_1^5*t_2^15*t_3^9+178*t_0^7*t_1^5*t_2^14*t_3^10+220*t_0^7*t_1^5*t_2^13*t_3^11+236*t_0^7*t_1^5*t_2^12*t_3^12+220*t_0^7*t_1^5*t_2^11*t_3^13+178*t_0^7*t_1^5*t_2^10*t_3^14+123*t_0^7*t_1^5*t_2^9*t_3^15+72*t_0^7*t_1^5*t_2^8*t_3^16+34*t_0^7*t_1^5*t_2^7*t_3^17+12*t_0^7*t_1^5*t_2^6*t_3^18+3*t_0^7*t_1^5*t_2^5*t_3^19+t_0^6*t_1^6*t_2^20*t_3^4+5*t_0^6*t_1^6*t_2^19*t_3^5+18*t_0^6*t_1^6*t_2^18*t_3^6+45*t_0^6*t_1^6*t_2^17*t_3^7+93*t_0^6*t_1^6*t_2^16*t_3^8+153*t_0^6*t_1^6*t_2^15*t_3^9+222*t_0^6*t_1^6*t_2^14*t_3^10+270*t_0^6*t_1^6*t_2^13*t_3^11+291*t_0^6*t_1^6*t_2^12*t_3^12+270*t_0^6*t_1^6*t_2^11*t_3^13+222*t_0^6*t_1^6*t_2^10*t_3^14+153*t_0^6*t_1^6*t_2^9*t_3^15+93*t_0^6*t_1^6*t_2^8*t_3^16+45*t_0^6*t_1^6*t_2^7*t_3^17+18*t_0^6*t_1^6*t_2^6*t_3^18+5*t_0^6*t_1^6*t_2^5*t_3^19+t_0^6*t_1^6*t_2^4*t_3^20+3*t_0^5*t_1^7*t_2^19*t_3^5+12*t_0^5*t_1^7*t_2^18*t_3^6+34*t_0^5*t_1^7*t_2^17*t_3^7+72*t_0^5*t_1^7*t_2^16*t_3^8+123*t_0^5*t_1^7*t_2^15*t_3^9+178*t_0^5*t_1^7*t_2^14*t_3^10+220*t_0^5*t_1^7*t_2^13*t_3^11+236*t_0^5*t_1^7*t_2^12*t_3^12+220*t_0^5*t_1^7*t_2^11*t_3^13+178*t_0^5*t_1^7*t_2^10*t_3^14+123*t_0^5*t_1^7*t_2^9*t_3^15+72*t_0^5*t_1^7*t_2^8*t_3^16+34*t_0^5*t_1^7*t_2^7*t_3^17+12*t_0^5*t_1^7*t_2^6*t_3^18+3*t_0^5*t_1^7*t_2^5*t_3^19+t_0^4*t_1^8*t_2^19*t_3^5+5*t_0^4*t_1^8*t_2^18*t_3^6+15*t_0^4*t_1^8*t_2^17*t_3^7+36*t_0^4*t_1^8*t_2^16*t_3^8+63*t_0^4*t_1^8*t_2^15*t_3^9+96*t_0^4*t_1^8*t_2^14*t_3^10+119*t_0^4*t_1^8*t_2^13*t_3^11+130*t_0^4*t_1^8*t_2^12*t_3^12+119*t_0^4*t_1^8*t_2^11*t_3^13+96*t_0^4*t_1^8*t_2^10*t_3^14+63*t_0^4*t_1^8*t_2^9*t_3^15+36*t_0^4*t_1^8*t_2^8*t_3^16+15*t_0^4*t_1^8*t_2^7*t_3^17+5*t_0^4*t_1^8*t_2^6*t_3^18+t_0^4*t_1^8*t_2^5*t_3^19+3*t_0^3*t_1^9*t_2^17*t_3^7+9*t_0^3*t_1^9*t_2^16*t_3^8+18*t_0^3*t_1^9*t_2^15*t_3^9+30*t_0^3*t_1^9*t_2^14*t_3^10+39*t_0^3*t_1^9*t_2^13*t_3^11+42*t_0^3*t_1^9*t_2^12*t_3^12+39*t_0^3*t_1^9*t_2^11*t_3^13+30*t_0^3*t_1^9*t_2^10*t_3^14+18*t_0^3*t_1^9*t_2^9*t_3^15+9*t_0^3*t_1^9*t_2^8*t_3^16+3*t_0^3*t_1^9*t_2^7*t_3^17+t_0^2*t_1^10*t_2^16*t_3^8+2*t_0^2*t_1^10*t_2^15*t_3^9+5*t_0^2*t_1^10*t_2^14*t_3^10+6*t_0^2*t_1^10*t_2^13*t_3^11+7*t_0^2*t_1^10*t_2^12*t_3^12+6*t_0^2*t_1^10*t_2^11*t_3^13+5*t_0^2*t_1^10*t_2^10*t_3^14+2*t_0^2*t_1^10*t_2^9*t_3^15+t_0^2*t_1^10*t_2^8*t_3^16};
end;