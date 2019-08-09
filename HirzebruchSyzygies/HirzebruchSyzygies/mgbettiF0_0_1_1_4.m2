A := QQ[t_0,t_1,t_2,t_3];
--mb stands for Multigraded Betti numbers
mb0114 = new HashTable from {(5,2) => 0, (6,1) => t_0^4*t_1^3*t_2^18*t_3^11+2*t_0^4*t_1^3*t_2^17*t_3^12+3*t_0^4*t_1^3*t_2^16*t_3^13+4*t_0^4*t_1^3*t_2^15*t_3^14+4*t_0^4*t_1^3*t_2^14*t_3^15+3*t_0^4*t_1^3*t_2^13*t_3^16+2*t_0^4*t_1^3*t_2^12*t_3^17+t_0^4*t_1^3*t_2^11*t_3^18+t_0^3*t_1^4*t_2^18*t_3^11+2*t_0^3*t_1^4*t_2^17*t_3^12+3*t_0^3*t_1^4*t_2^16*t_3^13+4*t_0^3*t_1^4*t_2^15*t_3^14+4*t_0^3*t_1^4*t_2^14*t_3^15+3*t_0^3*t_1^4*t_2^13*t_3^16+2*t_0^3*t_1^4*t_2^12*t_3^17+t_0^3*t_1^4*t_2^11*t_3^18, (7,0) => 0, (6,2) => 0, (7,1) => t_0^4*t_1^4*t_2^19*t_3^14+t_0^4*t_1^4*t_2^18*t_3^15+t_0^4*t_1^4*t_2^17*t_3^16+t_0^4*t_1^4*t_2^16*t_3^17+t_0^4*t_1^4*t_2^15*t_3^18+t_0^4*t_1^4*t_2^14*t_3^19, (8,1) => 0, (7,2) => 0, (9,1) => 0, (0,0) => t_2+t_3, (1,0) => t_0*t_2^4*t_3+t_0*t_2^3*t_3^2+t_0*t_2^2*t_3^3+t_0*t_2*t_3^4+t_1*t_2^4*t_3+t_1*t_2^3*t_3^2+t_1*t_2^2*t_3^3+t_1*t_2*t_3^4, (2,0) => 0, (1,1) => 0, (3,0) => 0, (2,1) => t_0^3*t_2^8*t_3^5+t_0^3*t_2^7*t_3^6+t_0^3*t_2^6*t_3^7+t_0^3*t_2^5*t_3^8+t_0^2*t_1*t_2^10*t_3^3+2*t_0^2*t_1*t_2^9*t_3^4+4*t_0^2*t_1*t_2^8*t_3^5+5*t_0^2*t_1*t_2^7*t_3^6+5*t_0^2*t_1*t_2^6*t_3^7+4*t_0^2*t_1*t_2^5*t_3^8+2*t_0^2*t_1*t_2^4*t_3^9+t_0^2*t_1*t_2^3*t_3^10+t_0*t_1^2*t_2^10*t_3^3+2*t_0*t_1^2*t_2^9*t_3^4+4*t_0*t_1^2*t_2^8*t_3^5+5*t_0*t_1^2*t_2^7*t_3^6+5*t_0*t_1^2*t_2^6*t_3^7+4*t_0*t_1^2*t_2^5*t_3^8+2*t_0*t_1^2*t_2^4*t_3^9+t_0*t_1^2*t_2^3*t_3^10+t_1^3*t_2^8*t_3^5+t_1^3*t_2^7*t_3^6+t_1^3*t_2^6*t_3^7+t_1^3*t_2^5*t_3^8, (2,2) => 0, (4,0) => 0, (3,1) => t_0^4*t_2^9*t_3^8+t_0^4*t_2^8*t_3^9+t_0^3*t_1*t_2^12*t_3^5+3*t_0^3*t_1*t_2^11*t_3^6+5*t_0^3*t_1*t_2^10*t_3^7+7*t_0^3*t_1*t_2^9*t_3^8+7*t_0^3*t_1*t_2^8*t_3^9+5*t_0^3*t_1*t_2^7*t_3^10+3*t_0^3*t_1*t_2^6*t_3^11+t_0^3*t_1*t_2^5*t_3^12+t_0^2*t_1^2*t_2^13*t_3^4+3*t_0^2*t_1^2*t_2^12*t_3^5+7*t_0^2*t_1^2*t_2^11*t_3^6+11*t_0^2*t_1^2*t_2^10*t_3^7+14*t_0^2*t_1^2*t_2^9*t_3^8+14*t_0^2*t_1^2*t_2^8*t_3^9+11*t_0^2*t_1^2*t_2^7*t_3^10+7*t_0^2*t_1^2*t_2^6*t_3^11+3*t_0^2*t_1^2*t_2^5*t_3^12+t_0^2*t_1^2*t_2^4*t_3^13+t_0*t_1^3*t_2^12*t_3^5+3*t_0*t_1^3*t_2^11*t_3^6+5*t_0*t_1^3*t_2^10*t_3^7+7*t_0*t_1^3*t_2^9*t_3^8+7*t_0*t_1^3*t_2^8*t_3^9+5*t_0*t_1^3*t_2^7*t_3^10+3*t_0*t_1^3*t_2^6*t_3^11+t_0*t_1^3*t_2^5*t_3^12+t_1^4*t_2^9*t_3^8+t_1^4*t_2^8*t_3^9, (3,2) => 0, (5,0) => 0, (4,1) => t_0^4*t_1*t_2^13*t_3^8+2*t_0^4*t_1*t_2^12*t_3^9+3*t_0^4*t_1*t_2^11*t_3^10+3*t_0^4*t_1*t_2^10*t_3^11+2*t_0^4*t_1*t_2^9*t_3^12+t_0^4*t_1*t_2^8*t_3^13+t_0^3*t_1^2*t_2^15*t_3^6+3*t_0^3*t_1^2*t_2^14*t_3^7+7*t_0^3*t_1^2*t_2^13*t_3^8+11*t_0^3*t_1^2*t_2^12*t_3^9+14*t_0^3*t_1^2*t_2^11*t_3^10+14*t_0^3*t_1^2*t_2^10*t_3^11+11*t_0^3*t_1^2*t_2^9*t_3^12+7*t_0^3*t_1^2*t_2^8*t_3^13+3*t_0^3*t_1^2*t_2^7*t_3^14+t_0^3*t_1^2*t_2^6*t_3^15+t_0^2*t_1^3*t_2^15*t_3^6+3*t_0^2*t_1^3*t_2^14*t_3^7+7*t_0^2*t_1^3*t_2^13*t_3^8+11*t_0^2*t_1^3*t_2^12*t_3^9+14*t_0^2*t_1^3*t_2^11*t_3^10+14*t_0^2*t_1^3*t_2^10*t_3^11+11*t_0^2*t_1^3*t_2^9*t_3^12+7*t_0^2*t_1^3*t_2^8*t_3^13+3*t_0^2*t_1^3*t_2^7*t_3^14+t_0^2*t_1^3*t_2^6*t_3^15+t_0*t_1^4*t_2^13*t_3^8+2*t_0*t_1^4*t_2^12*t_3^9+3*t_0*t_1^4*t_2^11*t_3^10+3*t_0*t_1^4*t_2^10*t_3^11+2*t_0*t_1^4*t_2^9*t_3^12+t_0*t_1^4*t_2^8*t_3^13, (4,2) => 0, (5,1) => t_0^4*t_1^2*t_2^16*t_3^9+2*t_0^4*t_1^2*t_2^15*t_3^10+4*t_0^4*t_1^2*t_2^14*t_3^11+5*t_0^4*t_1^2*t_2^13*t_3^12+5*t_0^4*t_1^2*t_2^12*t_3^13+4*t_0^4*t_1^2*t_2^11*t_3^14+2*t_0^4*t_1^2*t_2^10*t_3^15+t_0^4*t_1^2*t_2^9*t_3^16+t_0^3*t_1^3*t_2^17*t_3^8+3*t_0^3*t_1^3*t_2^16*t_3^9+6*t_0^3*t_1^3*t_2^15*t_3^10+10*t_0^3*t_1^3*t_2^14*t_3^11+12*t_0^3*t_1^3*t_2^13*t_3^12+12*t_0^3*t_1^3*t_2^12*t_3^13+10*t_0^3*t_1^3*t_2^11*t_3^14+6*t_0^3*t_1^3*t_2^10*t_3^15+3*t_0^3*t_1^3*t_2^9*t_3^16+t_0^3*t_1^3*t_2^8*t_3^17+t_0^2*t_1^4*t_2^16*t_3^9+2*t_0^2*t_1^4*t_2^15*t_3^10+4*t_0^2*t_1^4*t_2^14*t_3^11+5*t_0^2*t_1^4*t_2^13*t_3^12+5*t_0^2*t_1^4*t_2^12*t_3^13+4*t_0^2*t_1^4*t_2^11*t_3^14+2*t_0^2*t_1^4*t_2^10*t_3^15+t_0^2*t_1^4*t_2^9*t_3^16, (6,0) => 0};
end;