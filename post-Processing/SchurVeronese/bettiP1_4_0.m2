--This file computes Betti tables for P^2 for d = 4 and b = 0
A := QQ[t_0,t_1,t_2];
--tb stands for Total Betti numbers
tb40 = new HashTable from {(0,0) => 1, (1,0) => 0, (0,1) => 0, (1,1) => 6, (2,0) => 0, (3,0) => 0, (2,1) => 8, (3,1) => 3};
--mb stands for Multigraded Betti numbers
mb40 = new HashTable from {(0,0) => 1, (0,1) => 0, (1,0) => 0, (1,1) => t_0^6*t_1^2+t_0^5*t_1^3+2*t_0^4*t_1^4+t_0^3*t_1^5+t_0^2*t_1^6, (2,0) => 0, (2,1) => t_0^8*t_1^4+2*t_0^7*t_1^5+2*t_0^6*t_1^6+2*t_0^5*t_1^7+t_0^4*t_1^8, (3,0) => 0, (3,1) => t_0^9*t_1^7+t_0^8*t_1^8+t_0^7*t_1^9};
--sb represents the betti numbers as sums of Schur functors
sb40 = new HashTable from {(0,0) => {({0,0},1)}, (1,0) => {}, (0,1) => {}, (1,1) => {({6,2},1)}, (2,0) => {}, (3,0) => {}, (2,1) => {({8,4},1)}, (3,1) => {({9,7},1)}};
--dw encodes the dominant weights in each entry
dw40 = new HashTable from {(0,0) => {{0,0}}, (1,0) => {}, (0,1) => {}, (1,1) => {{6,2}}, (2,0) => {}, (3,0) => {}, (2,1) => {{8,4}}, (3,1) => {{9,7}}};
--lw encodes the lex leading weight in each entry
lw40 = new HashTable from {(0,0) => {0,0}, (1,0) => {}, (0,1) => {}, (1,1) => {6,2}, (2,0) => {}, (3,0) => {}, (2,1) => {8,4}, (3,1) => {9,7}};
--nr encodes the number of disctinct reprsentations in each entry
nr40 = new HashTable from {(0,0) => 1, (1,0) => 0, (0,1) => 0, (1,1) => 1, (2,0) => 0, (3,0) => 0, (2,1) => 1, (3,1) => 1};
--nrm encodes the number of representations with multiplicity in each entry
nrm40 = new HashTable from {(0,0) => 1, (1,0) => 0, (0,1) => 0, (1,1) => 1, (2,0) => 0, (3,0) => 0, (2,1) => 1, (3,1) => 1};
--er encodes the errors in the computed multigraded Hilbert series via our Schur method in each entry
er40 = new HashTable from {(0,0) => 0, (1,0) => 0, (0,1) => 0, (1,1) => 0, (2,0) => 0, (3,0) => 0, (2,1) => 0, (3,1) => 0};
--bs encodes the Boij-Soederberg coefficients each entry
bs40 = {24/1};
end;