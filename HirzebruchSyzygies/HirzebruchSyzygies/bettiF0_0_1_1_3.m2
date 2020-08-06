A := QQ[t_0,t_1,t_2,t_3];
--tb stands for Total Betti numbers
tb0113 = new HashTable from {(5,2) => 0, (6,1) => 0, (7,1) => 0, (0,0) => 2/1, (1,0) => 6/1, (1,1) => 0, (2,0) => 0, (2,1) => 20/1, (3,0) => 0, (2,2) => 0, (3,1) => 30/1, (4,0) => 0, (4,1) => 18/1, (5,0) => 0, (3,2) => 0, (5,1) => 4/1, (4,2) => 0};
--sb represents the betti numbers as sums of Schur functors
sb0113 = new HashTable from {(5,2) => {}, (6,1) => {}, (7,1) => {}, (0,0) => {({0,0,1,0},1/1)}, (1,0) => {({1,0,3,1},1/1)}, (1,1) => {}, (2,0) => {}, (2,1) => {({3,0,5,5},1/1),({2,1,7,3},1/1),({2,1,6,4},1/1)}, (3,0) => {}, (2,2) => {}, (3,1) => {({3,1,8,5},1/1),({3,1,7,6},1/1),({2,2,9,4},1/1),({2,2,8,5},1/1),({2,2,7,6},1/1)}, (4,0) => {}, (4,1) => {({3,2,10,6},1/1),({3,2,9,7},1/1),({3,2,8,8},1/1)}, (5,0) => {}, (3,2) => {}, (5,1) => {({3,3,11,8},1/1)}, (4,2) => {}};
--dw stands for dominant weights
dw0113 = new HashTable from {(5,2) => {}, (6,1) => {}, (7,1) => {}, (0,0) => {({0,0,1,0},1/1)}, (1,0) => {({1,0,3,1},1/1)}, (1,1) => {}, (2,0) => {}, (2,1) => {({3,0,5,5},1/1),({2,1,7,3},1/1)}, (3,0) => {}, (2,2) => {}, (3,1) => {({3,1,8,5},1/1),({2,2,9,4},1/1)}, (4,0) => {}, (3,2) => {}, (5,0) => {}, (4,1) => {({3,2,10,6},1/1)}, (4,2) => {}, (5,1) => {({3,3,11,8},1/1)}};
--dmw stands for dominant monomial weights
dmw0113 = new HashTable from {(5,2) => {{3,4,11,11}}, (6,1) => {}, (7,1) => {}, (0,0) => {}, (1,0) => {}, (1,1) => {}, (2,0) => {}, (2,1) => {{3,0,5,5}}, (3,0) => {}, (2,2) => {}, (3,1) => {{3,1,8,5}}, (4,0) => {}, (3,2) => {}, (5,0) => {}, (4,1) => {{3,2,10,6}}, (4,2) => {{3,3,8,11}}, (5,1) => {{3,3,11,8}}};
end;