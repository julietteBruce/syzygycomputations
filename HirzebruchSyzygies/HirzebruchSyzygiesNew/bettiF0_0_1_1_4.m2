A := QQ[t_0,t_1,t_2,t_3];
--tb stands for Total Betti numbers
tb0114 = new HashTable from {(7,0) => 0, (6,1) => 40/1, (5,2) => 0, (7,1) => 6/1, (6,2) => 0, (7,2) => 0, (8,1) => 0, (9,1) => 0, (0,0) => 2/1, (1,0) => 8/1, (1,1) => 0, (2,0) => 0, (2,1) => 56/1, (3,0) => 0, (3,1) => 140/1, (4,0) => 0, (2,2) => 0, (4,1) => 168/1, (5,0) => 0, (3,2) => 0, (6,0) => 0, (5,1) => 112/1, (4,2) => 0};
--sb represents the betti numbers as sums of Schur functors
sb0114 = new HashTable from {(7,0) => {}, (6,1) => {({4,3,18,11},1/1),({4,3,17,12},1/1),({4,3,16,13},1/1),({4,3,15,14},1/1)}, (5,2) => {}, (7,1) => {({4,4,19,14},1/1)}, (6,2) => {}, (7,2) => {}, (8,1) => {}, (9,1) => {}, (0,0) => {({0,0,1,0},1/1)}, (1,0) => {({1,0,4,1},1/1)}, (1,1) => {}, (2,0) => {}, (2,1) => {({3,0,8,5},1/1),({2,1,10,3},1/1),({2,1,9,4},1/1),({2,1,8,5},1/1),({2,1,7,6},1/1)}, (3,0) => {}, (3,1) => {({4,0,9,8},1/1),({3,1,12,5},1/1),({3,1,11,6},2/1),({3,1,10,7},2/1),({3,1,9,8},1/1),({2,2,13,4},1/1),({2,2,12,5},1/1),({2,2,11,6},2/1),({2,2,10,7},2/1),({2,2,9,8},1/1)}, (4,0) => {}, (2,2) => {}, (4,1) => {({4,1,13,8},1/1),({4,1,12,9},1/1),({4,1,11,10},1/1),({3,2,15,6},1/1),({3,2,14,7},2/1),({3,2,13,8},3/1),({3,2,12,9},3/1),({3,2,11,10},2/1)}, (5,0) => {}, (3,2) => {}, (6,0) => {}, (5,1) => {({4,2,16,9},1/1),({4,2,15,10},1/1),({4,2,14,11},2/1),({4,2,13,12},1/1),({3,3,17,8},1/1),({3,3,16,9},1/1),({3,3,15,10},2/1),({3,3,14,11},2/1),({3,3,13,12},1/1)}, (4,2) => {}};
--dw stands for dominant weights
dw0114 = new HashTable from {(7,0) => {}, (6,1) => {({4,3,18,11},1/1)}, (5,2) => {}, (6,2) => {}, (7,1) => {({4,4,19,14},1/1)}, (7,2) => {}, (8,1) => {}, (9,1) => {}, (0,0) => {({0,0,1,0},1/1)}, (1,0) => {({1,0,4,1},1/1)}, (1,1) => {}, (2,0) => {}, (3,0) => {}, (2,1) => {({3,0,8,5},1/1),({2,1,10,3},1/1)}, (2,2) => {}, (4,0) => {}, (3,1) => {({4,0,9,8},1/1),({3,1,12,5},1/1),({2,2,13,4},1/1)}, (3,2) => {}, (5,0) => {}, (4,1) => {({4,1,13,8},1/1),({3,2,15,6},1/1)}, (4,2) => {}, (5,1) => {({4,2,16,9},1/1),({3,3,17,8},1/1)}, (6,0) => {}};
--dmw stands for dominant monomial weights
dmw0114 = new HashTable from {(7,0) => {}, (6,1) => {{5,2,14,15},{4,3,15,14}}, (5,2) => {}, (6,2) => {}, (7,1) => {{5,3,15,18}}, (7,2) => {}, (8,1) => {}, (9,1) => {}, (0,0) => {}, (1,0) => {}, (1,1) => {}, (2,0) => {}, (3,0) => {}, (2,1) => {{3,0,8,5}}, (2,2) => {}, (4,0) => {}, (3,1) => {{4,0,9,8},{3,1,11,6}}, (3,2) => {}, (5,0) => {}, (4,1) => {{5,0,9,12},{3,2,13,8},{4,1,12,9}}, (4,2) => {}, (5,1) => {{5,1,12,13},{4,2,14,11}}, (6,0) => {}};
end;