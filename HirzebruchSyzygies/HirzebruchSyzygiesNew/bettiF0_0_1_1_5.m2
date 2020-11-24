A := QQ[t_0,t_1,t_2,t_3];
--tb stands for Total Betti numbers
tb0115 = new HashTable from {(5,2) => 0, (6,1) => 600/1, (7,0) => 0, (6,2) => 0, (7,1) => 270/1, (8,0) => 0, (8,1) => 70/1, (9,0) => 0, (7,2) => 0, (9,1) => 8/1, (8,2) => 0, (9,2) => 0, (10,1) => 0, (11,1) => 0, (0,0) => 2/1, (1,0) => 10/1, (2,0) => 0, (1,1) => 0, (3,0) => 0, (2,1) => 120/1, (2,2) => 0, (3,1) => 420/1, (4,0) => 0, (3,2) => 0, (4,1) => 756/1, (5,0) => 0, (4,2) => 0, (5,1) => 840/1, (6,0) => 0};
--sb represents the betti numbers as sums of Schur functors
sb0115 = new HashTable from {(5,2) => {}, (6,1) => {({5,2,23,13},1/1),({5,2,22,14},1/1),({5,2,21,15},2/1),({5,2,20,16},2/1),({5,2,19,17},2/1),({4,3,25,11},1/1),({4,3,24,12},2/1),({4,3,23,13},3/1),({4,3,22,14},5/1),({4,3,21,15},6/1),({4,3,20,16},5/1),({4,3,19,17},4/1),({4,3,18,18},2/1)}, (7,0) => {}, (6,2) => {}, (7,1) => {({5,3,26,15},1/1),({5,3,25,16},1/1),({5,3,24,17},2/1),({5,3,23,18},2/1),({5,3,22,19},2/1),({5,3,21,20},1/1),({4,4,27,14},1/1),({4,4,26,15},1/1),({4,4,25,16},2/1),({4,4,24,17},2/1),({4,4,23,18},3/1),({4,4,22,19},2/1),({4,4,21,20},1/1)}, (8,0) => {}, (8,1) => {({5,4,28,18},1/1),({5,4,27,19},1/1),({5,4,26,20},1/1),({5,4,25,21},1/1),({5,4,24,22},1/1)}, (9,0) => {}, (7,2) => {}, (9,1) => {({5,5,29,22},1/1)}, (8,2) => {}, (9,2) => {}, (10,1) => {}, (11,1) => {}, (0,0) => {({0,0,1,0},1/1)}, (1,0) => {({1,0,5,1},1/1)}, (2,0) => {}, (1,1) => {}, (3,0) => {}, (2,1) => {({3,0,11,5},1/1),({3,0,9,7},1/1),({2,1,13,3},1/1),({2,1,12,4},1/1),({2,1,11,5},1/1),({2,1,10,6},2/1),({2,1,9,7},1/1)}, (2,2) => {}, (3,1) => {({4,0,13,8},1/1),({4,0,12,9},1/1),({3,1,16,5},1/1),({3,1,15,6},2/1),({3,1,14,7},3/1),({3,1,13,8},3/1),({3,1,12,9},3/1),({3,1,11,10},2/1),({2,2,17,4},1/1),({2,2,16,5},1/1),({2,2,15,6},2/1),({2,2,14,7},3/1),({2,2,13,8},3/1),({2,2,12,9},2/1),({2,2,11,10},2/1)}, (4,0) => {}, (3,2) => {}, (4,1) => {({5,0,14,12},1/1),({4,1,18,8},1/1),({4,1,17,9},2/1),({4,1,16,10},3/1),({4,1,15,11},3/1),({4,1,14,12},2/1),({4,1,13,13},1/1),({3,2,20,6},1/1),({3,2,19,7},2/1),({3,2,18,8},4/1),({3,2,17,9},5/1),({3,2,16,10},7/1),({3,2,15,11},6/1),({3,2,14,12},5/1),({3,2,13,13},1/1)}, (5,0) => {}, (4,2) => {}, (5,1) => {({5,1,19,12},1/1),({5,1,18,13},1/1),({5,1,17,14},1/1),({5,1,16,15},1/1),({4,2,22,9},1/1),({4,2,21,10},2/1),({4,2,20,11},4/1),({4,2,19,12},5/1),({4,2,18,13},6/1),({4,2,17,14},5/1),({4,2,16,15},3/1),({3,3,23,8},1/1),({3,3,22,9},1/1),({3,3,21,10},3/1),({3,3,20,11},4/1),({3,3,19,12},5/1),({3,3,18,13},5/1),({3,3,17,14},5/1),({3,3,16,15},2/1)}, (6,0) => {}};
--dw stands for dominant weights
dw0115 = new HashTable from {(7,0) => {}, (6,1) => {({5,2,23,13},1/1),({4,3,25,11},1/1)}, (5,2) => {}, (6,2) => {}, (7,1) => {({5,3,26,15},1/1),({4,4,27,14},1/1)}, (8,0) => {}, (7,2) => {}, (9,0) => {}, (8,1) => {({5,4,28,18},1/1)}, (8,2) => {}, (9,1) => {({5,5,29,22},1/1)}, (9,2) => {}, (10,1) => {}, (11,1) => {}, (0,0) => {({0,0,1,0},1/1)}, (1,0) => {({1,0,5,1},1/1)}, (2,0) => {}, (1,1) => {}, (3,0) => {}, (2,1) => {({3,0,11,5},1/1),({2,1,13,3},1/1)}, (2,2) => {}, (3,1) => {({4,0,13,8},1/1),({3,1,16,5},1/1),({2,2,17,4},1/1)}, (4,0) => {}, (5,0) => {}, (4,1) => {({5,0,14,12},1/1),({4,1,18,8},1/1),({3,2,20,6},1/1)}, (3,2) => {}, (6,0) => {}, (5,1) => {({5,1,19,12},1/1),({4,2,22,9},1/1),({3,3,23,8},1/1)}, (4,2) => {}};
--dmw stands for dominant monomial weights
dmw0115 = new HashTable from {(7,0) => {}, (6,1) => {{6,1,18,18},{4,3,22,14},{5,2,21,15}}, (5,2) => {}, (6,2) => {}, (7,1) => {{6,2,21,20},{5,3,23,18}}, (8,0) => {}, (7,2) => {}, (9,0) => {}, (8,1) => {{6,3,23,23},{5,4,24,22}}, (8,2) => {}, (9,1) => {{6,4,24,27}}, (9,2) => {}, (10,1) => {}, (11,1) => {}, (0,0) => {}, (1,0) => {}, (2,0) => {}, (1,1) => {}, (3,0) => {}, (2,1) => {{3,0,11,5}}, (2,2) => {}, (3,1) => {{4,0,13,8},{3,1,15,6}}, (4,0) => {}, (5,0) => {}, (4,1) => {{5,0,14,12},{3,2,18,8},{4,1,17,9}}, (3,2) => {}, (6,0) => {}, (5,1) => {{5,1,18,13},{6,0,14,17},{4,2,20,11}}, (4,2) => {}};
end;