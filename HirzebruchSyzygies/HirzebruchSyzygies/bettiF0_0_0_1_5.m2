A := QQ[t_0,t_1,t_2,t_3];
--tb stands for Total Betti numbers
tb0015 = new HashTable from {(5,2) => 0, (6,1) => 720/1, (7,0) => 0, (6,2) => 0, (7,1) => 315/1, (8,0) => 0, (8,1) => 80/1, (9,0) => 0, (7,2) => 0, (9,1) => 9/1, (8,2) => 0, (9,2) => 0, (10,1) => 0, (11,1) => 0, (0,0) => 1/1, (1,0) => 0, (2,0) => 0, (1,1) => 45/1, (3,0) => 0, (2,1) => 240/1, (2,2) => 0, (3,1) => 630/1, (4,0) => 0, (3,2) => 0, (4,1) => 1008/1, (5,0) => 0, (4,2) => 0, (5,1) => 1050/1, (6,0) => 0};
--sb represents the betti numbers as sums of Schur functors
sb0015 = new HashTable from {(5,2) => {}, (6,1) => {({5,2,23,12},1/1),({5,2,22,13},1/1),({5,2,21,14},2/1),({5,2,20,15},2/1),({5,2,19,16},2/1),({5,2,18,17},1/1),({4,3,25,10},1/1),({4,3,24,11},2/1),({4,3,23,12},3/1),({4,3,22,13},5/1),({4,3,21,14},6/1),({4,3,20,15},6/1),({4,3,19,16},5/1),({4,3,18,17},3/1)}, (7,0) => {}, (6,2) => {}, (7,1) => {({5,3,26,14},1/1),({5,3,25,15},1/1),({5,3,24,16},2/1),({5,3,23,17},2/1),({5,3,22,18},2/1),({5,3,21,19},1/1),({5,3,20,20},1/1),({4,4,27,13},1/1),({4,4,26,14},1/1),({4,4,25,15},2/1),({4,4,24,16},2/1),({4,4,23,17},3/1),({4,4,22,18},2/1),({4,4,21,19},2/1)}, (8,0) => {}, (8,1) => {({5,4,28,17},1/1),({5,4,27,18},1/1),({5,4,26,19},1/1),({5,4,25,20},1/1),({5,4,24,21},1/1)}, (9,0) => {}, (7,2) => {}, (9,1) => {({5,5,29,21},1/1)}, (8,2) => {}, (9,2) => {}, (10,1) => {}, (11,1) => {}, (0,0) => {({0,0,0,0},1/1)}, (1,0) => {}, (2,0) => {}, (1,1) => {({2,0,8,2},1/1),({2,0,6,4},1/1),({1,1,9,1},1/1),({1,1,7,3},1/1),({1,1,5,5},1/1)}, (3,0) => {}, (2,1) => {({3,0,11,4},1/1),({3,0,10,5},1/1),({3,0,9,6},1/1),({3,0,8,7},1/1),({2,1,13,2},1/1),({2,1,12,3},2/1),({2,1,11,4},2/1),({2,1,10,5},3/1),({2,1,9,6},3/1),({2,1,8,7},1/1)}, (2,2) => {}, (3,1) => {({4,0,13,7},1/1),({4,0,12,8},1/1),({4,0,11,9},1/1),({3,1,16,4},1/1),({3,1,15,5},2/1),({3,1,14,6},4/1),({3,1,13,7},4/1),({3,1,12,8},5/1),({3,1,11,9},3/1),({3,1,10,10},2/1),({2,2,17,3},1/1),({2,2,16,4},1/1),({2,2,15,5},3/1),({2,2,14,6},3/1),({2,2,13,7},5/1),({2,2,12,8},3/1),({2,2,11,9},4/1)}, (4,0) => {}, (3,2) => {}, (4,1) => {({5,0,14,11},1/1),({4,1,18,7},1/1),({4,1,17,8},2/1),({4,1,16,9},3/1),({4,1,15,10},4/1),({4,1,14,11},3/1),({4,1,13,12},2/1),({3,2,20,5},1/1),({3,2,19,6},2/1),({3,2,18,7},4/1),({3,2,17,8},6/1),({3,2,16,9},8/1),({3,2,15,10},8/1),({3,2,14,11},7/1),({3,2,13,12},4/1)}, (5,0) => {}, (4,2) => {}, (5,1) => {({5,1,19,11},1/1),({5,1,18,12},1/1),({5,1,17,13},1/1),({5,1,16,14},1/1),({5,1,15,15},1/1),({4,2,22,8},1/1),({4,2,21,9},2/1),({4,2,20,10},4/1),({4,2,19,11},5/1),({4,2,18,12},7/1),({4,2,17,13},6/1),({4,2,16,14},5/1),({4,2,15,15},1/1),({3,3,23,7},1/1),({3,3,22,8},1/1),({3,3,21,9},3/1),({3,3,20,10},4/1),({3,3,19,11},6/1),({3,3,18,12},5/1),({3,3,17,13},7/1),({3,3,16,14},3/1),({3,3,15,15},2/1)}, (6,0) => {}};
--dw stands for dominant weights
dw0015 = new HashTable from {(7,0) => {}, (6,1) => {({5,2,23,12},1/1),({4,3,25,10},1/1)}, (5,2) => {}, (6,2) => {}, (7,1) => {({5,3,26,14},1/1),({4,4,27,13},1/1)}, (8,0) => {}, (7,2) => {}, (9,0) => {}, (8,1) => {({5,4,28,17},1/1)}, (8,2) => {}, (9,1) => {({5,5,29,21},1/1)}, (9,2) => {}, (10,1) => {}, (11,1) => {}, (0,0) => {({0,0,0,0},1/1)}, (1,0) => {}, (2,0) => {}, (1,1) => {({2,0,8,2},1/1),({1,1,9,1},1/1)}, (3,0) => {}, (2,1) => {({3,0,11,4},1/1),({2,1,13,2},1/1)}, (2,2) => {}, (3,1) => {({4,0,13,7},1/1),({3,1,16,4},1/1),({2,2,17,3},1/1)}, (4,0) => {}, (5,0) => {}, (4,1) => {({5,0,14,11},1/1),({4,1,18,7},1/1),({3,2,20,5},1/1)}, (3,2) => {}, (6,0) => {}, (5,1) => {({5,1,19,11},1/1),({4,2,22,8},1/1),({3,3,23,7},1/1)}, (4,2) => {}};
--dmw stands for dominant monomial weights
dmw0015 = new HashTable from {(7,0) => {}, (6,1) => {{4,3,25,10},{5,2,23,12}}, (5,2) => {}, (6,2) => {}, (7,1) => {{5,3,26,14},{4,4,27,13}}, (8,0) => {}, (7,2) => {}, (9,0) => {}, (8,1) => {{5,4,28,17}}, (8,2) => {{5,5,24,26}}, (9,1) => {{5,5,29,21}}, (9,2) => {{5,6,29,26}}, (10,1) => {}, (11,1) => {}, (0,0) => {}, (1,0) => {}, (2,0) => {}, (1,1) => {{2,0,8,2}}, (3,0) => {}, (2,1) => {{2,1,13,2},{3,0,11,4}}, (2,2) => {}, (3,1) => {{3,1,16,4},{4,0,13,7},{2,2,17,3}}, (4,0) => {}, (5,0) => {}, (4,1) => {{5,0,14,11},{4,1,18,7},{3,2,20,5}}, (3,2) => {}, (6,0) => {}, (5,1) => {{4,2,22,8},{5,1,19,11},{3,3,23,7}}, (4,2) => {}};
end;