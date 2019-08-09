A := QQ[t_0,t_1,t_2,t_3];
--tb stands for Total Betti numbers
tb1023 = new HashTable from {(5,2) => 0, (6,1) => 768/1, (7,0) => 0, (6,2) => 0, (7,1) => 342/1, (8,0) => 0, (8,1) => 88/1, (9,0) => 0, (7,2) => 0, (9,1) => 10/1, (8,2) => 0, (9,2) => 0, (10,1) => 0, (11,1) => 0, (0,0) => 2/1, (1,0) => 8/1, (2,0) => 0, (1,1) => 18/1, (3,0) => 0, (2,1) => 192/1, (2,2) => 0, (3,1) => 588/1, (4,0) => 0, (3,2) => 0, (4,1) => 1008/1, (5,0) => 0, (4,2) => 0, (5,1) => 1092/1, (6,0) => 0};
--sb represents the betti numbers as sums of Schur functors
sb1023 = new HashTable from {(5,2) => {}, (6,1) => {({11,4,11,10},1/1),({10,5,14,7},1/1),({10,5,13,8},2/1),({10,5,12,9},3/1),({10,5,11,10},2/1),({9,6,15,6},1/1),({9,6,14,7},2/1),({9,6,13,8},5/1),({9,6,12,9},6/1),({9,6,11,10},4/1),({8,7,15,6},1/1),({8,7,14,7},3/1),({8,7,13,8},5/1),({8,7,12,9},5/1),({8,7,11,10},4/1)}, (7,0) => {}, (6,2) => {}, (7,1) => {({11,6,14,10},1/1),({11,6,13,11},1/1),({11,6,12,12},1/1),({10,7,16,8},1/1),({10,7,15,9},2/1),({10,7,14,10},3/1),({10,7,13,11},3/1),({10,7,12,12},1/1),({9,8,16,8},1/1),({9,8,15,9},2/1),({9,8,14,10},3/1),({9,8,13,11},3/1),({9,8,12,12},1/1)}, (8,0) => {}, (8,1) => {({11,8,16,11},1/1),({11,8,15,12},1/1),({11,8,14,13},1/1),({10,9,17,10},1/1),({10,9,16,11},1/1),({10,9,15,12},1/1),({10,9,14,13},1/1)}, (9,0) => {}, (7,2) => {}, (9,1) => {({11,10,17,13},1/1)}, (8,2) => {}, (9,2) => {}, (10,1) => {}, (11,1) => {}, (0,0) => {({1,0,0,0},1/1)}, (1,0) => {({2,1,3,0},1/1)}, (2,0) => {}, (1,1) => {({5,0,4,2},1/1)}, (3,0) => {}, (2,1) => {({7,0,5,4},1/1),({6,1,7,2},1/1),({6,1,6,3},1/1),({6,1,5,4},1/1),({5,2,7,2},1/1),({5,2,6,3},2/1),({5,2,5,4},1/1),({4,3,8,1},1/1),({4,3,7,2},1/1),({4,3,6,3},1/1),({4,3,5,4},1/1)}, (2,2) => {}, (3,1) => {({8,1,8,4},1/1),({8,1,7,5},1/1),({7,2,9,3},1/1),({7,2,8,4},2/1),({7,2,7,5},3/1),({7,2,6,6},1/1),({6,3,10,2},1/1),({6,3,9,3},3/1),({6,3,8,4},4/1),({6,3,7,5},4/1),({6,3,6,6},2/1),({5,4,10,2},1/1),({5,4,9,3},2/1),({5,4,8,4},4/1),({5,4,7,5},3/1),({5,4,6,6},1/1)}, (4,0) => {}, (3,2) => {}, (4,1) => {({9,2,10,5},1/1),({9,2,9,6},1/1),({9,2,8,7},1/1),({8,3,11,4},1/1),({8,3,10,5},3/1),({8,3,9,6},4/1),({8,3,8,7},3/1),({7,4,12,3},1/1),({7,4,11,4},3/1),({7,4,10,5},6/1),({7,4,9,6},7/1),({7,4,8,7},5/1),({6,5,12,3},1/1),({6,5,11,4},3/1),({6,5,10,5},5/1),({6,5,9,6},6/1),({6,5,8,7},4/1)}, (5,0) => {}, (4,2) => {}, (5,1) => {({10,3,11,7},1/1),({10,3,10,8},1/1),({9,4,13,5},1/1),({9,4,12,6},2/1),({9,4,11,7},4/1),({9,4,10,8},4/1),({9,4,9,9},2/1),({8,5,13,5},2/1),({8,5,12,6},5/1),({8,5,11,7},7/1),({8,5,10,8},7/1),({8,5,9,9},3/1),({7,6,14,4},1/1),({7,6,13,5},2/1),({7,6,12,6},5/1),({7,6,11,7},7/1),({7,6,10,8},6/1),({7,6,9,9},2/1)}, (6,0) => {}};
end;