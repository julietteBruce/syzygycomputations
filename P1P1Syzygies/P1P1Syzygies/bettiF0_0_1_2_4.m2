A := QQ[t_0,t_1,t_2,t_3];
--tb stands for Total Betti numbers
tb0124 = new HashTable from {(5,2) => 0, (7,0) => 0, (6,1) => 7920/1, (7,1) => 6666/1, (8,0) => 0, (6,2) => 0, (7,2) => 0, (9,0) => 0, (8,1) => 3916/1, (8,2) => 0, (10,0) => 0, (9,1) => 1518/1, (10,1) => 328/1, (11,0) => 0, (9,2) => 0, (11,1) => 10/1, (12,0) => 0, (10,2) => 0, (12,1) => 0, (11,2) => 12/1, (12,2) => 2/1, (13,2) => 0, (0,0) => 2/1, (1,0) => 12/1, (2,0) => 0, (1,1) => 10/1, (3,0) => 0, (2,1) => 328/1, (2,2) => 0, (4,0) => 0, (3,1) => 1518/1, (3,2) => 0, (5,0) => 0, (4,1) => 3916/1, (4,2) => 0, (6,0) => 0, (5,1) => 6666/1};
--sb represents the betti numbers as sums of Schur functors
sb0124 = new HashTable from {(5,2) => {}, (7,0) => {}, (6,1) => {({11,3,18,11},2/1),({11,3,17,12},2/1),({11,3,16,13},2/1),({11,3,15,14},2/1),({10,4,20,9},2/1),({10,4,19,10},4/1),({10,4,18,11},8/1),({10,4,17,12},10/1),({10,4,16,13},10/1),({10,4,15,14},6/1),({9,5,21,8},2/1),({9,5,20,9},6/1),({9,5,19,10},12/1),({9,5,18,11},18/1),({9,5,17,12},22/1),({9,5,16,13},20/1),({9,5,15,14},12/1),({8,6,22,7},2/1),({8,6,21,8},4/1),({8,6,20,9},10/1),({8,6,19,10},16/1),({8,6,18,11},24/1),({8,6,17,12},26/1),({8,6,16,13},24/1),({8,6,15,14},14/1),({7,7,21,8},2/1),({7,7,20,9},4/1),({7,7,19,10},8/1),({7,7,18,11},10/1),({7,7,17,12},12/1),({7,7,16,13},10/1),({7,7,15,14},6/1)}, (7,1) => {({12,4,20,13},1/1),({12,4,19,14},2/1),({12,4,18,15},2/1),({12,4,17,16},1/1),({11,5,22,11},1/1),({11,5,21,12},4/1),({11,5,20,13},6/1),({11,5,19,14},8/1),({11,5,18,15},9/1),({11,5,17,16},5/1),({10,6,23,10},2/1),({10,6,22,11},5/1),({10,6,21,12},10/1),({10,6,20,13},16/1),({10,6,19,14},19/1),({10,6,18,15},17/1),({10,6,17,16},11/1),({9,7,24,9},1/1),({9,7,23,10},3/1),({9,7,22,11},8/1),({9,7,21,12},14/1),({9,7,20,13},20/1),({9,7,19,14},23/1),({9,7,18,15},20/1),({9,7,17,16},12/1),({8,8,24,9},1/1),({8,8,23,10},2/1),({8,8,22,11},4/1),({8,8,21,12},7/1),({8,8,20,13},9/1),({8,8,19,14},11/1),({8,8,18,15},10/1),({8,8,17,16},5/1)}, (8,0) => {}, (6,2) => {}, (7,2) => {}, (9,0) => {}, (8,1) => {({13,5,21,16},1/1),({13,5,20,17},1/1),({13,5,19,18},1/1),({12,6,23,14},2/1),({12,6,22,15},3/1),({12,6,21,16},5/1),({12,6,20,17},5/1),({12,6,19,18},3/1),({11,7,25,12},1/1),({11,7,24,13},3/1),({11,7,23,14},6/1),({11,7,22,15},10/1),({11,7,21,16},12/1),({11,7,20,17},11/1),({11,7,19,18},7/1),({10,8,25,12},2/1),({10,8,24,13},4/1),({10,8,23,14},9/1),({10,8,22,15},12/1),({10,8,21,16},15/1),({10,8,20,17},13/1),({10,8,19,18},8/1),({9,9,26,11},1/1),({9,9,25,12},1/1),({9,9,24,13},3/1),({9,9,23,14},5/1),({9,9,22,15},7/1),({9,9,21,16},7/1),({9,9,20,17},7/1),({9,9,19,18},3/1)}, (8,2) => {}, (10,0) => {}, (9,1) => {({14,6,21,20},1/1),({13,7,24,17},1/1),({13,7,23,18},2/1),({13,7,22,19},2/1),({13,7,21,20},1/1),({12,8,26,15},1/1),({12,8,25,16},2/1),({12,8,24,17},4/1),({12,8,23,18},5/1),({12,8,22,19},5/1),({12,8,21,20},3/1),({11,9,27,14},1/1),({11,9,26,15},2/1),({11,9,25,16},4/1),({11,9,24,17},6/1),({11,9,23,18},7/1),({11,9,22,19},6/1),({11,9,21,20},4/1),({10,10,26,15},1/1),({10,10,25,16},2/1),({10,10,24,17},3/1),({10,10,23,18},3/1),({10,10,22,19},3/1),({10,10,21,20},2/1)}, (10,1) => {({14,8,24,21},1/1),({13,9,26,19},1/1),({13,9,25,20},1/1),({13,9,24,21},1/1),({13,9,23,22},1/1),({12,10,28,17},1/1),({12,10,27,18},1/1),({12,10,26,19},2/1),({12,10,25,20},2/1),({12,10,24,21},2/1),({12,10,23,22},1/1),({11,11,27,18},1/1),({11,11,25,20},1/1),({11,11,24,21},1/1)}, (11,0) => {}, (9,2) => {}, (11,1) => {({12,12,29,20},1/1)}, (12,0) => {}, (10,2) => {}, (12,1) => {}, (11,2) => {({14,12,28,25},1/1)}, (12,2) => {({14,14,29,28},1/1)}, (13,2) => {}, (0,0) => {({0,0,1,0},1/1)}, (1,0) => {({2,0,4,1},1/1)}, (2,0) => {}, (1,1) => {({2,2,9,0},1/1)}, (3,0) => {}, (2,1) => {({6,0,8,5},1/1),({5,1,10,3},1/1),({5,1,9,4},1/1),({5,1,8,5},1/1),({5,1,7,6},1/1),({4,2,12,1},1/1),({4,2,11,2},1/1),({4,2,10,3},2/1),({4,2,9,4},2/1),({4,2,8,5},2/1),({4,2,7,6},1/1),({3,3,11,2},1/1),({3,3,9,4},1/1),({3,3,8,5},1/1)}, (2,2) => {}, (4,0) => {}, (3,1) => {({8,0,9,8},1/1),({7,1,12,5},1/1),({7,1,11,6},2/1),({7,1,10,7},2/1),({7,1,9,8},1/1),({6,2,14,3},1/1),({6,2,13,4},2/1),({6,2,12,5},4/1),({6,2,11,6},5/1),({6,2,10,7},5/1),({6,2,9,8},3/1),({5,3,15,2},1/1),({5,3,14,3},2/1),({5,3,13,4},4/1),({5,3,12,5},6/1),({5,3,11,6},7/1),({5,3,10,7},6/1),({5,3,9,8},4/1),({4,4,14,3},1/1),({4,4,13,4},2/1),({4,4,12,5},3/1),({4,4,11,6},3/1),({4,4,10,7},3/1),({4,4,9,8},2/1)}, (3,2) => {}, (5,0) => {}, (4,1) => {({9,1,13,8},1/1),({9,1,12,9},1/1),({9,1,11,10},1/1),({8,2,15,6},2/1),({8,2,14,7},3/1),({8,2,13,8},5/1),({8,2,12,9},5/1),({8,2,11,10},3/1),({7,3,17,4},1/1),({7,3,16,5},3/1),({7,3,15,6},6/1),({7,3,14,7},10/1),({7,3,13,8},12/1),({7,3,12,9},11/1),({7,3,11,10},7/1),({6,4,17,4},2/1),({6,4,16,5},4/1),({6,4,15,6},9/1),({6,4,14,7},12/1),({6,4,13,8},15/1),({6,4,12,9},13/1),({6,4,11,10},8/1),({5,5,18,3},1/1),({5,5,17,4},1/1),({5,5,16,5},3/1),({5,5,15,6},5/1),({5,5,14,7},7/1),({5,5,13,8},7/1),({5,5,12,9},7/1),({5,5,11,10},3/1)}, (4,2) => {}, (6,0) => {}, (5,1) => {({10,2,16,9},1/1),({10,2,15,10},2/1),({10,2,14,11},2/1),({10,2,13,12},1/1),({9,3,18,7},1/1),({9,3,17,8},4/1),({9,3,16,9},6/1),({9,3,15,10},8/1),({9,3,14,11},9/1),({9,3,13,12},5/1),({8,4,19,6},2/1),({8,4,18,7},5/1),({8,4,17,8},10/1),({8,4,16,9},16/1),({8,4,15,10},19/1),({8,4,14,11},17/1),({8,4,13,12},11/1),({7,5,20,5},1/1),({7,5,19,6},3/1),({7,5,18,7},8/1),({7,5,17,8},14/1),({7,5,16,9},20/1),({7,5,15,10},23/1),({7,5,14,11},20/1),({7,5,13,12},12/1),({6,6,20,5},1/1),({6,6,19,6},2/1),({6,6,18,7},4/1),({6,6,17,8},7/1),({6,6,16,9},9/1),({6,6,15,10},11/1),({6,6,14,11},10/1),({6,6,13,12},5/1)}};
end;