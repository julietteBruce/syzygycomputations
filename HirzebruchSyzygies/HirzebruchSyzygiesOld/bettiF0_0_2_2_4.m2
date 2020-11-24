A := QQ[t_0,t_1,t_2,t_3];
--tb stands for Total Betti numbers
tb0224 = new HashTable from {(6,1) => 7920/1, (7,0) => 0, (5,2) => 0, (7,1) => 7095/1, (8,0) => 0, (6,2) => 0, (8,1) => 4488/1, (9,0) => 0, (7,2) => 0, (9,1) => 1947/1, (10,0) => 0, (8,2) => 0, (10,1) => 536/1, (11,0) => 0, (9,2) => 0, (11,1) => 75/1, (12,0) => 0, (10,2) => 0, (12,1) => 0, (11,2) => 0, (12,2) => 1/1, (13,2) => 0, (0,0) => 3/1, (1,0) => 24/1, (1,1) => 11/1, (2,0) => 66/1, (2,1) => 120/1, (3,0) => 0, (3,1) => 1089/1, (4,0) => 0, (2,2) => 0, (4,1) => 3344/1, (5,0) => 0, (3,2) => 0, (5,1) => 6237/1, (6,0) => 0, (4,2) => 0};
--sb represents the betti numbers as sums of Schur functors
sb0224 = new HashTable from {(6,1) => {({11,3,19,11},1/1),({11,3,18,12},2/1),({11,3,17,13},2/1),({11,3,16,14},2/1),({11,3,15,15},1/1),({10,4,21,9},1/1),({10,4,20,10},3/1),({10,4,19,11},6/1),({10,4,18,12},9/1),({10,4,17,13},10/1),({10,4,16,14},8/1),({10,4,15,15},3/1),({9,5,22,8},1/1),({9,5,21,9},4/1),({9,5,20,10},9/1),({9,5,19,11},15/1),({9,5,18,12},20/1),({9,5,17,13},21/1),({9,5,16,14},16/1),({9,5,15,15},6/1),({8,6,23,7},1/1),({8,6,22,8},3/1),({8,6,21,9},7/1),({8,6,20,10},13/1),({8,6,19,11},20/1),({8,6,18,12},25/1),({8,6,17,13},25/1),({8,6,16,14},19/1),({8,6,15,15},7/1),({7,7,22,8},1/1),({7,7,21,9},3/1),({7,7,20,10},6/1),({7,7,19,11},9/1),({7,7,18,12},11/1),({7,7,17,13},11/1),({7,7,16,14},8/1),({7,7,15,15},3/1)}, (7,0) => {}, (5,2) => {}, (7,1) => {({12,4,21,13},1/1),({12,4,20,14},1/1),({12,4,19,15},3/1),({12,4,18,16},1/1),({12,4,17,17},1/1),({11,5,22,12},3/1),({11,5,21,13},5/1),({11,5,20,14},8/1),({11,5,19,15},9/1),({11,5,18,16},8/1),({11,5,17,17},2/1),({10,6,24,10},1/1),({10,6,23,11},4/1),({10,6,22,12},7/1),({10,6,21,13},15/1),({10,6,20,14},18/1),({10,6,19,15},21/1),({10,6,18,16},14/1),({10,6,17,17},8/1),({9,7,24,10},2/1),({9,7,23,11},5/1),({9,7,22,12},12/1),({9,7,21,13},17/1),({9,7,20,14},24/1),({9,7,19,15},22/1),({9,7,18,16},19/1),({9,7,17,17},5/1),({8,8,25,9},1/1),({8,8,24,10},1/1),({8,8,23,11},4/1),({8,8,22,12},5/1),({8,8,21,13},10/1),({8,8,20,14},9/1),({8,8,19,15},14/1),({8,8,18,16},6/1),({8,8,17,17},5/1)}, (8,0) => {}, (6,2) => {}, (8,1) => {({13,5,22,16},1/1),({13,5,21,17},1/1),({13,5,20,18},1/1),({13,5,19,19},1/1),({12,6,24,14},1/1),({12,6,23,15},3/1),({12,6,22,16},5/1),({12,6,21,17},6/1),({12,6,20,18},5/1),({12,6,19,19},2/1),({11,7,25,13},2/1),({11,7,24,14},5/1),({11,7,23,15},9/1),({11,7,22,16},13/1),({11,7,21,17},14/1),({11,7,20,18},11/1),({11,7,19,19},4/1),({10,8,26,12},1/1),({10,8,25,13},3/1),({10,8,24,14},7/1),({10,8,23,15},12/1),({10,8,22,16},16/1),({10,8,21,17},17/1),({10,8,20,18},13/1),({10,8,19,19},5/1),({9,9,26,12},1/1),({9,9,25,13},2/1),({9,9,24,14},4/1),({9,9,23,15},7/1),({9,9,22,16},8/1),({9,9,21,17},8/1),({9,9,20,18},6/1),({9,9,19,19},2/1)}, (9,0) => {}, (7,2) => {}, (9,1) => {({14,6,22,20},1/1),({13,7,25,17},1/1),({13,7,24,18},2/1),({13,7,23,19},3/1),({13,7,22,20},2/1),({13,7,21,21},1/1),({12,8,26,16},2/1),({12,8,25,17},3/1),({12,8,24,18},7/1),({12,8,23,19},6/1),({12,8,22,20},7/1),({12,8,21,21},1/1),({11,9,27,15},2/1),({11,9,26,16},3/1),({11,9,25,17},7/1),({11,9,24,18},8/1),({11,9,23,19},10/1),({11,9,22,20},6/1),({11,9,21,21},4/1),({10,10,26,16},2/1),({10,10,25,17},2/1),({10,10,24,18},5/1),({10,10,23,19},3/1),({10,10,22,20},5/1)}, (10,0) => {}, (8,2) => {}, (10,1) => {({14,8,25,21},1/1),({14,8,24,22},1/1),({13,9,27,19},1/1),({13,9,26,20},2/1),({13,9,25,21},2/1),({13,9,24,22},2/1),({13,9,23,23},1/1),({12,10,28,18},1/1),({12,10,27,19},2/1),({12,10,26,20},3/1),({12,10,25,21},4/1),({12,10,24,22},3/1),({12,10,23,23},1/1),({11,11,28,18},1/1),({11,11,27,19},1/1),({11,11,26,20},1/1),({11,11,25,21},2/1),({11,11,24,22},1/1)}, (11,0) => {}, (9,2) => {}, (11,1) => {({14,10,27,23},1/1),({14,10,25,25},1/1),({13,11,28,22},1/1),({13,11,26,24},1/1),({12,12,29,21},1/1),({12,12,27,23},1/1),({12,12,25,25},1/1)}, (12,0) => {}, (10,2) => {}, (12,1) => {}, (11,2) => {}, (12,2) => {({14,14,29,29},1/1)}, (13,2) => {}, (0,0) => {({0,0,2,0},1/1)}, (1,0) => {({2,0,5,1},1/1),({2,0,4,2},1/1)}, (1,1) => {({2,2,10,0},1/1)}, (2,0) => {({4,0,7,3},1/1),({4,0,5,5},1/1),({3,1,8,2},1/1),({3,1,6,4},1/1),({2,2,7,3},1/1),({2,2,5,5},1/1)}, (2,1) => {({4,2,13,1},1/1),({4,2,12,2},1/1),({4,2,11,3},1/1),({4,2,10,4},1/1)}, (3,0) => {}, (3,1) => {({8,0,9,9},1/1),({7,1,12,6},1/1),({7,1,11,7},1/1),({7,1,10,8},1/1),({6,2,15,3},1/1),({6,2,14,4},1/1),({6,2,13,5},3/1),({6,2,12,6},2/1),({6,2,11,7},4/1),({6,2,10,8},1/1),({6,2,9,9},2/1),({5,3,16,2},1/1),({5,3,15,3},1/1),({5,3,14,4},3/1),({5,3,13,5},3/1),({5,3,12,6},5/1),({5,3,11,7},3/1),({5,3,10,8},4/1),({4,4,15,3},1/1),({4,4,14,4},1/1),({4,4,13,5},3/1),({4,4,12,6},1/1),({4,4,11,7},3/1),({4,4,9,9},2/1)}, (4,0) => {}, (2,2) => {}, (4,1) => {({9,1,13,9},1/1),({9,1,12,10},1/1),({8,2,16,6},1/1),({8,2,15,7},2/1),({8,2,14,8},3/1),({8,2,13,9},4/1),({8,2,12,10},3/1),({8,2,11,11},1/1),({7,3,18,4},1/1),({7,3,17,5},2/1),({7,3,16,6},4/1),({7,3,15,7},7/1),({7,3,14,8},9/1),({7,3,13,9},9/1),({7,3,12,10},7/1),({7,3,11,11},3/1),({6,4,18,4},1/1),({6,4,17,5},3/1),({6,4,16,6},6/1),({6,4,15,7},9/1),({6,4,14,8},11/1),({6,4,13,9},11/1),({6,4,12,10},8/1),({6,4,11,11},3/1),({5,5,19,3},1/1),({5,5,18,4},1/1),({5,5,17,5},2/1),({5,5,16,6},4/1),({5,5,15,7},5/1),({5,5,14,8},6/1),({5,5,13,9},6/1),({5,5,12,10},4/1),({5,5,11,11},1/1)}, (5,0) => {}, (3,2) => {}, (5,1) => {({10,2,16,10},2/1),({10,2,15,11},1/1),({10,2,14,12},2/1),({9,3,19,7},1/1),({9,3,18,8},2/1),({9,3,17,9},5/1),({9,3,16,10},6/1),({9,3,15,11},8/1),({9,3,14,12},6/1),({9,3,13,13},3/1),({8,4,20,6},1/1),({8,4,19,7},3/1),({8,4,18,8},8/1),({8,4,17,9},11/1),({8,4,16,10},17/1),({8,4,15,11},15/1),({8,4,14,12},14/1),({8,4,13,13},3/1),({7,5,21,5},1/1),({7,5,20,6},2/1),({7,5,19,7},6/1),({7,5,18,8},10/1),({7,5,17,9},17/1),({7,5,16,10},19/1),({7,5,15,11},21/1),({7,5,14,12},13/1),({7,5,13,13},7/1),({6,6,20,6},2/1),({6,6,19,7},2/1),({6,6,18,8},6/1),({6,6,17,9},6/1),({6,6,16,10},11/1),({6,6,15,11},7/1),({6,6,14,12},9/1)}, (6,0) => {}, (4,2) => {}};
--dw stands for dominant weights
dw0224 = new HashTable from {(5,2) => {}, (7,0) => {}, (6,1) => {({11,3,19,11},1/1)}, (7,1) => {({12,4,21,13},1/1)}, (8,0) => {}, (6,2) => {}, (7,2) => {}, (9,0) => {}, (8,1) => {({13,5,22,16},1/1)}, (8,2) => {}, (10,0) => {}, (9,1) => {({14,6,22,20},1/1)}, (10,1) => {({14,8,25,21},1/1)}, (11,0) => {}, (9,2) => {}, (11,1) => {({14,10,27,23},1/1)}, (12,0) => {}, (10,2) => {}, (12,1) => {}, (11,2) => {}, (12,2) => {({14,14,29,29},1/1)}, (13,2) => {}, (0,0) => {({0,0,2,0},1/1)}, (1,0) => {({2,0,5,1},1/1)}, (2,0) => {({4,0,7,3},1/1)}, (1,1) => {({2,2,10,0},1/1)}, (3,0) => {}, (2,1) => {({4,2,13,1},1/1)}, (2,2) => {}, (4,0) => {}, (3,1) => {({8,0,9,9},1/1)}, (3,2) => {}, (5,0) => {}, (4,1) => {({9,1,13,9},1/1)}, (4,2) => {}, (6,0) => {}, (5,1) => {({10,2,16,10},2/1)}};
end;