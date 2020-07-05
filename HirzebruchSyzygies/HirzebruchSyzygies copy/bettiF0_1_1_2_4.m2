A := QQ[t_0,t_1,t_2,t_3];
--tb stands for Total Betti numbers
tb1124 = new HashTable from {(6,1) => 7920/1, (7,0) => 0, (5,2) => 0, (7,1) => 7524/1, (8,0) => 0, (6,2) => 0, (8,1) => 5060/1, (9,0) => 0, (7,2) => 0, (9,1) => 2376/1, (10,0) => 0, (8,2) => 0, (10,1) => 744/1, (11,0) => 0, (9,2) => 0, (11,1) => 140/1, (12,0) => 0, (10,2) => 0, (12,1) => 12/1, (11,2) => 0, (12,2) => 0, (13,1) => 0, (14,1) => 0, (0,0) => 4/1, (1,0) => 36/1, (1,1) => 0, (2,0) => 120/1, (2,1) => 32/1, (3,0) => 120/1, (3,1) => 660/1, (4,0) => 0, (2,2) => 0, (4,1) => 2772/1, (5,0) => 0, (3,2) => 0, (5,1) => 5808/1, (6,0) => 0, (4,2) => 0};
--sb represents the betti numbers as sums of Schur functors
sb1124 = new HashTable from {(6,1) => {({12,3,18,11},1/1),({12,3,17,12},1/1),({12,3,16,13},1/1),({12,3,15,14},1/1),({11,4,20,9},1/1),({11,4,19,10},2/1),({11,4,18,11},5/1),({11,4,17,12},6/1),({11,4,16,13},6/1),({11,4,15,14},4/1),({10,5,21,8},1/1),({10,5,20,9},4/1),({10,5,19,10},8/1),({10,5,18,11},13/1),({10,5,17,12},16/1),({10,5,16,13},15/1),({10,5,15,14},9/1),({9,6,22,7},1/1),({9,6,21,8},3/1),({9,6,20,9},8/1),({9,6,19,10},14/1),({9,6,18,11},21/1),({9,6,17,12},24/1),({9,6,16,13},22/1),({9,6,15,14},13/1),({8,7,22,7},1/1),({8,7,21,8},3/1),({8,7,20,9},7/1),({8,7,19,10},12/1),({8,7,18,11},17/1),({8,7,17,12},19/1),({8,7,16,13},17/1),({8,7,15,14},10/1)}, (7,0) => {}, (5,2) => {}, (7,1) => {({13,4,19,14},1/1),({13,4,18,15},1/1),({12,5,22,11},1/1),({12,5,21,12},2/1),({12,5,20,13},4/1),({12,5,19,14},5/1),({12,5,18,15},6/1),({12,5,17,16},3/1),({11,6,23,10},1/1),({11,6,22,11},3/1),({11,6,21,12},8/1),({11,6,20,13},12/1),({11,6,19,14},15/1),({11,6,18,15},14/1),({11,6,17,16},9/1),({10,7,24,9},1/1),({10,7,23,10},3/1),({10,7,22,11},8/1),({10,7,21,12},14/1),({10,7,20,13},21/1),({10,7,19,14},24/1),({10,7,18,15},22/1),({10,7,17,16},13/1),({9,8,24,9},1/1),({9,8,23,10},3/1),({9,8,22,11},7/1),({9,8,21,12},12/1),({9,8,20,13},17/1),({9,8,19,14},20/1),({9,8,18,15},17/1),({9,8,17,16},10/1)}, (8,0) => {}, (6,2) => {}, (8,1) => {({14,5,19,18},1/1),({13,6,23,14},1/1),({13,6,22,15},2/1),({13,6,21,16},3/1),({13,6,20,17},3/1),({13,6,19,18},2/1),({12,7,25,12},1/1),({12,7,24,13},2/1),({12,7,23,14},5/1),({12,7,22,15},8/1),({12,7,21,16},11/1),({12,7,20,17},10/1),({12,7,19,18},6/1),({11,8,25,12},2/1),({11,8,24,13},5/1),({11,8,23,14},10/1),({11,8,22,15},15/1),({11,8,21,16},18/1),({11,8,20,17},16/1),({11,8,19,18},10/1),({10,9,26,11},1/1),({10,9,25,12},2/1),({10,9,24,13},5/1),({10,9,23,14},10/1),({10,9,22,15},13/1),({10,9,21,16},15/1),({10,9,20,17},14/1),({10,9,19,18},8/1)}, (9,0) => {}, (7,2) => {}, (9,1) => {({14,7,23,18},1/1),({14,7,22,19},1/1),({14,7,21,20},1/1),({13,8,26,15},1/1),({13,8,25,16},2/1),({13,8,24,17},4/1),({13,8,23,18},5/1),({13,8,22,19},5/1),({13,8,21,20},3/1),({12,9,27,14},1/1),({12,9,26,15},2/1),({12,9,25,16},5/1),({12,9,24,17},8/1),({12,9,23,18},10/1),({12,9,22,19},9/1),({12,9,21,20},6/1),({11,10,27,14},1/1),({11,10,26,15},3/1),({11,10,25,16},5/1),({11,10,24,17},8/1),({11,10,23,18},9/1),({11,10,22,19},8/1),({11,10,21,20},5/1)}, (10,0) => {}, (8,2) => {}, (10,1) => {({14,9,26,19},1/1),({14,9,25,20},1/1),({14,9,24,21},2/1),({14,9,23,22},1/1),({13,10,28,17},1/1),({13,10,27,18},2/1),({13,10,26,19},3/1),({13,10,25,20},4/1),({13,10,24,21},4/1),({13,10,23,22},2/1),({12,11,28,17},1/1),({12,11,27,18},2/1),({12,11,26,19},3/1),({12,11,25,20},4/1),({12,11,24,21},4/1),({12,11,23,22},2/1)}, (11,0) => {}, (9,2) => {}, (11,1) => {({14,11,28,21},1/1),({14,11,27,22},1/1),({14,11,26,23},1/1),({14,11,25,24},1/1),({13,12,29,20},1/1),({13,12,28,21},1/1),({13,12,27,22},1/1),({13,12,26,23},1/1),({13,12,25,24},1/1)}, (12,0) => {}, (10,2) => {}, (12,1) => {({14,13,29,24},1/1)}, (11,2) => {}, (12,2) => {}, (13,1) => {}, (14,1) => {}, (0,0) => {({1,0,1,0},1/1)}, (1,0) => {({3,0,4,1},1/1),({2,1,5,0},1/1),({2,1,4,1},1/1)}, (1,1) => {}, (2,0) => {({4,1,8,1},1/1),({4,1,7,2},1/1),({4,1,6,3},1/1),({4,1,5,4},1/1),({3,2,8,1},1/1),({3,2,7,2},1/1),({3,2,6,3},1/1),({3,2,5,4},1/1)}, (2,1) => {({7,0,8,5},1/1)}, (3,0) => {({5,2,11,2},1/1),({5,2,9,4},1/1),({5,2,8,5},1/1),({4,3,10,3},1/1),({4,3,9,4},1/1),({4,3,8,5},1/1),({4,3,7,6},1/1)}, (3,1) => {({9,0,9,8},1/1),({8,1,12,5},1/1),({8,1,11,6},1/1),({8,1,10,7},1/1),({8,1,9,8},1/1),({7,2,12,5},1/1),({7,2,11,6},2/1),({7,2,10,7},2/1),({7,2,9,8},1/1),({6,3,14,3},1/1),({6,3,13,4},1/1),({6,3,12,5},2/1),({6,3,11,6},2/1),({6,3,10,7},2/1),({6,3,9,8},1/1),({5,4,13,4},1/1),({5,4,12,5},1/1),({5,4,11,6},1/1),({5,4,10,7},1/1),({5,4,9,8},1/1)}, (4,0) => {}, (2,2) => {}, (4,1) => {({10,1,13,8},1/1),({10,1,12,9},1/1),({9,2,15,6},1/1),({9,2,14,7},1/1),({9,2,13,8},3/1),({9,2,12,9},3/1),({9,2,11,10},2/1),({8,3,16,5},1/1),({8,3,15,6},3/1),({8,3,14,7},5/1),({8,3,13,8},6/1),({8,3,12,9},6/1),({8,3,11,10},4/1),({7,4,17,4},1/1),({7,4,16,5},2/1),({7,4,15,6},5/1),({7,4,14,7},7/1),({7,4,13,8},9/1),({7,4,12,9},8/1),({7,4,11,10},5/1),({6,5,17,4},1/1),({6,5,16,5},2/1),({6,5,15,6},4/1),({6,5,14,7},6/1),({6,5,13,8},7/1),({6,5,12,9},6/1),({6,5,11,10},3/1)}, (5,0) => {}, (3,2) => {}, (5,1) => {({11,2,16,9},1/1),({11,2,15,10},1/1),({11,2,14,11},1/1),({11,2,13,12},1/1),({10,3,17,8},2/1),({10,3,16,9},3/1),({10,3,15,10},5/1),({10,3,14,11},5/1),({10,3,13,12},3/1),({9,4,19,6},1/1),({9,4,18,7},3/1),({9,4,17,8},6/1),({9,4,16,9},10/1),({9,4,15,10},12/1),({9,4,14,11},12/1),({9,4,13,12},7/1),({8,5,19,6},2/1),({8,5,18,7},5/1),({8,5,17,8},10/1),({8,5,16,9},15/1),({8,5,15,10},18/1),({8,5,14,11},15/1),({8,5,13,12},10/1),({7,6,20,5},1/1),({7,6,19,6},2/1),({7,6,18,7},5/1),({7,6,17,8},9/1),({7,6,16,9},12/1),({7,6,15,10},14/1),({7,6,14,11},13/1),({7,6,13,12},7/1)}, (6,0) => {}, (4,2) => {}}; --dw stands for dominant weights
