A := QQ[t_0,t_1,t_2,t_3];
--tb stands for Total Betti numbers
tb0324 = new HashTable from {(5,2) => 0, (7,0) => 0, (6,1) => 7920/1, (7,1) => 7524/1, (8,0) => 0, (6,2) => 0, (7,2) => 0, (9,0) => 0, (8,1) => 5060/1, (8,2) => 0, (10,0) => 0, (9,1) => 2376/1, (10,1) => 744/1, (11,0) => 0, (9,2) => 0, (11,1) => 140/1, (12,0) => 0, (10,2) => 0, (12,1) => 12/1, (11,2) => 0, (12,2) => 0, (13,1) => 0, (14,1) => 0, (0,0) => 4/1, (1,0) => 36/1, (2,0) => 132/1, (1,1) => 12/1, (3,0) => 220/1, (2,1) => 132/1, (2,2) => 0, (4,0) => 0, (3,1) => 660/1, (3,2) => 0, (5,0) => 0, (4,1) => 2772/1, (4,2) => 0, (6,0) => 0, (5,1) => 5808/1};
--sb represents the betti numbers as sums of Schur functors
sb0324 = new HashTable from {(5,2) => {}, (7,0) => {}, (6,1) => {({11,3,20,11},1/1),({11,3,19,12},1/1),({11,3,18,13},2/1),({11,3,17,14},2/1),({11,3,16,15},1/1),({10,4,22,9},1/1),({10,4,21,10},2/1),({10,4,20,11},5/1),({10,4,19,12},7/1),({10,4,18,13},9/1),({10,4,17,14},8/1),({10,4,16,15},5/1),({9,5,23,8},1/1),({9,5,22,9},3/1),({9,5,21,10},7/1),({9,5,20,11},12/1),({9,5,19,12},17/1),({9,5,18,13},19/1),({9,5,17,14},17/1),({9,5,16,15},10/1),({8,6,24,7},1/1),({8,6,23,8},2/1),({8,6,22,9},6/1),({8,6,21,10},10/1),({8,6,20,11},17/1),({8,6,19,12},21/1),({8,6,18,13},24/1),({8,6,17,14},20/1),({8,6,16,15},12/1),({7,7,23,8},1/1),({7,7,22,9},2/1),({7,7,21,10},5/1),({7,7,20,11},7/1),({7,7,19,12},10/1),({7,7,18,13},10/1),({7,7,17,14},9/1),({7,7,16,15},5/1)}, (7,1) => {({12,4,22,13},1/1),({12,4,21,14},1/1),({12,4,20,15},2/1),({12,4,19,16},2/1),({12,4,18,17},1/1),({11,5,23,12},2/1),({11,5,22,13},4/1),({11,5,21,14},7/1),({11,5,20,15},9/1),({11,5,19,16},8/1),({11,5,18,17},5/1),({10,6,25,10},1/1),({10,6,24,11},3/1),({10,6,23,12},6/1),({10,6,22,13},12/1),({10,6,21,14},17/1),({10,6,20,15},20/1),({10,6,19,16},18/1),({10,6,18,17},11/1),({9,7,25,10},1/1),({9,7,24,11},4/1),({9,7,23,12},9/1),({9,7,22,13},15/1),({9,7,21,14},21/1),({9,7,20,15},23/1),({9,7,19,16},21/1),({9,7,18,17},12/1),({8,8,26,9},1/1),({8,8,25,10},1/1),({8,8,24,11},3/1),({8,8,23,12},5/1),({8,8,22,13},8/1),({8,8,21,14},10/1),({8,8,20,15},12/1),({8,8,19,16},10/1),({8,8,18,17},6/1)}, (8,0) => {}, (6,2) => {}, (7,2) => {}, (9,0) => {}, (8,1) => {({13,5,23,16},1/1),({13,5,22,17},1/1),({13,5,21,18},1/1),({13,5,20,19},1/1),({12,6,25,14},1/1),({12,6,24,15},2/1),({12,6,23,16},5/1),({12,6,22,17},6/1),({12,6,21,18},6/1),({12,6,20,19},4/1),({11,7,26,13},1/1),({11,7,25,14},4/1),({11,7,24,15},8/1),({11,7,23,16},12/1),({11,7,22,17},15/1),({11,7,21,18},14/1),({11,7,20,19},8/1),({10,8,27,12},1/1),({10,8,26,13},2/1),({10,8,25,14},6/1),({10,8,24,15},10/1),({10,8,23,16},16/1),({10,8,22,17},18/1),({10,8,21,18},17/1),({10,8,20,19},10/1),({9,9,26,13},2/1),({9,9,25,14},3/1),({9,9,24,15},6/1),({9,9,23,16},8/1),({9,9,22,17},9/1),({9,9,21,18},7/1),({9,9,20,19},5/1)}, (8,2) => {}, (10,0) => {}, (9,1) => {({14,6,23,20},1/1),({13,7,26,17},1/1),({13,7,25,18},2/1),({13,7,24,19},3/1),({13,7,23,20},3/1),({13,7,22,21},2/1),({12,8,27,16},1/1),({12,8,26,17},3/1),({12,8,25,18},6/1),({12,8,24,19},8/1),({12,8,23,20},8/1),({12,8,22,21},5/1),({11,9,28,15},1/1),({11,9,27,16},3/1),({11,9,26,17},6/1),({11,9,25,18},9/1),({11,9,24,19},11/1),({11,9,23,20},10/1),({11,9,22,21},6/1),({10,10,27,16},1/1),({10,10,26,17},2/1),({10,10,25,18},4/1),({10,10,24,19},5/1),({10,10,23,20},5/1),({10,10,22,21},3/1)}, (10,1) => {({14,8,26,21},1/1),({14,8,25,22},1/1),({14,8,24,23},1/1),({13,9,28,19},1/1),({13,9,27,20},2/1),({13,9,26,21},3/1),({13,9,25,22},3/1),({13,9,24,23},2/1),({12,10,28,19},2/1),({12,10,27,20},3/1),({12,10,26,21},5/1),({12,10,25,22},5/1),({12,10,24,23},3/1),({11,11,29,18},1/1),({11,11,28,19},1/1),({11,11,27,20},2/1),({11,11,26,21},2/1),({11,11,25,22},2/1),({11,11,24,23},1/1)}, (11,0) => {}, (9,2) => {}, (11,1) => {({14,10,28,23},1/1),({14,10,27,24},1/1),({14,10,26,25},1/1),({13,11,29,22},1/1),({13,11,28,23},1/1),({13,11,27,24},1/1),({13,11,26,25},1/1),({12,12,29,22},1/1),({12,12,28,23},1/1),({12,12,27,24},1/1),({12,12,26,25},1/1)}, (12,0) => {}, (10,2) => {}, (12,1) => {({14,12,29,26},1/1)}, (11,2) => {}, (12,2) => {}, (13,1) => {}, (14,1) => {}, (0,0) => {({0,0,3,0},1/1)}, (1,0) => {({2,0,6,1},1/1),({2,0,5,2},1/1),({2,0,4,3},1/1)}, (2,0) => {({4,0,8,3},1/1),({4,0,7,4},1/1),({4,0,6,5},1/1),({3,1,9,2},1/1),({3,1,8,3},1/1),({3,1,7,4},1/1),({3,1,6,5},1/1),({2,2,8,3},1/1),({2,2,7,4},1/1),({2,2,6,5},1/1)}, (1,1) => {({2,2,11,0},1/1)}, (3,0) => {({6,0,9,6},1/1),({5,1,11,4},1/1),({5,1,10,5},1/1),({5,1,9,6},1/1),({5,1,8,7},1/1),({4,2,11,4},1/1),({4,2,10,5},1/1),({4,2,9,6},2/1),({4,2,8,7},1/1),({3,3,12,3},1/1),({3,3,10,5},1/1),({3,3,9,6},1/1)}, (2,1) => {({4,2,14,1},1/1),({4,2,13,2},1/1),({4,2,12,3},1/1),({4,2,11,4},1/1)}, (2,2) => {}, (4,0) => {}, (3,1) => {({6,2,16,3},1/1),({6,2,15,4},1/1),({6,2,14,5},2/1),({6,2,13,6},1/1),({6,2,12,7},1/1),({5,3,17,2},1/1),({5,3,16,3},1/1),({5,3,15,4},2/1),({5,3,14,5},2/1),({5,3,13,6},2/1),({5,3,12,7},1/1),({5,3,11,8},1/1),({4,4,16,3},1/1),({4,4,15,4},1/1),({4,4,14,5},2/1),({4,4,13,6},1/1),({4,4,12,7},1/1)}, (3,2) => {}, (5,0) => {}, (4,1) => {({9,1,13,10},1/1),({8,2,17,6},1/1),({8,2,16,7},1/1),({8,2,15,8},2/1),({8,2,14,9},2/1),({8,2,13,10},2/1),({8,2,12,11},1/1),({7,3,19,4},1/1),({7,3,18,5},2/1),({7,3,17,6},3/1),({7,3,16,7},5/1),({7,3,15,8},6/1),({7,3,14,9},6/1),({7,3,13,10},5/1),({7,3,12,11},3/1),({6,4,19,4},1/1),({6,4,18,5},2/1),({6,4,17,6},5/1),({6,4,16,7},6/1),({6,4,15,8},8/1),({6,4,14,9},7/1),({6,4,13,10},6/1),({6,4,12,11},3/1),({5,5,20,3},1/1),({5,5,19,4},1/1),({5,5,18,5},2/1),({5,5,17,6},3/1),({5,5,16,7},4/1),({5,5,15,8},4/1),({5,5,14,9},5/1),({5,5,13,10},3/1),({5,5,12,11},2/1)}, (4,2) => {}, (6,0) => {}, (5,1) => {({10,2,17,10},1/1),({10,2,16,11},1/1),({10,2,15,12},1/1),({10,2,14,13},1/1),({9,3,20,7},1/1),({9,3,19,8},2/1),({9,3,18,9},3/1),({9,3,17,10},5/1),({9,3,16,11},6/1),({9,3,15,12},5/1),({9,3,14,13},4/1),({8,4,21,6},1/1),({8,4,20,7},2/1),({8,4,19,8},6/1),({8,4,18,9},9/1),({8,4,17,10},12/1),({8,4,16,11},13/1),({8,4,15,12},12/1),({8,4,14,13},6/1),({7,5,22,5},1/1),({7,5,21,6},2/1),({7,5,20,7},5/1),({7,5,19,8},8/1),({7,5,18,9},13/1),({7,5,17,10},16/1),({7,5,16,11},17/1),({7,5,15,12},14/1),({7,5,14,13},8/1),({6,6,21,6},1/1),({6,6,20,7},2/1),({6,6,19,8},4/1),({6,6,18,9},5/1),({6,6,17,10},8/1),({6,6,16,11},7/1),({6,6,15,12},6/1),({6,6,14,13},4/1)}};
--dw stands for dominant weights
dw0324 = new HashTable from {(6,1) => {({11,3,20,11},1/1),({10,4,22,9},1/1),({9,5,23,8},1/1),({8,6,24,7},1/1)}, (7,0) => {}, (5,2) => {}, (7,1) => {({12,4,22,13},1/1),({11,5,23,12},2/1),({10,6,25,10},1/1),({8,8,26,9},1/1)}, (8,0) => {}, (6,2) => {}, (8,1) => {({13,5,23,16},1/1),({12,6,25,14},1/1),({11,7,26,13},1/1),({10,8,27,12},1/1)}, (9,0) => {}, (7,2) => {}, (9,1) => {({14,6,23,20},1/1),({13,7,26,17},1/1),({12,8,27,16},1/1),({11,9,28,15},1/1)}, (10,0) => {}, (8,2) => {}, (10,1) => {({14,8,26,21},1/1),({13,9,28,19},1/1),({11,11,29,18},1/1)}, (11,0) => {}, (9,2) => {}, (11,1) => {({14,10,28,23},1/1),({13,11,29,22},1/1)}, (12,0) => {}, (10,2) => {}, (12,1) => {({14,12,29,26},1/1)}, (11,2) => {}, (12,2) => {}, (13,1) => {}, (14,1) => {}, (0,0) => {({0,0,3,0},1/1)}, (1,0) => {({2,0,6,1},1/1)}, (1,1) => {({2,2,11,0},1/1)}, (2,0) => {({4,0,8,3},1/1),({3,1,9,2},1/1)}, (2,1) => {({4,2,14,1},1/1)}, (3,0) => {({6,0,9,6},1/1),({5,1,11,4},1/1),({3,3,12,3},1/1)}, (3,1) => {({6,2,16,3},1/1),({5,3,17,2},1/1)}, (4,0) => {}, (2,2) => {}, (4,1) => {({9,1,13,10},1/1),({8,2,17,6},1/1),({7,3,19,4},1/1),({5,5,20,3},1/1)}, (5,0) => {}, (3,2) => {}, (5,1) => {({10,2,17,10},1/1),({9,3,20,7},1/1),({8,4,21,6},1/1),({7,5,22,5},1/1)}, (6,0) => {}, (4,2) => {}};
end;