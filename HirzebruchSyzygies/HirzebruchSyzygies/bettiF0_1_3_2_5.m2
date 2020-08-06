A := QQ[t_0,t_1,t_2,t_3];
--tb stands for Total Betti numbers
tb1325 = new HashTable from {(6,1) => 9372/1, (7,0) => 792/1, (5,2) => 0, (7,1) => 25740/1, (8,0) => 0, (6,2) => 0, (8,1) => 37180/1, (9,0) => 0, (7,2) => 0, (10,0) => 0, (9,1) => 36036/1, (8,2) => 0, (11,0) => 0, (10,1) => 25116/1, (9,2) => 0, (12,0) => 0, (11,1) => 12740/1, (10,2) => 0, (13,0) => 0, (12,1) => 4620/1, (11,2) => 0, (14,0) => 0, (13,1) => 1140/1, (12,2) => 0, (15,0) => 0, (14,1) => 172/1, (13,2) => 0, (15,1) => 12/1, (14,2) => 0, (15,2) => 0, (16,1) => 0, (17,1) => 0, (0,0) => 8/1, (1,0) => 108/1, (1,1) => 0, (2,0) => 660/1, (2,1) => 0, (3,0) => 2380/1, (3,1) => 0, (4,0) => 5460/1, (2,2) => 0, (4,1) => 12/1, (5,0) => 7656/1, (3,2) => 0, (5,1) => 1012/1, (6,0) => 5016/1, (4,2) => 0};
--sb represents the betti numbers as sums of Schur functors
sb1325 = new HashTable from {(6,1) => {({13,2,23,15},1/1),({13,2,21,17},1/1),({13,2,19,19},1/1),({12,3,24,14},1/1),({12,3,23,15},2/1),({12,3,22,16},2/1),({12,3,21,17},2/1),({12,3,20,18},2/1),({12,3,19,19},1/1),({11,4,27,11},1/1),({11,4,26,12},2/1),({11,4,25,13},3/1),({11,4,24,14},4/1),({11,4,23,15},6/1),({11,4,22,16},6/1),({11,4,21,17},7/1),({11,4,20,18},4/1),({11,4,19,19},3/1),({10,5,27,11},1/1),({10,5,26,12},3/1),({10,5,25,13},6/1),({10,5,24,14},8/1),({10,5,23,15},10/1),({10,5,22,16},11/1),({10,5,21,17},10/1),({10,5,20,18},7/1),({10,5,19,19},2/1),({9,6,29,9},1/1),({9,6,28,10},1/1),({9,6,27,11},3/1),({9,6,26,12},4/1),({9,6,25,13},8/1),({9,6,24,14},10/1),({9,6,23,15},14/1),({9,6,22,16},13/1),({9,6,21,17},13/1),({9,6,20,18},7/1),({9,6,19,19},4/1),({8,7,28,10},1/1),({8,7,27,11},2/1),({8,7,26,12},3/1),({8,7,25,13},5/1),({8,7,24,14},7/1),({8,7,23,15},9/1),({8,7,22,16},9/1),({8,7,21,17},8/1),({8,7,20,18},6/1),({8,7,19,19},2/1)}, (7,0) => {({11,4,26,12},1/1),({11,4,24,14},1/1),({11,4,22,16},1/1),({11,4,20,18},1/1),({10,5,25,13},1/1),({10,5,24,14},1/1),({10,5,23,15},1/1),({10,5,22,16},1/1),({10,5,21,17},1/1),({10,5,20,18},1/1),({9,6,24,14},1/1),({9,6,23,15},1/1),({9,6,22,16},2/1),({9,6,21,17},1/1),({9,6,20,18},1/1),({8,7,23,15},1/1),({8,7,22,16},1/1),({8,7,21,17},1/1),({8,7,20,18},1/1)}, (5,2) => {}, (7,1) => {({14,3,26,17},1/1),({14,3,24,19},1/1),({14,3,23,20},1/1),({13,4,28,15},1/1),({13,4,27,16},2/1),({13,4,26,17},3/1),({13,4,25,18},4/1),({13,4,24,19},5/1),({13,4,23,20},4/1),({13,4,22,21},2/1),({12,5,30,13},1/1),({12,5,29,14},2/1),({12,5,28,15},5/1),({12,5,27,16},8/1),({12,5,26,17},12/1),({12,5,25,18},14/1),({12,5,24,19},15/1),({12,5,23,20},12/1),({12,5,22,21},7/1),({11,6,31,12},1/1),({11,6,30,13},3/1),({11,6,29,14},7/1),({11,6,28,15},12/1),({11,6,27,16},19/1),({11,6,26,17},26/1),({11,6,25,18},29/1),({11,6,24,19},29/1),({11,6,23,20},24/1),({11,6,22,21},13/1),({10,7,32,11},1/1),({10,7,31,12},2/1),({10,7,30,13},6/1),({10,7,29,14},11/1),({10,7,28,15},19/1),({10,7,27,16},27/1),({10,7,26,17},35/1),({10,7,25,18},39/1),({10,7,24,19},38/1),({10,7,23,20},30/1),({10,7,22,21},17/1),({9,8,32,11},1/1),({9,8,31,12},2/1),({9,8,30,13},5/1),({9,8,29,14},9/1),({9,8,28,15},15/1),({9,8,27,16},21/1),({9,8,26,17},26/1),({9,8,25,18},29/1),({9,8,24,19},28/1),({9,8,23,20},22/1),({9,8,22,21},12/1)}, (8,0) => {}, (6,2) => {}, (8,1) => {({15,4,28,20},1/1),({15,4,26,22},1/1),({15,4,24,24},1/1),({14,5,31,17},1/1),({14,5,30,18},2/1),({14,5,29,19},3/1),({14,5,28,20},5/1),({14,5,27,21},6/1),({14,5,26,22},5/1),({14,5,25,23},4/1),({14,5,24,24},2/1),({13,6,32,16},2/1),({13,6,31,17},4/1),({13,6,30,18},9/1),({13,6,29,19},13/1),({13,6,28,20},19/1),({13,6,27,21},19/1),({13,6,26,22},20/1),({13,6,25,23},12/1),({13,6,24,24},6/1),({12,7,34,14},1/1),({12,7,33,15},3/1),({12,7,32,16},7/1),({12,7,31,17},14/1),({12,7,30,18},23/1),({12,7,29,19},33/1),({12,7,28,20},41/1),({12,7,27,21},44/1),({12,7,26,22},40/1),({12,7,25,23},28/1),({12,7,24,24},10/1),({11,8,34,14},2/1),({11,8,33,15},5/1),({11,8,32,16},12/1),({11,8,31,17},21/1),({11,8,30,18},35/1),({11,8,29,19},45/1),({11,8,28,20},57/1),({11,8,27,21},58/1),({11,8,26,22},53/1),({11,8,25,23},35/1),({11,8,24,24},15/1),({10,9,35,13},1/1),({10,9,34,14},2/1),({10,9,33,15},5/1),({10,9,32,16},11/1),({10,9,31,17},18/1),({10,9,30,18},27/1),({10,9,29,19},37/1),({10,9,28,20},43/1),({10,9,27,21},44/1),({10,9,26,22},40/1),({10,9,25,23},27/1),({10,9,24,24},9/1)}, (9,0) => {}, (7,2) => {}, (10,0) => {}, (9,1) => {({16,5,29,24},1/1),({15,6,33,20},1/1),({15,6,32,21},2/1),({15,6,31,22},3/1),({15,6,30,23},4/1),({15,6,29,24},5/1),({15,6,28,25},4/1),({15,6,27,26},2/1),({14,7,35,18},1/1),({14,7,34,19},2/1),({14,7,33,20},6/1),({14,7,32,21},10/1),({14,7,31,22},15/1),({14,7,30,23},18/1),({14,7,29,24},19/1),({14,7,28,25},16/1),({14,7,27,26},9/1),({13,8,36,17},1/1),({13,8,35,18},4/1),({13,8,34,19},9/1),({13,8,33,20},17/1),({13,8,32,21},27/1),({13,8,31,22},37/1),({13,8,30,23},43/1),({13,8,29,24},43/1),({13,8,28,25},35/1),({13,8,27,26},20/1),({12,9,37,16},1/1),({12,9,36,17},3/1),({12,9,35,18},8/1),({12,9,34,19},16/1),({12,9,33,20},28/1),({12,9,32,21},41/1),({12,9,31,22},54/1),({12,9,30,23},61/1),({12,9,29,24},60/1),({12,9,28,25},48/1),({12,9,27,26},27/1),({11,10,37,16},1/1),({11,10,36,17},3/1),({11,10,35,18},7/1),({11,10,34,19},14/1),({11,10,33,20},23/1),({11,10,32,21},33/1),({11,10,31,22},42/1),({11,10,30,23},47/1),({11,10,29,24},46/1),({11,10,28,25},36/1),({11,10,27,26},20/1)}, (8,2) => {}, (11,0) => {}, (10,1) => {({17,6,29,29},1/1),({16,7,34,24},1/1),({16,7,33,25},2/1),({16,7,32,26},2/1),({16,7,31,27},2/1),({16,7,30,28},2/1),({16,7,29,29},1/1),({15,8,37,21},1/1),({15,8,36,22},2/1),({15,8,35,23},5/1),({15,8,34,24},7/1),({15,8,33,25},12/1),({15,8,32,26},12/1),({15,8,31,27},13/1),({15,8,30,28},8/1),({15,8,29,29},4/1),({14,9,38,20},1/1),({14,9,37,21},3/1),({14,9,36,22},8/1),({14,9,35,23},15/1),({14,9,34,24},22/1),({14,9,33,25},29/1),({14,9,32,26},33/1),({14,9,31,27},30/1),({14,9,30,28},21/1),({14,9,29,29},8/1),({13,10,39,19},1/1),({13,10,38,20},3/1),({13,10,37,21},8/1),({13,10,36,22},15/1),({13,10,35,23},26/1),({13,10,34,24},36/1),({13,10,33,25},46/1),({13,10,32,26},48/1),({13,10,31,27},45/1),({13,10,30,28},30/1),({13,10,29,29},12/1),({12,11,39,19},1/1),({12,11,38,20},3/1),({12,11,37,21},7/1),({12,11,36,22},13/1),({12,11,35,23},21/1),({12,11,34,24},30/1),({12,11,33,25},36/1),({12,11,32,26},38/1),({12,11,31,27},34/1),({12,11,30,28},24/1),({12,11,29,29},9/1)}, (9,2) => {}, (12,0) => {}, (11,1) => {({17,8,34,29},1/1),({17,8,33,30},1/1),({16,9,38,25},1/1),({16,9,37,26},2/1),({16,9,36,27},4/1),({16,9,35,28},5/1),({16,9,34,29},6/1),({16,9,33,30},5/1),({16,9,32,31},3/1),({15,10,40,23},1/1),({15,10,39,24},2/1),({15,10,38,25},5/1),({15,10,37,26},9/1),({15,10,36,27},14/1),({15,10,35,28},17/1),({15,10,34,29},18/1),({15,10,33,30},15/1),({15,10,32,31},9/1),({14,11,40,23},2/1),({14,11,39,24},5/1),({14,11,38,25},10/1),({14,11,37,26},17/1),({14,11,36,27},24/1),({14,11,35,28},28/1),({14,11,34,29},29/1),({14,11,33,30},24/1),({14,11,32,31},13/1),({13,12,41,22},1/1),({13,12,40,23},2/1),({13,12,39,24},5/1),({13,12,38,25},10/1),({13,12,37,26},15/1),({13,12,36,27},20/1),({13,12,35,28},24/1),({13,12,34,29},24/1),({13,12,33,30},19/1),({13,12,32,31},11/1)}, (10,2) => {}, (13,0) => {}, (12,1) => {({17,10,38,30},1/1),({17,10,37,31},1/1),({17,10,36,32},2/1),({17,10,35,33},1/1),({17,10,34,34},1/1),({16,11,41,27},1/1),({16,11,40,28},2/1),({16,11,39,29},4/1),({16,11,38,30},6/1),({16,11,37,31},7/1),({16,11,36,32},7/1),({16,11,35,33},5/1),({16,11,34,34},2/1),({15,12,42,26},1/1),({15,12,41,27},2/1),({15,12,40,28},5/1),({15,12,39,29},8/1),({15,12,38,30},12/1),({15,12,37,31},13/1),({15,12,36,32},13/1),({15,12,35,33},9/1),({15,12,34,34},4/1),({14,13,42,26},1/1),({14,13,41,27},3/1),({14,13,40,28},5/1),({14,13,39,29},8/1),({14,13,38,30},11/1),({14,13,37,31},12/1),({14,13,36,32},11/1),({14,13,35,33},8/1),({14,13,34,34},3/1)}, (11,2) => {}, (14,0) => {}, (13,1) => {({17,12,41,32},1/1),({17,12,40,33},1/1),({17,12,39,34},2/1),({17,12,38,35},2/1),({17,12,37,36},1/1),({16,13,43,30},1/1),({16,13,42,31},2/1),({16,13,41,32},3/1),({16,13,40,33},4/1),({16,13,39,34},5/1),({16,13,38,35},4/1),({16,13,37,36},2/1),({15,14,43,30},1/1),({15,14,42,31},2/1),({15,14,41,32},3/1),({15,14,40,33},4/1),({15,14,39,34},5/1),({15,14,38,35},4/1),({15,14,37,36},2/1)}, (12,2) => {}, (15,0) => {}, (14,1) => {({17,14,43,35},1/1),({17,14,42,36},1/1),({17,14,41,37},1/1),({17,14,40,38},1/1),({17,14,39,39},1/1),({16,15,44,34},1/1),({16,15,43,35},1/1),({16,15,42,36},1/1),({16,15,41,37},1/1),({16,15,40,38},1/1),({16,15,39,39},1/1)}, (13,2) => {}, (15,1) => {({17,16,44,39},1/1)}, (14,2) => {}, (15,2) => {}, (16,1) => {}, (17,1) => {}, (0,0) => {({1,0,3,0},1/1)}, (1,0) => {({3,0,7,1},1/1),({3,0,6,2},1/1),({3,0,5,3},1/1),({2,1,8,0},1/1),({2,1,7,1},1/1),({2,1,6,2},1/1),({2,1,5,3},1/1)}, (1,1) => {}, (2,0) => {({5,0,10,3},1/1),({5,0,9,4},1/1),({5,0,8,5},1/1),({5,0,7,6},1/1),({4,1,12,1},1/1),({4,1,11,2},2/1),({4,1,10,3},3/1),({4,1,9,4},3/1),({4,1,8,5},3/1),({4,1,7,6},2/1),({3,2,12,1},1/1),({3,2,11,2},2/1),({3,2,10,3},3/1),({3,2,9,4},3/1),({3,2,8,5},3/1),({3,2,7,6},2/1)}, (2,1) => {}, (3,0) => {({7,0,12,6},1/1),({7,0,10,8},1/1),({6,1,15,3},1/1),({6,1,14,4},2/1),({6,1,13,5},3/1),({6,1,12,6},4/1),({6,1,11,7},4/1),({6,1,10,8},3/1),({6,1,9,9},1/1),({5,2,16,2},1/1),({5,2,15,3},2/1),({5,2,14,4},5/1),({5,2,13,5},6/1),({5,2,12,6},8/1),({5,2,11,7},8/1),({5,2,10,8},6/1),({5,2,9,9},1/1),({4,3,16,2},1/1),({4,3,15,3},3/1),({4,3,14,4},4/1),({4,3,13,5},6/1),({4,3,12,6},8/1),({4,3,11,7},7/1),({4,3,10,8},5/1),({4,3,9,9},2/1)}, (3,1) => {}, (4,0) => {({8,1,17,6},1/1),({8,1,16,7},1/1),({8,1,15,8},2/1),({8,1,14,9},2/1),({8,1,13,10},2/1),({8,1,12,11},1/1),({7,2,19,4},1/1),({7,2,18,5},2/1),({7,2,17,6},4/1),({7,2,16,7},6/1),({7,2,15,8},8/1),({7,2,14,9},8/1),({7,2,13,10},7/1),({7,2,12,11},4/1),({6,3,19,4},2/1),({6,3,18,5},4/1),({6,3,17,6},7/1),({6,3,16,7},11/1),({6,3,15,8},14/1),({6,3,14,9},14/1),({6,3,13,10},12/1),({6,3,12,11},7/1),({5,4,20,3},1/1),({5,4,19,4},2/1),({5,4,18,5},4/1),({5,4,17,6},7/1),({5,4,16,7},10/1),({5,4,15,8},12/1),({5,4,14,9},12/1),({5,4,13,10},10/1),({5,4,12,11},6/1)}, (2,2) => {}, (4,1) => {({11,0,14,14},1/1)}, (5,0) => {({9,2,21,7},1/1),({9,2,20,8},1/1),({9,2,19,9},3/1),({9,2,18,10},2/1),({9,2,17,11},4/1),({9,2,16,12},2/1),({9,2,15,13},3/1),({8,3,22,6},1/1),({8,3,21,7},2/1),({8,3,20,8},4/1),({8,3,19,9},7/1),({8,3,18,10},9/1),({8,3,17,11},10/1),({8,3,16,12},10/1),({8,3,15,13},7/1),({8,3,14,14},2/1),({7,4,23,5},1/1),({7,4,22,6},2/1),({7,4,21,7},5/1),({7,4,20,8},8/1),({7,4,19,9},13/1),({7,4,18,10},15/1),({7,4,17,11},18/1),({7,4,16,12},15/1),({7,4,15,13},12/1),({7,4,14,14},3/1),({6,5,23,5},1/1),({6,5,22,6},2/1),({6,5,21,7},4/1),({6,5,20,8},7/1),({6,5,19,9},10/1),({6,5,18,10},13/1),({6,5,17,11},14/1),({6,5,16,12},13/1),({6,5,15,13},9/1),({6,5,14,14},3/1)}, (3,2) => {}, (5,1) => {({12,1,19,14},1/1),({11,2,19,14},1/1),({11,2,18,15},1/1),({10,3,23,10},1/1),({10,3,22,11},1/1),({10,3,21,12},1/1),({10,3,20,13},1/1),({10,3,19,14},1/1),({10,3,18,15},1/1),({10,3,17,16},1/1),({9,4,22,11},1/1),({9,4,21,12},1/1),({9,4,20,13},1/1),({9,4,19,14},1/1),({9,4,18,15},1/1),({9,4,17,16},1/1),({8,5,21,12},1/1),({8,5,20,13},1/1),({8,5,19,14},1/1),({8,5,18,15},1/1),({7,6,20,13},1/1),({7,6,19,14},1/1)}, (6,0) => {({10,3,24,9},1/1),({10,3,23,10},1/1),({10,3,22,11},2/1),({10,3,21,12},2/1),({10,3,20,13},2/1),({10,3,19,14},2/1),({10,3,18,15},2/1),({10,3,17,16},1/1),({9,4,24,9},1/1),({9,4,23,10},2/1),({9,4,22,11},4/1),({9,4,21,12},5/1),({9,4,20,13},6/1),({9,4,19,14},6/1),({9,4,18,15},5/1),({9,4,17,16},3/1),({8,5,26,7},1/1),({8,5,25,8},1/1),({8,5,24,9},2/1),({8,5,23,10},4/1),({8,5,22,11},6/1),({8,5,21,12},8/1),({8,5,20,13},11/1),({8,5,19,14},10/1),({8,5,18,15},8/1),({8,5,17,16},5/1),({7,6,25,8},1/1),({7,6,24,9},2/1),({7,6,23,10},3/1),({7,6,22,11},5/1),({7,6,21,12},7/1),({7,6,20,13},8/1),({7,6,19,14},8/1),({7,6,18,15},7/1),({7,6,17,16},3/1)}, (4,2) => {}};
--dw stands for dominant weights
dw1325 = new HashTable from {(5,2) => {}, (7,0) => {({11,4,26,12},1/1)}, (6,1) => {({13,2,23,15},1/1),({12,3,24,14},1/1),({11,4,27,11},1/1),({9,6,29,9},1/1)}, (7,1) => {({14,3,26,17},1/1),({13,4,28,15},1/1),({12,5,30,13},1/1),({11,6,31,12},1/1),({10,7,32,11},1/1)}, (8,0) => {}, (6,2) => {}, (7,2) => {}, (9,0) => {}, (8,1) => {({15,4,28,20},1/1),({14,5,31,17},1/1),({13,6,32,16},2/1),({12,7,34,14},1/1),({10,9,35,13},1/1)}, (8,2) => {}, (9,1) => {({16,5,29,24},1/1),({15,6,33,20},1/1),({14,7,35,18},1/1),({13,8,36,17},1/1),({12,9,37,16},1/1)}, (10,0) => {}, (11,0) => {}, (10,1) => {({17,6,29,29},1/1),({16,7,34,24},1/1),({15,8,37,21},1/1),({14,9,38,20},1/1),({13,10,39,19},1/1)}, (9,2) => {}, (12,0) => {}, (11,1) => {({17,8,34,29},1/1),({16,9,38,25},1/1),({15,10,40,23},1/1),({13,12,41,22},1/1)}, (10,2) => {}, (13,0) => {}, (12,1) => {({17,10,38,30},1/1),({16,11,41,27},1/1),({15,12,42,26},1/1)}, (11,2) => {}, (14,0) => {}, (13,1) => {({17,12,41,32},1/1),({16,13,43,30},1/1)}, (12,2) => {}, (13,2) => {}, (14,1) => {({17,14,43,35},1/1),({16,15,44,34},1/1)}, (15,0) => {}, (14,2) => {}, (15,1) => {({17,16,44,39},1/1)}, (16,1) => {}, (15,2) => {}, (17,1) => {}, (0,0) => {({1,0,3,0},1/1)}, (1,0) => {({3,0,7,1},1/1),({2,1,8,0},1/1)}, (2,0) => {({5,0,10,3},1/1),({4,1,12,1},1/1)}, (1,1) => {}, (3,0) => {({7,0,12,6},1/1),({6,1,15,3},1/1),({5,2,16,2},1/1)}, (2,1) => {}, (2,2) => {}, (4,0) => {({8,1,17,6},1/1),({7,2,19,4},1/1),({5,4,20,3},1/1)}, (3,1) => {}, (3,2) => {}, (5,0) => {({9,2,21,7},1/1),({8,3,22,6},1/1),({7,4,23,5},1/1)}, (4,1) => {({11,0,14,14},1/1)}, (4,2) => {}, (6,0) => {({10,3,24,9},1/1),({8,5,26,7},1/1)}, (5,1) => {({12,1,19,14},1/1),({10,3,23,10},1/1)}};
--dmw stands for dominant monomial weights
dmw1325 = new HashTable from {(5,2) => {}, (7,0) => {}, (6,1) => {{13,2,23,15},{12,3,24,14}}, (7,1) => {{13,4,28,15},{14,3,26,17}}, (8,0) => {}, (6,2) => {}, (7,2) => {}, (9,0) => {}, (8,1) => {{14,5,31,17},{15,4,28,20},{13,6,32,16}}, (8,2) => {}, (9,1) => {{16,5,29,24},{14,7,35,18},{15,6,33,20}}, (10,0) => {}, (11,0) => {}, (10,1) => {{14,9,38,20},{16,7,34,24},{17,6,29,29},{15,8,37,21}}, (9,2) => {}, (12,0) => {}, (11,1) => {{17,8,34,29},{15,10,40,23},{16,9,38,25}}, (10,2) => {}, (13,0) => {}, (12,1) => {{17,10,38,30},{15,12,42,26},{16,11,41,27}}, (11,2) => {}, (14,0) => {}, (13,1) => {{16,13,43,30},{17,12,41,32}}, (12,2) => {}, (13,2) => {{16,15,34,44}}, (14,1) => {{17,14,43,35},{16,15,44,34}}, (15,0) => {}, (14,2) => {{17,16,39,44}}, (15,1) => {{17,16,44,39}}, (16,1) => {}, (15,2) => {{17,18,44,44}}, (17,1) => {}, (0,0) => {}, (1,0) => {}, (2,0) => {}, (1,1) => {}, (3,0) => {}, (2,1) => {}, (2,2) => {}, (4,0) => {}, (3,1) => {}, (3,2) => {}, (5,0) => {}, (4,1) => {{11,0,14,14}}, (4,2) => {}, (6,0) => {}, (5,1) => {{12,1,19,14}}};
end;