A := QQ[t_0,t_1,t_2,t_3];
--tb stands for Total Betti numbers
tb0325 = new HashTable from {(6,1) => 46332/1, (7,0) => 0, (5,2) => 0, (7,1) => 65780/1, (8,0) => 0, (6,2) => 0, (8,1) => 70070/1, (9,0) => 0, (7,2) => 0, (10,0) => 0, (9,1) => 56628/1, (8,2) => 0, (11,0) => 0, (10,1) => 34580/1, (9,2) => 0, (12,0) => 0, (11,1) => 15652/1, (10,2) => 0, (13,0) => 0, (12,1) => 5040/1, (11,2) => 0, (14,0) => 0, (13,1) => 1060/1, (12,2) => 0, (15,0) => 0, (14,1) => 116/1, (13,2) => 0, (15,1) => 0, (14,2) => 0, (15,2) => 1/1, (16,2) => 0, (0,0) => 4/1, (1,0) => 45/1, (1,1) => 14/1, (2,0) => 210/1, (2,1) => 195/1, (3,0) => 455/1, (3,1) => 1260/1, (4,0) => 0, (2,2) => 0, (4,1) => 8008/1, (5,0) => 0, (3,2) => 0, (5,1) => 23660/1, (6,0) => 0, (4,2) => 0};
--sb represents the betti numbers as sums of Schur functors
sb0325 = new HashTable from {(6,1) => {({12,2,23,15},2/1),({12,2,22,16},1/1),({12,2,21,17},2/1),({12,2,20,18},1/1),({12,2,19,19},1/1),({11,3,27,11},1/1),({11,3,26,12},2/1),({11,3,25,13},4/1),({11,3,24,14},7/1),({11,3,23,15},9/1),({11,3,22,16},11/1),({11,3,21,17},10/1),({11,3,20,18},8/1),({11,3,19,19},3/1),({10,4,29,9},1/1),({10,4,28,10},2/1),({10,4,27,11},6/1),({10,4,26,12},11/1),({10,4,25,13},19/1),({10,4,24,14},25/1),({10,4,23,15},33/1),({10,4,22,16},33/1),({10,4,21,17},32/1),({10,4,20,18},20/1),({10,4,19,19},9/1),({9,5,30,8},1/1),({9,5,29,9},3/1),({9,5,28,10},8/1),({9,5,27,11},15/1),({9,5,26,12},27/1),({9,5,25,13},39/1),({9,5,24,14},53/1),({9,5,23,15},60/1),({9,5,22,16},64/1),({9,5,21,17},54/1),({9,5,20,18},39/1),({9,5,19,19},12/1),({8,6,31,7},1/1),({8,6,30,8},2/1),({8,6,29,9},6/1),({8,6,28,10},11/1),({8,6,27,11},22/1),({8,6,26,12},32/1),({8,6,25,13},49/1),({8,6,24,14},59/1),({8,6,23,15},71/1),({8,6,22,16},67/1),({8,6,21,17},62/1),({8,6,20,18},38/1),({8,6,19,19},17/1),({7,7,30,8},1/1),({7,7,29,9},2/1),({7,7,28,10},6/1),({7,7,27,11},9/1),({7,7,26,12},16/1),({7,7,25,13},20/1),({7,7,24,14},28/1),({7,7,23,15},29/1),({7,7,22,16},32/1),({7,7,21,17},24/1),({7,7,20,18},20/1),({7,7,19,19},4/1)}, (7,0) => {}, (5,2) => {}, (7,1) => {({13,3,27,16},1/1),({13,3,26,17},2/1),({13,3,25,18},2/1),({13,3,24,19},3/1),({13,3,23,20},3/1),({13,3,22,21},1/1),({12,4,30,13},1/1),({12,4,29,14},2/1),({12,4,28,15},5/1),({12,4,27,16},9/1),({12,4,26,17},13/1),({12,4,25,18},16/1),({12,4,24,19},17/1),({12,4,23,20},14/1),({12,4,22,21},8/1),({11,5,31,12},2/1),({11,5,30,13},5/1),({11,5,29,14},12/1),({11,5,28,15},21/1),({11,5,27,16},32/1),({11,5,26,17},42/1),({11,5,25,18},49/1),({11,5,24,19},48/1),({11,5,23,20},39/1),({11,5,22,21},22/1),({10,6,33,10},1/1),({10,6,32,11},3/1),({10,6,31,12},7/1),({10,6,30,13},16/1),({10,6,29,14},29/1),({10,6,28,15},46/1),({10,6,27,16},65/1),({10,6,26,17},82/1),({10,6,25,18},89/1),({10,6,24,19},86/1),({10,6,23,20},68/1),({10,6,22,21},37/1),({9,7,33,10},1/1),({9,7,32,11},4/1),({9,7,31,12},10/1),({9,7,30,13},20/1),({9,7,29,14},35/1),({9,7,28,15},53/1),({9,7,27,16},73/1),({9,7,26,17},88/1),({9,7,25,18},96/1),({9,7,24,19},90/1),({9,7,23,20},71/1),({9,7,22,21},39/1),({8,8,34,9},1/1),({8,8,33,10},1/1),({8,8,32,11},3/1),({8,8,31,12},6/1),({8,8,30,13},11/1),({8,8,29,14},17/1),({8,8,28,15},27/1),({8,8,27,16},34/1),({8,8,26,17},41/1),({8,8,25,18},44/1),({8,8,24,19},41/1),({8,8,23,20},31/1),({8,8,22,21},18/1)}, (8,0) => {}, (6,2) => {}, (8,1) => {({14,4,30,18},1/1),({14,4,29,19},1/1),({14,4,28,20},3/1),({14,4,27,21},3/1),({14,4,26,22},3/1),({14,4,25,23},2/1),({14,4,24,24},2/1),({13,5,32,16},1/1),({13,5,31,17},4/1),({13,5,30,18},7/1),({13,5,29,19},12/1),({13,5,28,20},16/1),({13,5,27,21},19/1),({13,5,26,22},17/1),({13,5,25,23},13/1),({13,5,24,24},4/1),({12,6,34,14},1/1),({12,6,33,15},3/1),({12,6,32,16},9/1),({12,6,31,17},16/1),({12,6,30,18},29/1),({12,6,29,19},40/1),({12,6,28,20},52/1),({12,6,27,21},53/1),({12,6,26,22},51/1),({12,6,25,23},33/1),({12,6,24,24},14/1),({11,7,35,13},1/1),({11,7,34,14},4/1),({11,7,33,15},11/1),({11,7,32,16},21/1),({11,7,31,17},39/1),({11,7,30,18},58/1),({11,7,29,19},79/1),({11,7,28,20},92/1),({11,7,27,21},99/1),({11,7,26,22},84/1),({11,7,25,23},61/1),({11,7,24,24},20/1),({10,8,36,12},1/1),({10,8,35,13},2/1),({10,8,34,14},7/1),({10,8,33,15},14/1),({10,8,32,16},29/1),({10,8,31,17},45/1),({10,8,30,18},69/1),({10,8,29,19},86/1),({10,8,28,20},104/1),({10,8,27,21},101/1),({10,8,26,22},93/1),({10,8,25,23},59/1),({10,8,24,24},25/1),({9,9,35,13},2/1),({9,9,34,14},3/1),({9,9,33,15},8/1),({9,9,32,16},13/1),({9,9,31,17},23/1),({9,9,30,18},29/1),({9,9,29,19},42/1),({9,9,28,20},43/1),({9,9,27,21},48/1),({9,9,26,22},38/1),({9,9,25,23},30/1),({9,9,24,24},6/1)}, (9,0) => {}, (7,2) => {}, (10,0) => {}, (9,1) => {({15,5,32,21},1/1),({15,5,31,22},1/1),({15,5,30,23},2/1),({15,5,29,24},3/1),({15,5,28,25},2/1),({15,5,27,26},1/1),({14,6,34,19},1/1),({14,6,33,20},4/1),({14,6,32,21},7/1),({14,6,31,22},11/1),({14,6,30,23},14/1),({14,6,29,24},15/1),({14,6,28,25},13/1),({14,6,27,26},7/1),({13,7,36,17},1/1),({13,7,35,18},4/1),({13,7,34,19},9/1),({13,7,33,20},17/1),({13,7,32,21},28/1),({13,7,31,22},38/1),({13,7,30,23},45/1),({13,7,29,24},45/1),({13,7,28,25},37/1),({13,7,27,26},21/1),({12,8,37,16},1/1),({12,8,36,17},4/1),({12,8,35,18},11/1),({12,8,34,19},22/1),({12,8,33,20},38/1),({12,8,32,21},56/1),({12,8,31,22},73/1),({12,8,30,23},82/1),({12,8,29,24},80/1),({12,8,28,25},64/1),({12,8,27,26},36/1),({11,9,38,15},1/1),({11,9,37,16},3/1),({11,9,36,17},8/1),({11,9,35,18},16/1),({11,9,34,19},30/1),({11,9,33,20},47/1),({11,9,32,21},66/1),({11,9,31,22},82/1),({11,9,30,23},91/1),({11,9,29,24},87/1),({11,9,28,25},69/1),({11,9,27,26},38/1),({10,10,37,16},1/1),({10,10,36,17},3/1),({10,10,35,18},7/1),({10,10,34,19},13/1),({10,10,33,20},21/1),({10,10,32,21},29/1),({10,10,31,22},36/1),({10,10,30,23},39/1),({10,10,29,24},38/1),({10,10,28,25},29/1),({10,10,27,26},16/1)}, (8,2) => {}, (11,0) => {}, (10,1) => {({16,6,33,25},1/1),({16,6,32,26},1/1),({16,6,31,27},1/1),({16,6,30,28},1/1),({16,6,29,29},1/1),({15,7,36,22},1/1),({15,7,35,23},2/1),({15,7,34,24},5/1),({15,7,33,25},7/1),({15,7,32,26},9/1),({15,7,31,27},8/1),({15,7,30,28},7/1),({15,7,29,29},2/1),({14,8,37,21},3/1),({14,8,36,22},6/1),({14,8,35,23},13/1),({14,8,34,24},19/1),({14,8,33,25},27/1),({14,8,32,26},28/1),({14,8,31,27},29/1),({14,8,30,28},18/1),({14,8,29,29},8/1),({13,9,39,19},1/1),({13,9,38,20},4/1),({13,9,37,21},8/1),({13,9,36,22},18/1),({13,9,35,23},29/1),({13,9,34,24},42/1),({13,9,33,25},51/1),({13,9,32,26},58/1),({13,9,31,27},49/1),({13,9,30,28},37/1),({13,9,29,29},12/1),({12,10,39,19},2/1),({12,10,38,20},5/1),({12,10,37,21},13/1),({12,10,36,22},22/1),({12,10,35,23},37/1),({12,10,34,24},48/1),({12,10,33,25},61/1),({12,10,32,26},60/1),({12,10,31,27},57/1),({12,10,30,28},36/1),({12,10,29,29},16/1),({11,11,40,18},1/1),({11,11,39,19},1/1),({11,11,38,20},4/1),({11,11,37,21},6/1),({11,11,36,22},12/1),({11,11,35,23},15/1),({11,11,34,24},25/1),({11,11,33,25},25/1),({11,11,32,26},30/1),({11,11,31,27},23/1),({11,11,30,28},20/1),({11,11,29,29},3/1)}, (9,2) => {}, (12,0) => {}, (11,1) => {({17,7,33,30},1/1),({16,8,37,26},1/1),({16,8,36,27},2/1),({16,8,35,28},3/1),({16,8,34,29},4/1),({16,8,33,30},3/1),({16,8,32,31},2/1),({15,9,39,24},1/1),({15,9,38,25},3/1),({15,9,37,26},6/1),({15,9,36,27},10/1),({15,9,35,28},13/1),({15,9,34,29},14/1),({15,9,33,30},12/1),({15,9,32,31},7/1),({14,10,40,23},2/1),({14,10,39,24},5/1),({14,10,38,25},10/1),({14,10,37,26},17/1),({14,10,36,27},24/1),({14,10,35,28},28/1),({14,10,34,29},29/1),({14,10,33,30},24/1),({14,10,32,31},13/1),({13,11,41,22},1/1),({13,11,40,23},3/1),({13,11,39,24},7/1),({13,11,38,25},13/1),({13,11,37,26},21/1),({13,11,36,27},28/1),({13,11,35,28},33/1),({13,11,34,29},33/1),({13,11,33,30},27/1),({13,11,32,31},15/1),({12,12,41,22},1/1),({12,12,40,23},2/1),({12,12,39,24},4/1),({12,12,38,25},8/1),({12,12,37,26},11/1),({12,12,36,27},14/1),({12,12,35,28},16/1),({12,12,34,29},15/1),({12,12,33,30},12/1),({12,12,32,31},7/1)}, (10,2) => {}, (13,0) => {}, (12,1) => {({17,9,37,31},1/1),({17,9,36,32},1/1),({17,9,35,33},1/1),({16,10,40,28},1/1),({16,10,39,29},2/1),({16,10,38,30},4/1),({16,10,37,31},4/1),({16,10,36,32},5/1),({16,10,35,33},3/1),({16,10,34,34},2/1),({15,11,41,27},2/1),({15,11,40,28},3/1),({15,11,39,29},7/1),({15,11,38,30},9/1),({15,11,37,31},12/1),({15,11,36,32},10/1),({15,11,35,33},9/1),({15,11,34,34},2/1),({14,12,42,26},2/1),({14,12,41,27},3/1),({14,12,40,28},7/1),({14,12,39,29},9/1),({14,12,38,30},14/1),({14,12,37,31},13/1),({14,12,36,32},14/1),({14,12,35,33},8/1),({14,12,34,34},5/1),({13,13,41,27},2/1),({13,13,40,28},2/1),({13,13,39,29},5/1),({13,13,38,30},5/1),({13,13,37,31},8/1),({13,13,36,32},5/1),({13,13,35,33},6/1)}, (11,2) => {}, (14,0) => {}, (13,1) => {({17,11,40,33},1/1),({17,11,39,34},1/1),({17,11,38,35},1/1),({17,11,37,36},1/1),({16,12,42,31},1/1),({16,12,41,32},2/1),({16,12,40,33},2/1),({16,12,39,34},3/1),({16,12,38,35},3/1),({16,12,37,36},1/1),({15,13,43,30},1/1),({15,13,42,31},2/1),({15,13,41,32},3/1),({15,13,40,33},4/1),({15,13,39,34},5/1),({15,13,38,35},4/1),({15,13,37,36},2/1),({14,14,43,30},1/1),({14,14,42,31},1/1),({14,14,41,32},1/1),({14,14,40,33},2/1),({14,14,39,34},2/1),({14,14,38,35},1/1),({14,14,37,36},1/1)}, (12,2) => {}, (15,0) => {}, (14,1) => {({17,13,42,36},1/1),({17,13,40,38},1/1),({16,14,43,35},1/1),({16,14,41,37},1/1),({16,14,39,39},1/1),({15,15,44,34},1/1),({15,15,42,36},1/1),({15,15,40,38},1/1)}, (13,2) => {}, (15,1) => {}, (14,2) => {}, (15,2) => {({17,17,44,44},1/1)}, (16,2) => {}, (0,0) => {({0,0,3,0},1/1)}, (1,0) => {({2,0,7,1},1/1),({2,0,6,2},1/1),({2,0,5,3},1/1)}, (1,1) => {({2,2,13,0},1/1)}, (2,0) => {({4,0,10,3},1/1),({4,0,9,4},1/1),({4,0,8,5},1/1),({4,0,7,6},1/1),({3,1,11,2},1/1),({3,1,10,3},1/1),({3,1,9,4},1/1),({3,1,8,5},1/1),({3,1,7,6},1/1),({2,2,10,3},1/1),({2,2,9,4},1/1),({2,2,8,5},1/1),({2,2,7,6},1/1)}, (2,1) => {({4,2,17,1},1/1),({4,2,16,2},1/1),({4,2,15,3},1/1),({4,2,14,4},1/1),({4,2,13,5},1/1)}, (3,0) => {({6,0,12,6},1/1),({6,0,10,8},1/1),({5,1,14,4},1/1),({5,1,13,5},1/1),({5,1,12,6},1/1),({5,1,11,7},2/1),({5,1,10,8},1/1),({4,2,14,4},1/1),({4,2,13,5},1/1),({4,2,12,6},2/1),({4,2,11,7},2/1),({4,2,10,8},2/1),({3,3,15,3},1/1),({3,3,13,5},1/1),({3,3,12,6},1/1),({3,3,11,7},1/1),({3,3,9,9},1/1)}, (3,1) => {({6,2,20,3},1/1),({6,2,19,4},1/1),({6,2,18,5},2/1),({6,2,17,6},2/1),({6,2,16,7},2/1),({6,2,15,8},1/1),({6,2,14,9},1/1),({5,3,21,2},1/1),({5,3,20,3},1/1),({5,3,19,4},2/1),({5,3,18,5},2/1),({5,3,17,6},3/1),({5,3,16,7},2/1),({5,3,15,8},2/1),({5,3,14,9},1/1),({5,3,13,10},1/1),({4,4,20,3},1/1),({4,4,19,4},1/1),({4,4,18,5},2/1),({4,4,17,6},2/1),({4,4,16,7},2/1),({4,4,15,8},1/1),({4,4,14,9},1/1)}, (4,0) => {}, (2,2) => {}, (4,1) => {({10,0,14,14},1/1),({9,1,18,10},1/1),({9,1,17,11},1/1),({9,1,16,12},1/1),({9,1,15,13},1/1),({8,2,22,6},1/1),({8,2,21,7},1/1),({8,2,20,8},3/1),({8,2,19,9},3/1),({8,2,18,10},5/1),({8,2,17,11},4/1),({8,2,16,12},5/1),({8,2,15,13},2/1),({8,2,14,14},2/1),({7,3,24,4},1/1),({7,3,23,5},2/1),({7,3,22,6},3/1),({7,3,21,7},6/1),({7,3,20,8},8/1),({7,3,19,9},10/1),({7,3,18,10},11/1),({7,3,17,11},12/1),({7,3,16,12},9/1),({7,3,15,13},7/1),({7,3,14,14},2/1),({6,4,24,4},1/1),({6,4,23,5},2/1),({6,4,22,6},5/1),({6,4,21,7},7/1),({6,4,20,8},11/1),({6,4,19,9},12/1),({6,4,18,10},15/1),({6,4,17,11},12/1),({6,4,16,12},12/1),({6,4,15,13},6/1),({6,4,14,14},4/1),({5,5,25,3},1/1),({5,5,24,4},1/1),({5,5,23,5},2/1),({5,5,22,6},3/1),({5,5,21,7},5/1),({5,5,20,8},5/1),({5,5,19,9},8/1),({5,5,18,10},6/1),({5,5,17,11},8/1),({5,5,16,12},5/1),({5,5,15,13},5/1)}, (5,0) => {}, (3,2) => {}, (5,1) => {({11,1,19,14},1/1),({11,1,18,15},1/1),({10,2,23,10},1/1),({10,2,22,11},2/1),({10,2,21,12},3/1),({10,2,20,13},4/1),({10,2,19,14},5/1),({10,2,18,15},4/1),({10,2,17,16},2/1),({9,3,26,7},1/1),({9,3,25,8},2/1),({9,3,24,9},4/1),({9,3,23,10},7/1),({9,3,22,11},11/1),({9,3,21,12},14/1),({9,3,20,13},17/1),({9,3,19,14},16/1),({9,3,18,15},13/1),({9,3,17,16},8/1),({8,4,27,6},1/1),({8,4,26,7},2/1),({8,4,25,8},6/1),({8,4,24,9},11/1),({8,4,23,10},17/1),({8,4,22,11},24/1),({8,4,21,12},31/1),({8,4,20,13},32/1),({8,4,19,14},31/1),({8,4,18,15},25/1),({8,4,17,16},13/1),({7,5,28,5},1/1),({7,5,27,6},2/1),({7,5,26,7},5/1),({7,5,25,8},9/1),({7,5,24,9},16/1),({7,5,23,10},23/1),({7,5,22,11},31/1),({7,5,21,12},36/1),({7,5,20,13},39/1),({7,5,19,14},36/1),({7,5,18,15},28/1),({7,5,17,16},15/1),({6,6,27,6},1/1),({6,6,26,7},2/1),({6,6,25,8},4/1),({6,6,24,9},6/1),({6,6,23,10},11/1),({6,6,22,11},13/1),({6,6,21,12},16/1),({6,6,20,13},17/1),({6,6,19,14},16/1),({6,6,18,15},11/1),({6,6,17,16},7/1)}, (6,0) => {}, (4,2) => {}};
--dw stands for dominant weights
dw0325 = new HashTable from {(5,2) => {}, (7,0) => {}, (6,1) => {({12,2,23,15},2/1),({11,3,27,11},1/1),({10,4,29,9},1/1),({9,5,30,8},1/1),({8,6,31,7},1/1)}, (7,1) => {({13,3,27,16},1/1),({12,4,30,13},1/1),({11,5,31,12},2/1),({10,6,33,10},1/1),({8,8,34,9},1/1)}, (8,0) => {}, (6,2) => {}, (7,2) => {}, (9,0) => {}, (8,1) => {({14,4,30,18},1/1),({13,5,32,16},1/1),({12,6,34,14},1/1),({11,7,35,13},1/1),({10,8,36,12},1/1)}, (8,2) => {}, (9,1) => {({15,5,32,21},1/1),({14,6,34,19},1/1),({13,7,36,17},1/1),({12,8,37,16},1/1),({11,9,38,15},1/1)}, (10,0) => {}, (11,0) => {}, (10,1) => {({16,6,33,25},1/1),({15,7,36,22},1/1),({14,8,37,21},3/1),({13,9,39,19},1/1),({11,11,40,18},1/1)}, (9,2) => {}, (12,0) => {}, (11,1) => {({17,7,33,30},1/1),({16,8,37,26},1/1),({15,9,39,24},1/1),({14,10,40,23},2/1),({13,11,41,22},1/1)}, (10,2) => {}, (13,0) => {}, (12,1) => {({17,9,37,31},1/1),({16,10,40,28},1/1),({15,11,41,27},2/1),({14,12,42,26},2/1)}, (11,2) => {}, (14,0) => {}, (13,1) => {({17,11,40,33},1/1),({16,12,42,31},1/1),({15,13,43,30},1/1)}, (12,2) => {}, (13,2) => {}, (14,1) => {({17,13,42,36},1/1),({16,14,43,35},1/1),({15,15,44,34},1/1)}, (15,0) => {}, (14,2) => {}, (15,1) => {}, (15,2) => {({17,17,44,44},1/1)}, (16,2) => {}, (0,0) => {({0,0,3,0},1/1)}, (1,0) => {({2,0,7,1},1/1)}, (2,0) => {({4,0,10,3},1/1),({3,1,11,2},1/1)}, (1,1) => {({2,2,13,0},1/1)}, (3,0) => {({6,0,12,6},1/1),({5,1,14,4},1/1),({3,3,15,3},1/1)}, (2,1) => {({4,2,17,1},1/1)}, (2,2) => {}, (4,0) => {}, (3,1) => {({6,2,20,3},1/1),({5,3,21,2},1/1)}, (3,2) => {}, (5,0) => {}, (4,1) => {({10,0,14,14},1/1),({9,1,18,10},1/1),({8,2,22,6},1/1),({7,3,24,4},1/1),({5,5,25,3},1/1)}, (4,2) => {}, (6,0) => {}, (5,1) => {({11,1,19,14},1/1),({10,2,23,10},1/1),({9,3,26,7},1/1),({8,4,27,6},1/1),({7,5,28,5},1/1)}};
--dmw stands for dominant monomial weights
dmw0325 = new HashTable from {(5,2) => {}, (7,0) => {}, (6,1) => {{12,2,23,15},{11,3,24,14}}, (7,1) => {{12,4,28,15},{13,3,26,17}}, (8,0) => {}, (6,2) => {}, (7,2) => {}, (9,0) => {}, (8,1) => {{13,5,31,17},{14,4,28,20},{12,6,32,16}}, (8,2) => {}, (9,1) => {{15,5,29,24},{14,6,33,20},{13,7,35,18}}, (10,0) => {}, (11,0) => {}, (10,1) => {{15,7,34,24},{13,9,38,20},{14,8,37,21}}, (9,2) => {{16,6,24,34}}, (12,0) => {}, (11,1) => {{14,10,40,23},{15,9,38,25}}, (10,2) => {{17,7,29,34}}, (13,0) => {}, (12,1) => {{14,12,42,26},{15,11,41,27}}, (11,2) => {{17,9,34,34}}, (14,0) => {}, (13,1) => {{15,13,43,30}}, (12,2) => {{17,11,38,35}}, (13,2) => {{17,13,41,37}}, (14,1) => {{15,15,44,34}}, (15,0) => {}, (14,2) => {{17,15,43,40}}, (15,1) => {}, (15,2) => {{17,17,44,44}}, (16,2) => {}, (0,0) => {}, (1,0) => {}, (2,0) => {}, (1,1) => {}, (3,0) => {}, (2,1) => {}, (2,2) => {}, (4,0) => {}, (3,1) => {}, (3,2) => {}, (5,0) => {}, (4,1) => {{10,0,14,14}}, (4,2) => {}, (6,0) => {}, (5,1) => {{11,1,19,14}}};
end;