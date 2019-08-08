A := QQ[t_0,t_1,t_2,t_3];
--tb stands for Total Betti numbers
tb0025 = new HashTable from {(6,1) => 56628/1, (7,0) => 0, (5,2) => 0, (7,1) => 70070/1, (8,0) => 0, (6,2) => 0, (8,1) => 65780/1, (9,0) => 0, (7,2) => 0, (10,0) => 0, (9,1) => 46332/1, (8,2) => 0, (11,0) => 0, (10,1) => 23660/1, (9,2) => 0, (12,0) => 0, (11,1) => 8008/1, (10,2) => 0, (13,0) => 0, (12,1) => 1260/1, (11,2) => 0, (14,0) => 0, (13,1) => 195/1, (12,2) => 455/1, (15,0) => 0, (14,1) => 14/1, (13,2) => 210/1, (15,1) => 0, (14,2) => 45/1, (15,2) => 4/1, (16,2) => 0, (0,0) => 1/1, (1,0) => 0, (1,1) => 116/1, (2,0) => 0, (2,1) => 1060/1, (3,0) => 0, (3,1) => 5040/1, (4,0) => 0, (2,2) => 0, (4,1) => 15652/1, (5,0) => 0, (3,2) => 0, (5,1) => 34580/1, (6,0) => 0, (4,2) => 0};
end;--sb represents the betti numbers as sums of Schur functors
sb0025 = new HashTable from {(6,1) => {({12,2,23,12},1/1),({12,2,22,13},1/1),({12,2,21,14},2/1),({12,2,20,15},3/1),({12,2,19,16},2/1),({12,2,18,17},1/1),({11,3,25,10},1/1),({11,3,24,11},4/1),({11,3,23,12},7/1),({11,3,22,13},11/1),({11,3,21,14},14/1),({11,3,20,15},15/1),({11,3,19,16},13/1),({11,3,18,17},7/1),({10,4,27,8},1/1),({10,4,26,9},4/1),({10,4,25,10},9/1),({10,4,24,11},17/1),({10,4,23,12},28/1),({10,4,22,13},38/1),({10,4,21,14},45/1),({10,4,20,15},45/1),({10,4,19,16},37/1),({10,4,18,17},21/1),({9,5,28,7},1/1),({9,5,27,8},4/1),({9,5,26,9},11/1),({9,5,25,10},22/1),({9,5,24,11},38/1),({9,5,23,12},56/1),({9,5,22,13},73/1),({9,5,21,14},82/1),({9,5,20,15},80/1),({9,5,19,16},64/1),({9,5,18,17},36/1),({8,6,29,6},1/1),({8,6,28,7},3/1),({8,6,27,8},8/1),({8,6,26,9},16/1),({8,6,25,10},30/1),({8,6,24,11},47/1),({8,6,23,12},66/1),({8,6,22,13},82/1),({8,6,21,14},91/1),({8,6,20,15},87/1),({8,6,19,16},69/1),({8,6,18,17},38/1),({7,7,28,7},1/1),({7,7,27,8},3/1),({7,7,26,9},7/1),({7,7,25,10},13/1),({7,7,24,11},21/1),({7,7,23,12},29/1),({7,7,22,13},36/1),({7,7,21,14},39/1),({7,7,20,15},38/1),({7,7,19,16},29/1),({7,7,18,17},16/1)}, (7,0) => {}, (5,2) => {}, (6,2) => {}, (8,0) => {}, (7,1) => {({13,3,26,14},1/1),({13,3,25,15},1/1),({13,3,24,16},3/1),({13,3,23,17},3/1),({13,3,22,18},3/1),({13,3,21,19},2/1),({13,3,20,20},2/1),({12,4,28,12},1/1),({12,4,27,13},4/1),({12,4,26,14},7/1),({12,4,25,15},12/1),({12,4,24,16},16/1),({12,4,23,17},19/1),({12,4,22,18},17/1),({12,4,21,19},13/1),({12,4,20,20},4/1),({11,5,30,10},1/1),({11,5,29,11},3/1),({11,5,28,12},9/1),({11,5,27,13},16/1),({11,5,26,14},29/1),({11,5,25,15},40/1),({11,5,24,16},52/1),({11,5,23,17},53/1),({11,5,22,18},51/1),({11,5,21,19},33/1),({11,5,20,20},14/1),({10,6,31,9},1/1),({10,6,30,10},4/1),({10,6,29,11},11/1),({10,6,28,12},21/1),({10,6,27,13},39/1),({10,6,26,14},58/1),({10,6,25,15},79/1),({10,6,24,16},92/1),({10,6,23,17},99/1),({10,6,22,18},84/1),({10,6,21,19},61/1),({10,6,20,20},20/1),({9,7,32,8},1/1),({9,7,31,9},2/1),({9,7,30,10},7/1),({9,7,29,11},14/1),({9,7,28,12},29/1),({9,7,27,13},45/1),({9,7,26,14},69/1),({9,7,25,15},86/1),({9,7,24,16},104/1),({9,7,23,17},101/1),({9,7,22,18},93/1),({9,7,21,19},59/1),({9,7,20,20},25/1),({8,8,31,9},2/1),({8,8,30,10},3/1),({8,8,29,11},8/1),({8,8,28,12},13/1),({8,8,27,13},23/1),({8,8,26,14},29/1),({8,8,25,15},42/1),({8,8,24,16},43/1),({8,8,23,17},48/1),({8,8,22,18},38/1),({8,8,21,19},30/1),({8,8,20,20},6/1)}, (8,1) => {({14,4,28,17},1/1),({14,4,27,18},2/1),({14,4,26,19},2/1),({14,4,25,20},3/1),({14,4,24,21},3/1),({14,4,23,22},1/1),({13,5,31,14},1/1),({13,5,30,15},2/1),({13,5,29,16},5/1),({13,5,28,17},9/1),({13,5,27,18},13/1),({13,5,26,19},16/1),({13,5,25,20},17/1),({13,5,24,21},14/1),({13,5,23,22},8/1),({12,6,32,13},2/1),({12,6,31,14},5/1),({12,6,30,15},12/1),({12,6,29,16},21/1),({12,6,28,17},32/1),({12,6,27,18},42/1),({12,6,26,19},49/1),({12,6,25,20},48/1),({12,6,24,21},39/1),({12,6,23,22},22/1),({11,7,34,11},1/1),({11,7,33,12},3/1),({11,7,32,13},7/1),({11,7,31,14},16/1),({11,7,30,15},29/1),({11,7,29,16},46/1),({11,7,28,17},65/1),({11,7,27,18},82/1),({11,7,26,19},89/1),({11,7,25,20},86/1),({11,7,24,21},68/1),({11,7,23,22},37/1),({10,8,34,11},1/1),({10,8,33,12},4/1),({10,8,32,13},10/1),({10,8,31,14},20/1),({10,8,30,15},35/1),({10,8,29,16},53/1),({10,8,28,17},73/1),({10,8,27,18},88/1),({10,8,26,19},96/1),({10,8,25,20},90/1),({10,8,24,21},71/1),({10,8,23,22},39/1),({9,9,35,10},1/1),({9,9,34,11},1/1),({9,9,33,12},3/1),({9,9,32,13},6/1),({9,9,31,14},11/1),({9,9,30,15},17/1),({9,9,29,16},27/1),({9,9,28,17},34/1),({9,9,27,18},41/1),({9,9,26,19},44/1),({9,9,25,20},41/1),({9,9,24,21},31/1),({9,9,23,22},18/1)}, (9,0) => {}, (7,2) => {}, (10,0) => {}, (9,1) => {({15,5,29,21},2/1),({15,5,28,22},1/1),({15,5,27,23},2/1),({15,5,26,24},1/1),({15,5,25,25},1/1),({14,6,33,17},1/1),({14,6,32,18},2/1),({14,6,31,19},4/1),({14,6,30,20},7/1),({14,6,29,21},9/1),({14,6,28,22},11/1),({14,6,27,23},10/1),({14,6,26,24},8/1),({14,6,25,25},3/1),({13,7,35,15},1/1),({13,7,34,16},2/1),({13,7,33,17},6/1),({13,7,32,18},11/1),({13,7,31,19},19/1),({13,7,30,20},25/1),({13,7,29,21},33/1),({13,7,28,22},33/1),({13,7,27,23},32/1),({13,7,26,24},20/1),({13,7,25,25},9/1),({12,8,36,14},1/1),({12,8,35,15},3/1),({12,8,34,16},8/1),({12,8,33,17},15/1),({12,8,32,18},27/1),({12,8,31,19},39/1),({12,8,30,20},53/1),({12,8,29,21},60/1),({12,8,28,22},64/1),({12,8,27,23},54/1),({12,8,26,24},39/1),({12,8,25,25},12/1),({11,9,37,13},1/1),({11,9,36,14},2/1),({11,9,35,15},6/1),({11,9,34,16},11/1),({11,9,33,17},22/1),({11,9,32,18},32/1),({11,9,31,19},49/1),({11,9,30,20},59/1),({11,9,29,21},71/1),({11,9,28,22},67/1),({11,9,27,23},62/1),({11,9,26,24},38/1),({11,9,25,25},17/1),({10,10,36,14},1/1),({10,10,35,15},2/1),({10,10,34,16},6/1),({10,10,33,17},9/1),({10,10,32,18},16/1),({10,10,31,19},20/1),({10,10,30,20},28/1),({10,10,29,21},29/1),({10,10,28,22},32/1),({10,10,27,23},24/1),({10,10,26,24},20/1),({10,10,25,25},4/1)}, (8,2) => {}, (9,2) => {}, (10,1) => {({16,6,30,25},1/1),({16,6,29,26},1/1),({15,7,34,21},1/1),({15,7,33,22},2/1),({15,7,32,23},3/1),({15,7,31,24},4/1),({15,7,30,25},5/1),({15,7,29,26},4/1),({15,7,28,27},2/1),({14,8,37,18},1/1),({14,8,36,19},2/1),({14,8,35,20},4/1),({14,8,34,21},7/1),({14,8,33,22},11/1),({14,8,32,23},14/1),({14,8,31,24},17/1),({14,8,30,25},16/1),({14,8,29,26},13/1),({14,8,28,27},8/1),({13,9,38,17},1/1),({13,9,37,18},2/1),({13,9,36,19},6/1),({13,9,35,20},11/1),({13,9,34,21},17/1),({13,9,33,22},24/1),({13,9,32,23},31/1),({13,9,31,24},32/1),({13,9,30,25},31/1),({13,9,29,26},25/1),({13,9,28,27},13/1),({12,10,39,16},1/1),({12,10,38,17},2/1),({12,10,37,18},5/1),({12,10,36,19},9/1),({12,10,35,20},16/1),({12,10,34,21},23/1),({12,10,33,22},31/1),({12,10,32,23},36/1),({12,10,31,24},39/1),({12,10,30,25},36/1),({12,10,29,26},28/1),({12,10,28,27},15/1),({11,11,38,17},1/1),({11,11,37,18},2/1),({11,11,36,19},4/1),({11,11,35,20},6/1),({11,11,34,21},11/1),({11,11,33,22},13/1),({11,11,32,23},16/1),({11,11,31,24},17/1),({11,11,30,25},16/1),({11,11,29,26},11/1),({11,11,28,27},7/1)}, (11,0) => {}, (10,2) => {}, (11,1) => {({17,7,30,30},1/1),({16,8,34,26},1/1),({16,8,33,27},1/1),({16,8,32,28},1/1),({16,8,31,29},1/1),({15,9,38,22},1/1),({15,9,37,23},1/1),({15,9,36,24},3/1),({15,9,35,25},3/1),({15,9,34,26},5/1),({15,9,33,27},4/1),({15,9,32,28},5/1),({15,9,31,29},2/1),({15,9,30,30},2/1),({14,10,40,20},1/1),({14,10,39,21},2/1),({14,10,38,22},3/1),({14,10,37,23},6/1),({14,10,36,24},8/1),({14,10,35,25},10/1),({14,10,34,26},11/1),({14,10,33,27},12/1),({14,10,32,28},9/1),({14,10,31,29},7/1),({14,10,30,30},2/1),({13,11,40,20},1/1),({13,11,39,21},2/1),({13,11,38,22},5/1),({13,11,37,23},7/1),({13,11,36,24},11/1),({13,11,35,25},12/1),({13,11,34,26},15/1),({13,11,33,27},12/1),({13,11,32,28},12/1),({13,11,31,29},6/1),({13,11,30,30},4/1),({12,12,41,19},1/1),({12,12,40,20},1/1),({12,12,39,21},2/1),({12,12,38,22},3/1),({12,12,37,23},5/1),({12,12,36,24},5/1),({12,12,35,25},8/1),({12,12,34,26},6/1),({12,12,33,27},8/1),({12,12,32,28},5/1),({12,12,31,29},5/1)}, (12,0) => {}, (11,2) => {}, (12,1) => {({15,11,41,24},1/1),({15,11,40,25},1/1),({15,11,39,26},2/1),({15,11,38,27},2/1),({15,11,37,28},2/1),({15,11,36,29},1/1),({15,11,35,30},1/1),({14,12,42,23},1/1),({14,12,41,24},1/1),({14,12,40,25},2/1),({14,12,39,26},2/1),({14,12,38,27},3/1),({14,12,37,28},2/1),({14,12,36,29},2/1),({14,12,35,30},1/1),({14,12,34,31},1/1),({13,13,41,24},1/1),({13,13,40,25},1/1),({13,13,39,26},2/1),({13,13,38,27},2/1),({13,13,37,28},2/1),({13,13,36,29},1/1),({13,13,35,30},1/1)}, (13,0) => {}, (12,2) => {({17,11,38,32},1/1),({17,11,36,34},1/1),({16,12,40,30},1/1),({16,12,39,31},1/1),({16,12,38,32},1/1),({16,12,37,33},2/1),({16,12,36,34},1/1),({15,13,40,30},1/1),({15,13,39,31},1/1),({15,13,38,32},2/1),({15,13,37,33},2/1),({15,13,36,34},2/1),({14,14,41,29},1/1),({14,14,39,31},1/1),({14,14,38,32},1/1),({14,14,37,33},1/1),({14,14,35,35},1/1)}, (13,1) => {({15,13,43,27},1/1),({15,13,42,28},1/1),({15,13,41,29},1/1),({15,13,40,30},1/1),({15,13,39,31},1/1)}, (14,0) => {}, (15,0) => {}, (14,1) => {({15,15,44,31},1/1)}, (13,2) => {({17,13,41,34},1/1),({17,13,40,35},1/1),({17,13,39,36},1/1),({17,13,38,37},1/1),({16,14,42,33},1/1),({16,14,41,34},1/1),({16,14,40,35},1/1),({16,14,39,36},1/1),({16,14,38,37},1/1),({15,15,41,34},1/1),({15,15,40,35},1/1),({15,15,39,36},1/1),({15,15,38,37},1/1)}, (15,1) => {}, (14,2) => {({17,15,43,37},1/1),({17,15,42,38},1/1),({17,15,41,39},1/1)}, (15,2) => {({17,17,44,41},1/1)}, (16,2) => {}, (0,0) => {({0,0,0,0},1/1)}, (1,0) => {}, (1,1) => {({4,0,8,2},1/1),({4,0,6,4},1/1),({3,1,9,1},1/1),({3,1,7,3},1/1),({3,1,5,5},1/1),({2,2,10,0},1/1),({2,2,8,2},1/1),({2,2,6,4},1/1)}, (2,0) => {}, (2,1) => {({6,0,11,4},1/1),({6,0,10,5},1/1),({6,0,9,6},1/1),({6,0,8,7},1/1),({5,1,13,2},1/1),({5,1,12,3},2/1),({5,1,11,4},2/1),({5,1,10,5},3/1),({5,1,9,6},3/1),({5,1,8,7},1/1),({4,2,14,1},1/1),({4,2,13,2},2/1),({4,2,12,3},3/1),({4,2,11,4},4/1),({4,2,10,5},5/1),({4,2,9,6},4/1),({4,2,8,7},2/1),({3,3,14,1},1/1),({3,3,13,2},1/1),({3,3,12,3},1/1),({3,3,11,4},2/1),({3,3,10,5},2/1),({3,3,9,6},1/1),({3,3,8,7},1/1)}, (3,0) => {}, (3,1) => {({8,0,13,7},1/1),({8,0,12,8},1/1),({8,0,11,9},1/1),({7,1,16,4},1/1),({7,1,15,5},2/1),({7,1,14,6},4/1),({7,1,13,7},4/1),({7,1,12,8},5/1),({7,1,11,9},3/1),({7,1,10,10},2/1),({6,2,17,3},2/1),({6,2,16,4},3/1),({6,2,15,5},7/1),({6,2,14,6},9/1),({6,2,13,7},12/1),({6,2,12,8},10/1),({6,2,11,9},9/1),({6,2,10,10},2/1),({5,3,18,2},2/1),({5,3,17,3},3/1),({5,3,16,4},7/1),({5,3,15,5},9/1),({5,3,14,6},14/1),({5,3,13,7},13/1),({5,3,12,8},14/1),({5,3,11,9},8/1),({5,3,10,10},5/1),({4,4,17,3},2/1),({4,4,16,4},2/1),({4,4,15,5},5/1),({4,4,14,6},5/1),({4,4,13,7},8/1),({4,4,12,8},5/1),({4,4,11,9},6/1)}, (4,0) => {}, (2,2) => {}, (4,1) => {({10,0,14,11},1/1),({9,1,18,7},1/1),({9,1,17,8},2/1),({9,1,16,9},3/1),({9,1,15,10},4/1),({9,1,14,11},3/1),({9,1,13,12},2/1),({8,2,20,5},1/1),({8,2,19,6},3/1),({8,2,18,7},6/1),({8,2,17,8},10/1),({8,2,16,9},13/1),({8,2,15,10},14/1),({8,2,14,11},12/1),({8,2,13,12},7/1),({7,3,21,4},2/1),({7,3,20,5},5/1),({7,3,19,6},10/1),({7,3,18,7},17/1),({7,3,17,8},24/1),({7,3,16,9},28/1),({7,3,15,10},29/1),({7,3,14,11},24/1),({7,3,13,12},13/1),({6,4,22,3},1/1),({6,4,21,4},3/1),({6,4,20,5},7/1),({6,4,19,6},13/1),({6,4,18,7},21/1),({6,4,17,8},28/1),({6,4,16,9},33/1),({6,4,15,10},33/1),({6,4,14,11},27/1),({6,4,13,12},15/1),({5,5,22,3},1/1),({5,5,21,4},2/1),({5,5,20,5},4/1),({5,5,19,6},8/1),({5,5,18,7},11/1),({5,5,17,8},14/1),({5,5,16,9},16/1),({5,5,15,10},15/1),({5,5,14,11},12/1),({5,5,13,12},7/1)}, (5,0) => {}, (3,2) => {}, (5,1) => {({11,1,19,11},1/1),({11,1,18,12},1/1),({11,1,17,13},1/1),({11,1,16,14},1/1),({11,1,15,15},1/1),({10,2,22,8},1/1),({10,2,21,9},2/1),({10,2,20,10},5/1),({10,2,19,11},7/1),({10,2,18,12},9/1),({10,2,17,13},8/1),({10,2,16,14},7/1),({10,2,15,15},2/1),({9,3,23,7},3/1),({9,3,22,8},6/1),({9,3,21,9},13/1),({9,3,20,10},19/1),({9,3,19,11},27/1),({9,3,18,12},28/1),({9,3,17,13},29/1),({9,3,16,14},18/1),({9,3,15,15},8/1),({8,4,25,5},1/1),({8,4,24,6},4/1),({8,4,23,7},8/1),({8,4,22,8},18/1),({8,4,21,9},29/1),({8,4,20,10},42/1),({8,4,19,11},51/1),({8,4,18,12},58/1),({8,4,17,13},49/1),({8,4,16,14},37/1),({8,4,15,15},12/1),({7,5,25,5},2/1),({7,5,24,6},5/1),({7,5,23,7},13/1),({7,5,22,8},22/1),({7,5,21,9},37/1),({7,5,20,10},48/1),({7,5,19,11},61/1),({7,5,18,12},60/1),({7,5,17,13},57/1),({7,5,16,14},36/1),({7,5,15,15},16/1),({6,6,26,4},1/1),({6,6,25,5},1/1),({6,6,24,6},4/1),({6,6,23,7},6/1),({6,6,22,8},12/1),({6,6,21,9},15/1),({6,6,20,10},25/1),({6,6,19,11},25/1),({6,6,18,12},30/1),({6,6,17,13},23/1),({6,6,16,14},20/1),({6,6,15,15},3/1)}, (6,0) => {}, (4,2) => {}};
end;