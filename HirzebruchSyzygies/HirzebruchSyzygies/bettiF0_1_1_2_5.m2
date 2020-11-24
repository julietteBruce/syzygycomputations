A := QQ[t_0,t_1,t_2,t_3];
--tb stands for Total Betti numbers
tb1125 = new HashTable from {(5,2) => 0, (7,0) => 0, (6,1) => 54340/1, (7,1) => 77220/1, (8,0) => 0, (6,2) => 0, (7,2) => 0, (9,0) => 0, (8,1) => 82940/1, (8,2) => 0, (9,1) => 68068/1, (10,0) => 0, (11,0) => 0, (10,1) => 42588/1, (9,2) => 0, (12,0) => 0, (11,1) => 20020/1, (10,2) => 0, (13,0) => 0, (12,1) => 6860/1, (11,2) => 0, (14,0) => 0, (13,1) => 1620/1, (12,2) => 0, (13,2) => 0, (14,1) => 236/1, (15,0) => 0, (14,2) => 0, (15,1) => 16/1, (16,1) => 0, (15,2) => 0, (17,1) => 0, (0,0) => 4/1, (1,0) => 44/1, (2,0) => 180/1, (1,1) => 0, (3,0) => 220/1, (2,1) => 80/1, (2,2) => 0, (4,0) => 0, (3,1) => 1820/1, (3,2) => 0, (5,0) => 0, (4,1) => 9828/1, (4,2) => 0, (6,0) => 0, (5,1) => 28028/1};
--sb represents the betti numbers as sums of Schur functors
sb1125 = new HashTable from {(5,2) => {}, (7,0) => {}, (6,1) => {({13,2,23,13},1/1),({13,2,22,14},1/1),({13,2,21,15},2/1),({13,2,20,16},1/1),({13,2,19,17},2/1),({12,3,25,11},1/1),({12,3,24,12},3/1),({12,3,23,13},5/1),({12,3,22,14},8/1),({12,3,21,15},10/1),({12,3,20,16},9/1),({12,3,19,17},7/1),({12,3,18,18},3/1),({11,4,27,9},1/1),({11,4,26,10},3/1),({11,4,25,11},8/1),({11,4,24,12},13/1),({11,4,23,13},22/1),({11,4,22,14},27/1),({11,4,21,15},32/1),({11,4,20,16},28/1),({11,4,19,17},22/1),({11,4,18,18},6/1),({10,5,28,8},1/1),({10,5,27,9},4/1),({10,5,26,10},10/1),({10,5,25,11},20/1),({10,5,24,12},33/1),({10,5,23,13},47/1),({10,5,22,14},58/1),({10,5,21,15},62/1),({10,5,20,16},56/1),({10,5,19,17},39/1),({10,5,18,18},14/1),({9,6,29,7},1/1),({9,6,28,8},3/1),({9,6,27,9},9/1),({9,6,26,10},17/1),({9,6,25,11},32/1),({9,6,24,12},48/1),({9,6,23,13},66/1),({9,6,22,14},77/1),({9,6,21,15},83/1),({9,6,20,16},70/1),({9,6,19,17},51/1),({9,6,18,18},17/1),({8,7,29,7},1/1),({8,7,28,8},3/1),({8,7,27,9},7/1),({8,7,26,10},15/1),({8,7,25,11},25/1),({8,7,24,12},37/1),({8,7,23,13},50/1),({8,7,22,14},58/1),({8,7,21,15},59/1),({8,7,20,16},53/1),({8,7,19,17},36/1),({8,7,18,18},12/1)}, (7,1) => {({14,3,26,15},1/1),({14,3,25,16},1/1),({14,3,24,17},2/1),({14,3,23,18},2/1),({14,3,22,19},2/1),({14,3,21,20},1/1),({13,4,28,13},1/1),({13,4,27,14},3/1),({13,4,26,15},6/1),({13,4,25,16},9/1),({13,4,24,17},12/1),({13,4,23,18},13/1),({13,4,22,19},11/1),({13,4,21,20},6/1),({12,5,30,11},1/1),({12,5,29,12},3/1),({12,5,28,13},8/1),({12,5,27,14},15/1),({12,5,26,15},25/1),({12,5,25,16},34/1),({12,5,24,17},41/1),({12,5,23,18},41/1),({12,5,22,19},34/1),({12,5,21,20},19/1),({11,6,31,10},1/1),({11,6,30,11},4/1),({11,6,29,12},11/1),({11,6,28,13},22/1),({11,6,27,14},38/1),({11,6,26,15},57/1),({11,6,25,16},74/1),({11,6,24,17},84/1),({11,6,23,18},82/1),({11,6,22,19},66/1),({11,6,21,20},37/1),({10,7,32,9},1/1),({10,7,31,10},3/1),({10,7,30,11},9/1),({10,7,29,12},19/1),({10,7,28,13},36/1),({10,7,27,14},57/1),({10,7,26,15},81/1),({10,7,25,16},101/1),({10,7,24,17},112/1),({10,7,23,18},107/1),({10,7,22,19},85/1),({10,7,21,20},47/1),({9,8,32,9},1/1),({9,8,31,10},3/1),({9,8,30,11},8/1),({9,8,29,12},16/1),({9,8,28,13},29/1),({9,8,27,14},45/1),({9,8,26,15},62/1),({9,8,25,16},76/1),({9,8,24,17},83/1),({9,8,23,18},79/1),({9,8,22,19},62/1),({9,8,21,20},34/1)}, (8,0) => {}, (6,2) => {}, (7,2) => {}, (9,0) => {}, (8,1) => {({15,4,28,18},1/1),({15,4,27,19},1/1),({15,4,26,20},2/1),({15,4,25,21},1/1),({15,4,24,22},2/1),({14,5,31,15},1/1),({14,5,30,16},2/1),({14,5,29,17},4/1),({14,5,28,18},8/1),({14,5,27,19},11/1),({14,5,26,20},12/1),({14,5,25,21},12/1),({14,5,24,22},9/1),({14,5,23,23},3/1),({13,6,32,14},2/1),({13,6,31,15},5/1),({13,6,30,16},12/1),({13,6,29,17},20/1),({13,6,28,18},31/1),({13,6,27,19},38/1),({13,6,26,20},44/1),({13,6,25,21},38/1),({13,6,24,22},29/1),({13,6,23,23},9/1),({12,7,34,12},1/1),({12,7,33,13},3/1),({12,7,32,14},8/1),({12,7,31,15},18/1),({12,7,30,16},33/1),({12,7,29,17},52/1),({12,7,28,18},71/1),({12,7,27,19},86/1),({12,7,26,20},91/1),({12,7,25,21},81/1),({12,7,24,22},56/1),({12,7,23,23},20/1),({11,8,34,12},2/1),({11,8,33,13},6/1),({11,8,32,14},15/1),({11,8,31,15},29/1),({11,8,30,16},51/1),({11,8,29,17},74/1),({11,8,28,18},101/1),({11,8,27,19},116/1),({11,8,26,20},121/1),({11,8,25,21},105/1),({11,8,24,22},74/1),({11,8,23,23},24/1),({10,9,35,11},1/1),({10,9,34,12},2/1),({10,9,33,13},6/1),({10,9,32,14},14/1),({10,9,31,15},25/1),({10,9,30,16},41/1),({10,9,29,17},60/1),({10,9,28,18},77/1),({10,9,27,19},89/1),({10,9,26,20},91/1),({10,9,25,21},78/1),({10,9,24,22},54/1),({10,9,23,23},20/1)}, (8,2) => {}, (9,1) => {({16,5,29,22},1/1),({16,5,28,23},1/1),({16,5,27,24},1/1),({15,6,33,18},1/1),({15,6,32,19},2/1),({15,6,31,20},4/1),({15,6,30,21},6/1),({15,6,29,22},8/1),({15,6,28,23},9/1),({15,6,27,24},8/1),({15,6,26,25},4/1),({14,7,35,16},1/1),({14,7,34,17},2/1),({14,7,33,18},6/1),({14,7,32,19},12/1),({14,7,31,20},20/1),({14,7,30,21},27/1),({14,7,29,22},33/1),({14,7,28,23},33/1),({14,7,27,24},27/1),({14,7,26,25},16/1),({13,8,36,15},1/1),({13,8,35,16},4/1),({13,8,34,17},10/1),({13,8,33,18},20/1),({13,8,32,19},34/1),({13,8,31,20},51/1),({13,8,30,21},66/1),({13,8,29,22},75/1),({13,8,28,23},73/1),({13,8,27,24},59/1),({13,8,26,25},33/1),({12,9,37,14},1/1),({12,9,36,15},3/1),({12,9,35,16},8/1),({12,9,34,17},18/1),({12,9,33,18},34/1),({12,9,32,19},53/1),({12,9,31,20},76/1),({12,9,30,21},95/1),({12,9,29,22},104/1),({12,9,28,23},100/1),({12,9,27,24},80/1),({12,9,26,25},43/1),({11,10,37,14},1/1),({11,10,36,15},3/1),({11,10,35,16},8/1),({11,10,34,17},16/1),({11,10,33,18},28/1),({11,10,32,19},44/1),({11,10,31,20},60/1),({11,10,30,21},73/1),({11,10,29,22},81/1),({11,10,28,23},76/1),({11,10,27,24},59/1),({11,10,26,25},33/1)}, (10,0) => {}, (11,0) => {}, (10,1) => {({17,6,29,27},1/1),({16,7,34,22},1/1),({16,7,33,23},2/1),({16,7,32,24},3/1),({16,7,31,25},4/1),({16,7,30,26},4/1),({16,7,29,27},3/1),({16,7,28,28},1/1),({15,8,37,19},1/1),({15,8,36,20},2/1),({15,8,35,21},5/1),({15,8,34,22},8/1),({15,8,33,23},14/1),({15,8,32,24},17/1),({15,8,31,25},21/1),({15,8,30,26},18/1),({15,8,29,27},14/1),({15,8,28,28},4/1),({14,9,38,18},1/1),({14,9,37,19},3/1),({14,9,36,20},8/1),({14,9,35,21},16/1),({14,9,34,22},26/1),({14,9,33,23},37/1),({14,9,32,24},46/1),({14,9,31,25},49/1),({14,9,30,26},44/1),({14,9,29,27},31/1),({14,9,28,28},11/1),({13,10,39,17},1/1),({13,10,38,18},3/1),({13,10,37,19},8/1),({13,10,36,20},16/1),({13,10,35,21},29/1),({13,10,34,22},43/1),({13,10,33,23},59/1),({13,10,32,24},69/1),({13,10,31,25},73/1),({13,10,30,26},63/1),({13,10,29,27},45/1),({13,10,28,28},15/1),({12,11,39,17},1/1),({12,11,38,18},3/1),({12,11,37,19},7/1),({12,11,36,20},14/1),({12,11,35,21},24/1),({12,11,34,22},36/1),({12,11,33,23},47/1),({12,11,32,24},55/1),({12,11,31,25},56/1),({12,11,30,26},49/1),({12,11,29,27},34/1),({12,11,28,28},12/1)}, (9,2) => {}, (12,0) => {}, (11,1) => {({17,8,34,27},1/1),({17,8,33,28},1/1),({17,8,32,29},1/1),({17,8,31,30},1/1),({16,9,38,23},1/1),({16,9,37,24},2/1),({16,9,36,25},4/1),({16,9,35,26},6/1),({16,9,34,27},8/1),({16,9,33,28},8/1),({16,9,32,29},7/1),({16,9,31,30},4/1),({15,10,40,21},1/1),({15,10,39,22},2/1),({15,10,38,23},5/1),({15,10,37,24},9/1),({15,10,36,25},15/1),({15,10,35,26},20/1),({15,10,34,27},24/1),({15,10,33,28},24/1),({15,10,32,29},20/1),({15,10,31,30},11/1),({14,11,40,21},2/1),({14,11,39,22},5/1),({14,11,38,23},10/1),({14,11,37,24},18/1),({14,11,36,25},27/1),({14,11,35,26},34/1),({14,11,34,27},39/1),({14,11,33,28},38/1),({14,11,32,29},30/1),({14,11,31,30},17/1),({13,12,41,20},1/1),({13,12,40,21},2/1),({13,12,39,22},5/1),({13,12,38,23},10/1),({13,12,37,24},16/1),({13,12,36,25},23/1),({13,12,35,26},29/1),({13,12,34,27},32/1),({13,12,33,28},31/1),({13,12,32,29},25/1),({13,12,31,30},13/1)}, (10,2) => {}, (13,0) => {}, (12,1) => {({17,10,38,28},1/1),({17,10,37,29},1/1),({17,10,36,30},2/1),({17,10,35,31},2/1),({17,10,34,32},2/1),({16,11,41,25},1/1),({16,11,40,26},2/1),({16,11,39,27},4/1),({16,11,38,28},6/1),({16,11,37,29},8/1),({16,11,36,30},9/1),({16,11,35,31},8/1),({16,11,34,32},6/1),({16,11,33,33},2/1),({15,12,42,24},1/1),({15,12,41,25},2/1),({15,12,40,26},5/1),({15,12,39,27},8/1),({15,12,38,28},12/1),({15,12,37,29},15/1),({15,12,36,30},17/1),({15,12,35,31},14/1),({15,12,34,32},11/1),({15,12,33,33},4/1),({14,13,42,24},1/1),({14,13,41,25},3/1),({14,13,40,26},5/1),({14,13,39,27},8/1),({14,13,38,28},12/1),({14,13,37,29},14/1),({14,13,36,30},14/1),({14,13,35,31},13/1),({14,13,34,32},9/1),({14,13,33,33},3/1)}, (11,2) => {}, (14,0) => {}, (13,1) => {({17,12,41,30},1/1),({17,12,40,31},1/1),({17,12,39,32},2/1),({17,12,38,33},2/1),({17,12,37,34},2/1),({17,12,36,35},1/1),({16,13,43,28},1/1),({16,13,42,29},2/1),({16,13,41,30},3/1),({16,13,40,31},4/1),({16,13,39,32},5/1),({16,13,38,33},5/1),({16,13,37,34},4/1),({16,13,36,35},2/1),({15,14,43,28},1/1),({15,14,42,29},2/1),({15,14,41,30},3/1),({15,14,40,31},4/1),({15,14,39,32},5/1),({15,14,38,33},5/1),({15,14,37,34},4/1),({15,14,36,35},2/1)}, (12,2) => {}, (13,2) => {}, (14,1) => {({17,14,43,33},1/1),({17,14,42,34},1/1),({17,14,41,35},1/1),({17,14,40,36},1/1),({17,14,39,37},1/1),({16,15,44,32},1/1),({16,15,43,33},1/1),({16,15,42,34},1/1),({16,15,41,35},1/1),({16,15,40,36},1/1),({16,15,39,37},1/1)}, (15,0) => {}, (14,2) => {}, (15,1) => {({17,16,44,37},1/1)}, (16,1) => {}, (15,2) => {}, (17,1) => {}, (0,0) => {({1,0,1,0},1/1)}, (1,0) => {({3,0,5,1},1/1),({2,1,6,0},1/1),({2,1,5,1},1/1)}, (2,0) => {({4,1,10,1},1/1),({4,1,9,2},1/1),({4,1,8,3},1/1),({4,1,7,4},1/1),({4,1,6,5},1/1),({3,2,10,1},1/1),({3,2,9,2},1/1),({3,2,8,3},1/1),({3,2,7,4},1/1),({3,2,6,5},1/1)}, (1,1) => {}, (3,0) => {({5,2,14,2},1/1),({5,2,12,4},1/1),({5,2,11,5},1/1),({5,2,10,6},1/1),({5,2,8,8},1/1),({4,3,13,3},1/1),({4,3,12,4},1/1),({4,3,11,5},1/1),({4,3,10,6},2/1),({4,3,9,7},1/1)}, (2,1) => {({7,0,11,5},1/1),({7,0,9,7},1/1)}, (2,2) => {}, (4,0) => {}, (3,1) => {({9,0,13,8},1/1),({9,0,12,9},1/1),({8,1,16,5},1/1),({8,1,15,6},1/1),({8,1,14,7},2/1),({8,1,13,8},2/1),({8,1,12,9},2/1),({8,1,11,10},1/1),({7,2,16,5},1/1),({7,2,15,6},2/1),({7,2,14,7},3/1),({7,2,13,8},3/1),({7,2,12,9},3/1),({7,2,11,10},2/1),({6,3,18,3},1/1),({6,3,17,4},1/1),({6,3,16,5},2/1),({6,3,15,6},3/1),({6,3,14,7},3/1),({6,3,13,8},3/1),({6,3,12,9},3/1),({6,3,11,10},1/1),({5,4,17,4},1/1),({5,4,16,5},1/1),({5,4,15,6},1/1),({5,4,14,7},2/1),({5,4,13,8},2/1),({5,4,12,9},1/1),({5,4,11,10},1/1)}, (3,2) => {}, (5,0) => {}, (4,1) => {({11,0,14,12},1/1),({10,1,18,8},1/1),({10,1,17,9},2/1),({10,1,16,10},2/1),({10,1,15,11},2/1),({10,1,14,12},2/1),({10,1,13,13},1/1),({9,2,20,6},1/1),({9,2,19,7},1/1),({9,2,18,8},4/1),({9,2,17,9},5/1),({9,2,16,10},8/1),({9,2,15,11},6/1),({9,2,14,12},6/1),({9,2,13,13},1/1),({8,3,21,5},1/1),({8,3,20,6},3/1),({8,3,19,7},6/1),({8,3,18,8},9/1),({8,3,17,9},12/1),({8,3,16,10},14/1),({8,3,15,11},13/1),({8,3,14,12},9/1),({8,3,13,13},3/1),({7,4,22,4},1/1),({7,4,21,5},2/1),({7,4,20,6},5/1),({7,4,19,7},8/1),({7,4,18,8},13/1),({7,4,17,9},15/1),({7,4,16,10},18/1),({7,4,15,11},15/1),({7,4,14,12},12/1),({7,4,13,13},3/1),({6,5,22,4},1/1),({6,5,21,5},2/1),({6,5,20,6},4/1),({6,5,19,7},7/1),({6,5,18,8},10/1),({6,5,17,9},12/1),({6,5,16,10},12/1),({6,5,15,11},11/1),({6,5,14,12},8/1),({6,5,13,13},3/1)}, (4,2) => {}, (6,0) => {}, (5,1) => {({12,1,19,12},1/1),({12,1,18,13},1/1),({12,1,17,14},1/1),({11,2,22,9},1/1),({11,2,21,10},2/1),({11,2,20,11},3/1),({11,2,19,12},5/1),({11,2,18,13},6/1),({11,2,17,14},5/1),({11,2,16,15},3/1),({10,3,23,8},2/1),({10,3,22,9},4/1),({10,3,21,10},9/1),({10,3,20,11},13/1),({10,3,19,12},16/1),({10,3,18,13},17/1),({10,3,17,14},15/1),({10,3,16,15},8/1),({9,4,25,6},1/1),({9,4,24,7},3/1),({9,4,23,8},7/1),({9,4,22,9},13/1),({9,4,21,10},21/1),({9,4,20,11},29/1),({9,4,19,12},34/1),({9,4,18,13},34/1),({9,4,17,14},28/1),({9,4,16,15},16/1),({8,5,25,6},2/1),({8,5,24,7},5/1),({8,5,23,8},11/1),({8,5,22,9},20/1),({8,5,21,10},30/1),({8,5,20,11},38/1),({8,5,19,12},45/1),({8,5,18,13},43/1),({8,5,17,14},34/1),({8,5,16,15},20/1),({7,6,26,5},1/1),({7,6,25,6},2/1),({7,6,24,7},5/1),({7,6,23,8},10/1),({7,6,22,9},16/1),({7,6,21,10},23/1),({7,6,20,11},30/1),({7,6,19,12},33/1),({7,6,18,13},32/1),({7,6,17,14},26/1),({7,6,16,15},13/1)}};
--dw stands for dominant weights
dw1125 = new HashTable from {(6,1) => {({13,2,23,13},1/1),({12,3,25,11},1/1),({11,4,27,9},1/1),({10,5,28,8},1/1),({9,6,29,7},1/1)}, (7,0) => {}, (5,2) => {}, (7,1) => {({14,3,26,15},1/1),({13,4,28,13},1/1),({12,5,30,11},1/1),({11,6,31,10},1/1),({10,7,32,9},1/1)}, (8,0) => {}, (6,2) => {}, (8,1) => {({15,4,28,18},1/1),({14,5,31,15},1/1),({13,6,32,14},2/1),({12,7,34,12},1/1),({10,9,35,11},1/1)}, (9,0) => {}, (7,2) => {}, (10,0) => {}, (9,1) => {({16,5,29,22},1/1),({15,6,33,18},1/1),({14,7,35,16},1/1),({13,8,36,15},1/1),({12,9,37,14},1/1)}, (8,2) => {}, (11,0) => {}, (10,1) => {({17,6,29,27},1/1),({16,7,34,22},1/1),({15,8,37,19},1/1),({14,9,38,18},1/1),({13,10,39,17},1/1)}, (9,2) => {}, (12,0) => {}, (11,1) => {({17,8,34,27},1/1),({16,9,38,23},1/1),({15,10,40,21},1/1),({13,12,41,20},1/1)}, (10,2) => {}, (13,0) => {}, (12,1) => {({17,10,38,28},1/1),({16,11,41,25},1/1),({15,12,42,24},1/1)}, (11,2) => {}, (14,0) => {}, (13,1) => {({17,12,41,30},1/1),({16,13,43,28},1/1)}, (12,2) => {}, (15,0) => {}, (14,1) => {({17,14,43,33},1/1),({16,15,44,32},1/1)}, (13,2) => {}, (15,1) => {({17,16,44,37},1/1)}, (14,2) => {}, (15,2) => {}, (16,1) => {}, (17,1) => {}, (0,0) => {({1,0,1,0},1/1)}, (1,0) => {({3,0,5,1},1/1),({2,1,6,0},1/1)}, (1,1) => {}, (2,0) => {({4,1,10,1},1/1)}, (2,1) => {({7,0,11,5},1/1)}, (3,0) => {({5,2,14,2},1/1)}, (3,1) => {({9,0,13,8},1/1),({8,1,16,5},1/1),({6,3,18,3},1/1)}, (4,0) => {}, (2,2) => {}, (4,1) => {({11,0,14,12},1/1),({10,1,18,8},1/1),({9,2,20,6},1/1),({8,3,21,5},1/1),({7,4,22,4},1/1)}, (5,0) => {}, (3,2) => {}, (5,1) => {({12,1,19,12},1/1),({11,2,22,9},1/1),({10,3,23,8},2/1),({9,4,25,6},1/1),({7,6,26,5},1/1)}, (6,0) => {}, (4,2) => {}};
end;