A := QQ[t_0,t_1,t_2,t_3];
--tb stands for Total Betti numbers
tb1334 = new HashTable from {(5,2) => 0, (7,0) => 792/1, (7,1) => 116688/1, (9,0) => 0, (7,2) => 0, (9,1) => 233376/1, (9,2) => 0, (11,0) => 0, (11,1) => 148512/1, (11,2) => 0, (13,0) => 0, (13,1) => 32640/1, (15,0) => 0, (13,2) => 0, (15,1) => 2040/1, (15,2) => 0, (17,0) => 0, (17,1) => 16/1, (17,2) => 0, (19,1) => 0, (0,0) => 8/1, (2,0) => 816/1, (2,1) => 16/1, (2,2) => 0, (4,0) => 8400/1, (4,1) => 1680/1, (4,2) => 0, (6,0) => 8580/1, (6,1) => 43224/1, (6,2) => 0, (8,0) => 0, (8,1) => 194480/1, (8,2) => 0, (10,0) => 0, (10,1) => 212160/1, (10,2) => 0, (12,0) => 0, (12,1) => 79968/1, (12,2) => 0, (14,0) => 0, (14,1) => 9792/1, (14,2) => 0, (16,0) => 0, (16,1) => 264/1, (16,2) => 0, (18,1) => 0, (1,0) => 120/1, (1,1) => 0, (3,0) => 3280/1, (3,1) => 240/1, (5,0) => 13104/1, (3,2) => 0, (5,1) => 8580/1};
--sb represents the betti numbers as sums of Schur functors
sb1334 = new HashTable from {(5,2) => {}, (7,0) => {({18,4,19,12},1/1),({17,5,18,13},1/1),({16,6,19,12},1/1),({16,6,18,13},1/1),({16,6,17,14},1/1),({15,7,18,13},1/1),({15,7,17,14},1/1),({15,7,16,15},1/1),({14,8,19,12},1/1),({14,8,18,13},1/1),({14,8,17,14},2/1),({14,8,16,15},1/1),({13,9,18,13},1/1),({13,9,17,14},1/1),({13,9,16,15},1/1),({12,10,19,12},1/1),({12,10,18,13},1/1),({12,10,17,14},1/1),({12,10,16,15},1/1)}, (7,1) => {({20,5,22,13},1/1),({20,5,21,14},1/1),({20,5,20,15},2/1),({20,5,19,16},1/1),({20,5,18,17},1/1),({19,6,24,11},1/1),({19,6,23,12},2/1),({19,6,22,13},4/1),({19,6,21,14},6/1),({19,6,20,15},7/1),({19,6,19,16},6/1),({19,6,18,17},4/1),({18,7,26,9},1/1),({18,7,25,10},2/1),({18,7,24,11},6/1),({18,7,23,12},10/1),({18,7,22,13},16/1),({18,7,21,14},20/1),({18,7,20,15},22/1),({18,7,19,16},18/1),({18,7,18,17},11/1),({17,8,27,8},1/1),({17,8,26,9},3/1),({17,8,25,10},7/1),({17,8,24,11},15/1),({17,8,23,12},25/1),({17,8,22,13},36/1),({17,8,21,14},44/1),({17,8,20,15},46/1),({17,8,19,16},38/1),({17,8,18,17},22/1),({16,9,28,7},1/1),({16,9,27,8},3/1),({16,9,26,9},9/1),({16,9,25,10},18/1),({16,9,24,11},33/1),({16,9,23,12},49/1),({16,9,22,13},67/1),({16,9,21,14},77/1),({16,9,20,15},78/1),({16,9,19,16},63/1),({16,9,18,17},36/1),({15,10,28,7},1/1),({15,10,27,8},4/1),({15,10,26,9},12/1),({15,10,25,10},25/1),({15,10,24,11},44/1),({15,10,23,12},66/1),({15,10,22,13},87/1),({15,10,21,14},99/1),({15,10,20,15},98/1),({15,10,19,16},78/1),({15,10,18,17},44/1),({14,11,29,6},1/1),({14,11,28,7},3/1),({14,11,27,8},7/1),({14,11,26,9},16/1),({14,11,25,10},29/1),({14,11,24,11},49/1),({14,11,23,12},70/1),({14,11,22,13},90/1),({14,11,21,14},100/1),({14,11,20,15},97/1),({14,11,19,16},77/1),({14,11,18,17},43/1),({13,12,28,7},1/1),({13,12,27,8},4/1),({13,12,26,9},9/1),({13,12,25,10},18/1),({13,12,24,11},30/1),({13,12,23,12},43/1),({13,12,22,13},55/1),({13,12,21,14},61/1),({13,12,20,15},59/1),({13,12,19,16},47/1),({13,12,18,17},26/1)}, (9,0) => {}, (7,2) => {}, (9,1) => {({24,7,23,20},1/1),({23,8,27,16},1/1),({23,8,26,17},2/1),({23,8,25,18},3/1),({23,8,24,19},4/1),({23,8,23,20},4/1),({23,8,22,21},2/1),({22,9,29,14},1/1),({22,9,28,15},2/1),({22,9,27,16},6/1),({22,9,26,17},10/1),({22,9,25,18},15/1),({22,9,24,19},17/1),({22,9,23,20},16/1),({22,9,22,21},9/1),({21,10,30,13},1/1),({21,10,29,14},5/1),({21,10,28,15},11/1),({21,10,27,16},22/1),({21,10,26,17},35/1),({21,10,25,18},46/1),({21,10,24,19},50/1),({21,10,23,20},44/1),({21,10,22,21},25/1),({20,11,31,12},1/1),({20,11,30,13},4/1),({20,11,29,14},13/1),({20,11,28,15},28/1),({20,11,27,16},52/1),({20,11,26,17},77/1),({20,11,25,18},99/1),({20,11,24,19},104/1),({20,11,23,20},89/1),({20,11,22,21},51/1),({19,12,32,11},1/1),({19,12,31,12},4/1),({19,12,30,13},12/1),({19,12,29,14},29/1),({19,12,28,15},57/1),({19,12,27,16},95/1),({19,12,26,17},135/1),({19,12,25,18},166/1),({19,12,24,19},171/1),({19,12,23,20},144/1),({19,12,22,21},82/1),({18,13,32,11},1/1),({18,13,31,12},6/1),({18,13,30,13},17/1),({18,13,29,14},41/1),({18,13,28,15},77/1),({18,13,27,16},126/1),({18,13,26,17},175/1),({18,13,25,18},211/1),({18,13,24,19},215/1),({18,13,23,20},179/1),({18,13,22,21},101/1),({17,14,33,10},1/1),({17,14,32,11},3/1),({17,14,31,12},9/1),({17,14,30,13},22/1),({17,14,29,14},46/1),({17,14,28,15},82/1),({17,14,27,16},129/1),({17,14,26,17},175/1),({17,14,25,18},207/1),({17,14,24,19},209/1),({17,14,23,20},172/1),({17,14,22,21},97/1),({16,15,32,11},1/1),({16,15,31,12},5/1),({16,15,30,13},13/1),({16,15,29,14},29/1),({16,15,28,15},51/1),({16,15,27,16},80/1),({16,15,26,17},108/1),({16,15,25,18},127/1),({16,15,24,19},128/1),({16,15,23,20},105/1),({16,15,22,21},59/1)}, (9,2) => {}, (11,0) => {}, (11,1) => {({26,11,30,21},1/1),({26,11,29,22},1/1),({26,11,28,23},2/1),({26,11,27,24},2/1),({26,11,26,25},1/1),({25,12,32,19},1/1),({25,12,31,20},2/1),({25,12,30,21},5/1),({25,12,29,22},8/1),({25,12,28,23},10/1),({25,12,27,24},9/1),({25,12,26,25},6/1),({24,13,33,18},1/1),({24,13,32,19},5/1),({24,13,31,20},11/1),({24,13,30,21},20/1),({24,13,29,22},28/1),({24,13,28,23},33/1),({24,13,27,24},29/1),({24,13,26,25},18/1),({23,14,34,17},1/1),({23,14,33,18},5/1),({23,14,32,19},14/1),({23,14,31,20},29/1),({23,14,30,21},49/1),({23,14,29,22},65/1),({23,14,28,23},73/1),({23,14,27,24},64/1),({23,14,26,25},37/1),({22,15,35,16},1/1),({22,15,34,17},5/1),({22,15,33,18},13/1),({22,15,32,19},32/1),({22,15,31,20},58/1),({22,15,30,21},90/1),({22,15,29,22},117/1),({22,15,28,23},126/1),({22,15,27,24},107/1),({22,15,26,25},63/1),({21,16,35,16},1/1),({21,16,34,17},7/1),({21,16,33,18},20/1),({21,16,32,19},45/1),({21,16,31,20},80/1),({21,16,30,21},121/1),({21,16,29,22},153/1),({21,16,28,23},163/1),({21,16,27,24},138/1),({21,16,26,25},80/1),({20,17,36,15},1/1),({20,17,35,16},3/1),({20,17,34,17},10/1),({20,17,33,18},25/1),({20,17,32,19},50/1),({20,17,31,20},85/1),({20,17,30,21},125/1),({20,17,29,22},153/1),({20,17,28,23},161/1),({20,17,27,24},136/1),({20,17,26,25},77/1),({19,18,35,16},2/1),({19,18,34,17},6/1),({19,18,33,18},15/1),({19,18,32,19},32/1),({19,18,31,20},54/1),({19,18,30,21},77/1),({19,18,29,22},96/1),({19,18,28,23},100/1),({19,18,27,24},83/1),({19,18,26,25},49/1)}, (11,2) => {}, (13,0) => {}, (13,1) => {({28,15,33,26},1/1),({28,15,32,27},1/1),({28,15,31,28},1/1),({28,15,30,29},1/1),({27,16,35,24},1/1),({27,16,34,25},2/1),({27,16,33,26},5/1),({27,16,32,27},6/1),({27,16,31,28},6/1),({27,16,30,29},4/1),({26,17,36,23},1/1),({26,17,35,24},4/1),({26,17,34,25},8/1),({26,17,33,26},14/1),({26,17,32,27},17/1),({26,17,31,28},16/1),({26,17,30,29},10/1),({25,18,37,22},1/1),({25,18,36,23},4/1),({25,18,35,24},11/1),({25,18,34,25},20/1),({25,18,33,26},30/1),({25,18,32,27},35/1),({25,18,31,28},32/1),({25,18,30,29},19/1),({24,19,37,22},2/1),({24,19,36,23},7/1),({24,19,35,24},17/1),({24,19,34,25},30/1),({24,19,33,26},43/1),({24,19,32,27},49/1),({24,19,31,28},44/1),({24,19,30,29},26/1),({23,20,38,21},1/1),({23,20,37,22},4/1),({23,20,36,23},10/1),({23,20,35,24},21/1),({23,20,34,25},34/1),({23,20,33,26},47/1),({23,20,32,27},52/1),({23,20,31,28},46/1),({23,20,30,29},27/1),({22,21,37,22},2/1),({22,21,36,23},6/1),({22,21,35,24},13/1),({22,21,34,25},22/1),({22,21,33,26},30/1),({22,21,32,27},33/1),({22,21,31,28},29/1),({22,21,30,29},17/1)}, (15,0) => {}, (13,2) => {}, (15,1) => {({29,20,36,31},1/1),({29,20,35,32},1/1),({29,20,34,33},1/1),({28,21,38,29},1/1),({28,21,37,30},2/1),({28,21,36,31},3/1),({28,21,35,32},3/1),({28,21,34,33},2/1),({27,22,38,29},2/1),({27,22,37,30},3/1),({27,22,36,31},5/1),({27,22,35,32},5/1),({27,22,34,33},3/1),({26,23,39,28},1/1),({26,23,38,29},3/1),({26,23,37,30},5/1),({26,23,36,31},7/1),({26,23,35,32},7/1),({26,23,34,33},4/1),({25,24,39,28},1/1),({25,24,38,29},2/1),({25,24,37,30},3/1),({25,24,36,31},4/1),({25,24,35,32},4/1),({25,24,34,33},2/1)}, (15,2) => {}, (17,0) => {}, (17,1) => {({29,26,39,36},1/1)}, (17,2) => {}, (19,1) => {}, (0,0) => {({1,0,3,0},1/1)}, (2,0) => {({7,0,8,3},1/1),({7,0,7,4},1/1),({7,0,6,5},1/1),({6,1,10,1},1/1),({6,1,9,2},2/1),({6,1,8,3},3/1),({6,1,7,4},3/1),({6,1,6,5},2/1),({5,2,10,1},1/1),({5,2,9,2},2/1),({5,2,8,3},3/1),({5,2,7,4},3/1),({5,2,6,5},2/1),({4,3,10,1},1/1),({4,3,9,2},2/1),({4,3,8,3},3/1),({4,3,7,4},3/1),({4,3,6,5},2/1)}, (2,1) => {({5,5,15,0},1/1)}, (2,2) => {}, (4,0) => {({12,1,13,6},1/1),({12,1,12,7},1/1),({12,1,11,8},1/1),({12,1,10,9},1/1),({11,2,15,4},1/1),({11,2,14,5},2/1),({11,2,13,6},4/1),({11,2,12,7},5/1),({11,2,11,8},5/1),({11,2,10,9},3/1),({10,3,15,4},2/1),({10,3,14,5},4/1),({10,3,13,6},8/1),({10,3,12,7},10/1),({10,3,11,8},10/1),({10,3,10,9},6/1),({9,4,16,3},1/1),({9,4,15,4},4/1),({9,4,14,5},8/1),({9,4,13,6},14/1),({9,4,12,7},17/1),({9,4,11,8},16/1),({9,4,10,9},10/1),({8,5,16,3},1/1),({8,5,15,4},4/1),({8,5,14,5},8/1),({8,5,13,6},14/1),({8,5,12,7},17/1),({8,5,11,8},16/1),({8,5,10,9},10/1),({7,6,16,3},1/1),({7,6,15,4},3/1),({7,6,14,5},6/1),({7,6,13,6},10/1),({7,6,12,7},12/1),({7,6,11,8},11/1),({7,6,10,9},7/1)}, (4,1) => {({11,5,20,3},1/1),({11,5,19,4},1/1),({11,5,18,5},2/1),({11,5,17,6},1/1),({11,5,16,7},1/1),({10,6,21,2},1/1),({10,6,20,3},1/1),({10,6,19,4},2/1),({10,6,18,5},2/1),({10,6,17,6},2/1),({10,6,16,7},1/1),({10,6,15,8},1/1),({9,7,20,3},1/1),({9,7,19,4},1/1),({9,7,18,5},2/1),({9,7,17,6},1/1),({9,7,16,7},1/1),({8,8,21,2},1/1),({8,8,20,3},1/1),({8,8,19,4},2/1),({8,8,18,5},2/1),({8,8,17,6},2/1),({8,8,16,7},1/1),({8,8,15,8},1/1)}, (4,2) => {}, (6,0) => {({16,3,18,9},1/1),({16,3,17,10},1/1),({16,3,16,11},1/1),({16,3,15,12},1/1),({15,4,18,9},1/1),({15,4,17,10},2/1),({15,4,16,11},3/1),({15,4,15,12},2/1),({15,4,14,13},1/1),({14,5,20,7},1/1),({14,5,19,8},1/1),({14,5,18,9},3/1),({14,5,17,10},5/1),({14,5,16,11},7/1),({14,5,15,12},6/1),({14,5,14,13},4/1),({13,6,20,7},1/1),({13,6,19,8},2/1),({13,6,18,9},5/1),({13,6,17,10},7/1),({13,6,16,11},10/1),({13,6,15,12},10/1),({13,6,14,13},6/1),({12,7,20,7},1/1),({12,7,19,8},3/1),({12,7,18,9},7/1),({12,7,17,10},10/1),({12,7,16,11},14/1),({12,7,15,12},13/1),({12,7,14,13},9/1),({11,8,20,7},1/1),({11,8,19,8},2/1),({11,8,18,9},6/1),({11,8,17,10},10/1),({11,8,16,11},13/1),({11,8,15,12},12/1),({11,8,14,13},8/1),({10,9,19,8},1/1),({10,9,18,9},4/1),({10,9,17,10},6/1),({10,9,16,11},9/1),({10,9,15,12},8/1),({10,9,14,13},4/1)}, (6,1) => {({18,4,20,11},1/1),({18,4,19,12},1/1),({18,4,18,13},1/1),({18,4,17,14},1/1),({17,5,21,10},1/1),({17,5,20,11},2/1),({17,5,19,12},3/1),({17,5,18,13},4/1),({17,5,17,14},3/1),({17,5,16,15},1/1),({16,6,24,7},1/1),({16,6,23,8},2/1),({16,6,22,9},4/1),({16,6,21,10},5/1),({16,6,20,11},8/1),({16,6,19,12},9/1),({16,6,18,13},10/1),({16,6,17,14},8/1),({16,6,16,15},4/1),({15,7,25,6},1/1),({15,7,24,7},2/1),({15,7,23,8},6/1),({15,7,22,9},9/1),({15,7,21,10},14/1),({15,7,20,11},17/1),({15,7,19,12},19/1),({15,7,18,13},18/1),({15,7,17,14},15/1),({15,7,16,15},8/1),({14,8,26,5},1/1),({14,8,25,6},2/1),({14,8,24,7},6/1),({14,8,23,8},10/1),({14,8,22,9},18/1),({14,8,21,10},23/1),({14,8,20,11},29/1),({14,8,19,12},29/1),({14,8,18,13},28/1),({14,8,17,14},21/1),({14,8,16,15},12/1),({13,9,26,5},1/1),({13,9,25,6},4/1),({13,9,24,7},8/1),({13,9,23,8},16/1),({13,9,22,9},23/1),({13,9,21,10},32/1),({13,9,20,11},35/1),({13,9,19,12},37/1),({13,9,18,13},32/1),({13,9,17,14},25/1),({13,9,16,15},13/1),({12,10,26,5},1/1),({12,10,25,6},2/1),({12,10,24,7},7/1),({12,10,23,8},11/1),({12,10,22,9},19/1),({12,10,21,10},23/1),({12,10,20,11},29/1),({12,10,19,12},28/1),({12,10,18,13},26/1),({12,10,17,14},18/1),({12,10,16,15},10/1),({11,11,27,4},1/1),({11,11,26,5},1/1),({11,11,25,6},3/1),({11,11,24,7},4/1),({11,11,23,8},8/1),({11,11,22,9},9/1),({11,11,21,10},14/1),({11,11,20,11},13/1),({11,11,19,12},15/1),({11,11,18,13},11/1),({11,11,17,14},9/1),({11,11,16,15},4/1)}, (6,2) => {}, (8,0) => {}, (8,1) => {({22,6,23,16},1/1),({22,6,22,17},1/1),({22,6,21,18},1/1),({22,6,20,19},1/1),({21,7,26,13},1/1),({21,7,25,14},2/1),({21,7,24,15},4/1),({21,7,23,16},6/1),({21,7,22,17},7/1),({21,7,21,18},6/1),({21,7,20,19},4/1),({20,8,27,12},2/1),({20,8,26,13},4/1),({20,8,25,14},10/1),({20,8,24,15},16/1),({20,8,23,16},22/1),({20,8,22,17},24/1),({20,8,21,18},22/1),({20,8,20,19},12/1),({19,9,29,10},1/1),({19,9,28,11},3/1),({19,9,27,12},7/1),({19,9,26,13},17/1),({19,9,25,14},29/1),({19,9,24,15},44/1),({19,9,23,16},56/1),({19,9,22,17},60/1),({19,9,21,18},50/1),({19,9,20,19},30/1),({18,10,29,10},2/1),({18,10,28,11},7/1),({18,10,27,12},19/1),({18,10,26,13},36/1),({18,10,25,14},62/1),({18,10,24,15},86/1),({18,10,23,16},107/1),({18,10,22,17},109/1),({18,10,21,18},92/1),({18,10,20,19},52/1),({17,11,30,9},2/1),({17,11,29,10},5/1),({17,11,28,11},15/1),({17,11,27,12},32/1),({17,11,26,13},60/1),({17,11,25,14},94/1),({17,11,24,15},132/1),({17,11,23,16},154/1),({17,11,22,17},158/1),({17,11,21,18},130/1),({17,11,20,19},73/1),({16,12,31,8},1/1),({16,12,30,9},2/1),({16,12,29,10},9/1),({16,12,28,11},20/1),({16,12,27,12},43/1),({16,12,26,13},74/1),({16,12,25,14},115/1),({16,12,24,15},150/1),({16,12,23,16},179/1),({16,12,22,17},176/1),({16,12,21,18},144/1),({16,12,20,19},82/1),({15,13,30,9},2/1),({15,13,29,10},6/1),({15,13,28,11},17/1),({15,13,27,12},33/1),({15,13,26,13},61/1),({15,13,25,14},90/1),({15,13,24,15},122/1),({15,13,23,16},139/1),({15,13,22,17},140/1),({15,13,21,18},112/1),({15,13,20,19},64/1),({14,14,31,8},1/1),({14,14,30,9},2/1),({14,14,29,10},5/1),({14,14,28,11},8/1),({14,14,27,12},18/1),({14,14,26,13},25/1),({14,14,25,14},40/1),({14,14,24,15},50/1),({14,14,23,16},58/1),({14,14,22,17},55/1),({14,14,21,18},47/1),({14,14,20,19},24/1)}, (8,2) => {}, (10,0) => {}, (10,1) => {({25,9,27,20},1/1),({25,9,26,21},1/1),({25,9,25,22},1/1),({25,9,24,23},1/1),({24,10,30,17},1/1),({24,10,29,18},2/1),({24,10,28,19},4/1),({24,10,27,20},6/1),({24,10,26,21},8/1),({24,10,25,22},7/1),({24,10,24,23},5/1),({23,11,31,16},1/1),({23,11,30,17},4/1),({23,11,29,18},10/1),({23,11,28,19},17/1),({23,11,27,20},25/1),({23,11,26,21},28/1),({23,11,25,22},25/1),({23,11,24,23},15/1),({22,12,32,15},2/1),({22,12,31,16},6/1),({22,12,30,17},16/1),({22,12,29,18},30/1),({22,12,28,19},49/1),({22,12,27,20},64/1),({22,12,26,21},71/1),({22,12,25,22},61/1),({22,12,24,23},36/1),({21,13,33,14},1/1),({21,13,32,15},5/1),({21,13,31,16},16/1),({21,13,30,17},35/1),({21,13,29,18},65/1),({21,13,28,19},97/1),({21,13,27,20},125/1),({21,13,26,21},132/1),({21,13,25,22},113/1),({21,13,24,23},65/1),({20,14,34,13},1/1),({20,14,33,14},3/1),({20,14,32,15},12/1),({20,14,31,16},28/1),({20,14,30,17},60/1),({20,14,29,18},101/1),({20,14,28,19},149/1),({20,14,27,20},183/1),({20,14,26,21},193/1),({20,14,25,22},161/1),({20,14,24,23},93/1),({19,15,34,13},1/1),({19,15,33,14},6/1),({19,15,32,15},16/1),({19,15,31,16},39/1),({19,15,30,17},74/1),({19,15,29,18},123/1),({19,15,28,19},172/1),({19,15,27,20},211/1),({19,15,26,21},216/1),({19,15,25,22},181/1),({19,15,24,23},103/1),({18,16,34,13},1/1),({18,16,33,14},4/1),({18,16,32,15},14/1),({18,16,31,16},31/1),({18,16,30,17},62/1),({18,16,29,18},98/1),({18,16,28,19},140/1),({18,16,27,20},167/1),({18,16,26,21},173/1),({18,16,25,22},142/1),({18,16,24,23},82/1),({17,17,35,12},1/1),({17,17,34,13},1/1),({17,17,33,14},4/1),({17,17,32,15},7/1),({17,17,31,16},16/1),({17,17,30,17},26/1),({17,17,29,18},43/1),({17,17,28,19},56/1),({17,17,27,20},69/1),({17,17,26,21},68/1),({17,17,25,22},57/1),({17,17,24,23},32/1)}, (10,2) => {}, (12,0) => {}, (12,1) => {({27,13,32,23},1/1),({27,13,31,24},1/1),({27,13,30,25},2/1),({27,13,29,26},2/1),({27,13,28,27},1/1),({26,14,33,22},2/1),({26,14,32,23},4/1),({26,14,31,24},7/1),({26,14,30,25},9/1),({26,14,29,26},9/1),({26,14,28,27},5/1),({25,15,35,20},1/1),({25,15,34,21},3/1),({25,15,33,22},8/1),({25,15,32,23},16/1),({25,15,31,24},23/1),({25,15,30,25},28/1),({25,15,29,26},25/1),({25,15,28,27},15/1),({24,16,35,20},3/1),({24,16,34,21},9/1),({24,16,33,22},21/1),({24,16,32,23},36/1),({24,16,31,24},51/1),({24,16,30,25},57/1),({24,16,29,26},52/1),({24,16,28,27},30/1),({23,17,36,19},2/1),({23,17,35,20},7/1),({23,17,34,21},19/1),({23,17,33,22},37/1),({23,17,32,23},61/1),({23,17,31,24},81/1),({23,17,30,25},90/1),({23,17,29,26},78/1),({23,17,28,27},46/1),({22,18,37,18},1/1),({22,18,36,19},3/1),({22,18,35,20},11/1),({22,18,34,21},25/1),({22,18,33,22},49/1),({22,18,32,23},75/1),({22,18,31,24},99/1),({22,18,30,25},106/1),({22,18,29,26},92/1),({22,18,28,27},53/1),({21,19,36,19},3/1),({21,19,35,20},9/1),({21,19,34,21},22/1),({21,19,33,22},40/1),({21,19,32,23},63/1),({21,19,31,24},80/1),({21,19,30,25},87/1),({21,19,29,26},74/1),({21,19,28,27},43/1),({20,20,37,18},1/1),({20,20,36,19},2/1),({20,20,35,20},6/1),({20,20,34,21},10/1),({20,20,33,22},19/1),({20,20,32,23},26/1),({20,20,31,24},34/1),({20,20,30,25},35/1),({20,20,29,26},30/1),({20,20,28,27},17/1)}, (12,2) => {}, (14,0) => {}, (14,1) => {({29,17,33,30},1/1),({28,18,36,27},1/1),({28,18,35,28},2/1),({28,18,34,29},3/1),({28,18,33,30},3/1),({28,18,32,31},2/1),({27,19,37,26},1/1),({27,19,36,27},3/1),({27,19,35,28},6/1),({27,19,34,29},8/1),({27,19,33,30},8/1),({27,19,32,31},5/1),({26,20,38,25},1/1),({26,20,37,26},3/1),({26,20,36,27},8/1),({26,20,35,28},12/1),({26,20,34,29},16/1),({26,20,33,30},15/1),({26,20,32,31},9/1),({25,21,38,25},2/1),({25,21,37,26},6/1),({25,21,36,27},11/1),({25,21,35,28},18/1),({25,21,34,29},21/1),({25,21,33,30},19/1),({25,21,32,31},12/1),({24,22,38,25},2/1),({24,22,37,26},5/1),({24,22,36,27},10/1),({24,22,35,28},15/1),({24,22,34,29},18/1),({24,22,33,30},16/1),({24,22,32,31},10/1),({23,23,39,24},1/1),({23,23,38,25},1/1),({23,23,37,26},3/1),({23,23,36,27},5/1),({23,23,35,28},7/1),({23,23,34,29},8/1),({23,23,33,30},8/1),({23,23,32,31},4/1)}, (14,2) => {}, (16,0) => {}, (16,1) => {({29,23,38,33},1/1),({29,23,37,34},1/1),({29,23,36,35},1/1),({28,24,39,32},1/1),({28,24,38,33},1/1),({28,24,37,34},1/1),({28,24,36,35},1/1),({27,25,39,32},1/1),({27,25,38,33},1/1),({27,25,37,34},1/1),({27,25,36,35},1/1),({26,26,39,32},1/1),({26,26,38,33},1/1),({26,26,37,34},1/1),({26,26,36,35},1/1)}, (16,2) => {}, (18,1) => {}, (1,0) => {({4,0,6,1},1/1),({4,0,5,2},1/1),({4,0,4,3},1/1),({3,1,7,0},1/1),({3,1,6,1},1/1),({3,1,5,2},1/1),({3,1,4,3},1/1)}, (1,1) => {}, (3,0) => {({10,0,9,6},1/1),({9,1,12,3},1/1),({9,1,11,4},2/1),({9,1,10,5},3/1),({9,1,9,6},3/1),({9,1,8,7},2/1),({8,2,13,2},1/1),({8,2,12,3},2/1),({8,2,11,4},5/1),({8,2,10,5},6/1),({8,2,9,6},6/1),({8,2,8,7},4/1),({7,3,13,2},1/1),({7,3,12,3},4/1),({7,3,11,4},6/1),({7,3,10,5},9/1),({7,3,9,6},9/1),({7,3,8,7},5/1),({6,4,13,2},2/1),({6,4,12,3},4/1),({6,4,11,4},7/1),({6,4,10,5},9/1),({6,4,9,6},9/1),({6,4,8,7},5/1),({5,5,12,3},1/1),({5,5,11,4},2/1),({5,5,10,5},3/1),({5,5,9,6},2/1),({5,5,8,7},2/1)}, (3,1) => {({8,5,18,1},1/1),({8,5,17,2},1/1),({8,5,16,3},1/1),({8,5,15,4},1/1)}, (5,0) => {({14,2,16,7},1/1),({14,2,15,8},1/1),({14,2,14,9},2/1),({14,2,13,10},1/1),({14,2,12,11},1/1),({13,3,17,6},1/1),({13,3,16,7},2/1),({13,3,15,8},4/1),({13,3,14,9},5/1),({13,3,13,10},5/1),({13,3,12,11},3/1),({12,4,18,5},1/1),({12,4,17,6},2/1),({12,4,16,7},6/1),({12,4,15,8},9/1),({12,4,14,9},12/1),({12,4,13,10},11/1),({12,4,12,11},7/1),({11,5,18,5},2/1),({11,5,17,6},5/1),({11,5,16,7},10/1),({11,5,15,8},16/1),({11,5,14,9},19/1),({11,5,13,10},19/1),({11,5,12,11},11/1),({10,6,18,5},2/1),({10,6,17,6},5/1),({10,6,16,7},12/1),({10,6,15,8},18/1),({10,6,14,9},23/1),({10,6,13,10},21/1),({10,6,12,11},13/1),({9,7,19,4},1/1),({9,7,18,5},2/1),({9,7,17,6},6/1),({9,7,16,7},11/1),({9,7,15,8},17/1),({9,7,14,9},20/1),({9,7,13,10},19/1),({9,7,12,11},11/1),({8,8,18,5},1/1),({8,8,17,6},1/1),({8,8,16,7},4/1),({8,8,15,8},5/1),({8,8,14,9},7/1),({8,8,13,10},6/1),({8,8,12,11},4/1)}, (3,2) => {}, (5,1) => {({16,3,17,10},1/1),({15,4,17,10},1/1),({15,4,16,11},1/1),({14,5,21,6},1/1),({14,5,20,7},1/1),({14,5,19,8},1/1),({14,5,18,9},1/1),({14,5,17,10},1/1),({14,5,16,11},1/1),({14,5,15,12},1/1),({13,6,23,4},1/1),({13,6,22,5},2/1),({13,6,21,6},3/1),({13,6,20,7},5/1),({13,6,19,8},5/1),({13,6,18,9},4/1),({13,6,17,10},3/1),({13,6,16,11},2/1),({13,6,15,12},1/1),({13,6,14,13},1/1),({12,7,23,4},1/1),({12,7,22,5},2/1),({12,7,21,6},4/1),({12,7,20,7},5/1),({12,7,19,8},6/1),({12,7,18,9},5/1),({12,7,17,10},3/1),({12,7,16,11},2/1),({12,7,15,12},1/1),({12,7,14,13},1/1),({11,8,24,3},1/1),({11,8,23,4},2/1),({11,8,22,5},4/1),({11,8,21,6},7/1),({11,8,20,7},8/1),({11,8,19,8},8/1),({11,8,18,9},8/1),({11,8,17,10},5/1),({11,8,16,11},3/1),({11,8,15,12},2/1),({10,9,23,4},1/1),({10,9,22,5},2/1),({10,9,21,6},3/1),({10,9,20,7},4/1),({10,9,19,8},4/1),({10,9,18,9},3/1),({10,9,17,10},3/1),({10,9,16,11},2/1)}};
--dw stands for dominant weights
dw1334 = new HashTable from {(5,2) => {}, (7,0) => {({18,4,19,12},1/1)}, (7,1) => {({20,5,22,13},1/1),({19,6,24,11},1/1),({18,7,26,9},1/1),({17,8,27,8},1/1),({16,9,28,7},1/1),({14,11,29,6},1/1)}, (7,2) => {}, (9,0) => {}, (9,1) => {({24,7,23,20},1/1),({23,8,27,16},1/1),({22,9,29,14},1/1),({21,10,30,13},1/1),({20,11,31,12},1/1),({19,12,32,11},1/1),({17,14,33,10},1/1)}, (9,2) => {}, (11,0) => {}, (11,1) => {({26,11,30,21},1/1),({25,12,32,19},1/1),({24,13,33,18},1/1),({23,14,34,17},1/1),({22,15,35,16},1/1),({20,17,36,15},1/1)}, (11,2) => {}, (13,0) => {}, (13,1) => {({28,15,33,26},1/1),({27,16,35,24},1/1),({26,17,36,23},1/1),({25,18,37,22},1/1),({23,20,38,21},1/1)}, (13,2) => {}, (15,0) => {}, (15,1) => {({29,20,36,31},1/1),({28,21,38,29},1/1),({26,23,39,28},1/1)}, (17,0) => {}, (15,2) => {}, (17,1) => {({29,26,39,36},1/1)}, (17,2) => {}, (19,1) => {}, (0,0) => {({1,0,3,0},1/1)}, (2,0) => {({7,0,8,3},1/1),({6,1,10,1},1/1)}, (2,1) => {({5,5,15,0},1/1)}, (2,2) => {}, (4,0) => {({12,1,13,6},1/1),({11,2,15,4},1/1),({9,4,16,3},1/1)}, (4,1) => {({11,5,20,3},1/1),({10,6,21,2},1/1)}, (4,2) => {}, (6,0) => {({16,3,18,9},1/1),({14,5,20,7},1/1)}, (6,1) => {({18,4,20,11},1/1),({17,5,21,10},1/1),({16,6,24,7},1/1),({15,7,25,6},1/1),({14,8,26,5},1/1),({11,11,27,4},1/1)}, (6,2) => {}, (8,0) => {}, (8,1) => {({22,6,23,16},1/1),({21,7,26,13},1/1),({20,8,27,12},2/1),({19,9,29,10},1/1),({17,11,30,9},2/1),({16,12,31,8},1/1)}, (8,2) => {}, (10,0) => {}, (10,1) => {({25,9,27,20},1/1),({24,10,30,17},1/1),({23,11,31,16},1/1),({22,12,32,15},2/1),({21,13,33,14},1/1),({20,14,34,13},1/1),({17,17,35,12},1/1)}, (10,2) => {}, (12,0) => {}, (12,1) => {({27,13,32,23},1/1),({26,14,33,22},2/1),({25,15,35,20},1/1),({23,17,36,19},2/1),({22,18,37,18},1/1)}, (12,2) => {}, (14,0) => {}, (14,1) => {({29,17,33,30},1/1),({28,18,36,27},1/1),({27,19,37,26},1/1),({26,20,38,25},1/1),({23,23,39,24},1/1)}, (14,2) => {}, (16,0) => {}, (16,1) => {({29,23,38,33},1/1),({28,24,39,32},1/1)}, (16,2) => {}, (18,1) => {}, (1,0) => {({4,0,6,1},1/1),({3,1,7,0},1/1)}, (1,1) => {}, (3,0) => {({10,0,9,6},1/1),({9,1,12,3},1/1),({8,2,13,2},1/1)}, (3,1) => {({8,5,18,1},1/1)}, (3,2) => {}, (5,0) => {({14,2,16,7},1/1),({13,3,17,6},1/1),({12,4,18,5},1/1),({9,7,19,4},1/1)}, (5,1) => {({16,3,17,10},1/1),({14,5,21,6},1/1),({13,6,23,4},1/1),({11,8,24,3},1/1)}};
end;