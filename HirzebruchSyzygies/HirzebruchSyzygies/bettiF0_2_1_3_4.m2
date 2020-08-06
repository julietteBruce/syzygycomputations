A := QQ[t_0,t_1,t_2,t_3];
--tb stands for Total Betti numbers
tb2134 = new HashTable from {(5,2) => 0, (7,0) => 0, (7,1) => 204204/1, (7,2) => 0, (9,0) => 0, (9,1) => 320892/1, (9,2) => 0, (11,0) => 0, (11,1) => 185640/1, (11,2) => 0, (13,0) => 0, (13,1) => 38760/1, (13,2) => 0, (15,0) => 0, (15,1) => 2346/1, (17,0) => 0, (15,2) => 0, (17,1) => 18/1, (17,2) => 0, (19,1) => 0, (0,0) => 6/1, (2,0) => 510/1, (2,1) => 48/1, (2,2) => 0, (4,0) => 2730/1, (4,1) => 6504/1, (4,2) => 0, (6,0) => 0, (6,1) => 106080/1, (6,2) => 0, (8,0) => 0, (8,1) => 291720/1, (8,2) => 0, (10,0) => 0, (10,1) => 275808/1, (10,2) => 0, (12,0) => 0, (12,1) => 97104/1, (12,2) => 0, (14,0) => 0, (14,1) => 11424/1, (14,2) => 0, (16,0) => 0, (16,1) => 300/1, (16,2) => 0, (18,1) => 0, (1,0) => 84/1, (1,1) => 0, (3,0) => 1680/1, (3,1) => 690/1, (3,2) => 0, (5,0) => 792/1, (5,1) => 37128/1};
--sb represents the betti numbers as sums of Schur functors
sb2134 = new HashTable from {(5,2) => {}, (7,0) => {}, (7,1) => {({22,4,19,14},1/1),({22,4,18,15},1/1),({21,5,22,11},1/1),({21,5,21,12},2/1),({21,5,20,13},3/1),({21,5,19,14},4/1),({21,5,18,15},5/1),({21,5,17,16},2/1),({20,6,23,10},1/1),({20,6,22,11},4/1),({20,6,21,12},8/1),({20,6,20,13},13/1),({20,6,19,14},16/1),({20,6,18,15},15/1),({20,6,17,16},9/1),({19,7,24,9},2/1),({19,7,23,10},5/1),({19,7,22,11},14/1),({19,7,21,12},25/1),({19,7,20,13},36/1),({19,7,19,14},41/1),({19,7,18,15},38/1),({19,7,17,16},22/1),({18,8,25,8},1/1),({18,8,24,9},6/1),({18,8,23,10},16/1),({18,8,22,11},33/1),({18,8,21,12},55/1),({18,8,20,13},75/1),({18,8,19,14},84/1),({18,8,18,15},74/1),({18,8,17,16},43/1),({17,9,26,7},1/1),({17,9,25,8},4/1),({17,9,24,9},14/1),({17,9,23,10},32/1),({17,9,22,11},62/1),({17,9,21,12},95/1),({17,9,20,13},126/1),({17,9,19,14},135/1),({17,9,18,15},117/1),({17,9,17,16},68/1),({16,10,26,7},2/1),({16,10,25,8},8/1),({16,10,24,9},23/1),({16,10,23,10},50/1),({16,10,22,11},89/1),({16,10,21,12},133/1),({16,10,20,13},168/1),({16,10,19,14},178/1),({16,10,18,15},152/1),({16,10,17,16},87/1),({15,11,26,7},3/1),({15,11,25,8},10/1),({15,11,24,9},28/1),({15,11,23,10},56/1),({15,11,22,11},99/1),({15,11,21,12},143/1),({15,11,20,13},178/1),({15,11,19,14},185/1),({15,11,18,15},157/1),({15,11,17,16},89/1),({14,12,27,6},1/1),({14,12,26,7},3/1),({14,12,25,8},10/1),({14,12,24,9},24/1),({14,12,23,10},48/1),({14,12,22,11},80/1),({14,12,21,12},114/1),({14,12,20,13},139/1),({14,12,19,14},145/1),({14,12,18,15},120/1),({14,12,17,16},68/1),({13,13,26,7},1/1),({13,13,25,8},3/1),({13,13,24,9},9/1),({13,13,23,10},17/1),({13,13,22,11},30/1),({13,13,21,12},42/1),({13,13,20,13},52/1),({13,13,19,14},53/1),({13,13,18,15},45/1),({13,13,17,16},25/1)}, (7,2) => {}, (9,0) => {}, (9,1) => {({25,7,23,18},1/1),({25,7,22,19},1/1),({24,8,26,15},1/1),({24,8,25,16},2/1),({24,8,24,17},3/1),({24,8,23,18},5/1),({24,8,22,19},5/1),({24,8,21,20},3/1),({23,9,27,14},2/1),({23,9,26,15},5/1),({23,9,25,16},11/1),({23,9,24,17},16/1),({23,9,23,18},20/1),({23,9,22,19},18/1),({23,9,21,20},11/1),({22,10,29,12},1/1),({22,10,28,13},3/1),({22,10,27,14},9/1),({22,10,26,15},20/1),({22,10,25,16},35/1),({22,10,24,17},49/1),({22,10,23,18},56/1),({22,10,22,19},50/1),({22,10,21,20},30/1),({21,11,29,12},3/1),({21,11,28,13},10/1),({21,11,27,14},27/1),({21,11,26,15},51/1),({21,11,25,16},83/1),({21,11,24,17},109/1),({21,11,23,18},121/1),({21,11,22,19},104/1),({21,11,21,20},61/1),({20,12,30,11},2/1),({20,12,29,12},8/1),({20,12,28,13},24/1),({20,12,27,14},54/1),({20,12,26,15},98/1),({20,12,25,16},149/1),({20,12,24,17},191/1),({20,12,23,18},204/1),({20,12,22,19},174/1),({20,12,21,20},101/1),({19,13,31,10},1/1),({19,13,30,11},4/1),({19,13,29,12},16/1),({19,13,28,13},40/1),({19,13,27,14},85/1),({19,13,26,15},145/1),({19,13,25,16},214/1),({19,13,24,17},264/1),({19,13,23,18},279/1),({19,13,22,19},234/1),({19,13,21,20},135/1),({18,14,31,10},1/1),({18,14,30,11},6/1),({18,14,29,12},20/1),({18,14,28,13},49/1),({18,14,27,14},97/1),({18,14,26,15},164/1),({18,14,25,16},234/1),({18,14,24,17},287/1),({18,14,23,18},297/1),({18,14,22,19},249/1),({18,14,21,20},142/1),({17,15,31,10},2/1),({17,15,30,11},6/1),({17,15,29,12},19/1),({17,15,28,13},42/1),({17,15,27,14},83/1),({17,15,26,15},133/1),({17,15,25,16},189/1),({17,15,24,17},227/1),({17,15,23,18},235/1),({17,15,22,19},194/1),({17,15,21,20},110/1),({16,16,30,11},2/1),({16,16,29,12},6/1),({16,16,28,13},16/1),({16,16,27,14},30/1),({16,16,26,15},51/1),({16,16,25,16},70/1),({16,16,24,17},86/1),({16,16,23,18},87/1),({16,16,22,19},73/1),({16,16,21,20},41/1)}, (9,2) => {}, (11,0) => {}, (11,1) => {({27,11,28,21},1/1),({27,11,27,22},1/1),({27,11,26,23},1/1),({27,11,25,24},1/1),({26,12,30,19},1/1),({26,12,29,20},3/1),({26,12,28,21},6/1),({26,12,27,22},7/1),({26,12,26,23},7/1),({26,12,25,24},5/1),({25,13,32,17},1/1),({25,13,31,18},3/1),({25,13,30,19},8/1),({25,13,29,20},14/1),({25,13,28,21},23/1),({25,13,27,22},26/1),({25,13,26,23},24/1),({25,13,25,24},15/1),({24,14,33,16},1/1),({24,14,32,17},4/1),({24,14,31,18},12/1),({24,14,30,19},25/1),({24,14,29,20},42/1),({24,14,28,21},59/1),({24,14,27,22},66/1),({24,14,26,23},58/1),({24,14,25,24},35/1),({23,15,34,15},1/1),({23,15,33,16},3/1),({23,15,32,17},12/1),({23,15,31,18},28/1),({23,15,30,19},55/1),({23,15,29,20},85/1),({23,15,28,21},114/1),({23,15,27,22},123/1),({23,15,26,23},107/1),({23,15,25,24},62/1),({22,16,34,15},2/1),({22,16,33,16},8/1),({22,16,32,17},22/1),({22,16,31,18},49/1),({22,16,30,19},88/1),({22,16,29,20},133/1),({22,16,28,21},169/1),({22,16,27,22},180/1),({22,16,26,23},153/1),({22,16,25,24},89/1),({21,17,34,15},3/1),({21,17,33,16},10/1),({21,17,32,17},29/1),({21,17,31,18},59/1),({21,17,30,19},105/1),({21,17,29,20},152/1),({21,17,28,21},192/1),({21,17,27,22},200/1),({21,17,26,23},170/1),({21,17,25,24},97/1),({20,18,35,14},1/1),({20,18,34,15},3/1),({20,18,33,16},11/1),({20,18,32,17},26/1),({20,18,31,18},53/1),({20,18,30,19},88/1),({20,18,29,20},128/1),({20,18,28,21},156/1),({20,18,27,22},163/1),({20,18,26,23},136/1),({20,18,25,24},78/1),({19,19,34,15},1/1),({19,19,33,16},3/1),({19,19,32,17},10/1),({19,19,31,18},19/1),({19,19,30,19},34/1),({19,19,29,20},47/1),({19,19,28,21},60/1),({19,19,27,22},61/1),({19,19,26,23},52/1),({19,19,25,24},29/1)}, (11,2) => {}, (13,0) => {}, (13,1) => {({29,15,29,28},1/1),({28,16,33,24},1/1),({28,16,32,25},2/1),({28,16,31,26},3/1),({28,16,30,27},3/1),({28,16,29,28},2/1),({27,17,35,22},1/1),({27,17,34,23},2/1),({27,17,33,24},6/1),({27,17,32,25},9/1),({27,17,31,26},12/1),({27,17,30,27},11/1),({27,17,29,28},7/1),({26,18,36,21},1/1),({26,18,35,22},4/1),({26,18,34,23},9/1),({26,18,33,24},16/1),({26,18,32,25},24/1),({26,18,31,26},28/1),({26,18,30,27},25/1),({26,18,29,28},15/1),({25,19,37,20},1/1),({25,19,36,21},3/1),({25,19,35,22},9/1),({25,19,34,23},18/1),({25,19,33,24},32/1),({25,19,32,25},42/1),({25,19,31,26},48/1),({25,19,30,27},42/1),({25,19,29,28},25/1),({24,20,37,20},1/1),({24,20,36,21},5/1),({24,20,35,22},12/1),({24,20,34,23},25/1),({24,20,33,24},40/1),({24,20,32,25},53/1),({24,20,31,26},58/1),({24,20,30,27},51/1),({24,20,29,28},29/1),({23,21,37,20},2/1),({23,21,36,21},5/1),({23,21,35,22},13/1),({23,21,34,23},23/1),({23,21,33,24},37/1),({23,21,32,25},46/1),({23,21,31,26},51/1),({23,21,30,27},43/1),({23,21,29,28},25/1),({22,22,36,21},2/1),({22,22,35,22},4/1),({22,22,34,23},9/1),({22,22,33,24},13/1),({22,22,32,25},18/1),({22,22,31,26},19/1),({22,22,30,27},17/1),({22,22,29,28},9/1)}, (13,2) => {}, (15,0) => {}, (15,1) => {({29,21,36,29},1/1),({29,21,35,30},1/1),({29,21,34,31},2/1),({29,21,33,32},1/1),({28,22,38,27},1/1),({28,22,37,28},2/1),({28,22,36,29},3/1),({28,22,35,30},4/1),({28,22,34,31},4/1),({28,22,33,32},2/1),({27,23,38,27},2/1),({27,23,37,28},3/1),({27,23,36,29},5/1),({27,23,35,30},6/1),({27,23,34,31},6/1),({27,23,33,32},3/1),({26,24,39,26},1/1),({26,24,38,27},2/1),({26,24,37,28},4/1),({26,24,36,29},5/1),({26,24,35,30},7/1),({26,24,34,31},6/1),({26,24,33,32},3/1),({25,25,38,27},1/1),({25,25,37,28},1/1),({25,25,36,29},2/1),({25,25,35,30},2/1),({25,25,34,31},2/1),({25,25,33,32},1/1)}, (17,0) => {}, (15,2) => {}, (17,1) => {({29,27,39,34},1/1)}, (17,2) => {}, (19,1) => {}, (0,0) => {({2,0,1,0},1/1)}, (2,0) => {({7,1,8,1},1/1),({7,1,7,2},1/1),({7,1,6,3},1/1),({7,1,5,4},1/1),({6,2,8,1},2/1),({6,2,7,2},2/1),({6,2,6,3},2/1),({6,2,5,4},2/1),({5,3,9,0},1/1),({5,3,8,1},2/1),({5,3,7,2},2/1),({5,3,6,3},2/1),({5,3,5,4},2/1),({4,4,8,1},1/1),({4,4,7,2},1/1),({4,4,6,3},1/1),({4,4,5,4},1/1)}, (2,1) => {({11,0,8,5},1/1)}, (2,2) => {}, (4,0) => {({11,3,13,4},1/1),({11,3,12,5},1/1),({11,3,11,6},1/1),({11,3,10,7},1/1),({11,3,9,8},1/1),({10,4,15,2},1/1),({10,4,14,3},1/1),({10,4,13,4},2/1),({10,4,12,5},3/1),({10,4,11,6},4/1),({10,4,10,7},3/1),({10,4,9,8},2/1),({9,5,15,2},1/1),({9,5,14,3},2/1),({9,5,13,4},4/1),({9,5,12,5},6/1),({9,5,11,6},7/1),({9,5,10,7},6/1),({9,5,9,8},4/1),({8,6,15,2},1/1),({8,6,14,3},2/1),({8,6,13,4},4/1),({8,6,12,5},6/1),({8,6,11,6},7/1),({8,6,10,7},6/1),({8,6,9,8},4/1),({7,7,14,3},1/1),({7,7,13,4},2/1),({7,7,12,5},3/1),({7,7,11,6},3/1),({7,7,10,7},3/1),({7,7,9,8},2/1)}, (4,1) => {({16,1,13,8},1/1),({16,1,12,9},1/1),({15,2,15,6},1/1),({15,2,14,7},1/1),({15,2,13,8},3/1),({15,2,12,9},3/1),({15,2,11,10},1/1),({14,3,16,5},1/1),({14,3,15,6},2/1),({14,3,14,7},3/1),({14,3,13,8},5/1),({14,3,12,9},5/1),({14,3,11,10},3/1),({13,4,16,5},1/1),({13,4,15,6},4/1),({13,4,14,7},5/1),({13,4,13,8},7/1),({13,4,12,9},7/1),({13,4,11,10},4/1),({12,5,16,5},1/1),({12,5,15,6},3/1),({12,5,14,7},5/1),({12,5,13,8},6/1),({12,5,12,9},6/1),({12,5,11,10},3/1),({11,6,18,3},1/1),({11,6,17,4},1/1),({11,6,16,5},2/1),({11,6,15,6},3/1),({11,6,14,7},4/1),({11,6,13,8},5/1),({11,6,12,9},4/1),({11,6,11,10},2/1),({10,7,17,4},1/1),({10,7,16,5},1/1),({10,7,15,6},2/1),({10,7,14,7},2/1),({10,7,13,8},2/1),({10,7,12,9},2/1),({10,7,11,10},1/1),({9,8,16,5},1/1),({9,8,15,6},1/1),({9,8,14,7},1/1),({9,8,13,8},1/1),({9,8,12,9},1/1),({9,8,11,10},1/1)}, (4,2) => {}, (6,0) => {}, (6,1) => {({20,3,18,11},1/1),({20,3,17,12},1/1),({20,3,16,13},1/1),({20,3,15,14},1/1),({19,4,20,9},1/1),({19,4,19,10},2/1),({19,4,18,11},4/1),({19,4,17,12},5/1),({19,4,16,13},5/1),({19,4,15,14},3/1),({18,5,21,8},1/1),({18,5,20,9},4/1),({18,5,19,10},8/1),({18,5,18,11},13/1),({18,5,17,12},16/1),({18,5,16,13},15/1),({18,5,15,14},9/1),({17,6,22,7},1/1),({17,6,21,8},4/1),({17,6,20,9},11/1),({17,6,19,10},20/1),({17,6,18,11},30/1),({17,6,17,12},35/1),({17,6,16,13},32/1),({17,6,15,14},19/1),({16,7,23,6},1/1),({16,7,22,7},4/1),({16,7,21,8},11/1),({16,7,20,9},24/1),({16,7,19,10},40/1),({16,7,18,11},55/1),({16,7,17,12},62/1),({16,7,16,13},55/1),({16,7,15,14},32/1),({15,8,23,6},2/1),({15,8,22,7},8/1),({15,8,21,8},19/1),({15,8,20,9},38/1),({15,8,19,10},60/1),({15,8,18,11},80/1),({15,8,17,12},87/1),({15,8,16,13},76/1),({15,8,15,14},44/1),({14,9,24,5},1/1),({14,9,23,6},4/1),({14,9,22,7},12/1),({14,9,21,8},27/1),({14,9,20,9},49/1),({14,9,19,10},74/1),({14,9,18,11},95/1),({14,9,17,12},101/1),({14,9,16,13},86/1),({14,9,15,14},50/1),({13,10,24,5},1/1),({13,10,23,6},4/1),({13,10,22,7},12/1),({13,10,21,8},26/1),({13,10,20,9},46/1),({13,10,19,10},68/1),({13,10,18,11},86/1),({13,10,17,12},90/1),({13,10,16,13},76/1),({13,10,15,14},44/1),({12,11,24,5},1/1),({12,11,23,6},3/1),({12,11,22,7},8/1),({12,11,21,8},17/1),({12,11,20,9},29/1),({12,11,19,10},42/1),({12,11,18,11},52/1),({12,11,17,12},54/1),({12,11,16,13},45/1),({12,11,15,14},26/1)}, (6,2) => {}, (8,0) => {}, (8,1) => {({24,5,19,18},1/1),({23,6,23,14},1/1),({23,6,22,15},2/1),({23,6,21,16},2/1),({23,6,20,17},2/1),({23,6,19,18},2/1),({22,7,25,12},1/1),({22,7,24,13},2/1),({22,7,23,14},5/1),({22,7,22,15},9/1),({22,7,21,16},11/1),({22,7,20,17},10/1),({22,7,19,18},7/1),({21,8,26,11},1/1),({21,8,25,12},5/1),({21,8,24,13},11/1),({21,8,23,14},20/1),({21,8,22,15},30/1),({21,8,21,16},35/1),({21,8,20,17},31/1),({21,8,19,18},19/1),({20,9,27,10},1/1),({20,9,26,11},5/1),({20,9,25,12},15/1),({20,9,24,13},31/1),({20,9,23,14},53/1),({20,9,22,15},72/1),({20,9,21,16},81/1),({20,9,20,17},71/1),({20,9,19,18},42/1),({19,10,28,9},1/1),({19,10,27,10},5/1),({19,10,26,11},15/1),({19,10,25,12},36/1),({19,10,24,13},68/1),({19,10,23,14},106/1),({19,10,22,15},138/1),({19,10,21,16},150/1),({19,10,20,17},129/1),({19,10,19,18},76/1),({18,11,28,9},2/1),({18,11,27,10},10/1),({18,11,26,11},28/1),({18,11,25,12},62/1),({18,11,24,13},110/1),({18,11,23,14},165/1),({18,11,22,15},209/1),({18,11,21,16},222/1),({18,11,20,17},188/1),({18,11,19,18},109/1),({17,12,29,8},1/1),({17,12,28,9},5/1),({17,12,27,10},16/1),({17,12,26,11},41/1),({17,12,25,12},84/1),({17,12,24,13},142/1),({17,12,23,14},206/1),({17,12,22,15},255/1),({17,12,21,16},266/1),({17,12,20,17},223/1),({17,12,19,18},128/1),({16,13,29,8},1/1),({16,13,28,9},5/1),({16,13,27,10},17/1),({16,13,26,11},41/1),({16,13,25,12},82/1),({16,13,24,13},136/1),({16,13,23,14},194/1),({16,13,22,15},237/1),({16,13,21,16},244/1),({16,13,20,17},203/1),({16,13,19,18},117/1),({15,14,29,8},1/1),({15,14,28,9},4/1),({15,14,27,10},12/1),({15,14,26,11},27/1),({15,14,25,12},52/1),({15,14,24,13},85/1),({15,14,23,14},120/1),({15,14,22,15},145/1),({15,14,21,16},148/1),({15,14,20,17},123/1),({15,14,19,18},71/1)}, (8,2) => {}, (10,0) => {}, (10,1) => {({26,9,26,19},1/1),({26,9,25,20},1/1),({26,9,24,21},1/1),({26,9,23,22},1/1),({25,10,28,17},1/1),({25,10,27,18},3/1),({25,10,26,19},5/1),({25,10,25,20},7/1),({25,10,24,21},7/1),({25,10,23,22},4/1),({24,11,30,15},1/1),({24,11,29,16},3/1),({24,11,28,17},8/1),({24,11,27,18},15/1),({24,11,26,19},22/1),({24,11,25,20},26/1),({24,11,24,21},25/1),({24,11,23,22},14/1),({23,12,31,14},1/1),({23,12,30,15},4/1),({23,12,29,16},12/1),({23,12,28,17},26/1),({23,12,27,18},44/1),({23,12,26,19},61/1),({23,12,25,20},69/1),({23,12,24,21},61/1),({23,12,23,22},36/1),({22,13,32,13},1/1),({22,13,31,14},4/1),({22,13,30,15},14/1),({22,13,29,16},32/1),({22,13,28,17},62/1),({22,13,27,18},97/1),({22,13,26,19},127/1),({22,13,25,20},138/1),({22,13,24,21},120/1),({22,13,23,22},69/1),({21,14,32,13},2/1),({21,14,31,14},9/1),({21,14,30,15},27/1),({21,14,29,16},59/1),({21,14,28,17},106/1),({21,14,27,18},160/1),({21,14,26,19},203/1),({21,14,25,20},216/1),({21,14,24,21},184/1),({21,14,23,22},106/1),({20,15,33,12},1/1),({20,15,32,13},5/1),({20,15,31,14},16/1),({20,15,30,15},41/1),({20,15,29,16},83/1),({20,15,28,17},143/1),({20,15,27,18},207/1),({20,15,26,19},258/1),({20,15,25,20},269/1),({20,15,24,21},226/1),({20,15,23,22},130/1),({19,16,33,12},1/1),({19,16,32,13},5/1),({19,16,31,14},17/1),({19,16,30,15},42/1),({19,16,29,16},84/1),({19,16,28,17},140/1),({19,16,27,18},200/1),({19,16,26,19},245/1),({19,16,25,20},253/1),({19,16,24,21},212/1),({19,16,23,22},121/1),({18,17,33,12},1/1),({18,17,32,13},4/1),({18,17,31,14},12/1),({18,17,30,15},28/1),({18,17,29,16},54/1),({18,17,28,17},89/1),({18,17,27,18},125/1),({18,17,26,19},152/1),({18,17,25,20},156/1),({18,17,24,21},130/1),({18,17,23,22},74/1)}, (10,2) => {}, (12,0) => {}, (12,1) => {({28,13,29,24},1/1),({28,13,28,25},1/1),({27,14,32,21},1/1),({27,14,31,22},2/1),({27,14,30,23},4/1),({27,14,29,24},6/1),({27,14,28,25},6/1),({27,14,27,26},3/1),({26,15,33,20},2/1),({26,15,32,21},5/1),({26,15,31,22},11/1),({26,15,30,23},16/1),({26,15,29,24},20/1),({26,15,28,25},19/1),({26,15,27,26},11/1),({25,16,35,18},1/1),({25,16,34,19},3/1),({25,16,33,20},8/1),({25,16,32,21},18/1),({25,16,31,22},31/1),({25,16,30,23},43/1),({25,16,29,24},49/1),({25,16,28,25},44/1),({25,16,27,26},26/1),({24,17,35,18},2/1),({24,17,34,19},7/1),({24,17,33,20},18/1),({24,17,32,21},36/1),({24,17,31,22},59/1),({24,17,30,23},78/1),({24,17,29,24},87/1),({24,17,28,25},75/1),({24,17,27,26},44/1),({23,18,36,17},1/1),({23,18,35,18},4/1),({23,18,34,19},13/1),({23,18,33,20},30/1),({23,18,32,21},55/1),({23,18,31,22},84/1),({23,18,30,23},109/1),({23,18,29,24},117/1),({23,18,28,25},100/1),({23,18,27,26},58/1),({22,19,36,17},1/1),({22,19,35,18},5/1),({22,19,34,19},14/1),({22,19,33,20},32/1),({22,19,32,21},57/1),({22,19,31,22},86/1),({22,19,30,23},109/1),({22,19,29,24},116/1),({22,19,28,25},98/1),({22,19,27,26},57/1),({21,20,36,17},1/1),({21,20,35,18},4/1),({21,20,34,19},10/1),({21,20,33,20},22/1),({21,20,32,21},38/1),({21,20,31,22},56/1),({21,20,30,23},70/1),({21,20,29,24},73/1),({21,20,28,25},62/1),({21,20,27,26},36/1)}, (12,2) => {}, (14,0) => {}, (14,1) => {({29,18,33,28},1/1),({29,18,32,29},1/1),({29,18,31,30},1/1),({28,19,36,25},1/1),({28,19,35,26},2/1),({28,19,34,27},4/1),({28,19,33,28},5/1),({28,19,32,29},5/1),({28,19,31,30},3/1),({27,20,37,24},1/1),({27,20,36,25},3/1),({27,20,35,26},6/1),({27,20,34,27},10/1),({27,20,33,28},12/1),({27,20,32,29},11/1),({27,20,31,30},7/1),({26,21,38,23},1/1),({26,21,37,24},3/1),({26,21,36,25},7/1),({26,21,35,26},12/1),({26,21,34,27},18/1),({26,21,33,28},20/1),({26,21,32,29},18/1),({26,21,31,30},11/1),({25,22,38,23},1/1),({25,22,37,24},4/1),({25,22,36,25},8/1),({25,22,35,26},14/1),({25,22,34,27},20/1),({25,22,33,28},22/1),({25,22,32,29},19/1),({25,22,31,30},12/1),({24,23,38,23},1/1),({24,23,37,24},3/1),({24,23,36,25},6/1),({24,23,35,26},10/1),({24,23,34,27},14/1),({24,23,33,28},15/1),({24,23,32,29},13/1),({24,23,31,30},8/1)}, (14,2) => {}, (16,0) => {}, (16,1) => {({29,24,38,31},1/1),({29,24,37,32},1/1),({29,24,36,33},1/1),({29,24,35,34},1/1),({28,25,39,30},1/1),({28,25,38,31},1/1),({28,25,37,32},1/1),({28,25,36,33},1/1),({28,25,35,34},1/1),({27,26,39,30},1/1),({27,26,38,31},1/1),({27,26,37,32},1/1),({27,26,36,33},1/1),({27,26,35,34},1/1)}, (16,2) => {}, (18,1) => {}, (1,0) => {({5,0,4,1},1/1),({4,1,5,0},1/1),({4,1,4,1},1/1),({3,2,5,0},1/1),({3,2,4,1},1/1)}, (1,1) => {}, (3,0) => {({9,2,11,2},1/1),({9,2,10,3},1/1),({9,2,9,4},2/1),({9,2,8,5},2/1),({9,2,7,6},1/1),({8,3,12,1},1/1),({8,3,11,2},2/1),({8,3,10,3},3/1),({8,3,9,4},4/1),({8,3,8,5},4/1),({8,3,7,6},2/1),({7,4,12,1},1/1),({7,4,11,2},3/1),({7,4,10,3},4/1),({7,4,9,4},6/1),({7,4,8,5},6/1),({7,4,7,6},3/1),({6,5,12,1},1/1),({6,5,11,2},2/1),({6,5,10,3},3/1),({6,5,9,4},4/1),({6,5,8,5},4/1),({6,5,7,6},2/1)}, (3,1) => {({14,0,9,8},1/1),({13,1,12,5},1/1),({13,1,11,6},1/1),({13,1,10,7},1/1),({13,1,9,8},1/1),({12,2,12,5},1/1),({12,2,11,6},1/1),({12,2,10,7},1/1),({12,2,9,8},1/1),({11,3,12,5},1/1),({11,3,11,6},1/1),({11,3,10,7},1/1),({11,3,9,8},1/1)}, (3,2) => {}, (5,0) => {({11,6,18,3},1/1),({11,6,16,5},1/1),({11,6,15,6},1/1),({11,6,14,7},1/1),({11,6,13,8},1/1),({11,6,12,9},1/1),({10,7,17,4},1/1),({10,7,16,5},1/1),({10,7,15,6},2/1),({10,7,14,7},2/1),({10,7,13,8},2/1),({10,7,12,9},2/1),({10,7,11,10},1/1),({9,8,16,5},1/1),({9,8,15,6},1/1),({9,8,14,7},2/1),({9,8,13,8},2/1),({9,8,12,9},2/1),({9,8,11,10},1/1)}, (5,1) => {({18,2,16,9},1/1),({18,2,15,10},1/1),({18,2,14,11},1/1),({18,2,13,12},1/1),({17,3,17,8},2/1),({17,3,16,9},3/1),({17,3,15,10},4/1),({17,3,14,11},4/1),({17,3,13,12},3/1),({16,4,19,6},1/1),({16,4,18,7},2/1),({16,4,17,8},5/1),({16,4,16,9},9/1),({16,4,15,10},11/1),({16,4,14,11},10/1),({16,4,13,12},7/1),({15,5,19,6},2/1),({15,5,18,7},5/1),({15,5,17,8},11/1),({15,5,16,9},17/1),({15,5,15,10},21/1),({15,5,14,11},18/1),({15,5,13,12},12/1),({14,6,20,5},1/1),({14,6,19,6},5/1),({14,6,18,7},10/1),({14,6,17,8},18/1),({14,6,16,9},26/1),({14,6,15,10},30/1),({14,6,14,11},27/1),({14,6,13,12},16/1),({13,7,21,4},1/1),({13,7,20,5},2/1),({13,7,19,6},7/1),({13,7,18,7},14/1),({13,7,17,8},24/1),({13,7,16,9},32/1),({13,7,15,10},36/1),({13,7,14,11},31/1),({13,7,13,12},19/1),({12,8,21,4},1/1),({12,8,20,5},3/1),({12,8,19,6},8/1),({12,8,18,7},15/1),({12,8,17,8},24/1),({12,8,16,9},31/1),({12,8,15,10},34/1),({12,8,14,11},29/1),({12,8,13,12},17/1),({11,9,21,4},1/1),({11,9,20,5},3/1),({11,9,19,6},7/1),({11,9,18,7},12/1),({11,9,17,8},19/1),({11,9,16,9},23/1),({11,9,15,10},25/1),({11,9,14,11},21/1),({11,9,13,12},12/1),({10,10,20,5},1/1),({10,10,19,6},2/1),({10,10,18,7},4/1),({10,10,17,8},6/1),({10,10,16,9},8/1),({10,10,15,10},9/1),({10,10,14,11},7/1),({10,10,13,12},4/1)}};
--dw stands for dominant weights
dw2134 = new HashTable from {(5,2) => {}, (7,0) => {}, (7,1) => {({22,4,19,14},1/1),({21,5,22,11},1/1),({20,6,23,10},1/1),({19,7,24,9},2/1),({18,8,25,8},1/1),({17,9,26,7},1/1),({14,12,27,6},1/1)}, (9,0) => {}, (7,2) => {}, (9,1) => {({25,7,23,18},1/1),({24,8,26,15},1/1),({23,9,27,14},2/1),({22,10,29,12},1/1),({20,12,30,11},2/1),({19,13,31,10},1/1)}, (9,2) => {}, (11,0) => {}, (11,1) => {({27,11,28,21},1/1),({26,12,30,19},1/1),({25,13,32,17},1/1),({24,14,33,16},1/1),({23,15,34,15},1/1),({20,18,35,14},1/1)}, (11,2) => {}, (13,0) => {}, (13,1) => {({29,15,29,28},1/1),({28,16,33,24},1/1),({27,17,35,22},1/1),({26,18,36,21},1/1),({25,19,37,20},1/1)}, (15,0) => {}, (13,2) => {}, (15,1) => {({29,21,36,29},1/1),({28,22,38,27},1/1),({26,24,39,26},1/1)}, (15,2) => {}, (17,0) => {}, (17,1) => {({29,27,39,34},1/1)}, (17,2) => {}, (19,1) => {}, (0,0) => {({2,0,1,0},1/1)}, (2,0) => {({7,1,8,1},1/1),({5,3,9,0},1/1)}, (2,1) => {({11,0,8,5},1/1)}, (2,2) => {}, (4,0) => {({11,3,13,4},1/1),({10,4,15,2},1/1)}, (4,1) => {({16,1,13,8},1/1),({15,2,15,6},1/1),({14,3,16,5},1/1),({11,6,18,3},1/1)}, (4,2) => {}, (6,0) => {}, (6,1) => {({20,3,18,11},1/1),({19,4,20,9},1/1),({18,5,21,8},1/1),({17,6,22,7},1/1),({16,7,23,6},1/1),({14,9,24,5},1/1)}, (6,2) => {}, (8,0) => {}, (8,1) => {({24,5,19,18},1/1),({23,6,23,14},1/1),({22,7,25,12},1/1),({21,8,26,11},1/1),({20,9,27,10},1/1),({19,10,28,9},1/1),({17,12,29,8},1/1)}, (8,2) => {}, (10,0) => {}, (10,1) => {({26,9,26,19},1/1),({25,10,28,17},1/1),({24,11,30,15},1/1),({23,12,31,14},1/1),({22,13,32,13},1/1),({20,15,33,12},1/1)}, (10,2) => {}, (12,0) => {}, (12,1) => {({28,13,29,24},1/1),({27,14,32,21},1/1),({26,15,33,20},2/1),({25,16,35,18},1/1),({23,18,36,17},1/1)}, (12,2) => {}, (14,0) => {}, (14,1) => {({29,18,33,28},1/1),({28,19,36,25},1/1),({27,20,37,24},1/1),({26,21,38,23},1/1)}, (14,2) => {}, (16,0) => {}, (16,1) => {({29,24,38,31},1/1),({28,25,39,30},1/1)}, (16,2) => {}, (18,1) => {}, (1,0) => {({5,0,4,1},1/1),({4,1,5,0},1/1)}, (1,1) => {}, (3,0) => {({9,2,11,2},1/1),({8,3,12,1},1/1)}, (3,1) => {({14,0,9,8},1/1),({13,1,12,5},1/1)}, (5,0) => {({11,6,18,3},1/1)}, (3,2) => {}, (5,1) => {({18,2,16,9},1/1),({17,3,17,8},2/1),({16,4,19,6},1/1),({14,6,20,5},1/1),({13,7,21,4},1/1)}};
end;