A := QQ[t_0,t_1,t_2,t_3];
--tb stands for Total Betti numbers
tb0134 = new HashTable from {(5,2) => 0, (7,0) => 0, (7,1) => 251940/1, (7,2) => 0, (9,0) => 0, (9,1) => 301444/1, (9,2) => 0, (11,0) => 0, (11,1) => 132600/1, (11,2) => 0, (13,0) => 0, (13,1) => 16728/1, (13,2) => 0, (15,0) => 0, (15,1) => 58/1, (17,0) => 0, (15,2) => 240/1, (17,1) => 0, (17,2) => 4/1, (0,0) => 2/1, (2,0) => 0, (2,1) => 1020/1, (2,2) => 0, (4,0) => 0, (4,1) => 27744/1, (4,2) => 0, (6,0) => 0, (6,1) => 159120/1, (6,2) => 0, (8,0) => 0, (8,1) => 311168/1, (8,2) => 0, (10,0) => 0, (10,1) => 228072/1, (10,2) => 0, (12,0) => 0, (12,1) => 57120/1, (12,2) => 0, (14,0) => 0, (14,1) => 2448/1, (14,2) => 364/1, (16,0) => 0, (16,1) => 0, (16,2) => 50/1, (18,2) => 0, (1,0) => 16/1, (1,1) => 30/1, (3,0) => 0, (3,1) => 6936/1, (3,2) => 0, (5,0) => 0, (5,1) => 77112/1};
--sb represents the betti numbers as sums of Schur functors
sb0134 = new HashTable from {(5,2) => {}, (7,0) => {}, (7,1) => {({20,4,20,13},1/1),({20,4,19,14},2/1),({20,4,18,15},2/1),({20,4,17,16},1/1),({19,5,22,11},2/1),({19,5,21,12},5/1),({19,5,20,13},8/1),({19,5,19,14},10/1),({19,5,18,15},11/1),({19,5,17,16},6/1),({18,6,24,9},1/1),({18,6,23,10},4/1),({18,6,22,11},11/1),({18,6,21,12},21/1),({18,6,20,13},31/1),({18,6,19,14},37/1),({18,6,18,15},34/1),({18,6,17,16},20/1),({17,7,25,8},1/1),({17,7,24,9},6/1),({17,7,23,10},16/1),({17,7,22,11},35/1),({17,7,21,12},58/1),({17,7,20,13},80/1),({17,7,19,14},89/1),({17,7,18,15},79/1),({17,7,17,16},46/1),({16,8,26,7},1/1),({16,8,25,8},5/1),({16,8,24,9},17/1),({16,8,23,10},40/1),({16,8,22,11},75/1),({16,8,21,12},117/1),({16,8,20,13},152/1),({16,8,19,14},165/1),({16,8,18,15},142/1),({16,8,17,16},82/1),({15,9,26,7},3/1),({15,9,25,8},11/1),({15,9,24,9},32/1),({15,9,23,10},67/1),({15,9,22,11},120/1),({15,9,21,12},177/1),({15,9,20,13},224/1),({15,9,19,14},235/1),({15,9,18,15},200/1),({15,9,17,16},115/1),({14,10,27,6},1/1),({14,10,26,7},5/1),({14,10,25,8},17/1),({14,10,24,9},42/1),({14,10,23,10},85/1),({14,10,22,11},143/1),({14,10,21,12},206/1),({14,10,20,13},253/1),({14,10,19,14},263/1),({14,10,18,15},220/1),({14,10,17,16},126/1),({13,11,27,6},1/1),({13,11,26,7},5/1),({13,11,25,8},15/1),({13,11,24,9},37/1),({13,11,23,10},71/1),({13,11,22,11},119/1),({13,11,21,12},167/1),({13,11,20,13},204/1),({13,11,19,14},209/1),({13,11,18,15},175/1),({13,11,17,16},99/1),({12,12,27,6},1/1),({12,12,26,7},2/1),({12,12,25,8},7/1),({12,12,24,9},15/1),({12,12,23,10},30/1),({12,12,22,11},47/1),({12,12,21,12},66/1),({12,12,20,13},79/1),({12,12,19,14},82/1),({12,12,18,15},67/1),({12,12,17,16},38/1)}, (7,2) => {}, (9,0) => {}, (9,1) => {({24,6,21,20},1/1),({23,7,25,16},1/1),({23,7,24,17},2/1),({23,7,23,18},3/1),({23,7,22,19},3/1),({23,7,21,20},2/1),({22,8,27,14},1/1),({22,8,26,15},3/1),({22,8,25,16},7/1),({22,8,24,17},12/1),({22,8,23,18},15/1),({22,8,22,19},14/1),({22,8,21,20},9/1),({21,9,28,13},1/1),({21,9,27,14},7/1),({21,9,26,15},15/1),({21,9,25,16},28/1),({21,9,24,17},41/1),({21,9,23,18},48/1),({21,9,22,19},42/1),({21,9,21,20},26/1),({20,10,29,12},2/1),({20,10,28,13},8/1),({20,10,27,14},21/1),({20,10,26,15},44/1),({20,10,25,16},73/1),({20,10,24,17},98/1),({20,10,23,18},109/1),({20,10,22,19},96/1),({20,10,21,20},56/1),({19,11,30,11},1/1),({19,11,29,12},7/1),({19,11,28,13},21/1),({19,11,27,14},50/1),({19,11,26,15},91/1),({19,11,25,16},141/1),({19,11,24,17},180/1),({19,11,23,18},195/1),({19,11,22,19},166/1),({19,11,21,20},97/1),({18,12,30,11},4/1),({18,12,29,12},14/1),({18,12,28,13},38/1),({18,12,27,14},81/1),({18,12,26,15},141/1),({18,12,25,16},206/1),({18,12,24,17},260/1),({18,12,23,18},272/1),({18,12,22,19},229/1),({18,12,21,20},133/1),({17,13,31,10},2/1),({17,13,30,11},6/1),({17,13,29,12},21/1),({17,13,28,13},50/1),({17,13,27,14},100/1),({17,13,26,15},165/1),({17,13,25,16},238/1),({17,13,24,17},288/1),({17,13,23,18},300/1),({17,13,22,19},250/1),({17,13,21,20},142/1),({16,14,31,10},1/1),({16,14,30,11},6/1),({16,14,29,12},18/1),({16,14,28,13},43/1),({16,14,27,14},82/1),({16,14,26,15},136/1),({16,14,25,16},190/1),({16,14,24,17},231/1),({16,14,23,18},235/1),({16,14,22,19},196/1),({16,14,21,20},112/1),({15,15,31,10},1/1),({15,15,30,11},3/1),({15,15,29,12},9/1),({15,15,28,13},17/1),({15,15,27,14},35/1),({15,15,26,15},53/1),({15,15,25,16},75/1),({15,15,24,17},90/1),({15,15,23,18},92/1),({15,15,22,19},74/1),({15,15,21,20},44/1)}, (9,2) => {}, (11,0) => {}, (11,1) => {({26,10,28,21},1/1),({26,10,27,22},1/1),({26,10,26,23},1/1),({26,10,25,24},1/1),({25,11,30,19},1/1),({25,11,29,20},2/1),({25,11,28,21},5/1),({25,11,27,22},6/1),({25,11,26,23},6/1),({25,11,25,24},4/1),({24,12,31,18},1/1),({24,12,30,19},5/1),({24,12,29,20},11/1),({24,12,28,21},17/1),({24,12,27,22},21/1),({24,12,26,23},20/1),({24,12,25,24},12/1),({23,13,32,17},2/1),({23,13,31,18},7/1),({23,13,30,19},17/1),({23,13,29,20},29/1),({23,13,28,21},44/1),({23,13,27,22},50/1),({23,13,26,23},45/1),({23,13,25,24},27/1),({22,14,33,16},2/1),({22,14,32,17},7/1),({22,14,31,18},19/1),({22,14,30,19},38/1),({22,14,29,20},62/1),({22,14,28,21},83/1),({22,14,27,22},92/1),({22,14,26,23},80/1),({22,14,25,24},47/1),({21,15,34,15},1/1),({21,15,33,16},4/1),({21,15,32,17},15/1),({21,15,31,18},32/1),({21,15,30,19},62/1),({21,15,29,20},94/1),({21,15,28,21},122/1),({21,15,27,22},130/1),({21,15,26,23},113/1),({21,15,25,24},64/1),({20,16,34,15},2/1),({20,16,33,16},8/1),({20,16,32,17},20/1),({20,16,31,18},44/1),({20,16,30,19},75/1),({20,16,29,20},111/1),({20,16,28,21},139/1),({20,16,27,22},146/1),({20,16,26,23},122/1),({20,16,25,24},72/1),({19,17,34,15},2/1),({19,17,33,16},6/1),({19,17,32,17},18/1),({19,17,31,18},36/1),({19,17,30,19},63/1),({19,17,29,20},89/1),({19,17,28,21},112/1),({19,17,27,22},115/1),({19,17,26,23},97/1),({19,17,25,24},55/1),({18,18,35,14},1/1),({18,18,34,15},1/1),({18,18,33,16},4/1),({18,18,32,17},8/1),({18,18,31,18},16/1),({18,18,30,19},25/1),({18,18,29,20},38/1),({18,18,28,21},43/1),({18,18,27,22},46/1),({18,18,26,23},38/1),({18,18,25,24},21/1)}, (11,2) => {}, (13,0) => {}, (13,1) => {({28,14,31,26},1/1),({28,14,30,27},1/1),({27,15,33,24},1/1),({27,15,32,25},1/1),({27,15,31,26},3/1),({27,15,30,27},3/1),({27,15,29,28},2/1),({26,16,34,23},1/1),({26,16,33,24},3/1),({26,16,32,25},5/1),({26,16,31,26},7/1),({26,16,30,27},7/1),({26,16,29,28},4/1),({25,17,35,22},1/1),({25,17,34,23},3/1),({25,17,33,24},8/1),({25,17,32,25},11/1),({25,17,31,26},14/1),({25,17,30,27},13/1),({25,17,29,28},8/1),({24,18,36,21},1/1),({24,18,35,22},3/1),({24,18,34,23},7/1),({24,18,33,24},12/1),({24,18,32,25},18/1),({24,18,31,26},20/1),({24,18,30,27},18/1),({24,18,29,28},10/1),({23,19,37,20},1/1),({23,19,36,21},2/1),({23,19,35,22},6/1),({23,19,34,23},10/1),({23,19,33,24},17/1),({23,19,32,25},21/1),({23,19,31,26},24/1),({23,19,30,27},20/1),({23,19,29,28},12/1),({22,20,36,21},2/1),({22,20,35,22},4/1),({22,20,34,23},9/1),({22,20,33,24},13/1),({22,20,32,25},17/1),({22,20,31,26},18/1),({22,20,30,27},16/1),({22,20,29,28},9/1),({21,21,37,20},1/1),({21,21,36,21},1/1),({21,21,35,22},3/1),({21,21,34,23},4/1),({21,21,33,24},7/1),({21,21,32,25},7/1),({21,21,31,26},8/1),({21,21,30,27},6/1),({21,21,29,28},4/1)}, (13,2) => {}, (15,0) => {}, (15,1) => {({29,19,34,31},1/1),({24,24,39,26},1/1)}, (17,0) => {}, (15,2) => {({28,23,38,31},1/1),({28,23,37,32},1/1),({28,23,36,33},1/1),({28,23,35,34},1/1),({27,24,38,31},1/1),({27,24,37,32},1/1),({27,24,36,33},1/1),({27,24,35,34},1/1),({26,25,38,31},1/1),({26,25,37,32},1/1),({26,25,36,33},1/1),({26,25,35,34},1/1)}, (17,1) => {}, (17,2) => {({29,28,39,38},1/1)}, (0,0) => {({0,0,1,0},1/1)}, (2,0) => {}, (2,1) => {({9,0,8,5},1/1),({8,1,10,3},1/1),({8,1,9,4},1/1),({8,1,8,5},1/1),({8,1,7,6},1/1),({7,2,12,1},1/1),({7,2,11,2},1/1),({7,2,10,3},2/1),({7,2,9,4},2/1),({7,2,8,5},2/1),({7,2,7,6},1/1),({6,3,12,1},1/1),({6,3,11,2},2/1),({6,3,10,3},2/1),({6,3,9,4},3/1),({6,3,8,5},3/1),({6,3,7,6},1/1),({5,4,13,0},1/1),({5,4,12,1},1/1),({5,4,11,2},1/1),({5,4,10,3},2/1),({5,4,9,4},2/1),({5,4,8,5},1/1),({5,4,7,6},1/1)}, (2,2) => {}, (4,0) => {}, (4,1) => {({14,1,13,8},1/1),({14,1,12,9},1/1),({14,1,11,10},1/1),({13,2,15,6},2/1),({13,2,14,7},3/1),({13,2,13,8},5/1),({13,2,12,9},5/1),({13,2,11,10},3/1),({12,3,17,4},1/1),({12,3,16,5},3/1),({12,3,15,6},7/1),({12,3,14,7},11/1),({12,3,13,8},14/1),({12,3,12,9},13/1),({12,3,11,10},8/1),({11,4,18,3},1/1),({11,4,17,4},4/1),({11,4,16,5},9/1),({11,4,15,6},17/1),({11,4,14,7},24/1),({11,4,13,8},28/1),({11,4,12,9},25/1),({11,4,11,10},15/1),({10,5,19,2},1/1),({10,5,18,3},3/1),({10,5,17,4},8/1),({10,5,16,5},16/1),({10,5,15,6},27/1),({10,5,14,7},36/1),({10,5,13,8},40/1),({10,5,12,9},35/1),({10,5,11,10},20/1),({9,6,19,2},1/1),({9,6,18,3},4/1),({9,6,17,4},10/1),({9,6,16,5},19/1),({9,6,15,6},30/1),({9,6,14,7},38/1),({9,6,13,8},42/1),({9,6,12,9},36/1),({9,6,11,10},21/1),({8,7,19,2},1/1),({8,7,18,3},3/1),({8,7,17,4},7/1),({8,7,16,5},13/1),({8,7,15,6},20/1),({8,7,14,7},25/1),({8,7,13,8},27/1),({8,7,12,9},23/1),({8,7,11,10},13/1)}, (4,2) => {}, (6,0) => {}, (6,1) => {({18,3,18,11},2/1),({18,3,17,12},2/1),({18,3,16,13},2/1),({18,3,15,14},2/1),({17,4,20,9},2/1),({17,4,19,10},5/1),({17,4,18,11},9/1),({17,4,17,12},11/1),({17,4,16,13},11/1),({17,4,15,14},7/1),({16,5,22,7},1/1),({16,5,21,8},4/1),({16,5,20,9},11/1),({16,5,19,10},20/1),({16,5,18,11},31/1),({16,5,17,12},36/1),({16,5,16,13},33/1),({16,5,15,14},20/1),({15,6,23,6},1/1),({15,6,22,7},5/1),({15,6,21,8},15/1),({15,6,20,9},31/1),({15,6,19,10},51/1),({15,6,18,11},72/1),({15,6,17,12},80/1),({15,6,16,13},70/1),({15,6,15,14},42/1),({14,7,24,5},1/1),({14,7,23,6},4/1),({14,7,22,7},14/1),({14,7,21,8},32/1),({14,7,20,9},62/1),({14,7,19,10},96/1),({14,7,18,11},125/1),({14,7,17,12},135/1),({14,7,16,13},117/1),({14,7,15,14},67/1),({13,8,24,5},2/1),({13,8,23,6},8/1),({13,8,22,7},23/1),({13,8,21,8},49/1),({13,8,20,9},88/1),({13,8,19,10},131/1),({13,8,18,11},166/1),({13,8,17,12},175/1),({13,8,16,13},149/1),({13,8,15,14},86/1),({12,9,24,5},3/1),({12,9,23,6},10/1),({12,9,22,7},26/1),({12,9,21,8},54/1),({12,9,20,9},92/1),({12,9,19,10},132/1),({12,9,18,11},165/1),({12,9,17,12},171/1),({12,9,16,13},143/1),({12,9,15,14},83/1),({11,10,25,4},1/1),({11,10,24,5},2/1),({11,10,23,6},7/1),({11,10,22,7},18/1),({11,10,21,8},35/1),({11,10,20,9},59/1),({11,10,19,10},85/1),({11,10,18,11},103/1),({11,10,17,12},107/1),({11,10,16,13},90/1),({11,10,15,14},50/1)}, (6,2) => {}, (8,0) => {}, (8,1) => {({22,5,21,16},1/1),({22,5,20,17},1/1),({22,5,19,18},1/1),({21,6,24,13},1/1),({21,6,23,14},3/1),({21,6,22,15},5/1),({21,6,21,16},7/1),({21,6,20,17},7/1),({21,6,19,18},4/1),({20,7,25,12},3/1),({20,7,24,13},7/1),({20,7,23,14},15/1),({20,7,22,15},23/1),({20,7,21,16},28/1),({20,7,20,17},25/1),({20,7,19,18},16/1),({19,8,27,10},1/1),({19,8,26,11},4/1),({19,8,25,12},13/1),({19,8,24,13},28/1),({19,8,23,14},48/1),({19,8,22,15},67/1),({19,8,21,16},76/1),({19,8,20,17},67/1),({19,8,19,18},40/1),({18,9,27,10},4/1),({18,9,26,11},14/1),({18,9,25,12},35/1),({18,9,24,13},67/1),({18,9,23,14},108/1),({18,9,22,15},140/1),({18,9,21,16},154/1),({18,9,20,17},133/1),({18,9,19,18},77/1),({17,10,28,9},3/1),({17,10,27,10},11/1),({17,10,26,11},31/1),({17,10,25,12},69/1),({17,10,24,13},122/1),({17,10,23,14},182/1),({17,10,22,15},231/1),({17,10,21,16},245/1),({17,10,20,17},207/1),({17,10,19,18},121/1),({16,11,29,8},1/1),({16,11,28,9},5/1),({16,11,27,10},19/1),({16,11,26,11},47/1),({16,11,25,12},97/1),({16,11,24,13},164/1),({16,11,23,14},238/1),({16,11,22,15},293/1),({16,11,21,16},306/1),({16,11,20,17},256/1),({16,11,19,18},147/1),({15,12,29,8},2/1),({15,12,28,9},7/1),({15,12,27,10},22/1),({15,12,26,11},52/1),({15,12,25,12},100/1),({15,12,24,13},165/1),({15,12,23,14},234/1),({15,12,22,15},282/1),({15,12,21,16},291/1),({15,12,20,17},243/1),({15,12,19,18},137/1),({14,13,29,8},1/1),({14,13,28,9},5/1),({14,13,27,10},15/1),({14,13,26,11},33/1),({14,13,25,12},65/1),({14,13,24,13},104/1),({14,13,23,14},145/1),({14,13,22,15},176/1),({14,13,21,16},179/1),({14,13,20,17},147/1),({14,13,19,18},85/1)}, (8,2) => {}, (10,0) => {}, (10,1) => {({25,8,25,20},1/1),({25,8,24,21},1/1),({24,9,28,17},1/1),({24,9,27,18},2/1),({24,9,26,19},4/1),({24,9,25,20},5/1),({24,9,24,21},6/1),({24,9,23,22},3/1),({23,10,29,16},1/1),({23,10,28,17},5/1),({23,10,27,18},11/1),({23,10,26,19},17/1),({23,10,25,20},21/1),({23,10,24,21},20/1),({23,10,23,22},12/1),({22,11,30,15},3/1),({22,11,29,16},8/1),({22,11,28,17},20/1),({22,11,27,18},35/1),({22,11,26,19},50/1),({22,11,25,20},57/1),({22,11,24,21},52/1),({22,11,23,22},30/1),({21,12,31,14},2/1),({21,12,30,15},9/1),({21,12,29,16},24/1),({21,12,28,17},48/1),({21,12,27,18},78/1),({21,12,26,19},105/1),({21,12,25,20},116/1),({21,12,24,21},101/1),({21,12,23,22},59/1),({20,13,32,13},2/1),({20,13,31,14},7/1),({20,13,30,15},22/1),({20,13,29,16},48/1),({20,13,28,17},89/1),({20,13,27,18},134/1),({20,13,26,19},173/1),({20,13,25,20},184/1),({20,13,24,21},158/1),({20,13,23,22},91/1),({19,14,32,13},3/1),({19,14,31,14},12/1),({19,14,30,15},33/1),({19,14,29,16},69/1),({19,14,28,17},120/1),({19,14,27,18},176/1),({19,14,26,19},219/1),({19,14,25,20},229/1),({19,14,24,21},194/1),({19,14,23,22},111/1),({18,15,33,12},1/1),({18,15,32,13},5/1),({18,15,31,14},15/1),({18,15,30,15},37/1),({18,15,29,16},72/1),({18,15,28,17},122/1),({18,15,27,18},173/1),({18,15,26,19},212/1),({18,15,25,20},219/1),({18,15,24,21},183/1),({18,15,23,22},104/1),({17,16,33,12},1/1),({17,16,32,13},3/1),({17,16,31,14},10/1),({17,16,30,15},24/1),({17,16,29,16},47/1),({17,16,28,17},77/1),({17,16,27,18},108/1),({17,16,26,19},131/1),({17,16,25,20},135/1),({17,16,24,21},112/1),({17,16,23,22},63/1)}, (10,2) => {}, (12,0) => {}, (12,1) => {({27,12,30,23},1/1),({27,12,29,24},1/1),({27,12,28,25},1/1),({27,12,27,26},1/1),({26,13,31,22},2/1),({26,13,30,23},3/1),({26,13,29,24},5/1),({26,13,28,25},5/1),({26,13,27,26},3/1),({25,14,33,20},1/1),({25,14,32,21},3/1),({25,14,31,22},7/1),({25,14,30,23},12/1),({25,14,29,24},15/1),({25,14,28,25},14/1),({25,14,27,26},9/1),({24,15,33,20},4/1),({24,15,32,21},9/1),({24,15,31,22},18/1),({24,15,30,23},27/1),({24,15,29,24},32/1),({24,15,28,25},28/1),({24,15,27,26},18/1),({23,16,35,18},1/1),({23,16,34,19},4/1),({23,16,33,20},10/1),({23,16,32,21},21/1),({23,16,31,22},35/1),({23,16,30,23},47/1),({23,16,29,24},53/1),({23,16,28,25},47/1),({23,16,27,26},27/1),({22,17,35,18},2/1),({22,17,34,19},6/1),({22,17,33,20},16/1),({22,17,32,21},30/1),({22,17,31,22},48/1),({22,17,30,23},62/1),({22,17,29,24},68/1),({22,17,28,25},58/1),({22,17,27,26},34/1),({21,18,36,17},1/1),({21,18,35,18},3/1),({21,18,34,19},8/1),({21,18,33,20},19/1),({21,18,32,21},33/1),({21,18,31,22},49/1),({21,18,30,23},63/1),({21,18,29,24},66/1),({21,18,28,25},55/1),({21,18,27,26},33/1),({20,19,35,18},2/1),({20,19,34,19},6/1),({20,19,33,20},12/1),({20,19,32,21},21/1),({20,19,31,22},32/1),({20,19,30,23},38/1),({20,19,29,24},41/1),({20,19,28,25},35/1),({20,19,27,26},19/1)}, (12,2) => {}, (14,0) => {}, (14,1) => {({29,16,31,30},1/1),({28,17,34,27},1/1),({28,17,33,28},1/1),({28,17,32,29},1/1),({28,17,31,30},1/1),({27,18,34,27},1/1),({27,18,33,28},2/1),({27,18,32,29},2/1),({27,18,31,30},1/1),({26,19,36,25},1/1),({26,19,35,26},1/1),({26,19,34,27},3/1),({26,19,33,28},3/1),({26,19,32,29},3/1),({26,19,31,30},2/1),({25,20,36,25},1/1),({25,20,35,26},2/1),({25,20,34,27},3/1),({25,20,33,28},3/1),({25,20,32,29},3/1),({25,20,31,30},2/1),({24,21,38,23},1/1),({24,21,37,24},1/1),({24,21,36,25},2/1),({24,21,35,26},3/1),({24,21,34,27},3/1),({24,21,33,28},3/1),({24,21,32,29},3/1),({24,21,31,30},1/1),({23,22,37,24},1/1),({23,22,36,25},1/1),({23,22,35,26},1/1),({23,22,34,27},2/1),({23,22,33,28},2/1),({23,22,32,29},1/1),({23,22,31,30},1/1)}, (14,2) => {({27,21,37,28},1/1),({27,21,35,30},1/1),({27,21,34,31},1/1),({26,22,36,29},1/1),({26,22,35,30},1/1),({26,22,34,31},1/1),({26,22,33,32},1/1),({25,23,37,28},1/1),({25,23,36,29},1/1),({25,23,35,30},2/1),({25,23,34,31},2/1),({25,23,33,32},1/1),({24,24,34,31},1/1)}, (16,0) => {}, (16,1) => {}, (16,2) => {({29,25,38,35},1/1),({28,26,39,34},1/1),({28,26,38,35},1/1)}, (18,2) => {}, (1,0) => {({3,0,4,1},1/1)}, (1,1) => {({4,2,9,0},1/1)}, (3,0) => {}, (3,1) => {({12,0,9,8},1/1),({11,1,12,5},1/1),({11,1,11,6},2/1),({11,1,10,7},2/1),({11,1,9,8},1/1),({10,2,14,3},1/1),({10,2,13,4},2/1),({10,2,12,5},4/1),({10,2,11,6},5/1),({10,2,10,7},5/1),({10,2,9,8},3/1),({9,3,15,2},1/1),({9,3,14,3},3/1),({9,3,13,4},5/1),({9,3,12,5},9/1),({9,3,11,6},10/1),({9,3,10,7},9/1),({9,3,9,8},6/1),({8,4,16,1},1/1),({8,4,15,2},2/1),({8,4,14,3},5/1),({8,4,13,4},9/1),({8,4,12,5},12/1),({8,4,11,6},13/1),({8,4,10,7},12/1),({8,4,9,8},7/1),({7,5,16,1},1/1),({7,5,15,2},3/1),({7,5,14,3},5/1),({7,5,13,4},8/1),({7,5,12,5},11/1),({7,5,11,6},12/1),({7,5,10,7},10/1),({7,5,9,8},6/1),({6,6,15,2},1/1),({6,6,14,3},2/1),({6,6,13,4},3/1),({6,6,12,5},5/1),({6,6,11,6},5/1),({6,6,10,7},4/1),({6,6,9,8},3/1)}, (3,2) => {}, (5,0) => {}, (5,1) => {({16,2,16,9},1/1),({16,2,15,10},2/1),({16,2,14,11},2/1),({16,2,13,12},1/1),({15,3,18,7},1/1),({15,3,17,8},4/1),({15,3,16,9},6/1),({15,3,15,10},9/1),({15,3,14,11},9/1),({15,3,13,12},5/1),({14,4,19,6},3/1),({14,4,18,7},7/1),({14,4,17,8},14/1),({14,4,16,9},22/1),({14,4,15,10},26/1),({14,4,14,11},24/1),({14,4,13,12},15/1),({13,5,21,4},1/1),({13,5,20,5},3/1),({13,5,19,6},9/1),({13,5,18,7},20/1),({13,5,17,8},35/1),({13,5,16,9},48/1),({13,5,15,10},55/1),({13,5,14,11},48/1),({13,5,13,12},29/1),({12,6,21,4},2/1),({12,6,20,5},8/1),({12,6,19,6},19/1),({12,6,18,7},37/1),({12,6,17,8},59/1),({12,6,16,9},77/1),({12,6,15,10},85/1),({12,6,14,11},74/1),({12,6,13,12},42/1),({11,7,22,3},1/1),({11,7,21,4},4/1),({11,7,20,5},11/1),({11,7,19,6},27/1),({11,7,18,7},48/1),({11,7,17,8},73/1),({11,7,16,9},94/1),({11,7,15,10},100/1),({11,7,14,11},84/1),({11,7,13,12},50/1),({10,8,22,3},1/1),({10,8,21,4},4/1),({10,8,20,5},11/1),({10,8,19,6},24/1),({10,8,18,7},42/1),({10,8,17,8},62/1),({10,8,16,9},78/1),({10,8,15,10},82/1),({10,8,14,11},69/1),({10,8,13,12},40/1),({9,9,21,4},2/1),({9,9,20,5},5/1),({9,9,19,6},10/1),({9,9,18,7},17/1),({9,9,17,8},26/1),({9,9,16,9},30/1),({9,9,15,10},32/1),({9,9,14,11},28/1),({9,9,13,12},15/1)}};
--dw stands for dominant weights
dw0134 = new HashTable from {(5,2) => {}, (7,0) => {}, (7,1) => {({20,4,20,13},1/1),({19,5,22,11},2/1),({18,6,24,9},1/1),({17,7,25,8},1/1),({16,8,26,7},1/1),({14,10,27,6},1/1)}, (9,0) => {}, (7,2) => {}, (9,1) => {({24,6,21,20},1/1),({23,7,25,16},1/1),({22,8,27,14},1/1),({21,9,28,13},1/1),({20,10,29,12},2/1),({19,11,30,11},1/1),({17,13,31,10},2/1)}, (9,2) => {}, (11,0) => {}, (11,1) => {({26,10,28,21},1/1),({25,11,30,19},1/1),({24,12,31,18},1/1),({23,13,32,17},2/1),({22,14,33,16},2/1),({21,15,34,15},1/1),({18,18,35,14},1/1)}, (11,2) => {}, (13,0) => {}, (13,1) => {({28,14,31,26},1/1),({27,15,33,24},1/1),({26,16,34,23},1/1),({25,17,35,22},1/1),({24,18,36,21},1/1),({23,19,37,20},1/1)}, (15,0) => {}, (13,2) => {}, (15,1) => {({29,19,34,31},1/1),({24,24,39,26},1/1)}, (15,2) => {({28,23,38,31},1/1)}, (17,0) => {}, (17,1) => {}, (17,2) => {({29,28,39,38},1/1)}, (0,0) => {({0,0,1,0},1/1)}, (2,0) => {}, (2,1) => {({9,0,8,5},1/1),({8,1,10,3},1/1),({7,2,12,1},1/1),({5,4,13,0},1/1)}, (4,0) => {}, (2,2) => {}, (4,1) => {({14,1,13,8},1/1),({13,2,15,6},2/1),({12,3,17,4},1/1),({11,4,18,3},1/1),({10,5,19,2},1/1)}, (4,2) => {}, (6,0) => {}, (6,1) => {({18,3,18,11},2/1),({17,4,20,9},2/1),({16,5,22,7},1/1),({15,6,23,6},1/1),({14,7,24,5},1/1),({11,10,25,4},1/1)}, (6,2) => {}, (8,0) => {}, (8,1) => {({22,5,21,16},1/1),({21,6,24,13},1/1),({20,7,25,12},3/1),({19,8,27,10},1/1),({17,10,28,9},3/1),({16,11,29,8},1/1)}, (8,2) => {}, (10,0) => {}, (10,1) => {({25,8,25,20},1/1),({24,9,28,17},1/1),({23,10,29,16},1/1),({22,11,30,15},3/1),({21,12,31,14},2/1),({20,13,32,13},2/1),({18,15,33,12},1/1)}, (10,2) => {}, (12,0) => {}, (12,1) => {({27,12,30,23},1/1),({26,13,31,22},2/1),({25,14,33,20},1/1),({23,16,35,18},1/1),({21,18,36,17},1/1)}, (12,2) => {}, (14,0) => {}, (14,1) => {({29,16,31,30},1/1),({28,17,34,27},1/1),({26,19,36,25},1/1),({24,21,38,23},1/1)}, (14,2) => {({27,21,37,28},1/1)}, (16,0) => {}, (16,1) => {}, (16,2) => {({29,25,38,35},1/1),({28,26,39,34},1/1)}, (18,2) => {}, (1,0) => {({3,0,4,1},1/1)}, (1,1) => {({4,2,9,0},1/1)}, (3,0) => {}, (3,1) => {({12,0,9,8},1/1),({11,1,12,5},1/1),({10,2,14,3},1/1),({9,3,15,2},1/1),({8,4,16,1},1/1)}, (5,0) => {}, (3,2) => {}, (5,1) => {({16,2,16,9},1/1),({15,3,18,7},1/1),({14,4,19,6},3/1),({13,5,21,4},1/1),({11,7,22,3},1/1)}};
--dmw stands for dominant monomial weights
dmw0134 = {infinity};
end;