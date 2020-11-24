A := QQ[t_0,t_1,t_2,t_3];
--tb stands for Total Betti numbers
tb1234 = new HashTable from {(5,2) => 0, (7,0) => 0, (7,1) => 172380/1, (9,0) => 0, (7,2) => 0, (9,1) => 272272/1, (9,2) => 0, (11,0) => 0, (11,1) => 153816/1, (11,2) => 0, (13,0) => 0, (13,1) => 30192/1, (15,0) => 0, (13,2) => 0, (15,1) => 1530/1, (15,2) => 0, (17,0) => 0, (17,1) => 0, (17,2) => 1/1, (0,0) => 6/1, (2,0) => 528/1, (2,1) => 15/1, (4,0) => 3094/1, (2,2) => 0, (4,1) => 3939/1, (4,2) => 0, (6,0) => 0, (6,1) => 87516/1, (6,2) => 0, (8,0) => 0, (8,1) => 247962/1, (8,2) => 0, (10,0) => 0, (10,1) => 232050/1, (10,2) => 0, (12,0) => 0, (12,1) => 78540/1, (12,2) => 0, (14,0) => 0, (14,1) => 8364/1, (14,2) => 0, (16,0) => 0, (16,1) => 147/1, (16,2) => 0, (18,2) => 0, (1,0) => 85/1, (1,1) => 0, (3,0) => 1800/1, (3,1) => 238/1, (5,0) => 1287/1, (3,2) => 0, (5,1) => 28560/1};
--sb represents the betti numbers as sums of Schur functors
sb1234 = new HashTable from {(5,2) => {}, (7,0) => {}, (7,1) => {({21,4,19,15},1/1),({20,5,22,12},1/1),({20,5,21,13},2/1),({20,5,20,14},3/1),({20,5,19,15},4/1),({20,5,18,16},3/1),({20,5,17,17},1/1),({19,6,23,11},2/1),({19,6,22,12},5/1),({19,6,21,13},10/1),({19,6,20,14},13/1),({19,6,19,15},15/1),({19,6,18,16},11/1),({19,6,17,17},5/1),({18,7,25,9},1/1),({18,7,24,10},4/1),({18,7,23,11},10/1),({18,7,22,12},20/1),({18,7,21,13},31/1),({18,7,20,14},39/1),({18,7,19,15},40/1),({18,7,18,16},30/1),({18,7,17,17},11/1),({17,8,26,8},1/1),({17,8,25,9},4/1),({17,8,24,10},12/1),({17,8,23,11},27/1),({17,8,22,12},46/1),({17,8,21,13},68/1),({17,8,20,14},81/1),({17,8,19,15},81/1),({17,8,18,16},58/1),({17,8,17,17},23/1),({16,9,27,7},1/1),({16,9,26,8},4/1),({16,9,25,9},12/1),({16,9,24,10},28/1),({16,9,23,11},53/1),({16,9,22,12},85/1),({16,9,21,13},116/1),({16,9,20,14},134/1),({16,9,19,15},128/1),({16,9,18,16},93/1),({16,9,17,17},34/1),({15,10,27,7},1/1),({15,10,26,8},6/1),({15,10,25,9},18/1),({15,10,24,10},40/1),({15,10,23,11},73/1),({15,10,22,12},113/1),({15,10,21,13},150/1),({15,10,20,14},168/1),({15,10,19,15},159/1),({15,10,18,16},113/1),({15,10,17,17},42/1),({14,11,28,6},1/1),({14,11,27,7},3/1),({14,11,26,8},9/1),({14,11,25,9},22/1),({14,11,24,10},44/1),({14,11,23,11},77/1),({14,11,22,12},115/1),({14,11,21,13},148/1),({14,11,20,14},164/1),({14,11,19,15},152/1),({14,11,18,16},108/1),({14,11,17,17},39/1),({13,12,27,7},2/1),({13,12,26,8},5/1),({13,12,25,9},14/1),({13,12,24,10},28/1),({13,12,23,11},49/1),({13,12,22,12},71/1),({13,12,21,13},92/1),({13,12,20,14},100/1),({13,12,19,15},94/1),({13,12,18,16},65/1),({13,12,17,17},24/1)}, (9,0) => {}, (7,2) => {}, (9,1) => {({24,7,23,19},1/1),({24,7,22,20},1/1),({23,8,26,16},2/1),({23,8,25,17},3/1),({23,8,24,18},5/1),({23,8,23,19},6/1),({23,8,22,20},5/1),({23,8,21,21},1/1),({22,9,28,14},1/1),({22,9,27,15},4/1),({22,9,26,16},9/1),({22,9,25,17},16/1),({22,9,24,18},22/1),({22,9,23,19},23/1),({22,9,22,20},18/1),({22,9,21,21},7/1),({21,10,29,13},2/1),({21,10,28,14},7/1),({21,10,27,15},17/1),({21,10,26,16},33/1),({21,10,25,17},50/1),({21,10,24,18},63/1),({21,10,23,19},63/1),({21,10,22,20},48/1),({21,10,21,21},17/1),({20,11,30,12},1/1),({20,11,29,13},7/1),({20,11,28,14},20/1),({20,11,27,15},43/1),({20,11,26,16},76/1),({20,11,25,17},109/1),({20,11,24,18},130/1),({20,11,23,19},128/1),({20,11,22,20},94/1),({20,11,21,21},34/1),({19,12,31,11},1/1),({19,12,30,12},6/1),({19,12,29,13},17/1),({19,12,28,14},43/1),({19,12,27,15},83/1),({19,12,26,16},134/1),({19,12,25,17},183/1),({19,12,24,18},214/1),({19,12,23,19},202/1),({19,12,22,20},149/1),({19,12,21,21},54/1),({18,13,31,11},2/1),({18,13,30,12},9/1),({18,13,29,13},27/1),({18,13,28,14},61/1),({18,13,27,15},113/1),({18,13,26,16},177/1),({18,13,25,17},236/1),({18,13,24,18},268/1),({18,13,23,19},253/1),({18,13,22,20},182/1),({18,13,21,21},66/1),({17,14,32,10},1/1),({17,14,31,11},4/1),({17,14,30,12},12/1),({17,14,29,13},32/1),({17,14,28,14},66/1),({17,14,27,15},116/1),({17,14,26,16},178/1),({17,14,25,17},231/1),({17,14,24,18},258/1),({17,14,23,19},242/1),({17,14,22,20},173/1),({17,14,21,21},61/1),({16,15,31,11},2/1),({16,15,30,12},8/1),({16,15,29,13},20/1),({16,15,28,14},42/1),({16,15,27,15},74/1),({16,15,26,16},110/1),({16,15,25,17},143/1),({16,15,24,18},160/1),({16,15,23,19},147/1),({16,15,22,20},105/1),({16,15,21,21},39/1)}, (9,2) => {}, (11,0) => {}, (11,1) => {({26,11,29,21},1/1),({26,11,28,22},2/1),({26,11,27,23},2/1),({26,11,26,24},2/1),({26,11,25,25},1/1),({25,12,31,19},1/1),({25,12,30,20},3/1),({25,12,29,21},7/1),({25,12,28,22},10/1),({25,12,27,23},11/1),({25,12,26,24},9/1),({25,12,25,25},4/1),({24,13,32,18},2/1),({24,13,31,19},7/1),({24,13,30,20},15/1),({24,13,29,21},25/1),({24,13,28,22},34/1),({24,13,27,23},35/1),({24,13,26,24},27/1),({24,13,25,25},11/1),({23,14,33,17},2/1),({23,14,32,18},8/1),({23,14,31,19},21/1),({23,14,30,20},39/1),({23,14,29,21},61/1),({23,14,28,22},76/1),({23,14,27,23},77/1),({23,14,26,24},57/1),({23,14,25,25},22/1),({22,15,34,16},2/1),({22,15,33,17},7/1),({22,15,32,18},20/1),({22,15,31,19},44/1),({22,15,30,20},76/1),({22,15,29,21},110/1),({22,15,28,22},132/1),({22,15,27,23},129/1),({22,15,26,24},95/1),({22,15,25,25},35/1),({21,16,34,16},3/1),({21,16,33,17},12/1),({21,16,32,18},31/1),({21,16,31,19},63/1),({21,16,30,20},105/1),({21,16,29,21},147/1),({21,16,28,22},172/1),({21,16,27,23},166/1),({21,16,26,24},121/1),({21,16,25,25},45/1),({20,17,35,15},1/1),({20,17,34,16},5/1),({20,17,33,17},15/1),({20,17,32,18},36/1),({20,17,31,19},68/1),({20,17,30,20},109/1),({20,17,29,21},148/1),({20,17,28,22},170/1),({20,17,27,23},162/1),({20,17,26,24},117/1),({20,17,25,25},43/1),({19,18,35,15},1/1),({19,18,34,16},3/1),({19,18,33,17},10/1),({19,18,32,18},23/1),({19,18,31,19},44/1),({19,18,30,20},69/1),({19,18,29,21},93/1),({19,18,28,22},106/1),({19,18,27,23},101/1),({19,18,26,24},72/1),({19,18,25,25},27/1)}, (11,2) => {}, (13,0) => {}, (13,1) => {({28,15,32,26},1/1),({28,15,31,27},1/1),({28,15,30,28},1/1),({28,15,29,29},1/1),({27,16,34,24},1/1),({27,16,33,25},3/1),({27,16,32,26},5/1),({27,16,31,27},6/1),({27,16,30,28},5/1),({27,16,29,29},2/1),({26,17,35,23},2/1),({26,17,34,24},5/1),({26,17,33,25},10/1),({26,17,32,26},15/1),({26,17,31,27},16/1),({26,17,30,28},13/1),({26,17,29,29},5/1),({25,18,36,22},2/1),({25,18,35,23},6/1),({25,18,34,24},14/1),({25,18,33,25},23/1),({25,18,32,26},31/1),({25,18,31,27},32/1),({25,18,30,28},25/1),({25,18,29,29},9/1),({24,19,37,21},1/1),({24,19,36,22},4/1),({24,19,35,23},11/1),({24,19,34,24},22/1),({24,19,33,25},35/1),({24,19,32,26},44/1),({24,19,31,27},45/1),({24,19,30,28},34/1),({24,19,29,29},12/1),({23,20,37,21},2/1),({23,20,36,22},6/1),({23,20,35,23},14/1),({23,20,34,24},26/1),({23,20,33,25},38/1),({23,20,32,26},47/1),({23,20,31,27},47/1),({23,20,30,28},35/1),({23,20,29,29},13/1),({22,21,37,21},1/1),({22,21,36,22},4/1),({22,21,35,23},9/1),({22,21,34,24},17/1),({22,21,33,25},25/1),({22,21,32,26},30/1),({22,21,31,27},30/1),({22,21,30,28},22/1),({22,21,29,29},8/1)}, (15,0) => {}, (13,2) => {}, (15,1) => {({29,20,35,31},1/1),({29,20,34,32},1/1),({28,21,37,29},1/1),({28,21,36,30},2/1),({28,21,35,31},2/1),({28,21,34,32},2/1),({28,21,33,33},1/1),({27,22,38,28},1/1),({27,22,37,29},2/1),({27,22,36,30},3/1),({27,22,35,31},4/1),({27,22,34,32},3/1),({27,22,33,33},1/1),({26,23,38,28},2/1),({26,23,37,29},3/1),({26,23,36,30},4/1),({26,23,35,31},6/1),({26,23,34,32},4/1),({26,23,33,33},1/1),({25,24,39,27},1/1),({25,24,38,28},1/1),({25,24,37,29},2/1),({25,24,36,30},3/1),({25,24,35,31},3/1),({25,24,34,32},2/1),({25,24,33,33},1/1)}, (15,2) => {}, (17,0) => {}, (17,1) => {}, (17,2) => {({29,29,39,39},1/1)}, (0,0) => {({1,0,2,0},1/1)}, (2,0) => {({7,0,7,3},1/1),({7,0,5,5},1/1),({6,1,9,1},1/1),({6,1,8,2},2/1),({6,1,7,3},2/1),({6,1,6,4},2/1),({6,1,5,5},1/1),({5,2,9,1},1/1),({5,2,8,2},2/1),({5,2,7,3},2/1),({5,2,6,4},2/1),({5,2,5,5},1/1),({4,3,9,1},1/1),({4,3,8,2},2/1),({4,3,7,3},2/1),({4,3,6,4},2/1),({4,3,5,5},1/1)}, (2,1) => {({5,5,14,0},1/1)}, (4,0) => {({11,2,14,4},1/1),({11,2,13,5},1/1),({11,2,12,6},2/1),({11,2,11,7},1/1),({11,2,10,8},2/1),({10,3,14,4},1/1),({10,3,13,5},2/1),({10,3,12,6},3/1),({10,3,11,7},4/1),({10,3,10,8},3/1),({10,3,9,9},1/1),({9,4,15,3},1/1),({9,4,14,4},2/1),({9,4,13,5},4/1),({9,4,12,6},7/1),({9,4,11,7},7/1),({9,4,10,8},6/1),({9,4,9,9},2/1),({8,5,15,3},1/1),({8,5,14,4},2/1),({8,5,13,5},4/1),({8,5,12,6},7/1),({8,5,11,7},7/1),({8,5,10,8},6/1),({8,5,9,9},3/1),({7,6,14,4},2/1),({7,6,13,5},3/1),({7,6,12,6},5/1),({7,6,11,7},5/1),({7,6,10,8},5/1),({7,6,9,9},1/1)}, (2,2) => {}, (4,1) => {({15,1,13,9},1/1),({14,2,13,9},1/1),({14,2,12,10},1/1),({13,3,16,6},1/1),({13,3,15,7},1/1),({13,3,14,8},1/1),({13,3,13,9},2/1),({13,3,12,10},1/1),({13,3,11,11},1/1),({12,4,16,6},1/1),({12,4,15,7},2/1),({12,4,14,8},2/1),({12,4,13,9},2/1),({12,4,12,10},2/1),({11,5,19,3},1/1),({11,5,18,4},1/1),({11,5,17,5},2/1),({11,5,16,6},2/1),({11,5,15,7},3/1),({11,5,14,8},3/1),({11,5,13,9},3/1),({11,5,12,10},1/1),({11,5,11,11},1/1),({10,6,20,2},1/1),({10,6,19,3},1/1),({10,6,18,4},3/1),({10,6,17,5},3/1),({10,6,16,6},4/1),({10,6,15,7},3/1),({10,6,14,8},4/1),({10,6,13,9},2/1),({10,6,12,10},2/1),({9,7,19,3},1/1),({9,7,18,4},1/1),({9,7,17,5},3/1),({9,7,16,6},2/1),({9,7,15,7},3/1),({9,7,14,8},1/1),({9,7,13,9},2/1),({9,7,12,10},1/1),({9,7,11,11},1/1),({8,8,20,2},1/1),({8,8,19,3},1/1),({8,8,18,4},2/1),({8,8,17,5},2/1),({8,8,16,6},3/1),({8,8,15,7},1/1),({8,8,14,8},2/1),({8,8,12,10},1/1)}, (4,2) => {}, (6,0) => {}, (6,1) => {({19,3,18,12},1/1),({19,3,16,14},1/1),({18,4,20,10},1/1),({18,4,19,11},2/1),({18,4,18,12},3/1),({18,4,17,13},4/1),({18,4,16,14},3/1),({18,4,15,15},1/1),({17,5,21,9},1/1),({17,5,20,10},5/1),({17,5,19,11},8/1),({17,5,18,12},12/1),({17,5,17,13},12/1),({17,5,16,14},11/1),({17,5,15,15},2/1),({16,6,23,7},1/1),({16,6,22,8},3/1),({16,6,21,9},8/1),({16,6,20,10},14/1),({16,6,19,11},23/1),({16,6,18,12},29/1),({16,6,17,13},30/1),({16,6,16,14},22/1),({16,6,15,15},9/1),({15,7,24,6},1/1),({15,7,23,7},3/1),({15,7,22,8},10/1),({15,7,21,9},18/1),({15,7,20,10},34/1),({15,7,19,11},45/1),({15,7,18,12},57/1),({15,7,17,13},52/1),({15,7,16,14},42/1),({15,7,15,15},12/1),({14,8,25,5},1/1),({14,8,24,6},2/1),({14,8,23,7},8/1),({14,8,22,8},17/1),({14,8,21,9},33/1),({14,8,20,10},52/1),({14,8,19,11},71/1),({14,8,18,12},79/1),({14,8,17,13},78/1),({14,8,16,14},54/1),({14,8,15,15},21/1),({13,9,25,5},1/1),({13,9,24,6},5/1),({13,9,23,7},11/1),({13,9,22,8},25/1),({13,9,21,9},42/1),({13,9,20,10},65/1),({13,9,19,11},81/1),({13,9,18,12},95/1),({13,9,17,13},82/1),({13,9,16,14},64/1),({13,9,15,15},20/1),({12,10,25,5},1/1),({12,10,24,6},3/1),({12,10,23,7},10/1),({12,10,22,8},19/1),({12,10,21,9},36/1),({12,10,20,10},50/1),({12,10,19,11},68/1),({12,10,18,12},71/1),({12,10,17,13},69/1),({12,10,16,14},45/1),({12,10,15,15},19/1),({11,11,26,4},1/1),({11,11,25,5},1/1),({11,11,24,6},3/1),({11,11,23,7},5/1),({11,11,22,8},11/1),({11,11,21,9},14/1),({11,11,20,10},25/1),({11,11,19,11},26/1),({11,11,18,12},32/1),({11,11,17,13},25/1),({11,11,16,14},22/1),({11,11,15,15},3/1)}, (6,2) => {}, (8,0) => {}, (8,1) => {({23,5,19,19},1/1),({22,6,23,15},1/1),({22,6,22,16},2/1),({22,6,21,17},2/1),({22,6,20,18},2/1),({22,6,19,19},1/1),({21,7,25,13},2/1),({21,7,24,14},3/1),({21,7,23,15},8/1),({21,7,22,16},10/1),({21,7,21,17},12/1),({21,7,20,18},8/1),({21,7,19,19},5/1),({20,8,26,12},3/1),({20,8,25,13},8/1),({20,8,24,14},17/1),({20,8,23,15},27/1),({20,8,22,16},36/1),({20,8,21,17},36/1),({20,8,20,18},28/1),({20,8,19,19},10/1),({19,9,28,10},1/1),({19,9,27,11},4/1),({19,9,26,12},11/1),({19,9,25,13},27/1),({19,9,24,14},46/1),({19,9,23,15},71/1),({19,9,22,16},83/1),({19,9,21,17},86/1),({19,9,20,18},60/1),({19,9,19,19},26/1),({18,10,28,10},3/1),({18,10,27,11},11/1),({18,10,26,12},29/1),({18,10,25,13},57/1),({18,10,24,14},96/1),({18,10,23,15},131/1),({18,10,22,16},156/1),({18,10,21,17},148/1),({18,10,20,18},110/1),({18,10,19,19},39/1),({17,11,29,9},2/1),({17,11,28,10},7/1),({17,11,27,11},22/1),({17,11,26,12},48/1),({17,11,25,13},93/1),({17,11,24,14},142/1),({17,11,23,15},196/1),({17,11,22,16},218/1),({17,11,21,17},212/1),({17,11,20,18},147/1),({17,11,19,19},59/1),({16,12,30,8},1/1),({16,12,29,9},3/1),({16,12,28,10},12/1),({16,12,27,11},29/1),({16,12,26,12},64/1),({16,12,25,13},110/1),({16,12,24,14},170/1),({16,12,23,15},218/1),({16,12,22,16},249/1),({16,12,21,17},227/1),({16,12,20,18},167/1),({16,12,19,19},57/1),({15,13,29,9},3/1),({15,13,28,10},9/1),({15,13,27,11},26/1),({15,13,26,12},50/1),({15,13,25,13},93/1),({15,13,24,14},133/1),({15,13,23,15},179/1),({15,13,22,16},191/1),({15,13,21,17},184/1),({15,13,20,18},124/1),({15,13,19,19},52/1),({14,14,30,8},1/1),({14,14,29,9},2/1),({14,14,28,10},6/1),({14,14,27,11},11/1),({14,14,26,12},24/1),({14,14,25,13},36/1),({14,14,24,14},57/1),({14,14,23,15},68/1),({14,14,22,16},80/1),({14,14,21,17},67/1),({14,14,20,18},54/1),({14,14,19,19},14/1)}, (8,2) => {}, (10,0) => {}, (10,1) => {({25,9,26,20},2/1),({25,9,25,21},1/1),({25,9,24,22},2/1),({24,10,29,17},1/1),({24,10,28,18},3/1),({24,10,27,19},6/1),({24,10,26,20},8/1),({24,10,25,21},10/1),({24,10,24,22},8/1),({24,10,23,23},3/1),({23,11,30,16},2/1),({23,11,29,17},6/1),({23,11,28,18},15/1),({23,11,27,19},23/1),({23,11,26,20},33/1),({23,11,25,21},32/1),({23,11,24,22},27/1),({23,11,23,23},8/1),({22,12,31,15},3/1),({22,12,30,16},9/1),({22,12,29,17},23/1),({22,12,28,18},42/1),({22,12,27,19},65/1),({22,12,26,20},79/1),({22,12,25,21},81/1),({22,12,24,22},59/1),({22,12,23,23},23/1),({21,13,32,14},2/1),({21,13,31,15},8/1),({21,13,30,16},25/1),({21,13,29,17},50/1),({21,13,28,18},90/1),({21,13,27,19},124/1),({21,13,26,20},153/1),({21,13,25,21},144/1),({21,13,24,22},111/1),({21,13,23,23},36/1),({20,14,33,13},1/1),({20,14,32,14},5/1),({20,14,31,15},18/1),({20,14,30,16},42/1),({20,14,29,17},85/1),({20,14,28,18},136/1),({20,14,27,19},190/1),({20,14,26,20},217/1),({20,14,25,21},211/1),({20,14,24,22},150/1),({20,14,23,23},58/1),({19,15,33,13},2/1),({19,15,32,14},9/1),({19,15,31,15},24/1),({19,15,30,16},57/1),({19,15,29,17},102/1),({19,15,28,18},164/1),({19,15,27,19},214/1),({19,15,26,20},250/1),({19,15,25,21},229/1),({19,15,24,22},172/1),({19,15,23,23},57/1),({18,16,33,13},2/1),({18,16,32,14},7/1),({18,16,31,15},22/1),({18,16,30,16},46/1),({18,16,29,17},87/1),({18,16,28,18},130/1),({18,16,27,19},177/1),({18,16,26,20},194/1),({18,16,25,21},188/1),({18,16,24,22},129/1),({18,16,23,23},52/1),({17,17,34,12},1/1),({17,17,33,13},1/1),({17,17,32,14},5/1),({17,17,31,15},9/1),({17,17,30,16},22/1),({17,17,29,17},34/1),({17,17,28,18},56/1),({17,17,27,19},66/1),({17,17,26,20},82/1),({17,17,25,21},68/1),({17,17,24,22},56/1),({17,17,23,23},14/1)}, (10,2) => {}, (12,0) => {}, (12,1) => {({27,13,31,23},1/1),({27,13,30,24},1/1),({27,13,29,25},3/1),({27,13,28,26},1/1),({27,13,27,27},1/1),({26,14,32,22},3/1),({26,14,31,23},5/1),({26,14,30,24},8/1),({26,14,29,25},10/1),({26,14,28,26},8/1),({26,14,27,27},2/1),({25,15,34,20},1/1),({25,15,33,21},5/1),({25,15,32,22},10/1),({25,15,31,23},20/1),({25,15,30,24},25/1),({25,15,29,25},29/1),({25,15,28,26},20/1),({25,15,27,27},10/1),({24,16,35,19},1/1),({24,16,34,20},5/1),({24,16,33,21},13/1),({24,16,32,22},28/1),({24,16,31,23},43/1),({24,16,30,24},56/1),({24,16,29,25},56/1),({24,16,28,26},44/1),({24,16,27,27},15/1),({23,17,35,19},4/1),({23,17,34,20},11/1),({23,17,33,21},27/1),({23,17,32,22},47/1),({23,17,31,23},73/1),({23,17,30,24},85/1),({23,17,29,25},90/1),({23,17,28,26},62/1),({23,17,27,27},26/1),({22,18,36,18},2/1),({22,18,35,19},5/1),({22,18,34,20},17/1),({22,18,33,21},35/1),({22,18,32,22},62/1),({22,18,31,23},86/1),({22,18,30,24},107/1),({22,18,29,25},99/1),({22,18,28,26},77/1),({22,18,27,27},26/1),({21,19,36,18},1/1),({21,19,35,19},6/1),({21,19,34,20},14/1),({21,19,33,21},32/1),({21,19,32,22},50/1),({21,19,31,23},75/1),({21,19,30,24},83/1),({21,19,29,25},86/1),({21,19,28,26},57/1),({21,19,27,27},26/1),({20,20,36,18},1/1),({20,20,35,19},3/1),({20,20,34,20},8/1),({20,20,33,21},12/1),({20,20,32,22},24/1),({20,20,31,23},28/1),({20,20,30,24},36/1),({20,20,29,25},31/1),({20,20,28,26},27/1),({20,20,27,27},5/1)}, (12,2) => {}, (14,0) => {}, (14,1) => {({29,17,32,30},1/1),({28,18,35,27},1/1),({28,18,34,28},2/1),({28,18,33,29},3/1),({28,18,32,30},2/1),({28,18,31,31},1/1),({27,19,36,26},2/1),({27,19,35,27},3/1),({27,19,34,28},7/1),({27,19,33,29},6/1),({27,19,32,30},7/1),({27,19,31,31},1/1),({26,20,37,25},2/1),({26,20,36,26},4/1),({26,20,35,27},9/1),({26,20,34,28},12/1),({26,20,33,29},14/1),({26,20,32,30},10/1),({26,20,31,31},5/1),({25,21,38,24},1/1),({25,21,37,25},3/1),({25,21,36,26},8/1),({25,21,35,27},12/1),({25,21,34,28},18/1),({25,21,33,29},16/1),({25,21,32,30},15/1),({25,21,31,31},4/1),({24,22,38,24},1/1),({24,22,37,25},4/1),({24,22,36,26},6/1),({24,22,35,27},12/1),({24,22,34,28},14/1),({24,22,33,29},16/1),({24,22,32,30},10/1),({24,22,31,31},6/1),({23,23,38,24},1/1),({23,23,37,25},1/1),({23,23,36,26},4/1),({23,23,35,27},4/1),({23,23,34,28},8/1),({23,23,33,29},5/1),({23,23,32,30},7/1)}, (14,2) => {}, (16,0) => {}, (16,1) => {({29,23,37,33},1/1),({29,23,35,35},1/1),({28,24,38,32},1/1),({28,24,36,34},1/1),({27,25,39,31},1/1),({27,25,37,33},1/1),({27,25,35,35},1/1),({26,26,38,32},1/1),({26,26,36,34},1/1)}, (16,2) => {}, (18,2) => {}, (1,0) => {({4,0,5,1},1/1),({4,0,4,2},1/1),({3,1,6,0},1/1),({3,1,5,1},1/1),({3,1,4,2},1/1)}, (1,1) => {}, (3,0) => {({9,1,11,3},1/1),({9,1,10,4},1/1),({9,1,9,5},2/1),({9,1,8,6},1/1),({9,1,7,7},1/1),({8,2,12,2},1/1),({8,2,11,3},2/1),({8,2,10,4},3/1),({8,2,9,5},4/1),({8,2,8,6},3/1),({8,2,7,7},1/1),({7,3,12,2},1/1),({7,3,11,3},3/1),({7,3,10,4},4/1),({7,3,9,5},6/1),({7,3,8,6},4/1),({7,3,7,7},2/1),({6,4,12,2},2/1),({6,4,11,3},3/1),({6,4,10,4},5/1),({6,4,9,5},6/1),({6,4,8,6},5/1),({6,4,7,7},1/1),({5,5,11,3},1/1),({5,5,10,4},1/1),({5,5,9,5},2/1),({5,5,8,6},1/1),({5,5,7,7},1/1)}, (3,1) => {({13,0,9,9},1/1),({8,5,17,1},1/1),({8,5,16,2},1/1),({8,5,15,3},1/1),({8,5,14,4},1/1)}, (5,0) => {({13,3,16,6},1/1),({13,3,14,8},1/1),({13,3,12,10},1/1),({12,4,15,7},1/1),({12,4,14,8},1/1),({12,4,13,9},1/1),({12,4,12,10},1/1),({11,5,16,6},1/1),({11,5,15,7},1/1),({11,5,14,8},3/1),({11,5,13,9},2/1),({11,5,12,10},3/1),({10,6,15,7},1/1),({10,6,14,8},2/1),({10,6,13,9},3/1),({10,6,12,10},2/1),({10,6,11,11},1/1),({9,7,16,6},1/1),({9,7,15,7},1/1),({9,7,14,8},3/1),({9,7,13,9},2/1),({9,7,12,10},4/1),({8,8,13,9},1/1),({8,8,11,11},1/1)}, (3,2) => {}, (5,1) => {({17,2,16,10},1/1),({17,2,14,12},1/1),({16,3,17,9},1/1),({16,3,16,10},2/1),({16,3,15,11},2/1),({16,3,14,12},2/1),({16,3,13,13},1/1),({15,4,19,7},1/1),({15,4,18,8},2/1),({15,4,17,9},4/1),({15,4,16,10},7/1),({15,4,15,11},6/1),({15,4,14,12},6/1),({15,4,13,13},2/1),({14,5,20,6},1/1),({14,5,19,7},3/1),({14,5,18,8},6/1),({14,5,17,9},10/1),({14,5,16,10},13/1),({14,5,15,11},13/1),({14,5,14,12},10/1),({14,5,13,13},4/1),({13,6,22,4},1/1),({13,6,21,5},2/1),({13,6,20,6},5/1),({13,6,19,7},9/1),({13,6,18,8},15/1),({13,6,17,9},19/1),({13,6,16,10},24/1),({13,6,15,11},21/1),({13,6,14,12},17/1),({13,6,13,13},5/1),({12,7,22,4},1/1),({12,7,21,5},3/1),({12,7,20,6},7/1),({12,7,19,7},13/1),({12,7,18,8},20/1),({12,7,17,9},26/1),({12,7,16,10},29/1),({12,7,15,11},27/1),({12,7,14,12},19/1),({12,7,13,13},7/1),({11,8,23,3},1/1),({11,8,22,4},2/1),({11,8,21,5},5/1),({11,8,20,6},10/1),({11,8,19,7},16/1),({11,8,18,8},23/1),({11,8,17,9},28/1),({11,8,16,10},30/1),({11,8,15,11},26/1),({11,8,14,12},19/1),({11,8,13,13},6/1),({10,9,22,4},1/1),({10,9,21,5},3/1),({10,9,20,6},6/1),({10,9,19,7},10/1),({10,9,18,8},14/1),({10,9,17,9},17/1),({10,9,16,10},18/1),({10,9,15,11},16/1),({10,9,14,12},11/1),({10,9,13,13},4/1)}};
--dw stands for dominant weights
<<<<<<< HEAD
dw1234 = new HashTable from {(5,2) => {}, (7,0) => {}, (7,1) => {({21,4,19,15},1/1),({20,5,22,12},1/1),({19,6,23,11},2/1),({18,7,25,9},1/1),({17,8,26,8},1/1),({16,9,27,7},1/1),({14,11,28,6},1/1)}, (9,0) => {}, (7,2) => {}, (9,1) => {({24,7,23,19},1/1),({23,8,26,16},2/1),({22,9,28,14},1/1),({21,10,29,13},2/1),({20,11,30,12},1/1),({19,12,31,11},1/1),({17,14,32,10},1/1)}, (9,2) => {}, (11,0) => {}, (11,1) => {({26,11,29,21},1/1),({25,12,31,19},1/1),({24,13,32,18},2/1),({23,14,33,17},2/1),({22,15,34,16},2/1),({20,17,35,15},1/1)}, (11,2) => {}, (13,0) => {}, (13,1) => {({28,15,32,26},1/1),({27,16,34,24},1/1),({26,17,35,23},2/1),({25,18,36,22},2/1),({24,19,37,21},1/1)}, (15,0) => {}, (13,2) => {}, (15,1) => {({29,20,35,31},1/1),({28,21,37,29},1/1),({27,22,38,28},1/1),({25,24,39,27},1/1)}, (15,2) => {}, (17,0) => {}, (17,1) => {}, (17,2) => {({29,29,39,39},1/1)}, (0,0) => {({1,0,2,0},1/1)}, (2,0) => {({7,0,7,3},1/1),({6,1,9,1},1/1)}, (2,1) => {({5,5,14,0},1/1)}, (4,0) => {({11,2,14,4},1/1),({9,4,15,3},1/1)}, (2,2) => {}, (4,1) => {({15,1,13,9},1/1),({13,3,16,6},1/1),({11,5,19,3},1/1),({10,6,20,2},1/1)}, (4,2) => {}, (6,0) => {}, (6,1) => {({19,3,18,12},1/1),({18,4,20,10},1/1),({17,5,21,9},1/1),({16,6,23,7},1/1),({15,7,24,6},1/1),({14,8,25,5},1/1),({11,11,26,4},1/1)}, (6,2) => {}, (8,0) => {}, (8,1) => {({23,5,19,19},1/1),({22,6,23,15},1/1),({21,7,25,13},2/1),({20,8,26,12},3/1),({19,9,28,10},1/1),({17,11,29,9},2/1),({16,12,30,8},1/1)}, (8,2) => {}, (10,0) => {}, (10,1) => {({25,9,26,20},2/1),({24,10,29,17},1/1),({23,11,30,16},2/1),({22,12,31,15},3/1),({21,13,32,14},2/1),({20,14,33,13},1/1),({17,17,34,12},1/1)}, (10,2) => {}, (12,0) => {}, (12,1) => {({27,13,31,23},1/1),({26,14,32,22},3/1),({25,15,34,20},1/1),({24,16,35,19},1/1),({22,18,36,18},2/1)}, (12,2) => {}, (14,0) => {}, (14,1) => {({29,17,32,30},1/1),({28,18,35,27},1/1),({27,19,36,26},2/1),({26,20,37,25},2/1),({25,21,38,24},1/1)}, (14,2) => {}, (16,0) => {}, (16,1) => {({29,23,37,33},1/1),({28,24,38,32},1/1),({27,25,39,31},1/1)}, (16,2) => {}, (18,2) => {}, (1,0) => {({4,0,5,1},1/1),({3,1,6,0},1/1)}, (1,1) => {}, (3,0) => {({9,1,11,3},1/1),({8,2,12,2},1/1)}, (3,1) => {({13,0,9,9},1/1),({8,5,17,1},1/1)}, (5,0) => {({13,3,16,6},1/1)}, (3,2) => {}, (5,1) => {({17,2,16,10},1/1),({16,3,17,9},1/1),({15,4,19,7},1/1),({14,5,20,6},1/1),({13,6,22,4},1/1),({11,8,23,3},1/1)}};
--dmw stands for dominant monomial weights
dmw1234 = {infinity};
=======
dw1234 = new HashTable from {(5,2) => {}, (7,0) => {}, (7,1) => {({21,4,19,15},1/1),({20,5,22,12},1/1),({19,6,23,11},2/1),({18,7,25,9},1/1),({17,8,26,8},1/1),({16,9,27,7},1/1),({14,11,28,6},1/1)}, (7,2) => {}, (9,0) => {}, (9,1) => {({24,7,23,19},1/1),({23,8,26,16},2/1),({22,9,28,14},1/1),({21,10,29,13},2/1),({20,11,30,12},1/1),({19,12,31,11},1/1),({17,14,32,10},1/1)}, (9,2) => {}, (11,0) => {}, (11,1) => {({26,11,29,21},1/1),({25,12,31,19},1/1),({24,13,32,18},2/1),({23,14,33,17},2/1),({22,15,34,16},2/1),({20,17,35,15},1/1)}, (11,2) => {}, (13,0) => {}, (13,1) => {({28,15,32,26},1/1),({27,16,34,24},1/1),({26,17,35,23},2/1),({25,18,36,22},2/1),({24,19,37,21},1/1)}, (13,2) => {}, (15,0) => {}, (15,1) => {({29,20,35,31},1/1),({28,21,37,29},1/1),({27,22,38,28},1/1),({25,24,39,27},1/1)}, (17,0) => {}, (15,2) => {}, (17,1) => {}, (17,2) => {({29,29,39,39},1/1)}, (0,0) => {({1,0,2,0},1/1)}, (2,0) => {({7,0,7,3},1/1),({6,1,9,1},1/1)}, (2,1) => {({5,5,14,0},1/1)}, (2,2) => {}, (4,0) => {({11,2,14,4},1/1),({9,4,15,3},1/1)}, (4,1) => {({15,1,13,9},1/1),({13,3,16,6},1/1),({11,5,19,3},1/1),({10,6,20,2},1/1)}, (4,2) => {}, (6,0) => {}, (6,1) => {({19,3,18,12},1/1),({18,4,20,10},1/1),({17,5,21,9},1/1),({16,6,23,7},1/1),({15,7,24,6},1/1),({14,8,25,5},1/1),({11,11,26,4},1/1)}, (6,2) => {}, (8,0) => {}, (8,1) => {({23,5,19,19},1/1),({22,6,23,15},1/1),({21,7,25,13},2/1),({20,8,26,12},3/1),({19,9,28,10},1/1),({17,11,29,9},2/1),({16,12,30,8},1/1)}, (8,2) => {}, (10,0) => {}, (10,1) => {({25,9,26,20},2/1),({24,10,29,17},1/1),({23,11,30,16},2/1),({22,12,31,15},3/1),({21,13,32,14},2/1),({20,14,33,13},1/1),({17,17,34,12},1/1)}, (10,2) => {}, (12,0) => {}, (12,1) => {({27,13,31,23},1/1),({26,14,32,22},3/1),({25,15,34,20},1/1),({24,16,35,19},1/1),({22,18,36,18},2/1)}, (12,2) => {}, (14,0) => {}, (14,1) => {({29,17,32,30},1/1),({28,18,35,27},1/1),({27,19,36,26},2/1),({26,20,37,25},2/1),({25,21,38,24},1/1)}, (14,2) => {}, (16,0) => {}, (16,1) => {({29,23,37,33},1/1),({28,24,38,32},1/1),({27,25,39,31},1/1)}, (16,2) => {}, (18,2) => {}, (1,0) => {({4,0,5,1},1/1),({3,1,6,0},1/1)}, (1,1) => {}, (3,0) => {({9,1,11,3},1/1),({8,2,12,2},1/1)}, (3,1) => {({13,0,9,9},1/1),({8,5,17,1},1/1)}, (3,2) => {}, (5,0) => {({13,3,16,6},1/1)}, (5,1) => {({17,2,16,10},1/1),({16,3,17,9},1/1),({15,4,19,7},1/1),({14,5,20,6},1/1),({13,6,22,4},1/1),({11,8,23,3},1/1)}};
>>>>>>> parent of 67e86b7... updating data
end;