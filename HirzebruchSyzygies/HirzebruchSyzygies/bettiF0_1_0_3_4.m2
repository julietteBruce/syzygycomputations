A := QQ[t_0,t_1,t_2,t_3];
--tb stands for Total Betti numbers
tb1034 = new HashTable from {(5,2) => 0, (7,0) => 0, (7,1) => 283764/1, (9,0) => 0, (7,2) => 0, (9,1) => 350064/1, (9,2) => 0, (11,0) => 0, (11,1) => 164424/1, (11,2) => 0, (13,0) => 0, (13,1) => 25296/1, (15,0) => 0, (13,2) => 0, (15,1) => 510/1, (15,2) => 120/1, (17,0) => 0, (17,1) => 0, (17,2) => 3/1, (0,0) => 2/1, (2,0) => 0, (2,1) => 1173/1, (4,0) => 0, (2,2) => 0, (4,1) => 30804/1, (4,2) => 0, (6,0) => 0, (6,1) => 177684/1, (6,2) => 0, (8,0) => 0, (8,1) => 354926/1, (8,2) => 0, (10,0) => 0, (10,1) => 271830/1, (10,2) => 0, (12,0) => 0, (12,1) => 75684/1, (12,2) => 0, (14,0) => 0, (14,1) => 5508/1, (14,2) => 0, (16,0) => 0, (16,1) => 33/1, (16,2) => 32/1, (18,2) => 0, (1,0) => 15/1, (1,1) => 48/1, (3,0) => 0, (3,1) => 7752/1, (5,0) => 0, (3,2) => 0, (5,1) => 85680/1};
--sb represents the betti numbers as sums of Schur functors
sb1034 = new HashTable from {(5,2) => {}, (7,0) => {}, (7,1) => {({21,4,19,13},1/1),({21,4,18,14},1/1),({21,4,17,15},1/1),({20,5,22,10},1/1),({20,5,21,11},2/1),({20,5,20,12},4/1),({20,5,19,13},6/1),({20,5,18,14},7/1),({20,5,17,15},6/1),({20,5,16,16},2/1),({19,6,23,9},1/1),({19,6,22,10},4/1),({19,6,21,11},11/1),({19,6,20,12},18/1),({19,6,19,13},25/1),({19,6,18,14},26/1),({19,6,17,15},21/1),({19,6,16,16},7/1),({18,7,24,8},2/1),({18,7,23,9},7/1),({18,7,22,10},18/1),({18,7,21,11},35/1),({18,7,20,12},54/1),({18,7,19,13},68/1),({18,7,18,14},69/1),({18,7,17,15},52/1),({18,7,16,16},19/1),({17,8,25,7},1/1),({17,8,24,8},7/1),({17,8,23,9},21/1),({17,8,22,10},45/1),({17,8,21,11},80/1),({17,8,20,12},115/1),({17,8,19,13},140/1),({17,8,18,14},136/1),({17,8,17,15},102/1),({17,8,16,16},36/1),({16,9,26,6},1/1),({16,9,25,7},5/1),({16,9,24,8},17/1),({16,9,23,9},42/1),({16,9,22,10},84/1),({16,9,21,11},138/1),({16,9,20,12},191/1),({16,9,19,13},223/1),({16,9,18,14},214/1),({16,9,17,15},156/1),({16,9,16,16},57/1),({15,10,26,6},2/1),({15,10,25,7},9/1),({15,10,24,8},27/1),({15,10,23,9},62/1),({15,10,22,10},116/1),({15,10,21,11},183/1),({15,10,20,12},245/1),({15,10,19,13},280/1),({15,10,18,14},264/1),({15,10,17,15},192/1),({15,10,16,16},69/1),({14,11,26,6},3/1),({14,11,25,7},11/1),({14,11,24,8},30/1),({14,11,23,9},65/1),({14,11,22,10},117/1),({14,11,21,11},180/1),({14,11,20,12},237/1),({14,11,19,13},267/1),({14,11,18,14},250/1),({14,11,17,15},179/1),({14,11,16,16},65/1),({13,12,27,5},1/1),({13,12,26,6},2/1),({13,12,25,7},8/1),({13,12,24,8},20/1),({13,12,23,9},43/1),({13,12,22,10},75/1),({13,12,21,11},114/1),({13,12,20,12},147/1),({13,12,19,13},166/1),({13,12,18,14},153/1),({13,12,17,15},110/1),({13,12,16,16},39/1)}, (9,0) => {}, (7,2) => {}, (9,1) => {({24,7,23,17},1/1),({24,7,22,18},1/1),({24,7,21,19},1/1),({24,7,20,20},1/1),({23,8,26,14},1/1),({23,8,25,15},2/1),({23,8,24,16},5/1),({23,8,23,17},7/1),({23,8,22,18},8/1),({23,8,21,19},6/1),({23,8,20,20},3/1),({22,9,27,13},2/1),({22,9,26,14},6/1),({22,9,25,15},13/1),({22,9,24,16},22/1),({22,9,23,17},29/1),({22,9,22,18},30/1),({22,9,21,19},23/1),({22,9,20,20},9/1),({21,10,29,11},1/1),({21,10,28,12},3/1),({21,10,27,13},10/1),({21,10,26,14},24/1),({21,10,25,15},43/1),({21,10,24,16},65/1),({21,10,23,17},80/1),({21,10,22,18},80/1),({21,10,21,19},59/1),({21,10,20,20},23/1),({20,11,29,11},3/1),({20,11,28,12},11/1),({20,11,27,13},29/1),({20,11,26,14},59/1),({20,11,25,15},99/1),({20,11,24,16},139/1),({20,11,23,17},164/1),({20,11,22,18},159/1),({20,11,21,19},116/1),({20,11,20,20},43/1),({19,12,30,10},2/1),({19,12,29,11},8/1),({19,12,28,12},25/1),({19,12,27,13},57/1),({19,12,26,14},108/1),({19,12,25,15},170/1),({19,12,24,16},230/1),({19,12,23,17},262/1),({19,12,22,18},250/1),({19,12,21,19},179/1),({19,12,20,20},67/1),({18,13,31,9},1/1),({18,13,30,10},4/1),({18,13,29,11},15/1),({18,13,28,12},39/1),({18,13,27,13},83/1),({18,13,26,14},148/1),({18,13,25,15},225/1),({18,13,24,16},295/1),({18,13,23,17},331/1),({18,13,22,18},309/1),({18,13,21,19},221/1),({18,13,20,20},81/1),({17,14,31,9},1/1),({17,14,30,10},5/1),({17,14,29,11},17/1),({17,14,28,12},42/1),({17,14,27,13},86/1),({17,14,26,14},149/1),({17,14,25,15},221/1),({17,14,24,16},285/1),({17,14,23,17},315/1),({17,14,22,18},292/1),({17,14,21,19},207/1),({17,14,20,20},76/1),({16,15,31,9},1/1),({16,15,30,10},4/1),({16,15,29,11},12/1),({16,15,28,12},28/1),({16,15,27,13},56/1),({16,15,26,14},95/1),({16,15,25,15},139/1),({16,15,24,16},178/1),({16,15,23,17},195/1),({16,15,22,18},179/1),({16,15,21,19},127/1),({16,15,20,20},46/1)}, (9,2) => {}, (11,0) => {}, (11,1) => {({26,11,28,20},1/1),({26,11,27,21},2/1),({26,11,26,22},2/1),({26,11,25,23},2/1),({26,11,24,24},1/1),({25,12,30,18},1/1),({25,12,29,19},4/1),({25,12,28,20},7/1),({25,12,27,21},10/1),({25,12,26,22},11/1),({25,12,25,23},9/1),({25,12,24,24},3/1),({24,13,32,16},1/1),({24,13,31,17},3/1),({24,13,30,18},8/1),({24,13,29,19},16/1),({24,13,28,20},26/1),({24,13,27,21},33/1),({24,13,26,22},34/1),({24,13,25,23},26/1),({24,13,24,24},9/1),({23,14,33,15},1/1),({23,14,32,16},4/1),({23,14,31,17},12/1),({23,14,30,18},25/1),({23,14,29,19},43/1),({23,14,28,20},62/1),({23,14,27,21},76/1),({23,14,26,22},73/1),({23,14,25,23},55/1),({23,14,24,24},20/1),({22,15,34,14},1/1),({22,15,33,15},3/1),({22,15,32,16},11/1),({22,15,31,17},26/1),({22,15,30,18},50/1),({22,15,29,19},82/1),({22,15,28,20},111/1),({22,15,27,21},128/1),({22,15,26,22},123/1),({22,15,25,23},89/1),({22,15,24,24},32/1),({21,16,34,14},2/1),({21,16,33,15},7/1),({21,16,32,16},19/1),({21,16,31,17},41/1),({21,16,30,18},74/1),({21,16,29,19},114/1),({21,16,28,20},150/1),({21,16,27,21},169/1),({21,16,26,22},158/1),({21,16,25,23},114/1),({21,16,24,24},41/1),({20,17,34,14},2/1),({20,17,33,15},8/1),({20,17,32,16},21/1),({20,17,31,17},44/1),({20,17,30,18},78/1),({20,17,29,19},115/1),({20,17,28,20},149/1),({20,17,27,21},166/1),({20,17,26,22},152/1),({20,17,25,23},109/1),({20,17,24,24},40/1),({19,18,35,13},1/1),({19,18,34,14},2/1),({19,18,33,15},6/1),({19,18,32,16},15/1),({19,18,31,17},30/1),({19,18,30,18},50/1),({19,18,29,19},75/1),({19,18,28,20},95/1),({19,18,27,21},104/1),({19,18,26,22},96/1),({19,18,25,23},68/1),({19,18,24,24},23/1)}, (11,2) => {}, (13,0) => {}, (13,1) => {({28,15,30,26},1/1),({28,15,29,27},1/1),({27,16,33,23},1/1),({27,16,32,24},2/1),({27,16,31,25},3/1),({27,16,30,26},4/1),({27,16,29,27},3/1),({27,16,28,28},1/1),({26,17,35,21},1/1),({26,17,34,22},2/1),({26,17,33,23},5/1),({26,17,32,24},8/1),({26,17,31,25},10/1),({26,17,30,26},11/1),({26,17,29,27},8/1),({26,17,28,28},3/1),({25,18,36,20},1/1),({25,18,35,21},3/1),({25,18,34,22},7/1),({25,18,33,23},12/1),({25,18,32,24},18/1),({25,18,31,25},21/1),({25,18,30,26},21/1),({25,18,29,27},15/1),({25,18,28,28},6/1),({24,19,37,19},1/1),({24,19,36,20},3/1),({24,19,35,21},7/1),({24,19,34,22},13/1),({24,19,33,23},21/1),({24,19,32,24},28/1),({24,19,31,25},32/1),({24,19,30,26},30/1),({24,19,29,27},21/1),({24,19,28,28},8/1),({23,20,37,19},1/1),({23,20,36,20},3/1),({23,20,35,21},8/1),({23,20,34,22},15/1),({23,20,33,23},23/1),({23,20,32,24},30/1),({23,20,31,25},33/1),({23,20,30,26},31/1),({23,20,29,27},22/1),({23,20,28,28},8/1),({22,21,37,19},1/1),({22,21,36,20},3/1),({22,21,35,21},6/1),({22,21,34,22},11/1),({22,21,33,23},16/1),({22,21,32,24},20/1),({22,21,31,25},22/1),({22,21,30,26},20/1),({22,21,29,27},14/1),({22,21,28,28},5/1)}, (15,0) => {}, (13,2) => {}, (15,1) => {({27,22,38,26},1/1),({27,22,37,27},1/1),({27,22,36,28},1/1),({27,22,35,29},1/1),({26,23,38,26},1/1),({26,23,37,27},1/1),({26,23,36,28},1/1),({26,23,35,29},1/1),({25,24,39,25},1/1),({25,24,38,26},1/1),({25,24,37,27},1/1),({25,24,36,28},1/1),({25,24,35,29},1/1)}, (15,2) => {({29,23,36,32},1/1),({29,23,34,34},1/1),({28,24,37,31},1/1),({28,24,35,33},1/1),({27,25,36,32},1/1),({27,25,34,34},1/1),({26,26,37,31},1/1),({26,26,35,33},1/1)}, (17,0) => {}, (17,1) => {}, (17,2) => {({29,29,39,37},1/1)}, (0,0) => {({1,0,0,0},1/1)}, (2,0) => {}, (2,1) => {({10,0,8,4},1/1),({10,0,7,5},1/1),({9,1,10,2},1/1),({9,1,9,3},1/1),({9,1,8,4},2/1),({9,1,7,5},1/1),({9,1,6,6},1/1),({8,2,10,2},1/1),({8,2,9,3},2/1),({8,2,8,4},2/1),({8,2,7,5},2/1),({8,2,6,6},1/1),({7,3,11,1},1/1),({7,3,10,2},2/1),({7,3,9,3},2/1),({7,3,8,4},4/1),({7,3,7,5},2/1),({7,3,6,6},1/1),({6,4,11,1},1/1),({6,4,10,2},1/1),({6,4,9,3},2/1),({6,4,8,4},2/1),({6,4,7,5},2/1),({5,5,12,0},1/1),({5,5,10,2},1/1),({5,5,9,3},1/1),({5,5,8,4},1/1),({5,5,6,6},1/1)}, (4,0) => {}, (2,2) => {}, (4,1) => {({15,1,13,7},1/1),({15,1,12,8},1/1),({15,1,11,9},1/1),({14,2,15,5},1/1),({14,2,14,6},2/1),({14,2,13,7},4/1),({14,2,12,8},5/1),({14,2,11,9},4/1),({14,2,10,10},2/1),({13,3,16,4},1/1),({13,3,15,5},4/1),({13,3,14,6},7/1),({13,3,13,7},12/1),({13,3,12,8},12/1),({13,3,11,9},11/1),({13,3,10,10},3/1),({12,4,17,3},1/1),({12,4,16,4},4/1),({12,4,15,5},10/1),({12,4,14,6},18/1),({12,4,13,7},24/1),({12,4,12,8},26/1),({12,4,11,9},20/1),({12,4,10,10},8/1),({11,5,17,3},3/1),({11,5,16,4},8/1),({11,5,15,5},18/1),({11,5,14,6},28/1),({11,5,13,7},38/1),({11,5,12,8},37/1),({11,5,11,9},30/1),({11,5,10,10},9/1),({10,6,18,2},2/1),({10,6,17,3},5/1),({10,6,16,4},13/1),({10,6,15,5},23/1),({10,6,14,6},36/1),({10,6,13,7},43/1),({10,6,12,8},45/1),({10,6,11,9},32/1),({10,6,10,10},13/1),({9,7,18,2},1/1),({9,7,17,3},5/1),({9,7,16,4},10/1),({9,7,15,5},21/1),({9,7,14,6},28/1),({9,7,13,7},37/1),({9,7,12,8},34/1),({9,7,11,9},28/1),({9,7,10,10},8/1),({8,8,18,2},1/1),({8,8,17,3},2/1),({8,8,16,4},5/1),({8,8,15,5},8/1),({8,8,14,6},13/1),({8,8,13,7},13/1),({8,8,12,8},15/1),({8,8,11,9},9/1),({8,8,10,10},5/1)}, (4,2) => {}, (6,0) => {}, (6,1) => {({19,3,18,10},1/1),({19,3,17,11},1/1),({19,3,16,12},2/1),({19,3,15,13},1/1),({19,3,14,14},1/1),({18,4,20,8},1/1),({18,4,19,9},2/1),({18,4,18,10},5/1),({18,4,17,11},8/1),({18,4,16,12},8/1),({18,4,15,13},7/1),({18,4,14,14},3/1),({17,5,21,7},1/1),({17,5,20,8},5/1),({17,5,19,9},11/1),({17,5,18,10},20/1),({17,5,17,11},26/1),({17,5,16,12},29/1),({17,5,15,13},21/1),({17,5,14,14},10/1),({16,6,22,6},1/1),({16,6,21,7},6/1),({16,6,20,8},16/1),({16,6,19,9},32/1),({16,6,18,10},50/1),({16,6,17,11},65/1),({16,6,16,12},65/1),({16,6,15,13},50/1),({16,6,14,14},18/1),({15,7,23,5},1/1),({15,7,22,6},6/1),({15,7,21,7},16/1),({15,7,20,8},39/1),({15,7,19,9},66/1),({15,7,18,10},100/1),({15,7,17,11},118/1),({15,7,16,12},120/1),({15,7,15,13},85/1),({15,7,14,14},35/1),({14,8,23,5},3/1),({14,8,22,6},11/1),({14,8,21,7},30/1),({14,8,20,8},61/1),({14,8,19,9},104/1),({14,8,18,10},145/1),({14,8,17,11},173/1),({14,8,16,12},165/1),({14,8,15,13},123/1),({14,8,14,14},43/1),({13,9,24,4},1/1),({13,9,23,5},5/1),({13,9,22,6},17/1),({13,9,21,7},39/1),({13,9,20,8},77/1),({13,9,19,9},121/1),({13,9,18,10},169/1),({13,9,17,11},190/1),({13,9,16,12},186/1),({13,9,15,13},130/1),({13,9,14,14},52/1),({12,10,24,4},1/1),({12,10,23,5},5/1),({12,10,22,6},14/1),({12,10,21,7},35/1),({12,10,20,8},63/1),({12,10,19,9},102/1),({12,10,18,10},133/1),({12,10,17,11},156/1),({12,10,16,12},143/1),({12,10,15,13},107/1),({12,10,14,14},36/1),({11,11,24,4},1/1),({11,11,23,5},2/1),({11,11,22,6},7/1),({11,11,21,7},13/1),({11,11,20,8},27/1),({11,11,19,9},38/1),({11,11,18,10},55/1),({11,11,17,11},57/1),({11,11,16,12},59/1),({11,11,15,13},37/1),({11,11,14,14},18/1)}, (6,2) => {}, (8,0) => {}, (8,1) => {({23,5,19,17},1/1),({22,6,23,13},1/1),({22,6,22,14},2/1),({22,6,21,15},3/1),({22,6,20,16},4/1),({22,6,19,17},3/1),({22,6,18,18},1/1),({21,7,25,11},1/1),({21,7,24,12},2/1),({21,7,23,13},7/1),({21,7,22,14},11/1),({21,7,21,15},17/1),({21,7,20,16},17/1),({21,7,19,17},15/1),({21,7,18,18},4/1),({20,8,26,10},1/1),({20,8,25,11},5/1),({20,8,24,12},14/1),({20,8,23,13},26/1),({20,8,22,14},42/1),({20,8,21,15},53/1),({20,8,20,16},54/1),({20,8,19,17},40/1),({20,8,18,18},16/1),({19,9,27,9},1/1),({19,9,26,10},6/1),({19,9,25,11},18/1),({19,9,24,12},38/1),({19,9,23,13},71/1),({19,9,22,14},100/1),({19,9,21,15},124/1),({19,9,20,16},119/1),({19,9,19,17},92/1),({19,9,18,18},30/1),({18,10,28,8},1/1),({18,10,27,9},5/1),({18,10,26,10},17/1),({18,10,25,11},42/1),({18,10,24,12},84/1),({18,10,23,13},136/1),({18,10,22,14},190/1),({18,10,21,15},219/1),({18,10,20,16},213/1),({18,10,19,17},153/1),({18,10,18,18},58/1),({17,11,28,8},2/1),({17,11,27,9},11/1),({17,11,26,10},30/1),({17,11,25,11},71/1),({17,11,24,12},129/1),({17,11,23,13},205/1),({17,11,22,14},270/1),({17,11,21,15},314/1),({17,11,20,16},289/1),({17,11,19,17},215/1),({17,11,18,18},74/1),({16,12,29,7},1/1),({16,12,28,8},5/1),({16,12,27,9},15/1),({16,12,26,10},42/1),({16,12,25,11},86/1),({16,12,24,12},155/1),({16,12,23,13},233/1),({16,12,22,14},308/1),({16,12,21,15},339/1),({16,12,20,16},323/1),({16,12,19,17},225/1),({16,12,18,18},85/1),({15,13,29,7},1/1),({15,13,28,8},4/1),({15,13,27,9},15/1),({15,13,26,10},35/1),({15,13,25,11},75/1),({15,13,24,12},125/1),({15,13,23,13},193/1),({15,13,22,14},242/1),({15,13,21,15},275/1),({15,13,20,16},247/1),({15,13,19,17},184/1),({15,13,18,18},60/1),({14,14,28,8},2/1),({14,14,27,9},6/1),({14,14,26,10},15/1),({14,14,25,11},28/1),({14,14,24,12},52/1),({14,14,23,13},71/1),({14,14,22,14},97/1),({14,14,21,15},102/1),({14,14,20,16},99/1),({14,14,19,17},64/1),({14,14,18,18},30/1)}, (8,2) => {}, (10,0) => {}, (10,1) => {({25,9,26,18},1/1),({25,9,25,19},1/1),({25,9,24,20},3/1),({25,9,23,21},1/1),({25,9,22,22},1/1),({24,10,28,16},1/1),({24,10,27,17},4/1),({24,10,26,18},7/1),({24,10,25,19},10/1),({24,10,24,20},11/1),({24,10,23,21},9/1),({24,10,22,22},3/1),({23,11,30,14},1/1),({23,11,29,15},3/1),({23,11,28,16},9/1),({23,11,27,17},17/1),({23,11,26,18},29/1),({23,11,25,19},35/1),({23,11,24,20},39/1),({23,11,23,21},27/1),({23,11,22,22},12/1),({22,12,31,13},1/1),({22,12,30,14},4/1),({22,12,29,15},13/1),({22,12,28,16},28/1),({22,12,27,17},51/1),({22,12,26,18},73/1),({22,12,25,19},90/1),({22,12,24,20},87/1),({22,12,23,21},66/1),({22,12,22,22},23/1),({21,13,32,12},1/1),({21,13,31,13},4/1),({21,13,30,14},14/1),({21,13,29,15},32/1),({21,13,28,16},65/1),({21,13,27,17},102/1),({21,13,26,18},145/1),({21,13,25,19},164/1),({21,13,24,20},163/1),({21,13,23,21},113/1),({21,13,22,22},46/1),({20,14,32,12},2/1),({20,14,31,13},9/1),({20,14,30,14},25/1),({20,14,29,15},57/1),({20,14,28,16},102/1),({20,14,27,17},161/1),({20,14,26,18},211/1),({20,14,25,19},242/1),({20,14,24,20},224/1),({20,14,23,21},165/1),({20,14,22,22},57/1),({19,15,33,11},1/1),({19,15,32,12},4/1),({19,15,31,13},13/1),({19,15,30,14},35/1),({19,15,29,15},70/1),({19,15,28,16},126/1),({19,15,27,17},185/1),({19,15,26,18},245/1),({19,15,25,19},266/1),({19,15,24,20},254/1),({19,15,23,21},174/1),({19,15,22,22},69/1),({18,16,33,11},1/1),({18,16,32,12},4/1),({18,16,31,13},13/1),({18,16,30,14},30/1),({18,16,29,15},62/1),({18,16,28,16},103/1),({18,16,27,17},156/1),({18,16,26,18},194/1),({18,16,25,19},219/1),({18,16,24,20},196/1),({18,16,23,21},145/1),({18,16,22,22},47/1),({17,17,32,12},2/1),({17,17,31,13},4/1),({17,17,30,14},13/1),({17,17,29,15},23/1),({17,17,28,16},43/1),({17,17,27,17},57/1),({17,17,26,18},79/1),({17,17,25,19},79/1),({17,17,24,20},81/1),({17,17,23,21},49/1),({17,17,22,22},24/1)}, (10,2) => {}, (12,0) => {}, (12,1) => {({27,13,29,23},2/1),({27,13,28,24},1/1),({27,13,27,25},2/1),({26,14,32,20},1/1),({26,14,31,21},2/1),({26,14,30,22},5/1),({26,14,29,23},7/1),({26,14,28,24},8/1),({26,14,27,25},6/1),({26,14,26,26},3/1),({25,15,33,19},2/1),({25,15,32,20},5/1),({25,15,31,21},11/1),({25,15,30,22},16/1),({25,15,29,23},23/1),({25,15,28,24},21/1),({25,15,27,25},19/1),({25,15,26,26},5/1),({24,16,35,17},1/1),({24,16,34,18},3/1),({24,16,33,19},7/1),({24,16,32,20},16/1),({24,16,31,21},27/1),({24,16,30,22},40/1),({24,16,29,23},47/1),({24,16,28,24},47/1),({24,16,27,25},33/1),({24,16,26,26},14/1),({23,17,35,17},2/1),({23,17,34,18},6/1),({23,17,33,19},16/1),({23,17,32,20},29/1),({23,17,31,21},49/1),({23,17,30,22},63/1),({23,17,29,23},77/1),({23,17,28,24},69/1),({23,17,27,25},54/1),({23,17,26,26},16/1),({22,18,36,16},1/1),({22,18,35,17},3/1),({22,18,34,18},10/1),({22,18,33,19},21/1),({22,18,32,20},40/1),({22,18,31,21},59/1),({22,18,30,22},81/1),({22,18,29,23},87/1),({22,18,28,24},85/1),({22,18,27,25},57/1),({22,18,26,26},24/1),({21,19,36,16},1/1),({21,19,35,17},4/1),({21,19,34,18},9/1),({21,19,33,19},21/1),({21,19,32,20},34/1),({21,19,31,21},54/1),({21,19,30,22},65/1),({21,19,29,23},77/1),({21,19,28,24},65/1),({21,19,27,25},52/1),({21,19,26,26},14/1),({20,20,35,17},1/1),({20,20,34,18},4/1),({20,20,33,19},7/1),({20,20,32,20},15/1),({20,20,31,21},19/1),({20,20,30,22},28/1),({20,20,29,23},26/1),({20,20,28,24},29/1),({20,20,27,25},16/1),({20,20,26,26},10/1)}, (12,2) => {}, (14,0) => {}, (14,1) => {({29,17,30,30},1/1),({28,18,33,27},1/1),({28,18,32,28},1/1),({28,18,31,29},1/1),({27,19,36,24},1/1),({27,19,35,25},1/1),({27,19,34,26},3/1),({27,19,33,27},2/1),({27,19,32,28},4/1),({27,19,31,29},1/1),({27,19,30,30},2/1),({26,20,37,23},1/1),({26,20,36,24},2/1),({26,20,35,25},4/1),({26,20,34,26},5/1),({26,20,33,27},7/1),({26,20,32,28},5/1),({26,20,31,29},5/1),({26,20,30,30},1/1),({25,21,38,22},1/1),({25,21,37,23},2/1),({25,21,36,24},4/1),({25,21,35,25},6/1),({25,21,34,26},9/1),({25,21,33,27},7/1),({25,21,32,28},9/1),({25,21,31,29},4/1),({25,21,30,30},3/1),({24,22,38,22},1/1),({24,22,37,23},3/1),({24,22,36,24},4/1),({24,22,35,25},7/1),({24,22,34,26},7/1),({24,22,33,27},9/1),({24,22,32,28},6/1),({24,22,31,29},6/1),({23,23,36,24},2/1),({23,23,35,25},1/1),({23,23,34,26},4/1),({23,23,33,27},2/1),({23,23,32,28},4/1),({23,23,30,30},3/1)}, (14,2) => {}, (16,0) => {}, (16,1) => {({27,25,39,29},1/1)}, (16,2) => {({29,26,38,34},1/1),({29,26,37,35},1/1)}, (18,2) => {}, (1,0) => {({3,1,4,0},1/1)}, (1,1) => {({7,0,6,2},1/1),({7,0,4,4},1/1)}, (3,0) => {}, (3,1) => {({13,0,9,7},1/1),({12,1,12,4},1/1),({12,1,11,5},2/1),({12,1,10,6},2/1),({12,1,9,7},2/1),({12,1,8,8},1/1),({11,2,13,3},1/1),({11,2,12,4},2/1),({11,2,11,5},5/1),({11,2,10,6},5/1),({11,2,9,7},5/1),({11,2,8,8},1/1),({10,3,14,2},1/1),({10,3,13,3},3/1),({10,3,12,4},6/1),({10,3,11,5},9/1),({10,3,10,6},10/1),({10,3,9,7},8/1),({10,3,8,8},3/1),({9,4,14,2},2/1),({9,4,13,3},5/1),({9,4,12,4},9/1),({9,4,11,5},13/1),({9,4,10,6},13/1),({9,4,9,7},11/1),({9,4,8,8},4/1),({8,5,15,1},1/1),({8,5,14,2},3/1),({8,5,13,3},6/1),({8,5,12,4},10/1),({8,5,11,5},13/1),({8,5,10,6},13/1),({8,5,9,7},10/1),({8,5,8,8},4/1),({7,6,15,1},1/1),({7,6,14,2},2/1),({7,6,13,3},4/1),({7,6,12,4},6/1),({7,6,11,5},9/1),({7,6,10,6},8/1),({7,6,9,7},7/1),({7,6,8,8},2/1)}, (5,0) => {}, (3,2) => {}, (5,1) => {({17,2,16,8},1/1),({17,2,15,9},1/1),({17,2,14,10},2/1),({17,2,13,11},1/1),({17,2,12,12},1/1),({16,3,17,7},2/1),({16,3,16,8},4/1),({16,3,15,9},6/1),({16,3,14,10},8/1),({16,3,13,11},6/1),({16,3,12,12},2/1),({15,4,19,5},1/1),({15,4,18,6},3/1),({15,4,17,7},8/1),({15,4,16,8},15/1),({15,4,15,9},20/1),({15,4,14,10},23/1),({15,4,13,11},17/1),({15,4,12,12},7/1),({14,5,19,5},3/1),({14,5,18,6},10/1),({14,5,17,7},21/1),({14,5,16,8},35/1),({14,5,15,9},46/1),({14,5,14,10},47/1),({14,5,13,11},36/1),({14,5,12,12},14/1),({13,6,20,4},3/1),({13,6,19,5},9/1),({13,6,18,6},22/1),({13,6,17,7},41/1),({13,6,16,8},63/1),({13,6,15,9},76/1),({13,6,14,10},79/1),({13,6,13,11},57/1),({13,6,12,12},22/1),({12,7,21,3},1/1),({12,7,20,4},5/1),({12,7,19,5},15/1),({12,7,18,6},33/1),({12,7,17,7},58/1),({12,7,16,8},84/1),({12,7,15,9},101/1),({12,7,14,10},99/1),({12,7,13,11},73/1),({12,7,12,12},27/1),({11,8,21,3},2/1),({11,8,20,4},6/1),({11,8,19,5},17/1),({11,8,18,6},36/1),({11,8,17,7},59/1),({11,8,16,8},84/1),({11,8,15,9},99/1),({11,8,14,10},95/1),({11,8,13,11},69/1),({11,8,12,12},27/1),({10,9,21,3},1/1),({10,9,20,4},5/1),({10,9,19,5},12/1),({10,9,18,6},23/1),({10,9,17,7},39/1),({10,9,16,8},53/1),({10,9,15,9},61/1),({10,9,14,10},60/1),({10,9,13,11},43/1),({10,9,12,12},15/1)}};
--dw stands for dominant weights
<<<<<<< HEAD
<<<<<<< HEAD
dw1034 = new HashTable from {(5,2) => {}, (7,0) => {}, (7,1) => {({21,4,19,13},1/1),({20,5,22,10},1/1),({19,6,23,9},1/1),({18,7,24,8},2/1),({17,8,25,7},1/1),({16,9,26,6},1/1),({13,12,27,5},1/1)}, (9,0) => {}, (7,2) => {}, (9,1) => {({24,7,23,17},1/1),({23,8,26,14},1/1),({22,9,27,13},2/1),({21,10,29,11},1/1),({19,12,30,10},2/1),({18,13,31,9},1/1)}, (9,2) => {}, (11,0) => {}, (11,1) => {({26,11,28,20},1/1),({25,12,30,18},1/1),({24,13,32,16},1/1),({23,14,33,15},1/1),({22,15,34,14},1/1),({19,18,35,13},1/1)}, (11,2) => {}, (13,0) => {}, (13,1) => {({28,15,30,26},1/1),({27,16,33,23},1/1),({26,17,35,21},1/1),({25,18,36,20},1/1),({24,19,37,19},1/1)}, (15,0) => {}, (13,2) => {}, (15,1) => {({27,22,38,26},1/1),({25,24,39,25},1/1)}, (15,2) => {({29,23,36,32},1/1),({28,24,37,31},1/1)}, (17,0) => {}, (17,1) => {}, (17,2) => {({29,29,39,37},1/1)}, (0,0) => {({1,0,0,0},1/1)}, (2,0) => {}, (2,1) => {({10,0,8,4},1/1),({9,1,10,2},1/1),({7,3,11,1},1/1),({5,5,12,0},1/1)}, (4,0) => {}, (2,2) => {}, (4,1) => {({15,1,13,7},1/1),({14,2,15,5},1/1),({13,3,16,4},1/1),({12,4,17,3},1/1),({10,6,18,2},2/1)}, (4,2) => {}, (6,0) => {}, (6,1) => {({19,3,18,10},1/1),({18,4,20,8},1/1),({17,5,21,7},1/1),({16,6,22,6},1/1),({15,7,23,5},1/1),({13,9,24,4},1/1)}, (6,2) => {}, (8,0) => {}, (8,1) => {({23,5,19,17},1/1),({22,6,23,13},1/1),({21,7,25,11},1/1),({20,8,26,10},1/1),({19,9,27,9},1/1),({18,10,28,8},1/1),({16,12,29,7},1/1)}, (8,2) => {}, (10,0) => {}, (10,1) => {({25,9,26,18},1/1),({24,10,28,16},1/1),({23,11,30,14},1/1),({22,12,31,13},1/1),({21,13,32,12},1/1),({19,15,33,11},1/1)}, (10,2) => {}, (12,0) => {}, (12,1) => {({27,13,29,23},2/1),({26,14,32,20},1/1),({25,15,33,19},2/1),({24,16,35,17},1/1),({22,18,36,16},1/1)}, (12,2) => {}, (14,0) => {}, (14,1) => {({29,17,30,30},1/1),({28,18,33,27},1/1),({27,19,36,24},1/1),({26,20,37,23},1/1),({25,21,38,22},1/1)}, (14,2) => {}, (16,0) => {}, (16,1) => {({27,25,39,29},1/1)}, (16,2) => {({29,26,38,34},1/1)}, (18,2) => {}, (1,0) => {({3,1,4,0},1/1)}, (1,1) => {({7,0,6,2},1/1)}, (3,0) => {}, (3,1) => {({13,0,9,7},1/1),({12,1,12,4},1/1),({11,2,13,3},1/1),({10,3,14,2},1/1),({8,5,15,1},1/1)}, (5,0) => {}, (3,2) => {}, (5,1) => {({17,2,16,8},1/1),({16,3,17,7},2/1),({15,4,19,5},1/1),({13,6,20,4},3/1),({12,7,21,3},1/1)}};
--dmw stands for dominant monomial weights
dmw1034 = {infinity};
=======
dw1034 = new HashTable from {(5,2) => {}, (7,0) => {}, (7,1) => {({21,4,19,13},1/1),({20,5,22,10},1/1),({19,6,23,9},1/1),({18,7,24,8},2/1),({17,8,25,7},1/1),({16,9,26,6},1/1),({13,12,27,5},1/1)}, (7,2) => {}, (9,0) => {}, (9,1) => {({24,7,23,17},1/1),({23,8,26,14},1/1),({22,9,27,13},2/1),({21,10,29,11},1/1),({19,12,30,10},2/1),({18,13,31,9},1/1)}, (9,2) => {}, (11,0) => {}, (11,1) => {({26,11,28,20},1/1),({25,12,30,18},1/1),({24,13,32,16},1/1),({23,14,33,15},1/1),({22,15,34,14},1/1),({19,18,35,13},1/1)}, (11,2) => {}, (13,0) => {}, (13,1) => {({28,15,30,26},1/1),({27,16,33,23},1/1),({26,17,35,21},1/1),({25,18,36,20},1/1),({24,19,37,19},1/1)}, (13,2) => {}, (15,0) => {}, (15,1) => {({27,22,38,26},1/1),({25,24,39,25},1/1)}, (17,0) => {}, (15,2) => {({29,23,36,32},1/1),({28,24,37,31},1/1)}, (17,1) => {}, (17,2) => {({29,29,39,37},1/1)}, (0,0) => {({1,0,0,0},1/1)}, (2,0) => {}, (2,1) => {({10,0,8,4},1/1),({9,1,10,2},1/1),({7,3,11,1},1/1),({5,5,12,0},1/1)}, (2,2) => {}, (4,0) => {}, (4,1) => {({15,1,13,7},1/1),({14,2,15,5},1/1),({13,3,16,4},1/1),({12,4,17,3},1/1),({10,6,18,2},2/1)}, (4,2) => {}, (6,0) => {}, (6,1) => {({19,3,18,10},1/1),({18,4,20,8},1/1),({17,5,21,7},1/1),({16,6,22,6},1/1),({15,7,23,5},1/1),({13,9,24,4},1/1)}, (6,2) => {}, (8,0) => {}, (8,1) => {({23,5,19,17},1/1),({22,6,23,13},1/1),({21,7,25,11},1/1),({20,8,26,10},1/1),({19,9,27,9},1/1),({18,10,28,8},1/1),({16,12,29,7},1/1)}, (8,2) => {}, (10,0) => {}, (10,1) => {({25,9,26,18},1/1),({24,10,28,16},1/1),({23,11,30,14},1/1),({22,12,31,13},1/1),({21,13,32,12},1/1),({19,15,33,11},1/1)}, (10,2) => {}, (12,0) => {}, (12,1) => {({27,13,29,23},2/1),({26,14,32,20},1/1),({25,15,33,19},2/1),({24,16,35,17},1/1),({22,18,36,16},1/1)}, (12,2) => {}, (14,0) => {}, (14,1) => {({29,17,30,30},1/1),({28,18,33,27},1/1),({27,19,36,24},1/1),({26,20,37,23},1/1),({25,21,38,22},1/1)}, (14,2) => {}, (16,0) => {}, (16,1) => {({27,25,39,29},1/1)}, (16,2) => {({29,26,38,34},1/1)}, (18,2) => {}, (1,0) => {({3,1,4,0},1/1)}, (1,1) => {({7,0,6,2},1/1)}, (3,0) => {}, (3,1) => {({13,0,9,7},1/1),({12,1,12,4},1/1),({11,2,13,3},1/1),({10,3,14,2},1/1),({8,5,15,1},1/1)}, (3,2) => {}, (5,0) => {}, (5,1) => {({17,2,16,8},1/1),({16,3,17,7},2/1),({15,4,19,5},1/1),({13,6,20,4},3/1),({12,7,21,3},1/1)}};
>>>>>>> parent of 67e86b7... updating data
=======
dw1034 = new HashTable from {(5,2) => {}, (7,0) => {}, (7,1) => {({21,4,19,13},1/1),({20,5,22,10},1/1),({19,6,23,9},1/1),({18,7,24,8},2/1),({17,8,25,7},1/1),({16,9,26,6},1/1),({13,12,27,5},1/1)}, (7,2) => {}, (9,0) => {}, (9,1) => {({24,7,23,17},1/1),({23,8,26,14},1/1),({22,9,27,13},2/1),({21,10,29,11},1/1),({19,12,30,10},2/1),({18,13,31,9},1/1)}, (9,2) => {}, (11,0) => {}, (11,1) => {({26,11,28,20},1/1),({25,12,30,18},1/1),({24,13,32,16},1/1),({23,14,33,15},1/1),({22,15,34,14},1/1),({19,18,35,13},1/1)}, (11,2) => {}, (13,0) => {}, (13,1) => {({28,15,30,26},1/1),({27,16,33,23},1/1),({26,17,35,21},1/1),({25,18,36,20},1/1),({24,19,37,19},1/1)}, (13,2) => {}, (15,0) => {}, (15,1) => {({27,22,38,26},1/1),({25,24,39,25},1/1)}, (17,0) => {}, (15,2) => {({29,23,36,32},1/1),({28,24,37,31},1/1)}, (17,1) => {}, (17,2) => {({29,29,39,37},1/1)}, (0,0) => {({1,0,0,0},1/1)}, (2,0) => {}, (2,1) => {({10,0,8,4},1/1),({9,1,10,2},1/1),({7,3,11,1},1/1),({5,5,12,0},1/1)}, (2,2) => {}, (4,0) => {}, (4,1) => {({15,1,13,7},1/1),({14,2,15,5},1/1),({13,3,16,4},1/1),({12,4,17,3},1/1),({10,6,18,2},2/1)}, (4,2) => {}, (6,0) => {}, (6,1) => {({19,3,18,10},1/1),({18,4,20,8},1/1),({17,5,21,7},1/1),({16,6,22,6},1/1),({15,7,23,5},1/1),({13,9,24,4},1/1)}, (6,2) => {}, (8,0) => {}, (8,1) => {({23,5,19,17},1/1),({22,6,23,13},1/1),({21,7,25,11},1/1),({20,8,26,10},1/1),({19,9,27,9},1/1),({18,10,28,8},1/1),({16,12,29,7},1/1)}, (8,2) => {}, (10,0) => {}, (10,1) => {({25,9,26,18},1/1),({24,10,28,16},1/1),({23,11,30,14},1/1),({22,12,31,13},1/1),({21,13,32,12},1/1),({19,15,33,11},1/1)}, (10,2) => {}, (12,0) => {}, (12,1) => {({27,13,29,23},2/1),({26,14,32,20},1/1),({25,15,33,19},2/1),({24,16,35,17},1/1),({22,18,36,16},1/1)}, (12,2) => {}, (14,0) => {}, (14,1) => {({29,17,30,30},1/1),({28,18,33,27},1/1),({27,19,36,24},1/1),({26,20,37,23},1/1),({25,21,38,22},1/1)}, (14,2) => {}, (16,0) => {}, (16,1) => {({27,25,39,29},1/1)}, (16,2) => {({29,26,38,34},1/1)}, (18,2) => {}, (1,0) => {({3,1,4,0},1/1)}, (1,1) => {({7,0,6,2},1/1)}, (3,0) => {}, (3,1) => {({13,0,9,7},1/1),({12,1,12,4},1/1),({11,2,13,3},1/1),({10,3,14,2},1/1),({8,5,15,1},1/1)}, (3,2) => {}, (5,0) => {}, (5,1) => {({17,2,16,8},1/1),({16,3,17,7},2/1),({15,4,19,5},1/1),({13,6,20,4},3/1),({12,7,21,3},1/1)}};
>>>>>>> parent of 67e86b7... updating data
end;