A := QQ[t_0,t_1,t_2,t_3];
--tb stands for Total Betti numbers
tb0234 = new HashTable from {(5,2) => 0, (7,0) => 0, (7,1) => 271830/1, (9,0) => 0, (7,2) => 0, (9,1) => 354926/1, (9,2) => 0, (11,0) => 0, (11,1) => 177684/1, (11,2) => 0, (13,0) => 0, (13,1) => 30804/1, (15,0) => 0, (13,2) => 0, (15,1) => 1173/1, (15,2) => 0, (17,0) => 0, (17,1) => 0, (17,2) => 2/1, (0,0) => 3/1, (2,0) => 120/1, (2,1) => 510/1, (4,0) => 0, (2,2) => 0, (4,1) => 25296/1, (4,2) => 0, (6,0) => 0, (6,1) => 164424/1, (6,2) => 0, (8,0) => 0, (8,1) => 350064/1, (8,2) => 0, (10,0) => 0, (10,1) => 283764/1, (10,2) => 0, (12,0) => 0, (12,1) => 85680/1, (12,2) => 0, (14,0) => 0, (14,1) => 7752/1, (14,2) => 0, (16,0) => 0, (16,1) => 48/1, (16,2) => 15/1, (18,2) => 0, (1,0) => 32/1, (1,1) => 33/1, (3,0) => 0, (3,1) => 5508/1, (5,0) => 0, (3,2) => 0, (5,1) => 75684/1};
--sb represents the betti numbers as sums of Schur functors
sb0234 = new HashTable from {(5,2) => {}, (7,0) => {}, (7,1) => {({20,4,21,13},1/1),({20,4,20,14},1/1),({20,4,19,15},3/1),({20,4,18,16},1/1),({20,4,17,17},1/1),({19,5,23,11},1/1),({19,5,22,12},4/1),({19,5,21,13},7/1),({19,5,20,14},10/1),({19,5,19,15},11/1),({19,5,18,16},9/1),({19,5,17,17},3/1),({18,6,25,9},1/1),({18,6,24,10},3/1),({18,6,23,11},9/1),({18,6,22,12},17/1),({18,6,21,13},29/1),({18,6,20,14},35/1),({18,6,19,15},39/1),({18,6,18,16},27/1),({18,6,17,17},12/1),({17,7,26,8},1/1),({17,7,25,9},4/1),({17,7,24,10},13/1),({17,7,23,11},28/1),({17,7,22,12},51/1),({17,7,21,13},73/1),({17,7,20,14},90/1),({17,7,19,15},87/1),({17,7,18,16},66/1),({17,7,17,17},23/1),({16,8,27,7},1/1),({16,8,26,8},4/1),({16,8,25,9},14/1),({16,8,24,10},32/1),({16,8,23,11},65/1),({16,8,22,12},102/1),({16,8,21,13},145/1),({16,8,20,14},164/1),({16,8,19,15},163/1),({16,8,18,16},113/1),({16,8,17,17},46/1),({15,9,27,7},2/1),({15,9,26,8},9/1),({15,9,25,9},25/1),({15,9,24,10},57/1),({15,9,23,11},102/1),({15,9,22,12},161/1),({15,9,21,13},211/1),({15,9,20,14},242/1),({15,9,19,15},224/1),({15,9,18,16},165/1),({15,9,17,17},57/1),({14,10,28,6},1/1),({14,10,27,7},4/1),({14,10,26,8},13/1),({14,10,25,9},35/1),({14,10,24,10},70/1),({14,10,23,11},126/1),({14,10,22,12},185/1),({14,10,21,13},245/1),({14,10,20,14},266/1),({14,10,19,15},254/1),({14,10,18,16},174/1),({14,10,17,17},69/1),({13,11,28,6},1/1),({13,11,27,7},4/1),({13,11,26,8},13/1),({13,11,25,9},30/1),({13,11,24,10},62/1),({13,11,23,11},103/1),({13,11,22,12},156/1),({13,11,21,13},194/1),({13,11,20,14},219/1),({13,11,19,15},196/1),({13,11,18,16},145/1),({13,11,17,17},47/1),({12,12,27,7},2/1),({12,12,26,8},4/1),({12,12,25,9},13/1),({12,12,24,10},23/1),({12,12,23,11},43/1),({12,12,22,12},57/1),({12,12,21,13},79/1),({12,12,20,14},79/1),({12,12,19,15},81/1),({12,12,18,16},49/1),({12,12,17,17},24/1)}, (9,0) => {}, (7,2) => {}, (9,1) => {({24,6,22,20},1/1),({23,7,26,16},1/1),({23,7,25,17},2/1),({23,7,24,18},3/1),({23,7,23,19},4/1),({23,7,22,20},3/1),({23,7,21,21},1/1),({22,8,28,14},1/1),({22,8,27,15},2/1),({22,8,26,16},7/1),({22,8,25,17},11/1),({22,8,24,18},17/1),({22,8,23,19},17/1),({22,8,22,20},15/1),({22,8,21,21},4/1),({21,9,29,13},1/1),({21,9,28,14},5/1),({21,9,27,15},14/1),({21,9,26,16},26/1),({21,9,25,17},42/1),({21,9,24,18},53/1),({21,9,23,19},54/1),({21,9,22,20},40/1),({21,9,21,21},16/1),({20,10,30,12},1/1),({20,10,29,13},6/1),({20,10,28,14},18/1),({20,10,27,15},38/1),({20,10,26,16},71/1),({20,10,25,17},100/1),({20,10,24,18},124/1),({20,10,23,19},119/1),({20,10,22,20},92/1),({20,10,21,21},30/1),({19,11,31,11},1/1),({19,11,30,12},5/1),({19,11,29,13},17/1),({19,11,28,14},42/1),({19,11,27,15},84/1),({19,11,26,16},136/1),({19,11,25,17},190/1),({19,11,24,18},219/1),({19,11,23,19},213/1),({19,11,22,20},153/1),({19,11,21,21},58/1),({18,12,31,11},2/1),({18,12,30,12},11/1),({18,12,29,13},30/1),({18,12,28,14},71/1),({18,12,27,15},129/1),({18,12,26,16},205/1),({18,12,25,17},270/1),({18,12,24,18},314/1),({18,12,23,19},289/1),({18,12,22,20},215/1),({18,12,21,21},74/1),({17,13,32,10},1/1),({17,13,31,11},5/1),({17,13,30,12},15/1),({17,13,29,13},42/1),({17,13,28,14},86/1),({17,13,27,15},155/1),({17,13,26,16},233/1),({17,13,25,17},308/1),({17,13,24,18},339/1),({17,13,23,19},323/1),({17,13,22,20},225/1),({17,13,21,21},85/1),({16,14,32,10},1/1),({16,14,31,11},4/1),({16,14,30,12},15/1),({16,14,29,13},35/1),({16,14,28,14},75/1),({16,14,27,15},125/1),({16,14,26,16},193/1),({16,14,25,17},242/1),({16,14,24,18},275/1),({16,14,23,19},247/1),({16,14,22,20},184/1),({16,14,21,21},60/1),({15,15,31,11},2/1),({15,15,30,12},6/1),({15,15,29,13},15/1),({15,15,28,14},28/1),({15,15,27,15},52/1),({15,15,26,16},71/1),({15,15,25,17},97/1),({15,15,24,18},102/1),({15,15,23,19},99/1),({15,15,22,20},64/1),({15,15,21,21},30/1)}, (9,2) => {}, (11,0) => {}, (11,1) => {({26,10,29,21},1/1),({26,10,28,22},1/1),({26,10,27,23},2/1),({26,10,26,24},1/1),({26,10,25,25},1/1),({25,11,31,19},1/1),({25,11,30,20},2/1),({25,11,29,21},5/1),({25,11,28,22},8/1),({25,11,27,23},8/1),({25,11,26,24},7/1),({25,11,25,25},3/1),({24,12,32,18},1/1),({24,12,31,19},5/1),({24,12,30,20},11/1),({24,12,29,21},20/1),({24,12,28,22},26/1),({24,12,27,23},29/1),({24,12,26,24},21/1),({24,12,25,25},10/1),({23,13,33,17},1/1),({23,13,32,18},6/1),({23,13,31,19},16/1),({23,13,30,20},32/1),({23,13,29,21},50/1),({23,13,28,22},65/1),({23,13,27,23},65/1),({23,13,26,24},50/1),({23,13,25,25},18/1),({22,14,34,16},1/1),({22,14,33,17},6/1),({22,14,32,18},16/1),({22,14,31,19},39/1),({22,14,30,20},66/1),({22,14,29,21},100/1),({22,14,28,22},118/1),({22,14,27,23},120/1),({22,14,26,24},85/1),({22,14,25,25},35/1),({21,15,34,16},3/1),({21,15,33,17},11/1),({21,15,32,18},30/1),({21,15,31,19},61/1),({21,15,30,20},104/1),({21,15,29,21},145/1),({21,15,28,22},173/1),({21,15,27,23},165/1),({21,15,26,24},123/1),({21,15,25,25},43/1),({20,16,35,15},1/1),({20,16,34,16},5/1),({20,16,33,17},17/1),({20,16,32,18},39/1),({20,16,31,19},77/1),({20,16,30,20},121/1),({20,16,29,21},169/1),({20,16,28,22},190/1),({20,16,27,23},186/1),({20,16,26,24},130/1),({20,16,25,25},52/1),({19,17,35,15},1/1),({19,17,34,16},5/1),({19,17,33,17},14/1),({19,17,32,18},35/1),({19,17,31,19},63/1),({19,17,30,20},102/1),({19,17,29,21},133/1),({19,17,28,22},156/1),({19,17,27,23},143/1),({19,17,26,24},107/1),({19,17,25,25},36/1),({18,18,35,15},1/1),({18,18,34,16},2/1),({18,18,33,17},7/1),({18,18,32,18},13/1),({18,18,31,19},27/1),({18,18,30,20},38/1),({18,18,29,21},55/1),({18,18,28,22},57/1),({18,18,27,23},59/1),({18,18,26,24},37/1),({18,18,25,25},18/1)}, (11,2) => {}, (13,0) => {}, (13,1) => {({28,14,32,26},1/1),({28,14,31,27},1/1),({28,14,30,28},1/1),({27,15,34,24},1/1),({27,15,33,25},2/1),({27,15,32,26},4/1),({27,15,31,27},5/1),({27,15,30,28},4/1),({27,15,29,29},2/1),({26,16,35,23},1/1),({26,16,34,24},4/1),({26,16,33,25},7/1),({26,16,32,26},12/1),({26,16,31,27},12/1),({26,16,30,28},11/1),({26,16,29,29},3/1),({25,17,36,22},1/1),({25,17,35,23},4/1),({25,17,34,24},10/1),({25,17,33,25},18/1),({25,17,32,26},24/1),({25,17,31,27},26/1),({25,17,30,28},20/1),({25,17,29,29},8/1),({24,18,36,22},3/1),({24,18,35,23},8/1),({24,18,34,24},18/1),({24,18,33,25},28/1),({24,18,32,26},38/1),({24,18,31,27},37/1),({24,18,30,28},30/1),({24,18,29,29},9/1),({23,19,37,21},2/1),({23,19,36,22},5/1),({23,19,35,23},13/1),({23,19,34,24},23/1),({23,19,33,25},36/1),({23,19,32,26},43/1),({23,19,31,27},45/1),({23,19,30,28},32/1),({23,19,29,29},13/1),({22,20,37,21},1/1),({22,20,36,22},5/1),({22,20,35,23},10/1),({22,20,34,24},21/1),({22,20,33,25},28/1),({22,20,32,26},37/1),({22,20,31,27},34/1),({22,20,30,28},28/1),({22,20,29,29},8/1),({21,21,37,21},1/1),({21,21,36,22},2/1),({21,21,35,23},5/1),({21,21,34,24},8/1),({21,21,33,25},13/1),({21,21,32,26},13/1),({21,21,31,27},15/1),({21,21,30,28},9/1),({21,21,29,29},5/1)}, (15,0) => {}, (13,2) => {}, (15,1) => {({29,19,35,31},1/1),({29,19,34,32},1/1),({28,20,37,29},1/1),({28,20,36,30},1/1),({28,20,35,31},2/1),({28,20,34,32},1/1),({28,20,33,33},1/1),({27,21,37,29},1/1),({27,21,36,30},2/1),({27,21,35,31},2/1),({27,21,34,32},2/1),({27,21,33,33},1/1),({26,22,38,28},1/1),({26,22,37,29},2/1),({26,22,36,30},2/1),({26,22,35,31},4/1),({26,22,34,32},2/1),({26,22,33,33},1/1),({25,23,38,28},1/1),({25,23,37,29},1/1),({25,23,36,30},2/1),({25,23,35,31},2/1),({25,23,34,32},2/1),({24,24,39,27},1/1),({24,24,37,29},1/1),({24,24,36,30},1/1),({24,24,35,31},1/1),({24,24,33,33},1/1)}, (15,2) => {}, (17,0) => {}, (17,1) => {}, (17,2) => {({29,28,39,39},1/1)}, (0,0) => {({0,0,2,0},1/1)}, (2,0) => {({6,0,7,3},1/1),({6,0,5,5},1/1),({5,1,8,2},1/1),({5,1,6,4},1/1),({4,2,7,3},1/1),({4,2,5,5},1/1),({3,3,8,2},1/1),({3,3,6,4},1/1)}, (2,1) => {({7,2,13,1},1/1),({7,2,12,2},1/1),({7,2,11,3},1/1),({7,2,10,4},1/1),({6,3,13,1},1/1),({6,3,12,2},1/1),({6,3,11,3},1/1),({6,3,10,4},1/1),({5,4,14,0},1/1),({5,4,13,1},1/1),({5,4,12,2},1/1),({5,4,11,3},1/1),({5,4,10,4},1/1)}, (4,0) => {}, (2,2) => {}, (4,1) => {({14,1,13,9},1/1),({14,1,12,10},1/1),({13,2,16,6},1/1),({13,2,15,7},2/1),({13,2,14,8},3/1),({13,2,13,9},4/1),({13,2,12,10},3/1),({13,2,11,11},1/1),({12,3,18,4},1/1),({12,3,17,5},2/1),({12,3,16,6},5/1),({12,3,15,7},8/1),({12,3,14,8},10/1),({12,3,13,9},11/1),({12,3,12,10},8/1),({12,3,11,11},3/1),({11,4,19,3},1/1),({11,4,18,4},3/1),({11,4,17,5},7/1),({11,4,16,6},12/1),({11,4,15,7},18/1),({11,4,14,8},21/1),({11,4,13,9},21/1),({11,4,12,10},15/1),({11,4,11,11},6/1),({10,5,20,2},1/1),({10,5,19,3},3/1),({10,5,18,4},7/1),({10,5,17,5},13/1),({10,5,16,6},21/1),({10,5,15,7},28/1),({10,5,14,8},32/1),({10,5,13,9},30/1),({10,5,12,10},21/1),({10,5,11,11},8/1),({9,6,20,2},1/1),({9,6,19,3},3/1),({9,6,18,4},8/1),({9,6,17,5},15/1),({9,6,16,6},23/1),({9,6,15,7},30/1),({9,6,14,8},33/1),({9,6,13,9},31/1),({9,6,12,10},22/1),({9,6,11,11},8/1),({8,7,20,2},1/1),({8,7,19,3},3/1),({8,7,18,4},6/1),({8,7,17,5},11/1),({8,7,16,6},16/1),({8,7,15,7},20/1),({8,7,14,8},22/1),({8,7,13,9},20/1),({8,7,12,10},14/1),({8,7,11,11},5/1)}, (4,2) => {}, (6,0) => {}, (6,1) => {({18,3,19,11},1/1),({18,3,18,12},2/1),({18,3,17,13},2/1),({18,3,16,14},2/1),({18,3,15,15},1/1),({17,4,21,9},1/1),({17,4,20,10},4/1),({17,4,19,11},7/1),({17,4,18,12},10/1),({17,4,17,13},11/1),({17,4,16,14},9/1),({17,4,15,15},3/1),({16,5,23,7},1/1),({16,5,22,8},3/1),({16,5,21,9},8/1),({16,5,20,10},16/1),({16,5,19,11},26/1),({16,5,18,12},33/1),({16,5,17,13},34/1),({16,5,16,14},26/1),({16,5,15,15},9/1),({15,6,24,6},1/1),({15,6,23,7},4/1),({15,6,22,8},12/1),({15,6,21,9},25/1),({15,6,20,10},43/1),({15,6,19,11},62/1),({15,6,18,12},76/1),({15,6,17,13},73/1),({15,6,16,14},55/1),({15,6,15,15},20/1),({14,7,25,5},1/1),({14,7,24,6},3/1),({14,7,23,7},11/1),({14,7,22,8},26/1),({14,7,21,9},50/1),({14,7,20,10},82/1),({14,7,19,11},111/1),({14,7,18,12},128/1),({14,7,17,13},123/1),({14,7,16,14},89/1),({14,7,15,15},32/1),({13,8,25,5},2/1),({13,8,24,6},7/1),({13,8,23,7},19/1),({13,8,22,8},41/1),({13,8,21,9},74/1),({13,8,20,10},114/1),({13,8,19,11},150/1),({13,8,18,12},169/1),({13,8,17,13},158/1),({13,8,16,14},114/1),({13,8,15,15},41/1),({12,9,25,5},2/1),({12,9,24,6},8/1),({12,9,23,7},21/1),({12,9,22,8},44/1),({12,9,21,9},78/1),({12,9,20,10},115/1),({12,9,19,11},149/1),({12,9,18,12},166/1),({12,9,17,13},152/1),({12,9,16,14},109/1),({12,9,15,15},40/1),({11,10,26,4},1/1),({11,10,25,5},2/1),({11,10,24,6},6/1),({11,10,23,7},15/1),({11,10,22,8},30/1),({11,10,21,9},50/1),({11,10,20,10},75/1),({11,10,19,11},95/1),({11,10,18,12},104/1),({11,10,17,13},96/1),({11,10,16,14},68/1),({11,10,15,15},23/1)}, (6,2) => {}, (8,0) => {}, (8,1) => {({22,5,22,16},1/1),({22,5,21,17},1/1),({22,5,20,18},1/1),({22,5,19,19},1/1),({21,6,25,13},1/1),({21,6,24,14},2/1),({21,6,23,15},5/1),({21,6,22,16},7/1),({21,6,21,17},8/1),({21,6,20,18},6/1),({21,6,19,19},3/1),({20,7,26,12},2/1),({20,7,25,13},6/1),({20,7,24,14},13/1),({20,7,23,15},22/1),({20,7,22,16},29/1),({20,7,21,17},30/1),({20,7,20,18},23/1),({20,7,19,19},9/1),({19,8,28,10},1/1),({19,8,27,11},3/1),({19,8,26,12},10/1),({19,8,25,13},24/1),({19,8,24,14},43/1),({19,8,23,15},65/1),({19,8,22,16},80/1),({19,8,21,17},80/1),({19,8,20,18},59/1),({19,8,19,19},23/1),({18,9,28,10},3/1),({18,9,27,11},11/1),({18,9,26,12},29/1),({18,9,25,13},59/1),({18,9,24,14},99/1),({18,9,23,15},139/1),({18,9,22,16},164/1),({18,9,21,17},159/1),({18,9,20,18},116/1),({18,9,19,19},43/1),({17,10,29,9},2/1),({17,10,28,10},8/1),({17,10,27,11},25/1),({17,10,26,12},57/1),({17,10,25,13},108/1),({17,10,24,14},170/1),({17,10,23,15},230/1),({17,10,22,16},262/1),({17,10,21,17},250/1),({17,10,20,18},179/1),({17,10,19,19},67/1),({16,11,30,8},1/1),({16,11,29,9},4/1),({16,11,28,10},15/1),({16,11,27,11},39/1),({16,11,26,12},83/1),({16,11,25,13},148/1),({16,11,24,14},225/1),({16,11,23,15},295/1),({16,11,22,16},331/1),({16,11,21,17},309/1),({16,11,20,18},221/1),({16,11,19,19},81/1),({15,12,30,8},1/1),({15,12,29,9},5/1),({15,12,28,10},17/1),({15,12,27,11},42/1),({15,12,26,12},86/1),({15,12,25,13},149/1),({15,12,24,14},221/1),({15,12,23,15},285/1),({15,12,22,16},315/1),({15,12,21,17},292/1),({15,12,20,18},207/1),({15,12,19,19},76/1),({14,13,30,8},1/1),({14,13,29,9},4/1),({14,13,28,10},12/1),({14,13,27,11},28/1),({14,13,26,12},56/1),({14,13,25,13},95/1),({14,13,24,14},139/1),({14,13,23,15},178/1),({14,13,22,16},195/1),({14,13,21,17},179/1),({14,13,20,18},127/1),({14,13,19,19},46/1)}, (8,2) => {}, (10,0) => {}, (10,1) => {({25,8,26,20},1/1),({25,8,25,21},1/1),({25,8,24,22},1/1),({24,9,29,17},1/1),({24,9,28,18},2/1),({24,9,27,19},4/1),({24,9,26,20},6/1),({24,9,25,21},7/1),({24,9,24,22},6/1),({24,9,23,23},2/1),({23,10,30,16},1/1),({23,10,29,17},4/1),({23,10,28,18},11/1),({23,10,27,19},18/1),({23,10,26,20},25/1),({23,10,25,21},26/1),({23,10,24,22},21/1),({23,10,23,23},7/1),({22,11,31,15},2/1),({22,11,30,16},7/1),({22,11,29,17},18/1),({22,11,28,18},35/1),({22,11,27,19},54/1),({22,11,26,20},68/1),({22,11,25,21},69/1),({22,11,24,22},52/1),({22,11,23,23},19/1),({21,12,32,14},1/1),({21,12,31,15},7/1),({21,12,30,16},21/1),({21,12,29,17},45/1),({21,12,28,18},80/1),({21,12,27,19},115/1),({21,12,26,20},140/1),({21,12,25,21},136/1),({21,12,24,22},102/1),({21,12,23,23},36/1),({20,13,33,13},1/1),({20,13,32,14},5/1),({20,13,31,15},17/1),({20,13,30,16},42/1),({20,13,29,17},84/1),({20,13,28,18},138/1),({20,13,27,19},191/1),({20,13,26,20},223/1),({20,13,25,21},214/1),({20,13,24,22},156/1),({20,13,23,23},57/1),({19,14,33,13},2/1),({19,14,32,14},9/1),({19,14,31,15},27/1),({19,14,30,16},62/1),({19,14,29,17},116/1),({19,14,28,18},183/1),({19,14,27,19},245/1),({19,14,26,20},280/1),({19,14,25,21},264/1),({19,14,24,22},192/1),({19,14,23,23},69/1),({18,15,33,13},3/1),({18,15,32,14},11/1),({18,15,31,15},30/1),({18,15,30,16},65/1),({18,15,29,17},117/1),({18,15,28,18},180/1),({18,15,27,19},237/1),({18,15,26,20},267/1),({18,15,25,21},250/1),({18,15,24,22},179/1),({18,15,23,23},65/1),({17,16,34,12},1/1),({17,16,33,13},2/1),({17,16,32,14},8/1),({17,16,31,15},20/1),({17,16,30,16},43/1),({17,16,29,17},75/1),({17,16,28,18},114/1),({17,16,27,19},147/1),({17,16,26,20},166/1),({17,16,25,21},153/1),({17,16,24,22},110/1),({17,16,23,23},39/1)}, (10,2) => {}, (12,0) => {}, (12,1) => {({27,12,31,23},1/1),({27,12,30,24},1/1),({27,12,29,25},2/1),({27,12,28,26},1/1),({27,12,27,27},1/1),({26,13,32,22},2/1),({26,13,31,23},4/1),({26,13,30,24},6/1),({26,13,29,25},8/1),({26,13,28,26},6/1),({26,13,27,27},2/1),({25,14,34,20},1/1),({25,14,33,21},3/1),({25,14,32,22},8/1),({25,14,31,23},15/1),({25,14,30,24},20/1),({25,14,29,25},23/1),({25,14,28,26},17/1),({25,14,27,27},7/1),({24,15,34,20},3/1),({24,15,33,21},10/1),({24,15,32,22},21/1),({24,15,31,23},35/1),({24,15,30,24},46/1),({24,15,29,25},47/1),({24,15,28,26},36/1),({24,15,27,27},14/1),({23,16,35,19},3/1),({23,16,34,20},9/1),({23,16,33,21},22/1),({23,16,32,22},41/1),({23,16,31,23},63/1),({23,16,30,24},76/1),({23,16,29,25},79/1),({23,16,28,26},57/1),({23,16,27,27},22/1),({22,17,36,18},1/1),({22,17,35,19},5/1),({22,17,34,20},15/1),({22,17,33,21},33/1),({22,17,32,22},58/1),({22,17,31,23},84/1),({22,17,30,24},101/1),({22,17,29,25},99/1),({22,17,28,26},73/1),({22,17,27,27},27/1),({21,18,36,18},2/1),({21,18,35,19},6/1),({21,18,34,20},17/1),({21,18,33,21},36/1),({21,18,32,22},59/1),({21,18,31,23},84/1),({21,18,30,24},99/1),({21,18,29,25},95/1),({21,18,28,26},69/1),({21,18,27,27},27/1),({20,19,36,18},1/1),({20,19,35,19},5/1),({20,19,34,20},12/1),({20,19,33,21},23/1),({20,19,32,22},39/1),({20,19,31,23},53/1),({20,19,30,24},61/1),({20,19,29,25},60/1),({20,19,28,26},43/1),({20,19,27,27},15/1)}, (12,2) => {}, (14,0) => {}, (14,1) => {({29,16,32,30},1/1),({28,17,35,27},1/1),({28,17,34,28},2/1),({28,17,33,29},2/1),({28,17,32,30},2/1),({28,17,31,31},1/1),({27,18,36,26},1/1),({27,18,35,27},2/1),({27,18,34,28},5/1),({27,18,33,29},5/1),({27,18,32,30},5/1),({27,18,31,31},1/1),({26,19,37,25},1/1),({26,19,36,26},3/1),({26,19,35,27},6/1),({26,19,34,28},9/1),({26,19,33,29},10/1),({26,19,32,30},8/1),({26,19,31,31},3/1),({25,20,37,25},2/1),({25,20,36,26},5/1),({25,20,35,27},9/1),({25,20,34,28},13/1),({25,20,33,29},13/1),({25,20,32,30},11/1),({25,20,31,31},4/1),({24,21,38,24},1/1),({24,21,37,25},3/1),({24,21,36,26},6/1),({24,21,35,27},10/1),({24,21,34,28},13/1),({24,21,33,29},13/1),({24,21,32,30},10/1),({24,21,31,31},4/1),({23,22,38,24},1/1),({23,22,37,25},2/1),({23,22,36,26},4/1),({23,22,35,27},6/1),({23,22,34,28},9/1),({23,22,33,29},8/1),({23,22,32,30},7/1),({23,22,31,31},2/1)}, (14,2) => {}, (16,0) => {}, (16,1) => {({29,22,37,33},1/1),({29,22,35,35},1/1)}, (16,2) => {({28,26,39,35},1/1)}, (18,2) => {}, (1,0) => {({3,0,5,1},1/1),({3,0,4,2},1/1)}, (1,1) => {({4,2,10,0},1/1)}, (3,0) => {}, (3,1) => {({12,0,9,9},1/1),({11,1,12,6},1/1),({11,1,11,7},1/1),({11,1,10,8},1/1),({10,2,15,3},1/1),({10,2,14,4},1/1),({10,2,13,5},3/1),({10,2,12,6},2/1),({10,2,11,7},4/1),({10,2,10,8},1/1),({10,2,9,9},2/1),({9,3,16,2},1/1),({9,3,15,3},2/1),({9,3,14,4},4/1),({9,3,13,5},5/1),({9,3,12,6},7/1),({9,3,11,7},5/1),({9,3,10,8},5/1),({9,3,9,9},1/1),({8,4,17,1},1/1),({8,4,16,2},2/1),({8,4,15,3},4/1),({8,4,14,4},6/1),({8,4,13,5},9/1),({8,4,12,6},7/1),({8,4,11,7},9/1),({8,4,10,8},4/1),({8,4,9,9},3/1),({7,5,17,1},1/1),({7,5,16,2},3/1),({7,5,15,3},4/1),({7,5,14,4},7/1),({7,5,13,5},7/1),({7,5,12,6},9/1),({7,5,11,7},6/1),({7,5,10,8},6/1),({6,6,15,3},2/1),({6,6,14,4},1/1),({6,6,13,5},4/1),({6,6,12,6},2/1),({6,6,11,7},4/1),({6,6,9,9},3/1)}, (5,0) => {}, (3,2) => {}, (5,1) => {({16,2,16,10},2/1),({16,2,15,11},1/1),({16,2,14,12},2/1),({15,3,19,7},1/1),({15,3,18,8},2/1),({15,3,17,9},5/1),({15,3,16,10},7/1),({15,3,15,11},8/1),({15,3,14,12},6/1),({15,3,13,13},3/1),({14,4,20,6},2/1),({14,4,19,7},5/1),({14,4,18,8},11/1),({14,4,17,9},16/1),({14,4,16,10},23/1),({14,4,15,11},21/1),({14,4,14,12},19/1),({14,4,13,13},5/1),({13,5,22,4},1/1),({13,5,21,5},3/1),({13,5,20,6},7/1),({13,5,19,7},16/1),({13,5,18,8},27/1),({13,5,17,9},40/1),({13,5,16,10},47/1),({13,5,15,11},47/1),({13,5,14,12},33/1),({13,5,13,13},14/1),({12,6,22,4},2/1),({12,6,21,5},6/1),({12,6,20,6},16/1),({12,6,19,7},29/1),({12,6,18,8},49/1),({12,6,17,9},63/1),({12,6,16,10},77/1),({12,6,15,11},69/1),({12,6,14,12},54/1),({12,6,13,13},16/1),({11,7,23,3},1/1),({11,7,22,4},3/1),({11,7,21,5},10/1),({11,7,20,6},21/1),({11,7,19,7},40/1),({11,7,18,8},59/1),({11,7,17,9},81/1),({11,7,16,10},87/1),({11,7,15,11},85/1),({11,7,14,12},57/1),({11,7,13,13},24/1),({10,8,23,3},1/1),({10,8,22,4},4/1),({10,8,21,5},9/1),({10,8,20,6},21/1),({10,8,19,7},34/1),({10,8,18,8},54/1),({10,8,17,9},65/1),({10,8,16,10},77/1),({10,8,15,11},65/1),({10,8,14,12},52/1),({10,8,13,13},14/1),({9,9,22,4},1/1),({9,9,21,5},4/1),({9,9,20,6},7/1),({9,9,19,7},15/1),({9,9,18,8},19/1),({9,9,17,9},28/1),({9,9,16,10},26/1),({9,9,15,11},29/1),({9,9,14,12},16/1),({9,9,13,13},10/1)}};
--dw stands for dominant weights
<<<<<<< HEAD
<<<<<<< HEAD
dw0234 = new HashTable from {(5,2) => {}, (7,0) => {}, (7,1) => {({20,4,21,13},1/1),({19,5,23,11},1/1),({18,6,25,9},1/1),({17,7,26,8},1/1),({16,8,27,7},1/1),({14,10,28,6},1/1)}, (9,0) => {}, (7,2) => {}, (9,1) => {({24,6,22,20},1/1),({23,7,26,16},1/1),({22,8,28,14},1/1),({21,9,29,13},1/1),({20,10,30,12},1/1),({19,11,31,11},1/1),({17,13,32,10},1/1)}, (9,2) => {}, (11,0) => {}, (11,1) => {({26,10,29,21},1/1),({25,11,31,19},1/1),({24,12,32,18},1/1),({23,13,33,17},1/1),({22,14,34,16},1/1),({20,16,35,15},1/1)}, (11,2) => {}, (13,0) => {}, (13,1) => {({28,14,32,26},1/1),({27,15,34,24},1/1),({26,16,35,23},1/1),({25,17,36,22},1/1),({23,19,37,21},2/1)}, (15,0) => {}, (13,2) => {}, (15,1) => {({29,19,35,31},1/1),({28,20,37,29},1/1),({26,22,38,28},1/1),({24,24,39,27},1/1)}, (15,2) => {}, (17,0) => {}, (17,1) => {}, (17,2) => {({29,28,39,39},1/1)}, (0,0) => {({0,0,2,0},1/1)}, (2,0) => {({6,0,7,3},1/1),({5,1,8,2},1/1)}, (2,1) => {({7,2,13,1},1/1),({5,4,14,0},1/1)}, (4,0) => {}, (2,2) => {}, (4,1) => {({14,1,13,9},1/1),({13,2,16,6},1/1),({12,3,18,4},1/1),({11,4,19,3},1/1),({10,5,20,2},1/1)}, (4,2) => {}, (6,0) => {}, (6,1) => {({18,3,19,11},1/1),({17,4,21,9},1/1),({16,5,23,7},1/1),({15,6,24,6},1/1),({14,7,25,5},1/1),({11,10,26,4},1/1)}, (6,2) => {}, (8,0) => {}, (8,1) => {({22,5,22,16},1/1),({21,6,25,13},1/1),({20,7,26,12},2/1),({19,8,28,10},1/1),({17,10,29,9},2/1),({16,11,30,8},1/1)}, (8,2) => {}, (10,0) => {}, (10,1) => {({25,8,26,20},1/1),({24,9,29,17},1/1),({23,10,30,16},1/1),({22,11,31,15},2/1),({21,12,32,14},1/1),({20,13,33,13},1/1),({17,16,34,12},1/1)}, (10,2) => {}, (12,0) => {}, (12,1) => {({27,12,31,23},1/1),({26,13,32,22},2/1),({25,14,34,20},1/1),({23,16,35,19},3/1),({22,17,36,18},1/1)}, (12,2) => {}, (14,0) => {}, (14,1) => {({29,16,32,30},1/1),({28,17,35,27},1/1),({27,18,36,26},1/1),({26,19,37,25},1/1),({24,21,38,24},1/1)}, (14,2) => {}, (16,0) => {}, (16,1) => {({29,22,37,33},1/1)}, (16,2) => {({28,26,39,35},1/1)}, (18,2) => {}, (1,0) => {({3,0,5,1},1/1)}, (1,1) => {({4,2,10,0},1/1)}, (3,0) => {}, (3,1) => {({12,0,9,9},1/1),({11,1,12,6},1/1),({10,2,15,3},1/1),({9,3,16,2},1/1),({8,4,17,1},1/1)}, (5,0) => {}, (3,2) => {}, (5,1) => {({16,2,16,10},2/1),({15,3,19,7},1/1),({14,4,20,6},2/1),({13,5,22,4},1/1),({11,7,23,3},1/1)}};
--dmw stands for dominant monomial weights
dmw0234 = {infinity};
=======
dw0234 = new HashTable from {(5,2) => {}, (7,0) => {}, (7,1) => {({20,4,21,13},1/1),({19,5,23,11},1/1),({18,6,25,9},1/1),({17,7,26,8},1/1),({16,8,27,7},1/1),({14,10,28,6},1/1)}, (7,2) => {}, (9,0) => {}, (9,1) => {({24,6,22,20},1/1),({23,7,26,16},1/1),({22,8,28,14},1/1),({21,9,29,13},1/1),({20,10,30,12},1/1),({19,11,31,11},1/1),({17,13,32,10},1/1)}, (9,2) => {}, (11,0) => {}, (11,1) => {({26,10,29,21},1/1),({25,11,31,19},1/1),({24,12,32,18},1/1),({23,13,33,17},1/1),({22,14,34,16},1/1),({20,16,35,15},1/1)}, (11,2) => {}, (13,0) => {}, (13,1) => {({28,14,32,26},1/1),({27,15,34,24},1/1),({26,16,35,23},1/1),({25,17,36,22},1/1),({23,19,37,21},2/1)}, (13,2) => {}, (15,0) => {}, (15,1) => {({29,19,35,31},1/1),({28,20,37,29},1/1),({26,22,38,28},1/1),({24,24,39,27},1/1)}, (17,0) => {}, (15,2) => {}, (17,1) => {}, (17,2) => {({29,28,39,39},1/1)}, (0,0) => {({0,0,2,0},1/1)}, (2,0) => {({6,0,7,3},1/1),({5,1,8,2},1/1)}, (2,1) => {({7,2,13,1},1/1),({5,4,14,0},1/1)}, (2,2) => {}, (4,0) => {}, (4,1) => {({14,1,13,9},1/1),({13,2,16,6},1/1),({12,3,18,4},1/1),({11,4,19,3},1/1),({10,5,20,2},1/1)}, (4,2) => {}, (6,0) => {}, (6,1) => {({18,3,19,11},1/1),({17,4,21,9},1/1),({16,5,23,7},1/1),({15,6,24,6},1/1),({14,7,25,5},1/1),({11,10,26,4},1/1)}, (6,2) => {}, (8,0) => {}, (8,1) => {({22,5,22,16},1/1),({21,6,25,13},1/1),({20,7,26,12},2/1),({19,8,28,10},1/1),({17,10,29,9},2/1),({16,11,30,8},1/1)}, (8,2) => {}, (10,0) => {}, (10,1) => {({25,8,26,20},1/1),({24,9,29,17},1/1),({23,10,30,16},1/1),({22,11,31,15},2/1),({21,12,32,14},1/1),({20,13,33,13},1/1),({17,16,34,12},1/1)}, (10,2) => {}, (12,0) => {}, (12,1) => {({27,12,31,23},1/1),({26,13,32,22},2/1),({25,14,34,20},1/1),({23,16,35,19},3/1),({22,17,36,18},1/1)}, (12,2) => {}, (14,0) => {}, (14,1) => {({29,16,32,30},1/1),({28,17,35,27},1/1),({27,18,36,26},1/1),({26,19,37,25},1/1),({24,21,38,24},1/1)}, (14,2) => {}, (16,0) => {}, (16,1) => {({29,22,37,33},1/1)}, (16,2) => {({28,26,39,35},1/1)}, (18,2) => {}, (1,0) => {({3,0,5,1},1/1)}, (1,1) => {({4,2,10,0},1/1)}, (3,0) => {}, (3,1) => {({12,0,9,9},1/1),({11,1,12,6},1/1),({10,2,15,3},1/1),({9,3,16,2},1/1),({8,4,17,1},1/1)}, (3,2) => {}, (5,0) => {}, (5,1) => {({16,2,16,10},2/1),({15,3,19,7},1/1),({14,4,20,6},2/1),({13,5,22,4},1/1),({11,7,23,3},1/1)}};
>>>>>>> parent of 67e86b7... updating data
=======
dw0234 = new HashTable from {(5,2) => {}, (7,0) => {}, (7,1) => {({20,4,21,13},1/1),({19,5,23,11},1/1),({18,6,25,9},1/1),({17,7,26,8},1/1),({16,8,27,7},1/1),({14,10,28,6},1/1)}, (7,2) => {}, (9,0) => {}, (9,1) => {({24,6,22,20},1/1),({23,7,26,16},1/1),({22,8,28,14},1/1),({21,9,29,13},1/1),({20,10,30,12},1/1),({19,11,31,11},1/1),({17,13,32,10},1/1)}, (9,2) => {}, (11,0) => {}, (11,1) => {({26,10,29,21},1/1),({25,11,31,19},1/1),({24,12,32,18},1/1),({23,13,33,17},1/1),({22,14,34,16},1/1),({20,16,35,15},1/1)}, (11,2) => {}, (13,0) => {}, (13,1) => {({28,14,32,26},1/1),({27,15,34,24},1/1),({26,16,35,23},1/1),({25,17,36,22},1/1),({23,19,37,21},2/1)}, (13,2) => {}, (15,0) => {}, (15,1) => {({29,19,35,31},1/1),({28,20,37,29},1/1),({26,22,38,28},1/1),({24,24,39,27},1/1)}, (17,0) => {}, (15,2) => {}, (17,1) => {}, (17,2) => {({29,28,39,39},1/1)}, (0,0) => {({0,0,2,0},1/1)}, (2,0) => {({6,0,7,3},1/1),({5,1,8,2},1/1)}, (2,1) => {({7,2,13,1},1/1),({5,4,14,0},1/1)}, (2,2) => {}, (4,0) => {}, (4,1) => {({14,1,13,9},1/1),({13,2,16,6},1/1),({12,3,18,4},1/1),({11,4,19,3},1/1),({10,5,20,2},1/1)}, (4,2) => {}, (6,0) => {}, (6,1) => {({18,3,19,11},1/1),({17,4,21,9},1/1),({16,5,23,7},1/1),({15,6,24,6},1/1),({14,7,25,5},1/1),({11,10,26,4},1/1)}, (6,2) => {}, (8,0) => {}, (8,1) => {({22,5,22,16},1/1),({21,6,25,13},1/1),({20,7,26,12},2/1),({19,8,28,10},1/1),({17,10,29,9},2/1),({16,11,30,8},1/1)}, (8,2) => {}, (10,0) => {}, (10,1) => {({25,8,26,20},1/1),({24,9,29,17},1/1),({23,10,30,16},1/1),({22,11,31,15},2/1),({21,12,32,14},1/1),({20,13,33,13},1/1),({17,16,34,12},1/1)}, (10,2) => {}, (12,0) => {}, (12,1) => {({27,12,31,23},1/1),({26,13,32,22},2/1),({25,14,34,20},1/1),({23,16,35,19},3/1),({22,17,36,18},1/1)}, (12,2) => {}, (14,0) => {}, (14,1) => {({29,16,32,30},1/1),({28,17,35,27},1/1),({27,18,36,26},1/1),({26,19,37,25},1/1),({24,21,38,24},1/1)}, (14,2) => {}, (16,0) => {}, (16,1) => {({29,22,37,33},1/1)}, (16,2) => {({28,26,39,35},1/1)}, (18,2) => {}, (1,0) => {({3,0,5,1},1/1)}, (1,1) => {({4,2,10,0},1/1)}, (3,0) => {}, (3,1) => {({12,0,9,9},1/1),({11,1,12,6},1/1),({10,2,15,3},1/1),({9,3,16,2},1/1),({8,4,17,1},1/1)}, (3,2) => {}, (5,0) => {}, (5,1) => {({16,2,16,10},2/1),({15,3,19,7},1/1),({14,4,20,6},2/1),({13,5,22,4},1/1),({11,7,23,3},1/1)}};
>>>>>>> parent of 67e86b7... updating data
end;