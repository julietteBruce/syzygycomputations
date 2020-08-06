A := QQ[t_0,t_1,t_2,t_3];
--tb stands for Total Betti numbers
tb1226 = new HashTable from {(5,2) => 0, (7,0) => 0, (7,1) => 310284/1, (7,2) => 0, (9,0) => 0, (9,1) => 612612/1, (9,2) => 0, (11,0) => 0, (11,1) => 461448/1, (11,2) => 0, (13,0) => 0, (13,1) => 135864/1, (15,0) => 0, (13,2) => 0, (15,1) => 13770/1, (17,0) => 0, (15,2) => 0, (17,1) => 318/1, (17,2) => 0, (19,1) => 0, (0,0) => 6/1, (2,0) => 594/1, (2,1) => 0, (2,2) => 0, (4,0) => 3822/1, (4,1) => 5674/1, (4,2) => 0, (6,0) => 0, (6,1) => 143208/1, (6,2) => 0, (8,0) => 0, (8,1) => 495924/1, (8,2) => 0, (10,0) => 0, (10,1) => 596700/1, (10,2) => 0, (12,0) => 0, (12,1) => 282744/1, (12,2) => 0, (14,0) => 0, (14,1) => 50184/1, (16,0) => 0, (14,2) => 0, (16,1) => 2646/1, (18,0) => 0, (16,2) => 0, (18,1) => 18/1, (18,2) => 0, (20,1) => 0, (1,0) => 90/1, (1,1) => 0, (3,0) => 2142/1, (3,1) => 150/1, (3,2) => 0, (5,0) => 2002/1, (5,1) => 42840/1};
--sb represents the betti numbers as sums of Schur functors
sb1226 = new HashTable from {(5,2) => {}, (7,0) => {}, (7,1) => {({15,2,31,19},1/1),({15,2,30,20},1/1),({15,2,29,21},2/1),({15,2,28,22},1/1),({15,2,27,23},2/1),({15,2,26,24},1/1),({15,2,25,25},1/1),({14,3,34,16},1/1),({14,3,33,17},2/1),({14,3,32,18},4/1),({14,3,31,19},7/1),({14,3,30,20},10/1),({14,3,29,21},12/1),({14,3,28,22},12/1),({14,3,27,23},11/1),({14,3,26,24},8/1),({14,3,25,25},3/1),({13,4,36,14},1/1),({13,4,35,15},4/1),({13,4,34,16},8/1),({13,4,33,17},15/1),({13,4,32,18},23/1),({13,4,31,19},34/1),({13,4,30,20},41/1),({13,4,29,21},48/1),({13,4,28,22},46/1),({13,4,27,23},42/1),({13,4,26,24},26/1),({13,4,25,25},11/1),({12,5,38,12},1/1),({12,5,37,13},3/1),({12,5,36,14},8/1),({12,5,35,15},17/1),({12,5,34,16},31/1),({12,5,33,17},50/1),({12,5,32,18},72/1),({12,5,31,19},94/1),({12,5,30,20},112/1),({12,5,29,21},121/1),({12,5,28,22},117/1),({12,5,27,23},98/1),({12,5,26,24},65/1),({12,5,25,25},23/1),({11,6,39,11},1/1),({11,6,38,12},4/1),({11,6,37,13},11/1),({11,6,36,14},22/1),({11,6,35,15},42/1),({11,6,34,16},68/1),({11,6,33,17},104/1),({11,6,32,18},139/1),({11,6,31,19},177/1),({11,6,30,20},201/1),({11,6,29,21},215/1),({11,6,28,22},200/1),({11,6,27,23},169/1),({11,6,26,24},108/1),({11,6,25,25},41/1),({10,7,40,10},1/1),({10,7,39,11},3/1),({10,7,38,12},8/1),({10,7,37,13},18/1),({10,7,36,14},35/1),({10,7,35,15},60/1),({10,7,34,16},94/1),({10,7,33,17},135/1),({10,7,32,18},179/1),({10,7,31,19},219/1),({10,7,30,20},247/1),({10,7,29,21},256/1),({10,7,28,22},239/1),({10,7,27,23},196/1),({10,7,26,24},129/1),({10,7,25,25},45/1),({9,8,40,10},1/1),({9,8,39,11},3/1),({9,8,38,12},7/1),({9,8,37,13},15/1),({9,8,36,14},28/1),({9,8,35,15},47/1),({9,8,34,16},71/1),({9,8,33,17},100/1),({9,8,32,18},130/1),({9,8,31,19},158/1),({9,8,30,20},175/1),({9,8,29,21},181/1),({9,8,28,22},167/1),({9,8,27,23},138/1),({9,8,26,24},89/1),({9,8,25,25},32/1)}, (7,2) => {}, (9,0) => {}, (9,1) => {({17,4,38,24},1/1),({17,4,37,25},1/1),({17,4,36,26},2/1),({17,4,35,27},2/1),({17,4,34,28},3/1),({17,4,33,29},2/1),({17,4,32,30},2/1),({16,5,41,21},1/1),({16,5,40,22},2/1),({16,5,39,23},5/1),({16,5,38,24},9/1),({16,5,37,25},13/1),({16,5,36,26},17/1),({16,5,35,27},20/1),({16,5,34,28},21/1),({16,5,33,29},18/1),({16,5,32,30},12/1),({16,5,31,31},4/1),({15,6,43,19},1/1),({15,6,42,20},4/1),({15,6,41,21},9/1),({15,6,40,22},18/1),({15,6,39,23},30/1),({15,6,38,24},47/1),({15,6,37,25},62/1),({15,6,36,26},77/1),({15,6,35,27},83/1),({15,6,34,28},84/1),({15,6,33,29},69/1),({15,6,32,30},48/1),({15,6,31,31},15/1),({14,7,45,17},1/1),({14,7,44,18},3/1),({14,7,43,19},9/1),({14,7,42,20},20/1),({14,7,41,21},39/1),({14,7,40,22},66/1),({14,7,39,23},101/1),({14,7,38,24},141/1),({14,7,37,25},180/1),({14,7,36,26},210/1),({14,7,35,27},223/1),({14,7,34,28},213/1),({14,7,33,29},177/1),({14,7,32,30},118/1),({14,7,31,31},41/1),({13,8,46,16},1/1),({13,8,45,17},4/1),({13,8,44,18},12/1),({13,8,43,19},26/1),({13,8,42,20},52/1),({13,8,41,21},90/1),({13,8,40,22},144/1),({13,8,39,23},206/1),({13,8,38,24},277/1),({13,8,37,25},338/1),({13,8,36,26},387/1),({13,8,35,27},399/1),({13,8,34,28},378/1),({13,8,33,29},307/1),({13,8,32,30},206/1),({13,8,31,31},69/1),({12,9,47,15},1/1),({12,9,46,16},3/1),({12,9,45,17},9/1),({12,9,44,18},21/1),({12,9,43,19},43/1),({12,9,42,20},78/1),({12,9,41,21},129/1),({12,9,40,22},195/1),({12,9,39,23},273/1),({12,9,38,24},354/1),({12,9,37,25},426/1),({12,9,36,26},475/1),({12,9,35,27},487/1),({12,9,34,28},453/1),({12,9,33,29},369/1),({12,9,32,30},242/1),({12,9,31,31},84/1),({11,10,47,15},1/1),({11,10,46,16},3/1),({11,10,45,17},8/1),({11,10,44,18},18/1),({11,10,43,19},35/1),({11,10,42,20},62/1),({11,10,41,21},99/1),({11,10,40,22},147/1),({11,10,39,23},202/1),({11,10,38,24},259/1),({11,10,37,25},307/1),({11,10,36,26},340/1),({11,10,35,27},345/1),({11,10,34,28},321/1),({11,10,33,29},259/1),({11,10,32,30},170/1),({11,10,31,31},58/1)}, (9,2) => {}, (11,0) => {}, (11,1) => {({19,6,41,33},1/1),({19,6,40,34},1/1),({19,6,39,35},1/1),({18,7,46,28},1/1),({18,7,45,29},2/1),({18,7,44,30},4/1),({18,7,43,31},6/1),({18,7,42,32},8/1),({18,7,41,33},10/1),({18,7,40,34},11/1),({18,7,39,35},10/1),({18,7,38,36},6/1),({18,7,37,37},2/1),({17,8,49,25},1/1),({17,8,48,26},2/1),({17,8,47,27},5/1),({17,8,46,28},10/1),({17,8,45,29},19/1),({17,8,44,30},28/1),({17,8,43,31},40/1),({17,8,42,32},48/1),({17,8,41,33},55/1),({17,8,40,34},53/1),({17,8,39,35},47/1),({17,8,38,36},30/1),({17,8,37,37},12/1),({16,9,50,24},2/1),({16,9,49,25},5/1),({16,9,48,26},13/1),({16,9,47,27},26/1),({16,9,46,28},45/1),({16,9,45,29},70/1),({16,9,44,30},100/1),({16,9,43,31},129/1),({16,9,42,32},152/1),({16,9,41,33},163/1),({16,9,40,34},156/1),({16,9,39,35},131/1),({16,9,38,36},87/1),({16,9,37,37},31/1),({15,10,52,22},1/1),({15,10,51,23},3/1),({15,10,50,24},8/1),({15,10,49,25},19/1),({15,10,48,26},38/1),({15,10,47,27},68/1),({15,10,46,28},108/1),({15,10,45,29},159/1),({15,10,44,30},213/1),({15,10,43,31},266/1),({15,10,42,32},302/1),({15,10,41,33},317/1),({15,10,40,34},297/1),({15,10,39,35},247/1),({15,10,38,36},161/1),({15,10,37,37},58/1),({14,11,52,22},2/1),({14,11,51,23},6/1),({14,11,50,24},15/1),({14,11,49,25},32/1),({14,11,48,26},60/1),({14,11,47,27},100/1),({14,11,46,28},155/1),({14,11,45,29},219/1),({14,11,44,30},286/1),({14,11,43,31},347/1),({14,11,42,32},389/1),({14,11,41,33},401/1),({14,11,40,34},374/1),({14,11,39,35},306/1),({14,11,38,36},200/1),({14,11,37,37},70/1),({13,12,53,21},1/1),({13,12,52,22},2/1),({13,12,51,23},6/1),({13,12,50,24},14/1),({13,12,49,25},28/1),({13,12,48,26},49/1),({13,12,47,27},81/1),({13,12,46,28},120/1),({13,12,45,29},168/1),({13,12,44,30},215/1),({13,12,43,31},258/1),({13,12,42,32},285/1),({13,12,41,33},294/1),({13,12,40,34},270/1),({13,12,39,35},221/1),({13,12,38,36},143/1),({13,12,37,37},51/1)}, (11,2) => {}, (13,0) => {}, (13,1) => {({20,9,47,39},1/1),({20,9,46,40},1/1),({20,9,45,41},1/1),({20,9,44,42},1/1),({19,10,52,34},1/1),({19,10,51,35},2/1),({19,10,50,36},4/1),({19,10,49,37},6/1),({19,10,48,38},9/1),({19,10,47,39},10/1),({19,10,46,40},11/1),({19,10,45,41},9/1),({19,10,44,42},7/1),({19,10,43,43},2/1),({18,11,55,31},1/1),({18,11,54,32},2/1),({18,11,53,33},5/1),({18,11,52,34},9/1),({18,11,51,35},16/1),({18,11,50,36},24/1),({18,11,49,37},33/1),({18,11,48,38},41/1),({18,11,47,39},45/1),({18,11,46,40},45/1),({18,11,45,41},38/1),({18,11,44,42},26/1),({18,11,43,43},9/1),({17,12,56,30},1/1),({17,12,55,31},3/1),({17,12,54,32},8/1),({17,12,53,33},16/1),({17,12,52,34},29/1),({17,12,51,35},45/1),({17,12,50,36},65/1),({17,12,49,37},83/1),({17,12,48,38},99/1),({17,12,47,39},105/1),({17,12,46,40},102/1),({17,12,45,41},84/1),({17,12,44,42},57/1),({17,12,43,43},19/1),({16,13,57,29},1/1),({16,13,56,30},3/1),({16,13,55,31},8/1),({16,13,54,32},16/1),({16,13,53,33},30/1),({16,13,52,34},49/1),({16,13,51,35},73/1),({16,13,50,36},99/1),({16,13,49,37},124/1),({16,13,48,38},143/1),({16,13,47,39},150/1),({16,13,46,40},142/1),({16,13,45,41},117/1),({16,13,44,42},78/1),({16,13,43,43},27/1),({15,14,57,29},1/1),({15,14,56,30},3/1),({15,14,55,31},7/1),({15,14,54,32},14/1),({15,14,53,33},25/1),({15,14,52,34},41/1),({15,14,51,35},59/1),({15,14,50,36},79/1),({15,14,49,37},97/1),({15,14,48,38},111/1),({15,14,47,39},115/1),({15,14,46,40},108/1),({15,14,45,41},88/1),({15,14,44,42},59/1),({15,14,43,43},20/1)}, (15,0) => {}, (13,2) => {}, (15,1) => {({20,13,56,42},1/1),({20,13,55,43},1/1),({20,13,54,44},2/1),({20,13,53,45},3/1),({20,13,52,46},3/1),({20,13,51,47},3/1),({20,13,50,48},2/1),({20,13,49,49},1/1),({19,14,59,39},1/1),({19,14,58,40},2/1),({19,14,57,41},4/1),({19,14,56,42},6/1),({19,14,55,43},9/1),({19,14,54,44},11/1),({19,14,53,45},13/1),({19,14,52,46},12/1),({19,14,51,47},11/1),({19,14,50,48},7/1),({19,14,49,49},3/1),({18,15,60,38},1/1),({18,15,59,39},2/1),({18,15,58,40},5/1),({18,15,57,41},8/1),({18,15,56,42},12/1),({18,15,55,43},17/1),({18,15,54,44},21/1),({18,15,53,45},23/1),({18,15,52,46},22/1),({18,15,51,47},19/1),({18,15,50,48},13/1),({18,15,49,49},5/1),({17,16,60,38},1/1),({17,16,59,39},3/1),({17,16,58,40},5/1),({17,16,57,41},8/1),({17,16,56,42},12/1),({17,16,55,43},16/1),({17,16,54,44},19/1),({17,16,53,45},20/1),({17,16,52,46},19/1),({17,16,51,47},16/1),({17,16,50,48},11/1),({17,16,49,49},4/1)}, (17,0) => {}, (15,2) => {}, (17,1) => {({20,17,61,49},1/1),({20,17,60,50},1/1),({20,17,59,51},1/1),({20,17,58,52},1/1),({20,17,57,53},1/1),({20,17,56,54},1/1),({19,18,62,48},1/1),({19,18,61,49},1/1),({19,18,60,50},1/1),({19,18,59,51},1/1),({19,18,58,52},1/1),({19,18,57,53},1/1),({19,18,56,54},1/1)}, (17,2) => {}, (19,1) => {}, (0,0) => {({1,0,2,0},1/1)}, (2,0) => {({5,0,11,3},1/1),({5,0,9,5},1/1),({5,0,7,7},1/1),({4,1,13,1},1/1),({4,1,12,2},2/1),({4,1,11,3},2/1),({4,1,10,4},2/1),({4,1,9,5},2/1),({4,1,8,6},2/1),({4,1,7,7},1/1),({3,2,13,1},1/1),({3,2,12,2},2/1),({3,2,11,3},2/1),({3,2,10,4},2/1),({3,2,9,5},2/1),({3,2,8,6},2/1),({3,2,7,7},1/1)}, (2,1) => {}, (2,2) => {}, (4,0) => {({7,2,22,4},1/1),({7,2,21,5},1/1),({7,2,20,6},2/1),({7,2,19,7},2/1),({7,2,18,8},4/1),({7,2,17,9},3/1),({7,2,16,10},4/1),({7,2,15,11},2/1),({7,2,14,12},3/1),({6,3,22,4},1/1),({6,3,21,5},2/1),({6,3,20,6},3/1),({6,3,19,7},5/1),({6,3,18,8},7/1),({6,3,17,9},8/1),({6,3,16,10},8/1),({6,3,15,11},7/1),({6,3,14,12},5/1),({6,3,13,13},2/1),({5,4,23,3},1/1),({5,4,22,4},1/1),({5,4,21,5},2/1),({5,4,20,6},4/1),({5,4,19,7},5/1),({5,4,18,8},7/1),({5,4,17,9},7/1),({5,4,16,10},8/1),({5,4,15,11},6/1),({5,4,14,12},5/1),({5,4,13,13},1/1)}, (4,1) => {({11,0,19,13},1/1),({11,0,18,14},1/1),({10,1,23,9},1/1),({10,1,22,10},1/1),({10,1,21,11},2/1),({10,1,20,12},2/1),({10,1,19,13},3/1),({10,1,18,14},2/1),({10,1,17,15},2/1),({9,2,23,9},1/1),({9,2,22,10},2/1),({9,2,21,11},3/1),({9,2,20,12},4/1),({9,2,19,13},4/1),({9,2,18,14},4/1),({9,2,17,15},3/1),({9,2,16,16},1/1),({8,3,26,6},1/1),({8,3,25,7},1/1),({8,3,24,8},2/1),({8,3,23,9},3/1),({8,3,22,10},4/1),({8,3,21,11},5/1),({8,3,20,12},5/1),({8,3,19,13},5/1),({8,3,18,14},4/1),({8,3,17,15},3/1),({8,3,16,16},1/1),({7,4,25,7},1/1),({7,4,24,8},1/1),({7,4,23,9},2/1),({7,4,22,10},3/1),({7,4,21,11},3/1),({7,4,20,12},4/1),({7,4,19,13},4/1),({7,4,18,14},3/1),({7,4,17,15},2/1),({7,4,16,16},1/1),({6,5,24,8},1/1),({6,5,23,9},1/1),({6,5,22,10},1/1),({6,5,21,11},2/1),({6,5,20,12},2/1),({6,5,19,13},2/1),({6,5,18,14},2/1),({6,5,17,15},1/1)}, (4,2) => {}, (6,0) => {}, (6,1) => {({14,1,26,18},1/1),({14,1,25,19},1/1),({14,1,24,20},1/1),({13,2,30,14},1/1),({13,2,29,15},2/1),({13,2,28,16},3/1),({13,2,27,17},4/1),({13,2,26,18},6/1),({13,2,25,19},7/1),({13,2,24,20},6/1),({13,2,23,21},4/1),({13,2,22,22},1/1),({12,3,32,12},1/1),({12,3,31,13},2/1),({12,3,30,14},6/1),({12,3,29,15},10/1),({12,3,28,16},16/1),({12,3,27,17},19/1),({12,3,26,18},24/1),({12,3,25,19},23/1),({12,3,24,20},22/1),({12,3,23,21},13/1),({12,3,22,22},6/1),({11,4,34,10},1/1),({11,4,33,11},3/1),({11,4,32,12},7/1),({11,4,31,13},13/1),({11,4,30,14},22/1),({11,4,29,15},33/1),({11,4,28,16},45/1),({11,4,27,17},55/1),({11,4,26,18},60/1),({11,4,25,19},59/1),({11,4,24,20},50/1),({11,4,23,21},34/1),({11,4,22,22},12/1),({10,5,35,9},1/1),({10,5,34,10},3/1),({10,5,33,11},8/1),({10,5,32,12},17/1),({10,5,31,13},29/1),({10,5,30,14},46/1),({10,5,29,15},64/1),({10,5,28,16},84/1),({10,5,27,17},97/1),({10,5,26,18},105/1),({10,5,25,19},98/1),({10,5,24,20},84/1),({10,5,23,21},54/1),({10,5,22,22},21/1),({9,6,36,8},1/1),({9,6,35,9},3/1),({9,6,34,10},7/1),({9,6,33,11},14/1),({9,6,32,12},26/1),({9,6,31,13},42/1),({9,6,30,14},62/1),({9,6,29,15},84/1),({9,6,28,16},104/1),({9,6,27,17},119/1),({9,6,26,18},125/1),({9,6,25,19},118/1),({9,6,24,20},97/1),({9,6,23,21},64/1),({9,6,22,22},22/1),({8,7,36,8},1/1),({8,7,35,9},2/1),({8,7,34,10},6/1),({8,7,33,11},11/1),({8,7,32,12},20/1),({8,7,31,13},31/1),({8,7,30,14},46/1),({8,7,29,15},60/1),({8,7,28,16},75/1),({8,7,27,17},83/1),({8,7,26,18},88/1),({8,7,25,19},81/1),({8,7,24,20},68/1),({8,7,23,21},42/1),({8,7,22,22},16/1)}, (6,2) => {}, (8,0) => {}, (8,1) => {({16,3,35,21},1/1),({16,3,34,22},1/1),({16,3,33,23},2/1),({16,3,32,24},2/1),({16,3,31,25},3/1),({16,3,30,26},2/1),({16,3,29,27},2/1),({15,4,37,19},2/1),({15,4,36,20},4/1),({15,4,35,21},7/1),({15,4,34,22},11/1),({15,4,33,23},15/1),({15,4,32,24},18/1),({15,4,31,25},18/1),({15,4,30,26},16/1),({15,4,29,27},11/1),({15,4,28,28},4/1),({14,5,40,16},1/1),({14,5,39,17},3/1),({14,5,38,18},7/1),({14,5,37,19},15/1),({14,5,36,20},25/1),({14,5,35,21},39/1),({14,5,34,22},52/1),({14,5,33,23},66/1),({14,5,32,24},71/1),({14,5,31,25},72/1),({14,5,30,26},59/1),({14,5,29,27},42/1),({14,5,28,28},13/1),({13,6,41,15},2/1),({13,6,40,16},6/1),({13,6,39,17},15/1),({13,6,38,18},30/1),({13,6,37,19},52/1),({13,6,36,20},81/1),({13,6,35,21},115/1),({13,6,34,22},147/1),({13,6,33,23},173/1),({13,6,32,24},185/1),({13,6,31,25},177/1),({13,6,30,26},148/1),({13,6,29,27},98/1),({13,6,28,28},35/1),({12,7,43,13},1/1),({12,7,42,14},3/1),({12,7,41,15},9/1),({12,7,40,16},20/1),({12,7,39,17},41/1),({12,7,38,18},71/1),({12,7,37,19},115/1),({12,7,36,20},166/1),({12,7,35,21},225/1),({12,7,34,22},276/1),({12,7,33,23},317/1),({12,7,32,24},327/1),({12,7,31,25},312/1),({12,7,30,26},253/1),({12,7,29,27},170/1),({12,7,28,28},56/1),({11,8,43,13},2/1),({11,8,42,14},6/1),({11,8,41,15},15/1),({11,8,40,16},32/1),({11,8,39,17},60/1),({11,8,38,18},100/1),({11,8,37,19},154/1),({11,8,36,20},217/1),({11,8,35,21},283/1),({11,8,34,22},343/1),({11,8,33,23},384/1),({11,8,32,24},395/1),({11,8,31,25},368/1),({11,8,30,26},301/1),({11,8,29,27},197/1),({11,8,28,28},69/1),({10,9,44,12},1/1),({10,9,43,13},2/1),({10,9,42,14},6/1),({10,9,41,15},14/1),({10,9,40,16},27/1),({10,9,39,17},48/1),({10,9,38,18},78/1),({10,9,37,19},117/1),({10,9,36,20},161/1),({10,9,35,21},208/1),({10,9,34,22},247/1),({10,9,33,23},275/1),({10,9,32,24},280/1),({10,9,31,25},260/1),({10,9,30,26},210/1),({10,9,29,27},139/1),({10,9,28,28},47/1)}, (8,2) => {}, (10,0) => {}, (10,1) => {({18,5,40,28},1/1),({18,5,39,29},1/1),({18,5,38,30},2/1),({18,5,37,31},1/1),({18,5,36,32},2/1),({18,5,35,33},1/1),({18,5,34,34},1/1),({17,6,44,24},1/1),({17,6,43,25},2/1),({17,6,42,26},4/1),({17,6,41,27},7/1),({17,6,40,28},11/1),({17,6,39,29},15/1),({17,6,38,30},17/1),({17,6,37,31},17/1),({17,6,36,32},15/1),({17,6,35,33},11/1),({17,6,34,34},4/1),({16,7,46,22},1/1),({16,7,45,23},3/1),({16,7,44,24},8/1),({16,7,43,25},16/1),({16,7,42,26},28/1),({16,7,41,27},41/1),({16,7,40,28},57/1),({16,7,39,29},69/1),({16,7,38,30},78/1),({16,7,37,31},74/1),({16,7,36,32},65/1),({16,7,35,33},42/1),({16,7,34,34},17/1),({15,8,48,20},1/1),({15,8,47,21},3/1),({15,8,46,22},9/1),({15,8,45,23},19/1),({15,8,44,24},37/1),({15,8,43,25},63/1),({15,8,42,26},97/1),({15,8,41,27},135/1),({15,8,40,28},172/1),({15,8,39,29},201/1),({15,8,38,30},214/1),({15,8,37,31},205/1),({15,8,36,32},170/1),({15,8,35,33},113/1),({15,8,34,34},39/1),({14,9,49,19},1/1),({14,9,48,20},4/1),({14,9,47,21},11/1),({14,9,46,22},26/1),({14,9,45,23},50/1),({14,9,44,24},89/1),({14,9,43,25},139/1),({14,9,42,26},204/1),({14,9,41,27},270/1),({14,9,40,28},335/1),({14,9,39,29},377/1),({14,9,38,30},396/1),({14,9,37,31},369/1),({14,9,36,32},307/1),({14,9,35,33},198/1),({14,9,34,34},72/1),({13,10,50,18},1/1),({13,10,49,19},3/1),({13,10,48,20},9/1),({13,10,47,21},21/1),({13,10,46,22},43/1),({13,10,45,23},78/1),({13,10,44,24},129/1),({13,10,43,25},195/1),({13,10,42,26},273/1),({13,10,41,27},354/1),({13,10,40,28},427/1),({13,10,39,29},476/1),({13,10,38,30},488/1),({13,10,37,31},454/1),({13,10,36,32},370/1),({13,10,35,33},243/1),({13,10,34,34},84/1),({12,11,50,18},1/1),({12,11,49,19},3/1),({12,11,48,20},8/1),({12,11,47,21},18/1),({12,11,46,22},35/1),({12,11,45,23},62/1),({12,11,44,24},100/1),({12,11,43,25},148/1),({12,11,42,26},204/1),({12,11,41,27},260/1),({12,11,40,28},311/1),({12,11,39,29},343/1),({12,11,38,30},350/1),({12,11,37,31},322/1),({12,11,36,32},263/1),({12,11,35,33},171/1),({12,11,34,34},61/1)}, (10,2) => {}, (12,0) => {}, (12,1) => {({20,7,41,39},1/1),({19,8,47,33},1/1),({19,8,46,34},2/1),({19,8,45,35},3/1),({19,8,44,36},4/1),({19,8,43,37},4/1),({19,8,42,38},4/1),({19,8,41,39},3/1),({19,8,40,40},1/1),({18,9,51,29},1/1),({18,9,50,30},2/1),({18,9,49,31},5/1),({18,9,48,32},8/1),({18,9,47,33},14/1),({18,9,46,34},19/1),({18,9,45,35},26/1),({18,9,44,36},28/1),({18,9,43,37},30/1),({18,9,42,38},24/1),({18,9,41,39},18/1),({18,9,40,40},5/1),({17,10,53,27},1/1),({17,10,52,28},2/1),({17,10,51,29},6/1),({17,10,50,30},13/1),({17,10,49,31},24/1),({17,10,48,32},39/1),({17,10,47,33},57/1),({17,10,46,34},75/1),({17,10,45,35},90/1),({17,10,44,36},98/1),({17,10,43,37},95/1),({17,10,42,38},80/1),({17,10,41,39},53/1),({17,10,40,40},19/1),({16,11,54,26},1/1),({16,11,53,27},4/1),({16,11,52,28},10/1),({16,11,51,29},21/1),({16,11,50,30},39/1),({16,11,49,31},65/1),({16,11,48,32},97/1),({16,11,47,33},134/1),({16,11,46,34},167/1),({16,11,45,35},195/1),({16,11,44,36},204/1),({16,11,43,37},196/1),({16,11,42,38},160/1),({16,11,41,39},108/1),({16,11,40,40},36/1),({15,12,55,25},1/1),({15,12,54,26},3/1),({15,12,53,27},8/1),({15,12,52,28},18/1),({15,12,51,29},36/1),({15,12,50,30},62/1),({15,12,49,31},98/1),({15,12,48,32},141/1),({15,12,47,33},188/1),({15,12,46,34},231/1),({15,12,45,35},261/1),({15,12,44,36},271/1),({15,12,43,37},254/1),({15,12,42,38},209/1),({15,12,41,39},137/1),({15,12,40,40},48/1),({14,13,55,25},1/1),({14,13,54,26},3/1),({14,13,53,27},8/1),({14,13,52,28},16/1),({14,13,51,29},31/1),({14,13,50,30},51/1),({14,13,49,31},79/1),({14,13,48,32},110/1),({14,13,47,33},146/1),({14,13,46,34},175/1),({14,13,45,35},197/1),({14,13,44,36},201/1),({14,13,43,37},189/1),({14,13,42,38},153/1),({14,13,41,39},102/1),({14,13,40,40},34/1)}, (12,2) => {}, (14,0) => {}, (14,1) => {({20,11,52,40},1/1),({20,11,51,41},1/1),({20,11,50,42},2/1),({20,11,49,43},2/1),({20,11,48,44},3/1),({20,11,47,45},1/1),({20,11,46,46},1/1),({19,12,56,36},1/1),({19,12,55,37},2/1),({19,12,54,38},4/1),({19,12,53,39},7/1),({19,12,52,40},10/1),({19,12,51,41},13/1),({19,12,50,42},15/1),({19,12,49,43},15/1),({19,12,48,44},13/1),({19,12,47,45},9/1),({19,12,46,46},3/1),({18,13,58,34},1/1),({18,13,57,35},2/1),({18,13,56,36},5/1),({18,13,55,37},9/1),({18,13,54,38},16/1),({18,13,53,39},23/1),({18,13,52,40},32/1),({18,13,51,41},38/1),({18,13,50,42},43/1),({18,13,49,43},41/1),({18,13,48,44},36/1),({18,13,47,45},23/1),({18,13,46,46},9/1),({17,14,58,34},2/1),({17,14,57,35},5/1),({17,14,56,36},10/1),({17,14,55,37},18/1),({17,14,54,38},29/1),({17,14,53,39},41/1),({17,14,52,40},53/1),({17,14,51,41},62/1),({17,14,50,42},66/1),({17,14,49,43},64/1),({17,14,48,44},53/1),({17,14,47,45},35/1),({17,14,46,46},12/1),({16,15,59,33},1/1),({16,15,58,34},2/1),({16,15,57,35},5/1),({16,15,56,36},10/1),({16,15,55,37},16/1),({16,15,54,38},25/1),({16,15,53,39},35/1),({16,15,52,40},44/1),({16,15,51,41},50/1),({16,15,50,42},54/1),({16,15,49,43},51/1),({16,15,48,44},43/1),({16,15,47,45},27/1),({16,15,46,46},10/1)}, (16,0) => {}, (14,2) => {}, (16,1) => {({20,15,59,45},1/1),({20,15,58,46},1/1),({20,15,57,47},2/1),({20,15,56,48},2/1),({20,15,55,49},3/1),({20,15,54,50},2/1),({20,15,53,51},2/1),({19,16,61,43},1/1),({19,16,60,44},2/1),({19,16,59,45},3/1),({19,16,58,46},4/1),({19,16,57,47},5/1),({19,16,56,48},6/1),({19,16,55,49},6/1),({19,16,54,50},5/1),({19,16,53,51},3/1),({19,16,52,52},1/1),({18,17,61,43},1/1),({18,17,60,44},2/1),({18,17,59,45},3/1),({18,17,58,46},4/1),({18,17,57,47},5/1),({18,17,56,48},6/1),({18,17,55,49},6/1),({18,17,54,50},5/1),({18,17,53,51},3/1),({18,17,52,52},1/1)}, (18,0) => {}, (16,2) => {}, (18,1) => {({20,19,62,54},1/1)}, (18,2) => {}, (20,1) => {}, (1,0) => {({3,0,7,1},1/1),({3,0,6,2},1/1),({2,1,8,0},1/1),({2,1,7,1},1/1),({2,1,6,2},1/1)}, (1,1) => {}, (3,0) => {({6,1,17,3},1/1),({6,1,16,4},1/1),({6,1,15,5},2/1),({6,1,14,6},2/1),({6,1,13,7},3/1),({6,1,12,8},2/1),({6,1,11,9},2/1),({5,2,18,2},1/1),({5,2,17,3},2/1),({5,2,16,4},3/1),({5,2,15,5},4/1),({5,2,14,6},5/1),({5,2,13,7},6/1),({5,2,12,8},5/1),({5,2,11,9},3/1),({5,2,10,10},1/1),({4,3,18,2},1/1),({4,3,17,3},2/1),({4,3,16,4},3/1),({4,3,15,5},4/1),({4,3,14,6},5/1),({4,3,13,7},6/1),({4,3,12,8},5/1),({4,3,11,9},3/1),({4,3,10,10},1/1)}, (3,1) => {({9,0,17,9},1/1),({9,0,15,11},1/1),({9,0,13,13},1/1)}, (3,2) => {}, (5,0) => {({8,3,26,6},1/1),({8,3,24,8},1/1),({8,3,23,9},1/1),({8,3,22,10},2/1),({8,3,21,11},1/1),({8,3,20,12},2/1),({8,3,19,13},1/1),({8,3,18,14},2/1),({8,3,16,16},1/1),({7,4,25,7},1/1),({7,4,24,8},1/1),({7,4,23,9},2/1),({7,4,22,10},3/1),({7,4,21,11},3/1),({7,4,20,12},4/1),({7,4,19,13},4/1),({7,4,18,14},3/1),({7,4,17,15},2/1),({7,4,16,16},1/1),({6,5,24,8},1/1),({6,5,23,9},1/1),({6,5,22,10},2/1),({6,5,21,11},3/1),({6,5,20,12},4/1),({6,5,19,13},3/1),({6,5,18,14},4/1),({6,5,17,15},2/1),({6,5,16,16},1/1)}, (5,1) => {({13,0,20,18},1/1),({12,1,25,13},1/1),({12,1,24,14},2/1),({12,1,23,15},2/1),({12,1,22,16},2/1),({12,1,21,17},2/1),({12,1,20,18},2/1),({12,1,19,19},1/1),({11,2,28,10},1/1),({11,2,27,11},1/1),({11,2,26,12},3/1),({11,2,25,13},4/1),({11,2,24,14},8/1),({11,2,23,15},8/1),({11,2,22,16},10/1),({11,2,21,17},7/1),({11,2,20,18},7/1),({11,2,19,19},1/1),({10,3,29,9},1/1),({10,3,28,10},3/1),({10,3,27,11},6/1),({10,3,26,12},10/1),({10,3,25,13},14/1),({10,3,24,14},18/1),({10,3,23,15},21/1),({10,3,22,16},21/1),({10,3,21,17},18/1),({10,3,20,18},12/1),({10,3,19,19},4/1),({9,4,31,7},1/1),({9,4,30,8},2/1),({9,4,29,9},4/1),({9,4,28,10},8/1),({9,4,27,11},13/1),({9,4,26,12},20/1),({9,4,25,13},25/1),({9,4,24,14},32/1),({9,4,23,15},33/1),({9,4,22,16},35/1),({9,4,21,17},27/1),({9,4,20,18},20/1),({9,4,19,19},5/1),({8,5,31,7},1/1),({8,5,30,8},3/1),({8,5,29,9},6/1),({8,5,28,10},11/1),({8,5,27,11},17/1),({8,5,26,12},24/1),({8,5,25,13},31/1),({8,5,24,14},36/1),({8,5,23,15},38/1),({8,5,22,16},36/1),({8,5,21,17},30/1),({8,5,20,18},20/1),({8,5,19,19},7/1),({7,6,32,6},1/1),({7,6,31,7},1/1),({7,6,30,8},3/1),({7,6,29,9},5/1),({7,6,28,10},9/1),({7,6,27,11},12/1),({7,6,26,12},18/1),({7,6,25,13},21/1),({7,6,24,14},26/1),({7,6,23,15},25/1),({7,6,22,16},25/1),({7,6,21,17},19/1),({7,6,20,18},15/1),({7,6,19,19},4/1)}};
--dw stands for dominant weights
dw1226 = new HashTable from {(5,2) => {}, (7,0) => {}, (7,1) => {({15,2,31,19},1/1),({14,3,34,16},1/1),({13,4,36,14},1/1),({12,5,38,12},1/1),({11,6,39,11},1/1),({10,7,40,10},1/1)}, (9,0) => {}, (7,2) => {}, (9,1) => {({17,4,38,24},1/1),({16,5,41,21},1/1),({15,6,43,19},1/1),({14,7,45,17},1/1),({13,8,46,16},1/1),({12,9,47,15},1/1)}, (9,2) => {}, (11,0) => {}, (11,1) => {({19,6,41,33},1/1),({18,7,46,28},1/1),({17,8,49,25},1/1),({16,9,50,24},2/1),({15,10,52,22},1/1),({13,12,53,21},1/1)}, (11,2) => {}, (13,0) => {}, (13,1) => {({20,9,47,39},1/1),({19,10,52,34},1/1),({18,11,55,31},1/1),({17,12,56,30},1/1),({16,13,57,29},1/1)}, (13,2) => {}, (15,0) => {}, (15,1) => {({20,13,56,42},1/1),({19,14,59,39},1/1),({18,15,60,38},1/1)}, (15,2) => {}, (17,0) => {}, (17,1) => {({20,17,61,49},1/1),({19,18,62,48},1/1)}, (17,2) => {}, (19,1) => {}, (0,0) => {({1,0,2,0},1/1)}, (2,0) => {({5,0,11,3},1/1),({4,1,13,1},1/1)}, (2,1) => {}, (2,2) => {}, (4,0) => {({7,2,22,4},1/1),({5,4,23,3},1/1)}, (4,1) => {({11,0,19,13},1/1),({10,1,23,9},1/1),({8,3,26,6},1/1)}, (4,2) => {}, (6,0) => {}, (6,1) => {({14,1,26,18},1/1),({13,2,30,14},1/1),({12,3,32,12},1/1),({11,4,34,10},1/1),({10,5,35,9},1/1),({9,6,36,8},1/1)}, (6,2) => {}, (8,0) => {}, (8,1) => {({16,3,35,21},1/1),({15,4,37,19},2/1),({14,5,40,16},1/1),({13,6,41,15},2/1),({12,7,43,13},1/1),({10,9,44,12},1/1)}, (8,2) => {}, (10,0) => {}, (10,1) => {({18,5,40,28},1/1),({17,6,44,24},1/1),({16,7,46,22},1/1),({15,8,48,20},1/1),({14,9,49,19},1/1),({13,10,50,18},1/1)}, (10,2) => {}, (12,0) => {}, (12,1) => {({20,7,41,39},1/1),({19,8,47,33},1/1),({18,9,51,29},1/1),({17,10,53,27},1/1),({16,11,54,26},1/1),({15,12,55,25},1/1)}, (12,2) => {}, (14,0) => {}, (14,1) => {({20,11,52,40},1/1),({19,12,56,36},1/1),({18,13,58,34},1/1),({16,15,59,33},1/1)}, (16,0) => {}, (14,2) => {}, (16,1) => {({20,15,59,45},1/1),({19,16,61,43},1/1)}, (18,0) => {}, (16,2) => {}, (18,1) => {({20,19,62,54},1/1)}, (18,2) => {}, (20,1) => {}, (1,0) => {({3,0,7,1},1/1),({2,1,8,0},1/1)}, (1,1) => {}, (3,0) => {({6,1,17,3},1/1),({5,2,18,2},1/1)}, (3,1) => {({9,0,17,9},1/1)}, (5,0) => {({8,3,26,6},1/1)}, (3,2) => {}, (5,1) => {({13,0,20,18},1/1),({12,1,25,13},1/1),({11,2,28,10},1/1),({10,3,29,9},1/1),({9,4,31,7},1/1),({7,6,32,6},1/1)}};
--dmw stands for dominant monomial weights
dmw1226 = new HashTable from {(5,2) => {}, (7,0) => {}, (7,1) => {{13,4,36,14},{12,5,38,12},{14,3,34,16},{15,2,31,19},{11,6,39,11}}, (9,0) => {}, (7,2) => {}, (9,1) => {{17,4,38,24},{13,8,46,16},{12,9,47,15},{15,6,43,19},{14,7,45,17},{16,5,41,21}}, (9,2) => {}, (11,0) => {}, (11,1) => {{17,8,49,25},{19,6,41,33},{18,7,46,28},{16,9,50,24},{15,10,52,22},{13,12,53,21}}, (11,2) => {}, (13,0) => {}, (13,1) => {{20,9,47,39},{16,13,57,29},{18,11,55,31},{17,12,56,30},{19,10,52,34}}, (13,2) => {}, (15,0) => {}, (15,1) => {{18,15,60,38},{19,14,59,39},{20,13,56,42}}, (15,2) => {}, (17,0) => {}, (17,1) => {{19,18,62,48},{20,17,61,49}}, (17,2) => {{20,19,56,60}}, (19,1) => {}, (0,0) => {}, (2,0) => {}, (2,1) => {}, (2,2) => {}, (4,0) => {}, (4,1) => {{11,0,19,13},{10,1,23,9}}, (4,2) => {}, (6,0) => {}, (6,1) => {{12,3,32,12},{14,1,26,18},{13,2,30,14},{11,4,34,10}}, (6,2) => {}, (8,0) => {}, (8,1) => {{16,3,35,21},{12,7,43,13},{14,5,40,16},{13,6,41,15},{15,4,37,19}}, (8,2) => {}, (10,0) => {}, (10,1) => {{16,7,46,22},{17,6,44,24},{15,8,48,20},{18,5,40,28},{14,9,49,19},{13,10,50,18}}, (10,2) => {}, (12,0) => {}, (12,1) => {{19,8,47,33},{20,7,41,39},{16,11,54,26},{15,12,55,25},{18,9,51,29},{17,10,53,27}}, (12,2) => {}, (14,0) => {}, (14,1) => {{20,11,52,40},{19,12,56,36},{18,13,58,34},{16,15,59,33}}, (16,0) => {}, (14,2) => {}, (16,1) => {{20,15,59,45},{19,16,61,43}}, (18,0) => {}, (16,2) => {{19,18,50,60}}, (18,1) => {{20,19,62,54}}, (18,2) => {{20,21,62,60}}, (20,1) => {}, (1,0) => {}, (1,1) => {}, (3,0) => {}, (3,1) => {{9,0,17,9}}, (5,0) => {}, (3,2) => {}, (5,1) => {{13,0,20,18},{11,2,28,10},{10,3,29,9},{12,1,25,13}}};
end;