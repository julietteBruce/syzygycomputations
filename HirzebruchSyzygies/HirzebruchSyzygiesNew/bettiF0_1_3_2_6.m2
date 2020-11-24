A := QQ[t_0,t_1,t_2,t_3];
--tb stands for Total Betti numbers
tb1326 = new HashTable from {(5,2) => 0, (7,0) => 3432/1, (7,1) => 159120/1, (7,2) => 0, (9,0) => 0, (9,1) => 427856/1, (9,2) => 0, (11,0) => 0, (11,1) => 360672/1, (11,2) => 0, (13,0) => 0, (13,1) => 112608/1, (15,0) => 0, (13,2) => 0, (15,1) => 11832/1, (17,0) => 0, (15,2) => 0, (17,1) => 280/1, (17,2) => 0, (19,1) => 0, (0,0) => 8/1, (2,0) => 936/1, (2,1) => 0, (2,2) => 0, (4,0) => 11424/1, (4,1) => 72/1, (4,2) => 0, (6,0) => 16016/1, (6,1) => 45864/1, (6,2) => 0, (8,0) => 0, (8,1) => 311168/1, (8,2) => 0, (10,0) => 0, (10,1) => 445536/1, (10,2) => 0, (12,0) => 0, (12,1) => 228480/1, (12,2) => 0, (14,0) => 0, (14,1) => 42432/1, (16,0) => 0, (14,2) => 0, (16,1) => 2304/1, (18,0) => 0, (16,2) => 0, (18,1) => 16/1, (18,2) => 0, (20,1) => 0, (1,0) => 128/1, (1,1) => 0, (3,0) => 4080/1, (3,1) => 0, (3,2) => 0, (5,0) => 19656/1, (5,1) => 4592/1};
--sb represents the betti numbers as sums of Schur functors
sb1326 = new HashTable from {(5,2) => {}, (7,0) => {({11,4,33,12},1/1),({11,4,31,14},1/1),({11,4,30,15},1/1),({11,4,29,16},1/1),({11,4,28,17},1/1),({11,4,27,18},2/1),({11,4,26,19},1/1),({11,4,25,20},1/1),({11,4,24,21},1/1),({10,5,32,13},1/1),({10,5,31,14},1/1),({10,5,30,15},2/1),({10,5,29,16},2/1),({10,5,28,17},3/1),({10,5,27,18},3/1),({10,5,26,19},3/1),({10,5,25,20},3/1),({10,5,24,21},2/1),({10,5,23,22},1/1),({9,6,31,14},1/1),({9,6,30,15},1/1),({9,6,29,16},3/1),({9,6,28,17},3/1),({9,6,27,18},4/1),({9,6,26,19},4/1),({9,6,25,20},4/1),({9,6,24,21},3/1),({9,6,23,22},2/1),({8,7,30,15},1/1),({8,7,29,16},1/1),({8,7,28,17},2/1),({8,7,27,18},3/1),({8,7,26,19},3/1),({8,7,25,20},3/1),({8,7,24,21},3/1),({8,7,23,22},1/1)}, (7,1) => {({15,2,31,20},1/1),({15,2,30,21},1/1),({15,2,29,22},1/1),({15,2,28,23},1/1),({15,2,27,24},1/1),({15,2,26,25},1/1),({14,3,34,17},1/1),({14,3,33,18},1/1),({14,3,32,19},3/1),({14,3,31,20},5/1),({14,3,30,21},7/1),({14,3,29,22},7/1),({14,3,28,23},7/1),({14,3,27,24},6/1),({14,3,26,25},4/1),({13,4,36,15},1/1),({13,4,35,16},3/1),({13,4,34,17},6/1),({13,4,33,18},10/1),({13,4,32,19},15/1),({13,4,31,20},21/1),({13,4,30,21},25/1),({13,4,29,22},27/1),({13,4,28,23},25/1),({13,4,27,24},20/1),({13,4,26,25},11/1),({12,5,38,13},1/1),({12,5,37,14},2/1),({12,5,36,15},6/1),({12,5,35,16},11/1),({12,5,34,17},21/1),({12,5,33,18},31/1),({12,5,32,19},44/1),({12,5,31,20},54/1),({12,5,30,21},63/1),({12,5,29,22},64/1),({12,5,28,23},59/1),({12,5,27,24},44/1),({12,5,26,25},24/1),({11,6,39,12},1/1),({11,6,38,13},3/1),({11,6,37,14},8/1),({11,6,36,15},15/1),({11,6,35,16},27/1),({11,6,34,17},43/1),({11,6,33,18},62/1),({11,6,32,19},81/1),({11,6,31,20},98/1),({11,6,30,21},108/1),({11,6,29,22},109/1),({11,6,28,23},97/1),({11,6,27,24},73/1),({11,6,26,25},39/1),({10,7,40,11},1/1),({10,7,39,12},2/1),({10,7,38,13},6/1),({10,7,37,14},12/1),({10,7,36,15},23/1),({10,7,35,16},37/1),({10,7,34,17},57/1),({10,7,33,18},78/1),({10,7,32,19},101/1),({10,7,31,20},118/1),({10,7,30,21},129/1),({10,7,29,22},126/1),({10,7,28,23},112/1),({10,7,27,24},83/1),({10,7,26,25},45/1),({9,8,40,11},1/1),({9,8,39,12},2/1),({9,8,38,13},5/1),({9,8,37,14},10/1),({9,8,36,15},18/1),({9,8,35,16},29/1),({9,8,34,17},42/1),({9,8,33,18},57/1),({9,8,32,19},72/1),({9,8,31,20},84/1),({9,8,30,21},90/1),({9,8,29,22},88/1),({9,8,28,23},77/1),({9,8,27,24},58/1),({9,8,26,25},31/1)}, (7,2) => {}, (9,0) => {}, (9,1) => {({17,4,38,25},1/1),({17,4,37,26},1/1),({17,4,36,27},1/1),({17,4,35,28},2/1),({17,4,34,29},2/1),({17,4,33,30},1/1),({17,4,32,31},1/1),({16,5,41,22},1/1),({16,5,40,23},2/1),({16,5,39,24},5/1),({16,5,38,25},7/1),({16,5,37,26},11/1),({16,5,36,27},13/1),({16,5,35,28},15/1),({16,5,34,29},14/1),({16,5,33,30},11/1),({16,5,32,31},6/1),({15,6,43,20},1/1),({15,6,42,21},4/1),({15,6,41,22},8/1),({15,6,40,23},16/1),({15,6,39,24},26/1),({15,6,38,25},38/1),({15,6,37,26},49/1),({15,6,36,27},58/1),({15,6,35,28},60/1),({15,6,34,29},56/1),({15,6,33,30},43/1),({15,6,32,31},23/1),({14,7,45,18},1/1),({14,7,44,19},3/1),({14,7,43,20},9/1),({14,7,42,21},18/1),({14,7,41,22},35/1),({14,7,40,23},56/1),({14,7,39,24},84/1),({14,7,38,25},113/1),({14,7,37,26},140/1),({14,7,36,27},155/1),({14,7,35,28},160/1),({14,7,34,29},142/1),({14,7,33,30},108/1),({14,7,32,31},59/1),({13,8,46,17},1/1),({13,8,45,18},4/1),({13,8,44,19},11/1),({13,8,43,20},24/1),({13,8,42,21},46/1),({13,8,41,22},78/1),({13,8,40,23},120/1),({13,8,39,24},168/1),({13,8,38,25},218/1),({13,8,37,26},259/1),({13,8,36,27},284/1),({13,8,35,28},282/1),({13,8,34,29},251/1),({13,8,33,30},188/1),({13,8,32,31},101/1),({12,9,47,16},1/1),({12,9,46,17},3/1),({12,9,45,18},9/1),({12,9,44,19},19/1),({12,9,43,20},39/1),({12,9,42,21},68/1),({12,9,41,22},110/1),({12,9,40,23},161/1),({12,9,39,24},221/1),({12,9,38,25},276/1),({12,9,37,26},324/1),({12,9,36,27},348/1),({12,9,35,28},342/1),({12,9,34,29},300/1),({12,9,33,30},225/1),({12,9,32,31},119/1),({11,10,47,16},1/1),({11,10,46,17},3/1),({11,10,45,18},7/1),({11,10,44,19},17/1),({11,10,43,20},31/1),({11,10,42,21},53/1),({11,10,41,22},84/1),({11,10,40,23},121/1),({11,10,39,24},161/1),({11,10,38,25},202/1),({11,10,37,26},232/1),({11,10,36,27},247/1),({11,10,35,28},243/1),({11,10,34,29},212/1),({11,10,33,30},156/1),({11,10,32,31},85/1)}, (9,2) => {}, (11,0) => {}, (11,1) => {({19,6,41,34},1/1),({19,6,40,35},1/1),({18,7,46,29},1/1),({18,7,45,30},2/1),({18,7,44,31},4/1),({18,7,43,32},5/1),({18,7,42,33},7/1),({18,7,41,34},8/1),({18,7,40,35},9/1),({18,7,39,36},6/1),({18,7,38,37},3/1),({17,8,49,26},1/1),({17,8,48,27},2/1),({17,8,47,28},5/1),({17,8,46,29},10/1),({17,8,45,30},18/1),({17,8,44,31},26/1),({17,8,43,32},35/1),({17,8,42,33},41/1),({17,8,41,34},44/1),({17,8,40,35},41/1),({17,8,39,36},32/1),({17,8,38,37},17/1),({16,9,50,25},2/1),({16,9,49,26},5/1),({16,9,48,27},13/1),({16,9,47,28},25/1),({16,9,46,29},43/1),({16,9,45,30},64/1),({16,9,44,31},91/1),({16,9,43,32},112/1),({16,9,42,33},128/1),({16,9,41,34},131/1),({16,9,40,35},119/1),({16,9,39,36},90/1),({16,9,38,37},50/1),({15,10,52,23},1/1),({15,10,51,24},3/1),({15,10,50,25},8/1),({15,10,49,26},19/1),({15,10,48,27},37/1),({15,10,47,28},65/1),({15,10,46,29},101/1),({15,10,45,30},145/1),({15,10,44,31},190/1),({15,10,43,32},230/1),({15,10,42,33},253/1),({15,10,41,34},254/1),({15,10,40,35},226/1),({15,10,39,36},171/1),({15,10,38,37},92/1),({14,11,52,23},2/1),({14,11,51,24},6/1),({14,11,50,25},15/1),({14,11,49,26},31/1),({14,11,48,27},58/1),({14,11,47,28},94/1),({14,11,46,29},144/1),({14,11,45,30},198/1),({14,11,44,31},253/1),({14,11,43,32},298/1),({14,11,42,33},325/1),({14,11,41,34},320/1),({14,11,40,35},284/1),({14,11,39,36},212/1),({14,11,38,37},113/1),({13,12,53,22},1/1),({13,12,52,23},2/1),({13,12,51,24},6/1),({13,12,50,25},14/1),({13,12,49,26},27/1),({13,12,48,27},47/1),({13,12,47,28},76/1),({13,12,46,29},111/1),({13,12,45,30},151/1),({13,12,44,31},190/1),({13,12,43,32},221/1),({13,12,42,33},237/1),({13,12,41,34},235/1),({13,12,40,35},205/1),({13,12,39,36},152/1),({13,12,38,37},82/1)}, (11,2) => {}, (13,0) => {}, (13,1) => {({20,9,47,40},1/1),({20,9,46,41},1/1),({20,9,45,42},1/1),({19,10,52,35},1/1),({19,10,51,36},2/1),({19,10,50,37},4/1),({19,10,49,38},6/1),({19,10,48,39},8/1),({19,10,47,40},9/1),({19,10,46,41},9/1),({19,10,45,42},7/1),({19,10,44,43},4/1),({18,11,55,32},1/1),({18,11,54,33},2/1),({18,11,53,34},5/1),({18,11,52,35},9/1),({18,11,51,36},16/1),({18,11,50,37},23/1),({18,11,49,38},32/1),({18,11,48,39},37/1),({18,11,47,40},40/1),({18,11,46,41},37/1),({18,11,45,42},29/1),({18,11,44,43},16/1),({17,12,56,31},1/1),({17,12,55,32},3/1),({17,12,54,33},8/1),({17,12,53,34},16/1),({17,12,52,35},29/1),({17,12,51,36},44/1),({17,12,50,37},62/1),({17,12,49,38},78/1),({17,12,48,39},90/1),({17,12,47,40},92/1),({17,12,46,41},84/1),({17,12,45,42},64/1),({17,12,44,43},35/1),({16,13,57,30},1/1),({16,13,56,31},3/1),({16,13,55,32},8/1),({16,13,54,33},16/1),({16,13,53,34},30/1),({16,13,52,35},48/1),({16,13,51,36},71/1),({16,13,50,37},94/1),({16,13,49,38},116/1),({16,13,48,39},129/1),({16,13,47,40},131/1),({16,13,46,41},117/1),({16,13,45,42},89/1),({16,13,44,43},48/1),({15,14,57,30},1/1),({15,14,56,31},3/1),({15,14,55,32},7/1),({15,14,54,33},14/1),({15,14,53,34},25/1),({15,14,52,35},40/1),({15,14,51,36},57/1),({15,14,50,37},75/1),({15,14,49,38},90/1),({15,14,48,39},100/1),({15,14,47,40},100/1),({15,14,46,41},89/1),({15,14,45,42},67/1),({15,14,44,43},36/1)}, (15,0) => {}, (13,2) => {}, (15,1) => {({20,13,56,43},1/1),({20,13,55,44},1/1),({20,13,54,45},2/1),({20,13,53,46},3/1),({20,13,52,47},3/1),({20,13,51,48},2/1),({20,13,50,49},2/1),({19,14,59,40},1/1),({19,14,58,41},2/1),({19,14,57,42},4/1),({19,14,56,43},6/1),({19,14,55,44},9/1),({19,14,54,45},11/1),({19,14,53,46},12/1),({19,14,52,47},11/1),({19,14,51,48},9/1),({19,14,50,49},5/1),({18,15,60,39},1/1),({18,15,59,40},2/1),({18,15,58,41},5/1),({18,15,57,42},8/1),({18,15,56,43},12/1),({18,15,55,44},17/1),({18,15,54,45},21/1),({18,15,53,46},21/1),({18,15,52,47},20/1),({18,15,51,48},16/1),({18,15,50,49},9/1),({17,16,60,39},1/1),({17,16,59,40},3/1),({17,16,58,41},5/1),({17,16,57,42},8/1),({17,16,56,43},12/1),({17,16,55,44},16/1),({17,16,54,45},18/1),({17,16,53,46},19/1),({17,16,52,47},17/1),({17,16,51,48},13/1),({17,16,50,49},8/1)}, (17,0) => {}, (15,2) => {}, (17,1) => {({20,17,61,50},1/1),({20,17,60,51},1/1),({20,17,59,52},1/1),({20,17,58,53},1/1),({20,17,57,54},1/1),({20,17,56,55},1/1),({19,18,62,49},1/1),({19,18,61,50},1/1),({19,18,60,51},1/1),({19,18,59,52},1/1),({19,18,58,53},1/1),({19,18,57,54},1/1),({19,18,56,55},1/1)}, (17,2) => {}, (19,1) => {}, (0,0) => {({1,0,3,0},1/1)}, (2,0) => {({5,0,12,3},1/1),({5,0,11,4},1/1),({5,0,10,5},1/1),({5,0,9,6},1/1),({5,0,8,7},1/1),({4,1,14,1},1/1),({4,1,13,2},2/1),({4,1,12,3},3/1),({4,1,11,4},3/1),({4,1,10,5},3/1),({4,1,9,6},3/1),({4,1,8,7},2/1),({3,2,14,1},1/1),({3,2,13,2},2/1),({3,2,12,3},3/1),({3,2,11,4},3/1),({3,2,10,5},3/1),({3,2,9,6},3/1),({3,2,8,7},2/1)}, (2,1) => {}, (2,2) => {}, (4,0) => {({8,1,21,6},1/1),({8,1,20,7},1/1),({8,1,19,8},2/1),({8,1,18,9},3/1),({8,1,17,10},3/1),({8,1,16,11},3/1),({8,1,15,12},3/1),({8,1,14,13},1/1),({7,2,23,4},1/1),({7,2,22,5},2/1),({7,2,21,6},4/1),({7,2,20,7},6/1),({7,2,19,8},9/1),({7,2,18,9},11/1),({7,2,17,10},12/1),({7,2,16,11},11/1),({7,2,15,12},9/1),({7,2,14,13},5/1),({6,3,23,4},2/1),({6,3,22,5},4/1),({6,3,21,6},7/1),({6,3,20,7},11/1),({6,3,19,8},16/1),({6,3,18,9},19/1),({6,3,17,10},21/1),({6,3,16,11},19/1),({6,3,15,12},15/1),({6,3,14,13},9/1),({5,4,24,3},1/1),({5,4,23,4},2/1),({5,4,22,5},4/1),({5,4,21,6},7/1),({5,4,20,7},10/1),({5,4,19,8},14/1),({5,4,18,9},17/1),({5,4,17,10},17/1),({5,4,16,11},16/1),({5,4,15,12},13/1),({5,4,14,13},7/1)}, (4,1) => {({11,0,19,14},1/1)}, (4,2) => {}, (6,0) => {({10,3,30,9},1/1),({10,3,29,10},1/1),({10,3,28,11},2/1),({10,3,27,12},3/1),({10,3,26,13},3/1),({10,3,25,14},4/1),({10,3,24,15},5/1),({10,3,23,16},4/1),({10,3,22,17},4/1),({10,3,21,18},3/1),({10,3,20,19},1/1),({9,4,30,9},1/1),({9,4,29,10},2/1),({9,4,28,11},4/1),({9,4,27,12},6/1),({9,4,26,13},9/1),({9,4,25,14},11/1),({9,4,24,15},13/1),({9,4,23,16},13/1),({9,4,22,17},12/1),({9,4,21,18},9/1),({9,4,20,19},5/1),({8,5,32,7},1/1),({8,5,31,8},1/1),({8,5,30,9},2/1),({8,5,29,10},4/1),({8,5,28,11},7/1),({8,5,27,12},10/1),({8,5,26,13},15/1),({8,5,25,14},18/1),({8,5,24,15},21/1),({8,5,23,16},22/1),({8,5,22,17},20/1),({8,5,21,18},15/1),({8,5,20,19},9/1),({7,6,31,8},1/1),({7,6,30,9},2/1),({7,6,29,10},3/1),({7,6,28,11},6/1),({7,6,27,12},9/1),({7,6,26,13},12/1),({7,6,25,14},15/1),({7,6,24,15},18/1),({7,6,23,16},17/1),({7,6,22,17},16/1),({7,6,21,18},13/1),({7,6,20,19},6/1)}, (6,1) => {({14,1,26,19},1/1),({14,1,25,20},1/1),({13,2,30,15},1/1),({13,2,29,16},1/1),({13,2,28,17},2/1),({13,2,27,18},2/1),({13,2,26,19},4/1),({13,2,25,20},4/1),({13,2,24,21},3/1),({13,2,23,22},1/1),({12,3,31,14},1/1),({12,3,30,15},3/1),({12,3,29,16},5/1),({12,3,28,17},7/1),({12,3,27,18},9/1),({12,3,26,19},10/1),({12,3,25,20},10/1),({12,3,24,21},8/1),({12,3,23,22},4/1),({11,4,34,11},1/1),({11,4,33,12},2/1),({11,4,32,13},4/1),({11,4,31,14},6/1),({11,4,30,15},10/1),({11,4,29,16},14/1),({11,4,28,17},19/1),({11,4,27,18},21/1),({11,4,26,19},23/1),({11,4,25,20},20/1),({11,4,24,21},16/1),({11,4,23,22},9/1),({10,5,34,11},1/1),({10,5,33,12},3/1),({10,5,32,13},7/1),({10,5,31,14},11/1),({10,5,30,15},17/1),({10,5,29,16},23/1),({10,5,28,17},29/1),({10,5,27,18},32/1),({10,5,26,19},32/1),({10,5,25,20},29/1),({10,5,24,21},22/1),({10,5,23,22},12/1),({9,6,36,9},1/1),({9,6,35,10},1/1),({9,6,34,11},3/1),({9,6,33,12},5/1),({9,6,32,13},9/1),({9,6,31,14},14/1),({9,6,30,15},21/1),({9,6,29,16},26/1),({9,6,28,17},32/1),({9,6,27,18},35/1),({9,6,26,19},35/1),({9,6,25,20},31/1),({9,6,24,21},24/1),({9,6,23,22},12/1),({8,7,35,10},1/1),({8,7,34,11},2/1),({8,7,33,12},3/1),({8,7,32,13},6/1),({8,7,31,14},9/1),({8,7,30,15},13/1),({8,7,29,16},17/1),({8,7,28,17},20/1),({8,7,27,18},22/1),({8,7,26,19},23/1),({8,7,25,20},20/1),({8,7,24,21},14/1),({8,7,23,22},8/1)}, (6,2) => {}, (8,0) => {}, (8,1) => {({16,3,35,22},1/1),({16,3,34,23},1/1),({16,3,33,24},1/1),({16,3,32,25},2/1),({16,3,31,26},2/1),({16,3,30,27},1/1),({16,3,29,28},1/1),({15,4,37,20},2/1),({15,4,36,21},3/1),({15,4,35,22},6/1),({15,4,34,23},8/1),({15,4,33,24},11/1),({15,4,32,25},12/1),({15,4,31,26},12/1),({15,4,30,27},9/1),({15,4,29,28},5/1),({14,5,40,17},1/1),({14,5,39,18},3/1),({14,5,38,19},6/1),({14,5,37,20},12/1),({14,5,36,21},20/1),({14,5,35,22},29/1),({14,5,34,23},38/1),({14,5,33,24},45/1),({14,5,32,25},47/1),({14,5,31,26},44/1),({14,5,30,27},34/1),({14,5,29,28},18/1),({13,6,41,16},2/1),({13,6,40,17},5/1),({13,6,39,18},13/1),({13,6,38,19},24/1),({13,6,37,20},41/1),({13,6,36,21},60/1),({13,6,35,22},84/1),({13,6,34,23},102/1),({13,6,33,24},116/1),({13,6,32,25},118/1),({13,6,31,26},107/1),({13,6,30,27},80/1),({13,6,29,28},45/1),({12,7,43,14},1/1),({12,7,42,15},3/1),({12,7,41,16},8/1),({12,7,40,17},17/1),({12,7,39,18},33/1),({12,7,38,19},56/1),({12,7,37,20},87/1),({12,7,36,21},122/1),({12,7,35,22},159/1),({12,7,34,23},190/1),({12,7,33,24},208/1),({12,7,32,25},207/1),({12,7,31,26},184/1),({12,7,30,27},138/1),({12,7,29,28},74/1),({11,8,43,14},2/1),({11,8,42,15},5/1),({11,8,41,16},13/1),({11,8,40,17},26/1),({11,8,39,18},48/1),({11,8,38,19},76/1),({11,8,37,20},115/1),({11,8,36,21},156/1),({11,8,35,22},198/1),({11,8,34,23},231/1),({11,8,33,24},251/1),({11,8,32,25},245/1),({11,8,31,26},217/1),({11,8,30,27},162/1),({11,8,29,28},86/1),({10,9,44,13},1/1),({10,9,43,14},2/1),({10,9,42,15},5/1),({10,9,41,16},12/1),({10,9,40,17},22/1),({10,9,39,18},37/1),({10,9,38,19},60/1),({10,9,37,20},86/1),({10,9,36,21},115/1),({10,9,35,22},145/1),({10,9,34,23},166/1),({10,9,33,24},177/1),({10,9,32,25},175/1),({10,9,31,26},152/1),({10,9,30,27},112/1),({10,9,29,28},61/1)}, (8,2) => {}, (10,0) => {}, (10,1) => {({18,5,40,29},1/1),({18,5,39,30},1/1),({18,5,38,31},1/1),({18,5,37,32},1/1),({18,5,36,33},1/1),({18,5,35,34},1/1),({17,6,44,25},1/1),({17,6,43,26},2/1),({17,6,42,27},4/1),({17,6,41,28},6/1),({17,6,40,29},10/1),({17,6,39,30},12/1),({17,6,38,31},13/1),({17,6,37,32},12/1),({17,6,36,33},10/1),({17,6,35,34},6/1),({16,7,46,23},1/1),({16,7,45,24},3/1),({16,7,44,25},8/1),({16,7,43,26},15/1),({16,7,42,27},25/1),({16,7,41,28},36/1),({16,7,40,29},48/1),({16,7,39,30},56/1),({16,7,38,31},59/1),({16,7,37,32},54/1),({16,7,36,33},42/1),({16,7,35,34},23/1),({15,8,48,21},1/1),({15,8,47,22},3/1),({15,8,46,23},9/1),({15,8,45,24},18/1),({15,8,44,25},35/1),({15,8,43,26},57/1),({15,8,42,27},86/1),({15,8,41,28},115/1),({15,8,40,29},143/1),({15,8,39,30},160/1),({15,8,38,31},164/1),({15,8,37,32},147/1),({15,8,36,33},112/1),({15,8,35,34},60/1),({14,9,49,20},1/1),({14,9,48,21},4/1),({14,9,47,22},11/1),({14,9,46,23},25/1),({14,9,45,24},47/1),({14,9,44,25},81/1),({14,9,43,26},125/1),({14,9,42,27},177/1),({14,9,41,28},229/1),({14,9,40,29},274/1),({14,9,39,30},300/1),({14,9,38,31},300/1),({14,9,37,32},266/1),({14,9,36,33},200/1),({14,9,35,34},107/1),({13,10,50,19},1/1),({13,10,49,20},3/1),({13,10,48,21},9/1),({13,10,47,22},20/1),({13,10,46,23},41/1),({13,10,45,24},72/1),({13,10,44,25},117/1),({13,10,43,26},172/1),({13,10,42,27},236/1),({13,10,41,28},297/1),({13,10,40,29},349/1),({13,10,39,30},375/1),({13,10,38,31},370/1),({13,10,37,32},325/1),({13,10,36,33},243/1),({13,10,35,34},130/1),({12,11,50,19},1/1),({12,11,49,20},3/1),({12,11,48,21},8/1),({12,11,47,22},17/1),({12,11,46,23},33/1),({12,11,45,24},57/1),({12,11,44,25},90/1),({12,11,43,26},130/1),({12,11,42,27},175/1),({12,11,41,28},218/1),({12,11,40,29},253/1),({12,11,39,30},270/1),({12,11,38,31},264/1),({12,11,37,32},231/1),({12,11,36,33},172/1),({12,11,35,34},92/1)}, (10,2) => {}, (12,0) => {}, (12,1) => {({20,7,41,40},1/1),({19,8,47,34},1/1),({19,8,46,35},2/1),({19,8,45,36},3/1),({19,8,44,37},3/1),({19,8,43,38},3/1),({19,8,42,39},3/1),({19,8,41,40},2/1),({18,9,51,30},1/1),({18,9,50,31},2/1),({18,9,49,32},5/1),({18,9,48,33},8/1),({18,9,47,34},13/1),({18,9,46,35},18/1),({18,9,45,36},23/1),({18,9,44,37},24/1),({18,9,43,38},23/1),({18,9,42,39},18/1),({18,9,41,40},10/1),({17,10,53,28},1/1),({17,10,52,29},2/1),({17,10,51,30},6/1),({17,10,50,31},13/1),({17,10,49,32},24/1),({17,10,48,33},37/1),({17,10,47,34},54/1),({17,10,46,35},68/1),({17,10,45,36},79/1),({17,10,44,37},83/1),({17,10,43,38},76/1),({17,10,42,39},57/1),({17,10,41,40},32/1),({16,11,54,27},1/1),({16,11,53,28},4/1),({16,11,52,29},10/1),({16,11,51,30},21/1),({16,11,50,31},38/1),({16,11,49,32},63/1),({16,11,48,33},92/1),({16,11,47,34},124/1),({16,11,46,35},151/1),({16,11,45,36},170/1),({16,11,44,37},172/1),({16,11,43,38},155/1),({16,11,42,39},117/1),({16,11,41,40},63/1),({15,12,55,26},1/1),({15,12,54,27},3/1),({15,12,53,28},8/1),({15,12,52,29},18/1),({15,12,51,30},36/1),({15,12,50,31},60/1),({15,12,49,32},94/1),({15,12,48,33},133/1),({15,12,47,34},173/1),({15,12,46,35},207/1),({15,12,45,36},228/1),({15,12,44,37},226/1),({15,12,43,38},202/1),({15,12,42,39},152/1),({15,12,41,40},81/1),({14,13,55,26},1/1),({14,13,54,27},3/1),({14,13,53,28},8/1),({14,13,52,29},16/1),({14,13,51,30},30/1),({14,13,50,31},50/1),({14,13,49,32},75/1),({14,13,48,33},103/1),({14,13,47,34},134/1),({14,13,46,35},157/1),({14,13,45,36},170/1),({14,13,44,37},169/1),({14,13,43,38},149/1),({14,13,42,39},111/1),({14,13,41,40},61/1)}, (12,2) => {}, (14,0) => {}, (14,1) => {({20,11,52,41},1/1),({20,11,51,42},1/1),({20,11,50,43},2/1),({20,11,49,44},2/1),({20,11,48,45},2/1),({20,11,47,46},1/1),({19,12,56,37},1/1),({19,12,55,38},2/1),({19,12,54,39},4/1),({19,12,53,40},7/1),({19,12,52,41},10/1),({19,12,51,42},12/1),({19,12,50,43},14/1),({19,12,49,44},13/1),({19,12,48,45},10/1),({19,12,47,46},6/1),({18,13,58,35},1/1),({18,13,57,36},2/1),({18,13,56,37},5/1),({18,13,55,38},9/1),({18,13,54,39},16/1),({18,13,53,40},23/1),({18,13,52,41},31/1),({18,13,51,42},36/1),({18,13,50,43},39/1),({18,13,49,44},36/1),({18,13,48,45},28/1),({18,13,47,46},15/1),({17,14,58,35},2/1),({17,14,57,36},5/1),({17,14,56,37},10/1),({17,14,55,38},18/1),({17,14,54,39},29/1),({17,14,53,40},40/1),({17,14,52,41},51/1),({17,14,51,42},58/1),({17,14,50,43},60/1),({17,14,49,44},55/1),({17,14,48,45},42/1),({17,14,47,46},22/1),({16,15,59,34},1/1),({16,15,58,35},2/1),({16,15,57,36},5/1),({16,15,56,37},10/1),({16,15,55,38},16/1),({16,15,54,39},25/1),({16,15,53,40},34/1),({16,15,52,41},42/1),({16,15,51,42},47/1),({16,15,50,43},49/1),({16,15,49,44},44/1),({16,15,48,45},33/1),({16,15,47,46},18/1)}, (16,0) => {}, (14,2) => {}, (16,1) => {({20,15,59,46},1/1),({20,15,58,47},1/1),({20,15,57,48},2/1),({20,15,56,49},2/1),({20,15,55,50},3/1),({20,15,54,51},2/1),({20,15,53,52},1/1),({19,16,61,44},1/1),({19,16,60,45},2/1),({19,16,59,46},3/1),({19,16,58,47},4/1),({19,16,57,48},5/1),({19,16,56,49},6/1),({19,16,55,50},6/1),({19,16,54,51},4/1),({19,16,53,52},2/1),({18,17,61,44},1/1),({18,17,60,45},2/1),({18,17,59,46},3/1),({18,17,58,47},4/1),({18,17,57,48},5/1),({18,17,56,49},6/1),({18,17,55,50},6/1),({18,17,54,51},4/1),({18,17,53,52},2/1)}, (18,0) => {}, (16,2) => {}, (18,1) => {({20,19,62,55},1/1)}, (18,2) => {}, (20,1) => {}, (1,0) => {({3,0,8,1},1/1),({3,0,7,2},1/1),({3,0,6,3},1/1),({2,1,9,0},1/1),({2,1,8,1},1/1),({2,1,7,2},1/1),({2,1,6,3},1/1)}, (1,1) => {}, (3,0) => {({7,0,15,6},1/1),({7,0,13,8},1/1),({7,0,12,9},1/1),({6,1,18,3},1/1),({6,1,17,4},2/1),({6,1,16,5},3/1),({6,1,15,6},4/1),({6,1,14,7},5/1),({6,1,13,8},5/1),({6,1,12,9},4/1),({6,1,11,10},2/1),({5,2,19,2},1/1),({5,2,18,3},2/1),({5,2,17,4},5/1),({5,2,16,5},6/1),({5,2,15,6},8/1),({5,2,14,7},10/1),({5,2,13,8},10/1),({5,2,12,9},7/1),({5,2,11,10},4/1),({4,3,19,2},1/1),({4,3,18,3},3/1),({4,3,17,4},4/1),({4,3,16,5},6/1),({4,3,15,6},8/1),({4,3,14,7},9/1),({4,3,13,8},9/1),({4,3,12,9},7/1),({4,3,11,10},3/1)}, (3,1) => {}, (3,2) => {}, (5,0) => {({9,2,26,7},1/1),({9,2,25,8},1/1),({9,2,24,9},3/1),({9,2,23,10},3/1),({9,2,22,11},5/1),({9,2,21,12},5/1),({9,2,20,13},6/1),({9,2,19,14},5/1),({9,2,18,15},4/1),({9,2,17,16},2/1),({8,3,27,6},1/1),({8,3,26,7},2/1),({8,3,25,8},4/1),({8,3,24,9},7/1),({8,3,23,10},11/1),({8,3,22,11},14/1),({8,3,21,12},17/1),({8,3,20,13},18/1),({8,3,19,14},17/1),({8,3,18,15},13/1),({8,3,17,16},7/1),({7,4,28,5},1/1),({7,4,27,6},2/1),({7,4,26,7},5/1),({7,4,25,8},8/1),({7,4,24,9},14/1),({7,4,23,10},19/1),({7,4,22,11},25/1),({7,4,21,12},28/1),({7,4,20,13},30/1),({7,4,19,14},27/1),({7,4,18,15},21/1),({7,4,17,16},11/1),({6,5,28,5},1/1),({6,5,27,6},2/1),({6,5,26,7},4/1),({6,5,25,8},7/1),({6,5,24,9},11/1),({6,5,23,10},16/1),({6,5,22,11},20/1),({6,5,21,12},23/1),({6,5,20,13},24/1),({6,5,19,14},22/1),({6,5,18,15},17/1),({6,5,17,16},9/1)}, (5,1) => {({13,0,20,19},1/1),({12,1,25,14},1/1),({12,1,24,15},1/1),({12,1,23,16},1/1),({12,1,22,17},1/1),({12,1,21,18},1/1),({12,1,20,19},1/1),({11,2,25,14},1/1),({11,2,24,15},2/1),({11,2,23,16},2/1),({11,2,22,17},2/1),({11,2,21,18},2/1),({11,2,20,19},1/1),({10,3,29,10},1/1),({10,3,28,11},1/1),({10,3,27,12},2/1),({10,3,26,13},2/1),({10,3,25,14},3/1),({10,3,24,15},3/1),({10,3,23,16},4/1),({10,3,22,17},3/1),({10,3,21,18},2/1),({10,3,20,19},1/1),({9,4,28,11},1/1),({9,4,27,12},1/1),({9,4,26,13},2/1),({9,4,25,14},2/1),({9,4,24,15},3/1),({9,4,23,16},3/1),({9,4,22,17},3/1),({9,4,21,18},2/1),({9,4,20,19},1/1),({8,5,27,12},1/1),({8,5,26,13},1/1),({8,5,25,14},2/1),({8,5,24,15},2/1),({8,5,23,16},2/1),({8,5,22,17},2/1),({8,5,21,18},2/1),({8,5,20,19},1/1),({7,6,26,13},1/1),({7,6,25,14},1/1),({7,6,24,15},1/1),({7,6,23,16},1/1),({7,6,22,17},1/1),({7,6,21,18},1/1),({7,6,20,19},1/1)}};
--dw stands for dominant weights
dw1326 = new HashTable from {(5,2) => {}, (7,0) => {({11,4,33,12},1/1)}, (7,1) => {({15,2,31,20},1/1),({14,3,34,17},1/1),({13,4,36,15},1/1),({12,5,38,13},1/1),({11,6,39,12},1/1),({10,7,40,11},1/1)}, (9,0) => {}, (7,2) => {}, (9,1) => {({17,4,38,25},1/1),({16,5,41,22},1/1),({15,6,43,20},1/1),({14,7,45,18},1/1),({13,8,46,17},1/1),({12,9,47,16},1/1)}, (9,2) => {}, (11,0) => {}, (11,1) => {({19,6,41,34},1/1),({18,7,46,29},1/1),({17,8,49,26},1/1),({16,9,50,25},2/1),({15,10,52,23},1/1),({13,12,53,22},1/1)}, (11,2) => {}, (13,0) => {}, (13,1) => {({20,9,47,40},1/1),({19,10,52,35},1/1),({18,11,55,32},1/1),({17,12,56,31},1/1),({16,13,57,30},1/1)}, (13,2) => {}, (15,0) => {}, (15,1) => {({20,13,56,43},1/1),({19,14,59,40},1/1),({18,15,60,39},1/1)}, (15,2) => {}, (17,0) => {}, (17,1) => {({20,17,61,50},1/1),({19,18,62,49},1/1)}, (17,2) => {}, (19,1) => {}, (0,0) => {({1,0,3,0},1/1)}, (2,0) => {({5,0,12,3},1/1),({4,1,14,1},1/1)}, (2,1) => {}, (2,2) => {}, (4,0) => {({8,1,21,6},1/1),({7,2,23,4},1/1),({5,4,24,3},1/1)}, (4,1) => {({11,0,19,14},1/1)}, (4,2) => {}, (6,0) => {({10,3,30,9},1/1),({8,5,32,7},1/1)}, (6,1) => {({14,1,26,19},1/1),({13,2,30,15},1/1),({12,3,31,14},1/1),({11,4,34,11},1/1),({9,6,36,9},1/1)}, (6,2) => {}, (8,0) => {}, (8,1) => {({16,3,35,22},1/1),({15,4,37,20},2/1),({14,5,40,17},1/1),({13,6,41,16},2/1),({12,7,43,14},1/1),({10,9,44,13},1/1)}, (8,2) => {}, (10,0) => {}, (10,1) => {({18,5,40,29},1/1),({17,6,44,25},1/1),({16,7,46,23},1/1),({15,8,48,21},1/1),({14,9,49,20},1/1),({13,10,50,19},1/1)}, (10,2) => {}, (12,0) => {}, (12,1) => {({20,7,41,40},1/1),({19,8,47,34},1/1),({18,9,51,30},1/1),({17,10,53,28},1/1),({16,11,54,27},1/1),({15,12,55,26},1/1)}, (12,2) => {}, (14,0) => {}, (14,1) => {({20,11,52,41},1/1),({19,12,56,37},1/1),({18,13,58,35},1/1),({16,15,59,34},1/1)}, (16,0) => {}, (14,2) => {}, (16,1) => {({20,15,59,46},1/1),({19,16,61,44},1/1)}, (18,0) => {}, (16,2) => {}, (18,1) => {({20,19,62,55},1/1)}, (18,2) => {}, (20,1) => {}, (1,0) => {({3,0,8,1},1/1),({2,1,9,0},1/1)}, (1,1) => {}, (3,0) => {({7,0,15,6},1/1),({6,1,18,3},1/1),({5,2,19,2},1/1)}, (3,1) => {}, (5,0) => {({9,2,26,7},1/1),({8,3,27,6},1/1),({7,4,28,5},1/1)}, (3,2) => {}, (5,1) => {({13,0,20,19},1/1),({12,1,25,14},1/1),({10,3,29,10},1/1)}};
--dmw stands for dominant monomial weights
dmw1326 = new HashTable from {(5,2) => {}, (7,0) => {}, (7,1) => {{13,4,35,16},{15,2,31,20},{14,3,34,17},{16,1,26,25}}, (9,0) => {}, (7,2) => {}, (9,1) => {{17,4,38,25},{18,3,35,28},{15,6,42,21},{14,7,43,20},{16,5,40,23}}, (9,2) => {}, (11,0) => {}, (11,1) => {{17,8,47,28},{19,6,43,32},{16,9,48,27},{20,5,40,35},{18,7,45,30},{15,10,49,26}}, (11,2) => {}, (13,0) => {}, (13,1) => {{20,9,49,38},{22,7,41,46},{19,10,50,37},{21,8,46,41},{18,11,52,35},{16,13,53,34}}, (13,2) => {}, (15,0) => {}, (15,1) => {{20,13,54,45},{19,14,55,44},{21,12,53,46},{22,11,50,49}}, (15,2) => {}, (17,0) => {}, (17,1) => {{22,15,55,56},{21,16,56,55}}, (17,2) => {}, (19,1) => {}, (0,0) => {}, (2,0) => {}, (2,1) => {}, (2,2) => {}, (4,0) => {}, (4,1) => {{11,0,19,14}}, (4,2) => {}, (6,0) => {}, (6,1) => {{13,2,30,15},{14,1,26,19},{15,0,20,25}}, (6,2) => {}, (8,0) => {}, (8,1) => {{15,4,37,20},{16,3,35,22},{17,2,31,26},{14,5,39,18}}, (8,2) => {}, (10,0) => {}, (10,1) => {{18,5,40,29},{19,4,38,31},{16,7,44,25},{15,8,46,23},{17,6,43,26}}, (10,2) => {}, (12,0) => {}, (12,1) => {{16,11,51,30},{20,7,45,36},{18,9,49,32},{17,10,50,31},{19,8,47,34},{21,6,41,40}}, (12,2) => {}, (14,0) => {}, (14,1) => {{21,10,50,43},{22,9,46,47},{20,11,52,41},{19,12,53,40},{18,13,54,39}}, (16,0) => {}, (14,2) => {}, (16,1) => {{21,14,55,50},{22,13,53,52},{19,16,56,49}}, (18,0) => {}, (16,2) => {}, (18,1) => {{22,17,56,61}}, (18,2) => {{21,18,50,67}}, (20,1) => {}, (1,0) => {}, (1,1) => {}, (3,0) => {}, (3,1) => {}, (5,0) => {}, (3,2) => {}, (5,1) => {{13,0,20,19},{12,1,25,14}}};
end;