A := QQ[t_0,t_1,t_2,t_3];
--tb stands for Total Betti numbers
tb0426 = new HashTable from {(5,2) => 0, (7,0) => 0, (7,1) => 335478/1, (7,2) => 0, (9,0) => 0, (9,1) => 612612/1, (9,2) => 0, (11,0) => 0, (11,1) => 436254/1, (11,2) => 0, (13,0) => 0, (13,1) => 120360/1, (15,0) => 0, (13,2) => 0, (15,1) => 10863/1, (17,0) => 0, (15,2) => 0, (17,1) => 166/1, (17,2) => 0, (19,2) => 0, (0,0) => 5/1, (2,0) => 459/1, (2,1) => 288/1, (2,2) => 0, (4,0) => 3060/1, (4,1) => 11424/1, (4,2) => 0, (6,0) => 0, (6,1) => 166464/1, (6,2) => 0, (8,0) => 0, (8,1) => 512720/1, (8,2) => 0, (10,0) => 0, (10,1) => 579904/1, (10,2) => 0, (12,0) => 0, (12,1) => 259488/1, (12,2) => 0, (14,0) => 0, (14,1) => 42432/1, (16,0) => 0, (14,2) => 0, (16,1) => 1848/1, (18,0) => 0, (16,2) => 0, (18,1) => 0, (18,2) => 1/1, (1,0) => 72/1, (1,1) => 17/1, (3,0) => 1632/1, (3,1) => 2295/1, (3,2) => 0, (5,0) => 0, (5,1) => 58344/1};
--sb represents the betti numbers as sums of Schur functors
sb0426 = new HashTable from {(5,2) => {}, (7,0) => {}, (7,1) => {({14,2,31,21},2/1),({14,2,30,22},1/1),({14,2,29,23},2/1),({14,2,28,24},1/1),({14,2,27,25},2/1),({13,3,36,16},1/1),({13,3,35,17},2/1),({13,3,34,18},4/1),({13,3,33,19},6/1),({13,3,32,20},10/1),({13,3,31,21},12/1),({13,3,30,22},15/1),({13,3,29,23},14/1),({13,3,28,24},13/1),({13,3,27,25},9/1),({13,3,26,26},4/1),({12,4,39,13},1/1),({12,4,38,14},2/1),({12,4,37,15},5/1),({12,4,36,16},9/1),({12,4,35,17},18/1),({12,4,34,18},26/1),({12,4,33,19},39/1),({12,4,32,20},47/1),({12,4,31,21},60/1),({12,4,30,22},60/1),({12,4,29,23},62/1),({12,4,28,24},47/1),({12,4,27,25},36/1),({12,4,26,26},9/1),({11,5,40,12},2/1),({11,5,39,13},4/1),({11,5,38,14},11/1),({11,5,37,15},21/1),({11,5,36,16},37/1),({11,5,35,17},56/1),({11,5,34,18},84/1),({11,5,33,19},108/1),({11,5,32,20},135/1),({11,5,31,21},149/1),({11,5,30,22},158/1),({11,5,29,23},144/1),({11,5,28,24},123/1),({11,5,27,25},76/1),({11,5,26,26},29/1),({10,6,42,10},1/1),({10,6,41,11},3/1),({10,6,40,12},6/1),({10,6,39,13},15/1),({10,6,38,14},27/1),({10,6,37,15},49/1),({10,6,36,16},75/1),({10,6,35,17},115/1),({10,6,34,18},151/1),({10,6,33,19},199/1),({10,6,32,20},228/1),({10,6,31,21},258/1),({10,6,30,22},253/1),({10,6,29,23},243/1),({10,6,28,24},186/1),({10,6,27,25},131/1),({10,6,26,26},38/1),({9,7,42,10},1/1),({9,7,41,11},3/1),({9,7,40,12},9/1),({9,7,39,13},17/1),({9,7,38,14},34/1),({9,7,37,15},54/1),({9,7,36,16},88/1),({9,7,35,17},121/1),({9,7,34,18},167/1),({9,7,33,19},202/1),({9,7,32,20},243/1),({9,7,31,21},256/1),({9,7,30,22},265/1),({9,7,29,23},233/1),({9,7,28,24},197/1),({9,7,27,25},119/1),({9,7,26,26},49/1),({8,8,43,9},1/1),({8,8,42,10},1/1),({8,8,41,11},3/1),({8,8,40,12},5/1),({8,8,39,13},11/1),({8,8,38,14},16/1),({8,8,37,15},30/1),({8,8,36,16},39/1),({8,8,35,17},61/1),({8,8,34,18},73/1),({8,8,33,19},97/1),({8,8,32,20},102/1),({8,8,31,21},122/1),({8,8,30,22},109/1),({8,8,29,23},110/1),({8,8,28,24},77/1),({8,8,27,25},62/1),({8,8,26,26},11/1)}, (7,2) => {}, (9,0) => {}, (9,1) => {({16,4,40,24},1/1),({16,4,39,25},1/1),({16,4,38,26},3/1),({16,4,37,27},3/1),({16,4,36,28},5/1),({16,4,35,29},4/1),({16,4,34,30},5/1),({16,4,33,31},2/1),({16,4,32,32},2/1),({15,5,43,21},1/1),({15,5,42,22},2/1),({15,5,41,23},6/1),({15,5,40,24},10/1),({15,5,39,25},17/1),({15,5,38,26},23/1),({15,5,37,27},30/1),({15,5,36,28},33/1),({15,5,35,29},34/1),({15,5,34,30},28/1),({15,5,33,31},20/1),({15,5,32,32},6/1),({14,6,45,19},1/1),({14,6,44,20},4/1),({14,6,43,21},9/1),({14,6,42,22},20/1),({14,6,41,23},33/1),({14,6,40,24},55/1),({14,6,39,25},75/1),({14,6,38,26},101/1),({14,6,37,27},115/1),({14,6,36,28},128/1),({14,6,35,29},118/1),({14,6,34,30},104/1),({14,6,33,31},64/1),({14,6,32,32},27/1),({13,7,47,17},1/1),({13,7,46,18},3/1),({13,7,45,19},9/1),({13,7,44,20},19/1),({13,7,43,21},39/1),({13,7,42,22},66/1),({13,7,41,23},107/1),({13,7,40,24},152/1),({13,7,39,25},206/1),({13,7,38,26},250/1),({13,7,37,27},288/1),({13,7,36,28},295/1),({13,7,35,29},282/1),({13,7,34,30},227/1),({13,7,33,31},154/1),({13,7,32,32},50/1),({12,8,48,16},1/1),({12,8,47,17},3/1),({12,8,46,18},10/1),({12,8,45,19},21/1),({12,8,44,20},45/1),({12,8,43,21},77/1),({12,8,42,22},130/1),({12,8,41,23},189/1),({12,8,40,24},269/1),({12,8,39,25},338/1),({12,8,38,26},414/1),({12,8,37,27},449/1),({12,8,36,28},470/1),({12,8,35,29},424/1),({12,8,34,30},357/1),({12,8,33,31},222/1),({12,8,32,32},87/1),({11,9,49,15},1/1),({11,9,48,16},2/1),({11,9,47,17},7/1),({11,9,46,18},14/1),({11,9,45,19},31/1),({11,9,44,20},54/1),({11,9,43,21},95/1),({11,9,42,22},143/1),({11,9,41,23},213/1),({11,9,40,24},280/1),({11,9,39,25},361/1),({11,9,38,26},415/1),({11,9,37,27},464/1),({11,9,36,28},458/1),({11,9,35,29},432/1),({11,9,34,30},338/1),({11,9,33,31},231/1),({11,9,32,32},71/1),({10,10,48,16},1/1),({10,10,47,17},2/1),({10,10,46,18},7/1),({10,10,45,19},12/1),({10,10,44,20},26/1),({10,10,43,21},39/1),({10,10,42,22},66/1),({10,10,41,23},88/1),({10,10,40,24},126/1),({10,10,39,25},148/1),({10,10,38,26},184/1),({10,10,37,27},188/1),({10,10,36,28},203/1),({10,10,35,29},172/1),({10,10,34,30},153/1),({10,10,33,31},86/1),({10,10,32,32},41/1)}, (9,2) => {}, (11,0) => {}, (11,1) => {({18,6,45,31},1/1),({18,6,44,32},1/1),({18,6,43,33},2/1),({18,6,42,34},2/1),({18,6,41,35},4/1),({18,6,40,36},2/1),({18,6,39,37},2/1),({17,7,48,28},1/1),({17,7,47,29},2/1),({17,7,46,30},6/1),({17,7,45,31},10/1),({17,7,44,32},15/1),({17,7,43,33},19/1),({17,7,42,34},23/1),({17,7,41,35},23/1),({17,7,40,36},21/1),({17,7,39,37},13/1),({17,7,38,38},5/1),({16,8,50,26},1/1),({16,8,49,27},5/1),({16,8,48,28},10/1),({16,8,47,29},21/1),({16,8,46,30},33/1),({16,8,45,31},53/1),({16,8,44,32},68/1),({16,8,43,33},86/1),({16,8,42,34},90/1),({16,8,41,35},93/1),({16,8,40,36},74/1),({16,8,39,37},54/1),({16,8,38,38},15/1),({15,9,52,24},1/1),({15,9,51,25},3/1),({15,9,50,26},10/1),({15,9,49,27},20/1),({15,9,48,28},41/1),({15,9,47,29},68/1),({15,9,46,30},105/1),({15,9,45,31},143/1),({15,9,44,32},186/1),({15,9,43,33},212/1),({15,9,42,34},229/1),({15,9,41,35},214/1),({15,9,40,36},181/1),({15,9,39,37},117/1),({15,9,38,38},45/1),({14,10,53,23},1/1),({14,10,52,24},4/1),({14,10,51,25},12/1),({14,10,50,26},24/1),({14,10,49,27},50/1),({14,10,48,28},83/1),({14,10,47,29},134/1),({14,10,46,30},188/1),({14,10,45,31},256/1),({14,10,44,32},306/1),({14,10,43,33},356/1),({14,10,42,34},360/1),({14,10,41,35},347/1),({14,10,40,36},276/1),({14,10,39,37},191/1),({14,10,38,38},58/1),({13,11,54,22},1/1),({13,11,53,23},2/1),({13,11,52,24},7/1),({13,11,51,25},15/1),({13,11,50,26},33/1),({13,11,49,27},57/1),({13,11,48,28},98/1),({13,11,47,29},144/1),({13,11,46,30},208/1),({13,11,45,31},263/1),({13,11,44,32},324/1),({13,11,43,33},353/1),({13,11,42,34},371/1),({13,11,41,35},336/1),({13,11,40,36},283/1),({13,11,39,37},176/1),({13,11,38,38},69/1),({12,12,53,23},2/1),({12,12,52,24},3/1),({12,12,51,25},9/1),({12,12,50,26},15/1),({12,12,49,27},30/1),({12,12,48,28},43/1),({12,12,47,29},70/1),({12,12,46,30},88/1),({12,12,45,31},122/1),({12,12,44,32},136/1),({12,12,43,33},160/1),({12,12,42,34},152/1),({12,12,41,35},155/1),({12,12,40,36},112/1),({12,12,39,37},85/1),({12,12,38,38},20/1)}, (11,2) => {}, (13,0) => {}, (13,1) => {({20,8,46,42},1/1),({19,9,51,37},1/1),({19,9,50,38},2/1),({19,9,49,39},3/1),({19,9,48,40},4/1),({19,9,47,41},5/1),({19,9,46,42},4/1),({19,9,45,43},3/1),({19,9,44,44},1/1),({18,10,54,34},1/1),({18,10,53,35},2/1),({18,10,52,36},6/1),({18,10,51,37},9/1),({18,10,50,38},16/1),({18,10,49,39},19/1),({18,10,48,40},25/1),({18,10,47,41},23/1),({18,10,46,42},23/1),({18,10,45,43},13/1),({18,10,44,44},7/1),({17,11,55,33},3/1),({17,11,54,34},6/1),({17,11,53,35},14/1),({17,11,52,36},24/1),({17,11,51,37},38/1),({17,11,50,38},50/1),({17,11,49,39},65/1),({17,11,48,40},69/1),({17,11,47,41},70/1),({17,11,46,42},58/1),({17,11,45,43},41/1),({17,11,44,44},12/1),({16,12,57,31},1/1),({16,12,56,32},4/1),({16,12,55,33},8/1),({16,12,54,34},19/1),({16,12,53,35},32/1),({16,12,52,36},53/1),({16,12,51,37},73/1),({16,12,50,38},100/1),({16,12,49,39},113/1),({16,12,48,40},128/1),({16,12,47,41},118/1),({16,12,46,42},104/1),({16,12,45,43},64/1),({16,12,44,44},28/1),({15,13,57,31},2/1),({15,13,56,32},5/1),({15,13,55,33},13/1),({15,13,54,34},23/1),({15,13,53,35},42/1),({15,13,52,36},61/1),({15,13,51,37},88/1),({15,13,50,38},107/1),({15,13,49,39},129/1),({15,13,48,40},131/1),({15,13,47,41},130/1),({15,13,46,42},102/1),({15,13,45,43},73/1),({15,13,44,44},21/1),({14,14,58,30},1/1),({14,14,57,31},1/1),({14,14,56,32},4/1),({14,14,55,33},6/1),({14,14,54,34},13/1),({14,14,53,35},17/1),({14,14,52,36},31/1),({14,14,51,37},36/1),({14,14,50,38},51/1),({14,14,49,39},53/1),({14,14,48,40},63/1),({14,14,47,41},51/1),({14,14,46,42},52/1),({14,14,45,43},26/1),({14,14,44,44},15/1)}, (15,0) => {}, (13,2) => {}, (15,1) => {({20,12,55,45},1/1),({20,12,54,46},1/1),({20,12,53,47},2/1),({20,12,52,48},1/1),({20,12,51,49},2/1),({19,13,58,42},1/1),({19,13,57,43},2/1),({19,13,56,44},4/1),({19,13,55,45},5/1),({19,13,54,46},7/1),({19,13,53,47},7/1),({19,13,52,48},7/1),({19,13,51,49},4/1),({19,13,50,50},2/1),({18,14,59,41},2/1),({18,14,58,42},3/1),({18,14,57,43},7/1),({18,14,56,44},9/1),({18,14,55,45},15/1),({18,14,54,46},15/1),({18,14,53,47},18/1),({18,14,52,48},13/1),({18,14,51,49},12/1),({18,14,50,50},2/1),({17,15,60,40},2/1),({17,15,59,41},3/1),({17,15,58,42},7/1),({17,15,57,43},9/1),({17,15,56,44},15/1),({17,15,55,45},17/1),({17,15,54,46},21/1),({17,15,53,47},18/1),({17,15,52,48},18/1),({17,15,51,49},10/1),({17,15,50,50},6/1),({16,16,59,41},2/1),({16,16,58,42},2/1),({16,16,57,43},5/1),({16,16,56,44},5/1),({16,16,55,45},10/1),({16,16,54,46},8/1),({16,16,53,47},11/1),({16,16,52,48},6/1),({16,16,51,49},8/1)}, (17,0) => {}, (15,2) => {}, (17,1) => {({20,16,60,52},1/1),({20,16,58,54},1/1),({20,16,56,56},1/1),({19,17,61,51},1/1),({19,17,59,53},1/1),({19,17,57,55},1/1),({18,18,62,50},1/1),({18,18,60,52},1/1),({18,18,58,54},1/1),({18,18,56,56},1/1)}, (17,2) => {}, (19,2) => {}, (0,0) => {({0,0,4,0},1/1)}, (2,0) => {({4,0,13,3},1/1),({4,0,12,4},1/1),({4,0,11,5},2/1),({4,0,10,6},1/1),({4,0,9,7},2/1),({3,1,14,2},1/1),({3,1,13,3},1/1),({3,1,12,4},2/1),({3,1,11,5},1/1),({3,1,10,6},2/1),({3,1,9,7},1/1),({3,1,8,8},1/1),({2,2,13,3},1/1),({2,2,12,4},1/1),({2,2,11,5},2/1),({2,2,10,6},1/1),({2,2,9,7},2/1)}, (2,1) => {({4,2,21,1},1/1),({4,2,20,2},1/1),({4,2,19,3},1/1),({4,2,18,4},1/1),({4,2,17,5},1/1),({4,2,16,6},1/1)}, (2,2) => {}, (4,0) => {({8,0,18,10},1/1),({8,0,16,12},1/1),({8,0,14,14},1/1),({7,1,21,7},1/1),({7,1,20,8},1/1),({7,1,19,9},2/1),({7,1,18,10},2/1),({7,1,17,11},3/1),({7,1,16,12},2/1),({7,1,15,13},2/1),({6,2,22,6},1/1),({6,2,21,7},1/1),({6,2,20,8},3/1),({6,2,19,9},3/1),({6,2,18,10},6/1),({6,2,17,11},4/1),({6,2,16,12},6/1),({6,2,15,13},2/1),({6,2,14,14},3/1),({5,3,23,5},1/1),({5,3,22,6},1/1),({5,3,21,7},3/1),({5,3,20,8},3/1),({5,3,19,9},6/1),({5,3,18,10},5/1),({5,3,17,11},7/1),({5,3,16,12},4/1),({5,3,15,13},5/1),({4,4,22,6},1/1),({4,4,20,8},2/1),({4,4,19,9},1/1),({4,4,18,10},4/1),({4,4,17,11},1/1),({4,4,16,12},4/1),({4,4,14,14},3/1)}, (4,1) => {({8,2,28,6},1/1),({8,2,27,7},1/1),({8,2,26,8},2/1),({8,2,25,9},3/1),({8,2,24,10},3/1),({8,2,23,11},3/1),({8,2,22,12},3/1),({8,2,21,13},2/1),({8,2,20,14},1/1),({8,2,19,15},1/1),({7,3,30,4},1/1),({7,3,29,5},2/1),({7,3,28,6},3/1),({7,3,27,7},5/1),({7,3,26,8},7/1),({7,3,25,9},8/1),({7,3,24,10},9/1),({7,3,23,11},9/1),({7,3,22,12},8/1),({7,3,21,13},7/1),({7,3,20,14},5/1),({7,3,19,15},3/1),({7,3,18,16},2/1),({7,3,17,17},1/1),({6,4,30,4},1/1),({6,4,29,5},2/1),({6,4,28,6},4/1),({6,4,27,7},6/1),({6,4,26,8},9/1),({6,4,25,9},11/1),({6,4,24,10},12/1),({6,4,23,11},12/1),({6,4,22,12},11/1),({6,4,21,13},9/1),({6,4,20,14},6/1),({6,4,19,15},4/1),({6,4,18,16},2/1),({6,4,17,17},1/1),({5,5,31,3},1/1),({5,5,30,4},1/1),({5,5,29,5},2/1),({5,5,28,6},3/1),({5,5,27,7},4/1),({5,5,26,8},5/1),({5,5,25,9},6/1),({5,5,24,10},6/1),({5,5,23,11},6/1),({5,5,22,12},6/1),({5,5,21,13},5/1),({5,5,20,14},4/1),({5,5,19,15},3/1),({5,5,18,16},2/1)}, (4,2) => {}, (6,0) => {}, (6,1) => {({13,1,26,20},1/1),({13,1,25,21},1/1),({12,2,31,15},1/1),({12,2,30,16},2/1),({12,2,29,17},3/1),({12,2,28,18},4/1),({12,2,27,19},5/1),({12,2,26,20},6/1),({12,2,25,21},5/1),({12,2,24,22},3/1),({12,2,23,23},1/1),({11,3,35,11},1/1),({11,3,34,12},2/1),({11,3,33,13},4/1),({11,3,32,14},7/1),({11,3,31,15},11/1),({11,3,30,16},16/1),({11,3,29,17},21/1),({11,3,28,18},25/1),({11,3,27,19},27/1),({11,3,26,20},26/1),({11,3,25,21},22/1),({11,3,24,22},15/1),({11,3,23,23},5/1),({10,4,37,9},1/1),({10,4,36,10},2/1),({10,4,35,11},5/1),({10,4,34,12},10/1),({10,4,33,13},18/1),({10,4,32,14},28/1),({10,4,31,15},40/1),({10,4,30,16},53/1),({10,4,29,17},65/1),({10,4,28,18},73/1),({10,4,27,19},75/1),({10,4,26,20},70/1),({10,4,25,21},57/1),({10,4,24,22},38/1),({10,4,23,23},13/1),({9,5,38,8},1/1),({9,5,37,9},3/1),({9,5,36,10},7/1),({9,5,35,11},14/1),({9,5,34,12},25/1),({9,5,33,13},40/1),({9,5,32,14},59/1),({9,5,31,15},80/1),({9,5,30,16},101/1),({9,5,29,17},119/1),({9,5,28,18},130/1),({9,5,27,19},131/1),({9,5,26,20},120/1),({9,5,25,21},97/1),({9,5,24,22},63/1),({9,5,23,23},22/1),({8,6,39,7},1/1),({8,6,38,8},2/1),({8,6,37,9},5/1),({8,6,36,10},10/1),({8,6,35,11},19/1),({8,6,34,12},31/1),({8,6,33,13},48/1),({8,6,32,14},68/1),({8,6,31,15},91/1),({8,6,30,16},112/1),({8,6,29,17},129/1),({8,6,28,18},139/1),({8,6,27,19},138/1),({8,6,26,20},126/1),({8,6,25,21},100/1),({8,6,24,22},65/1),({8,6,23,23},22/1),({7,7,38,8},1/1),({7,7,37,9},2/1),({7,7,36,10},5/1),({7,7,35,11},9/1),({7,7,34,12},15/1),({7,7,33,13},22/1),({7,7,32,14},31/1),({7,7,31,15},41/1),({7,7,30,16},50/1),({7,7,29,17},57/1),({7,7,28,18},61/1),({7,7,27,19},61/1),({7,7,26,20},55/1),({7,7,25,21},44/1),({7,7,24,22},28/1),({7,7,23,23},10/1)}, (6,2) => {}, (8,0) => {}, (8,1) => {({15,3,36,22},1/1),({15,3,35,23},2/1),({15,3,34,24},2/1),({15,3,33,25},3/1),({15,3,32,26},4/1),({15,3,31,27},3/1),({15,3,30,28},2/1),({15,3,29,29},1/1),({14,4,40,18},1/1),({14,4,39,19},2/1),({14,4,38,20},4/1),({14,4,37,21},8/1),({14,4,36,22},13/1),({14,4,35,23},18/1),({14,4,34,24},23/1),({14,4,33,25},26/1),({14,4,32,26},26/1),({14,4,31,27},23/1),({14,4,30,28},15/1),({14,4,29,29},5/1),({13,5,42,16},1/1),({13,5,41,17},3/1),({13,5,40,18},8/1),({13,5,39,19},16/1),({13,5,38,20},28/1),({13,5,37,21},44/1),({13,5,36,22},62/1),({13,5,35,23},80/1),({13,5,34,24},94/1),({13,5,33,25},101/1),({13,5,32,26},97/1),({13,5,31,27},81/1),({13,5,30,28},54/1),({13,5,29,29},19/1),({12,6,44,14},1/1),({12,6,43,15},3/1),({12,6,42,16},8/1),({12,6,41,17},17/1),({12,6,40,18},33/1),({12,6,39,19},57/1),({12,6,38,20},90/1),({12,6,37,21},128/1),({12,6,36,22},170/1),({12,6,35,23},209/1),({12,6,34,24},235/1),({12,6,33,25},244/1),({12,6,32,26},229/1),({12,6,31,27},187/1),({12,6,30,28},123/1),({12,6,29,29},44/1),({11,7,45,13},1/1),({11,7,44,14},3/1),({11,7,43,15},9/1),({11,7,42,16},20/1),({11,7,41,17},39/1),({11,7,40,18},69/1),({11,7,39,19},111/1),({11,7,38,20},164/1),({11,7,37,21},226/1),({11,7,36,22},289/1),({11,7,35,23},343/1),({11,7,34,24},380/1),({11,7,33,25},386/1),({11,7,32,26},356/1),({11,7,31,27},290/1),({11,7,30,28},189/1),({11,7,29,29},65/1),({10,8,46,12},1/1),({10,8,45,13},2/1),({10,8,44,14},6/1),({10,8,43,15},13/1),({10,8,42,16},27/1),({10,8,41,17},49/1),({10,8,40,18},82/1),({10,8,39,19},126/1),({10,8,38,20},181/1),({10,8,37,21},242/1),({10,8,36,22},303/1),({10,8,35,23},354/1),({10,8,34,24},385/1),({10,8,33,25},388/1),({10,8,32,26},355/1),({10,8,31,27},287/1),({10,8,30,28},186/1),({10,8,29,29},65/1),({9,9,45,13},1/1),({9,9,44,14},3/1),({9,9,43,15},6/1),({9,9,42,16},12/1),({9,9,41,17},23/1),({9,9,40,18},36/1),({9,9,39,19},56/1),({9,9,38,20},80/1),({9,9,37,21},105/1),({9,9,36,22},131/1),({9,9,35,23},153/1),({9,9,34,24},164/1),({9,9,33,25},165/1),({9,9,32,26},152/1),({9,9,31,27},120/1),({9,9,30,28},79/1),({9,9,29,29},28/1)}, (8,2) => {}, (10,0) => {}, (10,1) => {({17,5,43,27},1/1),({17,5,42,28},1/1),({17,5,41,29},2/1),({17,5,40,30},4/1),({17,5,39,31},4/1),({17,5,38,32},4/1),({17,5,37,33},4/1),({17,5,36,34},3/1),({17,5,35,35},1/1),({16,6,45,25},2/1),({16,6,44,26},5/1),({16,6,43,27},9/1),({16,6,42,28},15/1),({16,6,41,29},22/1),({16,6,40,30},28/1),({16,6,39,31},32/1),({16,6,38,32},32/1),({16,6,37,33},27/1),({16,6,36,34},19/1),({16,6,35,35},7/1),({15,7,48,22},1/1),({15,7,47,23},3/1),({15,7,46,24},8/1),({15,7,45,25},17/1),({15,7,44,26},31/1),({15,7,43,27},50/1),({15,7,42,28},72/1),({15,7,41,29},95/1),({15,7,40,30},113/1),({15,7,39,31},122/1),({15,7,38,32},118/1),({15,7,37,33},99/1),({15,7,36,34},66/1),({15,7,35,35},23/1),({14,8,49,21},2/1),({14,8,48,22},6/1),({14,8,47,23},16/1),({14,8,46,24},33/1),({14,8,45,25},59/1),({14,8,44,26},96/1),({14,8,43,27},143/1),({14,8,42,28},193/1),({14,8,41,29},240/1),({14,8,40,30},276/1),({14,8,39,31},289/1),({14,8,38,32},273/1),({14,8,37,33},226/1),({14,8,36,34},149/1),({14,8,35,35},51/1),({13,9,51,19},1/1),({13,9,50,20},3/1),({13,9,49,21},8/1),({13,9,48,22},19/1),({13,9,47,23},39/1),({13,9,46,24},71/1),({13,9,45,25},118/1),({13,9,44,26},180/1),({13,9,43,27},252/1),({13,9,42,28},329/1),({13,9,41,29},397/1),({13,9,40,30},443/1),({13,9,39,31},456/1),({13,9,38,32},425/1),({13,9,37,33},346/1),({13,9,36,34},227/1),({13,9,35,35},80/1),({12,10,51,19},1/1),({12,10,50,20},4/1),({12,10,49,21},11/1),({12,10,48,22},24/1),({12,10,47,23},47/1),({12,10,46,24},82/1),({12,10,45,25},132/1),({12,10,44,26},194/1),({12,10,43,27},266/1),({12,10,42,28},339/1),({12,10,41,29},403/1),({12,10,40,30},445/1),({12,10,39,31},452/1),({12,10,38,32},418/1),({12,10,37,33},339/1),({12,10,36,34},222/1),({12,10,35,35},77/1),({11,11,52,18},1/1),({11,11,51,19},1/1),({11,11,50,20},3/1),({11,11,49,21},7/1),({11,11,48,22},13/1),({11,11,47,23},23/1),({11,11,46,24},40/1),({11,11,45,25},61/1),({11,11,44,26},88/1),({11,11,43,27},120/1),({11,11,42,28},150/1),({11,11,41,29},176/1),({11,11,40,30},194/1),({11,11,39,31},195/1),({11,11,38,32},178/1),({11,11,37,33},146/1),({11,11,36,34},94/1),({11,11,35,35},32/1)}, (10,2) => {}, (12,0) => {}, (12,1) => {({19,7,46,36},1/1),({19,7,45,37},1/1),({19,7,44,38},1/1),({19,7,43,39},1/1),({19,7,42,40},1/1),({19,7,41,41},1/1),({18,8,50,32},1/1),({18,8,49,33},2/1),({18,8,48,34},4/1),({18,8,47,35},7/1),({18,8,46,36},10/1),({18,8,45,37},12/1),({18,8,44,38},12/1),({18,8,43,39},11/1),({18,8,42,40},8/1),({18,8,41,41},3/1),({17,9,52,30},1/1),({17,9,51,31},4/1),({17,9,50,32},9/1),({17,9,49,33},17/1),({17,9,48,34},27/1),({17,9,47,35},38/1),({17,9,46,36},48/1),({17,9,45,37},54/1),({17,9,44,38},54/1),({17,9,43,39},46/1),({17,9,42,40},31/1),({17,9,41,41},11/1),({16,10,54,28},1/1),({16,10,53,29},4/1),({16,10,52,30},10/1),({16,10,51,31},20/1),({16,10,50,32},37/1),({16,10,49,33},59/1),({16,10,48,34},85/1),({16,10,47,35},111/1),({16,10,46,36},132/1),({16,10,45,37},143/1),({16,10,44,38},138/1),({16,10,43,39},116/1),({16,10,42,40},77/1),({16,10,41,41},27/1),({15,11,55,27},1/1),({15,11,54,28},4/1),({15,11,53,29},11/1),({15,11,52,30},24/1),({15,11,51,31},45/1),({15,11,50,32},75/1),({15,11,49,33},113/1),({15,11,48,34},155/1),({15,11,47,35},195/1),({15,11,46,36},225/1),({15,11,45,37},237/1),({15,11,44,38},225/1),({15,11,43,39},186/1),({15,11,42,40},123/1),({15,11,41,41},43/1),({14,12,56,26},1/1),({14,12,55,27},3/1),({14,12,54,28},8/1),({14,12,53,29},17/1),({14,12,52,30},33/1),({14,12,51,31},57/1),({14,12,50,32},89/1),({14,12,49,33},128/1),({14,12,48,34},170/1),({14,12,47,35},209/1),({14,12,46,36},237/1),({14,12,45,37},246/1),({14,12,44,38},231/1),({14,12,43,39},190/1),({14,12,42,40},125/1),({14,12,41,41},44/1),({13,13,55,27},1/1),({13,13,54,28},3/1),({13,13,53,29},7/1),({13,13,52,30},14/1),({13,13,51,31},25/1),({13,13,50,32},39/1),({13,13,49,33},56/1),({13,13,48,34},74/1),({13,13,47,35},91/1),({13,13,46,36},102/1),({13,13,45,37},105/1),({13,13,44,38},98/1),({13,13,43,39},80/1),({13,13,42,40},53/1),({13,13,41,41},18/1)}, (12,2) => {}, (14,0) => {}, (14,1) => {({20,10,51,43},1/1),({20,10,50,44},1/1),({20,10,49,45},1/1),({20,10,48,46},1/1),({19,11,55,39},1/1),({19,11,54,40},2/1),({19,11,53,41},4/1),({19,11,52,42},6/1),({19,11,51,43},7/1),({19,11,50,44},8/1),({19,11,49,45},7/1),({19,11,48,46},5/1),({19,11,47,47},2/1),({18,12,57,37},1/1),({18,12,56,38},3/1),({18,12,55,39},6/1),({18,12,54,40},11/1),({18,12,53,41},17/1),({18,12,52,42},22/1),({18,12,51,43},26/1),({18,12,50,44},27/1),({18,12,49,45},23/1),({18,12,48,46},16/1),({18,12,47,47},6/1),({17,13,58,36},2/1),({17,13,57,37},5/1),({17,13,56,38},10/1),({17,13,55,39},18/1),({17,13,54,40},28/1),({17,13,53,41},38/1),({17,13,52,42},47/1),({17,13,51,43},52/1),({17,13,50,44},51/1),({17,13,49,45},44/1),({17,13,48,46},29/1),({17,13,47,47},10/1),({16,14,59,35},1/1),({16,14,58,36},3/1),({16,14,57,37},7/1),({16,14,56,38},13/1),({16,14,55,39},22/1),({16,14,54,40},33/1),({16,14,53,41},44/1),({16,14,52,42},53/1),({16,14,51,43},58/1),({16,14,50,44},57/1),({16,14,49,45},48/1),({16,14,48,46},32/1),({16,14,47,47},11/1),({15,15,59,35},1/1),({15,15,58,36},2/1),({15,15,57,37},4/1),({15,15,56,38},8/1),({15,15,55,39},12/1),({15,15,54,40},17/1),({15,15,53,41},22/1),({15,15,52,42},25/1),({15,15,51,43},27/1),({15,15,50,44},26/1),({15,15,49,45},21/1),({15,15,48,46},14/1),({15,15,47,47},5/1)}, (16,0) => {}, (14,2) => {}, (16,1) => {({20,14,58,48},1/1),({20,14,57,49},1/1),({20,14,56,50},1/1),({20,14,55,51},2/1),({20,14,54,52},1/1),({19,15,60,46},1/1),({19,15,59,47},2/1),({19,15,58,48},2/1),({19,15,57,49},3/1),({19,15,56,50},4/1),({19,15,55,51},3/1),({19,15,54,52},2/1),({19,15,53,53},1/1),({18,16,61,45},1/1),({18,16,60,46},2/1),({18,16,59,47},3/1),({18,16,58,48},4/1),({18,16,57,49},5/1),({18,16,56,50},6/1),({18,16,55,51},5/1),({18,16,54,52},3/1),({18,16,53,53},1/1),({17,17,61,45},1/1),({17,17,60,46},1/1),({17,17,59,47},1/1),({17,17,58,48},2/1),({17,17,57,49},2/1),({17,17,56,50},2/1),({17,17,55,51},2/1),({17,17,54,52},1/1)}, (18,0) => {}, (16,2) => {}, (18,1) => {}, (18,2) => {({20,20,62,62},1/1)}, (1,0) => {({2,0,9,1},1/1),({2,0,8,2},1/1),({2,0,7,3},1/1),({2,0,6,4},1/1)}, (1,1) => {({2,2,16,0},1/1)}, (3,0) => {({6,0,16,6},1/1),({6,0,15,7},1/1),({6,0,14,8},1/1),({6,0,13,9},2/1),({6,0,12,10},1/1),({5,1,18,4},1/1),({5,1,17,5},2/1),({5,1,16,6},2/1),({5,1,15,7},3/1),({5,1,14,8},4/1),({5,1,13,9},3/1),({5,1,12,10},2/1),({5,1,11,11},1/1),({4,2,18,4},1/1),({4,2,17,5},2/1),({4,2,16,6},3/1),({4,2,15,7},4/1),({4,2,14,8},5/1),({4,2,13,9},5/1),({4,2,12,10},3/1),({4,2,11,11},1/1),({3,3,19,3},1/1),({3,3,18,4},1/1),({3,3,17,5},1/1),({3,3,16,6},2/1),({3,3,15,7},2/1),({3,3,14,8},2/1),({3,3,13,9},2/1),({3,3,12,10},1/1)}, (3,1) => {({6,2,25,3},1/1),({6,2,24,4},1/1),({6,2,23,5},2/1),({6,2,22,6},2/1),({6,2,21,7},3/1),({6,2,20,8},2/1),({6,2,19,9},2/1),({6,2,18,10},1/1),({6,2,17,11},1/1),({5,3,26,2},1/1),({5,3,25,3},1/1),({5,3,24,4},2/1),({5,3,23,5},2/1),({5,3,22,6},3/1),({5,3,21,7},3/1),({5,3,20,8},3/1),({5,3,19,9},2/1),({5,3,18,10},2/1),({5,3,17,11},1/1),({5,3,16,12},1/1),({4,4,25,3},1/1),({4,4,24,4},1/1),({4,4,23,5},2/1),({4,4,22,6},2/1),({4,4,21,7},3/1),({4,4,20,8},2/1),({4,4,19,9},2/1),({4,4,18,10},1/1),({4,4,17,11},1/1)}, (3,2) => {}, (5,0) => {}, (5,1) => {({12,0,20,20},1/1),({11,1,25,15},1/1),({11,1,24,16},1/1),({11,1,23,17},1/1),({11,1,22,18},1/1),({11,1,21,19},1/1),({10,2,30,10},1/1),({10,2,29,11},1/1),({10,2,28,12},3/1),({10,2,27,13},3/1),({10,2,26,14},6/1),({10,2,25,15},5/1),({10,2,24,16},8/1),({10,2,23,17},5/1),({10,2,22,18},7/1),({10,2,21,19},2/1),({10,2,20,20},3/1),({9,3,33,7},1/1),({9,3,32,8},2/1),({9,3,31,9},4/1),({9,3,30,10},6/1),({9,3,29,11},11/1),({9,3,28,12},13/1),({9,3,27,13},19/1),({9,3,26,14},21/1),({9,3,25,15},24/1),({9,3,24,16},23/1),({9,3,23,17},24/1),({9,3,22,18},16/1),({9,3,21,19},13/1),({9,3,20,20},3/1),({8,4,34,6},1/1),({8,4,33,7},2/1),({8,4,32,8},5/1),({8,4,31,9},9/1),({8,4,30,10},16/1),({8,4,29,11},21/1),({8,4,28,12},32/1),({8,4,27,13},36/1),({8,4,26,14},45/1),({8,4,25,15},44/1),({8,4,24,16},48/1),({8,4,23,17},37/1),({8,4,22,18},36/1),({8,4,21,19},17/1),({8,4,20,20},10/1),({7,5,35,5},1/1),({7,5,34,6},2/1),({7,5,33,7},5/1),({7,5,32,8},8/1),({7,5,31,9},15/1),({7,5,30,10},21/1),({7,5,29,11},31/1),({7,5,28,12},37/1),({7,5,27,13},48/1),({7,5,26,14},50/1),({7,5,25,15},57/1),({7,5,24,16},50/1),({7,5,23,17},50/1),({7,5,22,18},34/1),({7,5,21,19},27/1),({7,5,20,20},5/1),({6,6,34,6},1/1),({6,6,33,7},1/1),({6,6,32,8},4/1),({6,6,31,9},4/1),({6,6,30,10},10/1),({6,6,29,11},11/1),({6,6,28,12},18/1),({6,6,27,13},17/1),({6,6,26,14},26/1),({6,6,25,15},19/1),({6,6,24,16},26/1),({6,6,23,17},16/1),({6,6,22,18},19/1),({6,6,21,19},5/1),({6,6,20,20},9/1)}};
--dw stands for dominant weights
dw0426 = new HashTable from {(5,2) => {}, (7,0) => {}, (7,1) => {({14,2,31,21},2/1),({13,3,36,16},1/1),({12,4,39,13},1/1),({11,5,40,12},2/1),({10,6,42,10},1/1),({8,8,43,9},1/1)}, (9,0) => {}, (7,2) => {}, (9,1) => {({16,4,40,24},1/1),({15,5,43,21},1/1),({14,6,45,19},1/1),({13,7,47,17},1/1),({12,8,48,16},1/1),({11,9,49,15},1/1)}, (9,2) => {}, (11,0) => {}, (11,1) => {({18,6,45,31},1/1),({17,7,48,28},1/1),({16,8,50,26},1/1),({15,9,52,24},1/1),({14,10,53,23},1/1),({13,11,54,22},1/1)}, (11,2) => {}, (13,0) => {}, (13,1) => {({20,8,46,42},1/1),({19,9,51,37},1/1),({18,10,54,34},1/1),({17,11,55,33},3/1),({16,12,57,31},1/1),({14,14,58,30},1/1)}, (13,2) => {}, (15,0) => {}, (15,1) => {({20,12,55,45},1/1),({19,13,58,42},1/1),({18,14,59,41},2/1),({17,15,60,40},2/1)}, (15,2) => {}, (17,0) => {}, (17,1) => {({20,16,60,52},1/1),({19,17,61,51},1/1),({18,18,62,50},1/1)}, (17,2) => {}, (19,2) => {}, (0,0) => {({0,0,4,0},1/1)}, (2,0) => {({4,0,13,3},1/1),({3,1,14,2},1/1)}, (2,1) => {({4,2,21,1},1/1)}, (2,2) => {}, (4,0) => {({8,0,18,10},1/1),({7,1,21,7},1/1),({6,2,22,6},1/1),({5,3,23,5},1/1)}, (4,1) => {({8,2,28,6},1/1),({7,3,30,4},1/1),({5,5,31,3},1/1)}, (4,2) => {}, (6,0) => {}, (6,1) => {({13,1,26,20},1/1),({12,2,31,15},1/1),({11,3,35,11},1/1),({10,4,37,9},1/1),({9,5,38,8},1/1),({8,6,39,7},1/1)}, (6,2) => {}, (8,0) => {}, (8,1) => {({15,3,36,22},1/1),({14,4,40,18},1/1),({13,5,42,16},1/1),({12,6,44,14},1/1),({11,7,45,13},1/1),({10,8,46,12},1/1)}, (8,2) => {}, (10,0) => {}, (10,1) => {({17,5,43,27},1/1),({16,6,45,25},2/1),({15,7,48,22},1/1),({14,8,49,21},2/1),({13,9,51,19},1/1),({11,11,52,18},1/1)}, (10,2) => {}, (12,0) => {}, (12,1) => {({19,7,46,36},1/1),({18,8,50,32},1/1),({17,9,52,30},1/1),({16,10,54,28},1/1),({15,11,55,27},1/1),({14,12,56,26},1/1)}, (12,2) => {}, (14,0) => {}, (14,1) => {({20,10,51,43},1/1),({19,11,55,39},1/1),({18,12,57,37},1/1),({17,13,58,36},2/1),({16,14,59,35},1/1)}, (16,0) => {}, (14,2) => {}, (16,1) => {({20,14,58,48},1/1),({19,15,60,46},1/1),({18,16,61,45},1/1)}, (18,0) => {}, (16,2) => {}, (18,1) => {}, (18,2) => {({20,20,62,62},1/1)}, (1,0) => {({2,0,9,1},1/1)}, (1,1) => {({2,2,16,0},1/1)}, (3,0) => {({6,0,16,6},1/1),({5,1,18,4},1/1),({3,3,19,3},1/1)}, (3,1) => {({6,2,25,3},1/1),({5,3,26,2},1/1)}, (5,0) => {}, (3,2) => {}, (5,1) => {({12,0,20,20},1/1),({11,1,25,15},1/1),({10,2,30,10},1/1),({9,3,33,7},1/1),({8,4,34,6},1/1),({7,5,35,5},1/1)}};
--dmw stands for dominant monomial weights
dmw0426 = new HashTable from {(5,2) => {}, (7,0) => {}, (7,1) => {{15,1,26,26},{14,2,31,21}}, (9,0) => {}, (7,2) => {}, (9,1) => {{16,4,38,26},{15,5,40,24},{17,3,35,29}}, (9,2) => {}, (11,0) => {}, (11,1) => {{19,5,40,36},{16,8,47,29},{17,7,45,31},{18,6,43,33}}, (11,2) => {}, (13,0) => {}, (13,1) => {{19,9,49,39},{18,10,50,38},{20,8,46,42},{17,11,52,36}}, (13,2) => {{22,8,41,53}}, (15,0) => {}, (15,1) => {{19,13,54,46},{18,14,55,45},{20,12,53,47}}, (15,2) => {{22,12,50,56}}, (17,0) => {}, (17,1) => {{20,16,56,56}}, (17,2) => {{22,16,55,63}}, (19,2) => {}, (0,0) => {}, (2,0) => {}, (2,1) => {}, (2,2) => {}, (4,0) => {}, (4,1) => {}, (4,2) => {}, (6,0) => {}, (6,1) => {{13,1,26,20},{14,0,20,26}}, (6,2) => {}, (8,0) => {}, (8,1) => {{15,3,35,23},{16,2,31,27},{14,4,36,22}}, (8,2) => {}, (10,0) => {}, (10,1) => {{17,5,40,30},{18,4,38,32},{15,7,44,26},{16,6,43,27}}, (10,2) => {}, (12,0) => {}, (12,1) => {{16,10,50,32},{18,8,47,35},{19,7,45,37},{20,6,41,41},{17,9,49,33}}, (12,2) => {{21,7,35,53}}, (14,0) => {}, (14,1) => {{18,12,53,41},{19,11,52,42},{20,10,50,44},{17,13,54,40}}, (16,0) => {}, (14,2) => {{22,10,46,54}}, (16,1) => {{20,14,55,51},{18,16,56,50}}, (18,0) => {}, (16,2) => {{22,14,53,59}}, (18,1) => {}, (18,2) => {{22,18,56,68}}, (1,0) => {}, (1,1) => {}, (3,0) => {}, (3,1) => {}, (5,0) => {}, (3,2) => {}, (5,1) => {{12,0,20,20}}};
end;