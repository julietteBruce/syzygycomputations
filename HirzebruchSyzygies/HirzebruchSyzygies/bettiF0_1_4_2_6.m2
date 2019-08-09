A := QQ[t_0,t_1,t_2,t_3];
--tb stands for Total Betti numbers
tb1426 = new HashTable from {(5,2) => 0, (7,0) => 60346/1, (7,1) => 31122/1, (7,2) => 0, (9,0) => 2002/1, (9,1) => 243100/1, (9,2) => 0, (11,0) => 0, (11,1) => 259896/1, (11,2) => 0, (13,0) => 0, (13,1) => 89352/1, (15,0) => 0, (13,2) => 0, (15,1) => 9894/1, (17,0) => 0, (15,2) => 0, (17,1) => 242/1, (17,2) => 0, (19,1) => 0, (0,0) => 10/1, (2,0) => 1278/1, (2,1) => 0, (2,2) => 0, (4,0) => 19176/1, (4,1) => 0, (4,2) => 0, (6,0) => 65702/1, (6,1) => 2002/1, (6,2) => 0, (8,0) => 23166/1, (8,1) => 128414/1, (8,2) => 0, (10,0) => 0, (10,1) => 294372/1, (10,2) => 0, (12,0) => 0, (12,1) => 174216/1, (12,2) => 0, (14,0) => 0, (14,1) => 34680/1, (16,0) => 0, (14,2) => 0, (16,1) => 1962/1, (18,0) => 0, (16,2) => 0, (18,1) => 14/1, (18,2) => 0, (20,1) => 0, (1,0) => 166/1, (1,1) => 0, (3,0) => 6018/1, (3,1) => 0, (3,2) => 0, (5,0) => 42840/1, (5,1) => 14/1};
--sb represents the betti numbers as sums of Schur functors
sb1426 = new HashTable from {(5,2) => {}, (7,0) => {({12,3,33,13},1/1),({12,3,32,14},1/1),({12,3,31,15},3/1),({12,3,30,16},3/1),({12,3,29,17},5/1),({12,3,28,18},4/1),({12,3,27,19},6/1),({12,3,26,20},4/1),({12,3,25,21},5/1),({12,3,24,22},2/1),({12,3,23,23},2/1),({11,4,34,12},1/1),({11,4,33,13},2/1),({11,4,32,14},4/1),({11,4,31,15},8/1),({11,4,30,16},12/1),({11,4,29,17},15/1),({11,4,28,18},19/1),({11,4,27,19},20/1),({11,4,26,20},19/1),({11,4,25,21},17/1),({11,4,24,22},11/1),({11,4,23,23},4/1),({10,5,36,10},1/1),({10,5,35,11},2/1),({10,5,34,12},4/1),({10,5,33,13},8/1),({10,5,32,14},13/1),({10,5,31,15},21/1),({10,5,30,16},28/1),({10,5,29,17},38/1),({10,5,28,18},43/1),({10,5,27,19},48/1),({10,5,26,20},44/1),({10,5,25,21},39/1),({10,5,24,22},24/1),({10,5,23,23},10/1),({9,6,36,10},1/1),({9,6,35,11},3/1),({9,6,34,12},6/1),({9,6,33,13},12/1),({9,6,32,14},20/1),({9,6,31,15},29/1),({9,6,30,16},41/1),({9,6,29,17},52/1),({9,6,28,18},60/1),({9,6,27,19},64/1),({9,6,26,20},61/1),({9,6,25,21},50/1),({9,6,24,22},33/1),({9,6,23,23},12/1),({8,7,37,9},1/1),({8,7,36,10},1/1),({8,7,35,11},3/1),({8,7,34,12},6/1),({8,7,33,13},11/1),({8,7,32,14},16/1),({8,7,31,15},25/1),({8,7,30,16},32/1),({8,7,29,17},41/1),({8,7,28,18},46/1),({8,7,27,19},50/1),({8,7,26,20},45/1),({8,7,25,21},40/1),({8,7,24,22},24/1),({8,7,23,23},9/1)}, (7,1) => {({15,2,31,21},1/1),({15,2,29,23},1/1),({15,2,27,25},1/1),({14,3,32,20},1/1),({14,3,31,21},2/1),({14,3,30,22},2/1),({14,3,29,23},2/1),({14,3,28,24},2/1),({14,3,27,25},2/1),({14,3,26,26},1/1),({13,4,36,16},1/1),({13,4,35,17},2/1),({13,4,34,18},3/1),({13,4,33,19},4/1),({13,4,32,20},5/1),({13,4,31,21},7/1),({13,4,30,22},7/1),({13,4,29,23},8/1),({13,4,28,24},6/1),({13,4,27,25},6/1),({13,4,26,26},1/1),({12,5,36,16},1/1),({12,5,35,17},3/1),({12,5,34,18},6/1),({12,5,33,19},9/1),({12,5,32,20},11/1),({12,5,31,21},13/1),({12,5,30,22},14/1),({12,5,29,23},14/1),({12,5,28,24},12/1),({12,5,27,25},8/1),({12,5,26,26},3/1),({11,6,39,13},1/1),({11,6,38,14},1/1),({11,6,37,15},3/1),({11,6,36,16},3/1),({11,6,35,17},7/1),({11,6,34,18},9/1),({11,6,33,19},15/1),({11,6,32,20},17/1),({11,6,31,21},22/1),({11,6,30,22},21/1),({11,6,29,23},23/1),({11,6,28,24},16/1),({11,6,27,25},13/1),({11,6,26,26},2/1),({10,7,38,14},1/1),({10,7,37,15},2/1),({10,7,36,16},4/1),({10,7,35,17},6/1),({10,7,34,18},9/1),({10,7,33,19},13/1),({10,7,32,20},17/1),({10,7,31,21},21/1),({10,7,30,22},22/1),({10,7,29,23},21/1),({10,7,28,24},17/1),({10,7,27,25},11/1),({10,7,26,26},4/1),({9,8,37,15},1/1),({9,8,36,16},2/1),({9,8,35,17},4/1),({9,8,34,18},5/1),({9,8,33,19},8/1),({9,8,32,20},10/1),({9,8,31,21},13/1),({9,8,30,22},13/1),({9,8,29,23},13/1),({9,8,28,24},10/1),({9,8,27,25},8/1),({9,8,26,26},2/1)}, (7,2) => {}, (9,0) => {({14,5,38,20},1/1),({14,5,36,22},1/1),({14,5,34,24},1/1),({14,5,32,26},1/1),({14,5,30,28},1/1),({13,6,37,21},1/1),({13,6,36,22},1/1),({13,6,35,23},1/1),({13,6,34,24},1/1),({13,6,33,25},1/1),({13,6,32,26},1/1),({13,6,31,27},1/1),({13,6,30,28},1/1),({12,7,36,22},1/1),({12,7,35,23},1/1),({12,7,34,24},2/1),({12,7,33,25},1/1),({12,7,32,26},2/1),({12,7,31,27},1/1),({12,7,30,28},1/1),({11,8,35,23},1/1),({11,8,34,24},1/1),({11,8,33,25},2/1),({11,8,32,26},2/1),({11,8,31,27},1/1),({11,8,30,28},1/1),({10,9,34,24},1/1),({10,9,33,25},1/1),({10,9,32,26},1/1),({10,9,31,27},1/1),({10,9,30,28},1/1)}, (9,1) => {({17,4,38,26},1/1),({17,4,36,28},1/1),({17,4,35,29},1/1),({17,4,34,30},1/1),({17,4,32,32},1/1),({16,5,41,23},1/1),({16,5,40,24},2/1),({16,5,39,25},3/1),({16,5,38,26},5/1),({16,5,37,27},7/1),({16,5,36,28},8/1),({16,5,35,29},8/1),({16,5,34,30},7/1),({16,5,33,31},5/1),({16,5,32,32},2/1),({15,6,43,21},1/1),({15,6,42,22},3/1),({15,6,41,23},6/1),({15,6,40,24},12/1),({15,6,39,25},17/1),({15,6,38,26},25/1),({15,6,37,27},30/1),({15,6,36,28},35/1),({15,6,35,29},32/1),({15,6,34,30},30/1),({15,6,33,31},18/1),({15,6,32,32},8/1),({14,7,45,19},1/1),({14,7,44,20},3/1),({14,7,43,21},7/1),({14,7,42,22},14/1),({14,7,41,23},25/1),({14,7,40,24},39/1),({14,7,39,25},56/1),({14,7,38,26},73/1),({14,7,37,27},85/1),({14,7,36,28},92/1),({14,7,35,29},89/1),({14,7,34,30},73/1),({14,7,33,31},49/1),({14,7,32,32},18/1),({13,8,46,18},1/1),({13,8,45,19},3/1),({13,8,44,20},9/1),({13,8,43,21},18/1),({13,8,42,22},34/1),({13,8,41,23},54/1),({13,8,40,24},82/1),({13,8,39,25},109/1),({13,8,38,26},139/1),({13,8,37,27},156/1),({13,8,36,28},167/1),({13,8,35,29},155/1),({13,8,34,30},132/1),({13,8,33,31},83/1),({13,8,32,32},32/1),({12,9,47,17},1/1),({12,9,46,18},3/1),({12,9,45,19},7/1),({12,9,44,20},15/1),({12,9,43,21},29/1),({12,9,42,22},49/1),({12,9,41,23},76/1),({12,9,40,24},109/1),({12,9,39,25},143/1),({12,9,38,26},174/1),({12,9,37,27},197/1),({12,9,36,28},203/1),({12,9,35,29},189/1),({12,9,34,30},156/1),({12,9,33,31},102/1),({12,9,32,32},35/1),({11,10,47,17},1/1),({11,10,46,18},2/1),({11,10,45,19},6/1),({11,10,44,20},13/1),({11,10,43,21},22/1),({11,10,42,22},38/1),({11,10,41,23},58/1),({11,10,40,24},80/1),({11,10,39,25},104/1),({11,10,38,26},127/1),({11,10,37,27},139/1),({11,10,36,28},145/1),({11,10,35,29},134/1),({11,10,34,30},109/1),({11,10,33,31},71/1),({11,10,32,32},27/1)}, (9,2) => {}, (11,0) => {}, (11,1) => {({19,6,41,35},1/1),({18,7,46,30},1/1),({18,7,45,31},2/1),({18,7,44,32},3/1),({18,7,43,33},4/1),({18,7,42,34},5/1),({18,7,41,35},6/1),({18,7,40,36},5/1),({18,7,39,37},3/1),({18,7,38,38},1/1),({17,8,49,27},1/1),({17,8,48,28},2/1),({17,8,47,29},5/1),({17,8,46,30},9/1),({17,8,45,31},16/1),({17,8,44,32},21/1),({17,8,43,33},28/1),({17,8,42,34},30/1),({17,8,41,35},32/1),({17,8,40,36},26/1),({17,8,39,37},19/1),({17,8,38,38},5/1),({16,9,50,26},2/1),({16,9,49,27},5/1),({16,9,48,28},12/1),({16,9,47,29},23/1),({16,9,46,30},37/1),({16,9,45,31},55/1),({16,9,44,32},74/1),({16,9,43,33},88/1),({16,9,42,34},96/1),({16,9,41,35},94/1),({16,9,40,36},78/1),({16,9,39,37},53/1),({16,9,38,38},19/1),({15,10,52,24},1/1),({15,10,51,25},3/1),({15,10,50,26},8/1),({15,10,49,27},18/1),({15,10,48,28},34/1),({15,10,47,29},58/1),({15,10,46,30},87/1),({15,10,45,31},122/1),({15,10,44,32},154/1),({15,10,43,33},181/1),({15,10,42,34},190/1),({15,10,41,35},183/1),({15,10,40,36},150/1),({15,10,39,37},102/1),({15,10,38,38},34/1),({14,11,52,24},2/1),({14,11,51,25},6/1),({14,11,50,26},14/1),({14,11,49,27},29/1),({14,11,48,28},52/1),({14,11,47,29},83/1),({14,11,46,30},123/1),({14,11,45,31},165/1),({14,11,44,32},204/1),({14,11,43,33},234/1),({14,11,42,34},244/1),({14,11,41,35},230/1),({14,11,40,36},190/1),({14,11,39,37},125/1),({14,11,38,38},43/1),({13,12,53,23},1/1),({13,12,52,24},2/1),({13,12,51,25},6/1),({13,12,50,26},13/1),({13,12,49,27},25/1),({13,12,48,28},42/1),({13,12,47,29},67/1),({13,12,46,30},94/1),({13,12,45,31},126/1),({13,12,44,32},153/1),({13,12,43,33},173/1),({13,12,42,34},178/1),({13,12,41,35},170/1),({13,12,40,36},136/1),({13,12,39,37},91/1),({13,12,38,38},31/1)}, (11,2) => {}, (13,0) => {}, (13,1) => {({20,9,47,41},1/1),({20,9,46,42},1/1),({19,10,52,36},1/1),({19,10,51,37},2/1),({19,10,50,38},4/1),({19,10,49,39},5/1),({19,10,48,40},7/1),({19,10,47,41},7/1),({19,10,46,42},7/1),({19,10,45,43},4/1),({19,10,44,44},2/1),({18,11,55,33},1/1),({18,11,54,34},2/1),({18,11,53,35},5/1),({18,11,52,36},9/1),({18,11,51,37},15/1),({18,11,50,38},22/1),({18,11,49,39},28/1),({18,11,48,40},32/1),({18,11,47,41},32/1),({18,11,46,42},28/1),({18,11,45,43},19/1),({18,11,44,44},7/1),({17,12,56,32},1/1),({17,12,55,33},3/1),({17,12,54,34},8/1),({17,12,53,35},16/1),({17,12,52,36},28/1),({17,12,51,37},41/1),({17,12,50,38},57/1),({17,12,49,39},69/1),({17,12,48,40},77/1),({17,12,47,41},74/1),({17,12,46,42},64/1),({17,12,45,43},42/1),({17,12,44,44},16/1),({16,13,57,31},1/1),({16,13,56,32},3/1),({16,13,55,33},8/1),({16,13,54,34},16/1),({16,13,53,35},29/1),({16,13,52,36},46/1),({16,13,51,37},66/1),({16,13,50,38},86/1),({16,13,49,39},102/1),({16,13,48,40},110/1),({16,13,47,41},106/1),({16,13,46,42},89/1),({16,13,45,43},59/1),({16,13,44,44},21/1),({15,14,57,31},1/1),({15,14,56,32},3/1),({15,14,55,33},7/1),({15,14,54,34},14/1),({15,14,53,35},24/1),({15,14,52,36},38/1),({15,14,51,37},53/1),({15,14,50,38},68/1),({15,14,49,39},79/1),({15,14,48,40},85/1),({15,14,47,41},81/1),({15,14,46,42},68/1),({15,14,45,43},44/1),({15,14,44,44},16/1)}, (15,0) => {}, (13,2) => {}, (15,1) => {({20,13,56,44},1/1),({20,13,55,45},1/1),({20,13,54,46},2/1),({20,13,53,47},3/1),({20,13,52,48},2/1),({20,13,51,49},2/1),({20,13,50,50},1/1),({19,14,59,41},1/1),({19,14,58,42},2/1),({19,14,57,43},4/1),({19,14,56,44},6/1),({19,14,55,45},9/1),({19,14,54,46},10/1),({19,14,53,47},11/1),({19,14,52,48},9/1),({19,14,51,49},7/1),({19,14,50,50},2/1),({18,15,60,40},1/1),({18,15,59,41},2/1),({18,15,58,42},5/1),({18,15,57,43},8/1),({18,15,56,44},12/1),({18,15,55,45},17/1),({18,15,54,46},19/1),({18,15,53,47},19/1),({18,15,52,48},17/1),({18,15,51,49},12/1),({18,15,50,50},4/1),({17,16,60,40},1/1),({17,16,59,41},3/1),({17,16,58,42},5/1),({17,16,57,43},8/1),({17,16,56,44},12/1),({17,16,55,45},15/1),({17,16,54,46},17/1),({17,16,53,47},17/1),({17,16,52,48},14/1),({17,16,51,49},10/1),({17,16,50,50},4/1)}, (17,0) => {}, (15,2) => {}, (17,1) => {({20,17,61,51},1/1),({20,17,60,52},1/1),({20,17,59,53},1/1),({20,17,58,54},1/1),({20,17,57,55},1/1),({20,17,56,56},1/1),({19,18,62,50},1/1),({19,18,61,51},1/1),({19,18,60,52},1/1),({19,18,59,53},1/1),({19,18,58,54},1/1),({19,18,57,55},1/1),({19,18,56,56},1/1)}, (17,2) => {}, (19,1) => {}, (0,0) => {({1,0,4,0},1/1)}, (2,0) => {({5,0,13,3},1/1),({5,0,12,4},1/1),({5,0,11,5},2/1),({5,0,10,6},1/1),({5,0,9,7},2/1),({4,1,15,1},1/1),({4,1,14,2},2/1),({4,1,13,3},3/1),({4,1,12,4},4/1),({4,1,11,5},4/1),({4,1,10,6},4/1),({4,1,9,7},3/1),({4,1,8,8},1/1),({3,2,15,1},1/1),({3,2,14,2},2/1),({3,2,13,3},3/1),({3,2,12,4},4/1),({3,2,11,5},4/1),({3,2,10,6},4/1),({3,2,9,7},3/1),({3,2,8,8},1/1)}, (2,1) => {}, (2,2) => {}, (4,0) => {({9,0,18,10},1/1),({9,0,16,12},1/1),({9,0,14,14},1/1),({8,1,22,6},1/1),({8,1,21,7},2/1),({8,1,20,8},3/1),({8,1,19,9},5/1),({8,1,18,10},6/1),({8,1,17,11},6/1),({8,1,16,12},6/1),({8,1,15,13},4/1),({8,1,14,14},1/1),({7,2,24,4},1/1),({7,2,23,5},2/1),({7,2,22,6},5/1),({7,2,21,7},8/1),({7,2,20,8},13/1),({7,2,19,9},16/1),({7,2,18,10},20/1),({7,2,17,11},19/1),({7,2,16,12},18/1),({7,2,15,13},11/1),({7,2,14,14},5/1),({6,3,24,4},2/1),({6,3,23,5},5/1),({6,3,22,6},9/1),({6,3,21,7},15/1),({6,3,20,8},22/1),({6,3,19,9},28/1),({6,3,18,10},32/1),({6,3,17,11},32/1),({6,3,16,12},27/1),({6,3,15,13},19/1),({6,3,14,14},7/1),({5,4,25,3},1/1),({5,4,24,4},2/1),({5,4,23,5},5/1),({5,4,22,6},9/1),({5,4,21,7},13/1),({5,4,20,8},19/1),({5,4,19,9},24/1),({5,4,18,10},27/1),({5,4,17,11},25/1),({5,4,16,12},23/1),({5,4,15,13},15/1),({5,4,14,14},6/1)}, (4,1) => {}, (4,2) => {}, (6,0) => {({11,2,29,11},1/1),({11,2,28,12},1/1),({11,2,27,13},3/1),({11,2,26,14},3/1),({11,2,25,15},5/1),({11,2,24,16},4/1),({11,2,23,17},6/1),({11,2,22,18},3/1),({11,2,21,19},4/1),({10,3,31,9},1/1),({10,3,30,10},2/1),({10,3,29,11},4/1),({10,3,28,12},8/1),({10,3,27,13},12/1),({10,3,26,14},16/1),({10,3,25,15},21/1),({10,3,24,16},23/1),({10,3,23,17},22/1),({10,3,22,18},20/1),({10,3,21,19},13/1),({10,3,20,20},4/1),({9,4,32,8},1/1),({9,4,31,9},3/1),({9,4,30,10},7/1),({9,4,29,11},13/1),({9,4,28,12},21/1),({9,4,27,13},32/1),({9,4,26,14},41/1),({9,4,25,15},51/1),({9,4,24,16},53/1),({9,4,23,17},54/1),({9,4,22,18},43/1),({9,4,21,19},31/1),({9,4,20,20},9/1),({8,5,33,7},1/1),({8,5,32,8},3/1),({8,5,31,9},6/1),({8,5,30,10},12/1),({8,5,29,11},22/1),({8,5,28,12},33/1),({8,5,27,13},47/1),({8,5,26,14},61/1),({8,5,25,15},71/1),({8,5,24,16},77/1),({8,5,23,17},74/1),({8,5,22,18},61/1),({8,5,21,19},41/1),({8,5,20,20},15/1),({7,6,33,7},1/1),({7,6,32,8},2/1),({7,6,31,9},6/1),({7,6,30,10},10/1),({7,6,29,11},18/1),({7,6,28,12},27/1),({7,6,27,13},38/1),({7,6,26,14},46/1),({7,6,25,15},57/1),({7,6,24,16},58/1),({7,6,23,17},56/1),({7,6,22,18},46/1),({7,6,21,19},32/1),({7,6,20,20},9/1)}, (6,1) => {({14,1,26,20},1/1),({13,2,26,20},1/1),({13,2,25,21},1/1),({12,3,31,15},1/1),({12,3,30,16},1/1),({12,3,29,17},1/1),({12,3,28,18},1/1),({12,3,27,19},1/1),({12,3,26,20},1/1),({12,3,25,21},1/1),({12,3,24,22},1/1),({11,4,30,16},1/1),({11,4,29,17},1/1),({11,4,28,18},1/1),({11,4,27,19},1/1),({11,4,26,20},1/1),({11,4,25,21},1/1),({11,4,24,22},1/1),({11,4,23,23},1/1),({10,5,29,17},1/1),({10,5,28,18},1/1),({10,5,27,19},1/1),({10,5,26,20},1/1),({10,5,25,21},1/1),({10,5,24,22},1/1),({9,6,28,18},1/1),({9,6,27,19},1/1),({9,6,26,20},1/1),({9,6,25,21},1/1),({8,7,27,19},1/1),({8,7,26,20},1/1)}, (6,2) => {}, (8,0) => {({13,4,36,16},1/1),({13,4,35,17},1/1),({13,4,34,18},2/1),({13,4,33,19},2/1),({13,4,32,20},3/1),({13,4,31,21},2/1),({13,4,30,22},3/1),({13,4,29,23},2/1),({13,4,28,24},3/1),({13,4,27,25},1/1),({13,4,26,26},1/1),({12,5,36,16},1/1),({12,5,35,17},2/1),({12,5,34,18},4/1),({12,5,33,19},6/1),({12,5,32,20},7/1),({12,5,31,21},8/1),({12,5,30,22},8/1),({12,5,29,23},8/1),({12,5,28,24},7/1),({12,5,27,25},5/1),({12,5,26,26},2/1),({11,6,39,13},1/1),({11,6,38,14},1/1),({11,6,37,15},2/1),({11,6,36,16},3/1),({11,6,35,17},5/1),({11,6,34,18},8/1),({11,6,33,19},11/1),({11,6,32,20},15/1),({11,6,31,21},17/1),({11,6,30,22},19/1),({11,6,29,23},17/1),({11,6,28,24},15/1),({11,6,27,25},9/1),({11,6,26,26},4/1),({10,7,38,14},1/1),({10,7,37,15},2/1),({10,7,36,16},4/1),({10,7,35,17},6/1),({10,7,34,18},9/1),({10,7,33,19},13/1),({10,7,32,20},17/1),({10,7,31,21},21/1),({10,7,30,22},23/1),({10,7,29,23},22/1),({10,7,28,24},18/1),({10,7,27,25},12/1),({10,7,26,26},4/1),({9,8,37,15},1/1),({9,8,36,16},2/1),({9,8,35,17},4/1),({9,8,34,18},6/1),({9,8,33,19},9/1),({9,8,32,20},12/1),({9,8,31,21},14/1),({9,8,30,22},16/1),({9,8,29,23},15/1),({9,8,28,24},13/1),({9,8,27,25},8/1),({9,8,26,26},3/1)}, (8,1) => {({16,3,35,23},1/1),({16,3,33,25},1/1),({16,3,32,26},1/1),({16,3,31,27},1/1),({16,3,29,29},1/1),({15,4,37,21},1/1),({15,4,36,22},2/1),({15,4,35,23},3/1),({15,4,34,24},4/1),({15,4,33,25},5/1),({15,4,32,26},6/1),({15,4,31,27},5/1),({15,4,30,28},3/1),({15,4,29,29},1/1),({14,5,40,18},1/1),({14,5,39,19},2/1),({14,5,38,20},4/1),({14,5,37,21},7/1),({14,5,36,22},11/1),({14,5,35,23},15/1),({14,5,34,24},18/1),({14,5,33,25},21/1),({14,5,32,26},20/1),({14,5,31,27},19/1),({14,5,30,28},11/1),({14,5,29,29},5/1),({13,6,41,17},1/1),({13,6,40,18},3/1),({13,6,39,19},7/1),({13,6,38,20},13/1),({13,6,37,21},21/1),({13,6,36,22},30/1),({13,6,35,23},40/1),({13,6,34,24},46/1),({13,6,33,25},50/1),({13,6,32,26},49/1),({13,6,31,27},40/1),({13,6,30,28},28/1),({13,6,29,29},10/1),({12,7,43,15},1/1),({12,7,42,16},2/1),({12,7,41,17},5/1),({12,7,40,18},9/1),({12,7,39,19},18/1),({12,7,38,20},28/1),({12,7,37,21},43/1),({12,7,36,22},57/1),({12,7,35,23},74/1),({12,7,34,24},83/1),({12,7,33,25},89/1),({12,7,32,26},81/1),({12,7,31,27},70/1),({12,7,30,28},43/1),({12,7,29,29},18/1),({11,8,43,15},1/1),({11,8,42,16},3/1),({11,8,41,17},7/1),({11,8,40,18},14/1),({11,8,39,19},24/1),({11,8,38,20},37/1),({11,8,37,21},54/1),({11,8,36,22},71/1),({11,8,35,23},87/1),({11,8,34,24},99/1),({11,8,33,25},103/1),({11,8,32,26},96/1),({11,8,31,27},79/1),({11,8,30,28},52/1),({11,8,29,29},17/1),({10,9,44,14},1/1),({10,9,43,15},1/1),({10,9,42,16},3/1),({10,9,41,17},7/1),({10,9,40,18},11/1),({10,9,39,19},19/1),({10,9,38,20},29/1),({10,9,37,21},40/1),({10,9,36,22},52/1),({10,9,35,23},64/1),({10,9,34,24},69/1),({10,9,33,25},73/1),({10,9,32,26},68/1),({10,9,31,27},55/1),({10,9,30,28},35/1),({10,9,29,29},14/1)}, (8,2) => {}, (10,0) => {}, (10,1) => {({18,5,40,30},1/1),({18,5,38,32},1/1),({18,5,36,34},1/1),({17,6,44,26},1/1),({17,6,43,27},2/1),({17,6,42,28},3/1),({17,6,41,29},5/1),({17,6,40,30},7/1),({17,6,39,31},8/1),({17,6,38,32},8/1),({17,6,37,33},7/1),({17,6,36,34},5/1),({17,6,35,35},2/1),({16,7,46,24},1/1),({16,7,45,25},3/1),({16,7,44,26},7/1),({16,7,43,27},12/1),({16,7,42,28},20/1),({16,7,41,29},27/1),({16,7,40,30},35/1),({16,7,39,31},37/1),({16,7,38,32},39/1),({16,7,37,33},31/1),({16,7,36,34},23/1),({16,7,35,35},6/1),({15,8,48,22},1/1),({15,8,47,23},3/1),({15,8,46,24},8/1),({15,8,45,25},16/1),({15,8,44,26},29/1),({15,8,43,27},46/1),({15,8,42,28},66/1),({15,8,41,29},86/1),({15,8,40,30},102/1),({15,8,39,31},110/1),({15,8,38,32},106/1),({15,8,37,33},89/1),({15,8,36,34},59/1),({15,8,35,35},21/1),({14,9,49,21},1/1),({14,9,48,22},4/1),({14,9,47,23},10/1),({14,9,46,24},22/1),({14,9,45,25},39/1),({14,9,44,26},67/1),({14,9,43,27},98/1),({14,9,42,28},136/1),({14,9,41,29},168/1),({14,9,40,30},197/1),({14,9,39,31},204/1),({14,9,38,32},197/1),({14,9,37,33},159/1),({14,9,36,34},109/1),({14,9,35,35},35/1),({13,10,50,20},1/1),({13,10,49,21},3/1),({13,10,48,22},8/1),({13,10,47,23},18/1),({13,10,46,24},35/1),({13,10,45,25},60/1),({13,10,44,26},94/1),({13,10,43,27},135/1),({13,10,42,28},179/1),({13,10,41,29},219/1),({13,10,40,30},248/1),({13,10,39,31},257/1),({13,10,38,32},241/1),({13,10,37,33},198/1),({13,10,36,34},130/1),({13,10,35,35},46/1),({12,11,50,20},1/1),({12,11,49,21},3/1),({12,11,48,22},7/1),({12,11,47,23},15/1),({12,11,46,24},28/1),({12,11,45,25},47/1),({12,11,44,26},72/1),({12,11,43,27},101/1),({12,11,42,28},133/1),({12,11,41,29},160/1),({12,11,40,30},180/1),({12,11,39,31},184/1),({12,11,38,32},173/1),({12,11,37,33},140/1),({12,11,36,34},93/1),({12,11,35,35},31/1)}, (10,2) => {}, (12,0) => {}, (12,1) => {({20,7,41,41},1/1),({19,8,47,35},1/1),({19,8,46,36},2/1),({19,8,45,37},2/1),({19,8,44,38},2/1),({19,8,43,39},2/1),({19,8,42,40},2/1),({19,8,41,41},1/1),({18,9,51,31},1/1),({18,9,50,32},2/1),({18,9,49,33},5/1),({18,9,48,34},7/1),({18,9,47,35},12/1),({18,9,46,36},15/1),({18,9,45,37},19/1),({18,9,44,38},17/1),({18,9,43,39},17/1),({18,9,42,40},10/1),({18,9,41,41},5/1),({17,10,53,29},1/1),({17,10,52,30},2/1),({17,10,51,31},6/1),({17,10,50,32},13/1),({17,10,49,33},22/1),({17,10,48,34},34/1),({17,10,47,35},47/1),({17,10,46,36},57/1),({17,10,45,37},64/1),({17,10,44,38},64/1),({17,10,43,39},53/1),({17,10,42,40},36/1),({17,10,41,41},13/1),({16,11,54,28},1/1),({16,11,53,29},4/1),({16,11,52,30},10/1),({16,11,51,31},20/1),({16,11,50,32},36/1),({16,11,49,33},58/1),({16,11,48,34},82/1),({16,11,47,35},108/1),({16,11,46,36},126/1),({16,11,45,37},138/1),({16,11,44,38},131/1),({16,11,43,39},112/1),({16,11,42,40},72/1),({16,11,41,41},27/1),({15,12,55,27},1/1),({15,12,54,28},3/1),({15,12,53,29},8/1),({15,12,52,30},18/1),({15,12,51,31},34/1),({15,12,50,32},56/1),({15,12,49,33},86/1),({15,12,48,34},118/1),({15,12,47,35},149/1),({15,12,46,36},174/1),({15,12,45,37},183/1),({15,12,44,38},174/1),({15,12,43,39},145/1),({15,12,42,40},96/1),({15,12,41,41},33/1),({14,13,55,27},1/1),({14,13,54,28},3/1),({14,13,53,29},8/1),({14,13,52,30},15/1),({14,13,51,31},29/1),({14,13,50,32},46/1),({14,13,49,33},68/1),({14,13,48,34},91/1),({14,13,47,35},116/1),({14,13,46,36},130/1),({14,13,45,37},138/1),({14,13,44,38},129/1),({14,13,43,39},107/1),({14,13,42,40},70/1),({14,13,41,41},27/1)}, (12,2) => {}, (14,0) => {}, (14,1) => {({20,11,52,42},1/1),({20,11,51,43},1/1),({20,11,50,44},2/1),({20,11,49,45},1/1),({20,11,48,46},2/1),({19,12,56,38},1/1),({19,12,55,39},2/1),({19,12,54,40},4/1),({19,12,53,41},7/1),({19,12,52,42},9/1),({19,12,51,43},11/1),({19,12,50,44},12/1),({19,12,49,45},10/1),({19,12,48,46},7/1),({19,12,47,47},3/1),({18,13,58,36},1/1),({18,13,57,37},2/1),({18,13,56,38},5/1),({18,13,55,39},9/1),({18,13,54,40},16/1),({18,13,53,41},22/1),({18,13,52,42},29/1),({18,13,51,43},32/1),({18,13,50,44},34/1),({18,13,49,45},28/1),({18,13,48,46},20/1),({18,13,47,47},6/1),({17,14,58,36},2/1),({17,14,57,37},5/1),({17,14,56,38},10/1),({17,14,55,39},18/1),({17,14,54,40},28/1),({17,14,53,41},38/1),({17,14,52,42},47/1),({17,14,51,43},52/1),({17,14,50,44},51/1),({17,14,49,45},44/1),({17,14,48,46},29/1),({17,14,47,47},10/1),({16,15,59,35},1/1),({16,15,58,36},2/1),({16,15,57,37},5/1),({16,15,56,38},10/1),({16,15,55,39},16/1),({16,15,54,40},24/1),({16,15,53,41},32/1),({16,15,52,42},39/1),({16,15,51,43},42/1),({16,15,50,44},42/1),({16,15,49,45},34/1),({16,15,48,46},24/1),({16,15,47,47},8/1)}, (16,0) => {}, (14,2) => {}, (16,1) => {({20,15,59,47},1/1),({20,15,58,48},1/1),({20,15,57,49},2/1),({20,15,56,50},2/1),({20,15,55,51},3/1),({20,15,54,52},1/1),({20,15,53,53},1/1),({19,16,61,45},1/1),({19,16,60,46},2/1),({19,16,59,47},3/1),({19,16,58,48},4/1),({19,16,57,49},5/1),({19,16,56,50},6/1),({19,16,55,51},5/1),({19,16,54,52},3/1),({19,16,53,53},1/1),({18,17,61,45},1/1),({18,17,60,46},2/1),({18,17,59,47},3/1),({18,17,58,48},4/1),({18,17,57,49},5/1),({18,17,56,50},6/1),({18,17,55,51},5/1),({18,17,54,52},3/1),({18,17,53,53},1/1)}, (18,0) => {}, (16,2) => {}, (18,1) => {({20,19,62,56},1/1)}, (18,2) => {}, (20,1) => {}, (1,0) => {({3,0,9,1},1/1),({3,0,8,2},1/1),({3,0,7,3},1/1),({3,0,6,4},1/1),({2,1,10,0},1/1),({2,1,9,1},1/1),({2,1,8,2},1/1),({2,1,7,3},1/1),({2,1,6,4},1/1)}, (1,1) => {}, (3,0) => {({7,0,16,6},1/1),({7,0,15,7},1/1),({7,0,14,8},1/1),({7,0,13,9},2/1),({7,0,12,10},1/1),({6,1,19,3},1/1),({6,1,18,4},2/1),({6,1,17,5},4/1),({6,1,16,6},5/1),({6,1,15,7},7/1),({6,1,14,8},7/1),({6,1,13,9},7/1),({6,1,12,10},4/1),({6,1,11,11},2/1),({5,2,20,2},1/1),({5,2,19,3},2/1),({5,2,18,4},5/1),({5,2,17,5},8/1),({5,2,16,6},10/1),({5,2,15,7},13/1),({5,2,14,8},14/1),({5,2,13,9},12/1),({5,2,12,10},8/1),({5,2,11,11},3/1),({4,3,20,2},1/1),({4,3,19,3},3/1),({4,3,18,4},5/1),({4,3,17,5},7/1),({4,3,16,6},10/1),({4,3,15,7},12/1),({4,3,14,8},12/1),({4,3,13,9},11/1),({4,3,12,10},7/1),({4,3,11,11},2/1)}, (3,1) => {}, (3,2) => {}, (5,0) => {({10,1,24,10},1/1),({10,1,23,11},1/1),({10,1,22,12},2/1),({10,1,21,13},2/1),({10,1,20,14},3/1),({10,1,19,15},2/1),({10,1,18,16},2/1),({9,2,27,7},1/1),({9,2,26,8},2/1),({9,2,25,9},4/1),({9,2,24,10},7/1),({9,2,23,11},10/1),({9,2,22,12},13/1),({9,2,21,13},15/1),({9,2,20,14},15/1),({9,2,19,15},13/1),({9,2,18,16},9/1),({9,2,17,17},3/1),({8,3,28,6},1/1),({8,3,27,7},3/1),({8,3,26,8},7/1),({8,3,25,9},12/1),({8,3,24,10},20/1),({8,3,23,11},27/1),({8,3,22,12},35/1),({8,3,21,13},38/1),({8,3,20,14},39/1),({8,3,19,15},32/1),({8,3,18,16},23/1),({8,3,17,17},7/1),({7,4,29,5},1/1),({7,4,28,6},3/1),({7,4,27,7},7/1),({7,4,26,8},13/1),({7,4,25,9},22/1),({7,4,24,10},33/1),({7,4,23,11},44/1),({7,4,22,12},53/1),({7,4,21,13},58/1),({7,4,20,14},57/1),({7,4,19,15},48/1),({7,4,18,16},32/1),({7,4,17,17},11/1),({6,5,29,5},1/1),({6,5,28,6},3/1),({6,5,27,7},6/1),({6,5,26,8},11/1),({6,5,25,9},18/1),({6,5,24,10},27/1),({6,5,23,11},35/1),({6,5,22,12},42/1),({6,5,21,13},45/1),({6,5,20,14},45/1),({6,5,19,15},37/1),({6,5,18,16},25/1),({6,5,17,17},8/1)}, (5,1) => {({13,0,20,20},1/1)}};
end;