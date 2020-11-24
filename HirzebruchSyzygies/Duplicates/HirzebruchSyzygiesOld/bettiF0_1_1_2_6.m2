A := QQ[t_0,t_1,t_2,t_3];
--tb stands for Total Betti numbers
tb1126 = new HashTable from {(5,2) => 0, (7,0) => 0, (7,1) => 461448/1, (7,2) => 0, (9,0) => 0, (9,1) => 797368/1, (9,2) => 0, (11,0) => 0, (11,1) => 562224/1, (11,2) => 0, (13,0) => 0, (13,1) => 159120/1, (15,0) => 0, (13,2) => 0, (15,1) => 15708/1, (17,0) => 0, (15,2) => 0, (17,1) => 356/1, (17,2) => 0, (19,1) => 0, (0,0) => 4/1, (2,0) => 252/1, (2,1) => 160/1, (2,2) => 0, (4,0) => 0, (4,1) => 26928/1, (4,2) => 0, (6,0) => 0, (6,1) => 243984/1, (6,2) => 0, (8,0) => 0, (8,1) => 680680/1, (8,2) => 0, (10,0) => 0, (10,1) => 747864/1, (10,2) => 0, (12,0) => 0, (12,1) => 337008/1, (12,2) => 0, (14,0) => 0, (14,1) => 57936/1, (16,0) => 0, (14,2) => 0, (16,1) => 2988/1, (18,0) => 0, (16,2) => 0, (18,1) => 20/1, (18,2) => 0, (20,1) => 0, (1,0) => 52/1, (1,1) => 0, (3,0) => 364/1, (3,1) => 4080/1, (3,2) => 0, (5,0) => 0, (5,1) => 97104/1};
--sb represents the betti numbers as sums of Schur functors
sb1126 = new HashTable from {(5,2) => {}, (7,0) => {}, (7,1) => {({15,2,31,18},1/1),({15,2,30,19},1/1),({15,2,29,20},2/1),({15,2,28,21},2/1),({15,2,27,22},2/1),({15,2,26,23},2/1),({15,2,25,24},1/1),({14,3,34,15},1/1),({14,3,33,16},2/1),({14,3,32,17},5/1),({14,3,31,18},8/1),({14,3,30,19},12/1),({14,3,29,20},15/1),({14,3,28,21},17/1),({14,3,27,22},16/1),({14,3,26,23},13/1),({14,3,25,24},7/1),({13,4,36,13},1/1),({13,4,35,14},4/1),({13,4,34,15},9/1),({13,4,33,16},17/1),({13,4,32,17},28/1),({13,4,31,18},42/1),({13,4,30,19},54/1),({13,4,29,20},64/1),({13,4,28,21},67/1),({13,4,27,22},63/1),({13,4,26,23},48/1),({13,4,25,24},26/1),({12,5,38,11},1/1),({12,5,37,12},3/1),({12,5,36,13},9/1),({12,5,35,14},19/1),({12,5,34,15},37/1),({12,5,33,16},60/1),({12,5,32,17},91/1),({12,5,31,18},122/1),({12,5,30,19},152/1),({12,5,29,20},170/1),({12,5,28,21},174/1),({12,5,27,22},156/1),({12,5,26,23},119/1),({12,5,25,24},64/1),({11,6,39,10},1/1),({11,6,38,11},4/1),({11,6,37,12},12/1),({11,6,36,13},25/1),({11,6,35,14},49/1),({11,6,34,15},83/1),({11,6,33,16},129/1),({11,6,32,17},181/1),({11,6,31,18},235/1),({11,6,30,19},280/1),({11,6,29,20},308/1),({11,6,28,21},306/1),({11,6,27,22},272/1),({11,6,26,23},204/1),({11,6,25,24},110/1),({10,7,40,9},1/1),({10,7,39,10},3/1),({10,7,38,11},9/1),({10,7,37,12},20/1),({10,7,36,13},41/1),({10,7,35,14},72/1),({10,7,34,15},117/1),({10,7,33,16},172/1),({10,7,32,17},236/1),({10,7,31,18},297/1),({10,7,30,19},348/1),({10,7,29,20},374/1),({10,7,28,21},369/1),({10,7,27,22},323/1),({10,7,26,23},242/1),({10,7,25,24},129/1),({9,8,40,9},1/1),({9,8,39,10},3/1),({9,8,38,11},8/1),({9,8,37,12},17/1),({9,8,36,13},33/1),({9,8,35,14},57/1),({9,8,34,15},89/1),({9,8,33,16},129/1),({9,8,32,17},173/1),({9,8,31,18},216/1),({9,8,30,19},249/1),({9,8,29,20},266/1),({9,8,28,21},260/1),({9,8,27,22},228/1),({9,8,26,23},169/1),({9,8,25,24},90/1)}, (7,2) => {}, (9,0) => {}, (9,1) => {({17,4,38,23},1/1),({17,4,37,24},1/1),({17,4,36,25},2/1),({17,4,35,26},3/1),({17,4,34,27},3/1),({17,4,33,28},3/1),({17,4,32,29},3/1),({17,4,31,30},1/1),({16,5,41,20},1/1),({16,5,40,21},2/1),({16,5,39,22},5/1),({16,5,38,23},9/1),({16,5,37,24},15/1),({16,5,36,25},19/1),({16,5,35,26},24/1),({16,5,34,27},26/1),({16,5,33,28},25/1),({16,5,32,29},19/1),({16,5,31,30},10/1),({15,6,43,18},1/1),({15,6,42,19},4/1),({15,6,41,20},9/1),({15,6,40,21},19/1),({15,6,39,22},32/1),({15,6,38,23},51/1),({15,6,37,24},71/1),({15,6,36,25},90/1),({15,6,35,26},102/1),({15,6,34,27},107/1),({15,6,33,28},97/1),({15,6,32,29},74/1),({15,6,31,30},40/1),({14,7,45,16},1/1),({14,7,44,17},3/1),({14,7,43,18},9/1),({14,7,42,19},20/1),({14,7,41,20},41/1),({14,7,40,21},70/1),({14,7,39,22},111/1),({14,7,38,23},158/1),({14,7,37,24},208/1),({14,7,36,25},250/1),({14,7,35,26},278/1),({14,7,34,27},276/1),({14,7,33,28},248/1),({14,7,32,29},187/1),({14,7,31,30},100/1),({13,8,46,15},1/1),({13,8,45,16},4/1),({13,8,44,17},12/1),({13,8,43,18},27/1),({13,8,42,19},54/1),({13,8,41,20},96/1),({13,8,40,21},156/1),({13,8,39,22},230/1),({13,8,38,23},315/1),({13,8,37,24},397/1),({13,8,36,25},466/1),({13,8,35,26},502/1),({13,8,34,27},495/1),({13,8,33,28},434/1),({13,8,32,29},325/1),({13,8,31,30},174/1),({12,9,47,14},1/1),({12,9,46,15},3/1),({12,9,45,16},9/1),({12,9,44,17},21/1),({12,9,43,18},45/1),({12,9,42,19},82/1),({12,9,41,20},139/1),({12,9,40,21},214/1),({12,9,39,22},307/1),({12,9,38,23},406/1),({12,9,37,24},504/1),({12,9,36,25},577/1),({12,9,35,26},614/1),({12,9,34,27},598/1),({12,9,33,28},522/1),({12,9,32,29},386/1),({12,9,31,30},207/1),({11,10,47,14},1/1),({11,10,46,15},3/1),({11,10,45,16},8/1),({11,10,44,17},19/1),({11,10,43,18},36/1),({11,10,42,19},66/1),({11,10,41,20},108/1),({11,10,40,21},162/1),({11,10,39,22},228/1),({11,10,38,23},300/1),({11,10,37,24},364/1),({11,10,36,25},415/1),({11,10,35,26},438/1),({11,10,34,27},423/1),({11,10,33,28},368/1),({11,10,32,29},273/1),({11,10,31,30},143/1)}, (9,2) => {}, (11,0) => {}, (11,1) => {({19,6,41,32},1/1),({19,6,40,33},1/1),({19,6,39,34},1/1),({19,6,38,35},1/1),({18,7,46,27},1/1),({18,7,45,28},2/1),({18,7,44,29},4/1),({18,7,43,30},6/1),({18,7,42,31},9/1),({18,7,41,32},11/1),({18,7,40,33},13/1),({18,7,39,34},12/1),({18,7,38,35},10/1),({18,7,37,36},5/1),({17,8,49,24},1/1),({17,8,48,25},2/1),({17,8,47,26},5/1),({17,8,46,27},10/1),({17,8,45,28},19/1),({17,8,44,29},29/1),({17,8,43,30},42/1),({17,8,42,31},53/1),({17,8,41,32},62/1),({17,8,40,33},64/1),({17,8,39,34},59/1),({17,8,38,35},45/1),({17,8,37,36},25/1),({16,9,50,23},2/1),({16,9,49,24},5/1),({16,9,48,25},13/1),({16,9,47,26},26/1),({16,9,46,27},46/1),({16,9,45,28},72/1),({16,9,44,29},106/1),({16,9,43,30},138/1),({16,9,42,31},169/1),({16,9,41,32},187/1),({16,9,40,33},188/1),({16,9,39,34},168/1),({16,9,38,35},128/1),({16,9,37,36},68/1),({15,10,52,21},1/1),({15,10,51,22},3/1),({15,10,50,23},8/1),({15,10,49,24},19/1),({15,10,48,25},38/1),({15,10,47,26},69/1),({15,10,46,27},111/1),({15,10,45,28},166/1),({15,10,44,29},227/1),({15,10,43,30},289/1),({15,10,42,31},338/1),({15,10,41,32},366/1),({15,10,40,33},360/1),({15,10,39,34},318/1),({15,10,38,35},237/1),({15,10,37,36},127/1),({14,11,52,21},2/1),({14,11,51,22},6/1),({14,11,50,23},15/1),({14,11,49,24},32/1),({14,11,48,25},61/1),({14,11,47,26},102/1),({14,11,46,27},161/1),({14,11,45,28},230/1),({14,11,44,29},307/1),({14,11,43,30},380/1),({14,11,42,31},438/1),({14,11,41,32},465/1),({14,11,40,33},455/1),({14,11,39,34},396/1),({14,11,38,35},294/1),({14,11,37,36},157/1),({13,12,53,20},1/1),({13,12,52,21},2/1),({13,12,51,22},6/1),({13,12,50,23},14/1),({13,12,49,24},28/1),({13,12,48,25},50/1),({13,12,47,26},83/1),({13,12,46,27},125/1),({13,12,45,28},177/1),({13,12,44,29},232/1),({13,12,43,30},283/1),({13,12,42,31},322/1),({13,12,41,32},342/1),({13,12,40,33},329/1),({13,12,39,34},286/1),({13,12,38,35},212/1),({13,12,37,36},112/1)}, (11,2) => {}, (13,0) => {}, (13,1) => {({20,9,47,38},1/1),({20,9,46,39},1/1),({20,9,45,40},1/1),({20,9,44,41},1/1),({20,9,43,42},1/1),({19,10,52,33},1/1),({19,10,51,34},2/1),({19,10,50,35},4/1),({19,10,49,36},6/1),({19,10,48,37},9/1),({19,10,47,38},11/1),({19,10,46,39},12/1),({19,10,45,40},11/1),({19,10,44,41},9/1),({19,10,43,42},5/1),({18,11,55,30},1/1),({18,11,54,31},2/1),({18,11,53,32},5/1),({18,11,52,33},9/1),({18,11,51,34},16/1),({18,11,50,35},24/1),({18,11,49,36},34/1),({18,11,48,37},42/1),({18,11,47,38},49/1),({18,11,46,39},50/1),({18,11,45,40},46/1),({18,11,44,41},35/1),({18,11,43,42},19/1),({17,12,56,29},1/1),({17,12,55,30},3/1),({17,12,54,31},8/1),({17,12,53,32},16/1),({17,12,52,33},29/1),({17,12,51,34},45/1),({17,12,50,35},66/1),({17,12,49,36},86/1),({17,12,48,37},104/1),({17,12,47,38},114/1),({17,12,46,39},115/1),({17,12,45,40},102/1),({17,12,44,41},77/1),({17,12,43,42},41/1),({16,13,57,28},1/1),({16,13,56,29},3/1),({16,13,55,30},8/1),({16,13,54,31},16/1),({16,13,53,32},30/1),({16,13,52,33},49/1),({16,13,51,34},74/1),({16,13,50,35},101/1),({16,13,49,36},129/1),({16,13,48,37},151/1),({16,13,47,38},164/1),({16,13,46,39},161/1),({16,13,45,40},142/1),({16,13,44,41},106/1),({16,13,43,42},57/1),({15,14,57,28},1/1),({15,14,56,29},3/1),({15,14,55,30},7/1),({15,14,54,31},14/1),({15,14,53,32},25/1),({15,14,52,33},41/1),({15,14,51,34},60/1),({15,14,50,35},81/1),({15,14,49,36},101/1),({15,14,48,37},118/1),({15,14,47,38},126/1),({15,14,46,39},123/1),({15,14,45,40},107/1),({15,14,44,41},80/1),({15,14,43,42},43/1)}, (15,0) => {}, (13,2) => {}, (15,1) => {({20,13,56,41},1/1),({20,13,55,42},1/1),({20,13,54,43},2/1),({20,13,53,44},3/1),({20,13,52,45},3/1),({20,13,51,46},3/1),({20,13,50,47},3/1),({20,13,49,48},1/1),({19,14,59,38},1/1),({19,14,58,39},2/1),({19,14,57,40},4/1),({19,14,56,41},6/1),({19,14,55,42},9/1),({19,14,54,43},11/1),({19,14,53,44},13/1),({19,14,52,45},13/1),({19,14,51,46},12/1),({19,14,50,47},9/1),({19,14,49,48},5/1),({18,15,60,37},1/1),({18,15,59,38},2/1),({18,15,58,39},5/1),({18,15,57,40},8/1),({18,15,56,41},12/1),({18,15,55,42},17/1),({18,15,54,43},21/1),({18,15,53,44},23/1),({18,15,52,45},24/1),({18,15,51,46},21/1),({18,15,50,47},16/1),({18,15,49,48},9/1),({17,16,60,37},1/1),({17,16,59,38},3/1),({17,16,58,39},5/1),({17,16,57,40},8/1),({17,16,56,41},12/1),({17,16,55,42},16/1),({17,16,54,43},19/1),({17,16,53,44},21/1),({17,16,52,45},20/1),({17,16,51,46},18/1),({17,16,50,47},14/1),({17,16,49,48},7/1)}, (17,0) => {}, (15,2) => {}, (17,1) => {({20,17,61,48},1/1),({20,17,60,49},1/1),({20,17,59,50},1/1),({20,17,58,51},1/1),({20,17,57,52},1/1),({20,17,56,53},1/1),({19,18,62,47},1/1),({19,18,61,48},1/1),({19,18,60,49},1/1),({19,18,59,50},1/1),({19,18,58,51},1/1),({19,18,57,52},1/1),({19,18,56,53},1/1)}, (17,2) => {}, (19,1) => {}, (0,0) => {({1,0,1,0},1/1)}, (2,0) => {({4,1,12,1},1/1),({4,1,11,2},1/1),({4,1,10,3},1/1),({4,1,9,4},1/1),({4,1,8,5},1/1),({4,1,7,6},1/1),({3,2,12,1},1/1),({3,2,11,2},1/1),({3,2,10,3},1/1),({3,2,9,4},1/1),({3,2,8,5},1/1),({3,2,7,6},1/1)}, (2,1) => {({7,0,14,5},1/1),({7,0,12,7},1/1),({7,0,11,8},1/1)}, (2,2) => {}, (4,0) => {}, (4,1) => {({11,0,19,12},1/1),({11,0,18,13},1/1),({11,0,17,14},1/1),({10,1,23,8},1/1),({10,1,22,9},2/1),({10,1,21,10},3/1),({10,1,20,11},4/1),({10,1,19,12},5/1),({10,1,18,13},5/1),({10,1,17,14},4/1),({10,1,16,15},2/1),({9,2,25,6},1/1),({9,2,24,7},1/1),({9,2,23,8},4/1),({9,2,22,9},6/1),({9,2,21,10},10/1),({9,2,20,11},12/1),({9,2,19,12},14/1),({9,2,18,13},13/1),({9,2,17,14},11/1),({9,2,16,15},6/1),({8,3,26,5},1/1),({8,3,25,6},3/1),({8,3,24,7},6/1),({8,3,23,8},10/1),({8,3,22,9},15/1),({8,3,21,10},20/1),({8,3,20,11},24/1),({8,3,19,12},25/1),({8,3,18,13},23/1),({8,3,17,14},18/1),({8,3,16,15},10/1),({7,4,27,4},1/1),({7,4,26,5},2/1),({7,4,25,6},5/1),({7,4,24,7},8/1),({7,4,23,8},14/1),({7,4,22,9},19/1),({7,4,21,10},25/1),({7,4,20,11},28/1),({7,4,19,12},30/1),({7,4,18,13},27/1),({7,4,17,14},21/1),({7,4,16,15},11/1),({6,5,27,4},1/1),({6,5,26,5},2/1),({6,5,25,6},4/1),({6,5,24,7},7/1),({6,5,23,8},11/1),({6,5,22,9},15/1),({6,5,21,10},18/1),({6,5,20,11},20/1),({6,5,19,12},21/1),({6,5,18,13},19/1),({6,5,17,14},14/1),({6,5,16,15},7/1)}, (4,2) => {}, (6,0) => {}, (6,1) => {({14,1,26,17},1/1),({14,1,25,18},1/1),({14,1,24,19},1/1),({14,1,23,20},1/1),({13,2,30,13},1/1),({13,2,29,14},2/1),({13,2,28,15},4/1),({13,2,27,16},5/1),({13,2,26,17},8/1),({13,2,25,18},9/1),({13,2,24,19},9/1),({13,2,23,20},7/1),({13,2,22,21},4/1),({12,3,32,11},1/1),({12,3,31,12},3/1),({12,3,30,13},7/1),({12,3,29,14},13/1),({12,3,28,15},21/1),({12,3,27,16},28/1),({12,3,26,17},34/1),({12,3,25,18},37/1),({12,3,24,19},35/1),({12,3,23,20},27/1),({12,3,22,21},15/1),({11,4,34,9},1/1),({11,4,33,10},3/1),({11,4,32,11},9/1),({11,4,31,12},16/1),({11,4,30,13},30/1),({11,4,29,14},46/1),({11,4,28,15},65/1),({11,4,27,16},82/1),({11,4,26,17},96/1),({11,4,25,18},97/1),({11,4,24,19},90/1),({11,4,23,20},69/1),({11,4,22,21},37/1),({10,5,35,8},1/1),({10,5,34,9},4/1),({10,5,33,10},10/1),({10,5,32,11},22/1),({10,5,31,12},40/1),({10,5,30,13},65/1),({10,5,29,14},95/1),({10,5,28,15},127/1),({10,5,27,16},155/1),({10,5,26,17},173/1),({10,5,25,18},174/1),({10,5,24,19},156/1),({10,5,23,20},118/1),({10,5,22,21},64/1),({9,6,36,7},1/1),({9,6,35,8},3/1),({9,6,34,9},9/1),({9,6,33,10},18/1),({9,6,32,11},35/1),({9,6,31,12},59/1),({9,6,30,13},91/1),({9,6,29,14},126/1),({9,6,28,15},165/1),({9,6,27,16},194/1),({9,6,26,17},213/1),({9,6,25,18},212/1),({9,6,24,19},188/1),({9,6,23,20},140/1),({9,6,22,21},76/1),({8,7,36,7},1/1),({8,7,35,8},3/1),({8,7,34,9},7/1),({8,7,33,10},15/1),({8,7,32,11},28/1),({8,7,31,12},45/1),({8,7,30,13},68/1),({8,7,29,14},94/1),({8,7,28,15},119/1),({8,7,27,16},140/1),({8,7,26,17},152/1),({8,7,25,18},149/1),({8,7,24,19},132/1),({8,7,23,20},99/1),({8,7,22,21},51/1)}, (6,2) => {}, (8,0) => {}, (8,1) => {({16,3,35,20},1/1),({16,3,34,21},1/1),({16,3,33,22},2/1),({16,3,32,23},3/1),({16,3,31,24},3/1),({16,3,30,25},3/1),({16,3,29,26},3/1),({16,3,28,27},1/1),({15,4,37,18},2/1),({15,4,36,19},4/1),({15,4,35,20},8/1),({15,4,34,21},12/1),({15,4,33,22},18/1),({15,4,32,23},22/1),({15,4,31,24},24/1),({15,4,30,25},22/1),({15,4,29,26},18/1),({15,4,28,27},10/1),({14,5,40,15},1/1),({14,5,39,16},3/1),({14,5,38,17},7/1),({14,5,37,18},16/1),({14,5,36,19},28/1),({14,5,35,20},44/1),({14,5,34,21},62/1),({14,5,33,22},80/1),({14,5,32,23},92/1),({14,5,31,24},96/1),({14,5,30,25},87/1),({14,5,29,26},67/1),({14,5,28,27},37/1),({13,6,41,14},2/1),({13,6,40,15},6/1),({13,6,39,16},16/1),({13,6,38,17},32/1),({13,6,37,18},58/1),({13,6,36,19},92/1),({13,6,35,20},136/1),({13,6,34,21},178/1),({13,6,33,22},218/1),({13,6,32,23},242/1),({13,6,31,24},244/1),({13,6,30,25},218/1),({13,6,29,26},166/1),({13,6,28,27},88/1),({12,7,43,12},1/1),({12,7,42,13},3/1),({12,7,41,14},9/1),({12,7,40,15},21/1),({12,7,39,16},44/1),({12,7,38,17},79/1),({12,7,37,18},130/1),({12,7,36,19},194/1),({12,7,35,20},269/1),({12,7,34,21},342/1),({12,7,33,22},403/1),({12,7,32,23},436/1),({12,7,31,24},432/1),({12,7,30,25},381/1),({12,7,29,26},285/1),({12,7,28,27},152/1),({11,8,43,12},2/1),({11,8,42,13},6/1),({11,8,41,14},16/1),({11,8,40,15},34/1),({11,8,39,16},66/1),({11,8,38,17},112/1),({11,8,37,18},178/1),({11,8,36,19},256/1),({11,8,35,20},344/1),({11,8,34,21},428/1),({11,8,33,22},496/1),({11,8,32,23},528/1),({11,8,31,24},518/1),({11,8,30,25},452/1),({11,8,29,26},336/1),({11,8,28,27},180/1),({10,9,44,11},1/1),({10,9,43,12},2/1),({10,9,42,13},6/1),({10,9,41,14},15/1),({10,9,40,15},29/1),({10,9,39,16},53/1),({10,9,38,17},89/1),({10,9,37,18},135/1),({10,9,36,19},192/1),({10,9,35,20},254/1),({10,9,34,21},310/1),({10,9,33,22},356/1),({10,9,32,23},378/1),({10,9,31,24},365/1),({10,9,30,25},318/1),({10,9,29,26},237/1),({10,9,28,27},125/1)}, (8,2) => {}, (10,0) => {}, (10,1) => {({18,5,40,27},1/1),({18,5,39,28},1/1),({18,5,38,29},2/1),({18,5,37,30},2/1),({18,5,36,31},2/1),({18,5,35,32},2/1),({18,5,34,33},1/1),({17,6,44,23},1/1),({17,6,43,24},2/1),({17,6,42,25},4/1),({17,6,41,26},7/1),({17,6,40,27},12/1),({17,6,39,28},16/1),({17,6,38,29},20/1),({17,6,37,30},21/1),({17,6,36,31},20/1),({17,6,35,32},16/1),({17,6,34,33},9/1),({16,7,46,21},1/1),({16,7,45,22},3/1),({16,7,44,23},8/1),({16,7,43,24},16/1),({16,7,42,25},29/1),({16,7,41,26},44/1),({16,7,40,27},62/1),({16,7,39,28},78/1),({16,7,38,29},91/1),({16,7,37,30},93/1),({16,7,36,31},85/1),({16,7,35,32},65/1),({16,7,34,33},36/1),({15,8,48,19},1/1),({15,8,47,20},3/1),({15,8,46,21},9/1),({15,8,45,22},19/1),({15,8,44,23},38/1),({15,8,43,24},65/1),({15,8,42,25},103/1),({15,8,41,26},146/1),({15,8,40,27},192/1),({15,8,39,28},230/1),({15,8,38,29},255/1),({15,8,37,30},255/1),({15,8,36,31},228/1),({15,8,35,32},171/1),({15,8,34,33},92/1),({14,9,49,18},1/1),({14,9,48,19},4/1),({14,9,47,20},11/1),({14,9,46,21},26/1),({14,9,45,22},51/1),({14,9,44,23},92/1),({14,9,43,24},147/1),({14,9,42,25},218/1),({14,9,41,26},297/1),({14,9,40,27},376/1),({14,9,39,28},438/1),({14,9,38,29},473/1),({14,9,37,30},465/1),({14,9,36,31},410/1),({14,9,35,32},305/1),({14,9,34,33},163/1),({13,10,50,17},1/1),({13,10,49,18},3/1),({13,10,48,19},9/1),({13,10,47,20},21/1),({13,10,46,21},44/1),({13,10,45,22},80/1),({13,10,44,23},135/1),({13,10,43,24},207/1),({13,10,42,25},296/1),({13,10,41,26},391/1),({13,10,40,27},484/1),({13,10,39,28},554/1),({13,10,38,29},589/1),({13,10,37,30},572/1),({13,10,36,31},499/1),({13,10,35,32},370/1),({13,10,34,33},197/1),({12,11,50,17},1/1),({12,11,49,18},3/1),({12,11,48,19},8/1),({12,11,47,20},18/1),({12,11,46,21},36/1),({12,11,45,22},64/1),({12,11,44,23},105/1),({12,11,43,24},158/1),({12,11,42,25},222/1),({12,11,41,26},289/1),({12,11,40,27},353/1),({12,11,39,28},401/1),({12,11,38,29},423/1),({12,11,37,30},408/1),({12,11,36,31},354/1),({12,11,35,32},262/1),({12,11,34,33},140/1)}, (10,2) => {}, (12,0) => {}, (12,1) => {({20,7,41,38},1/1),({19,8,47,32},1/1),({19,8,46,33},2/1),({19,8,45,34},3/1),({19,8,44,35},4/1),({19,8,43,36},5/1),({19,8,42,37},5/1),({19,8,41,38},4/1),({19,8,40,39},2/1),({18,9,51,28},1/1),({18,9,50,29},2/1),({18,9,49,30},5/1),({18,9,48,31},8/1),({18,9,47,32},14/1),({18,9,46,33},20/1),({18,9,45,34},27/1),({18,9,44,35},31/1),({18,9,43,36},34/1),({18,9,42,37},31/1),({18,9,41,38},24/1),({18,9,40,39},13/1),({17,10,53,26},1/1),({17,10,52,27},2/1),({17,10,51,28},6/1),({17,10,50,29},13/1),({17,10,49,30},24/1),({17,10,48,31},39/1),({17,10,47,32},59/1),({17,10,46,33},78/1),({17,10,45,34},97/1),({17,10,44,35},109/1),({17,10,43,36},110/1),({17,10,42,37},99/1),({17,10,41,38},76/1),({17,10,40,39},40/1),({16,11,54,25},1/1),({16,11,53,26},4/1),({16,11,52,27},10/1),({16,11,51,28},21/1),({16,11,50,29},39/1),({16,11,49,30},66/1),({16,11,48,31},99/1),({16,11,47,32},139/1),({16,11,46,33},177/1),({16,11,45,34},211/1),({16,11,44,35},229/1),({16,11,43,36},228/1),({16,11,42,37},201/1),({16,11,41,38},151/1),({16,11,40,39},81/1),({15,12,55,24},1/1),({15,12,54,25},3/1),({15,12,53,26},8/1),({15,12,52,27},18/1),({15,12,51,28},36/1),({15,12,50,29},62/1),({15,12,49,30},100/1),({15,12,48,31},145/1),({15,12,47,32},196/1),({15,12,46,33},246/1),({15,12,45,34},285/1),({15,12,44,35},304/1),({15,12,43,36},299/1),({15,12,42,37},261/1),({15,12,41,38},194/1),({15,12,40,39},104/1),({14,13,55,24},1/1),({14,13,54,25},3/1),({14,13,53,26},8/1),({14,13,52,27},16/1),({14,13,51,28},31/1),({14,13,50,29},52/1),({14,13,49,30},80/1),({14,13,48,31},114/1),({14,13,47,32},153/1),({14,13,46,33},187/1),({14,13,45,34},215/1),({14,13,44,35},228/1),({14,13,43,36},221/1),({14,13,42,37},193/1),({14,13,41,38},144/1),({14,13,40,39},75/1)}, (12,2) => {}, (14,0) => {}, (14,1) => {({20,11,52,39},1/1),({20,11,51,40},1/1),({20,11,50,41},2/1),({20,11,49,42},2/1),({20,11,48,43},3/1),({20,11,47,44},2/1),({20,11,46,45},1/1),({19,12,56,35},1/1),({19,12,55,36},2/1),({19,12,54,37},4/1),({19,12,53,38},7/1),({19,12,52,39},10/1),({19,12,51,40},13/1),({19,12,50,41},16/1),({19,12,49,42},16/1),({19,12,48,43},15/1),({19,12,47,44},12/1),({19,12,46,45},6/1),({18,13,58,33},1/1),({18,13,57,34},2/1),({18,13,56,35},5/1),({18,13,55,36},9/1),({18,13,54,37},16/1),({18,13,53,38},23/1),({18,13,52,39},32/1),({18,13,51,40},39/1),({18,13,50,41},45/1),({18,13,49,42},45/1),({18,13,48,43},41/1),({18,13,47,44},31/1),({18,13,46,45},17/1),({17,14,58,33},2/1),({17,14,57,34},5/1),({17,14,56,35},10/1),({17,14,55,36},18/1),({17,14,54,37},29/1),({17,14,53,38},41/1),({17,14,52,39},54/1),({17,14,51,40},64/1),({17,14,50,41},70/1),({17,14,49,42},70/1),({17,14,48,43},62/1),({17,14,47,44},46/1),({17,14,46,45},25/1),({16,15,59,32},1/1),({16,15,58,33},2/1),({16,15,57,34},5/1),({16,15,56,35},10/1),({16,15,55,36},16/1),({16,15,54,37},25/1),({16,15,53,38},35/1),({16,15,52,39},45/1),({16,15,51,40},52/1),({16,15,50,41},57/1),({16,15,49,42},56/1),({16,15,48,43},50/1),({16,15,47,44},37/1),({16,15,46,45},19/1)}, (16,0) => {}, (14,2) => {}, (16,1) => {({20,15,59,44},1/1),({20,15,58,45},1/1),({20,15,57,46},2/1),({20,15,56,47},2/1),({20,15,55,48},3/1),({20,15,54,49},2/1),({20,15,53,50},2/1),({20,15,52,51},1/1),({19,16,61,42},1/1),({19,16,60,43},2/1),({19,16,59,44},3/1),({19,16,58,45},4/1),({19,16,57,46},5/1),({19,16,56,47},6/1),({19,16,55,48},6/1),({19,16,54,49},5/1),({19,16,53,50},4/1),({19,16,52,51},2/1),({18,17,61,42},1/1),({18,17,60,43},2/1),({18,17,59,44},3/1),({18,17,58,45},4/1),({18,17,57,46},5/1),({18,17,56,47},6/1),({18,17,55,48},6/1),({18,17,54,49},5/1),({18,17,53,50},4/1),({18,17,52,51},2/1)}, (18,0) => {}, (16,2) => {}, (18,1) => {({20,19,62,53},1/1)}, (18,2) => {}, (20,1) => {}, (1,0) => {({3,0,6,1},1/1),({2,1,7,0},1/1),({2,1,6,1},1/1)}, (1,1) => {}, (3,0) => {({5,2,17,2},1/1),({5,2,15,4},1/1),({5,2,14,5},1/1),({5,2,13,6},1/1),({5,2,12,7},1/1),({5,2,11,8},1/1),({4,3,16,3},1/1),({4,3,15,4},1/1),({4,3,14,5},1/1),({4,3,13,6},2/1),({4,3,12,7},2/1),({4,3,11,8},1/1),({4,3,10,9},1/1)}, (3,1) => {({9,0,17,8},1/1),({9,0,16,9},1/1),({9,0,15,10},1/1),({9,0,14,11},1/1),({9,0,13,12},1/1),({8,1,20,5},1/1),({8,1,19,6},1/1),({8,1,18,7},2/1),({8,1,17,8},3/1),({8,1,16,9},3/1),({8,1,15,10},3/1),({8,1,14,11},3/1),({8,1,13,12},1/1),({7,2,20,5},1/1),({7,2,19,6},2/1),({7,2,18,7},3/1),({7,2,17,8},4/1),({7,2,16,9},5/1),({7,2,15,10},5/1),({7,2,14,11},4/1),({7,2,13,12},2/1),({6,3,22,3},1/1),({6,3,21,4},1/1),({6,3,20,5},2/1),({6,3,19,6},3/1),({6,3,18,7},4/1),({6,3,17,8},4/1),({6,3,16,9},5/1),({6,3,15,10},4/1),({6,3,14,11},3/1),({6,3,13,12},2/1),({5,4,21,4},1/1),({5,4,20,5},1/1),({5,4,19,6},1/1),({5,4,18,7},2/1),({5,4,17,8},3/1),({5,4,16,9},2/1),({5,4,15,10},2/1),({5,4,14,11},2/1),({5,4,13,12},1/1)}, (3,2) => {}, (5,0) => {}, (5,1) => {({13,0,20,17},1/1),({12,1,25,12},1/1),({12,1,24,13},2/1),({12,1,23,14},3/1),({12,1,22,15},3/1),({12,1,21,16},3/1),({12,1,20,17},3/1),({12,1,19,18},2/1),({11,2,28,9},1/1),({11,2,27,10},2/1),({11,2,26,11},4/1),({11,2,25,12},7/1),({11,2,24,13},11/1),({11,2,23,14},14/1),({11,2,22,15},16/1),({11,2,21,16},15/1),({11,2,20,17},12/1),({11,2,19,18},7/1),({10,3,29,8},2/1),({10,3,28,9},4/1),({10,3,27,10},10/1),({10,3,26,11},17/1),({10,3,25,12},25/1),({10,3,24,13},33/1),({10,3,23,14},41/1),({10,3,22,15},42/1),({10,3,21,16},40/1),({10,3,20,17},31/1),({10,3,19,18},16/1),({9,4,31,6},1/1),({9,4,30,7},3/1),({9,4,29,8},7/1),({9,4,28,9},14/1),({9,4,27,10},24/1),({9,4,26,11},38/1),({9,4,25,12},52/1),({9,4,24,13},66/1),({9,4,23,14},75/1),({9,4,22,15},78/1),({9,4,21,16},71/1),({9,4,20,17},54/1),({9,4,19,18},29/1),({8,5,31,6},2/1),({8,5,30,7},5/1),({8,5,29,8},11/1),({8,5,28,9},21/1),({8,5,27,10},35/1),({8,5,26,11},50/1),({8,5,25,12},69/1),({8,5,24,13},83/1),({8,5,23,14},93/1),({8,5,22,15},94/1),({8,5,21,16},84/1),({8,5,20,17},63/1),({8,5,19,18},35/1),({7,6,32,5},1/1),({7,6,31,6},2/1),({7,6,30,7},5/1),({7,6,29,8},10/1),({7,6,28,9},17/1),({7,6,27,10},27/1),({7,6,26,11},39/1),({7,6,25,12},50/1),({7,6,24,13},61/1),({7,6,23,14},68/1),({7,6,22,15},66/1),({7,6,21,16},59/1),({7,6,20,17},46/1),({7,6,19,18},24/1)}};
--dw stands for dominant weights
dw1126 = new HashTable from {(5,2) => {}, (7,0) => {}, (7,1) => {({15,2,31,18},1/1)}, (9,0) => {}, (7,2) => {}, (9,1) => {({17,4,38,23},1/1)}, (9,2) => {}, (11,0) => {}, (11,1) => {({19,6,41,32},1/1)}, (11,2) => {}, (13,0) => {}, (13,1) => {({20,9,47,38},1/1)}, (13,2) => {}, (15,0) => {}, (15,1) => {({20,13,56,41},1/1)}, (15,2) => {}, (17,0) => {}, (17,1) => {({20,17,61,48},1/1)}, (17,2) => {}, (19,1) => {}, (0,0) => {({1,0,1,0},1/1)}, (2,0) => {({4,1,12,1},1/1)}, (2,1) => {({7,0,14,5},1/1)}, (2,2) => {}, (4,0) => {}, (4,1) => {({11,0,19,12},1/1)}, (4,2) => {}, (6,0) => {}, (6,1) => {({14,1,26,17},1/1)}, (6,2) => {}, (8,0) => {}, (8,1) => {({16,3,35,20},1/1)}, (8,2) => {}, (10,0) => {}, (10,1) => {({18,5,40,27},1/1)}, (10,2) => {}, (12,0) => {}, (12,1) => {({20,7,41,38},1/1)}, (12,2) => {}, (14,0) => {}, (14,1) => {({20,11,52,39},1/1)}, (16,0) => {}, (14,2) => {}, (16,1) => {({20,15,59,44},1/1)}, (18,0) => {}, (16,2) => {}, (18,1) => {({20,19,62,53},1/1)}, (18,2) => {}, (20,1) => {}, (1,0) => {({3,0,6,1},1/1)}, (1,1) => {}, (3,0) => {({5,2,17,2},1/1)}, (3,1) => {({9,0,17,8},1/1)}, (5,0) => {}, (3,2) => {}, (5,1) => {({13,0,20,17},1/1)}};
end;