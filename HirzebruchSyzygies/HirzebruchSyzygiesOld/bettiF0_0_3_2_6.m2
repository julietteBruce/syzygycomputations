A := QQ[t_0,t_1,t_2,t_3];
--tb stands for Total Betti numbers
tb0326 = new HashTable from {(5,2) => 0, (7,0) => 0, (7,1) => 360672/1, (7,2) => 0, (9,0) => 0, (9,1) => 612612/1, (9,2) => 0, (11,0) => 0, (11,1) => 411060/1, (11,2) => 0, (13,0) => 0, (13,1) => 104856/1, (15,0) => 0, (13,2) => 0, (15,1) => 7956/1, (17,0) => 0, (15,2) => 0, (17,1) => 14/1, (17,2) => 18/1, (19,2) => 0, (0,0) => 4/1, (2,0) => 306/1, (2,1) => 270/1, (2,2) => 0, (4,0) => 0, (4,1) => 19176/1, (4,2) => 0, (6,0) => 0, (6,1) => 189720/1, (6,2) => 0, (8,0) => 0, (8,1) => 529516/1, (8,2) => 0, (10,0) => 0, (10,1) => 563108/1, (10,2) => 0, (12,0) => 0, (12,1) => 236232/1, (12,2) => 0, (14,0) => 0, (14,1) => 34680/1, (16,0) => 0, (14,2) => 0, (16,1) => 1050/1, (18,0) => 0, (16,2) => 0, (18,1) => 0, (18,2) => 2/1, (1,0) => 54/1, (1,1) => 16/1, (3,0) => 816/1, (3,1) => 2142/1, (3,2) => 0, (5,0) => 0, (5,1) => 73848/1};
--sb represents the betti numbers as sums of Schur functors
sb0326 = new HashTable from {(5,2) => {}, (7,0) => {}, (7,1) => {({14,2,31,20},1/1),({14,2,30,21},2/1),({14,2,29,22},2/1),({14,2,28,23},2/1),({14,2,27,24},2/1),({14,2,26,25},1/1),({13,3,35,16},1/1),({13,3,34,17},3/1),({13,3,33,18},5/1),({13,3,32,19},9/1),({13,3,31,20},13/1),({13,3,30,21},16/1),({13,3,29,22},18/1),({13,3,28,23},17/1),({13,3,27,24},14/1),({13,3,26,25},8/1),({12,4,38,13},1/1),({12,4,37,14},2/1),({12,4,36,15},6/1),({12,4,35,16},12/1),({12,4,34,17},22/1),({12,4,33,18},34/1),({12,4,32,19},48/1),({12,4,31,20},61/1),({12,4,30,21},71/1),({12,4,29,22},73/1),({12,4,28,23},67/1),({12,4,27,24},51/1),({12,4,26,25},28/1),({11,5,39,12},2/1),({11,5,38,13},5/1),({11,5,37,14},13/1),({11,5,36,15},26/1),({11,5,35,16},45/1),({11,5,34,17},71/1),({11,5,33,18},103/1),({11,5,32,19},134/1),({11,5,31,20},162/1),({11,5,30,21},179/1),({11,5,29,22},179/1),({11,5,28,23},160/1),({11,5,27,24},121/1),({11,5,26,25},64/1),({10,6,41,10},1/1),({10,6,40,11},3/1),({10,6,39,12},7/1),({10,6,38,13},17/1),({10,6,37,14},33/1),({10,6,36,15},58/1),({10,6,35,16},93/1),({10,6,34,17},138/1),({10,6,33,18},186/1),({10,6,32,19},236/1),({10,6,31,20},274/1),({10,6,30,21},295/1),({10,6,29,22},290/1),({10,6,28,23},255/1),({10,6,27,24},189/1),({10,6,26,25},102/1),({9,7,41,10},1/1),({9,7,40,11},4/1),({9,7,39,12},10/1),({9,7,38,13},21/1),({9,7,37,14},40/1),({9,7,36,15},67/1),({9,7,35,16},105/1),({9,7,34,17},149/1),({9,7,33,18},199/1),({9,7,32,19},245/1),({9,7,31,20},282/1),({9,7,30,21},299/1),({9,7,29,22},291/1),({9,7,28,23},253/1),({9,7,27,24},188/1),({9,7,26,25},100/1),({8,8,42,9},1/1),({8,8,41,10},1/1),({8,8,40,11},3/1),({8,8,39,12},6/1),({8,8,38,13},12/1),({8,8,37,14},20/1),({8,8,36,15},34/1),({8,8,35,16},49/1),({8,8,34,17},70/1),({8,8,33,18},90/1),({8,8,32,19},110/1),({8,8,31,20},124/1),({8,8,30,21},133/1),({8,8,29,22},126/1),({8,8,28,23},110/1),({8,8,27,24},82/1),({8,8,26,25},43/1)}, (7,2) => {}, (9,0) => {}, (9,1) => {({16,4,39,24},1/1),({16,4,38,25},2/1),({16,4,37,26},3/1),({16,4,36,27},4/1),({16,4,35,28},5/1),({16,4,34,29},5/1),({16,4,33,30},4/1),({16,4,32,31},2/1),({15,5,42,21},1/1),({15,5,41,22},3/1),({15,5,40,23},7/1),({15,5,39,24},13/1),({15,5,38,25},20/1),({15,5,37,26},27/1),({15,5,36,27},33/1),({15,5,35,28},36/1),({15,5,34,29},34/1),({15,5,33,30},26/1),({15,5,32,31},14/1),({14,6,44,19},1/1),({14,6,43,20},5/1),({14,6,42,21},12/1),({14,6,41,22},24/1),({14,6,40,23},41/1),({14,6,39,24},64/1),({14,6,38,25},89/1),({14,6,37,26},112/1),({14,6,36,27},127/1),({14,6,35,28},131/1),({14,6,34,29},119/1),({14,6,33,30},91/1),({14,6,32,31},49/1),({13,7,46,17},1/1),({13,7,45,18},4/1),({13,7,44,19},11/1),({13,7,43,20},24/1),({13,7,42,21},47/1),({13,7,41,22},81/1),({13,7,40,23},126/1),({13,7,39,24},178/1),({13,7,38,25},232/1),({13,7,37,26},278/1),({13,7,36,27},306/1),({13,7,35,28},305/1),({13,7,34,29},271/1),({13,7,33,30},204/1),({13,7,32,31},110/1),({12,8,47,16},1/1),({12,8,46,17},4/1),({12,8,45,18},12/1),({12,8,44,19},27/1),({12,8,43,20},54/1),({12,8,42,21},95/1),({12,8,41,22},153/1),({12,8,40,23},224/1),({12,8,39,24},305/1),({12,8,38,25},383/1),({12,8,37,26},447/1),({12,8,36,27},480/1),({12,8,35,28},472/1),({12,8,34,29},414/1),({12,8,33,30},309/1),({12,8,32,31},165/1),({11,9,48,15},1/1),({11,9,47,16},3/1),({11,9,46,17},8/1),({11,9,45,18},18/1),({11,9,44,19},37/1),({11,9,43,20},67/1),({11,9,42,21},112/1),({11,9,41,22},171/1),({11,9,40,23},244/1),({11,9,39,24},322/1),({11,9,38,25},397/1),({11,9,37,26},454/1),({11,9,36,27},482/1),({11,9,35,28},468/1),({11,9,34,29},408/1),({11,9,33,30},302/1),({11,9,32,31},161/1),({10,10,47,16},1/1),({10,10,46,17},3/1),({10,10,45,18},8/1),({10,10,44,19},16/1),({10,10,43,20},30/1),({10,10,42,21},49/1),({10,10,41,22},75/1),({10,10,40,23},105/1),({10,10,39,24},139/1),({10,10,38,25},169/1),({10,10,37,26},193/1),({10,10,36,27},203/1),({10,10,35,28},198/1),({10,10,34,29},172/1),({10,10,33,30},127/1),({10,10,32,31},67/1)}, (9,2) => {}, (11,0) => {}, (11,1) => {({18,6,44,31},1/1),({18,6,43,32},1/1),({18,6,42,33},2/1),({18,6,41,34},3/1),({18,6,40,35},3/1),({18,6,39,36},2/1),({18,6,38,37},1/1),({17,7,47,28},1/1),({17,7,46,29},3/1),({17,7,45,30},7/1),({17,7,44,31},11/1),({17,7,43,32},16/1),({17,7,42,33},20/1),({17,7,41,34},22/1),({17,7,40,35},22/1),({17,7,39,36},17/1),({17,7,38,37},9/1),({16,8,49,26},2/1),({16,8,48,27},6/1),({16,8,47,28},13/1),({16,8,46,29},24/1),({16,8,45,30},39/1),({16,8,44,31},57/1),({16,8,43,32},73/1),({16,8,42,33},85/1),({16,8,41,34},89/1),({16,8,40,35},82/1),({16,8,39,36},63/1),({16,8,38,37},34/1),({15,9,51,24},1/1),({15,9,50,25},5/1),({15,9,49,26},12/1),({15,9,48,27},26/1),({15,9,47,28},49/1),({15,9,46,29},79/1),({15,9,45,30},116/1),({15,9,44,31},156/1),({15,9,43,32},191/1),({15,9,42,33},213/1),({15,9,41,34},216/1),({15,9,40,35},193/1),({15,9,39,36},146/1),({15,9,38,37},80/1),({14,10,52,23},2/1),({14,10,51,24},6/1),({14,10,50,25},15/1),({14,10,49,26},32/1),({14,10,48,27},60/1),({14,10,47,28},99/1),({14,10,46,29},151/1),({14,10,45,30},210/1),({14,10,44,31},269/1),({14,10,43,32},319/1),({14,10,42,33},348/1),({14,10,41,34},344/1),({14,10,40,35},305/1),({14,10,39,36},229/1),({14,10,38,37},122/1),({13,11,53,22},1/1),({13,11,52,23},3/1),({13,11,51,24},9/1),({13,11,50,25},20/1),({13,11,49,26},40/1),({13,11,48,27},70/1),({13,11,47,28},112/1),({13,11,46,29},165/1),({13,11,45,30},224/1),({13,11,44,31},281/1),({13,11,43,32},327/1),({13,11,42,33},351/1),({13,11,41,34},345/1),({13,11,40,35},302/1),({13,11,39,36},225/1),({13,11,38,37},120/1),({12,12,53,22},1/1),({12,12,52,23},2/1),({12,12,51,24},5/1),({12,12,50,25},11/1),({12,12,49,26},20/1),({12,12,48,27},34/1),({12,12,47,28},53/1),({12,12,46,29},75/1),({12,12,45,30},100/1),({12,12,44,31},125/1),({12,12,43,32},143/1),({12,12,42,33},152/1),({12,12,41,34},150/1),({12,12,40,35},131/1),({12,12,39,36},96/1),({12,12,38,37},52/1)}, (11,2) => {}, (13,0) => {}, (13,1) => {({20,8,45,42},1/1),({19,9,50,37},1/1),({19,9,49,38},2/1),({19,9,48,39},3/1),({19,9,47,40},4/1),({19,9,46,41},4/1),({19,9,45,42},3/1),({19,9,44,43},2/1),({18,10,53,34},1/1),({18,10,52,35},3/1),({18,10,51,36},6/1),({18,10,50,37},10/1),({18,10,49,38},15/1),({18,10,48,39},19/1),({18,10,47,40},21/1),({18,10,46,41},20/1),({18,10,45,42},16/1),({18,10,44,43},9/1),({17,11,55,32},1/1),({17,11,54,33},4/1),({17,11,53,34},8/1),({17,11,52,35},16/1),({17,11,51,36},27/1),({17,11,50,37},38/1),({17,11,49,38},50/1),({17,11,48,39},59/1),({17,11,47,40},61/1),({17,11,46,41},57/1),({17,11,45,42},44/1),({17,11,44,43},23/1),({16,12,56,31},2/1),({16,12,55,32},5/1),({16,12,54,33},11/1),({16,12,53,34},22/1),({16,12,52,35},37/1),({16,12,51,36},55/1),({16,12,50,37},76/1),({16,12,49,38},94/1),({16,12,48,39},106/1),({16,12,47,40},109/1),({16,12,46,41},98/1),({16,12,45,42},74/1),({16,12,44,43},41/1),({15,13,57,30},1/1),({15,13,56,31},3/1),({15,13,55,32},8/1),({15,13,54,33},16/1),({15,13,53,34},29/1),({15,13,52,35},46/1),({15,13,51,36},66/1),({15,13,50,37},87/1),({15,13,49,38},105/1),({15,13,48,39},116/1),({15,13,47,40},116/1),({15,13,46,41},103/1),({15,13,45,42},78/1),({15,13,44,43},42/1),({14,14,57,30},1/1),({14,14,56,31},2/1),({14,14,55,32},4/1),({14,14,54,33},8/1),({14,14,53,34},13/1),({14,14,52,35},21/1),({14,14,51,36},30/1),({14,14,50,37},38/1),({14,14,49,38},46/1),({14,14,48,39},52/1),({14,14,47,40},50/1),({14,14,46,41},45/1),({14,14,45,42},35/1),({14,14,44,43},18/1)}, (15,0) => {}, (13,2) => {}, (15,1) => {({20,12,54,45},1/1),({20,12,53,46},1/1),({20,12,52,47},1/1),({20,12,51,48},1/1),({20,12,50,49},1/1),({19,13,57,42},1/1),({19,13,56,43},2/1),({19,13,55,44},3/1),({19,13,54,45},4/1),({19,13,53,46},5/1),({19,13,52,47},5/1),({19,13,51,48},4/1),({19,13,50,49},2/1),({18,14,59,40},1/1),({18,14,58,41},2/1),({18,14,57,42},4/1),({18,14,56,43},6/1),({18,14,55,44},9/1),({18,14,54,45},11/1),({18,14,53,46},12/1),({18,14,52,47},11/1),({18,14,51,48},9/1),({18,14,50,49},5/1),({17,15,60,39},1/1),({17,15,59,40},2/1),({17,15,58,41},4/1),({17,15,57,42},6/1),({17,15,56,43},9/1),({17,15,55,44},12/1),({17,15,54,45},14/1),({17,15,53,46},14/1),({17,15,52,47},13/1),({17,15,51,48},10/1),({17,15,50,49},6/1),({16,16,59,40},1/1),({16,16,58,41},2/1),({16,16,57,42},3/1),({16,16,56,43},4/1),({16,16,55,44},6/1),({16,16,54,45},7/1),({16,16,53,46},7/1),({16,16,52,47},6/1),({16,16,51,48},5/1),({16,16,50,49},3/1)}, (17,0) => {}, (15,2) => {}, (17,1) => {({18,18,62,49},1/1)}, (17,2) => {({20,18,61,56},1/1)}, (19,2) => {}, (0,0) => {({0,0,3,0},1/1)}, (2,0) => {({4,0,12,3},1/1),({4,0,11,4},1/1),({4,0,10,5},1/1),({4,0,9,6},1/1),({4,0,8,7},1/1),({3,1,13,2},1/1),({3,1,12,3},1/1),({3,1,11,4},1/1),({3,1,10,5},1/1),({3,1,9,6},1/1),({3,1,8,7},1/1),({2,2,12,3},1/1),({2,2,11,4},1/1),({2,2,10,5},1/1),({2,2,9,6},1/1),({2,2,8,7},1/1)}, (2,1) => {({4,2,20,1},1/1),({4,2,19,2},1/1),({4,2,18,3},1/1),({4,2,17,4},1/1),({4,2,16,5},1/1),({4,2,15,6},1/1)}, (2,2) => {}, (4,0) => {}, (4,1) => {({10,0,19,14},1/1),({9,1,23,10},1/1),({9,1,22,11},1/1),({9,1,21,12},2/1),({9,1,20,13},2/1),({9,1,19,14},2/1),({9,1,18,15},2/1),({9,1,17,16},1/1),({8,2,27,6},1/1),({8,2,26,7},1/1),({8,2,25,8},3/1),({8,2,24,9},4/1),({8,2,23,10},6/1),({8,2,22,11},7/1),({8,2,21,12},9/1),({8,2,20,13},8/1),({8,2,19,14},8/1),({8,2,18,15},6/1),({8,2,17,16},3/1),({7,3,29,4},1/1),({7,3,28,5},2/1),({7,3,27,6},3/1),({7,3,26,7},6/1),({7,3,25,8},9/1),({7,3,24,9},12/1),({7,3,23,10},15/1),({7,3,22,11},18/1),({7,3,21,12},19/1),({7,3,20,13},19/1),({7,3,19,14},16/1),({7,3,18,15},12/1),({7,3,17,16},7/1),({6,4,29,4},1/1),({6,4,28,5},2/1),({6,4,27,6},5/1),({6,4,26,7},7/1),({6,4,25,8},12/1),({6,4,24,9},15/1),({6,4,23,10},20/1),({6,4,22,11},21/1),({6,4,21,12},23/1),({6,4,20,13},21/1),({6,4,19,14},19/1),({6,4,18,15},13/1),({6,4,17,16},7/1),({5,5,30,3},1/1),({5,5,29,4},1/1),({5,5,28,5},2/1),({5,5,27,6},3/1),({5,5,26,7},5/1),({5,5,25,8},6/1),({5,5,24,9},9/1),({5,5,23,10},9/1),({5,5,22,11},11/1),({5,5,21,12},11/1),({5,5,20,13},11/1),({5,5,19,14},9/1),({5,5,18,15},7/1),({5,5,17,16},3/1)}, (4,2) => {}, (6,0) => {}, (6,1) => {({13,1,26,19},1/1),({13,1,25,20},1/1),({13,1,24,21},1/1),({12,2,30,15},2/1),({12,2,29,16},3/1),({12,2,28,17},5/1),({12,2,27,18},6/1),({12,2,26,19},8/1),({12,2,25,20},8/1),({12,2,24,21},6/1),({12,2,23,22},3/1),({11,3,34,11},1/1),({11,3,33,12},2/1),({11,3,32,13},5/1),({11,3,31,14},9/1),({11,3,30,15},15/1),({11,3,29,16},22/1),({11,3,28,17},29/1),({11,3,27,18},34/1),({11,3,26,19},36/1),({11,3,25,20},33/1),({11,3,24,21},26/1),({11,3,23,22},14/1),({10,4,36,9},1/1),({10,4,35,10},2/1),({10,4,34,11},6/1),({10,4,33,12},12/1),({10,4,32,13},23/1),({10,4,31,14},35/1),({10,4,30,15},53/1),({10,4,29,16},69/1),({10,4,28,17},86/1),({10,4,27,18},94/1),({10,4,26,19},96/1),({10,4,25,20},85/1),({10,4,24,21},65/1),({10,4,23,22},35/1),({9,5,37,8},1/1),({9,5,36,9},3/1),({9,5,35,10},8/1),({9,5,34,11},16/1),({9,5,33,12},30/1),({9,5,32,13},49/1),({9,5,31,14},74/1),({9,5,30,15},101/1),({9,5,29,16},129/1),({9,5,28,17},151/1),({9,5,27,18},164/1),({9,5,26,19},161/1),({9,5,25,20},142/1),({9,5,24,21},106/1),({9,5,23,22},57/1),({8,6,38,7},1/1),({8,6,37,8},2/1),({8,6,36,9},6/1),({8,6,35,10},11/1),({8,6,34,11},23/1),({8,6,33,12},37/1),({8,6,32,13},60/1),({8,6,31,14},84/1),({8,6,30,15},115/1),({8,6,29,16},140/1),({8,6,28,17},163/1),({8,6,27,18},171/1),({8,6,26,19},169/1),({8,6,25,20},146/1),({8,6,24,21},109/1),({8,6,23,22},57/1),({7,7,37,8},1/1),({7,7,36,9},2/1),({7,7,35,10},6/1),({7,7,34,11},10/1),({7,7,33,12},18/1),({7,7,32,13},26/1),({7,7,31,14},39/1),({7,7,30,15},50/1),({7,7,29,16},63/1),({7,7,28,17},70/1),({7,7,27,18},76/1),({7,7,26,19},72/1),({7,7,25,20},64/1),({7,7,24,21},46/1),({7,7,23,22},25/1)}, (6,2) => {}, (8,0) => {}, (8,1) => {({15,3,35,22},2/1),({15,3,34,23},2/1),({15,3,33,24},3/1),({15,3,32,25},4/1),({15,3,31,26},4/1),({15,3,30,27},3/1),({15,3,29,28},2/1),({14,4,39,18},1/1),({14,4,38,19},2/1),({14,4,37,20},6/1),({14,4,36,21},10/1),({14,4,35,22},16/1),({14,4,34,23},22/1),({14,4,33,24},28/1),({14,4,32,25},29/1),({14,4,31,26},28/1),({14,4,30,27},22/1),({14,4,29,28},12/1),({13,5,41,16},1/1),({13,5,40,17},4/1),({13,5,39,18},10/1),({13,5,38,19},20/1),({13,5,37,20},35/1),({13,5,36,21},54/1),({13,5,35,22},75/1),({13,5,34,23},94/1),({13,5,33,24},108/1),({13,5,32,25},111/1),({13,5,31,26},101/1),({13,5,30,27},77/1),({13,5,29,28},42/1),({12,6,43,14},1/1),({12,6,42,15},3/1),({12,6,41,16},10/1),({12,6,40,17},20/1),({12,6,39,18},41/1),({12,6,38,19},69/1),({12,6,37,20},109/1),({12,6,36,21},152/1),({12,6,35,22},201/1),({12,6,34,23},238/1),({12,6,33,24},263/1),({12,6,32,25},262/1),({12,6,31,26},234/1),({12,6,30,27},174/1),({12,6,29,28},95/1),({11,7,44,13},1/1),({11,7,43,14},4/1),({11,7,42,15},11/1),({11,7,41,16},24/1),({11,7,40,17},48/1),({11,7,39,18},84/1),({11,7,38,19},134/1),({11,7,37,20},196/1),({11,7,36,21},267/1),({11,7,35,22},334/1),({11,7,34,23},390/1),({11,7,33,24},418/1),({11,7,32,25},410/1),({11,7,31,26},360/1),({11,7,30,27},269/1),({11,7,29,28},142/1),({10,8,45,12},1/1),({10,8,44,13},2/1),({10,8,43,14},7/1),({10,8,42,15},15/1),({10,8,41,16},33/1),({10,8,40,17},58/1),({10,8,39,18},99/1),({10,8,38,19},149/1),({10,8,37,20},215/1),({10,8,36,21},280/1),({10,8,35,22},348/1),({10,8,34,23},395/1),({10,8,33,24},422/1),({10,8,32,25},407/1),({10,8,31,26},356/1),({10,8,30,27},263/1),({10,8,29,28},141/1),({9,9,44,13},2/1),({9,9,43,14},3/1),({9,9,42,15},8/1),({9,9,41,16},15/1),({9,9,40,17},28/1),({9,9,39,18},43/1),({9,9,38,19},69/1),({9,9,37,20},93/1),({9,9,36,21},124/1),({9,9,35,22},150/1),({9,9,34,23},172/1),({9,9,33,24},178/1),({9,9,32,25},176/1),({9,9,31,26},150/1),({9,9,30,27},111/1),({9,9,29,28},60/1)}, (8,2) => {}, (10,0) => {}, (10,1) => {({17,5,42,27},1/1),({17,5,41,28},1/1),({17,5,40,29},3/1),({17,5,39,30},4/1),({17,5,38,31},4/1),({17,5,37,32},4/1),({17,5,36,33},4/1),({17,5,35,34},2/1),({16,6,44,25},3/1),({16,6,43,26},6/1),({16,6,42,27},11/1),({16,6,41,28},17/1),({16,6,40,29},25/1),({16,6,39,30},30/1),({16,6,38,31},33/1),({16,6,37,32},30/1),({16,6,36,33},24/1),({16,6,35,34},14/1),({15,7,47,22},1/1),({15,7,46,23},4/1),({15,7,45,24},10/1),({15,7,44,25},21/1),({15,7,43,26},37/1),({15,7,42,27},58/1),({15,7,41,28},81/1),({15,7,40,29},103/1),({15,7,39,30},118/1),({15,7,38,31},122/1),({15,7,37,32},111/1),({15,7,36,33},85/1),({15,7,35,34},46/1),({14,8,48,21},3/1),({14,8,47,22},8/1),({14,8,46,23},21/1),({14,8,45,24},40/1),({14,8,44,25},72/1),({14,8,43,26},112/1),({14,8,42,27},164/1),({14,8,41,28},212/1),({14,8,40,29},258/1),({14,8,39,30},284/1),({14,8,38,31},286/1),({14,8,37,32},254/1),({14,8,36,33},193/1),({14,8,35,34},102/1),({13,9,50,19},1/1),({13,9,49,20},4/1),({13,9,48,21},10/1),({13,9,47,22},24/1),({13,9,46,23},48/1),({13,9,45,24},86/1),({13,9,44,25},139/1),({13,9,43,26},207/1),({13,9,42,27},282/1),({13,9,41,28},359/1),({13,9,40,29},419/1),({13,9,39,30},452/1),({13,9,38,31},446/1),({13,9,37,32},393/1),({13,9,36,33},292/1),({13,9,35,34},157/1),({12,10,50,19},2/1),({12,10,49,20},5/1),({12,10,48,21},15/1),({12,10,47,22},30/1),({12,10,46,23},59/1),({12,10,45,24},98/1),({12,10,44,25},156/1),({12,10,43,26},221/1),({12,10,42,27},298/1),({12,10,41,28},366/1),({12,10,40,29},426/1),({12,10,39,30},451/1),({12,10,38,31},442/1),({12,10,37,32},384/1),({12,10,36,33},287/1),({12,10,35,34},153/1),({11,11,51,18},1/1),({11,11,50,19},1/1),({11,11,49,20},4/1),({11,11,48,21},8/1),({11,11,47,22},16/1),({11,11,46,23},27/1),({11,11,45,24},48/1),({11,11,44,25},69/1),({11,11,43,26},101/1),({11,11,42,27},131/1),({11,11,41,28},162/1),({11,11,40,29},183/1),({11,11,39,30},197/1),({11,11,38,31},187/1),({11,11,37,32},165/1),({11,11,36,33},122/1),({11,11,35,34},64/1)}, (10,2) => {}, (12,0) => {}, (12,1) => {({19,7,45,36},1/1),({19,7,44,37},1/1),({19,7,43,38},1/1),({19,7,42,39},1/1),({19,7,41,40},1/1),({18,8,49,32},1/1),({18,8,48,33},2/1),({18,8,47,34},5/1),({18,8,46,35},7/1),({18,8,45,36},10/1),({18,8,44,37},11/1),({18,8,43,38},11/1),({18,8,42,39},9/1),({18,8,41,40},5/1),({17,9,51,30},2/1),({17,9,50,31},5/1),({17,9,49,32},11/1),({17,9,48,33},19/1),({17,9,47,34},29/1),({17,9,46,35},39/1),({17,9,45,36},47/1),({17,9,44,37},50/1),({17,9,43,38},47/1),({17,9,42,39},36/1),({17,9,41,40},20/1),({16,10,53,28},2/1),({16,10,52,29},5/1),({16,10,51,30},13/1),({16,10,50,31},24/1),({16,10,49,32},43/1),({16,10,48,33},64/1),({16,10,47,34},90/1),({16,10,46,35},111/1),({16,10,45,36},128/1),({16,10,44,37},130/1),({16,10,43,38},119/1),({16,10,42,39},90/1),({16,10,41,40},49/1),({15,11,54,27},2/1),({15,11,53,28},6/1),({15,11,52,29},15/1),({15,11,51,30},30/1),({15,11,50,31},54/1),({15,11,49,32},85/1),({15,11,48,33},123/1),({15,11,47,34},161/1),({15,11,46,35},195/1),({15,11,45,36},215/1),({15,11,44,37},216/1),({15,11,43,38},192/1),({15,11,42,39},145/1),({15,11,41,40},78/1),({14,12,55,26},2/1),({14,12,54,27},4/1),({14,12,53,28},11/1),({14,12,52,29},21/1),({14,12,51,30},41/1),({14,12,50,31},65/1),({14,12,49,32},100/1),({14,12,48,33},136/1),({14,12,47,34},177/1),({14,12,46,35},206/1),({14,12,45,36},226/1),({14,12,44,37},222/1),({14,12,43,38},198/1),({14,12,42,39},147/1),({14,12,41,40},79/1),({13,13,54,27},2/1),({13,13,53,28},4/1),({13,13,52,29},10/1),({13,13,51,30},17/1),({13,13,50,31},30/1),({13,13,49,32},43/1),({13,13,48,33},61/1),({13,13,47,34},76/1),({13,13,46,35},91/1),({13,13,45,36},96/1),({13,13,44,37},96/1),({13,13,43,38},83/1),({13,13,42,39},63/1),({13,13,41,40},33/1)}, (12,2) => {}, (14,0) => {}, (14,1) => {({20,10,50,43},1/1),({20,10,49,44},1/1),({20,10,48,45},1/1),({19,11,54,39},1/1),({19,11,53,40},2/1),({19,11,52,41},4/1),({19,11,51,42},5/1),({19,11,50,43},6/1),({19,11,49,44},6/1),({19,11,48,45},5/1),({19,11,47,46},3/1),({18,12,56,37},2/1),({18,12,55,38},3/1),({18,12,54,39},7/1),({18,12,53,40},11/1),({18,12,52,41},16/1),({18,12,51,42},19/1),({18,12,50,43},22/1),({18,12,49,44},20/1),({18,12,48,45},16/1),({18,12,47,46},9/1),({17,13,58,35},1/1),({17,13,57,36},3/1),({17,13,56,37},6/1),({17,13,55,38},12/1),({17,13,54,39},19/1),({17,13,53,40},27/1),({17,13,52,41},35/1),({17,13,51,42},41/1),({17,13,50,43},42/1),({17,13,49,44},39/1),({17,13,48,45},30/1),({17,13,47,46},16/1),({16,14,58,35},2/1),({16,14,57,36},4/1),({16,14,56,37},9/1),({16,14,55,38},14/1),({16,14,54,39},24/1),({16,14,53,40},32/1),({16,14,52,41},41/1),({16,14,51,42},45/1),({16,14,50,43},48/1),({16,14,49,44},43/1),({16,14,48,45},33/1),({16,14,47,46},17/1),({15,15,59,34},1/1),({15,15,58,35},1/1),({15,15,57,36},3/1),({15,15,56,37},5/1),({15,15,55,38},9/1),({15,15,54,39},12/1),({15,15,53,40},17/1),({15,15,52,41},19/1),({15,15,51,42},22/1),({15,15,50,43},22/1),({15,15,49,44},20/1),({15,15,48,45},14/1),({15,15,47,46},8/1)}, (16,0) => {}, (14,2) => {}, (16,1) => {({20,14,57,48},1/1),({20,14,55,50},1/1),({20,14,54,51},1/1),({19,15,59,46},1/1),({19,15,58,47},1/1),({19,15,57,48},1/1),({19,15,56,49},2/1),({19,15,55,50},2/1),({19,15,54,51},1/1),({19,15,53,52},1/1),({18,16,61,44},1/1),({18,16,60,45},1/1),({18,16,59,46},2/1),({18,16,58,47},2/1),({18,16,57,48},3/1),({18,16,56,49},3/1),({18,16,55,50},3/1),({18,16,54,51},2/1),({18,16,53,52},1/1),({17,17,60,45},1/1),({17,17,58,47},1/1),({17,17,57,48},1/1),({17,17,56,49},1/1),({17,17,55,50},1/1),({17,17,54,51},1/1)}, (18,0) => {}, (16,2) => {}, (18,1) => {}, (18,2) => {({20,20,62,61},1/1)}, (1,0) => {({2,0,8,1},1/1),({2,0,7,2},1/1),({2,0,6,3},1/1)}, (1,1) => {({2,2,15,0},1/1)}, (3,0) => {({6,0,15,6},1/1),({6,0,13,8},1/1),({6,0,12,9},1/1),({5,1,17,4},1/1),({5,1,16,5},1/1),({5,1,15,6},1/1),({5,1,14,7},2/1),({5,1,13,8},2/1),({5,1,12,9},1/1),({5,1,11,10},1/1),({4,2,17,4},1/1),({4,2,16,5},1/1),({4,2,15,6},2/1),({4,2,14,7},2/1),({4,2,13,8},3/1),({4,2,12,9},2/1),({4,2,11,10},1/1),({3,3,18,3},1/1),({3,3,16,5},1/1),({3,3,15,6},1/1),({3,3,14,7},1/1),({3,3,13,8},1/1),({3,3,12,9},1/1)}, (3,1) => {({6,2,24,3},1/1),({6,2,23,4},1/1),({6,2,22,5},2/1),({6,2,21,6},2/1),({6,2,20,7},3/1),({6,2,19,8},2/1),({6,2,18,9},2/1),({6,2,17,10},1/1),({6,2,16,11},1/1),({5,3,25,2},1/1),({5,3,24,3},1/1),({5,3,23,4},2/1),({5,3,22,5},2/1),({5,3,21,6},3/1),({5,3,20,7},3/1),({5,3,19,8},3/1),({5,3,18,9},2/1),({5,3,17,10},2/1),({5,3,16,11},1/1),({5,3,15,12},1/1),({4,4,24,3},1/1),({4,4,23,4},1/1),({4,4,22,5},2/1),({4,4,21,6},2/1),({4,4,20,7},3/1),({4,4,19,8},2/1),({4,4,18,9},2/1),({4,4,17,10},1/1),({4,4,16,11},1/1)}, (3,2) => {}, (5,0) => {}, (5,1) => {({12,0,20,19},1/1),({11,1,25,14},1/1),({11,1,24,15},2/1),({11,1,23,16},2/1),({11,1,22,17},2/1),({11,1,21,18},2/1),({11,1,20,19},1/1),({10,2,29,10},1/1),({10,2,28,11},2/1),({10,2,27,12},4/1),({10,2,26,13},6/1),({10,2,25,14},9/1),({10,2,24,15},11/1),({10,2,23,16},12/1),({10,2,22,17},11/1),({10,2,21,18},9/1),({10,2,20,19},5/1),({9,3,32,7},1/1),({9,3,31,8},2/1),({9,3,30,9},4/1),({9,3,29,10},8/1),({9,3,28,11},13/1),({9,3,27,12},19/1),({9,3,26,13},27/1),({9,3,25,14},32/1),({9,3,24,15},36/1),({9,3,23,16},38/1),({9,3,22,17},34/1),({9,3,21,18},25/1),({9,3,20,19},14/1),({8,4,33,6},1/1),({8,4,32,7},2/1),({8,4,31,8},6/1),({8,4,30,9},11/1),({8,4,29,10},19/1),({8,4,28,11},29/1),({8,4,27,12},42/1),({8,4,26,13},52/1),({8,4,25,14},63/1),({8,4,24,15},68/1),({8,4,23,16},67/1),({8,4,22,17},59/1),({8,4,21,18},45/1),({8,4,20,19},23/1),({7,5,34,5},1/1),({7,5,33,6},2/1),({7,5,32,7},5/1),({7,5,31,8},9/1),({7,5,30,9},17/1),({7,5,29,10},26/1),({7,5,28,11},38/1),({7,5,27,12},50/1),({7,5,26,13},63/1),({7,5,25,14},72/1),({7,5,24,15},77/1),({7,5,23,16},74/1),({7,5,22,17},65/1),({7,5,21,18},48/1),({7,5,20,19},26/1),({6,6,33,6},1/1),({6,6,32,7},2/1),({6,6,31,8},4/1),({6,6,30,9},6/1),({6,6,29,10},12/1),({6,6,28,11},16/1),({6,6,27,12},22/1),({6,6,26,13},27/1),({6,6,25,14},32/1),({6,6,24,15},32/1),({6,6,23,16},33/1),({6,6,22,17},27/1),({6,6,21,18},20/1),({6,6,20,19},12/1)}};
end;