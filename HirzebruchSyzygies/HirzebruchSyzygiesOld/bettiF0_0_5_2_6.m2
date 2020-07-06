A := QQ[t_0,t_1,t_2,t_3];
--tb stands for Total Betti numbers
tb0526 = new HashTable from {(5,2) => 0, (7,0) => 0, (7,1) => 310284/1, (7,2) => 0, (9,0) => 0, (9,1) => 612612/1, (9,2) => 0, (11,0) => 0, (11,1) => 461448/1, (11,2) => 0, (13,0) => 0, (13,1) => 135864/1, (15,0) => 0, (13,2) => 0, (15,1) => 13770/1, (17,0) => 0, (15,2) => 0, (17,1) => 318/1, (17,2) => 0, (19,1) => 0, (0,0) => 6/1, (2,0) => 612/1, (2,1) => 306/1, (2,2) => 0, (4,0) => 6120/1, (4,1) => 12240/1, (4,2) => 0, (6,0) => 0, (6,1) => 143208/1, (6,2) => 0, (8,0) => 0, (8,1) => 495924/1, (8,2) => 0, (10,0) => 0, (10,1) => 596700/1, (10,2) => 0, (12,0) => 0, (12,1) => 282744/1, (12,2) => 0, (14,0) => 0, (14,1) => 50184/1, (16,0) => 0, (14,2) => 0, (16,1) => 2646/1, (18,0) => 0, (16,2) => 0, (18,1) => 18/1, (18,2) => 0, (20,1) => 0, (1,0) => 90/1, (1,1) => 18/1, (3,0) => 2448/1, (3,1) => 2448/1, (3,2) => 0, (5,0) => 8568/1, (5,1) => 42840/1};
--sb represents the betti numbers as sums of Schur functors
sb0526 = new HashTable from {(5,2) => {}, (7,0) => {}, (7,1) => {({14,2,32,21},1/1),({14,2,31,22},1/1),({14,2,30,23},1/1),({14,2,29,24},1/1),({14,2,28,25},1/1),({14,2,27,26},1/1),({13,3,37,16},1/1),({13,3,36,17},2/1),({13,3,35,18},3/1),({13,3,34,19},5/1),({13,3,33,20},7/1),({13,3,32,21},9/1),({13,3,31,22},11/1),({13,3,30,23},11/1),({13,3,29,24},10/1),({13,3,28,25},8/1),({13,3,27,26},5/1),({12,4,40,13},1/1),({12,4,39,14},2/1),({12,4,38,15},5/1),({12,4,37,16},8/1),({12,4,36,17},15/1),({12,4,35,18},22/1),({12,4,34,19},31/1),({12,4,33,20},38/1),({12,4,32,21},46/1),({12,4,31,22},49/1),({12,4,30,23},49/1),({12,4,29,24},42/1),({12,4,28,25},32/1),({12,4,27,26},17/1),({11,5,41,12},2/1),({11,5,40,13},4/1),({11,5,39,14},10/1),({11,5,38,15},19/1),({11,5,37,16},32/1),({11,5,36,17},48/1),({11,5,35,18},69/1),({11,5,34,19},89/1),({11,5,33,20},109/1),({11,5,32,21},122/1),({11,5,31,22},128/1),({11,5,30,23},123/1),({11,5,29,24},107/1),({11,5,28,25},78/1),({11,5,27,26},41/1),({10,6,43,10},1/1),({10,6,42,11},3/1),({10,6,41,12},6/1),({10,6,40,13},14/1),({10,6,39,14},25/1),({10,6,38,15},43/1),({10,6,37,16},66/1),({10,6,36,17},97/1),({10,6,35,18},128/1),({10,6,34,19},164/1),({10,6,33,20},191/1),({10,6,32,21},212/1),({10,6,31,22},216/1),({10,6,30,23},206/1),({10,6,29,24},174/1),({10,6,28,25},128/1),({10,6,27,26},67/1),({9,7,43,10},1/1),({9,7,42,11},3/1),({9,7,41,12},8/1),({9,7,40,13},16/1),({9,7,39,14},30/1),({9,7,38,15},48/1),({9,7,37,16},75/1),({9,7,36,17},104/1),({9,7,35,18},139/1),({9,7,34,19},170/1),({9,7,33,20},200/1),({9,7,32,21},217/1),({9,7,31,22},222/1),({9,7,30,23},207/1),({9,7,29,24},177/1),({9,7,28,25},128/1),({9,7,27,26},68/1),({8,8,44,9},1/1),({8,8,43,10},1/1),({8,8,42,11},3/1),({8,8,41,12},5/1),({8,8,40,13},10/1),({8,8,39,14},15/1),({8,8,38,15},26/1),({8,8,37,16},35/1),({8,8,36,17},51/1),({8,8,35,18},64/1),({8,8,34,19},80/1),({8,8,33,20},89/1),({8,8,32,21},100/1),({8,8,31,22},98/1),({8,8,30,23},93/1),({8,8,29,24},77/1),({8,8,28,25},57/1),({8,8,27,26},30/1)}, (7,2) => {}, (9,0) => {}, (9,1) => {({16,4,41,24},1/1),({16,4,40,25},1/1),({16,4,39,26},2/1),({16,4,38,27},3/1),({16,4,37,28},4/1),({16,4,36,29},4/1),({16,4,35,30},4/1),({16,4,34,31},3/1),({16,4,33,32},2/1),({15,5,44,21},1/1),({15,5,43,22},2/1),({15,5,42,23},5/1),({15,5,41,24},9/1),({15,5,40,25},14/1),({15,5,39,26},20/1),({15,5,38,27},26/1),({15,5,37,28},30/1),({15,5,36,29},31/1),({15,5,35,30},28/1),({15,5,34,31},22/1),({15,5,33,32},12/1),({14,6,46,19},1/1),({14,6,45,20},4/1),({14,6,44,21},8/1),({14,6,43,22},17/1),({14,6,42,23},29/1),({14,6,41,24},47/1),({14,6,40,25},66/1),({14,6,39,26},87/1),({14,6,38,27},104/1),({14,6,37,28},116/1),({14,6,36,29},115/1),({14,6,35,30},103/1),({14,6,34,31},77/1),({14,6,33,32},42/1),({13,7,48,17},1/1),({13,7,47,18},3/1),({13,7,46,19},8/1),({13,7,45,20},17/1),({13,7,44,21},34/1),({13,7,43,22},58/1),({13,7,42,23},92/1),({13,7,41,24},133/1),({13,7,40,25},180/1),({13,7,39,26},224/1),({13,7,38,27},260/1),({13,7,37,28},277/1),({13,7,36,29},272/1),({13,7,35,30},238/1),({13,7,34,31},177/1),({13,7,33,32},94/1),({12,8,49,16},1/1),({12,8,48,17},3/1),({12,8,47,18},9/1),({12,8,46,19},19/1),({12,8,45,20},39/1),({12,8,44,21},68/1),({12,8,43,22},112/1),({12,8,42,23},166/1),({12,8,41,24},234/1),({12,8,40,25},302/1),({12,8,39,26},369/1),({12,8,38,27},416/1),({12,8,37,28},439/1),({12,8,36,29},422/1),({12,8,35,30},367/1),({12,8,34,31},270/1),({12,8,33,32},144/1),({11,9,50,15},1/1),({11,9,49,16},2/1),({11,9,48,17},6/1),({11,9,47,18},13/1),({11,9,46,19},27/1),({11,9,45,20},48/1),({11,9,44,21},82/1),({11,9,43,22},126/1),({11,9,42,23},185/1),({11,9,41,24},249/1),({11,9,40,25},319/1),({11,9,39,26},379/1),({11,9,38,27},425/1),({11,9,37,28},440/1),({11,9,36,29},422/1),({11,9,35,30},362/1),({11,9,34,31},267/1),({11,9,33,32},141/1),({10,10,49,16},1/1),({10,10,48,17},2/1),({10,10,47,18},6/1),({10,10,46,19},11/1),({10,10,45,20},22/1),({10,10,44,21},35/1),({10,10,43,22},56/1),({10,10,42,23},79/1),({10,10,41,24},109/1),({10,10,40,25},135/1),({10,10,39,26},163/1),({10,10,38,27},179/1),({10,10,37,28},188/1),({10,10,36,29},177/1),({10,10,35,30},153/1),({10,10,34,31},112/1),({10,10,33,32},60/1)}, (9,2) => {}, (11,0) => {}, (11,1) => {({18,6,46,31},1/1),({18,6,45,32},1/1),({18,6,44,33},2/1),({18,6,43,34},2/1),({18,6,42,35},3/1),({18,6,41,36},3/1),({18,6,40,37},2/1),({18,6,39,38},1/1),({17,7,49,28},1/1),({17,7,48,29},2/1),({17,7,47,30},5/1),({17,7,46,31},9/1),({17,7,45,32},14/1),({17,7,44,33},18/1),({17,7,43,34},22/1),({17,7,42,35},24/1),({17,7,41,36},22/1),({17,7,40,37},17/1),({17,7,39,38},9/1),({16,8,51,26},1/1),({16,8,50,27},4/1),({16,8,49,28},9/1),({16,8,48,29},18/1),({16,8,47,30},30/1),({16,8,46,31},47/1),({16,8,45,32},64/1),({16,8,44,33},81/1),({16,8,43,34},91/1),({16,8,42,35},94/1),({16,8,41,36},85/1),({16,8,40,37},65/1),({16,8,39,38},35/1),({15,9,53,24},1/1),({15,9,52,25},3/1),({15,9,51,26},8/1),({15,9,50,27},18/1),({15,9,49,28},35/1),({15,9,48,29},60/1),({15,9,47,30},94/1),({15,9,46,31},132/1),({15,9,45,32},173/1),({15,9,44,33},207/1),({15,9,43,34},228/1),({15,9,42,35},227/1),({15,9,41,36},202/1),({15,9,40,37},152/1),({15,9,39,38},82/1),({14,10,54,23},1/1),({14,10,53,24},3/1),({14,10,52,25},10/1),({14,10,51,26},21/1),({14,10,50,27},42/1),({14,10,49,28},73/1),({14,10,48,29},118/1),({14,10,47,30},171/1),({14,10,46,31},234/1),({14,10,45,32},293/1),({14,10,44,33},343/1),({14,10,43,34},368/1),({14,10,42,35},363/1),({14,10,41,36},318/1),({14,10,40,37},238/1),({14,10,39,38},127/1),({13,11,55,22},1/1),({13,11,54,23},2/1),({13,11,53,24},6/1),({13,11,52,25},13/1),({13,11,51,26},28/1),({13,11,50,27},50/1),({13,11,49,28},85/1),({13,11,48,29},130/1),({13,11,47,30},187/1),({13,11,46,31},247/1),({13,11,45,32},306/1),({13,11,44,33},350/1),({13,11,43,34},373/1),({13,11,42,35},362/1),({13,11,41,36},317/1),({13,11,40,37},234/1),({13,11,39,38},125/1),({12,12,54,23},1/1),({12,12,53,24},3/1),({12,12,52,25},7/1),({12,12,51,26},13/1),({12,12,50,27},25/1),({12,12,49,28},39/1),({12,12,48,29},60/1),({12,12,47,30},83/1),({12,12,46,31},110/1),({12,12,45,32},133/1),({12,12,44,33},153/1),({12,12,43,34},160/1),({12,12,42,35},157/1),({12,12,41,36},136/1),({12,12,40,37},101/1),({12,12,39,38},53/1)}, (11,2) => {}, (13,0) => {}, (13,1) => {({20,8,47,42},1/1),({19,9,52,37},1/1),({19,9,51,38},2/1),({19,9,50,39},3/1),({19,9,49,40},4/1),({19,9,48,41},5/1),({19,9,47,42},5/1),({19,9,46,43},4/1),({19,9,45,44},2/1),({18,10,55,34},1/1),({18,10,54,35},2/1),({18,10,53,36},5/1),({18,10,52,37},9/1),({18,10,51,38},15/1),({18,10,50,39},20/1),({18,10,49,40},25/1),({18,10,48,41},27/1),({18,10,47,42},26/1),({18,10,46,43},20/1),({18,10,45,44},11/1),({17,11,56,33},2/1),({17,11,55,34},5/1),({17,11,54,35},12/1),({17,11,53,36},22/1),({17,11,52,37},35/1),({17,11,51,38},50/1),({17,11,50,39},65/1),({17,11,49,40},75/1),({17,11,48,41},78/1),({17,11,47,42},71/1),({17,11,46,43},55/1),({17,11,45,44},30/1),({16,12,58,31},1/1),({16,12,57,32},3/1),({16,12,56,33},7/1),({16,12,55,34},16/1),({16,12,54,35},29/1),({16,12,53,36},48/1),({16,12,52,37},71/1),({16,12,51,38},97/1),({16,12,50,39},119/1),({16,12,49,40},135/1),({16,12,48,41},137/1),({16,12,47,42},124/1),({16,12,46,43},94/1),({16,12,45,44},51/1),({15,13,58,31},1/1),({15,13,57,32},4/1),({15,13,56,33},10/1),({15,13,55,34},20/1),({15,13,54,35},36/1),({15,13,53,36},57/1),({15,13,52,37},83/1),({15,13,51,38},108/1),({15,13,50,39},131/1),({15,13,49,40},144/1),({15,13,48,41},145/1),({15,13,47,42},129/1),({15,13,46,43},97/1),({15,13,45,44},52/1),({14,14,59,30},1/1),({14,14,58,31},1/1),({14,14,57,32},3/1),({14,14,56,33},6/1),({14,14,55,34},11/1),({14,14,54,35},17/1),({14,14,53,36},27/1),({14,14,52,37},37/1),({14,14,51,38},49/1),({14,14,50,39},58/1),({14,14,49,40},64/1),({14,14,48,41},64/1),({14,14,47,42},58/1),({14,14,46,43},43/1),({14,14,45,44},23/1)}, (15,0) => {}, (13,2) => {}, (15,1) => {({20,12,56,45},1/1),({20,12,55,46},1/1),({20,12,54,47},2/1),({20,12,53,48},2/1),({20,12,52,49},2/1),({20,12,51,50},1/1),({19,13,59,42},1/1),({19,13,58,43},2/1),({19,13,57,44},4/1),({19,13,56,45},6/1),({19,13,55,46},8/1),({19,13,54,47},9/1),({19,13,53,48},9/1),({19,13,52,49},7/1),({19,13,51,50},4/1),({18,14,60,41},1/1),({18,14,59,42},3/1),({18,14,58,43},6/1),({18,14,57,44},10/1),({18,14,56,45},15/1),({18,14,55,46},19/1),({18,14,54,47},21/1),({18,14,53,48},20/1),({18,14,52,49},16/1),({18,14,51,50},9/1),({17,15,61,40},1/1),({17,15,60,41},3/1),({17,15,59,42},6/1),({17,15,58,43},10/1),({17,15,57,44},15/1),({17,15,56,45},20/1),({17,15,55,46},24/1),({17,15,54,47},25/1),({17,15,53,48},23/1),({17,15,52,49},18/1),({17,15,51,50},10/1),({16,16,60,41},1/1),({16,16,59,42},2/1),({16,16,58,43},4/1),({16,16,57,44},6/1),({16,16,56,45},9/1),({16,16,55,46},11/1),({16,16,54,47},12/1),({16,16,53,48},11/1),({16,16,52,49},9/1),({16,16,51,50},5/1)}, (17,0) => {}, (15,2) => {}, (17,1) => {({20,16,61,52},1/1),({20,16,60,53},1/1),({20,16,59,54},1/1),({20,16,58,55},1/1),({20,16,57,56},1/1),({19,17,62,51},1/1),({19,17,61,52},1/1),({19,17,60,53},1/1),({19,17,59,54},1/1),({19,17,58,55},1/1),({19,17,57,56},1/1),({18,18,62,51},1/1),({18,18,61,52},1/1),({18,18,60,53},1/1),({18,18,59,54},1/1),({18,18,58,55},1/1),({18,18,57,56},1/1)}, (17,2) => {}, (19,1) => {}, (0,0) => {({0,0,5,0},1/1)}, (2,0) => {({4,0,14,3},1/1),({4,0,13,4},1/1),({4,0,12,5},2/1),({4,0,11,6},2/1),({4,0,10,7},2/1),({4,0,9,8},1/1),({3,1,15,2},1/1),({3,1,14,3},1/1),({3,1,13,4},2/1),({3,1,12,5},2/1),({3,1,11,6},2/1),({3,1,10,7},2/1),({3,1,9,8},1/1),({2,2,14,3},1/1),({2,2,13,4},1/1),({2,2,12,5},2/1),({2,2,11,6},2/1),({2,2,10,7},2/1),({2,2,9,8},1/1)}, (2,1) => {({4,2,22,1},1/1),({4,2,21,2},1/1),({4,2,20,3},1/1),({4,2,19,4},1/1),({4,2,18,5},1/1),({4,2,17,6},1/1)}, (2,2) => {}, (4,0) => {({8,0,19,10},1/1),({8,0,18,11},1/1),({8,0,17,12},1/1),({8,0,16,13},1/1),({8,0,15,14},1/1),({7,1,22,7},1/1),({7,1,21,8},2/1),({7,1,20,9},3/1),({7,1,19,10},4/1),({7,1,18,11},5/1),({7,1,17,12},5/1),({7,1,16,13},4/1),({7,1,15,14},2/1),({6,2,23,6},1/1),({6,2,22,7},2/1),({6,2,21,8},4/1),({6,2,20,9},6/1),({6,2,19,10},9/1),({6,2,18,11},10/1),({6,2,17,12},10/1),({6,2,16,13},8/1),({6,2,15,14},5/1),({5,3,24,5},1/1),({5,3,23,6},2/1),({5,3,22,7},4/1),({5,3,21,8},6/1),({5,3,20,9},9/1),({5,3,19,10},11/1),({5,3,18,11},12/1),({5,3,17,12},11/1),({5,3,16,13},9/1),({5,3,15,14},5/1),({4,4,23,6},1/1),({4,4,22,7},1/1),({4,4,21,8},2/1),({4,4,20,9},3/1),({4,4,19,10},5/1),({4,4,18,11},5/1),({4,4,17,12},5/1),({4,4,16,13},4/1),({4,4,15,14},3/1)}, (4,1) => {({8,2,29,6},1/1),({8,2,28,7},1/1),({8,2,27,8},2/1),({8,2,26,9},3/1),({8,2,25,10},3/1),({8,2,24,11},3/1),({8,2,23,12},3/1),({8,2,22,13},2/1),({8,2,21,14},1/1),({8,2,20,15},1/1),({7,3,31,4},1/1),({7,3,30,5},2/1),({7,3,29,6},3/1),({7,3,28,7},5/1),({7,3,27,8},7/1),({7,3,26,9},8/1),({7,3,25,10},9/1),({7,3,24,11},9/1),({7,3,23,12},8/1),({7,3,22,13},7/1),({7,3,21,14},5/1),({7,3,20,15},3/1),({7,3,19,16},2/1),({7,3,18,17},1/1),({6,4,31,4},1/1),({6,4,30,5},2/1),({6,4,29,6},4/1),({6,4,28,7},6/1),({6,4,27,8},9/1),({6,4,26,9},11/1),({6,4,25,10},12/1),({6,4,24,11},12/1),({6,4,23,12},11/1),({6,4,22,13},9/1),({6,4,21,14},6/1),({6,4,20,15},4/1),({6,4,19,16},2/1),({6,4,18,17},1/1),({5,5,32,3},1/1),({5,5,31,4},1/1),({5,5,30,5},2/1),({5,5,29,6},3/1),({5,5,28,7},4/1),({5,5,27,8},5/1),({5,5,26,9},6/1),({5,5,25,10},6/1),({5,5,24,11},6/1),({5,5,23,12},6/1),({5,5,22,13},5/1),({5,5,21,14},4/1),({5,5,20,15},3/1),({5,5,19,16},2/1),({5,5,18,17},1/1)}, (4,2) => {}, (6,0) => {}, (6,1) => {({13,1,26,21},1/1),({12,2,32,15},1/1),({12,2,31,16},1/1),({12,2,30,17},2/1),({12,2,29,18},2/1),({12,2,28,19},3/1),({12,2,27,20},3/1),({12,2,26,21},3/1),({12,2,25,22},2/1),({12,2,24,23},1/1),({11,3,36,11},1/1),({11,3,35,12},2/1),({11,3,34,13},4/1),({11,3,33,14},6/1),({11,3,32,15},9/1),({11,3,31,16},12/1),({11,3,30,17},15/1),({11,3,29,18},17/1),({11,3,28,19},18/1),({11,3,27,20},17/1),({11,3,26,21},15/1),({11,3,25,22},11/1),({11,3,24,23},6/1),({10,4,38,9},1/1),({10,4,37,10},2/1),({10,4,36,11},5/1),({10,4,35,12},9/1),({10,4,34,13},16/1),({10,4,33,14},23/1),({10,4,32,15},33/1),({10,4,31,16},40/1),({10,4,30,17},49/1),({10,4,29,18},52/1),({10,4,28,19},54/1),({10,4,27,20},49/1),({10,4,26,21},42/1),({10,4,25,22},30/1),({10,4,24,23},16/1),({9,5,39,8},1/1),({9,5,38,9},3/1),({9,5,37,10},7/1),({9,5,36,11},13/1),({9,5,35,12},23/1),({9,5,34,13},35/1),({9,5,33,14},50/1),({9,5,32,15},65/1),({9,5,31,16},80/1),({9,5,30,17},91/1),({9,5,29,18},98/1),({9,5,28,19},97/1),({9,5,27,20},90/1),({9,5,26,21},75/1),({9,5,25,22},54/1),({9,5,24,23},28/1),({8,6,40,7},1/1),({8,6,39,8},2/1),({8,6,38,9},5/1),({8,6,37,10},9/1),({8,6,36,11},18/1),({8,6,35,12},27/1),({8,6,34,13},42/1),({8,6,33,14},56/1),({8,6,32,15},75/1),({8,6,31,16},88/1),({8,6,30,17},101/1),({8,6,29,18},105/1),({8,6,28,19},106/1),({8,6,27,20},95/1),({8,6,26,21},80/1),({8,6,25,22},56/1),({8,6,24,23},30/1),({7,7,39,8},1/1),({7,7,38,9},2/1),({7,7,37,10},5/1),({7,7,36,11},8/1),({7,7,35,12},14/1),({7,7,34,13},19/1),({7,7,33,14},27/1),({7,7,32,15},33/1),({7,7,31,16},41/1),({7,7,30,17},44/1),({7,7,29,18},48/1),({7,7,28,19},46/1),({7,7,27,20},44/1),({7,7,26,21},35/1),({7,7,25,22},26/1),({7,7,24,23},13/1)}, (6,2) => {}, (8,0) => {}, (8,1) => {({15,3,37,22},1/1),({15,3,36,23},1/1),({15,3,35,24},2/1),({15,3,34,25},2/1),({15,3,33,26},3/1),({15,3,32,27},3/1),({15,3,31,28},2/1),({15,3,30,29},1/1),({14,4,41,18},1/1),({14,4,40,19},2/1),({14,4,39,20},4/1),({14,4,38,21},6/1),({14,4,37,22},11/1),({14,4,36,23},15/1),({14,4,35,24},19/1),({14,4,34,25},21/1),({14,4,33,26},23/1),({14,4,32,27},21/1),({14,4,31,28},16/1),({14,4,30,29},8/1),({13,5,43,16},1/1),({13,5,42,17},3/1),({13,5,41,18},7/1),({13,5,40,19},14/1),({13,5,39,20},24/1),({13,5,38,21},37/1),({13,5,37,22},52/1),({13,5,36,23},67/1),({13,5,35,24},80/1),({13,5,34,25},87/1),({13,5,33,26},87/1),({13,5,32,27},77/1),({13,5,31,28},58/1),({13,5,30,29},31/1),({12,6,45,14},1/1),({12,6,44,15},3/1),({12,6,43,16},8/1),({12,6,42,17},15/1),({12,6,41,18},30/1),({12,6,40,19},49/1),({12,6,39,20},78/1),({12,6,38,21},109/1),({12,6,37,22},146/1),({12,6,36,23},178/1),({12,6,35,24},206/1),({12,6,34,25},216/1),({12,6,33,26},211/1),({12,6,32,27},182/1),({12,6,31,28},136/1),({12,6,30,29},72/1),({11,7,46,13},1/1),({11,7,45,14},3/1),({11,7,44,15},8/1),({11,7,43,16},18/1),({11,7,42,17},35/1),({11,7,41,18},60/1),({11,7,40,19},96/1),({11,7,39,20},141/1),({11,7,38,21},194/1),({11,7,37,22},248/1),({11,7,36,23},298/1),({11,7,35,24},333/1),({11,7,34,25},348/1),({11,7,33,26},332/1),({11,7,32,27},286/1),({11,7,31,28},210/1),({11,7,30,29},112/1),({10,8,47,12},1/1),({10,8,46,13},2/1),({10,8,45,14},6/1),({10,8,44,15},12/1),({10,8,43,16},25/1),({10,8,42,17},43/1),({10,8,41,18},73/1),({10,8,40,19},109/1),({10,8,39,20},158/1),({10,8,38,21},208/1),({10,8,37,22},265/1),({10,8,36,23},309/1),({10,8,35,24},344/1),({10,8,34,25},351/1),({10,8,33,26},336/1),({10,8,32,27},286/1),({10,8,31,28},210/1),({10,8,30,29},110/1),({9,9,46,13},1/1),({9,9,45,14},2/1),({9,9,44,15},6/1),({9,9,43,16},10/1),({9,9,42,17},20/1),({9,9,41,18},31/1),({9,9,40,19},49/1),({9,9,39,20},67/1),({9,9,38,21},92/1),({9,9,37,22},112/1),({9,9,36,23},134/1),({9,9,35,24},145/1),({9,9,34,25},151/1),({9,9,33,26},141/1),({9,9,32,27},122/1),({9,9,31,28},88/1),({9,9,30,29},47/1)}, (8,2) => {}, (10,0) => {}, (10,1) => {({17,5,44,27},1/1),({17,5,43,28},1/1),({17,5,42,29},2/1),({17,5,41,30},3/1),({17,5,40,31},4/1),({17,5,39,32},4/1),({17,5,38,33},4/1),({17,5,37,34},3/1),({17,5,36,35},2/1),({16,6,46,25},2/1),({16,6,45,26},4/1),({16,6,44,27},8/1),({16,6,43,28},13/1),({16,6,42,29},20/1),({16,6,41,30},25/1),({16,6,40,31},30/1),({16,6,39,32},31/1),({16,6,38,33},29/1),({16,6,37,34},22/1),({16,6,36,35},12/1),({15,7,49,22},1/1),({15,7,48,23},3/1),({15,7,47,24},7/1),({15,7,46,25},15/1),({15,7,45,26},27/1),({15,7,44,27},44/1),({15,7,43,28},64/1),({15,7,42,29},86/1),({15,7,41,30},105/1),({15,7,40,31},117/1),({15,7,39,32},118/1),({15,7,38,33},106/1),({15,7,37,34},80/1),({15,7,36,35},43/1),({14,8,50,21},2/1),({14,8,49,22},5/1),({14,8,48,23},14/1),({14,8,47,24},28/1),({14,8,46,25},52/1),({14,8,45,26},83/1),({14,8,44,27},127/1),({14,8,43,28},172/1),({14,8,42,29},221/1),({14,8,41,30},258/1),({14,8,40,31},281/1),({14,8,39,32},276/1),({14,8,38,33},245/1),({14,8,37,34},182/1),({14,8,36,35},98/1),({13,9,52,19},1/1),({13,9,51,20},3/1),({13,9,50,21},7/1),({13,9,49,22},17/1),({13,9,48,23},34/1),({13,9,47,24},62/1),({13,9,46,25},103/1),({13,9,45,26},159/1),({13,9,44,27},225/1),({13,9,43,28},299/1),({13,9,42,29},367/1),({13,9,41,30},421/1),({13,9,40,31},447/1),({13,9,39,32},435/1),({13,9,38,33},378/1),({13,9,37,34},281/1),({13,9,36,35},150/1),({12,10,52,19},1/1),({12,10,51,20},3/1),({12,10,50,21},10/1),({12,10,49,22},20/1),({12,10,48,23},41/1),({12,10,47,24},70/1),({12,10,46,25},116/1),({12,10,45,26},170/1),({12,10,44,27},239/1),({12,10,43,28},307/1),({12,10,42,29},376/1),({12,10,41,30},422/1),({12,10,40,31},446/1),({12,10,39,32},428/1),({12,10,38,33},373/1),({12,10,37,34},274/1),({12,10,36,35},146/1),({11,11,53,18},1/1),({11,11,52,19},1/1),({11,11,51,20},3/1),({11,11,50,21},6/1),({11,11,49,22},12/1),({11,11,48,23},20/1),({11,11,47,24},36/1),({11,11,46,25},53/1),({11,11,45,26},80/1),({11,11,44,27},107/1),({11,11,43,28},139/1),({11,11,42,29},164/1),({11,11,41,30},187/1),({11,11,40,31},192/1),({11,11,39,32},186/1),({11,11,38,33},159/1),({11,11,37,34},118/1),({11,11,36,35},62/1)}, (10,2) => {}, (12,0) => {}, (12,1) => {({19,7,47,36},1/1),({19,7,46,37},1/1),({19,7,45,38},1/1),({19,7,44,39},1/1),({19,7,43,40},1/1),({19,7,42,41},1/1),({18,8,51,32},1/1),({18,8,50,33},2/1),({18,8,49,34},4/1),({18,8,48,35},6/1),({18,8,47,36},10/1),({18,8,46,37},12/1),({18,8,45,38},13/1),({18,8,44,39},12/1),({18,8,43,40},10/1),({18,8,42,41},6/1),({17,9,53,30},1/1),({17,9,52,31},3/1),({17,9,51,32},8/1),({17,9,50,33},15/1),({17,9,49,34},25/1),({17,9,48,35},36/1),({17,9,47,36},47/1),({17,9,46,37},55/1),({17,9,45,38},58/1),({17,9,44,39},53/1),({17,9,43,40},41/1),({17,9,42,41},22/1),({16,10,55,28},1/1),({16,10,54,29},3/1),({16,10,53,30},9/1),({16,10,52,31},17/1),({16,10,51,32},33/1),({16,10,50,33},53/1),({16,10,49,34},80/1),({16,10,48,35},106/1),({16,10,47,36},132/1),({16,10,46,37},147/1),({16,10,45,38},151/1),({16,10,44,39},135/1),({16,10,43,40},103/1),({16,10,42,41},55/1),({15,11,56,27},1/1),({15,11,55,28},3/1),({15,11,54,29},9/1),({15,11,53,30},20/1),({15,11,52,31},39/1),({15,11,51,32},66/1),({15,11,50,33},103/1),({15,11,49,34},145/1),({15,11,48,35},189/1),({15,11,47,36},225/1),({15,11,46,37},247/1),({15,11,45,38},246/1),({15,11,44,39},219/1),({15,11,43,40},164/1),({15,11,42,41},88/1),({14,12,57,26},1/1),({14,12,56,27},2/1),({14,12,55,28},7/1),({14,12,54,29},14/1),({14,12,53,30},29/1),({14,12,52,31},49/1),({14,12,51,32},81/1),({14,12,50,33},117/1),({14,12,49,34},162/1),({14,12,48,35},202/1),({14,12,47,36},240/1),({14,12,46,37},257/1),({14,12,45,38},255/1),({14,12,44,39},223/1),({14,12,43,40},168/1),({14,12,42,41},90/1),({13,13,56,27},1/1),({13,13,55,28},2/1),({13,13,54,29},6/1),({13,13,53,30},11/1),({13,13,52,31},22/1),({13,13,51,32},34/1),({13,13,50,33},52/1),({13,13,49,34},69/1),({13,13,48,35},89/1),({13,13,47,36},102/1),({13,13,46,37},111/1),({13,13,45,38},107/1),({13,13,44,39},95/1),({13,13,43,40},70/1),({13,13,42,41},38/1)}, (12,2) => {}, (14,0) => {}, (14,1) => {({20,10,52,43},1/1),({20,10,51,44},1/1),({20,10,50,45},1/1),({20,10,49,46},1/1),({20,10,48,47},1/1),({19,11,56,39},1/1),({19,11,55,40},2/1),({19,11,54,41},4/1),({19,11,53,42},6/1),({19,11,52,43},8/1),({19,11,51,44},9/1),({19,11,50,45},9/1),({19,11,49,46},7/1),({19,11,48,47},4/1),({18,12,58,37},1/1),({18,12,57,38},2/1),({18,12,56,39},6/1),({18,12,55,40},10/1),({18,12,54,41},17/1),({18,12,53,42},23/1),({18,12,52,43},29/1),({18,12,51,44},31/1),({18,12,50,45},30/1),({18,12,49,46},23/1),({18,12,48,47},13/1),({17,13,59,36},1/1),({17,13,58,37},4/1),({17,13,57,38},9/1),({17,13,56,39},16/1),({17,13,55,40},27/1),({17,13,54,41},39/1),({17,13,53,42},50/1),({17,13,52,43},58/1),({17,13,51,44},61/1),({17,13,50,45},56/1),({17,13,49,46},43/1),({17,13,48,47},23/1),({16,14,60,35},1/1),({16,14,59,36},2/1),({16,14,58,37},6/1),({16,14,57,38},11/1),({16,14,56,39},21/1),({16,14,55,40},31/1),({16,14,54,41},45/1),({16,14,53,42},56/1),({16,14,52,43},66/1),({16,14,51,44},67/1),({16,14,50,45},62/1),({16,14,49,46},47/1),({16,14,48,47},26/1),({15,15,59,36},2/1),({15,15,58,37},3/1),({15,15,57,38},7/1),({15,15,56,39},11/1),({15,15,55,40},17/1),({15,15,54,41},22/1),({15,15,53,42},28/1),({15,15,52,43},30/1),({15,15,51,44},31/1),({15,15,50,45},27/1),({15,15,49,46},21/1),({15,15,48,47},11/1)}, (16,0) => {}, (14,2) => {}, (16,1) => {({20,14,59,48},1/1),({20,14,58,49},1/1),({20,14,57,50},2/1),({20,14,56,51},2/1),({20,14,55,52},2/1),({20,14,54,53},1/1),({19,15,61,46},1/1),({19,15,60,47},2/1),({19,15,59,48},3/1),({19,15,58,49},4/1),({19,15,57,50},5/1),({19,15,56,51},5/1),({19,15,55,52},4/1),({19,15,54,53},2/1),({18,16,61,46},2/1),({18,16,60,47},3/1),({18,16,59,48},5/1),({18,16,58,49},6/1),({18,16,57,50},8/1),({18,16,56,51},8/1),({18,16,55,52},6/1),({18,16,54,53},3/1),({17,17,62,45},1/1),({17,17,61,46},1/1),({17,17,60,47},2/1),({17,17,59,48},2/1),({17,17,58,49},3/1),({17,17,57,50},3/1),({17,17,56,51},3/1),({17,17,55,52},2/1),({17,17,54,53},1/1)}, (18,0) => {}, (16,2) => {}, (18,1) => {({20,18,62,57},1/1)}, (18,2) => {}, (20,1) => {}, (1,0) => {({2,0,10,1},1/1),({2,0,9,2},1/1),({2,0,8,3},1/1),({2,0,7,4},1/1),({2,0,6,5},1/1)}, (1,1) => {({2,2,17,0},1/1)}, (3,0) => {({6,0,17,6},1/1),({6,0,16,7},1/1),({6,0,15,8},2/1),({6,0,14,9},2/1),({6,0,13,10},2/1),({6,0,12,11},1/1),({5,1,19,4},1/1),({5,1,18,5},2/1),({5,1,17,6},3/1),({5,1,16,7},4/1),({5,1,15,8},5/1),({5,1,14,9},5/1),({5,1,13,10},4/1),({5,1,12,11},2/1),({4,2,19,4},1/1),({4,2,18,5},2/1),({4,2,17,6},4/1),({4,2,16,7},5/1),({4,2,15,8},7/1),({4,2,14,9},7/1),({4,2,13,10},6/1),({4,2,12,11},3/1),({3,3,20,3},1/1),({3,3,19,4},1/1),({3,3,18,5},2/1),({3,3,17,6},2/1),({3,3,16,7},3/1),({3,3,15,8},3/1),({3,3,14,9},3/1),({3,3,13,10},2/1),({3,3,12,11},1/1)}, (3,1) => {({6,2,26,3},1/1),({6,2,25,4},1/1),({6,2,24,5},2/1),({6,2,23,6},2/1),({6,2,22,7},3/1),({6,2,21,8},2/1),({6,2,20,9},2/1),({6,2,19,10},1/1),({6,2,18,11},1/1),({5,3,27,2},1/1),({5,3,26,3},1/1),({5,3,25,4},2/1),({5,3,24,5},2/1),({5,3,23,6},3/1),({5,3,22,7},3/1),({5,3,21,8},3/1),({5,3,20,9},2/1),({5,3,19,10},2/1),({5,3,18,11},1/1),({5,3,17,12},1/1),({4,4,26,3},1/1),({4,4,25,4},1/1),({4,4,24,5},2/1),({4,4,23,6},2/1),({4,4,22,7},3/1),({4,4,21,8},2/1),({4,4,20,9},2/1),({4,4,19,10},1/1),({4,4,18,11},1/1)}, (3,2) => {}, (5,0) => {({10,0,20,15},1/1),({9,1,24,11},1/1),({9,1,23,12},1/1),({9,1,22,13},2/1),({9,1,21,14},2/1),({9,1,20,15},2/1),({9,1,19,16},2/1),({9,1,18,17},1/1),({8,2,26,9},1/1),({8,2,25,10},1/1),({8,2,24,11},3/1),({8,2,23,12},4/1),({8,2,22,13},6/1),({8,2,21,14},6/1),({8,2,20,15},7/1),({8,2,19,16},5/1),({8,2,18,17},3/1),({7,3,27,8},1/1),({7,3,26,9},2/1),({7,3,25,10},4/1),({7,3,24,11},6/1),({7,3,23,12},9/1),({7,3,22,13},11/1),({7,3,21,14},12/1),({7,3,20,15},11/1),({7,3,19,16},9/1),({7,3,18,17},5/1),({6,4,28,7},1/1),({6,4,27,8},1/1),({6,4,26,9},3/1),({6,4,25,10},4/1),({6,4,24,11},8/1),({6,4,23,12},9/1),({6,4,22,13},12/1),({6,4,21,14},12/1),({6,4,20,15},13/1),({6,4,19,16},9/1),({6,4,18,17},5/1),({5,5,27,8},1/1),({5,5,26,9},1/1),({5,5,25,10},3/1),({5,5,24,11},3/1),({5,5,23,12},5/1),({5,5,22,13},5/1),({5,5,21,14},6/1),({5,5,20,15},5/1),({5,5,19,16},4/1),({5,5,18,17},2/1)}, (5,1) => {({10,2,31,10},1/1),({10,2,30,11},1/1),({10,2,29,12},2/1),({10,2,28,13},2/1),({10,2,27,14},3/1),({10,2,26,15},2/1),({10,2,25,16},2/1),({10,2,24,17},1/1),({10,2,23,18},1/1),({9,3,34,7},1/1),({9,3,33,8},2/1),({9,3,32,9},4/1),({9,3,31,10},6/1),({9,3,30,11},9/1),({9,3,29,12},11/1),({9,3,28,13},13/1),({9,3,27,14},13/1),({9,3,26,15},13/1),({9,3,25,16},11/1),({9,3,24,17},9/1),({9,3,23,18},6/1),({9,3,22,19},4/1),({9,3,21,20},2/1),({8,4,35,6},1/1),({8,4,34,7},2/1),({8,4,33,8},5/1),({8,4,32,9},8/1),({8,4,31,10},14/1),({8,4,30,11},18/1),({8,4,29,12},24/1),({8,4,28,13},26/1),({8,4,27,14},29/1),({8,4,26,15},26/1),({8,4,25,16},24/1),({8,4,24,17},18/1),({8,4,23,18},14/1),({8,4,22,19},8/1),({8,4,21,20},4/1),({7,5,36,5},1/1),({7,5,35,6},2/1),({7,5,34,7},5/1),({7,5,33,8},8/1),({7,5,32,9},14/1),({7,5,31,10},19/1),({7,5,30,11},26/1),({7,5,29,12},30/1),({7,5,28,13},35/1),({7,5,27,14},35/1),({7,5,26,15},35/1),({7,5,25,16},30/1),({7,5,24,17},26/1),({7,5,23,18},19/1),({7,5,22,19},13/1),({7,5,21,20},6/1),({6,6,35,6},1/1),({6,6,34,7},1/1),({6,6,33,8},3/1),({6,6,32,9},4/1),({6,6,31,10},8/1),({6,6,30,11},9/1),({6,6,29,12},13/1),({6,6,28,13},13/1),({6,6,27,14},16/1),({6,6,26,15},13/1),({6,6,25,16},13/1),({6,6,24,17},9/1),({6,6,23,18},8/1),({6,6,22,19},4/1),({6,6,21,20},2/1)}};
end;