A := QQ[t_0,t_1,t_2,t_3];
--tb stands for Total Betti numbers
tb1425 = new HashTable from {(5,2) => 0, (7,0) => 14520/1, (6,1) => 220/1, (7,1) => 4158/1, (8,0) => 4158/1, (6,2) => 0, (7,2) => 0, (9,0) => 220/1, (8,1) => 14520/1, (8,2) => 0, (9,1) => 20020/1, (10,0) => 0, (11,0) => 0, (10,1) => 16380/1, (9,2) => 0, (12,0) => 0, (11,1) => 9100/1, (10,2) => 0, (13,0) => 0, (12,1) => 3500/1, (11,2) => 0, (14,0) => 0, (13,1) => 900/1, (12,2) => 0, (13,2) => 0, (14,1) => 140/1, (15,0) => 0, (14,2) => 0, (15,1) => 10/1, (16,1) => 0, (15,2) => 0, (17,1) => 0, (0,0) => 10/1, (1,0) => 140/1, (2,0) => 900/1, (1,1) => 0, (3,0) => 3500/1, (2,1) => 0, (2,2) => 0, (4,0) => 9100/1, (3,1) => 0, (3,2) => 0, (5,0) => 16380/1, (4,1) => 0, (4,2) => 0, (6,0) => 20020/1, (5,1) => 0};
--sb represents the betti numbers as sums of Schur functors
sb1425 = new HashTable from {(5,2) => {}, (7,0) => {({12,3,26,13},1/1),({12,3,25,14},1/1),({12,3,24,15},2/1),({12,3,23,16},2/1),({12,3,22,17},2/1),({12,3,21,18},1/1),({12,3,20,19},1/1),({11,4,27,12},1/1),({11,4,26,13},2/1),({11,4,25,14},4/1),({11,4,24,15},6/1),({11,4,23,16},8/1),({11,4,22,17},7/1),({11,4,21,18},6/1),({11,4,20,19},3/1),({10,5,29,10},1/1),({10,5,28,11},2/1),({10,5,27,12},4/1),({10,5,26,13},7/1),({10,5,25,14},11/1),({10,5,24,15},15/1),({10,5,23,16},17/1),({10,5,22,17},18/1),({10,5,21,18},14/1),({10,5,20,19},8/1),({9,6,29,10},1/1),({9,6,28,11},3/1),({9,6,27,12},6/1),({9,6,26,13},11/1),({9,6,25,14},16/1),({9,6,24,15},20/1),({9,6,23,16},24/1),({9,6,22,17},23/1),({9,6,21,18},19/1),({9,6,20,19},11/1),({8,7,30,9},1/1),({8,7,29,10},1/1),({8,7,28,11},3/1),({8,7,27,12},6/1),({8,7,26,13},9/1),({8,7,25,14},13/1),({8,7,24,15},17/1),({8,7,23,16},18/1),({8,7,22,17},18/1),({8,7,21,18},15/1),({8,7,20,19},8/1)}, (6,1) => {({12,3,24,15},1/1),({11,4,23,16},1/1),({10,5,22,17},1/1),({9,6,21,18},1/1),({8,7,20,19},1/1)}, (7,1) => {({13,4,28,16},1/1),({13,4,27,17},1/1),({13,4,26,18},1/1),({13,4,25,19},1/1),({13,4,24,20},1/1),({12,5,28,16},1/1),({12,5,27,17},2/1),({12,5,26,18},3/1),({12,5,25,19},3/1),({12,5,24,20},2/1),({12,5,23,21},1/1),({11,6,31,13},1/1),({11,6,30,14},1/1),({11,6,29,15},2/1),({11,6,28,16},2/1),({11,6,27,17},4/1),({11,6,26,18},5/1),({11,6,25,19},6/1),({11,6,24,20},5/1),({11,6,23,21},4/1),({11,6,22,22},1/1),({10,7,30,14},1/1),({10,7,29,15},2/1),({10,7,28,16},3/1),({10,7,27,17},4/1),({10,7,26,18},5/1),({10,7,25,19},6/1),({10,7,24,20},6/1),({10,7,23,21},5/1),({10,7,22,22},2/1),({9,8,29,15},1/1),({9,8,28,16},2/1),({9,8,27,17},3/1),({9,8,26,18},3/1),({9,8,25,19},4/1),({9,8,24,20},4/1),({9,8,23,21},3/1),({9,8,22,22},1/1)}, (8,0) => {({13,4,28,16},1/1),({13,4,27,17},1/1),({13,4,26,18},1/1),({13,4,25,19},1/1),({13,4,24,20},1/1),({12,5,28,16},1/1),({12,5,27,17},2/1),({12,5,26,18},3/1),({12,5,25,19},3/1),({12,5,24,20},2/1),({12,5,23,21},1/1),({11,6,31,13},1/1),({11,6,30,14},1/1),({11,6,29,15},2/1),({11,6,28,16},2/1),({11,6,27,17},4/1),({11,6,26,18},5/1),({11,6,25,19},6/1),({11,6,24,20},5/1),({11,6,23,21},4/1),({11,6,22,22},1/1),({10,7,30,14},1/1),({10,7,29,15},2/1),({10,7,28,16},3/1),({10,7,27,17},4/1),({10,7,26,18},5/1),({10,7,25,19},6/1),({10,7,24,20},6/1),({10,7,23,21},5/1),({10,7,22,22},2/1),({9,8,29,15},1/1),({9,8,28,16},2/1),({9,8,27,17},3/1),({9,8,26,18},3/1),({9,8,25,19},4/1),({9,8,24,20},4/1),({9,8,23,21},3/1),({9,8,22,22},1/1)}, (6,2) => {}, (7,2) => {}, (9,0) => {({14,5,29,20},1/1),({13,6,28,21},1/1),({12,7,27,22},1/1),({11,8,26,23},1/1),({10,9,25,24},1/1)}, (8,1) => {({14,5,31,18},1/1),({14,5,30,19},1/1),({14,5,29,20},2/1),({14,5,28,21},2/1),({14,5,27,22},2/1),({14,5,26,23},1/1),({14,5,25,24},1/1),({13,6,32,17},1/1),({13,6,31,18},2/1),({13,6,30,19},4/1),({13,6,29,20},6/1),({13,6,28,21},8/1),({13,6,27,22},7/1),({13,6,26,23},6/1),({13,6,25,24},3/1),({12,7,34,15},1/1),({12,7,33,16},2/1),({12,7,32,17},4/1),({12,7,31,18},7/1),({12,7,30,19},11/1),({12,7,29,20},15/1),({12,7,28,21},17/1),({12,7,27,22},18/1),({12,7,26,23},14/1),({12,7,25,24},8/1),({11,8,34,15},1/1),({11,8,33,16},3/1),({11,8,32,17},6/1),({11,8,31,18},11/1),({11,8,30,19},16/1),({11,8,29,20},20/1),({11,8,28,21},24/1),({11,8,27,22},23/1),({11,8,26,23},19/1),({11,8,25,24},11/1),({10,9,35,14},1/1),({10,9,34,15},1/1),({10,9,33,16},3/1),({10,9,32,17},6/1),({10,9,31,18},9/1),({10,9,30,19},13/1),({10,9,29,20},17/1),({10,9,28,21},18/1),({10,9,27,22},18/1),({10,9,26,23},15/1),({10,9,25,24},8/1)}, (8,2) => {}, (9,1) => {({15,6,33,21},1/1),({15,6,32,22},1/1),({15,6,31,23},2/1),({15,6,30,24},2/1),({15,6,29,25},2/1),({15,6,28,26},1/1),({15,6,27,27},1/1),({14,7,35,19},1/1),({14,7,34,20},2/1),({14,7,33,21},4/1),({14,7,32,22},7/1),({14,7,31,23},9/1),({14,7,30,24},10/1),({14,7,29,25},10/1),({14,7,28,26},7/1),({14,7,27,27},2/1),({13,8,36,18},1/1),({13,8,35,19},3/1),({13,8,34,20},7/1),({13,8,33,21},12/1),({13,8,32,22},18/1),({13,8,31,23},23/1),({13,8,30,24},25/1),({13,8,29,25},23/1),({13,8,28,26},16/1),({13,8,27,27},6/1),({12,9,37,17},1/1),({12,9,36,18},3/1),({12,9,35,19},6/1),({12,9,34,20},12/1),({12,9,33,21},20/1),({12,9,32,22},27/1),({12,9,31,23},34/1),({12,9,30,24},37/1),({12,9,29,25},32/1),({12,9,28,26},23/1),({12,9,27,27},9/1),({11,10,37,17},1/1),({11,10,36,18},2/1),({11,10,35,19},6/1),({11,10,34,20},10/1),({11,10,33,21},16/1),({11,10,32,22},22/1),({11,10,31,23},27/1),({11,10,30,24},27/1),({11,10,29,25},26/1),({11,10,28,26},17/1),({11,10,27,27},6/1)}, (10,0) => {}, (11,0) => {}, (10,1) => {({16,7,34,25},1/1),({16,7,33,26},1/1),({16,7,32,27},1/1),({16,7,31,28},1/1),({16,7,30,29},1/1),({15,8,37,22},1/1),({15,8,36,23},2/1),({15,8,35,24},4/1),({15,8,34,25},6/1),({15,8,33,26},8/1),({15,8,32,27},8/1),({15,8,31,28},7/1),({15,8,30,29},4/1),({14,9,38,21},1/1),({14,9,37,22},3/1),({14,9,36,23},7/1),({14,9,35,24},12/1),({14,9,34,25},17/1),({14,9,33,26},21/1),({14,9,32,27},22/1),({14,9,31,28},18/1),({14,9,30,29},10/1),({13,10,39,20},1/1),({13,10,38,21},3/1),({13,10,37,22},7/1),({13,10,36,23},13/1),({13,10,35,24},21/1),({13,10,34,25},28/1),({13,10,33,26},33/1),({13,10,32,27},33/1),({13,10,31,28},27/1),({13,10,30,29},15/1),({12,11,39,20},1/1),({12,11,38,21},3/1),({12,11,37,22},6/1),({12,11,36,23},11/1),({12,11,35,24},17/1),({12,11,34,25},23/1),({12,11,33,26},26/1),({12,11,32,27},26/1),({12,11,31,28},21/1),({12,11,30,29},12/1)}, (9,2) => {}, (12,0) => {}, (11,1) => {({17,8,34,30},1/1),({16,9,38,26},1/1),({16,9,37,27},2/1),({16,9,36,28},3/1),({16,9,35,29},4/1),({16,9,34,30},4/1),({16,9,33,31},3/1),({16,9,32,32},1/1),({15,10,40,24},1/1),({15,10,39,25},2/1),({15,10,38,26},5/1),({15,10,37,27},8/1),({15,10,36,28},12/1),({15,10,35,29},13/1),({15,10,34,30},13/1),({15,10,33,31},9/1),({15,10,32,32},4/1),({14,11,40,24},2/1),({14,11,39,25},5/1),({14,11,38,26},9/1),({14,11,37,27},15/1),({14,11,36,28},20/1),({14,11,35,29},22/1),({14,11,34,30},21/1),({14,11,33,31},15/1),({14,11,32,32},5/1),({13,12,41,23},1/1),({13,12,40,24},2/1),({13,12,39,25},5/1),({13,12,38,26},9/1),({13,12,37,27},13/1),({13,12,36,28},17/1),({13,12,35,29},19/1),({13,12,34,30},17/1),({13,12,33,31},12/1),({13,12,32,32},5/1)}, (10,2) => {}, (13,0) => {}, (12,1) => {({17,10,38,31},1/1),({17,10,37,32},1/1),({17,10,36,33},1/1),({17,10,35,34},1/1),({16,11,41,28},1/1),({16,11,40,29},2/1),({16,11,39,30},4/1),({16,11,38,31},5/1),({16,11,37,32},6/1),({16,11,36,33},5/1),({16,11,35,34},3/1),({15,12,42,27},1/1),({15,12,41,28},2/1),({15,12,40,29},5/1),({15,12,39,30},8/1),({15,12,38,31},10/1),({15,12,37,32},11/1),({15,12,36,33},10/1),({15,12,35,34},5/1),({14,13,42,27},1/1),({14,13,41,28},3/1),({14,13,40,29},5/1),({14,13,39,30},7/1),({14,13,38,31},10/1),({14,13,37,32},10/1),({14,13,36,33},8/1),({14,13,35,34},5/1)}, (11,2) => {}, (14,0) => {}, (13,1) => {({17,12,41,33},1/1),({17,12,40,34},1/1),({17,12,39,35},2/1),({17,12,38,36},1/1),({17,12,37,37},1/1),({16,13,43,31},1/1),({16,13,42,32},2/1),({16,13,41,33},3/1),({16,13,40,34},4/1),({16,13,39,35},4/1),({16,13,38,36},3/1),({16,13,37,37},1/1),({15,14,43,31},1/1),({15,14,42,32},2/1),({15,14,41,33},3/1),({15,14,40,34},4/1),({15,14,39,35},4/1),({15,14,38,36},3/1),({15,14,37,37},1/1)}, (12,2) => {}, (13,2) => {}, (14,1) => {({17,14,43,36},1/1),({17,14,42,37},1/1),({17,14,41,38},1/1),({17,14,40,39},1/1),({16,15,44,35},1/1),({16,15,43,36},1/1),({16,15,42,37},1/1),({16,15,41,38},1/1),({16,15,40,39},1/1)}, (15,0) => {}, (14,2) => {}, (15,1) => {({17,16,44,40},1/1)}, (16,1) => {}, (15,2) => {}, (17,1) => {}, (0,0) => {({1,0,4,0},1/1)}, (1,0) => {({3,0,8,1},1/1),({3,0,7,2},1/1),({3,0,6,3},1/1),({3,0,5,4},1/1),({2,1,9,0},1/1),({2,1,8,1},1/1),({2,1,7,2},1/1),({2,1,6,3},1/1),({2,1,5,4},1/1)}, (2,0) => {({5,0,11,3},1/1),({5,0,10,4},1/1),({5,0,9,5},2/1),({5,0,8,6},1/1),({5,0,7,7},1/1),({4,1,13,1},1/1),({4,1,12,2},2/1),({4,1,11,3},3/1),({4,1,10,4},4/1),({4,1,9,5},4/1),({4,1,8,6},3/1),({4,1,7,7},1/1),({3,2,13,1},1/1),({3,2,12,2},2/1),({3,2,11,3},3/1),({3,2,10,4},4/1),({3,2,9,5},4/1),({3,2,8,6},3/1),({3,2,7,7},1/1)}, (1,1) => {}, (3,0) => {({7,0,13,6},1/1),({7,0,12,7},1/1),({7,0,11,8},1/1),({7,0,10,9},1/1),({6,1,16,3},1/1),({6,1,15,4},2/1),({6,1,14,5},4/1),({6,1,13,6},5/1),({6,1,12,7},6/1),({6,1,11,8},5/1),({6,1,10,9},3/1),({5,2,17,2},1/1),({5,2,16,3},2/1),({5,2,15,4},5/1),({5,2,14,5},8/1),({5,2,13,6},10/1),({5,2,12,7},11/1),({5,2,11,8},10/1),({5,2,10,9},5/1),({4,3,17,2},1/1),({4,3,16,3},3/1),({4,3,15,4},5/1),({4,3,14,5},7/1),({4,3,13,6},10/1),({4,3,12,7},10/1),({4,3,11,8},8/1),({4,3,10,9},5/1)}, (2,1) => {}, (2,2) => {}, (4,0) => {({9,0,14,10},1/1),({8,1,18,6},1/1),({8,1,17,7},2/1),({8,1,16,8},3/1),({8,1,15,9},4/1),({8,1,14,10},4/1),({8,1,13,11},3/1),({8,1,12,12},1/1),({7,2,20,4},1/1),({7,2,19,5},2/1),({7,2,18,6},5/1),({7,2,17,7},8/1),({7,2,16,8},12/1),({7,2,15,9},13/1),({7,2,14,10},13/1),({7,2,13,11},9/1),({7,2,12,12},4/1),({6,3,20,4},2/1),({6,3,19,5},5/1),({6,3,18,6},9/1),({6,3,17,7},15/1),({6,3,16,8},20/1),({6,3,15,9},22/1),({6,3,14,10},21/1),({6,3,13,11},15/1),({6,3,12,12},5/1),({5,4,21,3},1/1),({5,4,20,4},2/1),({5,4,19,5},5/1),({5,4,18,6},9/1),({5,4,17,7},13/1),({5,4,16,8},17/1),({5,4,15,9},19/1),({5,4,14,10},17/1),({5,4,13,11},12/1),({5,4,12,12},5/1)}, (3,1) => {}, (3,2) => {}, (5,0) => {({10,1,19,10},1/1),({10,1,18,11},1/1),({10,1,17,12},1/1),({10,1,16,13},1/1),({10,1,15,14},1/1),({9,2,22,7},1/1),({9,2,21,8},2/1),({9,2,20,9},4/1),({9,2,19,10},6/1),({9,2,18,11},8/1),({9,2,17,12},8/1),({9,2,16,13},7/1),({9,2,15,14},4/1),({8,3,23,6},1/1),({8,3,22,7},3/1),({8,3,21,8},7/1),({8,3,20,9},12/1),({8,3,19,10},17/1),({8,3,18,11},21/1),({8,3,17,12},22/1),({8,3,16,13},18/1),({8,3,15,14},10/1),({7,4,24,5},1/1),({7,4,23,6},3/1),({7,4,22,7},7/1),({7,4,21,8},13/1),({7,4,20,9},21/1),({7,4,19,10},28/1),({7,4,18,11},33/1),({7,4,17,12},33/1),({7,4,16,13},27/1),({7,4,15,14},15/1),({6,5,24,5},1/1),({6,5,23,6},3/1),({6,5,22,7},6/1),({6,5,21,8},11/1),({6,5,20,9},17/1),({6,5,19,10},23/1),({6,5,18,11},26/1),({6,5,17,12},26/1),({6,5,16,13},21/1),({6,5,15,14},12/1)}, (4,1) => {}, (4,2) => {}, (6,0) => {({11,2,23,11},1/1),({11,2,22,12},1/1),({11,2,21,13},2/1),({11,2,20,14},2/1),({11,2,19,15},2/1),({11,2,18,16},1/1),({11,2,17,17},1/1),({10,3,25,9},1/1),({10,3,24,10},2/1),({10,3,23,11},4/1),({10,3,22,12},7/1),({10,3,21,13},9/1),({10,3,20,14},10/1),({10,3,19,15},10/1),({10,3,18,16},7/1),({10,3,17,17},2/1),({9,4,26,8},1/1),({9,4,25,9},3/1),({9,4,24,10},7/1),({9,4,23,11},12/1),({9,4,22,12},18/1),({9,4,21,13},23/1),({9,4,20,14},25/1),({9,4,19,15},23/1),({9,4,18,16},16/1),({9,4,17,17},6/1),({8,5,27,7},1/1),({8,5,26,8},3/1),({8,5,25,9},6/1),({8,5,24,10},12/1),({8,5,23,11},20/1),({8,5,22,12},27/1),({8,5,21,13},34/1),({8,5,20,14},37/1),({8,5,19,15},32/1),({8,5,18,16},23/1),({8,5,17,17},9/1),({7,6,27,7},1/1),({7,6,26,8},2/1),({7,6,25,9},6/1),({7,6,24,10},10/1),({7,6,23,11},16/1),({7,6,22,12},22/1),({7,6,21,13},27/1),({7,6,20,14},27/1),({7,6,19,15},26/1),({7,6,18,16},17/1),({7,6,17,17},6/1)}, (5,1) => {}};
--dw stands for dominant weights
dw1425 = new HashTable from {(6,1) => {({12,3,24,15},1/1)}, (7,0) => {({12,3,26,13},1/1)}, (5,2) => {}, (7,1) => {({13,4,28,16},1/1)}, (8,0) => {({13,4,28,16},1/1)}, (6,2) => {}, (8,1) => {({14,5,31,18},1/1)}, (9,0) => {({14,5,29,20},1/1)}, (7,2) => {}, (10,0) => {}, (9,1) => {({15,6,33,21},1/1)}, (8,2) => {}, (11,0) => {}, (10,1) => {({16,7,34,25},1/1)}, (9,2) => {}, (12,0) => {}, (11,1) => {({17,8,34,30},1/1)}, (10,2) => {}, (13,0) => {}, (12,1) => {({17,10,38,31},1/1)}, (11,2) => {}, (14,0) => {}, (13,1) => {({17,12,41,33},1/1)}, (12,2) => {}, (15,0) => {}, (14,1) => {({17,14,43,36},1/1)}, (13,2) => {}, (15,1) => {({17,16,44,40},1/1)}, (14,2) => {}, (15,2) => {}, (16,1) => {}, (17,1) => {}, (0,0) => {({1,0,4,0},1/1)}, (1,0) => {({3,0,8,1},1/1)}, (1,1) => {}, (2,0) => {({5,0,11,3},1/1)}, (2,1) => {}, (3,0) => {({7,0,13,6},1/1)}, (3,1) => {}, (4,0) => {({9,0,14,10},1/1)}, (2,2) => {}, (4,1) => {}, (5,0) => {({10,1,19,10},1/1)}, (3,2) => {}, (5,1) => {}, (6,0) => {({11,2,23,11},1/1)}, (4,2) => {}};
--dmw stands for dominant monomial weights
dmw1425 = new HashTable from {(6,1) => {{12,3,20,18}}, (7,0) => {}, (5,2) => {}, (7,1) => {{13,4,25,18}}, (8,0) => {}, (6,2) => {}, (8,1) => {{13,6,30,18}}, (9,0) => {}, (7,2) => {}, (10,0) => {}, (9,1) => {{14,7,27,26}}, (8,2) => {}, (11,0) => {}, (10,1) => {{14,9,32,26}}, (9,2) => {}, (12,0) => {}, (11,1) => {{15,10,34,29}}, (10,2) => {}, (13,0) => {}, (12,1) => {{16,11,37,31}}, (11,2) => {}, (14,0) => {}, (13,1) => {{16,13,34,39}}, (12,2) => {}, (15,0) => {}, (14,1) => {{16,15,41,37}}, (13,2) => {}, (15,1) => {{17,16,44,39}}, (14,2) => {}, (15,2) => {{16,15,34,45}}, (16,1) => {}, (17,1) => {}, (0,0) => {}, (1,0) => {}, (1,1) => {}, (2,0) => {}, (2,1) => {}, (3,0) => {}, (3,1) => {}, (4,0) => {}, (2,2) => {}, (4,1) => {{11,0,14,14}}, (5,0) => {}, (3,2) => {}, (5,1) => {{12,1,16,17}}, (6,0) => {}, (4,2) => {}};
