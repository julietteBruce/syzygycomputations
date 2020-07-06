A := QQ[t_0,t_1,t_2,t_3];
--tb stands for Total Betti numbers
tb0233 = new HashTable from {(6,1) => 20592/1, (7,0) => 0, (5,2) => 0, (7,1) => 21879/1, (8,0) => 0, (6,2) => 0, (8,1) => 17160/1, (9,0) => 0, (7,2) => 0, (9,1) => 9867/1, (10,0) => 0, (8,2) => 0, (10,1) => 4056/1, (11,0) => 0, (9,2) => 0, (12,0) => 0, (11,1) => 1131/1, (10,2) => 0, (13,0) => 0, (12,1) => 192/1, (11,2) => 0, (13,1) => 15/1, (12,2) => 0, (13,2) => 0, (14,1) => 0, (15,1) => 0, (0,0) => 3/1, (1,0) => 24/1, (1,1) => 27/1, (2,0) => 66/1, (2,1) => 312/1, (3,0) => 0, (3,1) => 2145/1, (4,0) => 0, (2,2) => 0, (4,1) => 6864/1, (5,0) => 0, (3,2) => 0, (5,1) => 14157/1, (6,0) => 0, (4,2) => 0};
--sb represents the betti numbers as sums of Schur functors
sb0233 = new HashTable from {(6,1) => {({17,4,14,9},1/1),({17,4,13,10},1/1),({17,4,12,11},1/1),({16,5,16,7},1/1),({16,5,15,8},2/1),({16,5,14,9},5/1),({16,5,13,10},5/1),({16,5,12,11},4/1),({15,6,17,6},1/1),({15,6,16,7},4/1),({15,6,15,8},9/1),({15,6,14,9},14/1),({15,6,13,10},15/1),({15,6,12,11},10/1),({14,7,18,5},1/1),({14,7,17,6},3/1),({14,7,16,7},10/1),({14,7,15,8},19/1),({14,7,14,9},28/1),({14,7,13,10},29/1),({14,7,12,11},19/1),({13,8,18,5},2/1),({13,8,17,6},7/1),({13,8,16,7},17/1),({13,8,15,8},30/1),({13,8,14,9},41/1),({13,8,13,10},41/1),({13,8,12,11},26/1),({12,9,18,5},2/1),({12,9,17,6},8/1),({12,9,16,7},19/1),({12,9,15,8},32/1),({12,9,14,9},43/1),({12,9,13,10},41/1),({12,9,12,11},26/1),({11,10,19,4},1/1),({11,10,18,5},2/1),({11,10,17,6},6/1),({11,10,16,7},13/1),({11,10,15,8},22/1),({11,10,14,9},28/1),({11,10,13,10},27/1),({11,10,12,11},17/1)}, (7,0) => {}, (5,2) => {}, (7,1) => {({19,5,14,12},1/1),({18,6,17,9},1/1),({18,6,16,10},2/1),({18,6,15,11},3/1),({18,6,14,12},3/1),({18,6,13,13},1/1),({17,7,18,8},1/1),({17,7,17,9},3/1),({17,7,16,10},8/1),({17,7,15,11},10/1),({17,7,14,12},11/1),({17,7,13,13},3/1),({16,8,19,7},1/1),({16,8,18,8},4/1),({16,8,17,9},11/1),({16,8,16,10},19/1),({16,8,15,11},25/1),({16,8,14,12},21/1),({16,8,13,13},9/1),({15,9,19,7},2/1),({15,9,18,8},9/1),({15,9,17,9},19/1),({15,9,16,10},34/1),({15,9,15,11},39/1),({15,9,14,12},36/1),({15,9,13,13},12/1),({14,10,20,6},1/1),({14,10,19,7},4/1),({14,10,18,8},12/1),({14,10,17,9},27/1),({14,10,16,10},41/1),({14,10,15,11},50/1),({14,10,14,12},40/1),({14,10,13,13},17/1),({13,11,20,6},1/1),({13,11,19,7},4/1),({13,11,18,8},12/1),({13,11,17,9},23/1),({13,11,16,10},37/1),({13,11,15,11},40/1),({13,11,14,12},36/1),({13,11,13,13},11/1),({12,12,19,7},2/1),({12,12,18,8},4/1),({12,12,17,9},10/1),({12,12,16,10},13/1),({12,12,15,11},18/1),({12,12,14,12},11/1),({12,12,13,13},7/1)}, (8,0) => {}, (6,2) => {}, (8,1) => {({20,7,17,12},1/1),({20,7,16,13},1/1),({20,7,15,14},1/1),({19,8,19,10},1/1),({19,8,18,11},2/1),({19,8,17,12},4/1),({19,8,16,13},5/1),({19,8,15,14},4/1),({18,9,19,10},3/1),({18,9,18,11},7/1),({18,9,17,12},13/1),({18,9,16,13},14/1),({18,9,15,14},10/1),({17,10,20,9},2/1),({17,10,19,10},7/1),({17,10,18,11},16/1),({17,10,17,12},25/1),({17,10,16,13},27/1),({17,10,15,14},18/1),({16,11,21,8},1/1),({16,11,20,9},4/1),({16,11,19,10},13/1),({16,11,18,11},25/1),({16,11,17,12},37/1),({16,11,16,13},38/1),({16,11,15,14},25/1),({15,12,21,8},1/1),({15,12,20,9},5/1),({15,12,19,10},14/1),({15,12,18,11},27/1),({15,12,17,12},38/1),({15,12,16,13},39/1),({15,12,15,14},25/1),({14,13,21,8},1/1),({14,13,20,9},4/1),({14,13,19,10},10/1),({14,13,18,11},18/1),({14,13,17,12},25/1),({14,13,16,13},25/1),({14,13,15,14},16/1)}, (9,0) => {}, (7,2) => {}, (9,1) => {({21,9,19,13},1/1),({21,9,18,14},1/1),({21,9,17,15},2/1),({20,10,20,12},1/1),({20,10,19,13},3/1),({20,10,18,14},5/1),({20,10,17,15},5/1),({20,10,16,16},2/1),({19,11,21,11},1/1),({19,11,20,12},4/1),({19,11,19,13},9/1),({19,11,18,14},12/1),({19,11,17,15},13/1),({19,11,16,16},4/1),({18,12,21,11},2/1),({18,12,20,12},8/1),({18,12,19,13},16/1),({18,12,18,14},22/1),({18,12,17,15},19/1),({18,12,16,16},8/1),({17,13,22,10},1/1),({17,13,21,11},5/1),({17,13,20,12},11/1),({17,13,19,13},22/1),({17,13,18,14},26/1),({17,13,17,15},25/1),({17,13,16,16},8/1),({16,14,22,10},1/1),({16,14,21,11},4/1),({16,14,20,12},11/1),({16,14,19,13},18/1),({16,14,18,14},24/1),({16,14,17,15},19/1),({16,14,16,16},9/1),({15,15,21,11},2/1),({15,15,20,12},4/1),({15,15,19,13},8/1),({15,15,18,14},8/1),({15,15,17,15},9/1),({15,15,16,16},2/1)}, (10,0) => {}, (8,2) => {}, (10,1) => {({22,11,20,15},1/1),({22,11,19,16},1/1),({22,11,18,17},1/1),({21,12,21,14},1/1),({21,12,20,15},3/1),({21,12,19,16},4/1),({21,12,18,17},3/1),({20,13,22,13},1/1),({20,13,21,14},3/1),({20,13,20,15},7/1),({20,13,19,16},8/1),({20,13,18,17},6/1),({19,14,22,13},2/1),({19,14,21,14},6/1),({19,14,20,15},11/1),({19,14,19,16},13/1),({19,14,18,17},9/1),({18,15,22,13},3/1),({18,15,21,14},7/1),({18,15,20,15},12/1),({18,15,19,16},13/1),({18,15,18,17},9/1),({17,16,23,12},1/1),({17,16,22,13},2/1),({17,16,21,14},5/1),({17,16,20,15},8/1),({17,16,19,16},9/1),({17,16,18,17},6/1)}, (11,0) => {}, (9,2) => {}, (12,0) => {}, (11,1) => {({23,13,20,18},1/1),({22,14,22,16},1/1),({22,14,21,17},2/1),({22,14,20,18},2/1),({22,14,19,19},1/1),({21,15,22,16},2/1),({21,15,21,17},3/1),({21,15,20,18},4/1),({21,15,19,19},1/1),({20,16,23,15},1/1),({20,16,22,16},3/1),({20,16,21,17},5/1),({20,16,20,18},5/1),({20,16,19,19},2/1),({19,17,23,15},1/1),({19,17,22,16},3/1),({19,17,21,17},4/1),({19,17,20,18},5/1),({19,17,19,19},1/1),({18,18,23,15},1/1),({18,18,22,16},1/1),({18,18,21,17},2/1),({18,18,20,18},1/1),({18,18,19,19},1/1)}, (10,2) => {}, (13,0) => {}, (12,1) => {({23,16,22,19},1/1),({23,16,21,20},1/1),({22,17,23,18},1/1),({22,17,22,19},1/1),({22,17,21,20},1/1),({21,18,23,18},1/1),({21,18,22,19},1/1),({21,18,21,20},1/1),({20,19,23,18},1/1),({20,19,22,19},1/1),({20,19,21,20},1/1)}, (11,2) => {}, (13,1) => {({23,19,23,21},1/1)}, (12,2) => {}, (13,2) => {}, (14,1) => {}, (15,1) => {}, (0,0) => {({0,0,2,0},1/1)}, (1,0) => {({3,0,4,1},1/1),({3,0,3,2},1/1)}, (1,1) => {({4,2,8,0},1/1)}, (2,0) => {({6,0,5,3},1/1),({5,1,6,2},1/1),({5,1,4,4},1/1),({4,2,5,3},1/1),({3,3,6,2},1/1),({3,3,4,4},1/1)}, (2,1) => {({7,2,10,1},1/1),({7,2,9,2},1/1),({7,2,8,3},1/1),({6,3,10,1},1/1),({6,3,9,2},1/1),({6,3,8,3},1/1),({5,4,11,0},1/1),({5,4,10,1},1/1),({5,4,9,2},1/1),({5,4,8,3},1/1)}, (3,0) => {}, (3,1) => {({11,1,8,6},1/1),({10,2,11,3},1/1),({10,2,10,4},1/1),({10,2,9,5},2/1),({10,2,8,6},1/1),({10,2,7,7},1/1),({9,3,12,2},1/1),({9,3,11,3},2/1),({9,3,10,4},4/1),({9,3,9,5},3/1),({9,3,8,6},4/1),({8,4,13,1},1/1),({8,4,12,2},2/1),({8,4,11,3},4/1),({8,4,10,4},5/1),({8,4,9,5},6/1),({8,4,8,6},3/1),({8,4,7,7},2/1),({7,5,13,1},1/1),({7,5,12,2},3/1),({7,5,11,3},4/1),({7,5,10,4},6/1),({7,5,9,5},4/1),({7,5,8,6},5/1),({6,6,11,3},2/1),({6,6,10,4},1/1),({6,6,9,5},3/1),({6,6,7,7},2/1)}, (4,0) => {}, (2,2) => {}, (4,1) => {({13,2,11,6},1/1),({13,2,10,7},1/1),({13,2,9,8},1/1),({12,3,13,4},1/1),({12,3,12,5},2/1),({12,3,11,6},4/1),({12,3,10,7},4/1),({12,3,9,8},3/1),({11,4,14,3},1/1),({11,4,13,4},3/1),({11,4,12,5},6/1),({11,4,11,6},9/1),({11,4,10,7},9/1),({11,4,9,8},6/1),({10,5,15,2},1/1),({10,5,14,3},3/1),({10,5,13,4},7/1),({10,5,12,5},11/1),({10,5,11,6},15/1),({10,5,10,7},14/1),({10,5,9,8},9/1),({9,6,15,2},1/1),({9,6,14,3},3/1),({9,6,13,4},8/1),({9,6,12,5},13/1),({9,6,11,6},16/1),({9,6,10,7},15/1),({9,6,9,8},9/1),({8,7,15,2},1/1),({8,7,14,3},3/1),({8,7,13,4},6/1),({8,7,12,5},9/1),({8,7,11,6},11/1),({8,7,10,7},10/1),({8,7,9,8},6/1)}, (5,0) => {}, (3,2) => {}, (5,1) => {({15,3,13,7},1/1),({15,3,12,8},1/1),({15,3,11,9},2/1),({14,4,14,6},2/1),({14,4,13,7},4/1),({14,4,12,8},6/1),({14,4,11,9},5/1),({14,4,10,10},2/1),({13,5,16,4},1/1),({13,5,15,5},3/1),({13,5,14,6},6/1),({13,5,13,7},12/1),({13,5,12,8},14/1),({13,5,11,9},14/1),({13,5,10,10},4/1),({12,6,16,4},2/1),({12,6,15,5},6/1),({12,6,14,6},14/1),({12,6,13,7},21/1),({12,6,12,8},26/1),({12,6,11,9},20/1),({12,6,10,10},9/1),({11,7,17,3},1/1),({11,7,16,4},3/1),({11,7,15,5},10/1),({11,7,14,6},18/1),({11,7,13,7},29/1),({11,7,12,8},30/1),({11,7,11,9},28/1),({11,7,10,10},8/1),({10,8,17,3},1/1),({10,8,16,4},4/1),({10,8,15,5},9/1),({10,8,14,6},18/1),({10,8,13,7},24/1),({10,8,12,8},29/1),({10,8,11,9},20/1),({10,8,10,10},10/1),({9,9,16,4},1/1),({9,9,15,5},4/1),({9,9,14,6},6/1),({9,9,13,7},11/1),({9,9,12,8},9/1),({9,9,11,9},10/1),({9,9,10,10},1/1)}, (6,0) => {}, (4,2) => {}};
--dw stands for dominant weights
dw0233 = new HashTable from {(5,2) => {}, (7,0) => {}, (6,1) => {({17,4,14,9},1/1)}, (7,1) => {({19,5,14,12},1/1)}, (8,0) => {}, (6,2) => {}, (7,2) => {}, (9,0) => {}, (8,1) => {({20,7,17,12},1/1)}, (8,2) => {}, (10,0) => {}, (9,1) => {({21,9,19,13},1/1)}, (10,1) => {({22,11,20,15},1/1)}, (11,0) => {}, (9,2) => {}, (12,0) => {}, (11,1) => {({23,13,20,18},1/1)}, (10,2) => {}, (13,0) => {}, (12,1) => {({23,16,22,19},1/1)}, (11,2) => {}, (13,1) => {({23,19,23,21},1/1)}, (12,2) => {}, (14,1) => {}, (13,2) => {}, (15,1) => {}, (0,0) => {({0,0,2,0},1/1)}, (1,0) => {({3,0,4,1},1/1)}, (2,0) => {({6,0,5,3},1/1)}, (1,1) => {({4,2,8,0},1/1)}, (3,0) => {}, (2,1) => {({7,2,10,1},1/1)}, (2,2) => {}, (4,0) => {}, (3,1) => {({11,1,8,6},1/1)}, (3,2) => {}, (5,0) => {}, (4,1) => {({13,2,11,6},1/1)}, (4,2) => {}, (6,0) => {}, (5,1) => {({15,3,13,7},1/1)}};
end;