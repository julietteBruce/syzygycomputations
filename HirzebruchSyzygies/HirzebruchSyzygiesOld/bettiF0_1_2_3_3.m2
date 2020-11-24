A := QQ[t_0,t_1,t_2,t_3];
--tb stands for Total Betti numbers
tb1233 = new HashTable from {(6,1) => 10296/1, (7,0) => 0, (5,2) => 0, (7,1) => 12870/1, (8,0) => 0, (6,2) => 0, (8,1) => 11154/1, (9,0) => 0, (7,2) => 0, (9,1) => 6864/1, (10,0) => 0, (8,2) => 0, (10,1) => 2964/1, (11,0) => 0, (9,2) => 0, (12,0) => 0, (11,1) => 858/1, (10,2) => 0, (13,0) => 0, (12,1) => 150/1, (11,2) => 0, (13,1) => 12/1, (12,2) => 0, (13,2) => 0, (14,1) => 0, (15,1) => 0, (0,0) => 6/1, (1,0) => 66/1, (1,1) => 0, (2,0) => 312/1, (2,1) => 12/1, (3,0) => 792/1, (3,1) => 132/1, (4,0) => 990/1, (2,2) => 0, (4,1) => 1110/1, (5,0) => 252/1, (3,2) => 0, (5,1) => 5148/1, (6,0) => 0, (4,2) => 0};
--sb represents the betti numbers as sums of Schur functors
sb1233 = new HashTable from {(6,1) => {({17,5,14,9},1/1),({17,5,13,10},1/1),({17,5,12,11},1/1),({16,6,16,7},1/1),({16,6,15,8},2/1),({16,6,14,9},4/1),({16,6,13,10},4/1),({16,6,12,11},3/1),({15,7,17,6},1/1),({15,7,16,7},3/1),({15,7,15,8},7/1),({15,7,14,9},10/1),({15,7,13,10},11/1),({15,7,12,11},7/1),({14,8,18,5},1/1),({14,8,17,6},2/1),({14,8,16,7},7/1),({14,8,15,8},12/1),({14,8,14,9},18/1),({14,8,13,10},18/1),({14,8,12,11},12/1),({13,9,18,5},1/1),({13,9,17,6},5/1),({13,9,16,7},10/1),({13,9,15,8},18/1),({13,9,14,9},23/1),({13,9,13,10},23/1),({13,9,12,11},14/1),({12,10,18,5},1/1),({12,10,17,6},3/1),({12,10,16,7},9/1),({12,10,15,8},14/1),({12,10,14,9},20/1),({12,10,13,10},18/1),({12,10,12,11},12/1),({11,11,19,4},1/1),({11,11,18,5},1/1),({11,11,17,6},3/1),({11,11,16,7},4/1),({11,11,15,8},8/1),({11,11,14,9},8/1),({11,11,13,10},9/1),({11,11,12,11},5/1)}, (7,0) => {}, (5,2) => {}, (7,1) => {({19,6,14,12},1/1),({18,7,17,9},1/1),({18,7,16,10},2/1),({18,7,15,11},3/1),({18,7,14,12},3/1),({18,7,13,13},1/1),({17,8,18,8},1/1),({17,8,17,9},3/1),({17,8,16,10},7/1),({17,8,15,11},9/1),({17,8,14,12},9/1),({17,8,13,13},3/1),({16,9,19,7},1/1),({16,9,18,8},4/1),({16,9,17,9},9/1),({16,9,16,10},16/1),({16,9,15,11},20/1),({16,9,14,12},17/1),({16,9,13,13},7/1),({15,10,19,7},1/1),({15,10,18,8},6/1),({15,10,17,9},14/1),({15,10,16,10},24/1),({15,10,15,11},28/1),({15,10,14,12},25/1),({15,10,13,13},9/1),({14,11,20,6},1/1),({14,11,19,7},3/1),({14,11,18,8},8/1),({14,11,17,9},17/1),({14,11,16,10},26/1),({14,11,15,11},31/1),({14,11,14,12},26/1),({14,11,13,13},10/1),({13,12,19,7},2/1),({13,12,18,8},5/1),({13,12,17,9},11/1),({13,12,16,10},17/1),({13,12,15,11},20/1),({13,12,14,12},16/1),({13,12,13,13},6/1)}, (8,0) => {}, (6,2) => {}, (8,1) => {({20,8,17,12},1/1),({20,8,16,13},1/1),({20,8,15,14},1/1),({19,9,19,10},1/1),({19,9,18,11},2/1),({19,9,17,12},4/1),({19,9,16,13},5/1),({19,9,15,14},4/1),({18,10,19,10},3/1),({18,10,18,11},7/1),({18,10,17,12},12/1),({18,10,16,13},13/1),({18,10,15,14},9/1),({17,11,20,9},2/1),({17,11,19,10},6/1),({17,11,18,11},14/1),({17,11,17,12},21/1),({17,11,16,13},23/1),({17,11,15,14},15/1),({16,12,21,8},1/1),({16,12,20,9},3/1),({16,12,19,10},10/1),({16,12,18,11},19/1),({16,12,17,12},29/1),({16,12,16,13},29/1),({16,12,15,14},19/1),({15,13,20,9},3/1),({15,13,19,10},8/1),({15,13,18,11},17/1),({15,13,17,12},23/1),({15,13,16,13},25/1),({15,13,15,14},16/1),({14,14,21,8},1/1),({14,14,20,9},2/1),({14,14,19,10},5/1),({14,14,18,11},7/1),({14,14,17,12},11/1),({14,14,16,13},10/1),({14,14,15,14},7/1)}, (9,0) => {}, (7,2) => {}, (9,1) => {({21,10,19,13},1/1),({21,10,18,14},1/1),({21,10,17,15},2/1),({20,11,20,12},1/1),({20,11,19,13},3/1),({20,11,18,14},5/1),({20,11,17,15},5/1),({20,11,16,16},2/1),({19,12,21,11},1/1),({19,12,20,12},4/1),({19,12,19,13},9/1),({19,12,18,14},12/1),({19,12,17,15},12/1),({19,12,16,16},4/1),({18,13,21,11},2/1),({18,13,20,12},7/1),({18,13,19,13},14/1),({18,13,18,14},19/1),({18,13,17,15},17/1),({18,13,16,16},7/1),({17,14,22,10},1/1),({17,14,21,11},4/1),({17,14,20,12},9/1),({17,14,19,13},17/1),({17,14,18,14},21/1),({17,14,17,15},19/1),({17,14,16,16},7/1),({16,15,21,11},2/1),({16,15,20,12},6/1),({16,15,19,13},11/1),({16,15,18,14},14/1),({16,15,17,15},12/1),({16,15,16,16},5/1)}, (10,0) => {}, (8,2) => {}, (10,1) => {({22,12,20,15},1/1),({22,12,19,16},1/1),({22,12,18,17},1/1),({21,13,21,14},1/1),({21,13,20,15},3/1),({21,13,19,16},4/1),({21,13,18,17},3/1),({20,14,22,13},1/1),({20,14,21,14},3/1),({20,14,20,15},7/1),({20,14,19,16},8/1),({20,14,18,17},6/1),({19,15,22,13},2/1),({19,15,21,14},6/1),({19,15,20,15},10/1),({19,15,19,16},12/1),({19,15,18,17},8/1),({18,16,22,13},2/1),({18,16,21,14},5/1),({18,16,20,15},9/1),({18,16,19,16},10/1),({18,16,18,17},7/1),({17,17,23,12},1/1),({17,17,22,13},1/1),({17,17,21,14},3/1),({17,17,20,15},4/1),({17,17,19,16},5/1),({17,17,18,17},3/1)}, (11,0) => {}, (9,2) => {}, (12,0) => {}, (11,1) => {({23,14,20,18},1/1),({22,15,22,16},1/1),({22,15,21,17},2/1),({22,15,20,18},2/1),({22,15,19,19},1/1),({21,16,22,16},2/1),({21,16,21,17},3/1),({21,16,20,18},4/1),({21,16,19,19},1/1),({20,17,23,15},1/1),({20,17,22,16},3/1),({20,17,21,17},5/1),({20,17,20,18},5/1),({20,17,19,19},2/1),({19,18,23,15},1/1),({19,18,22,16},2/1),({19,18,21,17},3/1),({19,18,20,18},3/1),({19,18,19,19},1/1)}, (10,2) => {}, (13,0) => {}, (12,1) => {({23,17,22,19},1/1),({23,17,21,20},1/1),({22,18,23,18},1/1),({22,18,22,19},1/1),({22,18,21,20},1/1),({21,19,23,18},1/1),({21,19,22,19},1/1),({21,19,21,20},1/1),({20,20,23,18},1/1),({20,20,22,19},1/1),({20,20,21,20},1/1)}, (11,2) => {}, (13,1) => {({23,20,23,21},1/1)}, (12,2) => {}, (13,2) => {}, (14,1) => {}, (15,1) => {}, (0,0) => {({1,0,2,0},1/1)}, (1,0) => {({4,0,4,1},1/1),({4,0,3,2},1/1),({3,1,5,0},1/1),({3,1,4,1},1/1),({3,1,3,2},1/1)}, (1,1) => {}, (2,0) => {({7,0,5,3},1/1),({6,1,7,1},1/1),({6,1,6,2},2/1),({6,1,5,3},2/1),({6,1,4,4},1/1),({5,2,7,1},1/1),({5,2,6,2},2/1),({5,2,5,3},2/1),({5,2,4,4},1/1),({4,3,7,1},1/1),({4,3,6,2},2/1),({4,3,5,3},2/1),({4,3,4,4},1/1)}, (2,1) => {({5,5,11,0},1/1)}, (3,0) => {({9,1,8,3},1/1),({9,1,7,4},1/1),({9,1,6,5},1/1),({8,2,9,2},1/1),({8,2,8,3},2/1),({8,2,7,4},3/1),({8,2,6,5},2/1),({7,3,9,2},1/1),({7,3,8,3},3/1),({7,3,7,4},4/1),({7,3,6,5},3/1),({6,4,9,2},2/1),({6,4,8,3},3/1),({6,4,7,4},5/1),({6,4,6,5},3/1),({5,5,8,3},1/1),({5,5,7,4},1/1),({5,5,6,5},1/1)}, (3,1) => {({8,5,13,1},1/1),({8,5,12,2},1/1),({8,5,11,3},1/1)}, (4,0) => {({11,2,10,4},1/1),({11,2,9,5},1/1),({11,2,8,6},1/1),({10,3,10,4},1/1),({10,3,9,5},2/1),({10,3,8,6},2/1),({10,3,7,7},1/1),({9,4,11,3},1/1),({9,4,10,4},2/1),({9,4,9,5},4/1),({9,4,8,6},4/1),({9,4,7,7},2/1),({8,5,11,3},1/1),({8,5,10,4},2/1),({8,5,9,5},4/1),({8,5,8,6},4/1),({8,5,7,7},2/1),({7,6,10,4},2/1),({7,6,9,5},3/1),({7,6,8,6},3/1),({7,6,7,7},1/1)}, (2,2) => {}, (4,1) => {({13,3,11,6},1/1),({12,4,11,6},1/1),({12,4,10,7},1/1),({11,5,14,3},1/1),({11,5,13,4},1/1),({11,5,12,5},1/1),({11,5,11,6},1/1),({11,5,10,7},1/1),({11,5,9,8},1/1),({10,6,15,2},1/1),({10,6,14,3},1/1),({10,6,13,4},3/1),({10,6,12,5},2/1),({10,6,11,6},2/1),({10,6,10,7},1/1),({10,6,9,8},1/1),({9,7,14,3},1/1),({9,7,13,4},1/1),({9,7,12,5},2/1),({9,7,11,6},1/1),({9,7,10,7},1/1),({8,8,15,2},1/1),({8,8,14,3},1/1),({8,8,13,4},2/1),({8,8,12,5},1/1),({8,8,11,6},2/1)}, (5,0) => {({13,3,11,6},1/1),({12,4,10,7},1/1),({11,5,11,6},1/1),({11,5,10,7},1/1),({11,5,9,8},1/1),({10,6,10,7},1/1),({10,6,9,8},1/1),({9,7,11,6},1/1),({9,7,10,7},1/1),({9,7,9,8},1/1)}, (3,2) => {}, (5,1) => {({15,4,13,7},1/1),({15,4,12,8},1/1),({15,4,11,9},1/1),({14,5,14,6},1/1),({14,5,13,7},2/1),({14,5,12,8},3/1),({14,5,11,9},3/1),({14,5,10,10},1/1),({13,6,16,4},1/1),({13,6,15,5},2/1),({13,6,14,6},4/1),({13,6,13,7},6/1),({13,6,12,8},7/1),({13,6,11,9},6/1),({13,6,10,10},2/1),({12,7,16,4},1/1),({12,7,15,5},3/1),({12,7,14,6},6/1),({12,7,13,7},9/1),({12,7,12,8},10/1),({12,7,11,9},8/1),({12,7,10,10},3/1),({11,8,17,3},1/1),({11,8,16,4},2/1),({11,8,15,5},5/1),({11,8,14,6},8/1),({11,8,13,7},11/1),({11,8,12,8},11/1),({11,8,11,9},9/1),({11,8,10,10},3/1),({10,9,16,4},1/1),({10,9,15,5},3/1),({10,9,14,6},5/1),({10,9,13,7},7/1),({10,9,12,8},7/1),({10,9,11,9},5/1),({10,9,10,10},2/1)}, (6,0) => {}, (4,2) => {}};
--dw stands for dominant weights
dw1233 = new HashTable from {(5,2) => {}, (7,0) => {}, (6,1) => {({17,5,14,9},1/1)}, (7,1) => {({19,6,14,12},1/1)}, (8,0) => {}, (6,2) => {}, (7,2) => {}, (9,0) => {}, (8,1) => {({20,8,17,12},1/1)}, (8,2) => {}, (10,0) => {}, (9,1) => {({21,10,19,13},1/1)}, (10,1) => {({22,12,20,15},1/1)}, (11,0) => {}, (9,2) => {}, (12,0) => {}, (11,1) => {({23,14,20,18},1/1)}, (10,2) => {}, (13,0) => {}, (12,1) => {({23,17,22,19},1/1)}, (11,2) => {}, (13,1) => {({23,20,23,21},1/1)}, (12,2) => {}, (14,1) => {}, (13,2) => {}, (15,1) => {}, (0,0) => {({1,0,2,0},1/1)}, (1,0) => {({4,0,4,1},1/1)}, (2,0) => {({7,0,5,3},1/1)}, (1,1) => {}, (3,0) => {({9,1,8,3},1/1)}, (2,1) => {({5,5,11,0},1/1)}, (2,2) => {}, (4,0) => {({11,2,10,4},1/1)}, (3,1) => {({8,5,13,1},1/1)}, (3,2) => {}, (5,0) => {({13,3,11,6},1/1)}, (4,1) => {({13,3,11,6},1/1)}, (4,2) => {}, (6,0) => {}, (5,1) => {({15,4,13,7},1/1)}};
end;