A := QQ[t_0,t_1,t_2,t_3];
--tb stands for Total Betti numbers
tb1133 = new HashTable from {(5,2) => 0, (7,0) => 0, (6,1) => 14157/1, (7,1) => 15444/1, (8,0) => 0, (6,2) => 0, (7,2) => 0, (9,0) => 0, (8,1) => 12155/1, (8,2) => 0, (10,0) => 0, (9,1) => 6864/1, (10,1) => 2691/1, (11,0) => 0, (9,2) => 0, (12,0) => 0, (11,1) => 676/1, (10,2) => 0, (13,0) => 0, (12,1) => 87/1, (11,2) => 0, (13,1) => 0, (12,2) => 0, (13,2) => 1/1, (14,2) => 0, (0,0) => 4/1, (1,0) => 39/1, (2,0) => 144/1, (1,1) => 0, (3,0) => 165/1, (2,1) => 22/1, (2,2) => 0, (4,0) => 0, (3,1) => 780/1, (3,2) => 0, (5,0) => 0, (4,1) => 3861/1, (4,2) => 0, (6,0) => 0, (5,1) => 9152/1};
--sb represents the betti numbers as sums of Schur functors
sb1133 = new HashTable from {(5,2) => {}, (7,0) => {}, (6,1) => {({18,4,11,11},1/1),({17,5,14,8},1/1),({17,5,13,9},2/1),({17,5,12,10},2/1),({17,5,11,11},1/1),({16,6,15,7},2/1),({16,6,14,8},4/1),({16,6,13,9},8/1),({16,6,12,10},6/1),({16,6,11,11},4/1),({15,7,16,6},2/1),({15,7,15,7},6/1),({15,7,14,8},13/1),({15,7,13,9},17/1),({15,7,12,10},16/1),({15,7,11,11},6/1),({14,8,17,5},1/1),({14,8,16,6},4/1),({14,8,15,7},13/1),({14,8,14,8},21/1),({14,8,13,9},30/1),({14,8,12,10},24/1),({14,8,11,11},12/1),({13,9,17,5},2/1),({13,9,16,6},8/1),({13,9,15,7},17/1),({13,9,14,8},30/1),({13,9,13,9},35/1),({13,9,12,10},31/1),({13,9,11,11},11/1),({12,10,17,5},2/1),({12,10,16,6},6/1),({12,10,15,7},16/1),({12,10,14,8},24/1),({12,10,13,9},31/1),({12,10,12,10},23/1),({12,10,11,11},12/1),({11,11,18,4},1/1),({11,11,17,5},1/1),({11,11,16,6},4/1),({11,11,15,7},6/1),({11,11,14,8},12/1),({11,11,13,9},11/1),({11,11,12,10},12/1),({11,11,11,11},2/1)}, (7,1) => {({19,6,14,11},1/1),({19,6,13,12},1/1),({18,7,16,9},2/1),({18,7,15,10},3/1),({18,7,14,11},5/1),({18,7,13,12},3/1),({17,8,17,8},2/1),({17,8,16,9},6/1),({17,8,15,10},11/1),({17,8,14,11},13/1),({17,8,13,12},9/1),({16,9,18,7},2/1),({16,9,17,8},6/1),({16,9,16,9},15/1),({16,9,15,10},23/1),({16,9,14,11},26/1),({16,9,13,12},17/1),({15,10,18,7},3/1),({15,10,17,8},11/1),({15,10,16,9},23/1),({15,10,15,10},34/1),({15,10,14,11},36/1),({15,10,13,12},24/1),({14,11,19,6},1/1),({14,11,18,7},5/1),({14,11,17,8},13/1),({14,11,16,9},26/1),({14,11,15,10},36/1),({14,11,14,11},38/1),({14,11,13,12},24/1),({13,12,19,6},1/1),({13,12,18,7},3/1),({13,12,17,8},9/1),({13,12,16,9},17/1),({13,12,15,10},24/1),({13,12,14,11},24/1),({13,12,13,12},15/1)}, (8,0) => {}, (6,2) => {}, (7,2) => {}, (9,0) => {}, (8,1) => {({20,8,16,12},2/1),({20,8,15,13},1/1),({20,8,14,14},1/1),({19,9,18,10},1/1),({19,9,17,11},4/1),({19,9,16,12},5/1),({19,9,15,13},6/1),({19,9,14,14},3/1),({18,10,19,9},1/1),({18,10,18,10},5/1),({18,10,17,11},10/1),({18,10,16,12},16/1),({18,10,15,13},13/1),({18,10,14,14},7/1),({17,11,19,9},4/1),({17,11,18,10},10/1),({17,11,17,11},20/1),({17,11,16,12},26/1),({17,11,15,13},24/1),({17,11,14,14},8/1),({16,12,20,8},2/1),({16,12,19,9},5/1),({16,12,18,10},16/1),({16,12,17,11},26/1),({16,12,16,12},35/1),({16,12,15,13},27/1),({16,12,14,14},14/1),({15,13,20,8},1/1),({15,13,19,9},6/1),({15,13,18,10},13/1),({15,13,17,11},24/1),({15,13,16,12},27/1),({15,13,15,13},26/1),({15,13,14,14},8/1),({14,14,20,8},1/1),({14,14,19,9},3/1),({14,14,18,10},7/1),({14,14,17,11},8/1),({14,14,16,12},14/1),({14,14,15,13},8/1),({14,14,14,14},5/1)}, (8,2) => {}, (10,0) => {}, (9,1) => {({21,10,18,13},1/1),({21,10,17,14},2/1),({21,10,16,15},1/1),({20,11,19,12},2/1),({20,11,18,13},4/1),({20,11,17,14},6/1),({20,11,16,15},4/1),({19,12,20,11},2/1),({19,12,19,12},6/1),({19,12,18,13},11/1),({19,12,17,14},13/1),({19,12,16,15},9/1),({18,13,21,10},1/1),({18,13,20,11},4/1),({18,13,19,12},11/1),({18,13,18,13},17/1),({18,13,17,14},20/1),({18,13,16,15},13/1),({17,14,21,10},2/1),({17,14,20,11},6/1),({17,14,19,12},13/1),({17,14,18,13},20/1),({17,14,17,14},21/1),({17,14,16,15},14/1),({16,15,21,10},1/1),({16,15,20,11},4/1),({16,15,19,12},9/1),({16,15,18,13},13/1),({16,15,17,14},14/1),({16,15,16,15},9/1)}, (10,1) => {({22,12,19,15},1/1),({22,12,18,16},1/1),({22,12,17,17},1/1),({21,13,20,14},2/1),({21,13,19,15},3/1),({21,13,18,16},4/1),({21,13,17,17},1/1),({20,14,21,13},2/1),({20,14,20,14},4/1),({20,14,19,15},8/1),({20,14,18,16},6/1),({20,14,17,17},4/1),({19,15,22,12},1/1),({19,15,21,13},3/1),({19,15,20,14},8/1),({19,15,19,15},10/1),({19,15,18,16},10/1),({19,15,17,17},3/1),({18,16,22,12},1/1),({18,16,21,13},4/1),({18,16,20,14},6/1),({18,16,19,15},10/1),({18,16,18,16},7/1),({18,16,17,17},5/1),({17,17,22,12},1/1),({17,17,21,13},1/1),({17,17,20,14},4/1),({17,17,19,15},3/1),({17,17,18,16},5/1)}, (11,0) => {}, (9,2) => {}, (12,0) => {}, (11,1) => {({23,14,19,18},1/1),({22,15,21,16},1/1),({22,15,20,17},2/1),({22,15,19,18},1/1),({21,16,22,15},1/1),({21,16,21,16},2/1),({21,16,20,17},3/1),({21,16,19,18},2/1),({20,17,22,15},2/1),({20,17,21,16},3/1),({20,17,20,17},4/1),({20,17,19,18},3/1),({19,18,23,14},1/1),({19,18,22,15},1/1),({19,18,21,16},2/1),({19,18,20,17},3/1),({19,18,19,18},1/1)}, (10,2) => {}, (13,0) => {}, (12,1) => {({23,17,21,19},1/1),({22,18,22,18},1/1),({22,18,20,20},1/1),({21,19,23,17},1/1),({21,19,21,19},1/1),({20,20,22,18},1/1),({20,20,20,20},1/1)}, (11,2) => {}, (13,1) => {}, (12,2) => {}, (13,2) => {({23,23,23,23},1/1)}, (14,2) => {}, (0,0) => {({1,0,1,0},1/1)}, (1,0) => {({4,0,3,1},1/1),({3,1,4,0},1/1),({3,1,3,1},1/1)}, (2,0) => {({6,1,6,1},1/1),({6,1,5,2},1/1),({6,1,4,3},1/1),({5,2,6,1},1/1),({5,2,5,2},1/1),({5,2,4,3},1/1),({4,3,6,1},1/1),({4,3,5,2},1/1),({4,3,4,3},1/1)}, (1,1) => {}, (3,0) => {({8,2,8,2},1/1),({8,2,6,4},1/1),({7,3,7,3},1/1),({7,3,6,4},1/1),({6,4,8,2},1/1),({6,4,7,3},1/1),({6,4,6,4},2/1),({5,5,5,5},1/1)}, (2,1) => {({10,0,5,5},1/1),({5,5,10,0},1/1)}, (2,2) => {}, (4,0) => {}, (3,1) => {({12,1,8,5},1/1),({11,2,8,5},1/1),({11,2,7,6},1/1),({10,3,10,3},1/1),({10,3,9,4},1/1),({10,3,8,5},2/1),({10,3,7,6},1/1),({9,4,10,3},1/1),({9,4,9,4},2/1),({9,4,8,5},2/1),({9,4,7,6},1/1),({8,5,12,1},1/1),({8,5,11,2},1/1),({8,5,10,3},2/1),({8,5,9,4},2/1),({8,5,8,5},2/1),({8,5,7,6},1/1),({7,6,11,2},1/1),({7,6,10,3},1/1),({7,6,9,4},1/1),({7,6,8,5},1/1),({7,6,7,6},1/1)}, (3,2) => {}, (5,0) => {}, (4,1) => {({14,2,10,6},1/1),({14,2,8,8},1/1),({13,3,11,5},1/1),({13,3,10,6},2/1),({13,3,9,7},2/1),({13,3,8,8},1/1),({12,4,12,4},1/1),({12,4,11,5},3/1),({12,4,10,6},6/1),({12,4,9,7},4/1),({12,4,8,8},3/1),({11,5,13,3},1/1),({11,5,12,4},3/1),({11,5,11,5},6/1),({11,5,10,6},8/1),({11,5,9,7},8/1),({11,5,8,8},3/1),({10,6,14,2},1/1),({10,6,13,3},2/1),({10,6,12,4},6/1),({10,6,11,5},8/1),({10,6,10,6},12/1),({10,6,9,7},8/1),({10,6,8,8},5/1),({9,7,13,3},2/1),({9,7,12,4},4/1),({9,7,11,5},8/1),({9,7,10,6},8/1),({9,7,9,7},8/1),({9,7,8,8},2/1),({8,8,14,2},1/1),({8,8,13,3},1/1),({8,8,12,4},3/1),({8,8,11,5},3/1),({8,8,10,6},5/1),({8,8,9,7},2/1),({8,8,8,8},2/1)}, (4,2) => {}, (6,0) => {}, (5,1) => {({16,3,11,8},1/1),({15,4,13,6},1/1),({15,4,12,7},2/1),({15,4,11,8},3/1),({15,4,10,9},2/1),({14,5,13,6},4/1),({14,5,12,7},6/1),({14,5,11,8},8/1),({14,5,10,9},6/1),({13,6,15,4},1/1),({13,6,14,5},4/1),({13,6,13,6},8/1),({13,6,12,7},14/1),({13,6,11,8},16/1),({13,6,10,9},10/1),({12,7,15,4},2/1),({12,7,14,5},6/1),({12,7,13,6},14/1),({12,7,12,7},20/1),({12,7,11,8},22/1),({12,7,10,9},14/1),({11,8,16,3},1/1),({11,8,15,4},3/1),({11,8,14,5},8/1),({11,8,13,6},16/1),({11,8,12,7},22/1),({11,8,11,8},22/1),({11,8,10,9},15/1),({10,9,15,4},2/1),({10,9,14,5},6/1),({10,9,13,6},10/1),({10,9,12,7},14/1),({10,9,11,8},15/1),({10,9,10,9},8/1)}};
--dw stands for dominant weights
dw1133 = new HashTable from {(6,1) => {({18,4,11,11},1/1),({17,5,14,8},1/1),({16,6,15,7},2/1),({15,7,16,6},2/1),({14,8,17,5},1/1),({11,11,18,4},1/1)}, (7,0) => {}, (5,2) => {}, (7,1) => {({19,6,14,11},1/1),({18,7,16,9},2/1),({17,8,17,8},2/1),({16,9,18,7},2/1),({14,11,19,6},1/1)}, (8,0) => {}, (6,2) => {}, (8,1) => {({20,8,16,12},2/1),({19,9,18,10},1/1),({18,10,19,9},1/1),({16,12,20,8},2/1)}, (9,0) => {}, (7,2) => {}, (9,1) => {({21,10,18,13},1/1),({20,11,19,12},2/1),({19,12,20,11},2/1),({18,13,21,10},1/1)}, (10,0) => {}, (8,2) => {}, (10,1) => {({22,12,19,15},1/1),({21,13,20,14},2/1),({20,14,21,13},2/1),({19,15,22,12},1/1)}, (11,0) => {}, (9,2) => {}, (12,0) => {}, (11,1) => {({23,14,19,18},1/1),({22,15,21,16},1/1),({21,16,22,15},1/1),({19,18,23,14},1/1)}, (10,2) => {}, (13,0) => {}, (12,1) => {({23,17,21,19},1/1),({22,18,22,18},1/1),({21,19,23,17},1/1)}, (11,2) => {}, (13,1) => {}, (12,2) => {}, (13,2) => {({23,23,23,23},1/1)}, (14,2) => {}, (0,0) => {({1,0,1,0},1/1)}, (1,0) => {({4,0,3,1},1/1),({3,1,4,0},1/1)}, (1,1) => {}, (2,0) => {({6,1,6,1},1/1)}, (2,1) => {({10,0,5,5},1/1),({5,5,10,0},1/1)}, (3,0) => {({8,2,8,2},1/1)}, (3,1) => {({12,1,8,5},1/1),({10,3,10,3},1/1),({8,5,12,1},1/1)}, (4,0) => {}, (2,2) => {}, (4,1) => {({14,2,10,6},1/1),({13,3,11,5},1/1),({12,4,12,4},1/1),({11,5,13,3},1/1),({10,6,14,2},1/1)}, (5,0) => {}, (3,2) => {}, (5,1) => {({16,3,11,8},1/1),({15,4,13,6},1/1),({13,6,15,4},1/1),({11,8,16,3},1/1)}, (6,0) => {}, (4,2) => {}};
end;