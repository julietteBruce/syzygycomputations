A := QQ[t_0,t_1,t_2,t_3];
--tb stands for Total Betti numbers
tb0033 = new HashTable from {(6,1) => 15444/1, (7,0) => 0, (5,2) => 0, (7,1) => 14157/1, (8,0) => 0, (6,2) => 0, (8,1) => 9152/1, (9,0) => 0, (7,2) => 0, (9,1) => 3861/1, (10,0) => 0, (8,2) => 0, (10,1) => 780/1, (11,0) => 0, (9,2) => 0, (12,0) => 0, (11,1) => 22/1, (10,2) => 165/1, (13,0) => 0, (12,1) => 0, (11,2) => 144/1, (13,1) => 0, (12,2) => 39/1, (13,2) => 4/1, (14,2) => 0, (0,0) => 1/1, (1,0) => 0, (1,1) => 87/1, (2,0) => 0, (2,1) => 676/1, (3,0) => 0, (3,1) => 2691/1, (4,0) => 0, (2,2) => 0, (4,1) => 6864/1, (5,0) => 0, (3,2) => 0, (5,1) => 12155/1, (6,0) => 0, (4,2) => 0};
--sb represents the betti numbers as sums of Schur functors
sb0033 = new HashTable from {(6,1) => {({17,4,12,9},1/1),({17,4,11,10},1/1),({16,5,14,7},2/1),({16,5,13,8},3/1),({16,5,12,9},5/1),({16,5,11,10},3/1),({15,6,15,6},2/1),({15,6,14,7},6/1),({15,6,13,8},11/1),({15,6,12,9},13/1),({15,6,11,10},9/1),({14,7,16,5},2/1),({14,7,15,6},6/1),({14,7,14,7},15/1),({14,7,13,8},23/1),({14,7,12,9},26/1),({14,7,11,10},17/1),({13,8,16,5},3/1),({13,8,15,6},11/1),({13,8,14,7},23/1),({13,8,13,8},34/1),({13,8,12,9},36/1),({13,8,11,10},24/1),({12,9,17,4},1/1),({12,9,16,5},5/1),({12,9,15,6},13/1),({12,9,14,7},26/1),({12,9,13,8},36/1),({12,9,12,9},38/1),({12,9,11,10},24/1),({11,10,17,4},1/1),({11,10,16,5},3/1),({11,10,15,6},9/1),({11,10,14,7},17/1),({11,10,13,8},24/1),({11,10,12,9},24/1),({11,10,11,10},15/1)}, (7,0) => {}, (5,2) => {}, (7,1) => {({19,5,12,12},1/1),({18,6,15,9},1/1),({18,6,14,10},2/1),({18,6,13,11},2/1),({18,6,12,12},1/1),({17,7,16,8},2/1),({17,7,15,9},4/1),({17,7,14,10},8/1),({17,7,13,11},6/1),({17,7,12,12},4/1),({16,8,17,7},2/1),({16,8,16,8},6/1),({16,8,15,9},13/1),({16,8,14,10},17/1),({16,8,13,11},16/1),({16,8,12,12},6/1),({15,9,18,6},1/1),({15,9,17,7},4/1),({15,9,16,8},13/1),({15,9,15,9},21/1),({15,9,14,10},30/1),({15,9,13,11},24/1),({15,9,12,12},12/1),({14,10,18,6},2/1),({14,10,17,7},8/1),({14,10,16,8},17/1),({14,10,15,9},30/1),({14,10,14,10},35/1),({14,10,13,11},31/1),({14,10,12,12},11/1),({13,11,18,6},2/1),({13,11,17,7},6/1),({13,11,16,8},16/1),({13,11,15,9},24/1),({13,11,14,10},31/1),({13,11,13,11},23/1),({13,11,12,12},12/1),({12,12,19,5},1/1),({12,12,18,6},1/1),({12,12,17,7},4/1),({12,12,16,8},6/1),({12,12,15,9},12/1),({12,12,14,10},11/1),({12,12,13,11},12/1),({12,12,12,12},2/1)}, (8,0) => {}, (6,2) => {}, (8,1) => {({20,7,15,12},1/1),({19,8,17,10},1/1),({19,8,16,11},2/1),({19,8,15,12},3/1),({19,8,14,13},2/1),({18,9,17,10},4/1),({18,9,16,11},6/1),({18,9,15,12},8/1),({18,9,14,13},6/1),({17,10,19,8},1/1),({17,10,18,9},4/1),({17,10,17,10},8/1),({17,10,16,11},14/1),({17,10,15,12},16/1),({17,10,14,13},10/1),({16,11,19,8},2/1),({16,11,18,9},6/1),({16,11,17,10},14/1),({16,11,16,11},20/1),({16,11,15,12},22/1),({16,11,14,13},14/1),({15,12,20,7},1/1),({15,12,19,8},3/1),({15,12,18,9},8/1),({15,12,17,10},16/1),({15,12,16,11},22/1),({15,12,15,12},22/1),({15,12,14,13},15/1),({14,13,19,8},2/1),({14,13,18,9},6/1),({14,13,17,10},10/1),({14,13,16,11},14/1),({14,13,15,12},15/1),({14,13,14,13},8/1)}, (9,0) => {}, (7,2) => {}, (9,1) => {({21,9,17,13},1/1),({21,9,15,15},1/1),({20,10,18,12},1/1),({20,10,17,13},2/1),({20,10,16,14},2/1),({20,10,15,15},1/1),({19,11,19,11},1/1),({19,11,18,12},3/1),({19,11,17,13},6/1),({19,11,16,14},4/1),({19,11,15,15},3/1),({18,12,20,10},1/1),({18,12,19,11},3/1),({18,12,18,12},6/1),({18,12,17,13},8/1),({18,12,16,14},8/1),({18,12,15,15},3/1),({17,13,21,9},1/1),({17,13,20,10},2/1),({17,13,19,11},6/1),({17,13,18,12},8/1),({17,13,17,13},12/1),({17,13,16,14},8/1),({17,13,15,15},5/1),({16,14,20,10},2/1),({16,14,19,11},4/1),({16,14,18,12},8/1),({16,14,17,13},8/1),({16,14,16,14},8/1),({16,14,15,15},2/1),({15,15,21,9},1/1),({15,15,20,10},1/1),({15,15,19,11},3/1),({15,15,18,12},3/1),({15,15,17,13},5/1),({15,15,16,14},2/1),({15,15,15,15},2/1)}, (10,0) => {}, (8,2) => {}, (10,1) => {({22,11,18,15},1/1),({21,12,18,15},1/1),({21,12,17,16},1/1),({20,13,20,13},1/1),({20,13,19,14},1/1),({20,13,18,15},2/1),({20,13,17,16},1/1),({19,14,20,13},1/1),({19,14,19,14},2/1),({19,14,18,15},2/1),({19,14,17,16},1/1),({18,15,22,11},1/1),({18,15,21,12},1/1),({18,15,20,13},2/1),({18,15,19,14},2/1),({18,15,18,15},2/1),({18,15,17,16},1/1),({17,16,21,12},1/1),({17,16,20,13},1/1),({17,16,19,14},1/1),({17,16,18,15},1/1),({17,16,17,16},1/1)}, (11,0) => {}, (9,2) => {}, (12,0) => {}, (11,1) => {({23,13,18,18},1/1),({18,18,23,13},1/1)}, (10,2) => {({21,15,21,15},1/1),({21,15,19,17},1/1),({20,16,20,16},1/1),({20,16,19,17},1/1),({19,17,21,15},1/1),({19,17,20,16},1/1),({19,17,19,17},2/1),({18,18,18,18},1/1)}, (13,0) => {}, (12,1) => {}, (11,2) => {({22,17,22,17},1/1),({22,17,21,18},1/1),({22,17,20,19},1/1),({21,18,22,17},1/1),({21,18,21,18},1/1),({21,18,20,19},1/1),({20,19,22,17},1/1),({20,19,21,18},1/1),({20,19,20,19},1/1)}, (13,1) => {}, (12,2) => {({23,19,22,20},1/1),({22,20,23,19},1/1),({22,20,22,20},1/1)}, (13,2) => {({23,22,23,22},1/1)}, (14,2) => {}, (0,0) => {({0,0,0,0},1/1)}, (1,0) => {}, (1,1) => {({6,0,4,2},1/1),({5,1,5,1},1/1),({5,1,3,3},1/1),({4,2,6,0},1/1),({4,2,4,2},1/1),({3,3,5,1},1/1),({3,3,3,3},1/1)}, (2,0) => {}, (2,1) => {({9,0,5,4},1/1),({8,1,7,2},1/1),({8,1,6,3},2/1),({8,1,5,4},1/1),({7,2,8,1},1/1),({7,2,7,2},2/1),({7,2,6,3},3/1),({7,2,5,4},2/1),({6,3,8,1},2/1),({6,3,7,2},3/1),({6,3,6,3},4/1),({6,3,5,4},3/1),({5,4,9,0},1/1),({5,4,8,1},1/1),({5,4,7,2},2/1),({5,4,6,3},3/1),({5,4,5,4},1/1)}, (3,0) => {}, (3,1) => {({11,1,8,4},1/1),({11,1,7,5},1/1),({11,1,6,6},1/1),({10,2,9,3},2/1),({10,2,8,4},3/1),({10,2,7,5},4/1),({10,2,6,6},1/1),({9,3,10,2},2/1),({9,3,9,3},4/1),({9,3,8,4},8/1),({9,3,7,5},6/1),({9,3,6,6},4/1),({8,4,11,1},1/1),({8,4,10,2},3/1),({8,4,9,3},8/1),({8,4,8,4},10/1),({8,4,7,5},10/1),({8,4,6,6},3/1),({7,5,11,1},1/1),({7,5,10,2},4/1),({7,5,9,3},6/1),({7,5,8,4},10/1),({7,5,7,5},7/1),({7,5,6,6},5/1),({6,6,11,1},1/1),({6,6,10,2},1/1),({6,6,9,3},4/1),({6,6,8,4},3/1),({6,6,7,5},5/1)}, (4,0) => {}, (2,2) => {}, (4,1) => {({13,2,10,5},1/1),({13,2,9,6},2/1),({13,2,8,7},1/1),({12,3,11,4},2/1),({12,3,10,5},4/1),({12,3,9,6},6/1),({12,3,8,7},4/1),({11,4,12,3},2/1),({11,4,11,4},6/1),({11,4,10,5},11/1),({11,4,9,6},13/1),({11,4,8,7},9/1),({10,5,13,2},1/1),({10,5,12,3},4/1),({10,5,11,4},11/1),({10,5,10,5},17/1),({10,5,9,6},20/1),({10,5,8,7},13/1),({9,6,13,2},2/1),({9,6,12,3},6/1),({9,6,11,4},13/1),({9,6,10,5},20/1),({9,6,9,6},21/1),({9,6,8,7},14/1),({8,7,13,2},1/1),({8,7,12,3},4/1),({8,7,11,4},9/1),({8,7,10,5},13/1),({8,7,9,6},14/1),({8,7,8,7},9/1)}, (5,0) => {}, (3,2) => {}, (5,1) => {({15,3,11,7},2/1),({15,3,10,8},1/1),({15,3,9,9},1/1),({14,4,13,5},1/1),({14,4,12,6},4/1),({14,4,11,7},5/1),({14,4,10,8},6/1),({14,4,9,9},3/1),({13,5,14,4},1/1),({13,5,13,5},5/1),({13,5,12,6},10/1),({13,5,11,7},16/1),({13,5,10,8},13/1),({13,5,9,9},7/1),({12,6,14,4},4/1),({12,6,13,5},10/1),({12,6,12,6},20/1),({12,6,11,7},26/1),({12,6,10,8},24/1),({12,6,9,9},8/1),({11,7,15,3},2/1),({11,7,14,4},5/1),({11,7,13,5},16/1),({11,7,12,6},26/1),({11,7,11,7},35/1),({11,7,10,8},27/1),({11,7,9,9},14/1),({10,8,15,3},1/1),({10,8,14,4},6/1),({10,8,13,5},13/1),({10,8,12,6},24/1),({10,8,11,7},27/1),({10,8,10,8},26/1),({10,8,9,9},8/1),({9,9,15,3},1/1),({9,9,14,4},3/1),({9,9,13,5},7/1),({9,9,12,6},8/1),({9,9,11,7},14/1),({9,9,10,8},8/1),({9,9,9,9},5/1)}, (6,0) => {}, (4,2) => {}};
--dw stands for dominant weights
dw0033 = new HashTable from {(5,2) => {}, (7,0) => {}, (6,1) => {({17,4,12,9},1/1)}, (7,1) => {({19,5,12,12},1/1)}, (8,0) => {}, (6,2) => {}, (7,2) => {}, (9,0) => {}, (8,1) => {({20,7,15,12},1/1)}, (8,2) => {}, (10,0) => {}, (9,1) => {({21,9,17,13},1/1)}, (10,1) => {({22,11,18,15},1/1)}, (11,0) => {}, (9,2) => {}, (12,0) => {}, (11,1) => {({23,13,18,18},1/1)}, (10,2) => {({21,15,21,15},1/1)}, (13,0) => {}, (12,1) => {}, (11,2) => {({22,17,22,17},1/1)}, (13,1) => {}, (12,2) => {({23,19,22,20},1/1)}, (13,2) => {({23,22,23,22},1/1)}, (14,2) => {}, (0,0) => {({0,0,0,0},1/1)}, (1,0) => {}, (2,0) => {}, (1,1) => {({6,0,4,2},1/1)}, (3,0) => {}, (2,1) => {({9,0,5,4},1/1)}, (2,2) => {}, (4,0) => {}, (3,1) => {({11,1,8,4},1/1)}, (3,2) => {}, (5,0) => {}, (4,1) => {({13,2,10,5},1/1)}, (4,2) => {}, (6,0) => {}, (5,1) => {({15,3,11,7},2/1)}};
end;