A := QQ[t_0,t_1,t_2,t_3];
--tb stands for Total Betti numbers
tb0133 = new HashTable from {(6,1) => 18018/1, (7,0) => 0, (5,2) => 0, (7,1) => 18018/1, (8,0) => 0, (6,2) => 0, (8,1) => 13156/1, (9,0) => 0, (7,2) => 0, (9,1) => 6864/1, (10,0) => 0, (8,2) => 0, (10,1) => 2418/1, (11,0) => 0, (9,2) => 0, (12,0) => 0, (11,1) => 494/1, (10,2) => 0, (13,0) => 0, (12,1) => 24/1, (11,2) => 0, (13,1) => 0, (12,2) => 12/1, (13,2) => 2/1, (14,2) => 0, (0,0) => 2/1, (1,0) => 12/1, (1,1) => 24/1, (2,0) => 0, (2,1) => 494/1, (3,0) => 0, (3,1) => 2418/1, (4,0) => 0, (2,2) => 0, (4,1) => 6864/1, (5,0) => 0, (3,2) => 0, (5,1) => 13156/1, (6,0) => 0, (4,2) => 0};
--sb represents the betti numbers as sums of Schur functors
sb0133 = new HashTable from {(6,1) => {({17,4,13,9},1/1),({17,4,12,10},1/1),({17,4,11,11},1/1),({16,5,15,7},1/1),({16,5,14,8},3/1),({16,5,13,9},5/1),({16,5,12,10},5/1),({16,5,11,11},2/1),({15,6,16,6},1/1),({15,6,15,7},5/1),({15,6,14,8},10/1),({15,6,13,9},15/1),({15,6,12,10},13/1),({15,6,11,11},6/1),({14,7,17,5},1/1),({14,7,16,6},4/1),({14,7,15,7},12/1),({14,7,14,8},22/1),({14,7,13,9},29/1),({14,7,12,10},26/1),({14,7,11,11},10/1),({13,8,17,5},2/1),({13,8,16,6},8/1),({13,8,15,7},20/1),({13,8,14,8},33/1),({13,8,13,9},42/1),({13,8,12,10},35/1),({13,8,11,11},15/1),({12,9,17,5},3/1),({12,9,16,6},10/1),({12,9,15,7},22/1),({12,9,14,8},36/1),({12,9,13,9},43/1),({12,9,12,10},36/1),({12,9,11,11},14/1),({11,10,18,4},1/1),({11,10,17,5},2/1),({11,10,16,6},7/1),({11,10,15,7},15/1),({11,10,14,8},24/1),({11,10,13,9},28/1),({11,10,12,10},23/1),({11,10,11,11},9/1)}, (7,0) => {}, (5,2) => {}, (7,1) => {({19,5,13,12},1/1),({18,6,16,9},1/1),({18,6,15,10},2/1),({18,6,14,11},3/1),({18,6,13,12},2/1),({17,7,17,8},1/1),({17,7,16,9},4/1),({17,7,15,10},8/1),({17,7,14,11},10/1),({17,7,13,12},7/1),({16,8,18,7},1/1),({16,8,17,8},5/1),({16,8,16,9},12/1),({16,8,15,10},20/1),({16,8,14,11},22/1),({16,8,13,12},15/1),({15,9,18,7},3/1),({15,9,17,8},10/1),({15,9,16,9},22/1),({15,9,15,10},33/1),({15,9,14,11},36/1),({15,9,13,12},24/1),({14,10,19,6},1/1),({14,10,18,7},5/1),({14,10,17,8},15/1),({14,10,16,9},29/1),({14,10,15,10},42/1),({14,10,14,11},43/1),({14,10,13,12},28/1),({13,11,19,6},1/1),({13,11,18,7},5/1),({13,11,17,8},13/1),({13,11,16,9},26/1),({13,11,15,10},35/1),({13,11,14,11},36/1),({13,11,13,12},23/1),({12,12,19,6},1/1),({12,12,18,7},2/1),({12,12,17,8},6/1),({12,12,16,9},10/1),({12,12,15,10},15/1),({12,12,14,11},14/1),({12,12,13,12},9/1)}, (8,0) => {}, (6,2) => {}, (8,1) => {({20,7,16,12},1/1),({20,7,15,13},1/1),({19,8,18,10},1/1),({19,8,17,11},2/1),({19,8,16,12},4/1),({19,8,15,13},4/1),({19,8,14,14},2/1),({18,9,18,10},3/1),({18,9,17,11},8/1),({18,9,16,12},11/1),({18,9,15,13},11/1),({18,9,14,14},5/1),({17,10,19,9},3/1),({17,10,18,10},8/1),({17,10,17,11},16/1),({17,10,16,12},23/1),({17,10,15,13},20/1),({17,10,14,14},8/1),({16,11,20,8},1/1),({16,11,19,9},5/1),({16,11,18,10},14/1),({16,11,17,11},25/1),({16,11,16,12},32/1),({16,11,15,13},28/1),({16,11,14,14},11/1),({15,12,20,8},2/1),({15,12,19,9},6/1),({15,12,18,10},16/1),({15,12,17,11},27/1),({15,12,16,12},33/1),({15,12,15,13},28/1),({15,12,14,14},12/1),({14,13,20,8},1/1),({14,13,19,9},5/1),({14,13,18,10},11/1),({14,13,17,11},17/1),({14,13,16,12},22/1),({14,13,15,13},18/1),({14,13,14,14},6/1)}, (9,0) => {}, (7,2) => {}, (9,1) => {({21,9,18,13},1/1),({21,9,17,14},1/1),({21,9,16,15},1/1),({20,10,19,12},1/1),({20,10,18,13},3/1),({20,10,17,14},4/1),({20,10,16,15},3/1),({19,11,20,11},1/1),({19,11,19,12},4/1),({19,11,18,13},8/1),({19,11,17,14},10/1),({19,11,16,15},7/1),({18,12,20,11},3/1),({18,12,19,12},8/1),({18,12,18,13},14/1),({18,12,17,14},16/1),({18,12,16,15},11/1),({17,13,21,10},2/1),({17,13,20,11},5/1),({17,13,19,12},12/1),({17,13,18,13},18/1),({17,13,17,14},20/1),({17,13,16,15},13/1),({16,14,21,10},1/1),({16,14,20,11},5/1),({16,14,19,12},10/1),({16,14,18,13},16/1),({16,14,17,14},16/1),({16,14,16,15},11/1),({15,15,21,10},1/1),({15,15,20,11},2/1),({15,15,19,12},5/1),({15,15,18,13},6/1),({15,15,17,14},7/1),({15,15,16,15},4/1)}, (10,0) => {}, (8,2) => {}, (10,1) => {({22,11,19,15},1/1),({22,11,18,16},1/1),({21,12,20,14},1/1),({21,12,19,15},2/1),({21,12,18,16},3/1),({21,12,17,17},1/1),({20,13,21,13},1/1),({20,13,20,14},3/1),({20,13,19,15},5/1),({20,13,18,16},5/1),({20,13,17,17},2/1),({19,14,21,13},2/1),({19,14,20,14},5/1),({19,14,19,15},8/1),({19,14,18,16},7/1),({19,14,17,17},3/1),({18,15,22,12},1/1),({18,15,21,13},3/1),({18,15,20,14},6/1),({18,15,19,15},8/1),({18,15,18,16},7/1),({18,15,17,17},3/1),({17,16,22,12},1/1),({17,16,21,13},2/1),({17,16,20,14},4/1),({17,16,19,15},5/1),({17,16,18,16},5/1),({17,16,17,17},2/1)}, (11,0) => {}, (9,2) => {}, (12,0) => {}, (11,1) => {({23,13,19,18},1/1),({22,14,21,16},1/1),({22,14,20,17},1/1),({22,14,19,18},1/1),({21,15,21,16},1/1),({21,15,20,17},2/1),({21,15,19,18},1/1),({20,16,22,15},1/1),({20,16,21,16},2/1),({20,16,20,17},2/1),({20,16,19,18},2/1),({19,17,22,15},1/1),({19,17,21,16},1/1),({19,17,20,17},2/1),({19,17,19,18},1/1),({18,18,23,14},1/1),({18,18,21,16},1/1),({18,18,20,17},1/1)}, (10,2) => {}, (13,0) => {}, (12,1) => {({23,16,21,19},1/1)}, (11,2) => {}, (13,1) => {}, (12,2) => {({22,20,23,20},1/1)}, (13,2) => {({23,22,23,23},1/1)}, (14,2) => {}, (0,0) => {({0,0,1,0},1/1)}, (1,0) => {({3,0,3,1},1/1)}, (1,1) => {({4,2,7,0},1/1)}, (2,0) => {}, (2,1) => {({9,0,5,5},1/1),({8,1,7,3},1/1),({8,1,6,4},1/1),({7,2,9,1},1/1),({7,2,8,2},1/1),({7,2,7,3},2/1),({7,2,6,4},1/1),({7,2,5,5},1/1),({6,3,9,1},1/1),({6,3,8,2},2/1),({6,3,7,3},2/1),({6,3,6,4},2/1),({6,3,5,5},1/1),({5,4,10,0},1/1),({5,4,9,1},1/1),({5,4,8,2},1/1),({5,4,7,3},2/1),({5,4,6,4},1/1)}, (3,0) => {}, (3,1) => {({11,1,8,5},1/1),({11,1,7,6},1/1),({10,2,10,3},1/1),({10,2,9,4},2/1),({10,2,8,5},3/1),({10,2,7,6},2/1),({9,3,11,2},1/1),({9,3,10,3},3/1),({9,3,9,4},5/1),({9,3,8,5},6/1),({9,3,7,6},4/1),({8,4,12,1},1/1),({8,4,11,2},2/1),({8,4,10,3},5/1),({8,4,9,4},8/1),({8,4,8,5},8/1),({8,4,7,6},5/1),({7,5,12,1},1/1),({7,5,11,2},3/1),({7,5,10,3},5/1),({7,5,9,4},7/1),({7,5,8,5},7/1),({7,5,7,6},5/1),({6,6,11,2},1/1),({6,6,10,3},2/1),({6,6,9,4},3/1),({6,6,8,5},3/1),({6,6,7,6},2/1)}, (4,0) => {}, (2,2) => {}, (4,1) => {({13,2,10,6},2/1),({13,2,9,7},1/1),({13,2,8,8},1/1),({12,3,12,4},1/1),({12,3,11,5},3/1),({12,3,10,6},5/1),({12,3,9,7},5/1),({12,3,8,8},2/1),({11,4,13,3},1/1),({11,4,12,4},4/1),({11,4,11,5},8/1),({11,4,10,6},12/1),({11,4,9,7},10/1),({11,4,8,8},5/1),({10,5,14,2},1/1),({10,5,13,3},3/1),({10,5,12,4},8/1),({10,5,11,5},14/1),({10,5,10,6},18/1),({10,5,9,7},16/1),({10,5,8,8},6/1),({9,6,14,2},1/1),({9,6,13,3},4/1),({9,6,12,4},10/1),({9,6,11,5},16/1),({9,6,10,6},20/1),({9,6,9,7},16/1),({9,6,8,8},7/1),({8,7,14,2},1/1),({8,7,13,3},3/1),({8,7,12,4},7/1),({8,7,11,5},11/1),({8,7,10,6},13/1),({8,7,9,7},11/1),({8,7,8,8},4/1)}, (5,0) => {}, (3,2) => {}, (5,1) => {({15,3,12,7},1/1),({15,3,11,8},2/1),({15,3,10,9},1/1),({14,4,13,6},3/1),({14,4,12,7},5/1),({14,4,11,8},6/1),({14,4,10,9},5/1),({13,5,15,4},1/1),({13,5,14,5},3/1),({13,5,13,6},8/1),({13,5,12,7},14/1),({13,5,11,8},16/1),({13,5,10,9},11/1),({12,6,15,4},2/1),({12,6,14,5},8/1),({12,6,13,6},16/1),({12,6,12,7},25/1),({12,6,11,8},27/1),({12,6,10,9},17/1),({11,7,16,3},1/1),({11,7,15,4},4/1),({11,7,14,5},11/1),({11,7,13,6},23/1),({11,7,12,7},32/1),({11,7,11,8},33/1),({11,7,10,9},22/1),({10,8,16,3},1/1),({10,8,15,4},4/1),({10,8,14,5},11/1),({10,8,13,6},20/1),({10,8,12,7},28/1),({10,8,11,8},28/1),({10,8,10,9},18/1),({9,9,15,4},2/1),({9,9,14,5},5/1),({9,9,13,6},8/1),({9,9,12,7},11/1),({9,9,11,8},12/1),({9,9,10,9},6/1)}, (6,0) => {}, (4,2) => {}};
end;