A := QQ[t_0,t_1,t_2,t_3];
--tb stands for Total Betti numbers
tb2233 = new HashTable from {(5,2) => 0, (7,0) => 1080/1, (6,1) => 1080/1, (7,1) => 4026/1, (8,0) => 165/1, (6,2) => 0, (7,2) => 0, (9,0) => 0, (8,1) => 5148/1, (8,2) => 0, (10,0) => 0, (9,1) => 3861/1, (10,1) => 1872/1, (11,0) => 0, (9,2) => 0, (12,0) => 0, (11,1) => 585/1, (10,2) => 0, (13,0) => 0, (12,1) => 108/1, (11,2) => 0, (13,1) => 9/1, (12,2) => 0, (14,1) => 0, (13,2) => 0, (15,1) => 0, (0,0) => 9/1, (1,0) => 108/1, (2,0) => 585/1, (1,1) => 0, (3,0) => 1872/1, (2,1) => 0, (2,2) => 0, (4,0) => 3861/1, (3,1) => 0, (3,2) => 0, (5,0) => 5148/1, (4,1) => 0, (4,2) => 0, (6,0) => 4026/1, (5,1) => 165/1};
--sb represents the betti numbers as sums of Schur functors
sb2233 = new HashTable from {(5,2) => {}, (7,0) => {({16,7,16,7},1/1),({16,7,15,8},1/1),({16,7,14,9},1/1),({15,8,16,7},1/1),({15,8,15,8},2/1),({15,8,14,9},2/1),({15,8,13,10},1/1),({14,9,16,7},1/1),({14,9,15,8},2/1),({14,9,14,9},3/1),({14,9,13,10},2/1),({14,9,12,11},1/1),({13,10,15,8},1/1),({13,10,14,9},2/1),({13,10,13,10},3/1),({13,10,12,11},2/1),({12,11,14,9},1/1),({12,11,13,10},2/1),({12,11,12,11},2/1)}, (6,1) => {({16,7,16,7},1/1),({16,7,15,8},1/1),({16,7,14,9},1/1),({15,8,16,7},1/1),({15,8,15,8},2/1),({15,8,14,9},2/1),({15,8,13,10},1/1),({14,9,16,7},1/1),({14,9,15,8},2/1),({14,9,14,9},3/1),({14,9,13,10},2/1),({14,9,12,11},1/1),({13,10,15,8},1/1),({13,10,14,9},2/1),({13,10,13,10},3/1),({13,10,12,11},2/1),({12,11,14,9},1/1),({12,11,13,10},2/1),({12,11,12,11},2/1)}, (7,1) => {({19,7,14,12},1/1),({18,8,17,9},1/1),({18,8,16,10},1/1),({18,8,15,11},2/1),({18,8,14,12},1/1),({18,8,13,13},1/1),({17,9,18,8},1/1),({17,9,17,9},2/1),({17,9,16,10},4/1),({17,9,15,11},4/1),({17,9,14,12},5/1),({17,9,13,13},1/1),({16,10,18,8},1/1),({16,10,17,9},4/1),({16,10,16,10},7/1),({16,10,15,11},9/1),({16,10,14,12},6/1),({16,10,13,13},4/1),({15,11,18,8},2/1),({15,11,17,9},4/1),({15,11,16,10},9/1),({15,11,15,11},10/1),({15,11,14,12},11/1),({15,11,13,13},2/1),({14,12,19,7},1/1),({14,12,18,8},1/1),({14,12,17,9},5/1),({14,12,16,10},6/1),({14,12,15,11},11/1),({14,12,14,12},7/1),({14,12,13,13},5/1),({13,13,18,8},1/1),({13,13,17,9},1/1),({13,13,16,10},4/1),({13,13,15,11},2/1),({13,13,14,12},5/1)}, (8,0) => {({17,9,17,9},1/1),({16,10,16,10},1/1),({15,11,15,11},1/1),({14,12,14,12},1/1),({13,13,13,13},1/1)}, (6,2) => {}, (7,2) => {}, (9,0) => {}, (8,1) => {({20,9,17,12},1/1),({20,9,16,13},1/1),({20,9,15,14},1/1),({19,10,19,10},1/1),({19,10,18,11},2/1),({19,10,17,12},3/1),({19,10,16,13},4/1),({19,10,15,14},3/1),({18,11,19,10},2/1),({18,11,18,11},5/1),({18,11,17,12},8/1),({18,11,16,13},9/1),({18,11,15,14},6/1),({17,12,20,9},1/1),({17,12,19,10},3/1),({17,12,18,11},8/1),({17,12,17,12},13/1),({17,12,16,13},14/1),({17,12,15,14},9/1),({16,13,20,9},1/1),({16,13,19,10},4/1),({16,13,18,11},9/1),({16,13,17,12},14/1),({16,13,16,13},15/1),({16,13,15,14},10/1),({15,14,20,9},1/1),({15,14,19,10},3/1),({15,14,18,11},6/1),({15,14,17,12},9/1),({15,14,16,13},10/1),({15,14,15,14},7/1)}, (8,2) => {}, (10,0) => {}, (9,1) => {({21,11,19,13},1/1),({21,11,18,14},1/1),({21,11,17,15},2/1),({20,12,20,12},1/1),({20,12,19,13},3/1),({20,12,18,14},5/1),({20,12,17,15},4/1),({20,12,16,16},2/1),({19,13,21,11},1/1),({19,13,20,12},3/1),({19,13,19,13},7/1),({19,13,18,14},9/1),({19,13,17,15},10/1),({19,13,16,16},3/1),({18,14,21,11},1/1),({18,14,20,12},5/1),({18,14,19,13},9/1),({18,14,18,14},14/1),({18,14,17,15},11/1),({18,14,16,16},6/1),({17,15,21,11},2/1),({17,15,20,12},4/1),({17,15,19,13},10/1),({17,15,18,14},11/1),({17,15,17,15},12/1),({17,15,16,16},3/1),({16,16,20,12},2/1),({16,16,19,13},3/1),({16,16,18,14},6/1),({16,16,17,15},3/1),({16,16,16,16},3/1)}, (10,1) => {({22,13,20,15},1/1),({22,13,19,16},1/1),({22,13,18,17},1/1),({21,14,21,14},1/1),({21,14,20,15},3/1),({21,14,19,16},4/1),({21,14,18,17},3/1),({20,15,22,13},1/1),({20,15,21,14},3/1),({20,15,20,15},6/1),({20,15,19,16},7/1),({20,15,18,17},5/1),({19,16,22,13},1/1),({19,16,21,14},4/1),({19,16,20,15},7/1),({19,16,19,16},9/1),({19,16,18,17},6/1),({18,17,22,13},1/1),({18,17,21,14},3/1),({18,17,20,15},5/1),({18,17,19,16},6/1),({18,17,18,17},4/1)}, (11,0) => {}, (9,2) => {}, (12,0) => {}, (11,1) => {({23,15,20,18},1/1),({22,16,22,16},1/1),({22,16,21,17},2/1),({22,16,20,18},2/1),({22,16,19,19},1/1),({21,17,22,16},2/1),({21,17,21,17},3/1),({21,17,20,18},4/1),({21,17,19,19},1/1),({20,18,23,15},1/1),({20,18,22,16},2/1),({20,18,21,17},4/1),({20,18,20,18},3/1),({20,18,19,19},2/1),({19,19,22,16},1/1),({19,19,21,17},1/1),({19,19,20,18},2/1)}, (10,2) => {}, (13,0) => {}, (12,1) => {({23,18,22,19},1/1),({23,18,21,20},1/1),({22,19,23,18},1/1),({22,19,22,19},1/1),({22,19,21,20},1/1),({21,20,23,18},1/1),({21,20,22,19},1/1),({21,20,21,20},1/1)}, (11,2) => {}, (13,1) => {({23,21,23,21},1/1)}, (12,2) => {}, (14,1) => {}, (13,2) => {}, (15,1) => {}, (0,0) => {({2,0,2,0},1/1)}, (1,0) => {({5,0,4,1},1/1),({5,0,3,2},1/1),({4,1,5,0},1/1),({4,1,4,1},1/1),({4,1,3,2},1/1),({3,2,5,0},1/1),({3,2,4,1},1/1),({3,2,3,2},1/1)}, (2,0) => {({8,0,5,3},1/1),({7,1,7,1},1/1),({7,1,6,2},2/1),({7,1,5,3},2/1),({7,1,4,4},1/1),({6,2,7,1},2/1),({6,2,6,2},3/1),({6,2,5,3},4/1),({6,2,4,4},1/1),({5,3,8,0},1/1),({5,3,7,1},2/1),({5,3,6,2},4/1),({5,3,5,3},3/1),({5,3,4,4},2/1),({4,4,7,1},1/1),({4,4,6,2},1/1),({4,4,5,3},2/1)}, (1,1) => {}, (3,0) => {({10,1,8,3},1/1),({10,1,7,4},1/1),({10,1,6,5},1/1),({9,2,9,2},1/1),({9,2,8,3},3/1),({9,2,7,4},4/1),({9,2,6,5},3/1),({8,3,10,1},1/1),({8,3,9,2},3/1),({8,3,8,3},6/1),({8,3,7,4},7/1),({8,3,6,5},5/1),({7,4,10,1},1/1),({7,4,9,2},4/1),({7,4,8,3},7/1),({7,4,7,4},9/1),({7,4,6,5},6/1),({6,5,10,1},1/1),({6,5,9,2},3/1),({6,5,8,3},5/1),({6,5,7,4},6/1),({6,5,6,5},4/1)}, (2,1) => {}, (2,2) => {}, (4,0) => {({12,2,10,4},1/1),({12,2,9,5},1/1),({12,2,8,6},2/1),({11,3,11,3},1/1),({11,3,10,4},3/1),({11,3,9,5},5/1),({11,3,8,6},4/1),({11,3,7,7},2/1),({10,4,12,2},1/1),({10,4,11,3},3/1),({10,4,10,4},7/1),({10,4,9,5},9/1),({10,4,8,6},10/1),({10,4,7,7},3/1),({9,5,12,2},1/1),({9,5,11,3},5/1),({9,5,10,4},9/1),({9,5,9,5},14/1),({9,5,8,6},11/1),({9,5,7,7},6/1),({8,6,12,2},2/1),({8,6,11,3},4/1),({8,6,10,4},10/1),({8,6,9,5},11/1),({8,6,8,6},12/1),({8,6,7,7},3/1),({7,7,11,3},2/1),({7,7,10,4},3/1),({7,7,9,5},6/1),({7,7,8,6},3/1),({7,7,7,7},3/1)}, (3,1) => {}, (3,2) => {}, (5,0) => {({14,3,11,6},1/1),({14,3,10,7},1/1),({14,3,9,8},1/1),({13,4,13,4},1/1),({13,4,12,5},2/1),({13,4,11,6},3/1),({13,4,10,7},4/1),({13,4,9,8},3/1),({12,5,13,4},2/1),({12,5,12,5},5/1),({12,5,11,6},8/1),({12,5,10,7},9/1),({12,5,9,8},6/1),({11,6,14,3},1/1),({11,6,13,4},3/1),({11,6,12,5},8/1),({11,6,11,6},13/1),({11,6,10,7},14/1),({11,6,9,8},9/1),({10,7,14,3},1/1),({10,7,13,4},4/1),({10,7,12,5},9/1),({10,7,11,6},14/1),({10,7,10,7},15/1),({10,7,9,8},10/1),({9,8,14,3},1/1),({9,8,13,4},3/1),({9,8,12,5},6/1),({9,8,11,6},9/1),({9,8,10,7},10/1),({9,8,9,8},7/1)}, (4,1) => {}, (4,2) => {}, (6,0) => {({16,4,11,9},1/1),({15,5,14,6},1/1),({15,5,13,7},1/1),({15,5,12,8},2/1),({15,5,11,9},1/1),({15,5,10,10},1/1),({14,6,15,5},1/1),({14,6,14,6},2/1),({14,6,13,7},4/1),({14,6,12,8},4/1),({14,6,11,9},5/1),({14,6,10,10},1/1),({13,7,15,5},1/1),({13,7,14,6},4/1),({13,7,13,7},7/1),({13,7,12,8},9/1),({13,7,11,9},6/1),({13,7,10,10},4/1),({12,8,15,5},2/1),({12,8,14,6},4/1),({12,8,13,7},9/1),({12,8,12,8},10/1),({12,8,11,9},11/1),({12,8,10,10},2/1),({11,9,16,4},1/1),({11,9,15,5},1/1),({11,9,14,6},5/1),({11,9,13,7},6/1),({11,9,12,8},11/1),({11,9,11,9},7/1),({11,9,10,10},5/1),({10,10,15,5},1/1),({10,10,14,6},1/1),({10,10,13,7},4/1),({10,10,12,8},2/1),({10,10,11,9},5/1)}, (5,1) => {({14,6,14,6},1/1),({13,7,13,7},1/1),({12,8,12,8},1/1),({11,9,11,9},1/1),({10,10,10,10},1/1)}};
end;
