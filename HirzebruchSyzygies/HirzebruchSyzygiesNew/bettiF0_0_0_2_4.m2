A := QQ[t_0,t_1,t_2,t_3];
--tb stands for Total Betti numbers
tb0024 = new HashTable from {(6,1) => 7920/1, (7,0) => 0, (5,2) => 0, (7,1) => 6237/1, (8,0) => 0, (6,2) => 0, (8,1) => 3344/1, (9,0) => 0, (7,2) => 0, (9,1) => 1089/1, (10,0) => 0, (8,2) => 0, (10,1) => 120/1, (11,0) => 0, (9,2) => 0, (11,1) => 11/1, (12,0) => 0, (10,2) => 66/1, (12,1) => 0, (11,2) => 24/1, (12,2) => 3/1, (13,2) => 0, (0,0) => 1/1, (1,0) => 0, (1,1) => 75/1, (2,0) => 0, (2,1) => 536/1, (3,0) => 0, (3,1) => 1947/1, (4,0) => 0, (2,2) => 0, (4,1) => 4488/1, (5,0) => 0, (3,2) => 0, (5,1) => 7095/1, (6,0) => 0, (4,2) => 0};
--sb represents the betti numbers as sums of Schur functors
sb0024 = new HashTable from {(6,1) => {({11,3,18,10},1/1),({11,3,17,11},2/1),({11,3,16,12},2/1),({11,3,15,13},2/1),({11,3,14,14},1/1),({10,4,20,8},1/1),({10,4,19,9},3/1),({10,4,18,10},6/1),({10,4,17,11},9/1),({10,4,16,12},10/1),({10,4,15,13},8/1),({10,4,14,14},3/1),({9,5,21,7},1/1),({9,5,20,8},4/1),({9,5,19,9},9/1),({9,5,18,10},15/1),({9,5,17,11},20/1),({9,5,16,12},21/1),({9,5,15,13},16/1),({9,5,14,14},6/1),({8,6,22,6},1/1),({8,6,21,7},3/1),({8,6,20,8},7/1),({8,6,19,9},13/1),({8,6,18,10},20/1),({8,6,17,11},25/1),({8,6,16,12},25/1),({8,6,15,13},19/1),({8,6,14,14},7/1),({7,7,21,7},1/1),({7,7,20,8},3/1),({7,7,19,9},6/1),({7,7,18,10},9/1),({7,7,17,11},11/1),({7,7,16,12},11/1),({7,7,15,13},8/1),({7,7,14,14},3/1)}, (7,0) => {}, (5,2) => {}, (7,1) => {({12,4,19,13},2/1),({12,4,18,14},1/1),({12,4,17,15},2/1),({11,5,22,10},1/1),({11,5,21,11},2/1),({11,5,20,12},5/1),({11,5,19,13},6/1),({11,5,18,14},8/1),({11,5,17,15},6/1),({11,5,16,16},3/1),({10,6,23,9},1/1),({10,6,22,10},3/1),({10,6,21,11},8/1),({10,6,20,12},11/1),({10,6,19,13},17/1),({10,6,18,14},15/1),({10,6,17,15},14/1),({10,6,16,16},3/1),({9,7,24,8},1/1),({9,7,23,9},2/1),({9,7,22,10},6/1),({9,7,21,11},10/1),({9,7,20,12},17/1),({9,7,19,13},19/1),({9,7,18,14},21/1),({9,7,17,15},13/1),({9,7,16,16},7/1),({8,8,23,9},2/1),({8,8,22,10},2/1),({8,8,21,11},6/1),({8,8,20,12},6/1),({8,8,19,13},11/1),({8,8,18,14},7/1),({8,8,17,15},9/1)}, (8,0) => {}, (6,2) => {}, (8,1) => {({13,5,20,16},1/1),({13,5,19,17},1/1),({12,6,23,13},1/1),({12,6,22,14},2/1),({12,6,21,15},3/1),({12,6,20,16},4/1),({12,6,19,17},3/1),({12,6,18,18},1/1),({11,7,25,11},1/1),({11,7,24,12},2/1),({11,7,23,13},4/1),({11,7,22,14},7/1),({11,7,21,15},9/1),({11,7,20,16},9/1),({11,7,19,17},7/1),({11,7,18,18},3/1),({10,8,25,11},1/1),({10,8,24,12},3/1),({10,8,23,13},6/1),({10,8,22,14},9/1),({10,8,21,15},11/1),({10,8,20,16},11/1),({10,8,19,17},8/1),({10,8,18,18},3/1),({9,9,26,10},1/1),({9,9,25,11},1/1),({9,9,24,12},2/1),({9,9,23,13},4/1),({9,9,22,14},5/1),({9,9,21,15},6/1),({9,9,20,16},6/1),({9,9,19,17},4/1),({9,9,18,18},1/1)}, (9,0) => {}, (7,2) => {}, (9,1) => {({14,6,20,20},1/1),({13,7,23,17},1/1),({13,7,22,18},1/1),({13,7,21,19},1/1),({12,8,26,14},1/1),({12,8,25,15},1/1),({12,8,24,16},3/1),({12,8,23,17},2/1),({12,8,22,18},4/1),({12,8,21,19},1/1),({12,8,20,20},2/1),({11,9,27,13},1/1),({11,9,26,14},1/1),({11,9,25,15},3/1),({11,9,24,16},3/1),({11,9,23,17},5/1),({11,9,22,18},3/1),({11,9,21,19},4/1),({10,10,26,14},1/1),({10,10,25,15},1/1),({10,10,24,16},3/1),({10,10,23,17},1/1),({10,10,22,18},3/1),({10,10,20,20},2/1)}, (10,0) => {}, (8,2) => {}, (10,1) => {({12,10,28,16},1/1),({12,10,27,17},1/1),({12,10,26,18},1/1),({12,10,25,19},1/1)}, (11,0) => {}, (9,2) => {}, (11,1) => {({12,12,29,19},1/1)}, (12,0) => {}, (10,2) => {({14,10,26,22},1/1),({14,10,24,24},1/1),({13,11,27,21},1/1),({13,11,25,23},1/1),({12,12,26,22},1/1),({12,12,24,24},1/1)}, (12,1) => {}, (11,2) => {({14,12,28,24},1/1),({14,12,27,25},1/1)}, (12,2) => {({14,14,29,27},1/1)}, (13,2) => {}, (0,0) => {({0,0,0,0},1/1)}, (1,0) => {}, (1,1) => {({4,0,6,2},1/1),({4,0,4,4},1/1),({3,1,7,1},1/1),({3,1,5,3},1/1),({2,2,8,0},1/1),({2,2,6,2},1/1),({2,2,4,4},1/1)}, (2,0) => {}, (2,1) => {({6,0,8,4},1/1),({6,0,7,5},1/1),({5,1,10,2},1/1),({5,1,9,3},2/1),({5,1,8,4},2/1),({5,1,7,5},2/1),({5,1,6,6},1/1),({4,2,11,1},1/1),({4,2,10,2},2/1),({4,2,9,3},3/1),({4,2,8,4},4/1),({4,2,7,5},3/1),({4,2,6,6},1/1),({3,3,11,1},1/1),({3,3,10,2},1/1),({3,3,9,3},1/1),({3,3,8,4},2/1),({3,3,7,5},1/1)}, (3,0) => {}, (3,1) => {({8,0,9,7},1/1),({7,1,12,4},1/1),({7,1,11,5},2/1),({7,1,10,6},3/1),({7,1,9,7},2/1),({7,1,8,8},1/1),({6,2,13,3},2/1),({6,2,12,4},3/1),({6,2,11,5},7/1),({6,2,10,6},6/1),({6,2,9,7},7/1),({6,2,8,8},1/1),({5,3,14,2},2/1),({5,3,13,3},3/1),({5,3,12,4},7/1),({5,3,11,5},8/1),({5,3,10,6},10/1),({5,3,9,7},6/1),({5,3,8,8},4/1),({4,4,13,3},2/1),({4,4,12,4},2/1),({4,4,11,5},5/1),({4,4,10,6},3/1),({4,4,9,7},5/1)}, (4,0) => {}, (2,2) => {}, (4,1) => {({9,1,13,7},1/1),({9,1,12,8},1/1),({9,1,11,9},1/1),({9,1,10,10},1/1),({8,2,15,5},1/1),({8,2,14,6},3/1),({8,2,13,7},5/1),({8,2,12,8},6/1),({8,2,11,9},5/1),({8,2,10,10},2/1),({7,3,16,4},2/1),({7,3,15,5},5/1),({7,3,14,6},9/1),({7,3,13,7},13/1),({7,3,12,8},14/1),({7,3,11,9},11/1),({7,3,10,10},4/1),({6,4,17,3},1/1),({6,4,16,4},3/1),({6,4,15,5},7/1),({6,4,14,6},12/1),({6,4,13,7},16/1),({6,4,12,8},17/1),({6,4,11,9},13/1),({6,4,10,10},5/1),({5,5,17,3},1/1),({5,5,16,4},2/1),({5,5,15,5},4/1),({5,5,14,6},7/1),({5,5,13,7},8/1),({5,5,12,8},8/1),({5,5,11,9},6/1),({5,5,10,10},2/1)}, (5,0) => {}, (3,2) => {}, (5,1) => {({10,2,16,8},1/1),({10,2,15,9},1/1),({10,2,14,10},3/1),({10,2,13,11},1/1),({10,2,12,12},1/1),({9,3,17,7},3/1),({9,3,16,8},5/1),({9,3,15,9},8/1),({9,3,14,10},9/1),({9,3,13,11},8/1),({9,3,12,12},2/1),({8,4,19,5},1/1),({8,4,18,6},4/1),({8,4,17,7},7/1),({8,4,16,8},15/1),({8,4,15,9},18/1),({8,4,14,10},21/1),({8,4,13,11},14/1),({8,4,12,12},8/1),({7,5,19,5},2/1),({7,5,18,6},5/1),({7,5,17,7},12/1),({7,5,16,8},17/1),({7,5,15,9},24/1),({7,5,14,10},22/1),({7,5,13,11},19/1),({7,5,12,12},5/1),({6,6,20,4},1/1),({6,6,19,5},1/1),({6,6,18,6},4/1),({6,6,17,7},5/1),({6,6,16,8},10/1),({6,6,15,9},9/1),({6,6,14,10},14/1),({6,6,13,11},6/1),({6,6,12,12},5/1)}, (6,0) => {}, (4,2) => {}};
--dw stands for dominant weights
dw0024 = new HashTable from {(5,2) => {}, (7,0) => {}, (6,1) => {({11,3,18,10},1/1),({10,4,20,8},1/1),({9,5,21,7},1/1),({8,6,22,6},1/1)}, (7,1) => {({12,4,19,13},2/1),({11,5,22,10},1/1),({10,6,23,9},1/1),({9,7,24,8},1/1)}, (8,0) => {}, (6,2) => {}, (7,2) => {}, (9,0) => {}, (8,1) => {({13,5,20,16},1/1),({12,6,23,13},1/1),({11,7,25,11},1/1),({9,9,26,10},1/1)}, (8,2) => {}, (10,0) => {}, (9,1) => {({14,6,20,20},1/1),({13,7,23,17},1/1),({12,8,26,14},1/1),({11,9,27,13},1/1)}, (10,1) => {({12,10,28,16},1/1)}, (11,0) => {}, (9,2) => {}, (11,1) => {({12,12,29,19},1/1)}, (12,0) => {}, (10,2) => {({14,10,26,22},1/1),({13,11,27,21},1/1)}, (12,1) => {}, (11,2) => {({14,12,28,24},1/1)}, (12,2) => {({14,14,29,27},1/1)}, (13,2) => {}, (0,0) => {({0,0,0,0},1/1)}, (1,0) => {}, (2,0) => {}, (1,1) => {({4,0,6,2},1/1),({3,1,7,1},1/1),({2,2,8,0},1/1)}, (3,0) => {}, (2,1) => {({6,0,8,4},1/1),({5,1,10,2},1/1),({4,2,11,1},1/1)}, (2,2) => {}, (4,0) => {}, (3,1) => {({8,0,9,7},1/1),({7,1,12,4},1/1),({6,2,13,3},2/1),({5,3,14,2},2/1)}, (3,2) => {}, (5,0) => {}, (4,1) => {({9,1,13,7},1/1),({8,2,15,5},1/1),({7,3,16,4},2/1),({6,4,17,3},1/1)}, (4,2) => {}, (6,0) => {}, (5,1) => {({10,2,16,8},1/1),({9,3,17,7},3/1),({8,4,19,5},1/1),({6,6,20,4},1/1)}};
--dmw stands for dominant monomial weights
dmw0024 = new HashTable from {(5,2) => {}, (7,0) => {}, (6,1) => {{10,4,19,9},{9,5,20,8},{11,3,18,10},{12,2,16,12}}, (7,1) => {{13,3,18,14},{12,4,19,13},{9,7,22,10},{11,5,21,11}}, (8,0) => {}, (6,2) => {}, (7,2) => {}, (9,0) => {}, (8,1) => {{14,4,19,17},{13,5,21,15},{11,7,23,13},{12,6,22,14}}, (8,2) => {{15,5,15,25}}, (10,0) => {}, (9,1) => {{14,6,22,18},{12,8,24,16},{13,7,23,17}}, (10,1) => {{14,8,24,20},{12,10,25,19}}, (11,0) => {}, (9,2) => {{16,6,19,25}}, (11,1) => {{14,10,25,23}}, (12,0) => {}, (10,2) => {{16,8,22,26}}, (12,1) => {}, (11,2) => {{16,10,24,28}}, (12,2) => {{16,12,25,31}}, (13,2) => {}, (0,0) => {}, (1,0) => {}, (2,0) => {}, (1,1) => {{4,0,6,2}}, (3,0) => {}, (2,1) => {{5,1,10,2},{6,0,8,4}}, (2,2) => {}, (4,0) => {}, (3,1) => {{7,1,12,4},{6,2,13,3},{8,0,9,7}}, (3,2) => {}, (5,0) => {}, (4,1) => {{10,0,9,11},{6,4,16,4},{8,2,15,5},{9,1,13,7}}, (4,2) => {}, (6,0) => {}, (5,1) => {{11,1,13,11},{9,3,17,7},{10,2,16,8},{8,4,18,6}}};
end;