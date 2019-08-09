A := QQ[t_0,t_1,t_2,t_3];
--tb stands for Total Betti numbers
tb0117 = new HashTable from {(6,1) => 17160/1, (7,0) => 0, (5,2) => 0, (7,1) => 18018/1, (8,0) => 0, (6,2) => 0, (8,1) => 14014/1, (9,0) => 0, (7,2) => 0, (9,1) => 8008/1, (10,0) => 0, (8,2) => 0, (10,1) => 3276/1, (11,0) => 0, (9,2) => 0, (12,0) => 0, (11,1) => 910/1, (10,2) => 0, (13,0) => 0, (12,1) => 154/1, (11,2) => 0, (13,1) => 12/1, (12,2) => 0, (13,2) => 0, (14,1) => 0, (15,1) => 0, (0,0) => 2/1, (1,0) => 14/1, (1,1) => 0, (2,0) => 0, (2,1) => 364/1, (3,0) => 0, (3,1) => 2002/1, (4,0) => 0, (2,2) => 0, (4,1) => 6006/1, (5,0) => 0, (3,2) => 0, (5,1) => 12012/1, (6,0) => 0, (4,2) => 0};
--sb represents the betti numbers as sums of Schur functors
sb0117 = new HashTable from {(6,1) => {({7,0,27,23},1/1),({6,1,33,17},1/1),({6,1,32,18},2/1),({6,1,31,19},3/1),({6,1,30,20},4/1),({6,1,29,21},5/1),({6,1,28,22},5/1),({6,1,27,23},4/1),({6,1,26,24},3/1),({6,1,25,25},1/1),({5,2,37,13},1/1),({5,2,36,14},2/1),({5,2,35,15},5/1),({5,2,34,16},8/1),({5,2,33,17},13/1),({5,2,32,18},17/1),({5,2,31,19},23/1),({5,2,30,20},25/1),({5,2,29,21},28/1),({5,2,28,22},25/1),({5,2,27,23},22/1),({5,2,26,24},13/1),({5,2,25,25},6/1),({4,3,39,11},1/1),({4,3,38,12},2/1),({4,3,37,13},4/1),({4,3,36,14},8/1),({4,3,35,15},13/1),({4,3,34,16},19/1),({4,3,33,17},27/1),({4,3,32,18},34/1),({4,3,31,19},40/1),({4,3,30,20},45/1),({4,3,29,21},45/1),({4,3,28,22},41/1),({4,3,27,23},34/1),({4,3,26,24},22/1),({4,3,25,25},7/1)}, (7,0) => {}, (5,2) => {}, (7,1) => {({7,1,34,23},1/1),({7,1,33,24},1/1),({7,1,32,25},1/1),({7,1,31,26},1/1),({7,1,30,27},1/1),({7,1,29,28},1/1),({6,2,39,18},1/1),({6,2,38,19},2/1),({6,2,37,20},4/1),({6,2,36,21},6/1),({6,2,35,22},9/1),({6,2,34,23},11/1),({6,2,33,24},13/1),({6,2,32,25},13/1),({6,2,31,26},12/1),({6,2,30,27},9/1),({6,2,29,28},5/1),({5,3,42,15},1/1),({5,3,41,16},2/1),({5,3,40,17},5/1),({5,3,39,18},8/1),({5,3,38,19},14/1),({5,3,37,20},20/1),({5,3,36,21},28/1),({5,3,35,22},34/1),({5,3,34,23},40/1),({5,3,33,24},42/1),({5,3,32,25},42/1),({5,3,31,26},36/1),({5,3,30,27},27/1),({5,3,29,28},14/1),({4,4,43,14},1/1),({4,4,42,15},1/1),({4,4,41,16},3/1),({4,4,40,17},5/1),({4,4,39,18},9/1),({4,4,38,19},12/1),({4,4,37,20},18/1),({4,4,36,21},22/1),({4,4,35,22},28/1),({4,4,34,23},30/1),({4,4,33,24},32/1),({4,4,32,25},30/1),({4,4,31,26},27/1),({4,4,30,27},19/1),({4,4,29,28},10/1)}, (8,0) => {}, (6,2) => {}, (8,1) => {({7,2,40,24},1/1),({7,2,39,25},1/1),({7,2,38,26},2/1),({7,2,37,27},2/1),({7,2,36,28},3/1),({7,2,35,29},3/1),({7,2,34,30},3/1),({7,2,33,31},1/1),({7,2,32,32},1/1),({6,3,44,20},1/1),({6,3,43,21},2/1),({6,3,42,22},4/1),({6,3,41,23},7/1),({6,3,40,24},10/1),({6,3,39,25},14/1),({6,3,38,26},18/1),({6,3,37,27},20/1),({6,3,36,28},21/1),({6,3,35,29},20/1),({6,3,34,30},16/1),({6,3,33,31},11/1),({6,3,32,32},4/1),({5,4,46,18},1/1),({5,4,45,19},2/1),({5,4,44,20},4/1),({5,4,43,21},7/1),({5,4,42,22},12/1),({5,4,41,23},17/1),({5,4,40,24},24/1),({5,4,39,25},29/1),({5,4,38,26},35/1),({5,4,37,27},38/1),({5,4,36,28},39/1),({5,4,35,29},34/1),({5,4,34,30},29/1),({5,4,33,31},18/1),({5,4,32,32},7/1)}, (9,0) => {}, (7,2) => {}, (9,1) => {({7,3,45,26},1/1),({7,3,44,27},1/1),({7,3,43,28},2/1),({7,3,42,29},3/1),({7,3,41,30},4/1),({7,3,40,31},4/1),({7,3,39,32},5/1),({7,3,38,33},4/1),({7,3,37,34},3/1),({7,3,36,35},2/1),({6,4,48,23},1/1),({6,4,47,24},2/1),({6,4,46,25},4/1),({6,4,45,26},6/1),({6,4,44,27},10/1),({6,4,43,28},13/1),({6,4,42,29},17/1),({6,4,41,30},19/1),({6,4,40,31},21/1),({6,4,39,32},20/1),({6,4,38,33},18/1),({6,4,37,34},13/1),({6,4,36,35},7/1),({5,5,49,22},1/1),({5,5,48,23},1/1),({5,5,47,24},3/1),({5,5,46,25},4/1),({5,5,45,26},7/1),({5,5,44,27},9/1),({5,5,43,28},13/1),({5,5,42,29},14/1),({5,5,41,30},17/1),({5,5,40,31},17/1),({5,5,39,32},17/1),({5,5,38,33},14/1),({5,5,37,34},11/1),({5,5,36,35},5/1)}, (10,0) => {}, (8,2) => {}, (10,1) => {({7,4,49,29},1/1),({7,4,48,30},1/1),({7,4,47,31},2/1),({7,4,46,32},3/1),({7,4,45,33},4/1),({7,4,44,34},4/1),({7,4,43,35},5/1),({7,4,42,36},4/1),({7,4,41,37},4/1),({7,4,40,38},2/1),({7,4,39,39},1/1),({6,5,51,27},1/1),({6,5,50,28},2/1),({6,5,49,29},3/1),({6,5,48,30},5/1),({6,5,47,31},7/1),({6,5,46,32},9/1),({6,5,45,33},11/1),({6,5,44,34},12/1),({6,5,43,35},12/1),({6,5,42,36},11/1),({6,5,41,37},9/1),({6,5,40,38},6/1),({6,5,39,39},2/1)}, (11,0) => {}, (9,2) => {}, (12,0) => {}, (11,1) => {({7,5,52,33},1/1),({7,5,51,34},1/1),({7,5,50,35},2/1),({7,5,49,36},2/1),({7,5,48,37},3/1),({7,5,47,38},3/1),({7,5,46,39},3/1),({7,5,45,40},2/1),({7,5,44,41},2/1),({7,5,43,42},1/1),({6,6,53,32},1/1),({6,6,52,33},1/1),({6,6,51,34},2/1),({6,6,50,35},2/1),({6,6,49,36},3/1),({6,6,48,37},3/1),({6,6,47,38},4/1),({6,6,46,39},3/1),({6,6,45,40},3/1),({6,6,44,41},2/1),({6,6,43,42},1/1)}, (10,2) => {}, (13,0) => {}, (12,1) => {({7,6,54,38},1/1),({7,6,53,39},1/1),({7,6,52,40},1/1),({7,6,51,41},1/1),({7,6,50,42},1/1),({7,6,49,43},1/1),({7,6,48,44},1/1)}, (11,2) => {}, (13,1) => {({7,7,55,44},1/1)}, (12,2) => {}, (13,2) => {}, (14,1) => {}, (15,1) => {}, (0,0) => {({0,0,1,0},1/1)}, (1,0) => {({1,0,7,1},1/1)}, (1,1) => {}, (2,0) => {}, (2,1) => {({3,0,17,5},1/1),({3,0,15,7},1/1),({3,0,14,8},1/1),({3,0,13,9},1/1),({3,0,11,11},1/1),({2,1,19,3},1/1),({2,1,18,4},1/1),({2,1,17,5},1/1),({2,1,16,6},2/1),({2,1,15,7},2/1),({2,1,14,8},2/1),({2,1,13,9},2/1),({2,1,12,10},1/1)}, (3,0) => {}, (3,1) => {({4,0,21,8},1/1),({4,0,20,9},1/1),({4,0,19,10},1/1),({4,0,18,11},2/1),({4,0,17,12},2/1),({4,0,16,13},1/1),({4,0,15,14},1/1),({3,1,24,5},1/1),({3,1,23,6},2/1),({3,1,22,7},3/1),({3,1,21,8},4/1),({3,1,20,9},6/1),({3,1,19,10},7/1),({3,1,18,11},7/1),({3,1,17,12},6/1),({3,1,16,13},5/1),({3,1,15,14},3/1),({2,2,25,4},1/1),({2,2,24,5},1/1),({2,2,23,6},2/1),({2,2,22,7},3/1),({2,2,21,8},4/1),({2,2,20,9},5/1),({2,2,19,10},6/1),({2,2,18,11},5/1),({2,2,17,12},5/1),({2,2,16,13},4/1),({2,2,15,14},2/1)}, (4,0) => {}, (2,2) => {}, (4,1) => {({5,0,24,12},1/1),({5,0,23,13},1/1),({5,0,22,14},2/1),({5,0,21,15},1/1),({5,0,20,16},2/1),({5,0,19,17},1/1),({5,0,18,18},1/1),({4,1,28,8},1/1),({4,1,27,9},2/1),({4,1,26,10},4/1),({4,1,25,11},6/1),({4,1,24,12},8/1),({4,1,23,13},10/1),({4,1,22,14},11/1),({4,1,21,15},11/1),({4,1,20,16},9/1),({4,1,19,17},6/1),({4,1,18,18},2/1),({3,2,30,6},1/1),({3,2,29,7},2/1),({3,2,28,8},4/1),({3,2,27,9},6/1),({3,2,26,10},10/1),({3,2,25,11},13/1),({3,2,24,12},17/1),({3,2,23,13},18/1),({3,2,22,14},20/1),({3,2,21,15},18/1),({3,2,20,16},16/1),({3,2,19,17},9/1),({3,2,18,18},4/1)}, (5,0) => {}, (3,2) => {}, (5,1) => {({6,0,26,17},1/1),({6,0,25,18},1/1),({6,0,24,19},1/1),({6,0,23,20},1/1),({5,1,31,12},1/1),({5,1,30,13},2/1),({5,1,29,14},4/1),({5,1,28,15},6/1),({5,1,27,16},8/1),({5,1,26,17},9/1),({5,1,25,18},10/1),({5,1,24,19},9/1),({5,1,23,20},7/1),({5,1,22,21},4/1),({4,2,34,9},1/1),({4,2,33,10},2/1),({4,2,32,11},5/1),({4,2,31,12},8/1),({4,2,30,13},13/1),({4,2,29,14},18/1),({4,2,28,15},24/1),({4,2,27,16},28/1),({4,2,26,17},31/1),({4,2,25,18},30/1),({4,2,24,19},27/1),({4,2,23,20},20/1),({4,2,22,21},11/1),({3,3,35,8},1/1),({3,3,34,9},1/1),({3,3,33,10},3/1),({3,3,32,11},5/1),({3,3,31,12},8/1),({3,3,30,13},11/1),({3,3,29,14},16/1),({3,3,28,15},18/1),({3,3,27,16},22/1),({3,3,26,17},23/1),({3,3,25,18},22/1),({3,3,24,19},19/1),({3,3,23,20},15/1),({3,3,22,21},7/1)}, (6,0) => {}, (4,2) => {}};
end;