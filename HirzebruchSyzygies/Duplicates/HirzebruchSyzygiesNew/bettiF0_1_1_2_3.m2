A := QQ[t_0,t_1,t_2,t_3];
--tb stands for Total Betti numbers
tb1123 = new HashTable from {(5,2) => 0, (6,1) => 528/1, (7,0) => 0, (6,2) => 0, (7,1) => 252/1, (8,0) => 0, (8,1) => 68/1, (9,0) => 0, (7,2) => 0, (9,1) => 8/1, (8,2) => 0, (9,2) => 0, (10,1) => 0, (11,1) => 0, (0,0) => 4/1, (1,0) => 28/1, (2,0) => 72/1, (1,1) => 0, (3,0) => 56/1, (2,1) => 8/1, (2,2) => 0, (3,1) => 168/1, (4,0) => 0, (3,2) => 0, (4,1) => 504/1, (5,0) => 0, (4,2) => 0, (5,1) => 672/1, (6,0) => 0};
--sb represents the betti numbers as sums of Schur functors
sb1123 = new HashTable from {(5,2) => {}, (6,1) => {({11,4,11,11},1/1),({10,5,14,8},1/1),({10,5,13,9},2/1),({10,5,12,10},2/1),({10,5,11,11},1/1),({9,6,15,7},1/1),({9,6,14,8},2/1),({9,6,13,9},5/1),({9,6,12,10},4/1),({9,6,11,11},2/1),({8,7,15,7},1/1),({8,7,14,8},3/1),({8,7,13,9},4/1),({8,7,12,10},4/1),({8,7,11,11},2/1)}, (7,0) => {}, (6,2) => {}, (7,1) => {({11,6,14,11},1/1),({11,6,13,12},1/1),({10,7,16,9},1/1),({10,7,15,10},2/1),({10,7,14,11},3/1),({10,7,13,12},2/1),({9,8,16,9},1/1),({9,8,15,10},2/1),({9,8,14,11},3/1),({9,8,13,12},2/1)}, (8,0) => {}, (8,1) => {({11,8,16,12},1/1),({11,8,15,13},1/1),({11,8,14,14},1/1),({10,9,17,11},1/1),({10,9,16,12},1/1),({10,9,15,13},1/1),({10,9,14,14},1/1)}, (9,0) => {}, (7,2) => {}, (9,1) => {({11,10,17,14},1/1)}, (8,2) => {}, (9,2) => {}, (10,1) => {}, (11,1) => {}, (0,0) => {({1,0,1,0},1/1)}, (1,0) => {({3,0,3,1},1/1),({2,1,4,0},1/1),({2,1,3,1},1/1)}, (2,0) => {({4,1,6,1},1/1),({4,1,5,2},1/1),({4,1,4,3},1/1),({3,2,6,1},1/1),({3,2,5,2},1/1),({3,2,4,3},1/1)}, (1,1) => {}, (3,0) => {({5,2,8,2},1/1),({5,2,6,4},1/1),({4,3,7,3},1/1),({4,3,6,4},1/1)}, (2,1) => {({7,0,5,5},1/1)}, (2,2) => {}, (3,1) => {({8,1,8,5},1/1),({7,2,8,5},1/1),({7,2,7,6},1/1),({6,3,10,3},1/1),({6,3,9,4},1/1),({6,3,8,5},1/1),({6,3,7,6},1/1),({5,4,9,4},1/1),({5,4,8,5},1/1)}, (4,0) => {}, (3,2) => {}, (4,1) => {({9,2,10,6},1/1),({9,2,8,8},1/1),({8,3,11,5},1/1),({8,3,10,6},2/1),({8,3,9,7},2/1),({8,3,8,8},1/1),({7,4,12,4},1/1),({7,4,11,5},2/1),({7,4,10,6},4/1),({7,4,9,7},3/1),({7,4,8,8},2/1),({6,5,12,4},1/1),({6,5,11,5},2/1),({6,5,10,6},3/1),({6,5,9,7},3/1),({6,5,8,8},1/1)}, (5,0) => {}, (4,2) => {}, (5,1) => {({10,3,11,8},1/1),({9,4,13,6},1/1),({9,4,12,7},2/1),({9,4,11,8},3/1),({9,4,10,9},2/1),({8,5,13,6},2/1),({8,5,12,7},4/1),({8,5,11,8},5/1),({8,5,10,9},4/1),({7,6,14,5},1/1),({7,6,13,6},2/1),({7,6,12,7},4/1),({7,6,11,8},5/1),({7,6,10,9},3/1)}, (6,0) => {}};
--dw stands for dominant weights
dw1123 = new HashTable from {(7,0) => {}, (6,1) => {({11,4,11,11},1/1),({10,5,14,8},1/1),({9,6,15,7},1/1)}, (5,2) => {}, (6,2) => {}, (7,1) => {({11,6,14,11},1/1),({10,7,16,9},1/1)}, (8,0) => {}, (7,2) => {}, (9,0) => {}, (8,1) => {({11,8,16,12},1/1),({10,9,17,11},1/1)}, (8,2) => {}, (9,1) => {({11,10,17,14},1/1)}, (9,2) => {}, (10,1) => {}, (11,1) => {}, (0,0) => {({1,0,1,0},1/1)}, (1,0) => {({3,0,3,1},1/1),({2,1,4,0},1/1)}, (2,0) => {({4,1,6,1},1/1)}, (1,1) => {}, (3,0) => {({5,2,8,2},1/1)}, (2,1) => {({7,0,5,5},1/1)}, (2,2) => {}, (3,1) => {({8,1,8,5},1/1),({6,3,10,3},1/1)}, (4,0) => {}, (5,0) => {}, (4,1) => {({9,2,10,6},1/1),({8,3,11,5},1/1),({7,4,12,4},1/1)}, (3,2) => {}, (6,0) => {}, (5,1) => {({10,3,11,8},1/1),({9,4,13,6},1/1),({7,6,14,5},1/1)}, (4,2) => {}};
--dmw stands for dominant monomial weights
dmw1123 = new HashTable from {(7,0) => {}, (6,1) => {{12,3,11,11},{10,5,13,9},{11,4,12,10}}, (5,2) => {}, (6,2) => {}, (7,1) => {{13,4,11,14},{12,5,13,12},{10,7,14,11}}, (8,0) => {}, (7,2) => {}, (9,0) => {}, (8,1) => {{12,7,14,14},{13,6,13,15}}, (8,2) => {}, (9,1) => {{13,8,14,17}}, (9,2) => {{12,9,11,20}}, (10,1) => {}, (11,1) => {}, (0,0) => {}, (1,0) => {}, (2,0) => {}, (1,1) => {}, (3,0) => {}, (2,1) => {{7,0,5,5}}, (2,2) => {}, (3,1) => {{9,0,5,8},{8,1,8,5}}, (4,0) => {}, (5,0) => {}, (4,1) => {{10,1,8,8},{9,2,10,6}}, (3,2) => {}, (6,0) => {}, (5,1) => {{11,2,10,9},{9,4,12,7},{10,3,11,8}}, (4,2) => {}};
end;