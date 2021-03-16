A := QQ[t_0,t_1,t_2,t_3];
--tb stands for Total Betti numbers
tb1022 = new HashTable from {(6,1) => 6/1, (5,2) => 0, (7,1) => 0, (6,2) => 0, (8,1) => 0, (0,0) => 2/1, (1,0) => 6/1, (1,1) => 6/1, (2,0) => 0, (2,1) => 50/1, (3,0) => 0, (2,2) => 0, (4,0) => 0, (3,1) => 90/1, (3,2) => 0, (5,0) => 0, (4,1) => 78/1, (4,2) => 0, (6,0) => 0, (5,1) => 34/1};
--sb represents the betti numbers as sums of Schur functors
sb1022 = new HashTable from {(6,1) => {({8,7,8,6},1/1)}, (5,2) => {}, (7,1) => {}, (6,2) => {}, (8,1) => {}, (0,0) => {({1,0,0,0},1/1)}, (1,0) => {({2,1,2,0},1/1)}, (1,1) => {({5,0,2,2},1/1)}, (2,0) => {}, (2,1) => {({6,1,4,2},1/1),({5,2,4,2},1/1),({5,2,3,3},1/1),({4,3,5,1},1/1),({4,3,4,2},1/1)}, (3,0) => {}, (2,2) => {}, (4,0) => {}, (3,1) => {({7,2,5,3},1/1),({6,3,6,2},1/1),({6,3,5,3},2/1),({6,3,4,4},1/1),({5,4,6,2},1/1),({5,4,5,3},2/1),({5,4,4,4},1/1)}, (3,2) => {}, (5,0) => {}, (4,1) => {({8,3,5,5},1/1),({7,4,7,3},1/1),({7,4,6,4},2/1),({7,4,5,5},1/1),({6,5,7,3},1/1),({6,5,6,4},2/1),({6,5,5,5},1/1)}, (4,2) => {}, (6,0) => {}, (5,1) => {({8,5,7,5},1/1),({8,5,6,6},1/1),({7,6,8,4},1/1),({7,6,7,5},1/1),({7,6,6,6},1/1)}};
end;