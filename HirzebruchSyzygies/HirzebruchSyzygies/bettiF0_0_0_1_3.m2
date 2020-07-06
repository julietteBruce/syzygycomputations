A := QQ[t_0,t_1,t_2,t_3];
--tb stands for Total Betti numbers
tb0013 = new HashTable from {(5,2) => 0, (6,1) => 0, (7,1) => 0, (0,0) => 1/1, (1,0) => 0, (1,1) => 15/1, (2,0) => 0, (2,1) => 40/1, (3,0) => 0, (2,2) => 0, (3,1) => 45/1, (4,0) => 0, (4,1) => 24/1, (5,0) => 0, (3,2) => 0, (5,1) => 5/1, (4,2) => 0};
--sb represents the betti numbers as sums of Schur functors
sb0013 = new HashTable from {(5,2) => {}, (6,1) => {}, (7,1) => {}, (0,0) => {({0,0,0,0},1/1)}, (1,0) => {}, (1,1) => {({2,0,4,2},1/1),({1,1,5,1},1/1),({1,1,3,3},1/1)}, (2,0) => {}, (2,1) => {({3,0,5,4},1/1),({2,1,7,2},1/1),({2,1,6,3},2/1),({2,1,5,4},1/1)}, (3,0) => {}, (2,2) => {}, (3,1) => {({3,1,8,4},1/1),({3,1,7,5},1/1),({3,1,6,6},1/1),({2,2,9,3},1/1),({2,2,8,4},1/1),({2,2,7,5},2/1)}, (4,0) => {}, (4,1) => {({3,2,10,5},1/1),({3,2,9,6},1/1),({3,2,8,7},1/1)}, (5,0) => {}, (3,2) => {}, (5,1) => {({3,3,11,7},1/1)}, (4,2) => {}};
--dw stands for dominant weights
dw0013 = new HashTable from {(5,2) => {}, (6,1) => {}, (7,1) => {}, (0,0) => {({0,0,0,0},1/1)}, (1,0) => {}, (1,1) => {({2,0,4,2},1/1)}, (2,0) => {}, (2,1) => {({3,0,5,4},1/1)}, (3,0) => {}, (2,2) => {}, (3,1) => {({3,1,8,4},1/1)}, (4,0) => {}, (3,2) => {}, (5,0) => {}, (4,1) => {({3,2,10,5},1/1)}, (4,2) => {}, (5,1) => {({3,3,11,7},1/1)}};
end;