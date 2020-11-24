A := QQ[t_0,t_1,t_2,t_3];
--tb stands for Total Betti numbers
tb0012 = new HashTable from {(0,0) => 1/1, (1,0) => 0, (1,1) => 6/1, (2,0) => 0, (2,1) => 8/1, (3,0) => 0, (3,1) => 3/1, (2,2) => 0, (3,2) => 0, (4,1) => 0, (5,1) => 0};
--sb represents the betti numbers as sums of Schur functors
sb0012 = new HashTable from {(0,0) => {({0,0,0,0},1/1)}, (1,0) => {}, (1,1) => {({2,0,2,2},1/1),({1,1,3,1},1/1)}, (2,0) => {}, (2,1) => {({2,1,4,2},1/1),({2,1,3,3},1/1)}, (3,0) => {}, (3,1) => {({2,2,5,3},1/1)}, (2,2) => {}, (3,2) => {}, (4,1) => {}, (5,1) => {}};
--dw stands for dominant weights
dw0012 = new HashTable from {(0,0) => {({0,0,0,0},1/1)}, (1,0) => {}, (1,1) => {({2,0,2,2},1/1),({1,1,3,1},1/1)}, (2,0) => {}, (2,1) => {({2,1,4,2},1/1)}, (3,0) => {}, (2,2) => {}, (3,1) => {({2,2,5,3},1/1)}, (4,1) => {}, (3,2) => {}, (5,1) => {}};
--dmw stands for dominant monomial weights
dmw0012 = new HashTable from {(0,0) => {}, (1,0) => {}, (1,1) => {{2,0,2,2}}, (2,0) => {}, (2,1) => {{3,0,2,4},{2,1,3,3}}, (3,0) => {}, (2,2) => {}, (3,1) => {{3,1,3,5}}, (4,1) => {}, (3,2) => {}, (5,1) => {}};
end;