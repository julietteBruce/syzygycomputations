A := QQ[t_0,t_1,t_2,t_3];
--tb stands for Total Betti numbers
tb1135 = new HashTable from {(5,2) => 0, (7,0) => 0, (7,1) => 1697688/1, (9,0) => 0, (7,2) => 0, (9,1) => 4739056/1, (11,0) => 0, (9,2) => 0, (11,1) => 5878600/1, (13,0) => 0, (11,2) => 0, (13,1) => 3333360/1, (15,0) => 0, (13,2) => 0, (15,1) => 817836/1, (17,0) => 0, (15,2) => 0, (17,1) => 71288/1, (19,0) => 0, (17,2) => 0, (19,1) => 756/1, (21,0) => 0, (19,2) => 190/1, (21,1) => 0, (21,2) => 3/1, (0,0) => 4/1, (2,0) => 360/1, (2,1) => 127/1, (2,2) => 0, (4,0) => 0, (4,1) => 52269/1, (4,2) => 0, (6,0) => 0, (6,1) => 721905/1, (6,2) => 0, (8,0) => 0, (8,1) => 3155710/1, (10,0) => 0, (8,2) => 0, (10,1) => 5819814/1, (12,0) => 0, (10,2) => 0, (12,1) => 4888282/1, (14,0) => 0, (12,2) => 0, (14,1) => 1846914/1, (16,0) => 0, (14,2) => 0, (16,1) => 281295/1, (18,0) => 0, (16,2) => 0, (18,1) => 11795/1, (20,0) => 0, (18,2) => 0, (20,1) => 39/1, (20,2) => 40/1, (22,2) => 0, (1,0) => 61/1, (1,1) => 0, (3,0) => 680/1, (3,1) => 6020/1, (3,2) => 0, (5,0) => 0, (5,1) => 233016/1};
--sb represents the betti numbers as sums of Schur functors
sb1135 = new HashTable from {(5,2) => {}, (7,0) => {}, (7,1) => {({22,3,26,15},1/1),({22,3,25,16},1/1),({22,3,24,17},2/1),({22,3,23,18},2/1),({22,3,22,19},2/1),({22,3,21,20},1/1),({21,4,28,13},1/1),({21,4,27,14},3/1),({21,4,26,15},6/1),({21,4,25,16},9/1),({21,4,24,17},12/1),({21,4,23,18},13/1),({21,4,22,19},11/1),({21,4,21,20},6/1),({20,5,30,11},1/1),({20,5,29,12},3/1),({20,5,28,13},9/1),({20,5,27,14},17/1),({20,5,26,15},30/1),({20,5,25,16},40/1),({20,5,24,17},49/1),({20,5,23,18},49/1),({20,5,22,19},41/1),({20,5,21,20},23/1),({19,6,31,10},1/1),({19,6,30,11},6/1),({19,6,29,12},17/1),({19,6,28,13},35/1),({19,6,27,14},61/1),({19,6,26,15},92/1),({19,6,25,16},121/1),({19,6,24,17},138/1),({19,6,23,18},135/1),({19,6,22,19},109/1),({19,6,21,20},61/1),({18,7,32,9},3/1),({18,7,31,10},9/1),({18,7,30,11},26/1),({18,7,29,12},54/1),({18,7,28,13},101/1),({18,7,27,14},158/1),({18,7,26,15},224/1),({18,7,25,16},277/1),({18,7,24,17},307/1),({18,7,23,18},292/1),({18,7,22,19},232/1),({18,7,21,20},128/1),({17,8,33,8},2/1),({17,8,32,9},9/1),({17,8,31,10},27/1),({17,8,30,11},62/1),({17,8,29,12},121/1),({17,8,28,13},206/1),({17,8,27,14},311/1),({17,8,26,15},418/1),({17,8,25,16},505/1),({17,8,24,17},543/1),({17,8,23,18},511/1),({17,8,22,19},399/1),({17,8,21,20},219/1),({16,9,34,7},2/1),({16,9,33,8},7/1),({16,9,32,9},23/1),({16,9,31,10},55/1),({16,9,30,11},116/1),({16,9,29,12},208/1),({16,9,28,13},337/1),({16,9,27,14},485/1),({16,9,26,15},636/1),({16,9,25,16},747/1),({16,9,24,17},791/1),({16,9,23,18},732/1),({16,9,22,19},568/1),({16,9,21,20},310/1),({15,10,34,7},3/1),({15,10,33,8},12/1),({15,10,32,9},34/1),({15,10,31,10},79/1),({15,10,30,11},156/1),({15,10,29,12},272/1),({15,10,28,13},424/1),({15,10,27,14},598/1),({15,10,26,15},766/1),({15,10,25,16},889/1),({15,10,24,17},927/1),({15,10,23,18},852/1),({15,10,22,19},656/1),({15,10,21,20},358/1),({14,11,35,6},1/1),({14,11,34,7},5/1),({14,11,33,8},15/1),({14,11,32,9},39/1),({14,11,31,10},83/1),({14,11,30,11},159/1),({14,11,29,12},267/1),({14,11,28,13},409/1),({14,11,27,14},565/1),({14,11,26,15},715/1),({14,11,25,16},819/1),({14,11,24,17},849/1),({14,11,23,18},774/1),({14,11,22,19},594/1),({14,11,21,20},322/1),({13,12,35,6},1/1),({13,12,34,7},3/1),({13,12,33,8},10/1),({13,12,32,9},25/1),({13,12,31,10},54/1),({13,12,30,11},100/1),({13,12,29,12},167/1),({13,12,28,13},252/1),({13,12,27,14},347/1),({13,12,26,15},434/1),({13,12,25,16},495/1),({13,12,24,17},510/1),({13,12,23,18},465/1),({13,12,22,19},355/1),({13,12,21,20},192/1)}, (9,0) => {}, (7,2) => {}, (9,1) => {({26,5,29,22},1/1),({26,5,28,23},1/1),({26,5,27,24},1/1),({25,6,33,18},1/1),({25,6,32,19},2/1),({25,6,31,20},4/1),({25,6,30,21},6/1),({25,6,29,22},8/1),({25,6,28,23},9/1),({25,6,27,24},8/1),({25,6,26,25},4/1),({24,7,35,16},1/1),({24,7,34,17},2/1),({24,7,33,18},8/1),({24,7,32,19},15/1),({24,7,31,20},26/1),({24,7,30,21},35/1),({24,7,29,22},43/1),({24,7,28,23},43/1),({24,7,27,24},36/1),({24,7,26,25},20/1),({23,8,36,15},2/1),({23,8,35,16},8/1),({23,8,34,17},20/1),({23,8,33,18},40/1),({23,8,32,19},69/1),({23,8,31,20},102/1),({23,8,30,21},132/1),({23,8,29,22},150/1),({23,8,28,23},146/1),({23,8,27,24},117/1),({23,8,26,25},66/1),({22,9,37,14},4/1),({22,9,36,15},13/1),({22,9,35,16},35/1),({22,9,34,17},73/1),({22,9,33,18},134/1),({22,9,32,19},207/1),({22,9,31,20},291/1),({22,9,30,21},358/1),({22,9,29,22},393/1),({22,9,28,23},373/1),({22,9,27,24},296/1),({22,9,26,25},162/1),({21,10,39,12},1/1),({21,10,38,13},5/1),({21,10,37,14},18/1),({21,10,36,15},48/1),({21,10,35,16},105/1),({21,10,34,17},198/1),({21,10,33,18},328/1),({21,10,32,19},485/1),({21,10,31,20},644/1),({21,10,30,21},770/1),({21,10,29,22},822/1),({21,10,28,23},768/1),({21,10,27,24},598/1),({21,10,26,25},328/1),({20,11,39,12},4/1),({20,11,38,13},16/1),({20,11,37,14},48/1),({20,11,36,15},111/1),({20,11,35,16},226/1),({20,11,34,17},395/1),({20,11,33,18},626/1),({20,11,32,19},888/1),({20,11,31,20},1147/1),({20,11,30,21},1335/1),({20,11,29,22},1403/1),({20,11,28,23},1291/1),({20,11,27,24},997/1),({20,11,26,25},544/1),({19,12,40,11},3/1),({19,12,39,12},12/1),({19,12,38,13},37/1),({19,12,37,14},93/1),({19,12,36,15},200/1),({19,12,35,16},374/1),({19,12,34,17},630/1),({19,12,33,18},957/1),({19,12,32,19},1322/1),({19,12,31,20},1667/1),({19,12,30,21},1912/1),({19,12,29,22},1974/1),({19,12,28,23},1804/1),({19,12,27,24},1383/1),({19,12,26,25},750/1),({18,13,41,10},1/1),({18,13,40,11},5/1),({18,13,39,12},20/1),({18,13,38,13},55/1),({18,13,37,14},131/1),({18,13,36,15},265/1),({18,13,35,16},482/1),({18,13,34,17},783/1),({18,13,33,18},1165/1),({18,13,32,19},1577/1),({18,13,31,20},1960/1),({18,13,30,21},2217/1),({18,13,29,22},2273/1),({18,13,28,23},2059/1),({18,13,27,24},1571/1),({18,13,26,25},850/1),({17,14,41,10},2/1),({17,14,40,11},7/1),({17,14,39,12},23/1),({17,14,38,13},61/1),({17,14,37,14},135/1),({17,14,36,15},265/1),({17,14,35,16},468/1),({17,14,34,17},746/1),({17,14,33,18},1088/1),({17,14,32,19},1459/1),({17,14,31,20},1790/1),({17,14,30,21},2010/1),({17,14,29,22},2047/1),({17,14,28,23},1846/1),({17,14,27,24},1401/1),({17,14,26,25},759/1),({16,15,41,10},1/1),({16,15,40,11},5/1),({16,15,39,12},16/1),({16,15,38,13},39/1),({16,15,37,14},87/1),({16,15,36,15},167/1),({16,15,35,16},290/1),({16,15,34,17},459/1),({16,15,33,18},665/1),({16,15,32,19},881/1),({16,15,31,20},1078/1),({16,15,30,21},1205/1),({16,15,29,22},1221/1),({16,15,28,23},1099/1),({16,15,27,24},835/1),({16,15,26,25},448/1)}, (11,0) => {}, (9,2) => {}, (11,1) => {({29,8,34,27},1/1),({29,8,33,28},1/1),({29,8,32,29},1/1),({29,8,31,30},1/1),({28,9,38,23},1/1),({28,9,37,24},2/1),({28,9,36,25},5/1),({28,9,35,26},7/1),({28,9,34,27},10/1),({28,9,33,28},10/1),({28,9,32,29},9/1),({28,9,31,30},5/1),({27,10,40,21},1/1),({27,10,39,22},4/1),({27,10,38,23},10/1),({27,10,37,24},19/1),({27,10,36,25},31/1),({27,10,35,26},43/1),({27,10,34,27},51/1),({27,10,33,28},51/1),({27,10,32,29},42/1),({27,10,31,30},24/1),({26,11,42,19},1/1),({26,11,41,20},3/1),({26,11,40,21},12/1),({26,11,39,22},26/1),({26,11,38,23},52/1),({26,11,37,24},86/1),({26,11,36,25},126/1),({26,11,35,26},159/1),({26,11,34,27},181/1),({26,11,33,28},174/1),({26,11,32,29},139/1),({26,11,31,30},78/1),({25,12,43,18},1/1),({25,12,42,19},6/1),({25,12,41,20},20/1),({25,12,40,21},48/1),({25,12,39,22},97/1),({25,12,38,23},170/1),({25,12,37,24},261/1),({25,12,36,25},357/1),({25,12,35,26},436/1),({25,12,34,27},473/1),({25,12,33,28},448/1),({25,12,32,29},352/1),({25,12,31,30},193/1),({24,13,44,17},2/1),({24,13,43,18},8/1),({24,13,42,19},27/1),({24,13,41,20},66/1),({24,13,40,21},141/1),({24,13,39,22},255/1),({24,13,38,23},416/1),({24,13,37,24},602/1),({24,13,36,25},793/1),({24,13,35,26},935/1),({24,13,34,27},993/1),({24,13,33,28},921/1),({24,13,32,29},716/1),({24,13,31,30},391/1),({23,14,45,16},1/1),({23,14,44,17},7/1),({23,14,43,18},26/1),({23,14,42,19},68/1),({23,14,41,20},152/1),({23,14,40,21},296/1),({23,14,39,22},509/1),({23,14,38,23},786/1),({23,14,37,24},1104/1),({23,14,36,25},1406/1),({23,14,35,26},1626/1),({23,14,34,27},1694/1),({23,14,33,28},1554/1),({23,14,32,29},1194/1),({23,14,31,30},652/1),({22,15,46,15},1/1),({22,15,45,16},5/1),({22,15,44,17},19/1),({22,15,43,18},53/1),({22,15,42,19},129/1),({22,15,41,20},265/1),({22,15,40,21},487/1),({22,15,39,22},800/1),({22,15,38,23},1200/1),({22,15,37,24},1634/1),({22,15,36,25},2043/1),({22,15,35,26},2322/1),({22,15,34,27},2388/1),({22,15,33,28},2169/1),({22,15,32,29},1659/1),({22,15,31,30},898/1),({21,16,46,15},2/1),({21,16,45,16},9/1),({21,16,44,17},30/1),({21,16,43,18},79/1),({21,16,42,19},178/1),({21,16,41,20},351/1),({21,16,40,21},620/1),({21,16,39,22},993/1),({21,16,38,23},1453/1),({21,16,37,24},1949/1),({21,16,36,25},2397/1),({21,16,35,26},2696/1),({21,16,34,27},2746/1),({21,16,33,28},2480/1),({21,16,32,29},1885/1),({21,16,31,30},1019/1),({20,17,46,15},3/1),({20,17,45,16},11/1),({20,17,44,17},34/1),({20,17,43,18},84/1),({20,17,42,19},182/1),({20,17,41,20},346/1),({20,17,40,21},600/1),({20,17,39,22},940/1),({20,17,38,23},1356/1),({20,17,37,24},1796/1),({20,17,36,25},2189/1),({20,17,35,26},2440/1),({20,17,34,27},2474/1),({20,17,33,28},2221/1),({20,17,32,29},1683/1),({20,17,31,30},910/1),({19,18,47,14},1/1),({19,18,46,15},2/1),({19,18,45,16},8/1),({19,18,44,17},23/1),({19,18,43,18},55/1),({19,18,42,19},116/1),({19,18,41,20},219/1),({19,18,40,21},371/1),({19,18,39,22},579/1),({19,18,38,23},827/1),({19,18,37,24},1086/1),({19,18,36,25},1317/1),({19,18,35,26},1465/1),({19,18,34,27},1475/1),({19,18,33,28},1323/1),({19,18,32,29},1002/1),({19,18,31,30},539/1)}, (13,0) => {}, (11,2) => {}, (13,1) => {({31,12,41,30},1/1),({31,12,40,31},2/1),({31,12,39,32},3/1),({31,12,38,33},3/1),({31,12,37,34},3/1),({31,12,36,35},2/1),({30,13,43,28},3/1),({30,13,42,29},6/1),({30,13,41,30},11/1),({30,13,40,31},16/1),({30,13,39,32},21/1),({30,13,38,33},21/1),({30,13,37,34},18/1),({30,13,36,35},10/1),({29,14,46,25},1/1),({29,14,45,26},3/1),({29,14,44,27},9/1),({29,14,43,28},20/1),({29,14,42,29},36/1),({29,14,41,30},55/1),({29,14,40,31},73/1),({29,14,39,32},85/1),({29,14,38,33},84/1),({29,14,37,34},68/1),({29,14,36,35},38/1),({28,15,47,24},2/1),({28,15,46,25},7/1),({28,15,45,26},20/1),({28,15,44,27},42/1),({28,15,43,28},79/1),({28,15,42,29},125/1),({28,15,41,30},178/1),({28,15,40,31},221/1),({28,15,39,32},246/1),({28,15,38,33},235/1),({28,15,37,34},187/1),({28,15,36,35},103/1),({27,16,49,22},1/1),({27,16,48,23},3/1),({27,16,47,24},11/1),({27,16,46,25},30/1),({27,16,45,26},67/1),({27,16,44,27},127/1),({27,16,43,28},213/1),({27,16,42,29},318/1),({27,16,41,30},426/1),({27,16,40,31},512/1),({27,16,39,32},549/1),({27,16,38,33},515/1),({27,16,37,34},402/1),({27,16,36,35},221/1),({26,17,49,22},3/1),({26,17,48,23},11/1),({26,17,47,24},33/1),({26,17,46,25},76/1),({26,17,45,26},155/1),({26,17,44,27},272/1),({26,17,43,28},433/1),({26,17,42,29},616/1),({26,17,41,30},799/1),({26,17,40,31},933/1),({26,17,39,32},982/1),({26,17,38,33},906/1),({26,17,37,34},701/1),({26,17,36,35},382/1),({25,18,50,21},2/1),({25,18,49,22},8/1),({25,18,48,23},26/1),({25,18,47,24},66/1),({25,18,46,25},142/1),({25,18,45,26},268/1),({25,18,44,27},452/1),({25,18,43,28},688/1),({25,18,42,29},953/1),({25,18,41,30},1203/1),({25,18,40,31},1381/1),({25,18,39,32},1429/1),({25,18,38,33},1306/1),({25,18,37,34},1001/1),({25,18,36,35},544/1),({24,19,51,20},1/1),({24,19,50,21},4/1),({24,19,49,22},15/1),({24,19,48,23},41/1),({24,19,47,24},97/1),({24,19,46,25},196/1),({24,19,45,26},356/1),({24,19,44,27},579/1),({24,19,43,28},862/1),({24,19,42,29},1167/1),({24,19,41,30},1451/1),({24,19,40,31},1642/1),({24,19,39,32},1684/1),({24,19,38,33},1525/1),({24,19,37,34},1164/1),({24,19,36,35},630/1),({23,20,51,20},1/1),({23,20,50,21},5/1),({23,20,49,22},17/1),({23,20,48,23},45/1),({23,20,47,24},101/1),({23,20,46,25},199/1),({23,20,45,26},351/1),({23,20,44,27},561/1),({23,20,43,28},819/1),({23,20,42,29},1096/1),({23,20,41,30},1346/1),({23,20,40,31},1512/1),({23,20,39,32},1538/1),({23,20,38,33},1387/1),({23,20,37,34},1054/1),({23,20,36,35},570/1),({22,21,51,20},1/1),({22,21,50,21},4/1),({22,21,49,22},12/1),({22,21,48,23},30/1),({22,21,47,24},66/1),({22,21,46,25},127/1),({22,21,45,26},221/1),({22,21,44,27},349/1),({22,21,43,28},505/1),({22,21,42,29},670/1),({22,21,41,30},819/1),({22,21,40,31},914/1),({22,21,39,32},927/1),({22,21,38,33},833/1),({22,21,37,34},632/1),({22,21,36,35},341/1)}, (15,0) => {}, (13,2) => {}, (15,1) => {({33,16,45,36},1/1),({33,16,44,37},2/1),({33,16,43,38},2/1),({33,16,42,39},2/1),({33,16,41,40},1/1),({32,17,48,33},1/1),({32,17,47,34},3/1),({32,17,46,35},6/1),({32,17,45,36},9/1),({32,17,44,37},13/1),({32,17,43,38},13/1),({32,17,42,39},11/1),({32,17,41,40},7/1),({31,18,50,31},1/1),({31,18,49,32},4/1),({31,18,48,33},10/1),({31,18,47,34},18/1),({31,18,46,35},29/1),({31,18,45,36},40/1),({31,18,44,37},47/1),({31,18,43,38},47/1),({31,18,42,39},39/1),({31,18,41,40},22/1),({30,19,52,29},1/1),({30,19,51,30},3/1),({30,19,50,31},9/1),({30,19,49,32},19/1),({30,19,48,33},38/1),({30,19,47,34},61/1),({30,19,46,35},89/1),({30,19,45,36},112/1),({30,19,44,37},126/1),({30,19,43,38},121/1),({30,19,42,39},97/1),({30,19,41,40},54/1),({29,20,53,28},1/1),({29,20,52,29},4/1),({29,20,51,30},12/1),({29,20,50,31},28/1),({29,20,49,32},54/1),({29,20,48,33},93/1),({29,20,47,34},142/1),({29,20,46,35},192/1),({29,20,45,36},233/1),({29,20,44,37},253/1),({29,20,43,38},238/1),({29,20,42,39},186/1),({29,20,41,40},103/1),({28,21,54,27},1/1),({28,21,53,28},3/1),({28,21,52,29},11/1),({28,21,51,30},27/1),({28,21,50,31},57/1),({28,21,49,32},104/1),({28,21,48,33},169/1),({28,21,47,34},243/1),({28,21,46,35},320/1),({28,21,45,36},377/1),({28,21,44,37},399/1),({28,21,43,38},370/1),({28,21,42,39},288/1),({28,21,41,40},156/1),({27,22,54,27},2/1),({27,22,53,28},7/1),({27,22,52,29},19/1),({27,22,51,30},43/1),({27,22,50,31},85/1),({27,22,49,32},147/1),({27,22,48,33},229/1),({27,22,47,34},322/1),({27,22,46,35},412/1),({27,22,45,36},477/1),({27,22,44,37},497/1),({27,22,43,38},456/1),({27,22,42,39},351/1),({27,22,41,40},191/1),({26,23,54,27},2/1),({26,23,53,28},8/1),({26,23,52,29},21/1),({26,23,51,30},46/1),({26,23,50,31},90/1),({26,23,49,32},150/1),({26,23,48,33},229/1),({26,23,47,34},317/1),({26,23,46,35},399/1),({26,23,45,36},455/1),({26,23,44,37},472/1),({26,23,43,38},429/1),({26,23,42,39},328/1),({26,23,41,40},179/1),({25,24,55,26},1/1),({25,24,54,27},2/1),({25,24,53,28},6/1),({25,24,52,29},15/1),({25,24,51,30},32/1),({25,24,50,31},58/1),({25,24,49,32},97/1),({25,24,48,33},146/1),({25,24,47,34},199/1),({25,24,46,35},249/1),({25,24,45,36},283/1),({25,24,44,37},289/1),({25,24,43,38},264/1),({25,24,42,39},202/1),({25,24,41,40},108/1)}, (17,0) => {}, (15,2) => {}, (17,1) => {({35,20,46,45},1/1),({34,21,50,41},1/1),({34,21,49,42},2/1),({34,21,48,43},2/1),({34,21,47,44},2/1),({34,21,46,45},1/1),({33,22,53,38},1/1),({33,22,52,39},2/1),({33,22,51,40},4/1),({33,22,50,41},6/1),({33,22,49,42},8/1),({33,22,48,43},8/1),({33,22,47,44},7/1),({33,22,46,45},4/1),({32,23,55,36},1/1),({32,23,54,37},2/1),({32,23,53,38},5/1),({32,23,52,39},9/1),({32,23,51,40},14/1),({32,23,50,41},18/1),({32,23,49,42},22/1),({32,23,48,43},21/1),({32,23,47,44},17/1),({32,23,46,45},10/1),({31,24,56,35},1/1),({31,24,55,36},3/1),({31,24,54,37},7/1),({31,24,53,38},13/1),({31,24,52,39},21/1),({31,24,51,40},30/1),({31,24,50,41},37/1),({31,24,49,42},41/1),({31,24,48,43},39/1),({31,24,47,44},31/1),({31,24,46,45},17/1),({30,25,57,34},1/1),({30,25,56,35},3/1),({30,25,55,36},7/1),({30,25,54,37},13/1),({30,25,53,38},23/1),({30,25,52,39},34/1),({30,25,51,40},46/1),({30,25,50,41},54/1),({30,25,49,42},58/1),({30,25,48,43},54/1),({30,25,47,44},42/1),({30,25,46,45},23/1),({29,26,57,34},1/1),({29,26,56,35},3/1),({29,26,55,36},8/1),({29,26,54,37},15/1),({29,26,53,38},25/1),({29,26,52,39},37/1),({29,26,51,40},48/1),({29,26,50,41},56/1),({29,26,49,42},60/1),({29,26,48,43},55/1),({29,26,47,44},42/1),({29,26,46,45},23/1),({28,27,57,34},1/1),({28,27,56,35},3/1),({28,27,55,36},6/1),({28,27,54,37},11/1),({28,27,53,38},18/1),({28,27,52,39},25/1),({28,27,51,40},32/1),({28,27,50,41},37/1),({28,27,49,42},38/1),({28,27,48,43},35/1),({28,27,47,44},27/1),({28,27,46,45},14/1)}, (19,0) => {}, (17,2) => {}, (19,1) => {({33,28,58,43},1/1),({33,28,57,44},1/1),({33,28,56,45},1/1),({33,28,55,46},1/1),({33,28,54,47},1/1),({32,29,58,43},1/1),({32,29,57,44},1/1),({32,29,56,45},1/1),({32,29,55,46},1/1),({32,29,54,47},1/1),({31,30,59,42},1/1),({31,30,58,43},1/1),({31,30,57,44},1/1),({31,30,56,45},1/1),({31,30,55,46},1/1),({31,30,54,47},1/1)}, (21,0) => {}, (19,2) => {({35,29,56,50},1/1),({35,29,54,52},1/1),({34,30,57,49},1/1),({34,30,55,51},1/1),({34,30,53,53},1/1),({33,31,56,50},1/1),({33,31,54,52},1/1),({32,32,57,49},1/1),({32,32,55,51},1/1),({32,32,53,53},1/1)}, (21,1) => {}, (21,2) => {({35,35,59,57},1/1)}, (0,0) => {({1,0,1,0},1/1)}, (2,0) => {({6,1,10,1},1/1),({6,1,9,2},1/1),({6,1,8,3},1/1),({6,1,7,4},1/1),({6,1,6,5},1/1),({5,2,10,1},1/1),({5,2,9,2},1/1),({5,2,8,3},1/1),({5,2,7,4},1/1),({5,2,6,5},1/1),({4,3,10,1},1/1),({4,3,9,2},1/1),({4,3,8,3},1/1),({4,3,7,4},1/1),({4,3,6,5},1/1)}, (2,1) => {({10,0,11,5},1/1),({10,0,9,7},1/1),({5,5,16,0},1/1)}, (2,2) => {}, (4,0) => {}, (4,1) => {({16,0,14,12},1/1),({15,1,18,8},1/1),({15,1,17,9},2/1),({15,1,16,10},2/1),({15,1,15,11},2/1),({15,1,14,12},2/1),({15,1,13,13},1/1),({14,2,20,6},1/1),({14,2,19,7},1/1),({14,2,18,8},4/1),({14,2,17,9},5/1),({14,2,16,10},8/1),({14,2,15,11},6/1),({14,2,14,12},6/1),({14,2,13,13},1/1),({13,3,21,5},1/1),({13,3,20,6},3/1),({13,3,19,7},6/1),({13,3,18,8},10/1),({13,3,17,9},14/1),({13,3,16,10},16/1),({13,3,15,11},15/1),({13,3,14,12},11/1),({13,3,13,13},4/1),({12,4,22,4},1/1),({12,4,21,5},3/1),({12,4,20,6},8/1),({12,4,19,7},13/1),({12,4,18,8},21/1),({12,4,17,9},25/1),({12,4,16,10},30/1),({12,4,15,11},25/1),({12,4,14,12},20/1),({12,4,13,13},5/1),({11,5,23,3},1/1),({11,5,22,4},3/1),({11,5,21,5},7/1),({11,5,20,6},13/1),({11,5,19,7},22/1),({11,5,18,8},30/1),({11,5,17,9},37/1),({11,5,16,10},38/1),({11,5,15,11},35/1),({11,5,14,12},24/1),({11,5,13,13},9/1),({10,6,24,2},1/1),({10,6,23,3},2/1),({10,6,22,4},6/1),({10,6,21,5},10/1),({10,6,20,6},19/1),({10,6,19,7},26/1),({10,6,18,8},37/1),({10,6,17,9},40/1),({10,6,16,10},45/1),({10,6,15,11},36/1),({10,6,14,12},28/1),({10,6,13,13},7/1),({9,7,23,3},2/1),({9,7,22,4},4/1),({9,7,21,5},9/1),({9,7,20,6},14/1),({9,7,19,7},22/1),({9,7,18,8},27/1),({9,7,17,9},33/1),({9,7,16,10},32/1),({9,7,15,11},30/1),({9,7,14,12},18/1),({9,7,13,13},8/1),({8,8,24,2},1/1),({8,8,23,3},1/1),({8,8,22,4},3/1),({8,8,21,5},4/1),({8,8,20,6},8/1),({8,8,19,7},9/1),({8,8,18,8},13/1),({8,8,17,9},12/1),({8,8,16,10},15/1),({8,8,15,11},10/1),({8,8,14,12},9/1),({8,8,13,13},1/1)}, (4,2) => {}, (6,0) => {}, (6,1) => {({20,2,23,13},1/1),({20,2,22,14},1/1),({20,2,21,15},2/1),({20,2,20,16},1/1),({20,2,19,17},2/1),({19,3,25,11},1/1),({19,3,24,12},3/1),({19,3,23,13},5/1),({19,3,22,14},8/1),({19,3,21,15},10/1),({19,3,20,16},9/1),({19,3,19,17},7/1),({19,3,18,18},3/1),({18,4,27,9},1/1),({18,4,26,10},3/1),({18,4,25,11},8/1),({18,4,24,12},14/1),({18,4,23,13},24/1),({18,4,22,14},30/1),({18,4,21,15},36/1),({18,4,20,16},31/1),({18,4,19,17},25/1),({18,4,18,18},7/1),({17,5,28,8},1/1),({17,5,27,9},5/1),({17,5,26,10},14/1),({17,5,25,11},28/1),({17,5,24,12},47/1),({17,5,23,13},68/1),({17,5,22,14},85/1),({17,5,21,15},91/1),({17,5,20,16},83/1),({17,5,19,17},58/1),({17,5,18,18},21/1),({16,6,29,7},2/1),({16,6,28,8},7/1),({16,6,27,9},20/1),({16,6,26,10},39/1),({16,6,25,11},73/1),({16,6,24,12},109/1),({16,6,23,13},151/1),({16,6,22,14},176/1),({16,6,21,15},189/1),({16,6,20,16},161/1),({16,6,19,17},117/1),({16,6,18,18},38/1),({15,7,30,6},2/1),({15,7,29,7},7/1),({15,7,28,8},20/1),({15,7,27,9},44/1),({15,7,26,10},84/1),({15,7,25,11},137/1),({15,7,24,12},201/1),({15,7,23,13},259/1),({15,7,22,14},301/1),({15,7,21,15},305/1),({15,7,20,16},268/1),({15,7,19,17},181/1),({15,7,18,18},66/1),({14,8,31,5},1/1),({14,8,30,6},4/1),({14,8,29,7},15/1),({14,8,28,8},34/1),({14,8,27,9},73/1),({14,8,26,10},127/1),({14,8,25,11},203/1),({14,8,24,12},280/1),({14,8,23,13},361/1),({14,8,22,14},401/1),({14,8,21,15},409/1),({14,8,20,16},346/1),({14,8,19,17},242/1),({14,8,18,18},79/1),({13,9,31,5},2/1),({13,9,30,6},8/1),({13,9,29,7},20/1),({13,9,28,8},47/1),({13,9,27,9},89/1),({13,9,26,10},152/1),({13,9,25,11},229/1),({13,9,24,12},317/1),({13,9,23,13},387/1),({13,9,22,14},436/1),({13,9,21,15},429/1),({13,9,20,16},369/1),({13,9,19,17},246/1),({13,9,18,18},91/1),({12,10,31,5},2/1),({12,10,30,6},6/1),({12,10,29,7},18/1),({12,10,28,8},38/1),({12,10,27,9},75/1),({12,10,26,10},121/1),({12,10,25,11},186/1),({12,10,24,12},246/1),({12,10,23,13},308/1),({12,10,22,14},333/1),({12,10,21,15},336/1),({12,10,20,16},278/1),({12,10,19,17},195/1),({12,10,18,18},62/1),({11,11,32,4},1/1),({11,11,31,5},1/1),({11,11,30,6},4/1),({11,11,29,7},8/1),({11,11,28,8},18/1),({11,11,27,9},29/1),({11,11,26,10},52/1),({11,11,25,11},71/1),({11,11,24,12},100/1),({11,11,23,13},116/1),({11,11,22,14},133/1),({11,11,21,15},123/1),({11,11,20,16},113/1),({11,11,19,17},68/1),({11,11,18,18},29/1)}, (6,2) => {}, (8,0) => {}, (8,1) => {({24,4,28,18},1/1),({24,4,27,19},1/1),({24,4,26,20},2/1),({24,4,25,21},1/1),({24,4,24,22},2/1),({23,5,31,15},1/1),({23,5,30,16},2/1),({23,5,29,17},4/1),({23,5,28,18},8/1),({23,5,27,19},11/1),({23,5,26,20},12/1),({23,5,25,21},12/1),({23,5,24,22},9/1),({23,5,23,23},3/1),({22,6,32,14},2/1),({22,6,31,15},6/1),({22,6,30,16},15/1),({22,6,29,17},25/1),({22,6,28,18},39/1),({22,6,27,19},48/1),({22,6,26,20},56/1),({22,6,25,21},48/1),({22,6,24,22},37/1),({22,6,23,23},11/1),({21,7,34,12},1/1),({21,7,33,13},5/1),({21,7,32,14},14/1),({21,7,31,15},32/1),({21,7,30,16},59/1),({21,7,29,17},94/1),({21,7,28,18},129/1),({21,7,27,19},157/1),({21,7,26,20},165/1),({21,7,25,21},148/1),({21,7,24,22},102/1),({21,7,23,23},37/1),({20,8,35,11},1/1),({20,8,34,12},8/1),({20,8,33,13},22/1),({20,8,32,14},53/1),({20,8,31,15},100/1),({20,8,30,16},171/1),({20,8,29,17},247/1),({20,8,28,18},330/1),({20,8,27,19},377/1),({20,8,26,20},393/1),({20,8,25,21},338/1),({20,8,24,22},238/1),({20,8,23,23},78/1),({19,9,36,10},2/1),({19,9,35,11},9/1),({19,9,34,12},26/1),({19,9,33,13},65/1),({19,9,32,14},132/1),({19,9,31,15},234/1),({19,9,30,16},366/1),({19,9,29,17},517/1),({19,9,28,18},650/1),({19,9,27,19},739/1),({19,9,26,20},739/1),({19,9,25,21},639/1),({19,9,24,22},432/1),({19,9,23,23},156/1),({18,10,37,9},1/1),({18,10,36,10},7/1),({18,10,35,11},23/1),({18,10,34,12},62/1),({18,10,33,13},132/1),({18,10,32,14},253/1),({18,10,31,15},418/1),({18,10,30,16},634/1),({18,10,29,17},853/1),({18,10,28,18},1058/1),({18,10,27,19},1164/1),({18,10,26,20},1161/1),({18,10,25,21},979/1),({18,10,24,22},671/1),({18,10,23,23},228/1),({17,11,37,9},4/1),({17,11,36,10},14/1),({17,11,35,11},42/1),({17,11,34,12},100/1),({17,11,33,13},205/1),({17,11,32,14},365/1),({17,11,31,15},594/1),({17,11,30,16},861/1),({17,11,29,17},1147/1),({17,11,28,18},1382/1),({17,11,27,19},1517/1),({17,11,26,20},1475/1),({17,11,25,21},1257/1),({17,11,24,22},836/1),({17,11,23,23},299/1),({16,12,38,8},2/1),({16,12,37,9},6/1),({16,12,36,10},22/1),({16,12,35,11},55/1),({16,12,34,12},125/1),({16,12,33,13},238/1),({16,12,32,14},420/1),({16,12,31,15},650/1),({16,12,30,16},939/1),({16,12,29,17},1215/1),({16,12,28,18},1459/1),({16,12,27,19},1567/1),({16,12,26,20},1536/1),({16,12,25,21},1275/1),({16,12,24,22},868/1),({16,12,23,23},294/1),({15,13,38,8},1/1),({15,13,37,9},6/1),({15,13,36,10},18/1),({15,13,35,11},48/1),({15,13,34,12},101/1),({15,13,33,13},197/1),({15,13,32,14},332/1),({15,13,31,15},520/1),({15,13,30,16},727/1),({15,13,29,17},949/1),({15,13,28,18},1114/1),({15,13,27,19},1210/1),({15,13,26,20},1156/1),({15,13,25,21},981/1),({15,13,24,22},644/1),({15,13,23,23},236/1),({14,14,38,8},1/1),({14,14,37,9},3/1),({14,14,36,10},9/1),({14,14,35,11},19/1),({14,14,34,12},44/1),({14,14,33,13},76/1),({14,14,32,14},133/1),({14,14,31,15},198/1),({14,14,30,16},285/1),({14,14,29,17},355/1),({14,14,28,18},432/1),({14,14,27,19},449/1),({14,14,26,20},445/1),({14,14,25,21},360/1),({14,14,24,22},253/1),({14,14,23,23},76/1)}, (10,0) => {}, (8,2) => {}, (10,1) => {({28,6,29,27},1/1),({27,7,34,22},1/1),({27,7,33,23},2/1),({27,7,32,24},3/1),({27,7,31,25},4/1),({27,7,30,26},4/1),({27,7,29,27},3/1),({27,7,28,28},1/1),({26,8,37,19},1/1),({26,8,36,20},2/1),({26,8,35,21},6/1),({26,8,34,22},10/1),({26,8,33,23},18/1),({26,8,32,24},22/1),({26,8,31,25},27/1),({26,8,30,26},23/1),({26,8,29,27},19/1),({26,8,28,28},5/1),({25,9,38,18},2/1),({25,9,37,19},7/1),({25,9,36,20},18/1),({25,9,35,21},34/1),({25,9,34,22},56/1),({25,9,33,23},79/1),({25,9,32,24},98/1),({25,9,31,25},104/1),({25,9,30,26},94/1),({25,9,29,27},65/1),({25,9,28,28},24/1),({24,10,40,16},1/1),({24,10,39,17},5/1),({24,10,38,18},15/1),({24,10,37,19},38/1),({24,10,36,20},73/1),({24,10,35,21},128/1),({24,10,34,22},188/1),({24,10,33,23},254/1),({24,10,32,24},293/1),({24,10,31,25},308/1),({24,10,30,26},265/1),({24,10,29,27},188/1),({24,10,28,28},62/1),({23,11,41,15},1/1),({23,11,40,16},7/1),({23,11,39,17},23/1),({23,11,38,18},58/1),({23,11,37,19},119/1),({23,11,36,20},215/1),({23,11,35,21},339/1),({23,11,34,22},482/1),({23,11,33,23},609/1),({23,11,32,24},695/1),({23,11,31,25},697/1),({23,11,30,26},605/1),({23,11,29,27},409/1),({23,11,28,28},147/1),({22,12,42,14},1/1),({22,12,41,15},8/1),({22,12,40,16},25/1),({22,12,39,17},68/1),({22,12,38,18},146/1),({22,12,37,19},281/1),({22,12,36,20},466/1),({22,12,35,21},708/1),({22,12,34,22},955/1),({22,12,33,23},1186/1),({22,12,32,24},1307/1),({22,12,31,25},1304/1),({22,12,30,26},1101/1),({22,12,29,27},755/1),({22,12,28,28},257/1),({21,13,43,13},1/1),({21,13,42,14},6/1),({21,13,41,15},21/1),({21,13,40,16},61/1),({21,13,39,17},140/1),({21,13,38,18},284/1),({21,13,37,19},503/1),({21,13,36,20},810/1),({21,13,35,21},1171/1),({21,13,34,22},1555/1),({21,13,33,23},1867/1),({21,13,32,24},2046/1),({21,13,31,25},1990/1),({21,13,30,26},1690/1),({21,13,29,27},1126/1),({21,13,28,28},403/1),({20,14,43,13},3/1),({20,14,42,14},12/1),({20,14,41,15},40/1),({20,14,40,16},99/1),({20,14,39,17},219/1),({20,14,38,18},414/1),({20,14,37,19},717/1),({20,14,36,20},1105/1),({20,14,35,21},1577/1),({20,14,34,22},2033/1),({20,14,33,23},2427/1),({20,14,32,24},2602/1),({20,14,31,25},2534/1),({20,14,30,26},2110/1),({20,14,29,27},1426/1),({20,14,28,28},487/1),({19,15,44,12},1/1),({19,15,43,13},5/1),({19,15,42,14},19/1),({19,15,41,15},52/1),({19,15,40,16},125/1),({19,15,39,17},256/1),({19,15,38,18},476/1),({19,15,37,19},788/1),({19,15,36,20},1204/1),({19,15,35,21},1670/1),({19,15,34,22},2146/1),({19,15,33,23},2507/1),({19,15,32,24},2693/1),({19,15,31,25},2575/1),({19,15,30,26},2165/1),({19,15,29,27},1428/1),({19,15,28,28},512/1),({18,16,44,12},1/1),({18,16,43,13},5/1),({18,16,42,14},16/1),({18,16,41,15},46/1),({18,16,40,16},103/1),({18,16,39,17},213/1),({18,16,38,18},379/1),({18,16,37,19},631/1),({18,16,36,20},937/1),({18,16,35,21},1305/1),({18,16,34,22},1642/1),({18,16,33,23},1931/1),({18,16,32,24},2036/1),({18,16,31,25},1967/1),({18,16,30,26},1619/1),({18,16,29,27},1096/1),({18,16,28,28},368/1),({17,17,44,12},1/1),({17,17,43,13},2/1),({17,17,42,14},8/1),({17,17,41,15},18/1),({17,17,40,16},44/1),({17,17,39,17},82/1),({17,17,38,18},152/1),({17,17,37,19},239/1),({17,17,36,20},366/1),({17,17,35,21},489/1),({17,17,34,22},632/1),({17,17,33,23},718/1),({17,17,32,24},779/1),({17,17,31,25},725/1),({17,17,30,26},622/1),({17,17,29,27},394/1),({17,17,28,28},154/1)}, (12,0) => {}, (10,2) => {}, (12,1) => {({30,10,38,28},1/1),({30,10,37,29},1/1),({30,10,36,30},3/1),({30,10,35,31},2/1),({30,10,34,32},2/1),({29,11,41,25},1/1),({29,11,40,26},3/1),({29,11,39,27},7/1),({29,11,38,28},11/1),({29,11,37,29},15/1),({29,11,36,30},17/1),({29,11,35,31},17/1),({29,11,34,32},12/1),({29,11,33,33},4/1),({28,12,43,23},1/1),({28,12,42,24},5/1),({28,12,41,25},11/1),({28,12,40,26},24/1),({28,12,39,27},39/1),({28,12,38,28},58/1),({28,12,37,29},70/1),({28,12,36,30},79/1),({28,12,35,31},68/1),({28,12,34,32},51/1),({28,12,33,33},16/1),({27,13,45,21},1/1),({27,13,44,22},3/1),({27,13,43,23},11/1),({27,13,42,24},27/1),({27,13,41,25},56/1),({27,13,40,26},95/1),({27,13,39,27},146/1),({27,13,38,28},194/1),({27,13,37,29},231/1),({27,13,36,30},237/1),({27,13,35,31},211/1),({27,13,34,32},144/1),({27,13,33,33},53/1),({26,14,46,20},1/1),({26,14,45,21},5/1),({26,14,44,22},19/1),({26,14,43,23},44/1),({26,14,42,24},95/1),({26,14,41,25},169/1),({26,14,40,26},273/1),({26,14,39,27},382/1),({26,14,38,28},495/1),({26,14,37,29},557/1),({26,14,36,30},570/1),({26,14,35,31},486/1),({26,14,34,32},339/1),({26,14,33,33},113/1),({25,15,47,19},1/1),({25,15,46,20},6/1),({25,15,45,21},21/1),({25,15,44,22},54/1),({25,15,43,23},121/1),({25,15,42,24},228/1),({25,15,41,25},386/1),({25,15,40,26},580/1),({25,15,39,27},795/1),({25,15,38,28},976/1),({25,15,37,29},1092/1),({25,15,36,30},1077/1),({25,15,35,31},924/1),({25,15,34,32},619/1),({25,15,33,33},224/1),({24,16,48,18},1/1),({24,16,47,19},5/1),({24,16,46,20},18/1),({24,16,45,21},50/1),({24,16,44,22},119/1),({24,16,43,23},236/1),({24,16,42,24},426/1),({24,16,41,25},677/1),({24,16,40,26},993/1),({24,16,39,27},1306/1),({24,16,38,28},1586/1),({24,16,37,29},1721/1),({24,16,36,30},1695/1),({24,16,35,31},1420/1),({24,16,34,32},966/1),({24,16,33,33},329/1),({23,17,48,18},2/1),({23,17,47,19},11/1),({23,17,46,20},33/1),({23,17,45,21},86/1),({23,17,44,22},185/1),({23,17,43,23},357/1),({23,17,42,24},608/1),({23,17,41,25},952/1),({23,17,40,26},1342/1),({23,17,39,27},1750/1),({23,17,38,28},2070/1),({23,17,37,29},2243/1),({23,17,36,30},2158/1),({23,17,35,31},1826/1),({23,17,34,32},1208/1),({23,17,33,33},433/1),({22,18,49,17},1/1),({22,18,48,18},5/1),({22,18,47,19},15/1),({22,18,46,20},46/1),({22,18,45,21},106/1),({22,18,44,22},223/1),({22,18,43,23},407/1),({22,18,42,24},687/1),({22,18,41,25},1033/1),({22,18,40,26},1453/1),({22,18,39,27},1846/1),({22,18,38,28},2180/1),({22,18,37,29},2315/1),({22,18,36,30},2245/1),({22,18,35,31},1854/1),({22,18,34,32},1254/1),({22,18,33,33},425/1),({21,19,49,17},1/1),({21,19,48,18},4/1),({21,19,47,19},15/1),({21,19,46,20},39/1),({21,19,45,21},92/1),({21,19,44,22},182/1),({21,19,43,23},335/1),({21,19,42,24},543/1),({21,19,41,25},824/1),({21,19,40,26},1126/1),({21,19,39,27},1442/1),({21,19,38,28},1666/1),({21,19,37,29},1788/1),({21,19,36,30},1694/1),({21,19,35,31},1429/1),({21,19,34,32},932/1),({21,19,33,33},341/1),({20,20,48,18},2/1),({20,20,47,19},6/1),({20,20,46,20},17/1),({20,20,45,21},35/1),({20,20,44,22},75/1),({20,20,43,23},127/1),({20,20,42,24},214/1),({20,20,41,25},311/1),({20,20,40,26},437/1),({20,20,39,27},537/1),({20,20,38,28},642/1),({20,20,37,29},661/1),({20,20,36,30},651/1),({20,20,35,31},523/1),({20,20,34,32},365/1),({20,20,33,33},110/1)}, (14,0) => {}, (12,2) => {}, (14,1) => {({32,14,43,33},2/1),({32,14,42,34},2/1),({32,14,41,35},3/1),({32,14,40,36},3/1),({32,14,39,37},3/1),({31,15,46,30},1/1),({31,15,45,31},3/1),({31,15,44,32},8/1),({31,15,43,33},12/1),({31,15,42,34},17/1),({31,15,41,35},19/1),({31,15,40,36},19/1),({31,15,39,37},13/1),({31,15,38,38},5/1),({30,16,48,28},1/1),({30,16,47,29},5/1),({30,16,46,30},11/1),({30,16,45,31},23/1),({30,16,44,32},37/1),({30,16,43,33},56/1),({30,16,42,34},67/1),({30,16,41,35},75/1),({30,16,40,36},65/1),({30,16,39,37},49/1),({30,16,38,38},15/1),({29,17,50,26},1/1),({29,17,49,27},3/1),({29,17,48,28},10/1),({29,17,47,29},23/1),({29,17,46,30},48/1),({29,17,45,31},80/1),({29,17,44,32},122/1),({29,17,43,33},162/1),({29,17,42,34},193/1),({29,17,41,35},197/1),({29,17,40,36},177/1),({29,17,39,37},119/1),({29,17,38,38},44/1),({28,18,51,25},1/1),({28,18,50,26},4/1),({28,18,49,27},14/1),({28,18,48,28},33/1),({28,18,47,29},71/1),({28,18,46,30},124/1),({28,18,45,31},201/1),({28,18,44,32},280/1),({28,18,43,33},362/1),({28,18,42,34},406/1),({28,18,41,35},417/1),({28,18,40,36},353/1),({28,18,39,37},248/1),({28,18,38,38},82/1),({27,19,52,24},1/1),({27,19,51,25},4/1),({27,19,50,26},14/1),({27,19,49,27},35/1),({27,19,48,28},78/1),({27,19,47,29},145/1),({27,19,46,30},247/1),({27,19,45,31},369/1),({27,19,44,32},508/1),({27,19,43,33},621/1),({27,19,42,34},697/1),({27,19,41,35},684/1),({27,19,40,36},591/1),({27,19,39,37},392/1),({27,19,38,38},144/1),({26,20,52,24},2/1),({26,20,51,25},9/1),({26,20,50,26},25/1),({26,20,49,27},63/1),({26,20,48,28},125/1),({26,20,47,29},229/1),({26,20,46,30},365/1),({26,20,45,31},539/1),({26,20,44,32},707/1),({26,20,43,33},866/1),({26,20,42,34},936/1),({26,20,41,35},926/1),({26,20,40,36},773/1),({26,20,39,37},531/1),({26,20,38,38},176/1),({25,21,53,23},1/1),({25,21,52,24},4/1),({25,21,51,25},13/1),({25,21,50,26},36/1),({25,21,49,27},78/1),({25,21,48,28},155/1),({25,21,47,29},266/1),({25,21,46,30},422/1),({25,21,45,31},597/1),({25,21,44,32},787/1),({25,21,43,33},928/1),({25,21,42,34},1015/1),({25,21,41,35},974/1),({25,21,40,36},829/1),({25,21,39,37},544/1),({25,21,38,38},201/1),({24,22,53,23},1/1),({24,22,52,24},4/1),({24,22,51,25},13/1),({24,22,50,26},31/1),({24,22,49,27},69/1),({24,22,48,28},128/1),({24,22,47,29},223/1),({24,22,46,30},338/1),({24,22,45,31},485/1),({24,22,44,32},616/1),({24,22,43,33},740/1),({24,22,42,34},782/1),({24,22,41,35},769/1),({24,22,40,36},629/1),({24,22,39,37},434/1),({24,22,38,38},140/1),({23,23,52,24},2/1),({23,23,51,25},4/1),({23,23,50,26},13/1),({23,23,49,27},26/1),({23,23,48,28},53/1),({23,23,47,29},83/1),({23,23,46,30},136/1),({23,23,45,31},181/1),({23,23,44,32},243/1),({23,23,43,33},274/1),({23,23,42,34},307/1),({23,23,41,35},280/1),({23,23,40,36},251/1),({23,23,39,37},151/1),({23,23,38,38},65/1)}, (16,0) => {}, (14,2) => {}, (16,1) => {({34,18,46,40},1/1),({34,18,45,41},1/1),({34,18,44,42},1/1),({33,19,49,37},2/1),({33,19,48,38},3/1),({33,19,47,39},5/1),({33,19,46,40},6/1),({33,19,45,41},7/1),({33,19,44,42},4/1),({33,19,43,43},2/1),({32,20,52,34},1/1),({32,20,51,35},2/1),({32,20,50,36},6/1),({32,20,49,37},10/1),({32,20,48,38},17/1),({32,20,47,39},20/1),({32,20,46,40},25/1),({32,20,45,41},21/1),({32,20,44,42},17/1),({32,20,43,43},4/1),({31,21,53,33},2/1),({31,21,52,34},5/1),({31,21,51,35},12/1),({31,21,50,36},21/1),({31,21,49,37},35/1),({31,21,48,38},47/1),({31,21,47,39},59/1),({31,21,46,40},60/1),({31,21,45,41},56/1),({31,21,44,42},37/1),({31,21,43,43},15/1),({30,22,55,31},1/1),({30,22,54,32},3/1),({30,22,53,33},7/1),({30,22,52,34},17/1),({30,22,51,35},31/1),({30,22,50,36},53/1),({30,22,49,37},75/1),({30,22,48,38},101/1),({30,22,47,39},113/1),({30,22,46,40},120/1),({30,22,45,41},100/1),({30,22,44,42},73/1),({30,22,43,43},22/1),({29,23,55,31},2/1),({29,23,54,32},6/1),({29,23,53,33},16/1),({29,23,52,34},31/1),({29,23,51,35},57/1),({29,23,50,36},86/1),({29,23,49,37},124/1),({29,23,48,38},152/1),({29,23,47,39},176/1),({29,23,46,40},171/1),({29,23,45,41},152/1),({29,23,44,42},98/1),({29,23,43,43},39/1),({28,24,56,30},1/1),({28,24,55,31},3/1),({28,24,54,32},10/1),({28,24,53,33},21/1),({28,24,52,34},43/1),({28,24,51,35},70/1),({28,24,50,36},110/1),({28,24,49,37},144/1),({28,24,48,38},183/1),({28,24,47,39},196/1),({28,24,46,40},201/1),({28,24,45,41},163/1),({28,24,44,42},118/1),({28,24,43,43},35/1),({27,25,56,30},1/1),({27,25,55,31},4/1),({27,25,54,32},9/1),({27,25,53,33},21/1),({27,25,52,34},37/1),({27,25,51,35},64/1),({27,25,50,36},90/1),({27,25,49,37},125/1),({27,25,48,38},146/1),({27,25,47,39},167/1),({27,25,46,40},155/1),({27,25,45,41},139/1),({27,25,44,42},86/1),({27,25,43,43},37/1),({26,26,55,31},1/1),({26,26,54,32},4/1),({26,26,53,33},7/1),({26,26,52,34},16/1),({26,26,51,35},23/1),({26,26,50,36},38/1),({26,26,49,37},45/1),({26,26,48,38},61/1),({26,26,47,39},59/1),({26,26,46,40},66/1),({26,26,45,41},47/1),({26,26,44,42},40/1),({26,26,43,43},7/1)}, (18,0) => {}, (16,2) => {}, (18,1) => {({35,23,50,46},1/1),({34,24,53,43},1/1),({34,24,52,44},1/1),({34,24,51,45},2/1),({34,24,50,46},1/1),({34,24,49,47},2/1),({33,25,56,40},1/1),({33,25,55,41},1/1),({33,25,54,42},3/1),({33,25,53,43},3/1),({33,25,52,44},5/1),({33,25,51,45},4/1),({33,25,50,46},5/1),({33,25,49,47},2/1),({33,25,48,48},2/1),({32,26,57,39},1/1),({32,26,56,40},2/1),({32,26,55,41},4/1),({32,26,54,42},5/1),({32,26,53,43},9/1),({32,26,52,44},8/1),({32,26,51,45},10/1),({32,26,50,46},7/1),({32,26,49,47},7/1),({31,27,58,38},1/1),({31,27,57,39},2/1),({31,27,56,40},4/1),({31,27,55,41},6/1),({31,27,54,42},10/1),({31,27,53,43},10/1),({31,27,52,44},13/1),({31,27,51,45},11/1),({31,27,50,46},11/1),({31,27,49,47},5/1),({31,27,48,48},4/1),({30,28,58,38},1/1),({30,28,57,39},3/1),({30,28,56,40},4/1),({30,28,55,41},7/1),({30,28,54,42},8/1),({30,28,53,43},12/1),({30,28,52,44},10/1),({30,28,51,45},12/1),({30,28,50,46},7/1),({30,28,49,47},8/1),({29,29,56,40},2/1),({29,29,55,41},1/1),({29,29,54,42},4/1),({29,29,53,43},3/1),({29,29,52,44},6/1),({29,29,51,45},2/1),({29,29,50,46},6/1),({29,29,48,48},3/1)}, (20,0) => {}, (18,2) => {}, (20,1) => {({33,31,59,47},1/1)}, (20,2) => {({35,32,58,53},1/1),({35,32,57,54},1/1)}, (22,2) => {}, (1,0) => {({4,0,5,1},1/1),({3,1,6,0},1/1),({3,1,5,1},1/1)}, (1,1) => {}, (3,0) => {({8,2,14,2},1/1),({8,2,12,4},1/1),({8,2,11,5},1/1),({8,2,10,6},1/1),({8,2,8,8},1/1),({7,3,13,3},1/1),({7,3,12,4},1/1),({7,3,11,5},1/1),({7,3,10,6},2/1),({7,3,9,7},1/1),({6,4,14,2},1/1),({6,4,13,3},1/1),({6,4,12,4},2/1),({6,4,11,5},2/1),({6,4,10,6},3/1),({6,4,9,7},1/1),({6,4,8,8},1/1),({5,5,11,5},1/1),({5,5,9,7},1/1)}, (3,1) => {({13,0,13,8},1/1),({13,0,12,9},1/1),({12,1,16,5},1/1),({12,1,15,6},1/1),({12,1,14,7},2/1),({12,1,13,8},2/1),({12,1,12,9},2/1),({12,1,11,10},1/1),({11,2,16,5},1/1),({11,2,15,6},2/1),({11,2,14,7},3/1),({11,2,13,8},3/1),({11,2,12,9},3/1),({11,2,11,10},2/1),({10,3,18,3},1/1),({10,3,17,4},1/1),({10,3,16,5},3/1),({10,3,15,6},4/1),({10,3,14,7},5/1),({10,3,13,8},5/1),({10,3,12,9},5/1),({10,3,11,10},2/1),({9,4,18,3},1/1),({9,4,17,4},2/1),({9,4,16,5},3/1),({9,4,15,6},4/1),({9,4,14,7},5/1),({9,4,13,8},5/1),({9,4,12,9},4/1),({9,4,11,10},2/1),({8,5,20,1},1/1),({8,5,19,2},1/1),({8,5,18,3},2/1),({8,5,17,4},3/1),({8,5,16,5},4/1),({8,5,15,6},4/1),({8,5,14,7},5/1),({8,5,13,8},4/1),({8,5,12,9},3/1),({8,5,11,10},2/1),({7,6,19,2},1/1),({7,6,18,3},1/1),({7,6,17,4},1/1),({7,6,16,5},2/1),({7,6,15,6},3/1),({7,6,14,7},2/1),({7,6,13,8},2/1),({7,6,12,9},2/1),({7,6,11,10},1/1)}, (3,2) => {}, (5,0) => {}, (5,1) => {({18,1,19,12},1/1),({18,1,18,13},1/1),({18,1,17,14},1/1),({17,2,22,9},1/1),({17,2,21,10},2/1),({17,2,20,11},3/1),({17,2,19,12},5/1),({17,2,18,13},6/1),({17,2,17,14},5/1),({17,2,16,15},3/1),({16,3,23,8},2/1),({16,3,22,9},4/1),({16,3,21,10},9/1),({16,3,20,11},13/1),({16,3,19,12},17/1),({16,3,18,13},18/1),({16,3,17,14},16/1),({16,3,16,15},8/1),({15,4,25,6},1/1),({15,4,24,7},3/1),({15,4,23,8},8/1),({15,4,22,9},16/1),({15,4,21,10},27/1),({15,4,20,11},37/1),({15,4,19,12},44/1),({15,4,18,13},45/1),({15,4,17,14},37/1),({15,4,16,15},21/1),({14,5,25,6},4/1),({14,5,24,7},9/1),({14,5,23,8},21/1),({14,5,22,9},37/1),({14,5,21,10},57/1),({14,5,20,11},73/1),({14,5,19,12},86/1),({14,5,18,13},83/1),({14,5,17,14},67/1),({14,5,16,15},38/1),({13,6,27,4},1/1),({13,6,26,5},4/1),({13,6,25,6},10/1),({13,6,24,7},23/1),({13,6,23,8},42/1),({13,6,22,9},67/1),({13,6,21,10},95/1),({13,6,20,11},120/1),({13,6,19,12},132/1),({13,6,18,13},127/1),({13,6,17,14},101/1),({13,6,16,15},55/1),({12,7,27,4},2/1),({12,7,26,5},6/1),({12,7,25,6},16/1),({12,7,24,7},32/1),({12,7,23,8},58/1),({12,7,22,9},88/1),({12,7,21,10},122/1),({12,7,20,11},148/1),({12,7,19,12},162/1),({12,7,18,13},152/1),({12,7,17,14},120/1),({12,7,16,15},66/1),({11,8,28,3},1/1),({11,8,27,4},3/1),({11,8,26,5},8/1),({11,8,25,6},19/1),({11,8,24,7},36/1),({11,8,23,8},60/1),({11,8,22,9},90/1),({11,8,21,10},119/1),({11,8,20,11},142/1),({11,8,19,12},153/1),({11,8,18,13},142/1),({11,8,17,14},110/1),({11,8,16,15},62/1),({10,9,27,4},2/1),({10,9,26,5},6/1),({10,9,25,6},12/1),({10,9,24,7},23/1),({10,9,23,8},39/1),({10,9,22,9},55/1),({10,9,21,10},74/1),({10,9,20,11},88/1),({10,9,19,12},92/1),({10,9,18,13},86/1),({10,9,17,14},68/1),({10,9,16,15},36/1)}};
end;