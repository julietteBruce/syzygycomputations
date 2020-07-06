A := QQ[t_0,t_1,t_2,t_3];
--tb stands for Total Betti numbers
tb0227 = new HashTable from {(5,2) => 0, (7,0) => 0, (7,1) => 1614354/1, (9,0) => 0, (7,2) => 0, (9,1) => 4300422/1, (11,0) => 0, (9,2) => 0, (11,1) => 5114382/1, (13,0) => 0, (11,2) => 0, (13,1) => 2748730/1, (15,0) => 0, (13,2) => 0, (15,1) => 613377/1, (17,0) => 0, (15,2) => 0, (17,1) => 40299/1, (19,0) => 0, (17,2) => 0, (19,1) => 357/1, (21,0) => 0, (19,2) => 420/1, (21,1) => 0, (21,2) => 4/1, (0,0) => 3/1, (2,0) => 210/1, (2,1) => 336/1, (2,2) => 0, (4,0) => 0, (4,1) => 59318/1, (4,2) => 0, (6,0) => 0, (6,1) => 709308/1, (6,2) => 0, (8,0) => 0, (8,1) => 2926380/1, (10,0) => 0, (8,2) => 0, (10,1) => 5173168/1, (12,0) => 0, (10,2) => 0, (12,1) => 4151196/1, (14,0) => 0, (12,2) => 0, (14,1) => 1465128/1, (16,0) => 0, (14,2) => 0, (16,1) => 192318/1, (18,0) => 0, (16,2) => 0, (18,1) => 3360/1, (20,0) => 0, (18,2) => 1330/1, (20,1) => 18/1, (20,2) => 63/1, (22,2) => 0, (1,0) => 42/1, (1,1) => 17/1, (3,0) => 0, (3,1) => 9135/1, (3,2) => 0, (5,0) => 0, (5,1) => 240597/1};
--sb represents the betti numbers as sums of Schur functors
sb0227 = new HashTable from {(5,2) => {}, (7,0) => {}, (7,1) => {({15,1,34,24},1/1),({15,1,33,25},1/1),({15,1,32,26},1/1),({15,1,31,27},1/1),({15,1,30,28},1/1),({14,2,39,19},1/1),({14,2,38,20},2/1),({14,2,37,21},5/1),({14,2,36,22},7/1),({14,2,35,23},10/1),({14,2,34,24},12/1),({14,2,33,25},14/1),({14,2,32,26},13/1),({14,2,31,27},12/1),({14,2,30,28},7/1),({14,2,29,29},3/1),({13,3,42,16},2/1),({13,3,41,17},4/1),({13,3,40,18},10/1),({13,3,39,19},17/1),({13,3,38,20},29/1),({13,3,37,21},40/1),({13,3,36,22},56/1),({13,3,35,23},66/1),({13,3,34,24},77/1),({13,3,33,25},77/1),({13,3,32,26},75/1),({13,3,31,27},58/1),({13,3,30,28},41/1),({13,3,29,29},12/1),({12,4,45,13},1/1),({12,4,44,14},3/1),({12,4,43,15},9/1),({12,4,42,16},18/1),({12,4,41,17},35/1),({12,4,40,18},58/1),({12,4,39,19},91/1),({12,4,38,20},127/1),({12,4,37,21},172/1),({12,4,36,22},210/1),({12,4,35,23},246/1),({12,4,34,24},262/1),({12,4,33,25},265/1),({12,4,32,26},237/1),({12,4,31,27},195/1),({12,4,30,28},122/1),({12,4,29,29},45/1),({11,5,46,12},3/1),({11,5,45,13},7/1),({11,5,44,14},18/1),({11,5,43,15},36/1),({11,5,42,16},68/1),({11,5,41,17},110/1),({11,5,40,18},174/1),({11,5,39,19},245/1),({11,5,38,20},335/1),({11,5,37,21},419/1),({11,5,36,22},505/1),({11,5,35,23},559/1),({11,5,34,24},597/1),({11,5,33,25},575/1),({11,5,32,26},522/1),({11,5,31,27},407/1),({11,5,30,28},270/1),({11,5,29,29},86/1),({10,6,48,10},1/1),({10,6,47,11},4/1),({10,6,46,12},9/1),({10,6,45,13},23/1),({10,6,44,14},44/1),({10,6,43,15},82/1),({10,6,42,16},135/1),({10,6,41,17},214/1),({10,6,40,18},307/1),({10,6,39,19},429/1),({10,6,38,20},551/1),({10,6,37,21},684/1),({10,6,36,22},789/1),({10,6,35,23},877/1),({10,6,34,24},897/1),({10,6,33,25},877/1),({10,6,32,26},768/1),({10,6,31,27},614/1),({10,6,30,28},384/1),({10,6,29,29},141/1),({9,7,48,10},2/1),({9,7,47,11},5/1),({9,7,46,12},14/1),({9,7,45,13},27/1),({9,7,44,14},55/1),({9,7,43,15},92/1),({9,7,42,16},154/1),({9,7,41,17},227/1),({9,7,40,18},331/1),({9,7,39,19},438/1),({9,7,38,20},569/1),({9,7,37,21},679/1),({9,7,36,22},792/1),({9,7,35,23},849/1),({9,7,34,24},885/1),({9,7,33,25},834/1),({9,7,32,26},750/1),({9,7,31,27},574/1),({9,7,30,28},380/1),({9,7,29,29},119/1),({8,8,49,9},1/1),({8,8,48,10},1/1),({8,8,47,11},4/1),({8,8,46,12},7/1),({8,8,45,13},16/1),({8,8,44,14},26/1),({8,8,43,15},48/1),({8,8,42,16},69/1),({8,8,41,17},108/1),({8,8,40,18},145/1),({8,8,39,19},200/1),({8,8,38,20},242/1),({8,8,37,21},303/1),({8,8,36,22},333/1),({8,8,35,23},373/1),({8,8,34,24},368/1),({8,8,33,25},365/1),({8,8,32,26},306/1),({8,8,31,27},256/1),({8,8,30,28},148/1),({8,8,29,29},62/1)}, (9,0) => {}, (7,2) => {}, (9,1) => {({17,3,45,27},1/1),({17,3,44,28},1/1),({17,3,43,29},3/1),({17,3,42,30},4/1),({17,3,41,31},5/1),({17,3,40,32},5/1),({17,3,39,33},6/1),({17,3,38,34},4/1),({17,3,37,35},4/1),({17,3,36,36},1/1),({16,4,48,24},2/1),({16,4,47,25},5/1),({16,4,46,26},10/1),({16,4,45,27},16/1),({16,4,44,28},25/1),({16,4,43,29},34/1),({16,4,42,30},44/1),({16,4,41,31},49/1),({16,4,40,32},52/1),({16,4,39,33},48/1),({16,4,38,34},41/1),({16,4,37,35},26/1),({16,4,36,36},10/1),({15,5,51,21},2/1),({15,5,50,22},5/1),({15,5,49,23},14/1),({15,5,48,24},26/1),({15,5,47,25},48/1),({15,5,46,26},74/1),({15,5,45,27},111/1),({15,5,44,28},147/1),({15,5,43,29},188/1),({15,5,42,30},216/1),({15,5,41,31},240/1),({15,5,40,32},236/1),({15,5,39,33},221/1),({15,5,38,34},173/1),({15,5,37,35},117/1),({15,5,36,36},37/1),({14,6,53,19},2/1),({14,6,52,20},8/1),({14,6,51,21},19/1),({14,6,50,22},41/1),({14,6,49,23},75/1),({14,6,48,24},129/1),({14,6,47,25},199/1),({14,6,46,26},291/1),({14,6,45,27},392/1),({14,6,44,28},502/1),({14,6,43,29},598/1),({14,6,42,30},678/1),({14,6,41,31},710/1),({14,6,40,32},702/1),({14,6,39,33},625/1),({14,6,38,34},501/1),({14,6,37,35},318/1),({14,6,36,36},114/1),({13,7,55,17},2/1),({13,7,54,18},6/1),({13,7,53,19},17/1),({13,7,52,20},36/1),({13,7,51,21},75/1),({13,7,50,22},133/1),({13,7,49,23},226/1),({13,7,48,24},346/1),({13,7,47,25},510/1),({13,7,46,26},694/1),({13,7,45,27},909/1),({13,7,44,28},1108/1),({13,7,43,29},1298/1),({13,7,42,30},1416/1),({13,7,41,31},1476/1),({13,7,40,32},1413/1),({13,7,39,33},1265/1),({13,7,38,34},984/1),({13,7,37,35},641/1),({13,7,36,36},211/1),({12,8,56,16},2/1),({12,8,55,17},6/1),({12,8,54,18},18/1),({12,8,53,19},39/1),({12,8,52,20},81/1),({12,8,51,21},146/1),({12,8,50,22},251/1),({12,8,49,23},392/1),({12,8,48,24},589/1),({12,8,47,25},821/1),({12,8,46,26},1101/1),({12,8,45,27},1384/1),({12,8,44,28},1674/1),({12,8,43,29},1902/1),({12,8,42,30},2071/1),({12,8,41,31},2103/1),({12,8,40,32},2023/1),({12,8,39,33},1767/1),({12,8,38,34},1397/1),({12,8,37,35},878/1),({12,8,36,36},312/1),({11,9,57,15},2/1),({11,9,56,16},4/1),({11,9,55,17},12/1),({11,9,54,18},25/1),({11,9,53,19},54/1),({11,9,52,20},97/1),({11,9,51,21},173/1),({11,9,50,22},274/1),({11,9,49,23},427/1),({11,9,48,24},609/1),({11,9,47,25},846/1),({11,9,46,26},1094/1),({11,9,45,27},1376/1),({11,9,44,28},1617/1),({11,9,43,29},1843/1),({11,9,42,30},1960/1),({11,9,41,31},2005/1),({11,9,40,32},1885/1),({11,9,39,33},1670/1),({11,9,38,34},1283/1),({11,9,37,35},836/1),({11,9,36,36},272/1),({10,10,56,16},2/1),({10,10,55,17},4/1),({10,10,54,18},12/1),({10,10,53,19},22/1),({10,10,52,20},45/1),({10,10,51,21},73/1),({10,10,50,22},123/1),({10,10,49,23},180/1),({10,10,48,24},267/1),({10,10,47,25},354/1),({10,10,46,26},471/1),({10,10,45,27},570/1),({10,10,44,28},688/1),({10,10,43,29},758/1),({10,10,42,30},827/1),({10,10,41,31},818/1),({10,10,40,32},794/1),({10,10,39,33},675/1),({10,10,38,34},545/1),({10,10,37,35},328/1),({10,10,36,36},128/1)}, (11,0) => {}, (9,2) => {}, (11,1) => {({19,5,52,34},2/1),({19,5,51,35},2/1),({19,5,50,36},4/1),({19,5,49,37},5/1),({19,5,48,38},7/1),({19,5,47,39},7/1),({19,5,46,40},8/1),({19,5,45,41},5/1),({19,5,44,42},5/1),({19,5,43,43},1/1),({18,6,56,30},1/1),({18,6,55,31},3/1),({18,6,54,32},7/1),({18,6,53,33},14/1),({18,6,52,34},22/1),({18,6,51,35},33/1),({18,6,50,36},44/1),({18,6,49,37},55/1),({18,6,48,38},62/1),({18,6,47,39},65/1),({18,6,46,40},60/1),({18,6,45,41},50/1),({18,6,44,42},32/1),({18,6,43,43},12/1),({17,7,58,28},3/1),({17,7,57,29},8/1),({17,7,56,30},19/1),({17,7,55,31},35/1),({17,7,54,32},62/1),({17,7,53,33},94/1),({17,7,52,34},139/1),({17,7,51,35},182/1),({17,7,50,36},230/1),({17,7,49,37},263/1),({17,7,48,38},289/1),({17,7,47,39},284/1),({17,7,46,40},264/1),({17,7,45,41},207/1),({17,7,44,42},139/1),({17,7,43,43},44/1),({16,8,61,25},1/1),({16,8,60,26},4/1),({16,8,59,27},12/1),({16,8,58,28},27/1),({16,8,57,29},55/1),({16,8,56,30},98/1),({16,8,55,31},164/1),({16,8,54,32},249/1),({16,8,53,33},358/1),({16,8,52,34},477/1),({16,8,51,35},605/1),({16,8,50,36},716/1),({16,8,49,37},806/1),({16,8,48,38},841/1),({16,8,47,39},826/1),({16,8,46,40},735/1),({16,8,45,41},587/1),({16,8,44,42},373/1),({16,8,43,43},133/1),({15,9,62,24},3/1),({15,9,61,25},8/1),({15,9,60,26},23/1),({15,9,59,27},48/1),({15,9,58,28},96/1),({15,9,57,29},167/1),({15,9,56,30},280/1),({15,9,55,31},423/1),({15,9,54,32},616/1),({15,9,53,33},831/1),({15,9,52,34},1079/1),({15,9,51,35},1307/1),({15,9,50,36},1523/1),({15,9,49,37},1654/1),({15,9,48,38},1717/1),({15,9,47,39},1640/1),({15,9,46,40},1463/1),({15,9,45,41},1137/1),({15,9,44,42},740/1),({15,9,43,43},244/1),({14,10,64,22},1/1),({14,10,63,23},4/1),({14,10,62,24},10/1),({14,10,61,25},26/1),({14,10,60,26},54/1),({14,10,59,27},106/1),({14,10,58,28},187/1),({14,10,57,29},313/1),({14,10,56,30},482/1),({14,10,55,31},713/1),({14,10,54,32},985/1),({14,10,53,33},1306/1),({14,10,52,34},1633/1),({14,10,51,35},1960/1),({14,10,50,36},2217/1),({14,10,49,37},2401/1),({14,10,48,38},2433/1),({14,10,47,39},2330/1),({14,10,46,40},2034/1),({14,10,45,41},1602/1),({14,10,44,42},1009/1),({14,10,43,43},356/1),({13,11,64,22},2/1),({13,11,63,23},5/1),({13,11,62,24},15/1),({13,11,61,25},32/1),({13,11,60,26},67/1),({13,11,59,27},120/1),({13,11,58,28},211/1),({13,11,57,29},332/1),({13,11,56,30},510/1),({13,11,55,31},723/1),({13,11,54,32},995/1),({13,11,53,33},1281/1),({13,11,52,34},1599/1),({13,11,51,35},1873/1),({13,11,50,36},2123/1),({13,11,49,37},2253/1),({13,11,48,38},2295/1),({13,11,47,39},2156/1),({13,11,46,40},1904/1),({13,11,45,41},1464/1),({13,11,44,42},949/1),({13,11,43,43},311/1),({12,12,65,21},1/1),({12,12,64,22},1/1),({12,12,63,23},4/1),({12,12,62,24},8/1),({12,12,61,25},18/1),({12,12,60,26},31/1),({12,12,59,27},60/1),({12,12,58,28},94/1),({12,12,57,29},154/1),({12,12,56,30},222/1),({12,12,55,31},322/1),({12,12,54,32},424/1),({12,12,53,33},558/1),({12,12,52,34},671/1),({12,12,51,35},802/1),({12,12,50,36},882/1),({12,12,49,37},955/1),({12,12,48,38},943/1),({12,12,47,39},911/1),({12,12,46,40},775/1),({12,12,45,41},622/1),({12,12,44,42},376/1),({12,12,43,43},144/1)}, (13,0) => {}, (11,2) => {}, (13,1) => {({21,7,57,43},1/1),({21,7,56,44},1/1),({21,7,55,45},3/1),({21,7,54,46},2/1),({21,7,53,47},3/1),({21,7,52,48},2/1),({21,7,51,49},2/1),({20,8,61,39},2/1),({20,8,60,40},4/1),({20,8,59,41},7/1),({20,8,58,42},12/1),({20,8,57,43},17/1),({20,8,56,44},22/1),({20,8,55,45},26/1),({20,8,54,46},28/1),({20,8,53,47},26/1),({20,8,52,48},23/1),({20,8,51,49},15/1),({20,8,50,50},5/1),({19,9,65,35},1/1),({19,9,64,36},2/1),({19,9,63,37},7/1),({19,9,62,38},13/1),({19,9,61,39},25/1),({19,9,60,40},40/1),({19,9,59,41},62/1),({19,9,58,42},82/1),({19,9,57,43},108/1),({19,9,56,44},125/1),({19,9,55,45},141/1),({19,9,54,46},139/1),({19,9,53,47},132/1),({19,9,52,48},102/1),({19,9,51,49},71/1),({19,9,50,50},22/1),({18,10,67,33},1/1),({18,10,66,34},4/1),({18,10,65,35},10/1),({18,10,64,36},23/1),({18,10,63,37},43/1),({18,10,62,38},75/1),({18,10,61,39},118/1),({18,10,60,40},174/1),({18,10,59,41},237/1),({18,10,58,42},307/1),({18,10,57,43},368/1),({18,10,56,44},419/1),({18,10,55,45},442/1),({18,10,54,46},438/1),({18,10,53,47},391/1),({18,10,52,48},315/1),({18,10,51,49},200/1),({18,10,50,50},71/1),({17,11,69,31},1/1),({17,11,68,32},3/1),({17,11,67,33},10/1),({17,11,66,34},21/1),({17,11,65,35},45/1),({17,11,64,36},81/1),({17,11,63,37},140/1),({17,11,62,38},215/1),({17,11,61,39},322/1),({17,11,60,40},439/1),({17,11,59,41},579/1),({17,11,58,42},709/1),({17,11,57,43},835/1),({17,11,56,44},911/1),({17,11,55,45},956/1),({17,11,54,46},914/1),({17,11,53,47},821/1),({17,11,52,48},639/1),({17,11,51,49},418/1),({17,11,50,50},135/1),({16,12,70,30},1/1),({16,12,69,31},4/1),({16,12,68,32},11/1),({16,12,67,33},24/1),({16,12,66,34},51/1),({16,12,65,35},93/1),({16,12,64,36},160/1),({16,12,63,37},253/1),({16,12,62,38},382/1),({16,12,61,39},535/1),({16,12,60,40},721/1),({16,12,59,41},910/1),({16,12,58,42},1102/1),({16,12,57,43},1257/1),({16,12,56,44},1371/1),({16,12,55,45},1393/1),({16,12,54,46},1342/1),({16,12,53,47},1175/1),({16,12,52,48},927/1),({16,12,51,49},584/1),({16,12,50,50},208/1),({15,13,71,29},1/1),({15,13,70,30},2/1),({15,13,69,31},7/1),({15,13,68,32},15/1),({15,13,67,33},34/1),({15,13,66,34},62/1),({15,13,65,35},112/1),({15,13,64,36},179/1),({15,13,63,37},281/1),({15,13,62,38},402/1),({15,13,61,39},562/1),({15,13,60,40},728/1),({15,13,59,41},919/1),({15,13,58,42},1082/1),({15,13,57,43},1236/1),({15,13,56,44},1315/1),({15,13,55,45},1349/1),({15,13,54,46},1268/1),({15,13,53,47},1125/1),({15,13,52,48},864/1),({15,13,51,49},564/1),({15,13,50,50},182/1),({14,14,70,30},2/1),({14,14,69,31},3/1),({14,14,68,32},8/1),({14,14,67,33},15/1),({14,14,66,34},30/1),({14,14,65,35},48/1),({14,14,64,36},83/1),({14,14,63,37},120/1),({14,14,62,38},179/1),({14,14,61,39},239/1),({14,14,60,40},318/1),({14,14,59,41},384/1),({14,14,58,42},467/1),({14,14,57,43},513/1),({14,14,56,44},560/1),({14,14,55,45},556/1),({14,14,54,46},539/1),({14,14,53,47},456/1),({14,14,52,48},372/1),({14,14,51,49},222/1),({14,14,50,50},86/1)}, (15,0) => {}, (13,2) => {}, (15,1) => {({23,9,58,56},1/1),({22,10,64,50},1/1),({22,10,63,51},2/1),({22,10,62,52},3/1),({22,10,61,53},3/1),({22,10,60,54},3/1),({22,10,59,55},3/1),({22,10,58,56},2/1),({22,10,57,57},1/1),({21,11,68,46},2/1),({21,11,67,47},3/1),({21,11,66,48},7/1),({21,11,65,49},9/1),({21,11,64,50},15/1),({21,11,63,51},18/1),({21,11,62,52},23/1),({21,11,61,53},22/1),({21,11,60,54},23/1),({21,11,59,55},17/1),({21,11,58,56},14/1),({21,11,57,57},3/1),({20,12,72,42},1/1),({20,12,71,43},2/1),({20,12,70,44},5/1),({20,12,69,45},10/1),({20,12,68,46},17/1),({20,12,67,47},28/1),({20,12,66,48},41/1),({20,12,65,49},56/1),({20,12,64,50},70/1),({20,12,63,51},83/1),({20,12,62,52},90/1),({20,12,61,53},92/1),({20,12,60,54},84/1),({20,12,59,55},68/1),({20,12,58,56},43/1),({20,12,57,57},16/1),({19,13,74,40},1/1),({19,13,73,41},2/1),({19,13,72,42},6/1),({19,13,71,43},12/1),({19,13,70,44},24/1),({19,13,69,45},39/1),({19,13,68,46},63/1),({19,13,67,47},89/1),({19,13,66,48},125/1),({19,13,65,49},156/1),({19,13,64,50},191/1),({19,13,63,51},211/1),({19,13,62,52},228/1),({19,13,61,53},219/1),({19,13,60,54},202/1),({19,13,59,55},155/1),({19,13,58,56},105/1),({19,13,57,57},32/1),({18,14,75,39},1/1),({18,14,74,40},3/1),({18,14,73,41},8/1),({18,14,72,42},16/1),({18,14,71,43},31/1),({18,14,70,44},52/1),({18,14,69,45},84/1),({18,14,68,46},122/1),({18,14,67,47},171/1),({18,14,66,48},221/1),({18,14,65,49},276/1),({18,14,64,50},319/1),({18,14,63,51},355/1),({18,14,62,52},364/1),({18,14,61,53},356/1),({18,14,60,54},312/1),({18,14,59,55},250/1),({18,14,58,56},156/1),({18,14,57,57},57/1),({17,15,76,38},1/1),({17,15,75,39},2/1),({17,15,74,40},6/1),({17,15,73,41},11/1),({17,15,72,42},23/1),({17,15,71,43},38/1),({17,15,70,44},65/1),({17,15,69,45},95/1),({17,15,68,46},140/1),({17,15,67,47},184/1),({17,15,66,48},241/1),({17,15,65,49},286/1),({17,15,64,50},336/1),({17,15,63,51},358/1),({17,15,62,52},375/1),({17,15,61,53},351/1),({17,15,60,54},318/1),({17,15,59,55},241/1),({17,15,58,56},162/1),({17,15,57,57},49/1),({16,16,75,39},1/1),({16,16,74,40},2/1),({16,16,73,41},6/1),({16,16,72,42},10/1),({16,16,71,43},19/1),({16,16,70,44},28/1),({16,16,69,45},45/1),({16,16,68,46},61/1),({16,16,67,47},85/1),({16,16,66,48},103/1),({16,16,65,49},129/1),({16,16,64,50},142/1),({16,16,63,51},160/1),({16,16,62,52},156/1),({16,16,61,53},156/1),({16,16,60,54},130/1),({16,16,59,55},109/1),({16,16,58,56},63/1),({16,16,57,57},27/1)}, (17,0) => {}, (15,2) => {}, (17,1) => {({23,13,69,59},1/1),({23,13,67,61},1/1),({23,13,65,63},1/1),({22,14,73,55},1/1),({22,14,72,56},1/1),({22,14,71,57},2/1),({22,14,70,58},3/1),({22,14,69,59},3/1),({22,14,68,60},4/1),({22,14,67,61},4/1),({22,14,66,62},3/1),({22,14,65,63},2/1),({22,14,64,64},1/1),({21,15,77,51},1/1),({21,15,76,52},1/1),({21,15,75,53},3/1),({21,15,74,54},4/1),({21,15,73,55},7/1),({21,15,72,56},8/1),({21,15,71,57},12/1),({21,15,70,58},12/1),({21,15,69,59},15/1),({21,15,68,60},13/1),({21,15,67,61},14/1),({21,15,66,62},9/1),({21,15,65,63},8/1),({21,15,64,64},1/1),({20,16,79,49},1/1),({20,16,78,50},2/1),({20,16,77,51},3/1),({20,16,76,52},6/1),({20,16,75,53},9/1),({20,16,74,54},13/1),({20,16,73,55},17/1),({20,16,72,56},22/1),({20,16,71,57},25/1),({20,16,70,58},29/1),({20,16,69,59},29/1),({20,16,68,60},29/1),({20,16,67,61},25/1),({20,16,66,62},21/1),({20,16,65,63},12/1),({20,16,64,64},5/1),({19,17,79,49},1/1),({19,17,78,50},2/1),({19,17,77,51},5/1),({19,17,76,52},7/1),({19,17,75,53},12/1),({19,17,74,54},16/1),({19,17,73,55},23/1),({19,17,72,56},26/1),({19,17,71,57},32/1),({19,17,70,58},32/1),({19,17,69,59},36/1),({19,17,68,60},31/1),({19,17,67,61},30/1),({19,17,66,62},20/1),({19,17,65,63},16/1),({19,17,64,64},3/1),({18,18,80,48},1/1),({18,18,79,49},1/1),({18,18,78,50},2/1),({18,18,77,51},3/1),({18,18,76,52},5/1),({18,18,75,53},6/1),({18,18,74,54},10/1),({18,18,73,55},10/1),({18,18,72,56},14/1),({18,18,71,57},14/1),({18,18,70,58},17/1),({18,18,69,59},15/1),({18,18,68,60},17/1),({18,18,67,61},12/1),({18,18,66,62},12/1),({18,18,65,63},5/1),({18,18,64,64},4/1)}, (19,0) => {}, (17,2) => {}, (19,1) => {({21,19,82,60},1/1),({21,19,81,61},1/1),({21,19,80,62},1/1),({21,19,79,63},1/1),({21,19,78,64},1/1),({21,19,77,65},1/1),({21,19,76,66},1/1)}, (21,0) => {}, (19,2) => {({23,19,80,69},1/1),({23,19,79,70},1/1),({23,19,78,71},1/1),({23,19,77,72},1/1),({23,19,76,73},1/1),({23,19,75,74},1/1),({22,20,81,68},1/1),({22,20,80,69},1/1),({22,20,79,70},1/1),({22,20,78,71},1/1),({22,20,77,72},1/1),({22,20,76,73},1/1),({22,20,75,74},1/1),({21,21,80,69},1/1),({21,21,79,70},1/1),({21,21,78,71},1/1),({21,21,77,72},1/1),({21,21,76,73},1/1),({21,21,75,74},1/1)}, (21,1) => {}, (21,2) => {({23,23,83,80},1/1)}, (0,0) => {({0,0,2,0},1/1)}, (2,0) => {({4,0,13,3},1/1),({4,0,11,5},1/1),({4,0,9,7},1/1),({3,1,14,2},1/1),({3,1,12,4},1/1),({3,1,10,6},1/1),({3,1,8,8},1/1),({2,2,13,3},1/1),({2,2,11,5},1/1),({2,2,9,7},1/1)}, (2,1) => {({4,2,22,1},1/1),({4,2,21,2},1/1),({4,2,20,3},1/1),({4,2,19,4},1/1),({4,2,18,5},1/1),({4,2,17,6},1/1),({4,2,16,7},1/1)}, (2,2) => {}, (4,0) => {}, (4,1) => {({10,0,24,13},1/1),({10,0,23,14},1/1),({10,0,22,15},1/1),({10,0,21,16},1/1),({10,0,20,17},1/1),({10,0,19,18},1/1),({9,1,28,9},1/1),({9,1,27,10},2/1),({9,1,26,11},3/1),({9,1,25,12},5/1),({9,1,24,13},6/1),({9,1,23,14},7/1),({9,1,22,15},8/1),({9,1,21,16},7/1),({9,1,20,17},5/1),({9,1,19,18},3/1),({8,2,31,6},1/1),({8,2,30,7},2/1),({8,2,29,8},4/1),({8,2,28,9},7/1),({8,2,27,10},11/1),({8,2,26,11},15/1),({8,2,25,12},20/1),({8,2,24,13},23/1),({8,2,23,14},25/1),({8,2,22,15},25/1),({8,2,21,16},22/1),({8,2,20,17},16/1),({8,2,19,18},9/1),({7,3,33,4},1/1),({7,3,32,5},2/1),({7,3,31,6},4/1),({7,3,30,7},8/1),({7,3,29,8},13/1),({7,3,28,9},19/1),({7,3,27,10},27/1),({7,3,26,11},35/1),({7,3,25,12},42/1),({7,3,24,13},47/1),({7,3,23,14},49/1),({7,3,22,15},47/1),({7,3,21,16},41/1),({7,3,20,17},30/1),({7,3,19,18},15/1),({6,4,33,4},1/1),({6,4,32,5},3/1),({6,4,31,6},6/1),({6,4,30,7},10/1),({6,4,29,8},16/1),({6,4,28,9},24/1),({6,4,27,10},33/1),({6,4,26,11},41/1),({6,4,25,12},48/1),({6,4,24,13},53/1),({6,4,23,14},55/1),({6,4,22,15},52/1),({6,4,21,16},44/1),({6,4,20,17},32/1),({6,4,19,18},17/1),({5,5,34,3},1/1),({5,5,33,4},1/1),({5,5,32,5},2/1),({5,5,31,6},4/1),({5,5,30,7},6/1),({5,5,29,8},9/1),({5,5,28,9},13/1),({5,5,27,10},16/1),({5,5,26,11},20/1),({5,5,25,12},23/1),({5,5,24,13},25/1),({5,5,23,14},25/1),({5,5,22,15},24/1),({5,5,21,16},20/1),({5,5,20,17},14/1),({5,5,19,18},8/1)}, (4,2) => {}, (6,0) => {}, (6,1) => {({14,0,27,24},1/1),({13,1,33,18},1/1),({13,1,32,19},2/1),({13,1,31,20},3/1),({13,1,30,21},4/1),({13,1,29,22},4/1),({13,1,28,23},4/1),({13,1,27,24},3/1),({13,1,26,25},2/1),({12,2,37,14},1/1),({12,2,36,15},3/1),({12,2,35,16},6/1),({12,2,34,17},10/1),({12,2,33,18},15/1),({12,2,32,19},21/1),({12,2,31,20},26/1),({12,2,30,21},29/1),({12,2,29,22},29/1),({12,2,28,23},26/1),({12,2,27,24},20/1),({12,2,26,25},11/1),({11,3,40,11},1/1),({11,3,39,12},3/1),({11,3,38,13},7/1),({11,3,37,14},14/1),({11,3,36,15},25/1),({11,3,35,16},39/1),({11,3,34,17},57/1),({11,3,33,18},76/1),({11,3,32,19},93/1),({11,3,31,20},107/1),({11,3,30,21},114/1),({11,3,29,22},110/1),({11,3,28,23},96/1),({11,3,27,24},71/1),({11,3,26,25},37/1),({10,4,42,9},1/1),({10,4,41,10},3/1),({10,4,40,11},8/1),({10,4,39,12},17/1),({10,4,38,13},32/1),({10,4,37,14},54/1),({10,4,36,15},84/1),({10,4,35,16},121/1),({10,4,34,17},163/1),({10,4,33,18},205/1),({10,4,32,19},242/1),({10,4,31,20},267/1),({10,4,30,21},275/1),({10,4,29,22},261/1),({10,4,28,23},223/1),({10,4,27,24},163/1),({10,4,26,25},86/1),({9,5,43,8},1/1),({9,5,42,9},4/1),({9,5,41,10},10/1),({9,5,40,11},21/1),({9,5,39,12},40/1),({9,5,38,13},69/1),({9,5,37,14},109/1),({9,5,36,15},160/1),({9,5,35,16},220/1),({9,5,34,17},285/1),({9,5,33,18},348/1),({9,5,32,19},400/1),({9,5,31,20},432/1),({9,5,30,21},437/1),({9,5,29,22},409/1),({9,5,28,23},346/1),({9,5,27,24},251/1),({9,5,26,25},132/1),({8,6,44,7},1/1),({8,6,43,8},3/1),({8,6,42,9},7/1),({8,6,41,10},15/1),({8,6,40,11},29/1),({8,6,39,12},51/1),({8,6,38,13},83/1),({8,6,37,14},125/1),({8,6,36,15},178/1),({8,6,35,16},238/1),({8,6,34,17},302/1),({8,6,33,18},362/1),({8,6,32,19},410/1),({8,6,31,20},438/1),({8,6,30,21},439/1),({8,6,29,22},407/1),({8,6,28,23},343/1),({8,6,27,24},248/1),({8,6,26,25},130/1),({7,7,43,8},1/1),({7,7,42,9},3/1),({7,7,41,10},7/1),({7,7,40,11},13/1),({7,7,39,12},23/1),({7,7,38,13},37/1),({7,7,37,14},56/1),({7,7,36,15},79/1),({7,7,35,16},105/1),({7,7,34,17},132/1),({7,7,33,18},157/1),({7,7,32,19},177/1),({7,7,31,20},188/1),({7,7,30,21},187/1),({7,7,29,22},173/1),({7,7,28,23},145/1),({7,7,27,24},105/1),({7,7,26,25},55/1)}, (6,2) => {}, (8,0) => {}, (8,1) => {({16,2,40,25},1/1),({16,2,39,26},1/1),({16,2,38,27},2/1),({16,2,37,28},3/1),({16,2,36,29},3/1),({16,2,35,30},3/1),({16,2,34,31},2/1),({16,2,33,32},1/1),({15,3,44,21},1/1),({15,3,43,22},3/1),({15,3,42,23},6/1),({15,3,41,24},11/1),({15,3,40,25},16/1),({15,3,39,26},22/1),({15,3,38,27},27/1),({15,3,37,28},30/1),({15,3,36,29},31/1),({15,3,35,30},28/1),({15,3,34,31},21/1),({15,3,33,32},11/1),({14,4,47,18},1/1),({14,4,46,19},4/1),({14,4,45,20},10/1),({14,4,44,21},19/1),({14,4,43,22},34/1),({14,4,42,23},54/1),({14,4,41,24},78/1),({14,4,40,25},104/1),({14,4,39,26},128/1),({14,4,38,27},147/1),({14,4,37,28},156/1),({14,4,36,29},151/1),({14,4,35,30},131/1),({14,4,34,31},97/1),({14,4,33,32},52/1),({13,5,49,16},2/1),({13,5,48,17},6/1),({13,5,47,18},15/1),({13,5,46,19},32/1),({13,5,45,20},59/1),({13,5,44,21},99/1),({13,5,43,22},153/1),({13,5,42,23},218/1),({13,5,41,24},292/1),({13,5,40,25},366/1),({13,5,39,26},430/1),({13,5,38,27},473/1),({13,5,37,28},486/1),({13,5,36,29},459/1),({13,5,35,30},392/1),({13,5,34,31},287/1),({13,5,33,32},151/1),({12,6,51,14},1/1),({12,6,50,15},5/1),({12,6,49,16},13/1),({12,6,48,17},30/1),({12,6,47,18},60/1),({12,6,46,19},108/1),({12,6,45,20},179/1),({12,6,44,21},276/1),({12,6,43,22},397/1),({12,6,42,23},538/1),({12,6,41,24},688/1),({12,6,40,25},831/1),({12,6,39,26},948/1),({12,6,38,27},1019/1),({12,6,37,28},1025/1),({12,6,36,29},955/1),({12,6,35,30},806/1),({12,6,34,31},584/1),({12,6,33,32},307/1),({11,7,52,13},2/1),({11,7,51,14},6/1),({11,7,50,15},15/1),({11,7,49,16},34/1),({11,7,48,17},68/1),({11,7,47,18},123/1),({11,7,46,19},206/1),({11,7,45,20},322/1),({11,7,44,21},473/1),({11,7,43,22},657/1),({11,7,42,23},862/1),({11,7,41,24},1073/1),({11,7,40,25},1269/1),({11,7,39,26},1422/1),({11,7,38,27},1504/1),({11,7,37,28},1495/1),({11,7,36,29},1380/1),({11,7,35,30},1156/1),({11,7,34,31},833/1),({11,7,33,32},436/1),({10,8,53,12},1/1),({10,8,52,13},3/1),({10,8,51,14},9/1),({10,8,50,15},21/1),({10,8,49,16},44/1),({10,8,48,17},82/1),({10,8,47,18},142/1),({10,8,46,19},229/1),({10,8,45,20},347/1),({10,8,44,21},496/1),({10,8,43,22},673/1),({10,8,42,23},868/1),({10,8,41,24},1065/1),({10,8,40,25},1242/1),({10,8,39,26},1376/1),({10,8,38,27},1444/1),({10,8,37,28},1426/1),({10,8,36,29},1308/1),({10,8,35,30},1090/1),({10,8,34,31},784/1),({10,8,33,32},410/1),({9,9,53,12},1/1),({9,9,52,13},2/1),({9,9,51,14},5/1),({9,9,50,15},11/1),({9,9,49,16},21/1),({9,9,48,17},38/1),({9,9,47,18},65/1),({9,9,46,19},102/1),({9,9,45,20},153/1),({9,9,44,21},217/1),({9,9,43,22},291/1),({9,9,42,23},372/1),({9,9,41,24},454/1),({9,9,40,25},526/1),({9,9,39,26},580/1),({9,9,38,27},606/1),({9,9,37,28},596/1),({9,9,36,29},545/1),({9,9,35,30},454/1),({9,9,34,31},325/1),({9,9,33,32},170/1)}, (10,0) => {}, (8,2) => {}, (10,1) => {({18,4,49,30},1/1),({18,4,48,31},2/1),({18,4,47,32},3/1),({18,4,46,33},5/1),({18,4,45,34},6/1),({18,4,44,35},7/1),({18,4,43,36},8/1),({18,4,42,37},7/1),({18,4,41,38},5/1),({18,4,40,39},3/1),({17,5,52,27},2/1),({17,5,51,28},5/1),({17,5,50,29},10/1),({17,5,49,30},18/1),({17,5,48,31},28/1),({17,5,47,32},39/1),({17,5,46,33},51/1),({17,5,45,34},60/1),({17,5,44,35},65/1),({17,5,43,36},65/1),({17,5,42,37},57/1),({17,5,41,38},42/1),({17,5,40,39},23/1),({16,6,55,24},1/1),({16,6,54,25},5/1),({16,6,53,26},13/1),({16,6,52,27},26/1),({16,6,51,28},48/1),({16,6,50,29},78/1),({16,6,49,30},117/1),({16,6,48,31},163/1),({16,6,47,32},210/1),({16,6,46,33},252/1),({16,6,45,34},283/1),({16,6,44,35},295/1),({16,6,43,36},282/1),({16,6,42,37},243/1),({16,6,41,38},179/1),({16,6,40,39},94/1),({15,7,57,22},2/1),({15,7,56,23},7/1),({15,7,55,24},18/1),({15,7,54,25},39/1),({15,7,53,26},74/1),({15,7,52,27},129/1),({15,7,51,28},205/1),({15,7,50,29},304/1),({15,7,49,30},421/1),({15,7,48,31},548/1),({15,7,47,32},672/1),({15,7,46,33},777/1),({15,7,45,34},843/1),({15,7,44,35},855/1),({15,7,43,36},802/1),({15,7,42,37},680/1),({15,7,41,38},494/1),({15,7,40,39},261/1),({14,8,59,20},1/1),({14,8,58,21},5/1),({14,8,57,22},14/1),({14,8,56,23},33/1),({14,8,55,24},69/1),({14,8,54,25},128/1),({14,8,53,26},219/1),({14,8,52,27},349/1),({14,8,51,28},519/1),({14,8,50,29},727/1),({14,8,49,30},965/1),({14,8,48,31},1210/1),({14,8,47,32},1439/1),({14,8,46,33},1622/1),({14,8,45,34},1724/1),({14,8,44,35},1719/1),({14,8,43,36},1593/1),({14,8,42,37},1336/1),({14,8,41,38},964/1),({14,8,40,39},507/1),({13,9,60,19},2/1),({13,9,59,20},6/1),({13,9,58,21},16/1),({13,9,57,22},37/1),({13,9,56,23},76/1),({13,9,55,24},141/1),({13,9,54,25},244/1),({13,9,53,26},392/1),({13,9,52,27},593/1),({13,9,51,28},849/1),({13,9,50,29},1152/1),({13,9,49,30},1483/1),({13,9,48,31},1821/1),({13,9,47,32},2124/1),({13,9,46,33},2353/1),({13,9,45,34},2469/1),({13,9,44,35},2438/1),({13,9,43,36},2235/1),({13,9,42,37},1865/1),({13,9,41,38},1341/1),({13,9,40,39},700/1),({12,10,61,18},1/1),({12,10,60,19},3/1),({12,10,59,20},9/1),({12,10,58,21},22/1),({12,10,57,22},47/1),({12,10,56,23},90/1),({12,10,55,24},161/1),({12,10,54,25},267/1),({12,10,53,26},417/1),({12,10,52,27},615/1),({12,10,51,28},861/1),({12,10,50,29},1146/1),({12,10,49,30},1456/1),({12,10,48,31},1761/1),({12,10,47,32},2032/1),({12,10,46,33},2231/1),({12,10,45,34},2323/1),({12,10,44,35},2277/1),({12,10,43,36},2080/1),({12,10,42,37},1727/1),({12,10,41,38},1238/1),({12,10,40,39},647/1),({11,11,61,18},1/1),({11,11,60,19},2/1),({11,11,59,20},5/1),({11,11,58,21},12/1),({11,11,57,22},23/1),({11,11,56,23},43/1),({11,11,55,24},75/1),({11,11,54,25},121/1),({11,11,53,26},185/1),({11,11,52,27},271/1),({11,11,51,28},373/1),({11,11,50,29},492/1),({11,11,49,30},621/1),({11,11,48,31},746/1),({11,11,47,32},854/1),({11,11,46,33},936/1),({11,11,45,34},968/1),({11,11,44,35},946/1),({11,11,43,36},864/1),({11,11,42,37},715/1),({11,11,41,38},510/1),({11,11,40,39},268/1)}, (12,0) => {}, (10,2) => {}, (12,1) => {({20,6,55,38},1/1),({20,6,54,39},2/1),({20,6,53,40},3/1),({20,6,52,41},4/1),({20,6,51,42},5/1),({20,6,50,43},5/1),({20,6,49,44},5/1),({20,6,48,45},4/1),({20,6,47,46},2/1),({19,7,59,34},1/1),({19,7,58,35},3/1),({19,7,57,36},7/1),({19,7,56,37},12/1),({19,7,55,38},19/1),({19,7,54,39},28/1),({19,7,53,40},37/1),({19,7,52,41},44/1),({19,7,51,42},48/1),({19,7,50,43},48/1),({19,7,49,44},43/1),({19,7,48,45},32/1),({19,7,47,46},17/1),({18,8,62,31},1/1),({18,8,61,32},3/1),({18,8,60,33},9/1),({18,8,59,34},19/1),({18,8,58,35},35/1),({18,8,57,36},58/1),({18,8,56,37},88/1),({18,8,55,38},123/1),({18,8,54,39},160/1),({18,8,53,40},193/1),({18,8,52,41},218/1),({18,8,51,42},228/1),({18,8,50,43},219/1),({18,8,49,44},189/1),({18,8,48,45},139/1),({18,8,47,46},74/1),({17,9,64,29},1/1),({17,9,63,30},5/1),({17,9,62,31},13/1),({17,9,61,32},29/1),({17,9,60,33},56/1),({17,9,59,34},98/1),({17,9,58,35},158/1),({17,9,57,36},236/1),({17,9,56,37},328/1),({17,9,55,38},430/1),({17,9,54,39},529/1),({17,9,53,40},613/1),({17,9,52,41},667/1),({17,9,51,42},679/1),({17,9,50,43},637/1),({17,9,49,44},541/1),({17,9,48,45},394/1),({17,9,47,46},207/1),({16,10,66,27},1/1),({16,10,65,28},4/1),({16,10,64,29},11/1),({16,10,63,30},26/1),({16,10,62,31},54/1),({16,10,61,32},101/1),({16,10,60,33},174/1),({16,10,59,34},277/1),({16,10,58,35},414/1),({16,10,57,36},582/1),({16,10,56,37},774/1),({16,10,55,38},973/1),({16,10,54,39},1160/1),({16,10,53,40},1308/1),({16,10,52,41},1393/1),({16,10,51,42},1391/1),({16,10,50,43},1289/1),({16,10,49,44},1082/1),({16,10,48,45},782/1),({16,10,47,46},410/1),({15,11,67,26},1/1),({15,11,66,27},4/1),({15,11,65,28},12/1),({15,11,64,29},28/1),({15,11,63,30},59/1),({15,11,62,31},111/1),({15,11,61,32},194/1),({15,11,60,33},313/1),({15,11,59,34},477/1),({15,11,58,35},684/1),({15,11,57,36},932/1),({15,11,56,37},1203/1),({15,11,55,38},1479/1),({15,11,54,39},1728/1),({15,11,53,40},1919/1),({15,11,52,41},2014/1),({15,11,51,42},1990/1),({15,11,50,43},1827/1),({15,11,49,44},1525/1),({15,11,48,45},1096/1),({15,11,47,46},573/1),({14,12,68,25},1/1),({14,12,67,26},3/1),({14,12,66,27},8/1),({14,12,65,28},18/1),({14,12,64,29},39/1),({14,12,63,30},74/1),({14,12,62,31},132/1),({14,12,61,32},218/1),({14,12,60,33},342/1),({14,12,59,34},504/1),({14,12,58,35},707/1),({14,12,57,36},941/1),({14,12,56,37},1197/1),({14,12,55,38},1449/1),({14,12,54,39},1674/1),({14,12,53,40},1837/1),({14,12,52,41},1915/1),({14,12,51,42},1878/1),({14,12,50,43},1715/1),({14,12,49,44},1424/1),({14,12,48,45},1022/1),({14,12,47,46},534/1),({13,13,67,26},1/1),({13,13,66,27},3/1),({13,13,65,28},8/1),({13,13,64,29},17/1),({13,13,63,30},33/1),({13,13,62,31},58/1),({13,13,61,32},96/1),({13,13,60,33},148/1),({13,13,59,34},218/1),({13,13,58,35},302/1),({13,13,57,36},401/1),({13,13,56,37},505/1),({13,13,55,38},610/1),({13,13,54,39},700/1),({13,13,53,40},767/1),({13,13,52,41},795/1),({13,13,51,42},779/1),({13,13,50,43},709/1),({13,13,49,44},589/1),({13,13,48,45},421/1),({13,13,47,46},220/1)}, (14,0) => {}, (12,2) => {}, (14,1) => {({22,8,58,49},1/1),({22,8,57,50},1/1),({22,8,56,51},1/1),({22,8,55,52},1/1),({21,9,63,44},1/1),({21,9,62,45},3/1),({21,9,61,46},5/1),({21,9,60,47},7/1),({21,9,59,48},9/1),({21,9,58,49},11/1),({21,9,57,50},12/1),({21,9,56,51},11/1),({21,9,55,52},8/1),({21,9,54,53},4/1),({20,10,67,40},1/1),({20,10,66,41},3/1),({20,10,65,42},6/1),({20,10,64,43},12/1),({20,10,63,44},20/1),({20,10,62,45},30/1),({20,10,61,46},42/1),({20,10,60,47},53/1),({20,10,59,48},62/1),({20,10,58,49},67/1),({20,10,57,50},66/1),({20,10,56,51},58/1),({20,10,55,52},43/1),({20,10,54,53},23/1),({19,11,70,37},1/1),({19,11,69,38},2/1),({19,11,68,39},6/1),({19,11,67,40},13/1),({19,11,66,41},25/1),({19,11,65,42},43/1),({19,11,64,43},68/1),({19,11,63,44},98/1),({19,11,62,45},134/1),({19,11,61,46},170/1),({19,11,60,47},202/1),({19,11,59,48},224/1),({19,11,58,49},232/1),({19,11,57,50},220/1),({19,11,56,51},189/1),({19,11,55,52},139/1),({19,11,54,53},73/1),({18,12,71,36},2/1),({18,12,70,37},5/1),({18,12,69,38},13/1),({18,12,68,39},27/1),({18,12,67,40},50/1),({18,12,66,41},84/1),({18,12,65,42},133/1),({18,12,64,43},193/1),({18,12,63,44},265/1),({18,12,62,45},341/1),({18,12,61,46},415/1),({18,12,60,47},476/1),({18,12,59,48},515/1),({18,12,58,49},519/1),({18,12,57,50},485/1),({18,12,56,51},410/1),({18,12,55,52},298/1),({18,12,54,53},157/1),({17,13,73,34},1/1),({17,13,72,35},3/1),({17,13,71,36},7/1),({17,13,70,37},17/1),({17,13,69,38},34/1),({17,13,68,39},62/1),({17,13,67,40},105/1),({17,13,66,41},166/1),({17,13,65,42},244/1),({17,13,64,43},342/1),({17,13,63,44},450/1),({17,13,62,45},563/1),({17,13,61,46},667/1),({17,13,60,47},750/1),({17,13,59,48},794/1),({17,13,58,49},792/1),({17,13,57,50},731/1),({17,13,56,51},613/1),({17,13,55,52},442/1),({17,13,54,53},232/1),({16,14,73,34},1/1),({16,14,72,35},4/1),({16,14,71,36},10/1),({16,14,70,37},21/1),({16,14,69,38},41/1),({16,14,68,39},72/1),({16,14,67,40},119/1),({16,14,66,41},181/1),({16,14,65,42},262/1),({16,14,64,43},357/1),({16,14,63,44},464/1),({16,14,62,45},570/1),({16,14,61,46},668/1),({16,14,60,47},740/1),({16,14,59,48},779/1),({16,14,58,49},769/1),({16,14,57,50},707/1),({16,14,56,51},589/1),({16,14,55,52},424/1),({16,14,54,53},221/1),({15,15,74,33},1/1),({15,15,73,34},1/1),({15,15,72,35},3/1),({15,15,71,36},6/1),({15,15,70,37},12/1),({15,15,69,38},21/1),({15,15,68,39},37/1),({15,15,67,40},56/1),({15,15,66,41},85/1),({15,15,65,42},119/1),({15,15,64,43},161/1),({15,15,63,44},204/1),({15,15,62,45},251/1),({15,15,61,46},289/1),({15,15,60,47},320/1),({15,15,59,48},334/1),({15,15,58,49},329/1),({15,15,57,50},300/1),({15,15,56,51},251/1),({15,15,55,52},179/1),({15,15,54,53},94/1)}, (16,0) => {}, (14,2) => {}, (16,1) => {({23,11,64,57},1/1),({23,11,63,58},1/1),({22,12,69,52},1/1),({22,12,68,53},2/1),({22,12,67,54},3/1),({22,12,66,55},4/1),({22,12,65,56},5/1),({22,12,64,57},5/1),({22,12,63,58},5/1),({22,12,62,59},4/1),({22,12,61,60},2/1),({21,13,73,48},1/1),({21,13,72,49},2/1),({21,13,71,50},4/1),({21,13,70,51},7/1),({21,13,69,52},11/1),({21,13,68,53},15/1),({21,13,67,54},20/1),({21,13,66,55},23/1),({21,13,65,56},25/1),({21,13,64,57},25/1),({21,13,63,58},22/1),({21,13,62,59},16/1),({21,13,61,60},9/1),({20,14,76,45},1/1),({20,14,75,46},2/1),({20,14,74,47},4/1),({20,14,73,48},8/1),({20,14,72,49},14/1),({20,14,71,50},21/1),({20,14,70,51},32/1),({20,14,69,52},42/1),({20,14,68,53},53/1),({20,14,67,54},64/1),({20,14,66,55},71/1),({20,14,65,56},72/1),({20,14,64,57},70/1),({20,14,63,58},59/1),({20,14,62,59},43/1),({20,14,61,60},24/1),({19,15,77,44},1/1),({19,15,76,45},2/1),({19,15,75,46},6/1),({19,15,74,47},11/1),({19,15,73,48},19/1),({19,15,72,49},31/1),({19,15,71,50},47/1),({19,15,70,51},63/1),({19,15,69,52},83/1),({19,15,68,53},101/1),({19,15,67,54},116/1),({19,15,66,55},125/1),({19,15,65,56},127/1),({19,15,64,57},117/1),({19,15,63,58},100/1),({19,15,62,59},73/1),({19,15,61,60},37/1),({18,16,78,43},1/1),({18,16,77,44},2/1),({18,16,76,45},5/1),({18,16,75,46},9/1),({18,16,74,47},17/1),({18,16,73,48},27/1),({18,16,72,49},41/1),({18,16,71,50},57/1),({18,16,70,51},77/1),({18,16,69,52},96/1),({18,16,68,53},115/1),({18,16,67,54},128/1),({18,16,66,55},137/1),({18,16,65,56},136/1),({18,16,64,57},126/1),({18,16,63,58},105/1),({18,16,62,59},76/1),({18,16,61,60},40/1),({17,17,77,44},1/1),({17,17,76,45},2/1),({17,17,75,46},4/1),({17,17,74,47},6/1),({17,17,73,48},12/1),({17,17,72,49},17/1),({17,17,71,50},25/1),({17,17,70,51},33/1),({17,17,69,52},42/1),({17,17,68,53},48/1),({17,17,67,54},57/1),({17,17,66,55},58/1),({17,17,65,56},58/1),({17,17,64,57},54/1),({17,17,63,58},45/1),({17,17,62,59},31/1),({17,17,61,60},18/1)}, (18,0) => {}, (16,2) => {}, (18,1) => {({21,17,80,55},1/1),({21,17,79,56},1/1),({21,17,78,57},2/1),({21,17,77,58},2/1),({21,17,76,59},3/1),({21,17,75,60},3/1),({21,17,74,61},3/1),({21,17,73,62},2/1),({21,17,72,63},2/1),({21,17,71,64},1/1),({21,17,70,65},1/1),({20,18,81,54},1/1),({20,18,80,55},1/1),({20,18,79,56},2/1),({20,18,78,57},2/1),({20,18,77,58},3/1),({20,18,76,59},3/1),({20,18,75,60},4/1),({20,18,74,61},3/1),({20,18,73,62},3/1),({20,18,72,63},2/1),({20,18,71,64},2/1),({20,18,70,65},1/1),({20,18,69,66},1/1),({19,19,80,55},1/1),({19,19,79,56},1/1),({19,19,78,57},2/1),({19,19,77,58},2/1),({19,19,76,59},3/1),({19,19,75,60},3/1),({19,19,74,61},3/1),({19,19,73,62},2/1),({19,19,72,63},2/1),({19,19,71,64},1/1),({19,19,70,65},1/1)}, (20,0) => {}, (18,2) => {({23,17,77,65},1/1),({23,17,75,67},1/1),({23,17,74,68},1/1),({23,17,73,69},1/1),({23,17,71,71},1/1),({22,18,79,63},1/1),({22,18,78,64},1/1),({22,18,77,65},1/1),({22,18,76,66},2/1),({22,18,75,67},2/1),({22,18,74,68},2/1),({22,18,73,69},2/1),({22,18,72,70},1/1),({21,19,79,63},1/1),({21,19,78,64},1/1),({21,19,77,65},2/1),({21,19,76,66},2/1),({21,19,75,67},3/1),({21,19,74,68},3/1),({21,19,73,69},3/1),({21,19,72,70},1/1),({21,19,71,71},1/1),({20,20,80,62},1/1),({20,20,78,64},1/1),({20,20,77,65},1/1),({20,20,76,66},1/1),({20,20,75,67},1/1),({20,20,74,68},2/1),({20,20,72,70},1/1)}, (20,1) => {({21,21,83,66},1/1)}, (20,2) => {({23,21,82,74},1/1),({23,21,81,75},1/1),({23,21,80,76},1/1)}, (22,2) => {}, (1,0) => {({2,0,8,1},1/1),({2,0,7,2},1/1)}, (1,1) => {({2,2,16,0},1/1)}, (3,0) => {}, (3,1) => {({8,0,21,9},1/1),({8,0,19,11},1/1),({8,0,18,12},1/1),({8,0,17,13},1/1),({8,0,15,15},1/1),({7,1,24,6},1/1),({7,1,23,7},1/1),({7,1,22,8},2/1),({7,1,21,9},2/1),({7,1,20,10},4/1),({7,1,19,11},3/1),({7,1,18,12},4/1),({7,1,17,13},2/1),({7,1,16,14},3/1),({6,2,27,3},1/1),({6,2,26,4},1/1),({6,2,25,5},3/1),({6,2,24,6},3/1),({6,2,23,7},6/1),({6,2,22,8},6/1),({6,2,21,9},9/1),({6,2,20,10},8/1),({6,2,19,11},10/1),({6,2,18,12},7/1),({6,2,17,13},8/1),({6,2,16,14},3/1),({6,2,15,15},3/1),({5,3,28,2},1/1),({5,3,27,3},1/1),({5,3,26,4},3/1),({5,3,25,5},3/1),({5,3,24,6},6/1),({5,3,23,7},6/1),({5,3,22,8},10/1),({5,3,21,9},9/1),({5,3,20,10},12/1),({5,3,19,11},9/1),({5,3,18,12},11/1),({5,3,17,13},6/1),({5,3,16,14},7/1),({4,4,27,3},1/1),({4,4,26,4},1/1),({4,4,25,5},3/1),({4,4,24,6},2/1),({4,4,23,7},5/1),({4,4,22,8},4/1),({4,4,21,9},7/1),({4,4,20,10},4/1),({4,4,19,11},7/1),({4,4,18,12},3/1),({4,4,17,13},6/1),({4,4,15,15},3/1)}, (3,2) => {}, (5,0) => {}, (5,1) => {({12,0,26,18},1/1),({12,0,25,19},1/1),({12,0,24,20},1/1),({11,1,31,13},1/1),({11,1,30,14},2/1),({11,1,29,15},4/1),({11,1,28,16},5/1),({11,1,27,17},7/1),({11,1,26,18},7/1),({11,1,25,19},8/1),({11,1,24,20},6/1),({11,1,23,21},5/1),({11,1,22,22},1/1),({10,2,34,10},2/1),({10,2,33,11},3/1),({10,2,32,12},7/1),({10,2,31,13},11/1),({10,2,30,14},18/1),({10,2,29,15},23/1),({10,2,28,16},31/1),({10,2,27,17},33/1),({10,2,26,18},37/1),({10,2,25,19},33/1),({10,2,24,20},29/1),({10,2,23,21},17/1),({10,2,22,22},8/1),({9,3,37,7},1/1),({9,3,36,8},2/1),({9,3,35,9},6/1),({9,3,34,10},10/1),({9,3,33,11},20/1),({9,3,32,12},30/1),({9,3,31,13},46/1),({9,3,30,14},59/1),({9,3,29,15},78/1),({9,3,28,16},88/1),({9,3,27,17},100/1),({9,3,26,18},96/1),({9,3,25,19},92/1),({9,3,24,20},70/1),({9,3,23,21},50/1),({9,3,22,22},14/1),({8,4,38,6},1/1),({8,4,37,7},3/1),({8,4,36,8},8/1),({8,4,35,9},14/1),({8,4,34,10},27/1),({8,4,33,11},42/1),({8,4,32,12},65/1),({8,4,31,13},87/1),({8,4,30,14},116/1),({8,4,29,15},137/1),({8,4,28,16},161/1),({8,4,27,17},166/1),({8,4,26,18},169/1),({8,4,25,19},147/1),({8,4,24,20},123/1),({8,4,23,21},74/1),({8,4,22,22},30/1),({7,5,39,5},1/1),({7,5,38,6},2/1),({7,5,37,7},6/1),({7,5,36,8},11/1),({7,5,35,9},22/1),({7,5,34,10},34/1),({7,5,33,11},55/1),({7,5,32,12},75/1),({7,5,31,13},105/1),({7,5,30,14},127/1),({7,5,29,15},156/1),({7,5,28,16},168/1),({7,5,27,17},183/1),({7,5,26,18},171/1),({7,5,25,19},160/1),({7,5,24,20},119/1),({7,5,23,21},84/1),({7,5,22,22},23/1),({6,6,38,6},2/1),({6,6,37,7},2/1),({6,6,36,8},6/1),({6,6,35,9},8/1),({6,6,34,10},17/1),({6,6,33,11},22/1),({6,6,32,12},36/1),({6,6,31,13},42/1),({6,6,30,14},59/1),({6,6,29,15},63/1),({6,6,28,16},78/1),({6,6,27,17},73/1),({6,6,26,18},80/1),({6,6,25,19},62/1),({6,6,24,20},58/1),({6,6,23,21},29/1),({6,6,22,22},17/1)}};
end;