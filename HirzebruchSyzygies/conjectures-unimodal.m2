uninstallPackage "HirzebruchSyzygies"
restart
installPackage "HirzebruchSyzygies"
needsPackage "BoijSoederberg"
needsPackage "Posets";



dataRange = value get "dataRange.m2"

--validData = {{{0,0},{1,1}},{{1,0},{1,1}},{{0,0},{1,2}},{{0,1},{1,2}},{{0,0},{1,3}},{{0,1},{1,3}},{{0,0},{1,4}},{{0,1},{1,4}},{{0,0},{1,5}},{{0,1},{1,5}},{{0,0},{1,6}},{{0,1},{1,6}},{{0,0},{1,7}},{{0,1},{1,7}},{{0,0},{1,8}},{{0,1},{1,8}},{{0,0},{2,2}},{{0,1},{2,2}},{{1,0},{2,
--       2}},{{0,0},{2,3}},{{0,1},{2,3}},{{1,0},{2,3}},{{0,0},{2,4}},{{0,1},{2,4}},{{0,2},{2,4}},{{0,3},{2,4}},{{1,0},{2,4}},{{1,1},{2,4}},{{1,2},{2,4}},{{1,3},{2,4}},{{0,0},{2,5}},{{0,1},{2,5}},{{0,2},{2,5}},{{0,3},{2,5}},{{0,4},{2,5}},{{1,0},{2,5}},{{1,1},{2,5}},{{1,2
--       },{2,5}},{{1,3},{2,5}},{{1,4},{2,5}},{{0,0},{2,6}},{{0,1},{2,6}},{{0,2},{2,6}},{{0,3},{2,6}},{{0,4},{2,6}},{{0,5},{2,6}},{{1,0},{2,6}},{{1,1},{2,6}},{{1,2},{2,6}},{{1,3},{2,6}},{{1,4},{2,6}},{{1,5},{2,6}},{{0,0},{2,7}},{{0,1},{2,7}},{{0,2},{2,7}},{{0,3},{2,7}},{{
--       0,4},{2,7}},{{0,5},{2,7}},{{0,6},{2,7}},{{1,0},{2,7}},{{1,1},{2,7}},{{1,2},{2,7}},{{1,3},{2,7}},{{1,4},{2,7}},{{1,5},{2,7}},{{1,6},{2,7}},{{0,0},{2,8}},{{1,4},{2,8}},{{0,7},{2,9}},{{0,8},{2,10}},{{0,0},{3,3}},{{0,1},{3,3}},{{0,2},{3,3}},{{1,1},{3,3}},{{1,2},{3,3
--       }},{{2,2},{3,3}},{{0,0},{3,4}},{{0,1},{3,4}},{{0,2},{3,4}},{{0,3},{3,4}},{{1,0},{3,4}},{{1,1},{3,4}},{{1,2},{3,4}},{{1,3},{3,4}},{{2,0},{3,4}},{{2,1},{3,4}},{{2,2},{3,4}},{{2,3},{3,4}},{{0,0},{3,5}},{{0,1},{3,5}},{{0,3},{3,5}},{{0,4},{3,5}},{{1,0},{3,5}},{{1,1},{
--       3,5}},{{1,2},{3,5}},{{1,3},{3,5}},{{1,4},{3,5}},{{2,0},{3,5}},{{2,1},{3,5}},{{2,2},{3,5}},{{2,3},{3,5}},{{2,4},{3,5}},{{1,4},{3,6}},{{2,2},{4,4}}}



--------------------------------------------------------------------
--------------------------------------------------------------------
----- INPUT: (a,q,F) where F is a method whose input is of the form
----- (a,B,D,q) and whose output are Boolean values
-----
----- OUPUT: By defualt this returns a list of all the unique values 
----- that are returned when testing the function F on all data in 
----- dataRange. If the option ShowFails=>true then the output is a 
----- list {B,D} such that F(a,B,D,q)==false. 
-----
----- DESCRIPTION: Given a function F which tests a conjecture for a
----- particular (a,B,D,q) this function will test the function against
----- all of the computed data that appears in the list dataRange. 
--------------------------------------------------------------------
--------------------------------------------------------------------
testConjecture = method(Options => {ShowFails => false});
testConjecture  (ZZ,ZZ,MethodFunction) := opts -> (a,q,F) ->(
    if opts.ShowFails == false then (
	return unique delete(,apply(dataRange,L->(
		F(a,L#0,L#1,q)
		))) 
	);
    if opts.ShowFails == true then ( 
	return unique delete(,apply(dataRange,L->(
		if F(a,L#0,L#1,q) != true then L
		)))  
	); 
    )

testConjecture  (ZZ,MethodFunction) := opts -> (a,F) ->(
    if opts.ShowFails == false then (
	return unique delete(,apply(dataRange,L->(
		F(a,L#0,L#1)
		))) 
	);
    if opts.ShowFails == true then ( 
	return unique delete(,apply(dataRange,L->(
		if F(a,L#0,L#1) != true then L
		)))  
	); 
    )



---- Check whether a list is unimodal -----
isUnimodal = method();
isUnimodal  (List) := (L1) ->(
    L := apply(L1,i->(
	(if i == infinity or i == -infinity then i
	else sub(i,RR))
	));
    signs := delete(0_RR,apply(#L-1,i->(L#(i+1)-L#(i))));
    signFlips := sum delete(,apply(toList(1..#signs-1),j->(
		if ((signs#(j-1))<=0 and (signs#(j))>0) or ((signs#(j-1))>=0 and (signs#(j))<0) then 1
		)));
    signFlips <= 1
    )

G1 = {1,2,3,4,5,4,3,2,1}
G2 = {1,2,3_QQ,4,5,5,5_RR,5,4,3,2,1}
G3 = {5,4,3,2,1,1,2,3,4,5,5,5,5}
G4 = {1,2,3,4_RR,3_QQ,4,3,2,1}


isUnimodal(G1)
isUnimodal(G2)
isUnimodal(G3)
isUnimodal(G4)	    	



---- Get the indices for the non-zero entries for the q-th row of the Betti table. 
getRow = method();
getRow  (ZZ,List,List) := (q,B,D) ->(
    H := totalBetti(0,B,D);
    sort delete(,apply(keys H, k->if k#1 == q and H#k != 0 then k))
    )


---- Tests whether all the q=0,1,2 rows of the Betti table is unimodal. 
totalBettiUnimodal = method();
totalBettiUnimodal  (List,List) := (B,D) ->(
    H := totalBetti(0,B,D);
    unique apply(3,q->(
	    Keys := getRow(q,B,D);
	    isUnimodal(apply(Keys,k->H#k))
	    ))
    )

-- tests all data
delete(,apply(dataRange,D->(
	print D;
	if totalBettiUnimodal(D#0,D#1) != {true} then D
	)))
--- One failure....
totalBettiTally(0,{2,3},{4,5})

	



---- Tests whether all the q=0,1,2 rows of the schur Betti table WITHOUT multiplicty is unimodal. 
schurBettiNoMulUnimodal = method();
schurBettiNoMulUnimodal  (List,List) := (B,D) ->(
    H := schurBetti(0,B,D);
    unique apply(3,q->(
	    Keys := getRow(q,B,D);
	    isUnimodal(apply(Keys,k->#(H#k)))
	    ))
    )

-- tests function correct-ness
schurBettiNoMulUnimodal({0,0},{2,5})
H = schurBetti(0,{0,0},{2,5})
Keys = getRow(1,{0,0},{2,5})
apply(Keys,k->#(H#k))

-- tests all data
delete(,apply(dataRange,D->(
	print D;
	if schurBettiNoMulUnimodal(D#0,D#1) != {true} then D
	)))
--- Many failures (87/193) including...{{0,0},{2,2}} {{0,0},{3,3}},{{0,0},{3,4}},{{0,0},{3,5}}
--- Example of failure...
H = schurBetti(0,{0,0},{3,3})
Keys = getRow(1,{0,0},{3,3})
apply(Keys,k->#(H#k))


---- Tests whether all the q=0,1,2 rows of the schur Betti table WITH multiplicty is unimodal. 
schurBettiMulUnimodal = method();
schurBettiMulUnimodal  (List,List) := (B,D) ->(
    H := schurBetti(0,B,D);
    unique apply(3,q->(
	    Keys := getRow(q,B,D);
	    L1 := apply(Keys,k->(
		    sum apply(H#k,v->v#1)
		    ));
	    isUnimodal L1
	    ))
    )

-- tests function correct-ness
schurBettiMulUnimodal({0,0},{2,5})
H = schurBetti(0,{0,0},{2,5})
Keys = getRow(1,{0,0},{2,5})
apply(Keys,k->(
	sum apply(H#k,v->v#1)
	))

-- tests all data
delete(,apply(dataRange,D->(
	print D;
	if schurBettiMulUnimodal(D#0,D#1) != {true} then D
	)))
--- No failures....


---- Tests whether all the q=0,1,2 rows of the MAX mulitplicity schur Betti is unimodal. 
schurBettiMaxMulUnimodal = method();
schurBettiMaxMulUnimodal  (List,List) := (B,D) ->(
    H := schurBetti(0,B,D);
    unique apply(3,q->(
	    Keys := getRow(q,B,D);
	    L1 := apply(Keys,k->(
		    max apply(H#k,v->v#1)
		    ));
	    isUnimodal L1
	    ))
    )

-- tests function correct-ness
schurBettiMaxMulUnimodal({0,0},{2,5})
H = schurBetti(0,{0,0},{2,5})
Keys = getRow(1,{0,0},{2,5})
apply(Keys,k->(
	max apply(H#k,v->v#1)
	))

-- tests all data
delete(,apply(dataRange,D->(
	print D;
	if schurBettiMaxMulUnimodal(D#0,D#1) != {true} then D
	)))
--- No failures....



multiBettiUnimodal = method();
multiBettiUnimodal (List,List) := (B,D) ->(
    H := multiBetti(0,B,D);
    apply(3,q->(
	    Keys := getRow(q,B,D);
	    multiDegs := flatten apply(Keys,k->(
		    if H#k == 1 then 1
		    else exponents(H#k)
		    ));
	    H1 := new MutableHashTable;
	    apply(Keys,k->(
		    if H#k == 1 then 1
		    else apply(listForm(H#k),j->(H1#{k,j#0} = j#1;))
		    ));
	    delete(,apply(multiDegs, e->(
		    L1 := apply(Keys, k->(
			    if isSubset({{k,(k#0)*e}},keys H1) == true then H1#{k,e}
			    else 0
			    ));
		    if isUnimodal(L1) == false then e
		    )))
	    ))
    )


testedTrue = unique {
    {{0, 0}, {1, 1}} , {{0, 0}, {1, 2}} , {{0, 0}, {1, 3}} , {{0, 0},
    {1, 4}} , {{0, 0}, {1, 5}} , {{0, 0}, {1, 6}} , {{0, 0}, {1, 7}} ,
    {{0, 0}, {1, 8}} , {{0, 0}, {2, 2}} , {{0, 0}, {2, 3}} , {{0, 0},
    {2, 4}} , {{0, 0}, {2, 5}} , {{0, 0}, {2, 6}} , {{0, 0}, {2, 7}} ,
    {{0, 0}, {2, 8}} , {{0, 0}, {3, 3}} , {{0, 0}, {3, 4}} , {{0, 0},
    {3, 5}} , {{0, 1}, {1, 2}} , {{0, 1}, {1, 3}} , {{0, 1}, {1, 4}} ,
    {{0, 1}, {1, 5}} , {{0, 1}, {1, 6}} , {{0, 1}, {1, 7}} , {{0, 1},
    {1, 8}} , {{0, 1}, {2, 2}} , {{0, 1}, {2, 3}} , {{0, 1}, {2, 4}} ,
    {{0, 1}, {2, 5}} , {{0, 1}, {2, 6}} , {{0, 1}, {2, 7}} , {{0, 1},
    {2, 8}} , {{0, 1}, {3, 3}} , {{0, 1}, {3, 4}} , {{0, 1}, {3, 5}} ,
    {{0, 1}, {4, 4}} , {{0, 2}, {2, 3}} , {{0, 2}, {2, 4}} , {{0, 2},
    {2, 5}} , {{0, 2}, {2, 6}} , {{0, 2}, {2, 7}} , {{0, 2}, {2, 8}} ,
    {{0, 2}, {3, 3}} , {{0, 2}, {3, 4}} , {{0, 2}, {3, 5}} , {{0, 2},
    {4, 4}} , {{0, 3}, {2, 4}} , {{0, 3}, {2, 5}} , {{0, 3}, {2, 6}} ,
    {{0, 3}, {2, 7}} , {{0, 3}, {2, 8}} , {{0, 3}, {3, 4}} , {{0, 3},
    {3, 5}} , {{0, 3}, {4, 4}} , {{0, 4}, {2, 5}} , {{0, 4}, {2, 6}} ,
    {{0, 4}, {2, 7}} , {{0, 4}, {2, 8}} , {{0, 4}, {3, 5}} , {{0, 4},
    {3, 7}} , {{0, 4}, {4, 6}} , {{0, 5}, {2, 6}} , {{0, 5}, {2, 7}} ,
    {{0, 5}, {2, 8}} , {{0, 5}, {3, 7}} , {{0, 5}, {3, 8}} , {{0, 5},
    {4, 6}},{{1, 0}, {1, 1}} , {{1, 0}, {2, 2}} , {{1, 0}, {2, 3}} ,
    {{1, 1}, {2, 3}} , {{1, 2}, {2, 3}} , {{1, 0}, {2, 4}} , {{1, 1},
    {2, 4}} , {{1, 2}, {2, 4}} , {{1, 3}, {2, 4}} , {{1, 0}, {2, 5}} ,
    {{1, 1}, {2, 5}} , {{1, 2}, {2, 5}} , {{1, 3}, {2, 5}} , {{1, 4},
    {2, 5}} , {{1, 0}, {2, 6}} , {{1, 1}, {2, 6}} , {{1, 2}, {2, 6}} ,
    {{1, 3}, {2, 6}} , {{1, 4}, {2, 6}} , {{1, 5}, {2, 6}} , {{0, 6},
    {2, 7}} , {{1, 0}, {2, 7}} , {{1, 1}, {2, 7}} , {{1, 2}, {2, 7}} ,
    {{1, 3}, {2, 7}} , {{1, 4}, {2, 7}} , {{1, 5}, {2, 7}} , {{1, 6},
    {2, 7}} , {{0, 6}, {2, 8}} , {{0, 7}, {2, 8}} , {{1, 0}, {2, 8}} ,
    {{1, 1}, {2, 8}} , {{1, 2}, {2, 8}} , {{1, 3}, {2, 8}} , {{1, 4},
    {2, 8}} , {{0, 7}, {2, 9}} , {{1, 0}, {2, 9}} , {{1, 1}, {2, 9}} ,
    {{1, 2}, {2, 9}} , {{1, 3}, {2, 9}} , {{0, 8}, {2, 10}} , {{1, 0},
    {2, 10}} , {{1, 1}, {2, 10}} , {{1, 2}, {2, 10}} , {{1, 3}, {2, 10}}
    , {{1, 0}, {2, 11}},{{1, 1}, {3, 3}} , {{1, 2}, {3, 3}} , {{2, 2},
    {3, 3}} , {{1, 0}, {3, 4}} , {{1, 1}, {3, 4}} , {{1, 2}, {3, 4}} ,
    {{1, 3}, {3, 4}} , {{2, 0}, {3, 4}} , {{2, 1}, {3, 4}} , {{2, 2},
    {3, 4}} , {{2, 3}, {3, 4}} , {{1, 0}, {3, 5}} , {{1, 1}, {3, 5}} ,
    {{1, 2}, {3, 5}} , {{1, 3}, {3, 5}}, {{1, 4}, {3, 5}}, 
    {{2, 0}, {3, 5}}, {{2, 1}, {3, 5}}, {{2, 2}, {3, 5}}, {{2, 3}, {3, 5}},
    {{2, 4}, {3, 5}}, {{1, 0}, {4, 4}}, {{1, 1}, {4, 4}}, {{1, 2}, {4, 4}},
    {{1, 3}, {4, 4}}, {{2, 0}, {4, 4}}, {{2, 1}, {4, 4}}, {{2, 2}, {4, 4}},
    {{3, 0}, {4, 4}}, {{3, 1}, {4, 4}}, {{1, 0}, {3, 6}}, {{1, 1}, {3, 6}},
    {{1, 2}, {3, 6}}, {{1, 3}, {3, 6}}, {{1, 4}, {3, 6}}, {{2, 0}, {3, 6}},
    {{2, 0}, {3, 6}}, {{2, 1}, {3, 6}}, {{2, 2}, {3, 6}}, {{1, 2}, {4, 5}},
    {{1, 3}, {4, 5}}, {{2, 0}, {4, 5}}, {{2, 1}, {4, 5}}, {{2, 2}, {4, 5}},
    {{2, 3}, {4, 5}}, {{3, 0}, {4, 5}}, {{3, 1}, {4, 5}}, {{0, 6}, {3, 7}}
};

testNotTerminating = {{{0, 5}, {4, 7}}}

dataRange2 = toList(set dataRange - set testedTrue - set testNotTerminating)

dataRange3 = apply(sort apply(dataRange2,D->({((D#1)#0)*((D#1)#1),D#0,D#1})),D->({D#1,D#2}));

i = 0;
delete(,apply(dataRange3,D->(
	print i;
	i = i+1;
	print D;
	time L1 := multiBettiUnimodal(D#0,D#1);
	print (L1 == {{},{},{}});
	if L1 != {{},{},{}} then {D,L1}
	)))

