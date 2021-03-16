uninstallPackage "HirzebruchSyzygies"
restart
installPackage "HirzebruchSyzygies"
needsPackage "BoijSoederberg"
needsPackage "Posets";



dataRange = value get "dataRange.m2"


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




--------------------------------------------------------------------
--------------------------------------------------------------------
----- INPUT: A list L
-----
----- OUPUT: True if the list L is unimodal and false if the list
----- L is not unimodal. 
-----
----- DESCRIPTION: This function tests whether a given list L is
----- is unimodal by counting the number of sign flips in the 
----- the list of first differences. 
--------------------------------------------------------------------
--------------------------------------------------------------------
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

-*
--- Tests for isUnimodal
G1 = {1,2,3,4,5,4,3,2,1}
isUnimodal(G1)
--
G2 = {1,2,3_QQ,4,5,5,5_RR,5,4,3,2,1}
isUnimodal(G2)
--
G3 = {5,4,3,2,1,1,2,3,4,5,5,5,5}
isUnimodal(G3)
--
G4 = {1,2,3,4_RR,3_QQ,4,3,2,1}
isUnimodal(G4)
*-  	



--------------------------------------------------------------------
--------------------------------------------------------------------
----- INPUT: (a,q,B,D)
-----
----- OUPUT: A list of of the indices (p,q) such that the 
----- K_{p,q}(0,B;D) != 0.    
-----
----- DESCRIPTION: This function gets the indices for the non-zero
----- entries in the q-th row of the Betti table for (a,B,D). 
--------------------------------------------------------------------
--------------------------------------------------------------------
getRow = method();
getRow  (ZZ,List,List,ZZ) := (a,B,D,q) ->(
    H := totalBetti(0,B,D);
    sort delete(,apply(keys H, k->if k#1 == q and H#k != 0 then k))
    )








--------------------------------------------------------------------
--------------------------------------------------------------------
-----
----- CONJECTURE 5.1.(1)
-----
----- INPUT: (a,B,D,q)
-----
----- OUPUT: True or False
-----
----- DESCRIPTION: This function tests whether our conjecture that
----- the standard graded Betti numbers in the q-th row form a 
----- unimodal sequence.
-----
----- CAVEAT: This conjecture is currently only made for a=0. If any 
----- other case is tested this function will just return true. 
-----
----- NOTE: There is currently one failure, D = {4,5} and B = {2,3}
----- but this seems like an error in the data. 
--------------------------------------------------------------------
--------------------------------------------------------------------
totalBettiUnimodal = method();
totalBettiUnimodal  (ZZ,List,List,ZZ) := (a,B,D,q) ->(
    if a != 0 then true
    else (
	H := totalBetti(a,B,D);
    	Keys := getRow(a,B,D,q);
    	isUnimodal(apply(Keys,k->H#k))
	)
    )

testConjecture(0,0,totalBettiUnimodal,ShowFails=>true)
--
testConjecture(0,1,totalBettiUnimodal,ShowFails=>true)
--- One failure....
totalBettiTally(0,{2,3},{4,5})
--
testConjecture(0,2,totalBettiUnimodal,ShowFails=>true)




--------------------------------------------------------------------
--------------------------------------------------------------------
-----
----- CONJECTURE 5.1.(2)
-----
----- INPUT: (a,B,D,q)
-----
----- OUPUT: True or False
-----
----- DESCRIPTION: This function tests whether our conjecture that
----- the number of distinct irreducible Schur functors counted
----- WITH MULTIPLICITY appearing in K_{p,q}(a,B;D) in the q-th row 
----- form a unimodal sequence.
-----
----- CAVEAT: This conjecture is currently only made for a=0. If any 
----- other case is tested this function will just return true. 
-----
----- NOTE: No Failures
--------------------------------------------------------------------
--------------------------------------------------------------------
schurBettiMulUnimodal = method();
schurBettiMulUnimodal  (ZZ,List,List,ZZ) := (a,B,D,q) ->(
    if a != 0 then true
    else (
	H := schurBetti(a,B,D);
    	Keys := getRow(a,B,D,q);
	L1 := apply(Keys,k->(
		    sum apply(H#k,v->v#1)
		    ));
	isUnimodal L1
	)
    )

testConjecture(0,0,schurBettiMulUnimodal,ShowFails=>true)
--
testConjecture(0,1,schurBettiMulUnimodal,ShowFails=>true)
--
testConjecture(0,2,schurBettiMulUnimodal,ShowFails=>true)
--


--------------------------------------------------------------------
--------------------------------------------------------------------
-----
----- CONJECTURE 5.1.(3)
-----
----- INPUT: (a,B,D,q)
-----
----- OUPUT: True or False
-----
----- DESCRIPTION: This function tests whether our conjecture that
----- the MAXIMUM MULTIPLICITY of a irreducible Schur functor appearing 
----- in K_{p,q}(a,B;D) in the q-th row  form a unimodal sequence.
-----
----- CAVEAT: This conjecture is currently only made for a=0. If any 
----- other case is tested this function will just return true. 
-----
----- NOTE: No Failures
--------------------------------------------------------------------
--------------------------------------------------------------------
schurBettiMaxMulUnimodal = method();
schurBettiMaxMulUnimodal  (ZZ,List,List,ZZ) := (a,B,D,q) ->(
    if a != 0 then true
    else (
	H := schurBetti(a,B,D);
    	Keys := getRow(a,B,D,q);
	L1 := apply(Keys,k->(
		    max apply(H#k,v->v#1)
		    ));
	    isUnimodal L1
	)
    )

testConjecture(0,0,schurBettiMulUnimodal,ShowFails=>true)
--
testConjecture(0,1,schurBettiMulUnimodal,ShowFails=>true)
--
testConjecture(0,2,schurBettiMulUnimodal,ShowFails=>true)
--



--------------------------------------------------------------------
--------------------------------------------------------------------
-----
----- 
-----
----- INPUT: (B,D)
-----
----- OUPUT: True or False
-----
----- DESCRIPTION: This function tests whether our conjecture that for
----- any fixed multidegree e the function i---> \beta_{i,i*e}(0,B;D)
----- is true. 
-----
----- CAVEAT: This conjecture is currently only made for a=0.
-----
----- NOTE: No Failures, but was unable to check all the data. 
--------------------------------------------------------------------
--------------------------------------------------------------------
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

--- Pairs {B,D} for which the conjecture already holds true.
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

--- Pairs {B,D} for which testing the conjecture has yet to terminate.
testNotTerminating = {{{0, 5}, {4, 7}}}

--- Data remaining to test that might terminate. 
dataRange2 = toList(set dataRange - set testedTrue - set testNotTerminating)
--- Sorted data to test to improve runtime. 
dataRange3 = apply(sort apply(dataRange2,D->({((D#1)#0)*((D#1)#1),D#0,D#1})),D->({D#1,D#2}));

--- Testing the remaining data. 
i = 0;
delete(,apply(dataRange3,D->(
	print i;
	i = i+1;
	print D;
	time L1 := multiBettiUnimodal(D#0,D#1);
	print (L1 == {{},{},{}});
	if L1 != {{},{},{}} then {D,L1}
	)))




--------------------------------------------------------------------
--------------------------------------------------------------------
-----
----- 
-----
----- INPUT: (a,B,D,q)
-----
----- OUPUT: True or False
-----
----- DESCRIPTION: This function tests whether our conjecture that
----- the number of distinct irreducible Schur functors appearing in
----- K_{p,q}(a,B;D) in the q-th row form a unimodal sequence.
-----
----- CAVEAT: This conjecture is currently only made for a=0. If any 
----- other case is tested this function will just return true. 
-----
----- NOTE: There are many many many failures of this for both 
----- q=0 and 1, but now for q=2.
--------------------------------------------------------------------
--------------------------------------------------------------------
schurBettiNoMulUnimodal = method();
schurBettiNoMulUnimodal  (ZZ,List,List,ZZ) := (a,B,D,q) ->(
    if a != 0 then true
    else (
	H := schurBetti(a,B,D);
    	Keys := getRow(a,B,D,q);
	isUnimodal(apply(Keys,k->#(H#k)))
	)
    )

F1 = testConjecture(0,0,schurBettiNoMulUnimodal,ShowFails=>true)
--
F2 = testConjecture(0,1,schurBettiNoMulUnimodal,ShowFails=>true)
--
F3 = testConjecture(0,2,schurBettiNoMulUnimodal,ShowFails=>true)
--

---- Total Number of Failures (as of 2/24/21)
#(unique F1|F2|F3)
---- Total Number of Tests (as of 2/24/21)
#dataRange

---- EXAMPLE 5.3 -- Failure of the conjecture for B={0,0} and D={3,4},{3,5}.
H = schurBetti(0,{0,0},{3,4})
Keys = getRow(0,{0,0},{3,4},1)
apply(Keys,k->#(H#k))
--
H = schurBetti(0,{0,0},{3,5})
Keys = getRow(0,{0,0},{3,5},1)
apply(Keys,k->#(H#k))









