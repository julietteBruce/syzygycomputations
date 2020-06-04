uninstallPackage "HirzebruchSyzygies"
restart
installPackage "HirzebruchSyzygies"
needsPackage "BoijSoederberg"

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


--------------------------------------------------------------------
--------------------------------------------------------------------
----- INPUT: (a,B,D,q)
-----
----- OUPUT: By default a sequence (p,q) such that K_{p,q}(a,B,D)
----- is non-zero, but K_{p+1,q}(a,B,D) is equal to zero. If
----- the option Value=>true is used the value of K_{p+1,q}(a,B,D)
----- is also returned 
-----
----- DESCRIPTION: This function is finds the last non-zero entry
----- in the q-th row of the Betti table for (a,B,D). This is done 
----- using our computational data, i.e. totalBetti(a,B,D) is loaded.
-----
----- CAVEAT: Note that if the request row is empty then the index
----- returned is -infinity.
--------------------------------------------------------------------
--------------------------------------------------------------------
lastNonzeroEntry = method(Options => {Value => false});
lastNonzeroEntry  (ZZ,List,List,ZZ) := opts -> (a,B,D,q) ->(
    H =  totalBetti(a,B,D);
    lastKey = max delete(,apply(keys H, k->(
	    if k#1 == q and H#k !=0 then k
	    ))); 
    if opts.Value == false then ( 
	return lastKey 
	);
    if opts.Value == true then ( 
	return {lastKey,H#lastKey} 
	); 
    )

--------------------------------------------------------------------
--------------------------------------------------------------------
----- INPUT: (a,B,D,q)
-----
----- OUPUT: The Schur functor decomposition of K_{p,q}(a,B,D) where
----- (p,q) is such that  K_{p,q}(a,B,D) is non-zero, but 
----- K_{p+1,q}(a,B,D) is equal to zero.
-----
----- DESCRIPTION: This function is finds the Schur functor decomposition
----- for the last non-zero entry in the q-th row of the Betti table for (a,B,D). 
-----
----- CAVEAT: If the row is empty then an error is returned. 
--------------------------------------------------------------------
--------------------------------------------------------------------
lastSchurEntry = method();
lastSchurEntry (ZZ,List,List,ZZ) := (a,B,D,q) ->(
    lastKey := lastNonzeroEntry(a,B,D,q);
    if (#lastKey) == 1 then error "row empty";
    if (#lastKey) == 2 then return (schurBetti(a,B,D))#(lastKey)
    )


--------------------------------------------------------------------
--------------------------------------------------------------------
----- INPUT: (a,B,D,q)
-----
----- OUPUT: True or False
-----
----- DESCRIPTION: This function tests whether our conjecture concerning 
----- the Schur functor decomposition for the last entry in the q=1
----- row is true for the input (a,B,D,q).
-----
----- CAVEAT: This conjecture is currently only made for a=0, 
----- B={0,0}, and q=1. If any other case is tested this function will
----- just return true. 
-----
----- NOTES: This conjecture seems to be false for the diagonal entries,
----- namely for D={2,2} and D = {3,3}.
--------------------------------------------------------------------
--------------------------------------------------------------------
lastSchurConjecture = method();
lastSchurConjecture (ZZ,List,List,ZZ) := (a,B,D,q) ->(
    if a != 0 or B != {0,0} or q != 1 then true
    else (
	A1 := binomial((D#0)+1,2)*binomial((D#1),1);
	A2 := binomial((D#0)+1,2)*binomial((D#1),1); 
	A3 := binomial((D#0)+1,1)*binomial((D#1)+1,2)-1;
	A4 := binomial((D#0)+1,1)*binomial((D#1),2)+1;
	conj := {A1,A2,A3,A4};
	{(conj,1)} == lastSchurEntry(a,B,D,q)
	)
    )

testConjecture(0,1,lastSchurConjecture,ShowFails=>true)



--------------------------------------------------------------------
--------------------------------------------------------------------
----- INPUT: (a,B,D,q)
-----
----- OUPUT: The Schur functor decomposition of K_{p,q}(a,B,D) where
----- (p,q) is such that  K_{p,q}(a,B,D) is non-zero, but 
----- K_{p+2,q}(a,B,D) is equal to zero.
-----
----- DESCRIPTION: This function is finds the Schur functor decomposition
----- for the second to last non-zero entry in the q-th row of the 
----- Betti table for (a,B,D). 
-----
----- CAVEAT: If the row has only one non-zero entry the list
----- {-infinity} is returned. If the row is empty then an error
----- is returned.  
--------------------------------------------------------------------
--------------------------------------------------------------------
secondLastSchurEntry = method();
secondLastSchurEntry (ZZ,List,List,ZZ) := (a,B,D,q) ->(
    lastKey = lastNonzeroEntry(a,B,D,q);
    if (#lastKey) == 1 then error "row empty";
    if (#lastKey) == 2 then (
	penKey = (lastKey#0-1,lastKey#1);
	if isSubset({penKey},(keys (schurBetti(0,B,D)))) == false then return {-infinity};
	if isSubset({penKey},(keys (schurBetti(0,B,D)))) == true then (
	    return (schurBetti(a,B,D))#(penKey)
	    );
	)
    )

--------------------------------------------------------------------
--------------------------------------------------------------------
----- INPUT: (a,B,D,q)
-----
----- OUPUT: True or False
-----
----- DESCRIPTION: This function tests whether our conjecture concerning 
----- the Schur functor decomposition for the last entry in the q=1
----- row is true for the input (a,B,D,q).
-----
----- CAVEAT: This conjecture is currently only made for a=0, 
----- B={0,0}, and q=1. If any other case is tested this function will
----- just return true. 
-----
----- NOTE: This conjecture is FALSE for D={2,2},{2,3},{3,3},{3,4}.
--------------------------------------------------------------------
--------------------------------------------------------------------
secondLastSchurConjecture = method();
secondLastSchurConjecture (ZZ,List,List,ZZ) := (a,B,D,q) ->(
    if a != 0 or B != {0,0} or q != 1 then true
    else (
	ans := secondLastSchurEntry(a,B,D,q);
	if ans == {-infinity} then return true;
	if ans != {-infinity} then (
	    conj = apply((D#1),i->(
		    A1 := binomial((D#0)+1,2)*binomial((D#1),1);
		    A2 := binomial((D#0)+1,2)*binomial((D#1),1)-binomial((D#0),1); 
		    A3 := binomial((D#0)+1,1)*binomial((D#1)+1,2)-(2+i);
		    A4 := binomial((D#0)+1,1)*binomial((D#1),2)-binomial((D#1),1)+(2+i);
		    ({A1,A2,A3,A4},1)));
	    return conj == ans
	    )
    )
)



testConjecture(0,1,secondLastSchurConjecture,ShowFails=>true)




----  Boij Soderberg coefficients

koszulDual = (p,q,a,B,D) ->(
    pDim := (D#0+1)*(D#1+1)-3;
    dualB := D-B-{2,2};
    (pDim - p, 2 - 1, a,dualB,D)
    )



--------------------------------------------------------------------
--------------------------------------------------------------------
----- INPUT: (a,B,D)
-----
----- OUPUT: a sequence of integers. 
-----
----- DESCRIPTION: returns the Boij-Soederberg coefficients, following 
----- the conventions in BEGY.
-----
----- CAVEAT: only works for a=0
--------------------------------------------------------------------
--------------------------------------------------------------------
BoijSoederbergCoefficients = method();
BoijSoederbergCoefficients (ZZ,List,List) := (a,B,D) ->(
    BT := totalBettiTally(a,B,D);
    degs := decomposeDegrees(BT, TableEntries => HerzogKuhl);
    return apply(degs, e->e#0);
    )



BoijSoederbergDegreeSequence = method();
BoijSoederbergDegreeSequence (ZZ,List,List) := (a,B,D) ->(
    BT := totalBettiTally(a,B,D);
    degs := decomposeDegrees(BT, TableEntries => HerzogKuhl);
    return apply(degs, e->e#1);
    )




--------------------------------------------------------------------
--------------------------------------------------------------------
----- CONJECTURE: If (a_1,...,a_k) are the  Boij-Soederberg coefficients 
----- for B={b1,b2}, D={d1,d2}, then (a_k,..., a_1) are the  Boij-Soederberg 
-----ccoefficients for B'={d1-b1-2, d2-b2-2} D'={d1,d2} (the Koszul Dual)
-----
----- INPUT: (B,D)
-----
----- OUPUT: True or False
-----
----- DESCRIPTION: This function tests whether the above conjecture 
----- is true
-----
--------------------------------------------------------------------
--------------------------------------------------------------------

BoijSoederbergCoefficientsKoszulDualConjecture = method();
BoijSoederbergCoefficientsKoszulDualConjecture (List,List) := (B,D) -> (
    BSC := BoijSoederbergCoefficients(0,B,D);
    BSCDual := BoijSoederbergCoefficients(0,D-B-{2,2},D);
    return (reverse BSCDual) == BSC
    )



--------------------------------------------------------------------
--------------------------------------------------------------------
----- CONJECTURE: Boij-Soederberg coefficients for B={0,b1}, D={2,d1}
----- consist of the following. set [n] = [0,1,...,n-1] 
----- The degree sequence is delta_j = [3(d1+1)] - {b1+1,3d1+1-j} for j = 0,..., d1-b1-2
----- and delta_j = [3(d1+1)] - {d1-j-1,2d1+b1+3} for j = d1-b1-1,...,d1-2.
----- the coefficient on delta_j is 2(3d1)! for j neq d1-b1-2 and 2(d1+2)(3d1)! for j=d1-b1-2

----- INPUT: b1,d1 integers
-----
----- OUPUT: True or False
-----
----- DESCRIPTION: This function tests whether the above conjecture 
----- is true

--------------------------------------------------------------------
--------------------------------------------------------------------




BoijSoederberg2dConjecture = method();
BoijSoederberg2dConjecture (ZZ,ZZ) := (b1,d1) -> (
    BT := totalBettiTally(0,{0,b1},{2,d1});
    degs := decomposeDegrees(BT, TableEntries => HerzogKuhl); 
    conjDegs1 := for j in (0..d1-b1-3) list (
	(2*(3*d1)!, sort(toList(set(0..3*(d1+1)-1) - set{b1+1,3*d1+1-j})))
	);
    conjDegs2 := {(2*(d1+2)*(3*d1)!, sort(toList(set(0..3*(d1+1)-1) - set{b1+1,2*d1+b1+3})))} ;
    conjDegs := conjDegs1|conjDegs2;
    if b1>0 then (
	conjDegs3 := for j in (d1-b1-1..d1-2) list (
	    (2*(3*d1)!, sort(toList(set(0..3*(d1+1)-1) - set{d1-j-1,2*d1+b1+3})))
	    );
	conjDegs = conjDegs|conjDegs3;
    );
    return conjDegs == degs    
    )





--------------------------------------------------------------------
--------------------------------------------------------------------
----- CONJECTURE: Boij-Soederberg coefficients for B={0,d1-1}, D={2,d1}
----- consist of the following. set [n] = [0,1,...,n-1] 
----- The degree sequence is [3d1 +2]-{d1-j} for j=0,...,d1-1
----- the coefficient on [3d1 +2]-{d1} is  2*(d1+1)*(3*d1)!
----- the coefficient on the remaining ones is 2*(3*d1)! 

----- INPUT: d1 integer
-----
----- OUPUT: True or False
-----
----- DESCRIPTION: This function tests whether the above conjecture 
----- is true

--------------------------------------------------------------------
--------------------------------------------------------------------

BoijSoederberg2dLastBConjecture = method();
BoijSoederberg2dLastBConjecture (ZZ) := d1 -> (
    BT := totalBettiTally(0,{0,d1-1},{2,d1});
    degs := decomposeDegrees(BT, TableEntries => HerzogKuhl); 
    conjDegs1 := {(2*(d1+1)*(3*d1)!, delete(d1, toList(0..3*d1+1)) ) };
    conjDegs2 := for j in (1..d1-1) list(
	(2*(3*d1)!, delete(d1-j, toList(0..3*d1+1)) )
	);
    conjDegs = conjDegs1|conjDegs2;
    return conjDegs == degs    
    )





