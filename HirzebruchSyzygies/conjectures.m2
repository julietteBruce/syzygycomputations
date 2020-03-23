uninstallPackage "HirzebruchSyzygies"
restart
installPackage "HirzebruchSyzygies"

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
