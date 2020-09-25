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


delete(,apply(keys H, k->(
	if k#1 == 1 then (k#0,H#k)
	)))

    g = openOut ("HirzebruchSyzygiesNew2/bettiF0_"|fileName(B,D)|".m2");
    g<< "A := QQ[t_0,t_1,t_2,t_3];";
    g<< endl;
    g<< "--tb stands for Total Betti numbers";
    g<< endl;
    
D1 = delete(,apply(dataRange, i->(if i#0 == {0,0} then i)))

g = openOut ("bettiNumbers-Q1");
apply(D1,D->(
	H := totalBetti(0,{0,0},D#1);
	L := delete(,apply(keys H, k->(
		    --if k#1 == 1 then (k#0,sub(H#k,ZZ))
		    if k#1 == 2 then (k#0,sub(H#k,ZZ))
		    )));
    	g<< "D="|toString((D#1));
    	g<< endl;
    	g<< toExternalString L;
    	g<< endl;
	g<< endl;
	))

close g
to--------------------------------------------------------------------
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

H=schurBetti(0,{0,0},{2,4})
M1 = new MutableHashTable;
L1 = delete(,apply(keys H, h-> if h#1==2 then (if H#h != {} then M1#h = H#h)
	    ));

H=schurBetti(0,{0,0},{2,5});
M2 = new MutableHashTable;
L2 = delete(,apply(keys H, h-> if h#1==2 then (if H#h != {} then M2#h = H#h)
	    ));

H=schurBetti(0,{0,0},{3,3});
M3 = new MutableHashTable;
L3 = delete(,apply(keys H, h-> if h#1==2 then (if H#h != {} then M3#h = H#h)
	    ));

(schurBetti(0,{0,0},{3,3}))#(13,2)    
(schurBetti(0,{0,0},{3,4}))#(17,2)
(schurBetti(0,{0,0},{3,5}))#(21,2)
(schurBetti(0,{0,0},{2,5}))#(15-1,2)
(schurBetti(0,{0,0},{2,6}))#(18-1,2)
(schurBetti(0,{0,0},{2,7}))#(21-1,2)

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
monomialWeightsSubsetConjecture = method();
monomialWeightsSubsetConjecture (ZZ,List,List) := (a,B,D) ->(
    if a != 0 or B != {0,0} then true
    else (
	H1 = repsWithoutMultiplicity(dominantWeightsBetti(a,B,D));
	H2 = monomialDominantWeights(a,B,D);
	{{D,B},delete(true,delete(,apply(keys H1, h->(
		if H2 === {infinity} or H1#h === {} then true
		else (
		    if isSubset(H2#h,H1#h) == false then h
		    ))
	    )
	))}
))

apply(dataRange,D->(
	print D;
	monomialWeightsSubsetConjecture(0,D#0,D#1)
	))

testConjecture(0,1,secondLastSchurConjecture,ShowFails=>true)


H = dominantWeightsBetti(0,{0,0},{2,4})
W = monomialDominantWeights(0,{0,0},{2,4})
M = new MutableHashTable
apply(keys H, h->(
	M#h = toList (set apply(H#h,i->i#0) - set W#h)
	))

totalBettiTally(0,{0,0},{2,4})

(dominantWeightsBetti(0,{0,0},{2,3}))#(9,1)
(dominantWeightsBetti(0,{0,0},{2,3}))#(8,1)
(dominantWeightsBetti(0,{0,0},{2,3}))#(7,1)
(dominantWeightsBetti(0,{0,0},{2,3}))#(6,1)
(dominantWeightsBetti(0,{0,0},{2,3}))#(5,1)
(dominantWeightsBetti(0,{0,0},{2,3}))#(4,1)
(dominantWeightsBetti(0,{0,0},{2,3}))#(3,1)
(dominantWeightsBetti(0,{0,0},{2,3}))#(2,1)
(dominantWeightsBetti(0,{0,0},{2,3}))#(1,1)


totalBettiTally(0,{0,0},{2,4})
(dominantWeightsBetti(0,{0,0},{2,4}))
(monomialDominantWeights(0,{0,0},{2,4}))

apply(1..11,i->(
	L2 = apply((dominantWeightsBetti(0,{0,0},{2,4}))#(i,1),k->k#0);
	L1 = (monomialDominantWeights(0,{0,0},{2,4}))#(i,1);
	isSubset(L1,L2)
	))

apply(10..12,i->(
	L2 = apply((dominantWeightsBetti(0,{0,0},{2,4}))#(i,2),k->k#0);
	L1 = (monomialDominantWeights(0,{0,0},{2,4}))#(i,2);
	L1==L2
	))

(dominantWeightsBetti(0,{0,0},{2,4}))#(12,2)
monomialDominantWeights
(dominantWeightsBetti(0,{0,0},{2,4}))#(11,1)
(monomialDominantWeights(0,{0,0},{2,4}))#(11,1)

(dominantWeightsBetti(0,{0,0},{2,4}))#(10,2)
(dominantWeightsBetti(0,{0,0},{2,4}))#(8,1)
(dominantWeightsBetti(0,{0,0},{2,4}))#(7,1)
(dominantWeightsBetti(0,{0,0},{2,4}))#(6,1)
(dominantWeightsBetti(0,{0,0},{2,4}))#(5,1)
(dominantWeightsBetti(0,{0,0},{2,4}))#(4,1)
(dominantWeightsBetti(0,{0,0},{2,4}))#(3,1)
(dominantWeightsBetti(0,{0,0},{2,4}))#(2,1)
(dominantWeightsBetti(0,{0,0},{2,4}))#(1,1)

totalBettiTally(0,{0,0},{2,5})

(dominantWeightsBetti(0,{0,0},{2,5}))#(14,1)
(dominantWeightsBetti(0,{0,0},{2,5}))#(13,1)
(dominantWeightsBetti(0,{0,0},{2,5}))#(12,1)
(dominantWeightsBetti(0,{0,0},{2,5}))#(11,1)
(dominantWeightsBetti(0,{0,0},{2,5}))#(10,1)
(dominantWeightsBetti(0,{0,0},{2,5}))#(9,1)
(dominantWeightsBetti(0,{0,0},{2,5}))#(8,1)
(dominantWeightsBetti(0,{0,0},{2,5}))#(7,1)
(dominantWeightsBetti(0,{0,0},{2,5}))#(6,1)
(dominantWeightsBetti(0,{0,0},{2,5}))#(5,1)


totalBettiTally(0,{0,0},{2,6})
(dominantWeightsBetti(0,{0,0},{2,6}))#(13,1)
(dominantWeightsBetti(0,{0,0},{2,6}))#(12,1)
(dominantWeightsBetti(0,{0,0},{2,6}))#(11,1)
(dominantWeightsBetti(0,{0,0},{2,6}))#(10,1)
(dominantWeightsBetti(0,{0,0},{2,6}))#(9,1)
(dominantWeightsBetti(0,{0,0},{2,6}))#(8,1)

totalBettiTally(0,{0,0},{3,3})

(dominantWeightsBetti(0,{0,0},{3,3}))#(13,2)
(dominantWeightsBetti(0,{0,0},{3,3}))#(12,2)
(dominantWeightsBetti(0,{0,0},{3,3}))#(11,2)
(dominantWeightsBetti(0,{0,0},{3,3}))#(10,2)
(dominantWeightsBetti(0,{0,0},{3,3}))#(9,1)
3
(dominantWeightsBetti(0,{0,0},{3,4}))#(17,2)
(dominantWeightsBetti(0,{0,0},{3,5}))#(21,2)

dominantWeightsBetti(0,{0,0},{3,5})
dominantWeightsBetti(0,{0,0},{2,6})


M = new MutableHashTable
apply(dataRange,L->(
	H := dominantWeightsBetti(0,(L#0),(L#1));
	M#{L#1,L#0} = delete(,apply(keys H, k->(
		if (H#k) != {} then (
		    if (max apply((H#k),i->i#1)) > 1 then k
		    ))));
	))
M2 = new MutableHashTable
apply(keys M, k->(
	B := k#1;
	D := k#0;
	if ((D#0)-(B#0))>=2 and ((D#1)-(B#1))>=2 then (
	    if M#k != {} then M2#k = M#k;
	    if M#k == {} then M2#k = false;
	    )))

M3 = new MutableHashTable

guess = (D,B)->(
    {(2*(D#1)-1-(B#1)+(B#0)*((D#1)+1),1)}
    )

apply(keys M2,k->(
	guess((k#0),(k#1)) == M2#k
	))

multiplicityDominantWeightsConjecture = method()
multiplicityDominantWeightsConjecture (ZZ,List,List) := (a,B,D) ->(
    if a != 0 then true
    else (
	H := dominantWeightsBetti(a,B,D);
	L := delete(,apply(keys H, k->(
		if (H#k) != {} then (
		    if (max apply((H#k),i->i#1)) > 1 then k
		    ))));
	return max L <= 2
	)
    )


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
----- CONJECTURE: For B={0,0}, D={d0,d1} d0<=d1 the degree sequence is:
----- (set [n] = [0,1,...,n-1]) 
----- [(d0+1)*(d1+1)]-{1,(d0+1)(d1+1)-d1-j} for j=0,...,(d0-1)(d1-2)
----- the coefficient on [3d1 +2]-{d1} is  2*(d1+1)*(3*d1)!
----- the coefficient on the remaining ones is 2*(3*d1)! 

----- INPUT: D a list of 2 integers D#0<=D#1
-----
----- OUPUT: True or False
-----
----- DESCRIPTION: This function tests whether the above conjecture 
----- is true

--------------------------------------------------------------------
--------------------------------------------------------------------

degSeqb00 = D-> for j in (0..(D#0-1)*(D#1-2)) list( sort(toList( set (0..(D#0+1)*(D#1+1)-1) - set{1,(D#0+1)*(D#1+1)-D#0-j} )))

BoijSoederbergb00DegSequenceConjecture = method();
BoijSoederbergb00DegSequenceConjecture (List) := D -> (
    d0=D#0; d1=D#1;
    BT := totalBettiTally(0,{0,0},{d0,d1});
    degs := decomposeDegrees(BT, TableEntries => HerzogKuhl); 
    return degSeqb00(D) == degs    
    )





--------------------------------------------------------------------
--------------------------------------------------------------------
----- CONJECTURE: Boij-Soederberg coefficients for B={0,b1}, D={2,d1}, b1<=d1-2
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







--------------------------------------------------------------------
--------------------------------------------------------------------
----- CONJECTURE: For B={1,b1}, D={2,d1} 0 <= b1 <= d1-2, the degree sequence is:
----- (set [n] = [0,1,...,n-1]) 
----- [3d1+2]-{2b1+2-j} for j=0,...,b1+1

----- INPUT: integers b1, d1

----- OUPUT: True or False
-----
----- DESCRIPTION: This function tests whether the above conjecture 
----- is true

--------------------------------------------------------------------
--------------------------------------------------------------------

degSeq1b12d1 = (b1,d1)-> for j in (0..b1+1) list( sort(toList( set (0..3*d1+1) - set{2*b1+2-j} )))

BoijSoederberg1b12d1DegSequenceConjecture = method();
BoijSoederberg1b12d1DegSequenceConjecture (ZZ,ZZ) := (b1,d1) -> ( 
    return degSeq1b12d1(b1,d1) == BoijSoederbergDegreeSequence(0,{1,b1},{2,d1})
    )




--------------------------------------------------------------------
--------------------------------------------------------------------
----- CONJECTURE: For B={1,b1}, D={2,d1} 0 <= b1 <= d1-2, the degree sequence is:
----- (set [n] = [0,1,...,n-1]) 
----- [3d1+2]-{2b1+2-j} for j=0,...,b1+1

----- INPUT: integers b1, d1

----- OUPUT: True or False
-----
----- DESCRIPTION: This function tests whether the above conjecture 
----- is true

--------------------------------------------------------------------
--------------------------------------------------------------------

degSeqLastb2d1 = d1-> for j in (0..d1-1) list( sort(toList( set (0..3*d1+1) - set{2*d1-j} )))

BoijSoederbergDegSeqLastb2d1Conjecture = method();
BoijSoederbergDegSeqLastb2d1Conjecture (ZZ) := d1 -> ( 
    return degSeqLastb2d1(d1) == BoijSoederbergDegreeSequence(0,{1,d1-1},{2,d1})
    )



--------------------------------------------------------------------
--------------------------------------------------------------------
----- CONJECTURE: For B={1,0}, D={2,d1} d1>=2 the BS-coefficients are
----- {2d1+2,2d1-2}(3d1)!
----- INPUT: integer d1

----- OUPUT: True or False
-----
----- DESCRIPTION: This function tests whether the above conjecture 
----- is true

--------------------------------------------------------------------
--------------------------------------------------------------------

coefficientb10d2  = d1-> {(2*d1+2)*(3*d1)!, (2*d1-2)*(3*d1)!} 

BoijSoederbergCoefficientb10d2Conjecture = method();
BoijSoederbergCoefficientb10d2Conjecture (ZZ) := d1 -> (
    return coefficientb10d2(d1) == BoijSoederbergCoefficients(0,{1,0},{2,d1})
    )

--------------------------------------------------------------------
--------------------------------------------------------------------
----- CONJECTURE: For B={1,1}, D={2,d1} d1>=3 the BS-coefficients are
----- {4d1(2d1+1)(2d1+2), 20(d1-1)d1(d1+1),8(d1-2)(d1-1)d1}(3*d1-2)! 

----- INPUT: integer d1

----- OUPUT: True or False
-----
----- DESCRIPTION: This function tests whether the above conjecture 
----- is true

--------------------------------------------------------------------
--------------------------------------------------------------------

coefficientb11d2 = d1->{(2*d1)*(2*d1+1)*(2*d1+2)*(3*d1-2)!, 20*(d1-1)*d1*(d1+1)*(3*d1-2)!,8*(d1-2)*(d1-1)*d1*(3*d1-2)! }

BoijSoederbergCoefficientb11d2Conjecture = method();
BoijSoederbergCoefficientb11d2Conjecture (ZZ) := d1 -> (
    return coefficientb11d2(d1) == BoijSoederbergCoefficients(0,{1,1},{2,d1})
    )



--------------------------------------------------------------------
--------------------------------------------------------------------
----- CONJECTURE: For B={1,2}, D={2,d1} d1>=4  the BS-coefficients are
----- {(8/3)d1(d1+1)(2d1-1)(2d1+1), (8/3)(8d1-7)(2*d1+1)(d1+1)d1, 
-----    (2/3)d1(67*d1-61)(d1+1)(d1-2),10(d1-3)(d1-2)(d1-1)d1}(3*d1-3)! 

----- INPUT: integer d1

----- OUPUT: True or False
-----
----- DESCRIPTION: This function tests whether the above conjecture 
----- is true

--------------------------------------------------------------------
--------------------------------------------------------------------

coefficientb12d2 = d1->{(8/3)*d1*(d1+1)*(2*d1-1)*(2*d1+1)*(3*d1-3)!, (8/3)*(8*d1-7)*(2*d1+1)*(d1+1)*d1*(3*d1-3)!, 
    (2/3)*d1*(67*d1-61)*(d1+1)*(d1-2)*(3*d1-3)!,10*(d1-3)*(d1-2)*(d1-1)*d1*(3*d1-3)! }

BoijSoederbergCoefficientb12d2Conjecture = method();
BoijSoederbergCoefficientb12d2Conjecture (ZZ) := d1 -> (
    return coefficientb12d2(d1) == BoijSoederbergCoefficients(0,{1,2},{2,d1})
    )



--------------------------------------------------------------------
--------------------------------------------------------------------
----- CONJECTURE: For B={2,0}, D={3,d1} d1>=3 the BS-coefficients are
----- {3(3d1+2)(d1+1), 6(d1-1)(d1+1), 9(d1-1)d1}(4d1)!


----- INPUT: integer d1

----- OUPUT: True or False
-----
----- DESCRIPTION: This function tests whether the above conjecture 
----- is true

--------------------------------------------------------------------
--------------------------------------------------------------------
coefficientb20d3 = d1 -> {3*(3*d1+2)*(d1+1)*(4*d1)!, 6*(d1-1)*(d1+1)*(4*d1)!, 9*(d1-1)*d1*(4*d1)!}

BoijSoederbergCoefficientb20d3Conjecture = method();
BoijSoederbergCoefficientb20d3Conjecture (ZZ) := d1 -> (
    return coefficientb20d3(d1) == BoijSoederbergCoefficients(0,{2,0},{3,d1})
    )


--------------------------------------------------------------------
--------------------------------------------------------------------
----- CONJECTURE: For B={3,0}, D={4,d1} d1>=4 the BS-coefficients are
----- {8(4d1+3)(2d1+1)(d1+1), 8(4d1+3)(d1-1)(d1+1), 44d1(d1-1)(d1+1), 12d1(5d1+1)(d1-1)}(5d1)!


----- INPUT: integer d1

----- OUPUT: True or False
-----
----- DESCRIPTION: This function tests whether the above conjecture 
----- is true

--------------------------------------------------------------------
--------------------------------------------------------------------
coefficientb30d4 = d1 -> {8*(4*d1+3)*(2*d1+1)*(d1+1)*(5*d1)!, 8*(4*d1+3)*(d1-1)*(d1+1)*(5*d1)!, 44*d1*(d1-1)*(d1+1)*(5*d1)!, 12*d1*(5*d1+1)*(d1-1)*(5*d1)!}

BoijSoederbergCoefficientb30d4Conjecture = method();
BoijSoederbergCoefficientb30d4Conjecture (ZZ) := d1 -> (
    return coefficientb20d3(d1) == BoijSoederbergCoefficients(0,{3,0},{4,d1})
    )

