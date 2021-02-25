uninstallPackage "HirzebruchSyzygies"
restart
installPackage "HirzebruchSyzygies"
needsPackage "BoijSoederberg"
needsPackage "Posets";



dataRange = value get "dataRange.m2"

--------------------------------------------------------------------
--------------------------------------------------------------------
----- INPUT: H, A has htable containing the Schur functor decompositions
-----
----- OUPUT: H, A hash table containing the Schur functor decompositions
----- without the multiplicities.
-----
----- DESCRIPTION: This function removes the multiplicity from the
----- the hash table representing the contianing the Schur functor
----- decompositions. 
-----
----- CAVEAT:  
--------------------------------------------------------------------
--------------------------------------------------------------------
deleteMulti = method();
deleteMulti (HashTable) := (H) ->(
    applyValues(H,v->(
	    if v == {} then ( {} )
	    else (
		apply(v,i->i#0)
		)
	    ))
    )

--------------------------------------------------------------------
--------------------------------------------------------------------
----- INPUT: (a,B,D)
-----
----- OUPUT: A list of all the redundant Schur functors appearing in
----- for B and D. 
-----
----- DESCRIPTION: A Schur functor is redundant if it appears in both
----- K_{p,q}(a,B;D) as well as K_{p-1,q+1}(a,B;D). This function finds
----- all of the redundant Schur functors for a given a, B, and D.
-----
----- CAVEAT: This currently only works for a == 0, otherwise an 
----- error is returned. 
--------------------------------------------------------------------
--------------------------------------------------------------------
redundantSchurFunctors = method();
redundantSchurFunctors (ZZ,List,List) := (a,B,D) ->(
    if a != 0 then error;
    H = schurBetti(0,B,D);
    H1 = deleteMulti(H);
    unique delete(,flatten apply(keys H1,k->(
		if isSubset({(k#0-1,k#1+1)},keys H1) == true then (
	    	    H1 = deleteMulti(H);
		    toList( (set H1#k)*(set H1#(k#0-1,k#1+1))) 
	    	    )
		)))
    )





--- This returns the list of pairs {B,D}, which 
--- have redundant Schur functors. 
L1 = delete(,apply(dataRange,D->(
	L1 := redundantSchurFunctors(0,D#0,D#1);
	if L1 != {} then D
	)))


--- Number of pairs {B,D} with redundant Schur functors.
#L1 -- 125 (as of 2/24/21)
--- Number of cases tested. 
#dataRange -- 193 (as of 2/24/21)







--- CONJECTURE 6.4
--- This finds all the D that have redundant Schur
--- functors when B = {0,0}.
---
--- For small values of D there is not redundant Schur 
--- functors, but once D is large enough there is. 
L2 = delete(,apply(dataRange,D->(
	    if D#0 == {0,0} then (
		L1 := redundantSchurFunctors(0,D#0,D#1);
		if L1 != {} then D
		)
	    )))
--- {{{0, 0}, {2, 5}}, {{0, 0}, {2, 6}}, {{0, 0}, {2, 7}}, {{0, 0}, {2, 8}}, 
--- {{0, 0}, {3, 4}}, {{0, 0}, {3, 5}}} (as of 2/24/21)

	


--- This returns a list of tuples {n,{B,D}}
--- where n is the number or redundant Schur
--- functors appearing in {B,D}. The list is
--- sorted so that the {n,{B,D}} with the 
--- largest n comes first. 

L3 = rsort delete(,apply(dataRange,D->(
	L1 := redundantSchurFunctors(0,D#0,D#1);
	{#L1,D}
	)))

--- Returns {n,{B,D}} where {B,D} has the
--- largest number of redundant Schur functors
--- and n is the number of redundant Schur functors.
max rsort delete(,apply(dataRange,D->(
	    L1 := redundantSchurFunctors(0,D#0,D#1);
	    {#L1,D}
	)))

--- {596, {{0, 8}, {2, 10}}} (as of 2/24/21)


    
    
    
--------------------------------------------------------------------
--------------------------------------------------------------------
----- INPUT: (a,B,D)
-----
----- OUPUT: A list of all the distinct Schur functors appearing in 
----- B and D considered without multiplicity
-----  
-----
----- DESCRIPTION: A Schur functor is redundant if it appears in both
----- K_{p,q}(a,B;D) as well as K_{p-1,q+1}(a,B;D). This function finds
----- all of the redundant Schur functors for a given a, B, and D.
-----
----- CAVEAT: This currently only works for a == 0, otherwise an 
----- error is returned. 
--------------------------------------------------------------------
--------------------------------------------------------------------
allSchurFunctors = method();
allSchurFunctors (ZZ,List,List) := (a,B,D) ->(
    if a != 0 then error;
    H = schurBetti(0,B,D);
    H1 = deleteMulti(H);
    unique flatten apply(keys H1, k -> H1#k)
    )
    
    
--- Finding all the distinct Schur functors
--- without multiplicity. 
H = deleteMulti(schurBetti(0,{0,8},{2,10}));
L4 = allSchurFunctors(0,{0,8},{2,10});
--- Total number of distinct Schur functors
--- without multiplicity. 
#L4 --7135 

    
    
    
--- This returns a list of tuples {n,{B,D}}
--- where n is the % of the distinct Schur
--- functors appearing in {B,D} (counted without mutlplicity),
--- which are redundant. The list is sorted so that 
--- the {n,{B,D}} with the largest n comes first. 
L5 = rsort delete(,apply(dataRange,D->(
	L1 := redundantSchurFunctors(0,D#0,D#1);
	L2 := allSchurFunctors(0,D#0,D#1);
	{#L1/#L2,D}
	)))

--- Returns {n,{B,D}} where {B,D} has the
--- largest number of redundant Schur functors
--- and n is the percentage of redundant Schur functors.
max rsort delete(,apply(dataRange,D->(
	L1 := redundantSchurFunctors(0,D#0,D#1);
	L2 := allSchurFunctors(0,D#0,D#1);
	{#L1/#L2,D}
	)))




--------------------------------------------------------------------
--------------------------------------------------------------------
----- INPUT: (a,B,D)
-----
----- OUPUT: A list of (p,q) for which every Schur functor appearing 
----- in K_{p,q}(a,B;D) is redundant
-----
----- DESCRIPTION: A Schur functor is redundant if it appears in both
----- K_{p,q}(a,B;D) as well as K_{p-1,q+1}(a,B;D). This function finds
----- all the (p,q) for which every Schur functor appearing 
----- in K_{p,q}(a,B;D) is redundant
-----
----- CAVEAT: This currently only works for a == 0, otherwise an 
----- error is returned. 
--------------------------------------------------------------------
--------------------------------------------------------------------
totallyRedundantSchurFunctors = method();
totallyRedundantSchurFunctors (ZZ,List,List) := (a,B,D) ->(
    if a != 0 then error;
    H = schurBetti(0,B,D);
    H1 = deleteMulti(H);
    unique delete(,flatten apply(keys H1,k->(
		if isSubset({(k#0-1,k#1+1)},keys H1) == true then (
	    	    H1 = deleteMulti(H);
		    if (H1#k) == (H1#(k#0-1,k#1+1)) then k
	    	    )
		)))
    )






--- This returns the list of pairs {L,{B,D}} consisting
--- of pairs {B,D} which have entries with totally
--- redundant Schur functors. 
L6 = delete(,apply(dataRange,D->(
	L1 := totallyRedundantSchurFunctors(0,D#0,D#1);
	if L1 != {} then {D,L1}
	)))

--- This is the example of totally redundant Schur
--- functors when B = {1,2} and D = {2,3}.
totallyRedundantSchurFunctors(0,{1,2},{2,3})
H = deleteMulti(schurBetti(0,{1,2},{2,3}));
H#(5,0)
H#(4,1)    
  
 
