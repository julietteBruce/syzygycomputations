needsPackage "Polyhedra"

--------------------------------------------------------------------
--------------------------------------------------------------------
----- INPUT: (q,D) 
-----
----- OUPUT: an integer
-----
----- DESCRIPTION: This function computes h^0(F_a, O(D)).
-----
----- WARNINGS: This code can probably be made faster if we 
----- actually write down the formula's for number of lattice points
--------------------------------------------------------------------
--------------------------------------------------------------------
h0 = (a,D)->(
    M := matrix {{0, 0, D#1, (a+1)*D#1},{0, -D#0, -D#0, 0}};
    P := convexHull(M);
    #latticePoints(P)
    )


--------------------------------------------------------------------
--------------------------------------------------------------------
----- INPUT: (q,D) 
-----
----- OUPUT: an integer
-----
----- DESCRIPTION: This function computes the number of boundary
----- lattice points for the polytope associated to F_a and O(D)
----- where F_a is the a-th Hirzebruch surface
-----
----- WARNINGS: This code can probably be made faster if we 
----- actually write down the formula's for number of lattice points
--------------------------------------------------------------------
--------------------------------------------------------------------
boundary = (a,D)->(
    M := matrix {{0, 0, D#1, (a+1)*D#1},{0, -D#0, -D#0, 0}};
    P := convexHull(M);
    #latticePoints(P) - #interiorLatticePoints(P)
    )

--------------------------------------------------------------------
--------------------------------------------------------------------
----- INPUT: (a,D) 
-----
----- OUPUT: an integer
-----
----- DESCRIPTION: This function computes the lattice width of the
----- polytope associated to F_a and O(D) where F_a is the a-th 
----- Hirzebruch surface
-----
----- WARNINGS: This code makes a few assumptions about the lattice
----- width of a polytope with no internal lattice points...
--------------------------------------------------------------------
--------------------------------------------------------------------

--- This uses the recursive relationship described below Definition
--- 1.5 of CCDL
latticeWidth = (a,D)->(
    M1 = matrix {{0, 0, D#1, (a+1)*D#1},{0, -D#0, -D#0, 0}};
    P = convexHull(M1);
    i = 0;
    while dim(P) > 1 do (
	L1 = interiorLatticePoints(P);
       	if L1 == {} then break (
	    M2 = vertices(P);
	    -- converts vertices to column vectors
	    L2 = apply(rank source M2,i->M2_{i});
	    -- finds pairwise differences between vertices
	    L3 = unique flatten apply(L2,i->(apply(L2,j->(i-j))));
	    -- computes the horizontal distance betewen horizontalpoints 
    	    L4 = delete(,flatten apply(L3,i->(
		    if i^{1} == 0 then flatten entries i^{0}
		    )));
	    LT = max L4;
	    )
	else (
	    N = matrix{{},{}};
	    M3 = last apply(L1,i->(N = N|i));
	    P = convexHull(M3);
	    i = i+1;
	    )
	);
    if dim(P) > 1 then sub(LT+2*i,ZZ)
    else 2*i
    )


-------------------------------------------------------------------
--------------------------------------------------------------------
----- INPUT: (q,B,D) 
-----
----- OUPUT: a boolean value
-----
----- DESCRIPTION: This function checks whether the pushforward of
----- O(B) by O(D) on F_{a} is Cohen-Macaulay. Here F_a is the a-th
----- Hirzebruch surface.
-----
----- WARNINGS: Everything is hardcoded to only handle a = 0.
--------------------------------------------------------------------
--------------------------------------------------------------------
isCM = (a,B,D)->(
    if a != 0 then return "error: a != 0";
    if ((D#0)/(D#1)*(B#1)-B#0 < 2) and ((D#1)/(D#0)*(B#0)-B#1 < 2) then true
    else false
    )


--------------------------------------------------------------------
--------------------------------------------------------------------
----- INPUT: (q,a,D,B) 
-----
----- OUPUT: an integer
-----
----- DESCRIPTION: This function gives the largest N such that we
----- are certain that the Koszul cohomology groups 
----- K_{p,q}(F_a, O(B), O(D))=0 for all p<N. Here F_a is the a-th
----- Hirzebruch surface. Note this lower bound is not necesarily sharp. 
-----
----- CAVEATS: This function will at time assume certain conditions
----- on the relationship between B and D in order to apply results 
----- from the literature. 
-----
----- WARNINGS: Everything is hardcoded to only handle a = 0.
--------------------------------------------------------------------
--------------------------------------------------------------------




lowerBound = (q,a,B,D)->(
--    if isCM(a,B,D) == false then return "error not CM";
    --This is crucial for the use of duality esp in the a = 0 and B != (0,0) cases
    if B == {0,0} then (
	if q == 0 then (0)
	else if q == 1 then (1)
	-- for q = 2 see [CCDL, Thm 1.4]
	else if q == 2 then (
	    boundary(a,D)-2
	    )
	)
    else (
	if a != 0 then return "error: a != 0 and B not zero";
    	-- for q = 0 see Prop 5.1 in EL
    	if q == 0 then (0) 
    	else if q == 1 then min{B#1+1,B#0+1}
    	-- for q = 2 apply duality (see Prop 3.5 in EL)
    	else if q == 2 then ((D_0+1)*(D_1+1)-3-(D_0-B_0-1)*(D_1-B_1-1)+1)
	)
    )

--------------------------------------------------------------------
--------------------------------------------------------------------
----- INPUT: (q,a,D,B) 
-----
----- OUPUT: an integer
-----
----- DESCRIPTION: This function gives the smallest N such that we
----- are certain that the Koszul cohomology groups 
----- K_{p,q}(F_a, O(B), O(D))=0 for all p>N. Here F_a is the a-th
----- Hirzebruch surface. Note this lower bound is not necesarily sharp. 
-----
----- CAVEATS: This function will at time assume certain conditions
----- on the relationship between B and D in order to apply results 
----- from the literature. 
-----  it also assumes a conjecture of Juliette's...
-----
----- WARNINGS: Everything is hardcoded to only handle a = 0.
--------------------------------------------------------------------
--------------------------------------------------------------------
upperBound = (q,a,B,D)->(
--    if isCM(a,B,D) == false then return "error not CM";
    if a == 0 then(
	pDim := (D#0+1)*(D#1+1)-3;
    	dualB := D-B-{2,2};
    	pDim - lowerBound(2-q,a,dualB,D)
    	)
    else(    
     if B ! = {0,0} then return "error: uppperBound code currently assumes  a = 0 or B = {0,0}";
     if q == 0 then (0)
     else if q == 1 then (
	 h0(a,D)-(latticeWidth(a,D)+1)
	 )
     else if q == 2 then (
	 h0(a,D)-3
	 )
     ) 
    )


--------------------------------------------------------------------
--------------------------------------------------------------------
----- INPUT: (a,B,D) 
-----
----- OUPUT: a list of p,q entries where there could be overlaps.
-----
--------------------------------------------------------------------
--------------------------------------------------------------------

relevantRange = (a,B,D) ->(
    L = apply(3,q->(toList apply((lowerBound(q,a,B,D)..upperBound(q,a,B,D)),i->(i))));
    R1 = delete(,apply(L#0, i->(if member(i-1,L#1)==true or member(i-2,L#2)==true then (i,0))));
    R2 = delete(,apply(L#1, i->(if member(i+1,L#0)==true or member(i-1,L#2)==true then (i,1))));
    R3 = delete(,apply(L#2, i->(if member(i+2,L#0)==true or member(i+1,L#1)==true then (i,2))));
    unique R1|R2|R3
    )


--------------------------------------------------------------------
--------------------------------------------------------------------
----- INPUT: (p,q,a,B,D) 
-----
----- OUPUT: (p',q',a,B',D) where K_{p',q'}(a,B'; D) is isomorphic to K_{p,q}(a,B;D).
-----        This will basically tell you the choice you have
-----        when computing these groups.
-----
--------------------------------------------------------------------
--------------------------------------------------------------------


koszulDual = (p,q,a,B,D) ->(
    pDim := (D#0+1)*(D#1+1)-3;
    dualB := D-B-{2,2};
    (pDim - p, 2 - 1, a dualB,D)
    )


end;

relevantRange(0,{1,2},{3,3})
relevantRange(0,{0,0},{3,3})





q
D
latticeWidth(a,D)
a = 0
D = {2,4}
B = {0,2}
lowerBound(0,a,B,D)
lowerBound(1,a,B,D)
lowerBound(2,a,B,D)
upperBound(0,a,B,D)
upperBound(1,a,B,D)
upperBound(2,a,B,D)



a = 0
D = {2,3}
B = {0,0}
lowerBound(0,a,B,D)
lowerBound(1,a,B,D)
lowerBound(2,a,B,D)

a = 0
D = {2,4}
B = {0,0}
lowerBound(0,a,B,D)
lowerBound(1,a,B,D)
lowerBound(2,a,B,D)

a = 0
D = {3,3}
B = {0,0}
lowerBound(0,a,B,D)
lowerBound(1,a,B,D)
lowerBound(2,a,B,D)


a = 0
D = {4,5}
B = {2,3}
upperBound(0,a,B,D)
upperBound(1,a,B,D)
upperBound(2,a,B,D)

a = 0
D = {2,3}
B = {0,0}
upperBound(0,a,B,D)
upperBound(1,a,B,D)
upperBound(2,a,B,D)

a = 0
D = {2,4}
B = {0,0}
upperBound(0,a,B,D)
upperBound(1,a,B,D)
upperBound(2,a,B,D)

a = 0
D = {3,3}
B = {0,0}
upperBound(0,a,B,D)
upperBound(1,a,B,D) --- This email doesnt seem right...
upperBound(2,a,B,D)




--load "relevantRange.m2"
B = {0,0}
D = {2,4}
relevantRange(0,B,D)

B = {0,0}
D = {2,3}
relevantRange(0,B,D)

B = {0,0}
D = {3,3}
relevantRange(0,B,D)
