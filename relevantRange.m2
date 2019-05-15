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

D = {2,3}
h0(0,D)
(D#0+1)*(D#1+1)

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

D = {4,4}
latticeWidth(0,D)

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
    if B == {0,0} then (
	if q == 0 then (0)
	else if q == 1 then (1)
	-- for q = 2 see [CCDL, Thm 1.4]
	else if q == 2 then (
	    boundary(a,D)-2
	    )
	)
    else (
	if a != 0 then return "error: a != 0";
    	-- for q = 0 see Prop 5.1 in EL
    	if q == 0 then (0) 
    	else if q == 1 then (1)
    	-- for q = 2 apply duality (see Prop 3.5 in EL)
    	else if q == 2 then (
	if ((B#0+D#0)< -1 and (B#0+D#0)>0) or ((B#0+D#0)>0 and (B#1+D#1)< -1) then return "error: intermediate cohomology"
	else h0(a,{-2,-2}-B+D)
	)
    ))

a = 0
D = {4,5}
B = {2,3}
lowerBound(0,a,B,D)
lowerBound(1,a,B,D)
lowerBound(2,a,B,D)

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

--------------------------------------------------------------------
--------------------------------------------------------------------
----- INPUT: (q,a,D,B) 
-----
----- OUPUT: an integer
-----
----- DESCRIPTION: This function gives the smallest N such that we
----- are certain that the Koszul cohomology groups 
----- K_{p,q}(F_a, O(B), O(D))=0 for all p>N. Here F_a is the a-th
s----- Hirzebruch surface. Note this lower bound is not necesarily sharp. 
-----
----- CAVEATS: This function will at time assume certain conditions
----- on the relationship between B and D in order to apply results 
----- from the literature. 
-----
----- WARNINGS: Everything is hardcoded to only handle a = 0.
--------------------------------------------------------------------
--------------------------------------------------------------------
upperBound = (q,a,B,D)->(
    if B == {0,0} then (
	if q == 0 then (0)
	else if q == 1 then (
	    h0(a,D)-(latticeWidth(a,D)+1)
	    )
	else if q == 2 then (
	    h0(a,D)-3
	    )
	)
    else(
	if a != 0 then return "error: a != 0";
    	-- for q = 0 see Prop 5.1 in EL
    	if q == 0 then ( 
	    if (B#0 > D#0) or (B#1 > D#1) then return "error: B>D"
	    else h0(a,B)
	    )  
    	else if q == 1 then (
	    h0(a,D)
	    )
    	-- for q = 2 apply duality (see Prop 3.5 in EL)
    	else if q == 2 then (
	    h0(a,D)
	    )
	))

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


relevantRange = (a,B,D) ->(
    L = apply(3,q->(toList apply((lowerBound(q,0,B,D)..upperBound(q,0,B,D)),i->(i))));
    R1 = delete(,apply(L#0, i->(if member(i-1,L#1)==true or member(i-2,L#2)==true then (i,0))));
    R2 = delete(,apply(L#1, i->(if member(i+1,L#0)==true or member(i-1,L#2)==true then (i,1))));
    R3 = delete(,apply(L#2, i->(if member(i+2,L#0)==true or member(i+1,L#1)==true then (i,2))));
    unique R1|R2|R3
    )

B = {0,0}
D = {2,4}
relevantRange(0,B,D)

B = {0,0}
D = {2,3}
relevantRange(0,B,D)

B = {0,0}
D = {3,3}
relevantRange(0,B,D)
