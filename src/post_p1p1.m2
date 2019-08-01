--As of now, this file has two parts.
--Part 1: code for manipulating schur functors
--Part 2: naiveBetti code
--



---------------
---------------
--Part 1: code for manipulating schur functors
---------------
---------------

needsPackage "SchurRings";
A = QQ[t_0,t_1,t_2,t_3,MonomialOrder => Lex];

--  Since we are using GL2 x GL2 representation throughout
--  we will use the phrase 2-partition for a pair of partitions
--  each of which has at most 2 parts.  In other words,
--  2-partitions are in bijection with irred reps of GL2 x GL2

--  Input:  a 2-partition represented as a list L with 4 entries
--          (hard coded!!!)
--  Output:  the hilbert series of the corresponding product of Schur functors
--          NOTE: the output ring A is assumed to be hardcoded
--                for reasons that we forgot. 
hilbSeries22 = K->(
    D := schurRing(s,2);
    pole1 := toE D_{K#0-K#1};
    phi1 := map(A,ring pole1, 
	{t_0+t_1,t_0*t_1}|apply(4,i->0));
    hilb1 := (t_0*t_1)^(K#1) * phi1(pole1);
    pole2 := toE D_{K#2-K#3};
    phi2 := map(A,ring pole2, 
	{t_2+t_3,t_2*t_3}|apply(4,i->0));
    hilb2 := (t_2*t_3)^(K#3) * phi2(pole2);
    hilb1*hilb2
    )


--  Input: a list of pairs (2-partition, multiplicity)
--  Output: the multigraded Hilbert series of the corresponding representation
--  Note:  this code was taking too long on large examples, so we created
--         some hacky multistep stuff below.
totalHilb22 = (K)->(
    if #K == 0 then return 0_A;
    sum apply(K, i-> i_1*hilbSeries22(i_0))  
    );


--  Note: totalHilb was a bottleneck step previously.  And there is code
--        to speed that up yielding ``newTotalHilb'' which could likely
--        be adapted if/when needed.

--  Input:  A multigraded hilbert series P
--  Output: P(1,1,...,1), i.e. the corresponding total Betti number
--  Caveat: doesn't seem to work w quotient rings
totalRank = P ->(
    phi := map(ZZ,ring P, apply(numgens ring P, i-> 1));
    phi(P)
    )


--  
--  Input:  A list of pairs  (2-partition, multiplicity)
--  Output: The dimension of the corresponding representation
schurRank = L ->(
    D := schurRing(QQ,s,2);
    sum apply(L,l-> l#1*(dim s_({l#0#0,l#0#1}))*(dim s_({l#0#2,l#0#3})))
    )

--  Input:  A polynomial
--  Output: the exponent vector of the leading term
expM = m->(
    if m == 0 then return {};
    (exponents m)#0
    );


--  This is a subroutine for taking a GL2 x GL2 representation
--  and removing the product of Schur functors corresponding
--  to its highest weight vector.
--  Input:  A multigraded Hilbert series P
--  Output:  The hilbert series of P - h, where
--  h is the Hilbert series of the Schur functor
--  of the highest weight vector of P.
decomposeHilb1 = P ->(
    if P == 0 then return "P is zero";
    m := expM(leadTerm P);
    lc := leadCoefficient P;
    h := lc*hilbSeries22(m);
    (P - h,m,lc)
    )



--  By iterating decomposeHilb1, this should decompose
--  a representation into Schur functors.
--  It spits out an error if we ever get a negative
--  leading coefficient in the Hilbert series.
--  This allows us to catch errors in our Hilbert series computation.
--  Input:  A multigraded Hilbert series P
--  Output:  A pair (L,P) where L is the corresponding
--        list of Schur functors and P is the error.
--       P should be 0 unless there was an error
--       our computation of multigraded numbers.
decomposeHilb = (P) -> (
    L := {};
    if P == 0 then return {{},0};
    while leadCoefficient P > 0 do(
      	  K := decomposeHilb1(P);
	  L = L|{(K#1,K#2)};
	  P = K#0;
	  );
    (L,P)
    )


--  Input:  the output of decomposeHilb.  namely a pair
--       where the first entry is a list of 2-partiions and
--       the second entry is a polynoimal
--  Output:  the output formatted in the way we want
--          for schurBetti output,
--        


------------------------
------------------------
--Part 2: naiveBetti code
------------------------
------------------------


--  Input: an integer d and a polynomial f
--  Output:  the degree d part of f
--  Caveat:  may be screwy if ring f has a funky grading.
degPart = (d,f)->(
    L = terms f;
    out = 0;
    scan(L,l-> if degree(l) == {d} then out = out+l);
    out
    )

--Checking degree
--Input: Two degrees

--Faster Multiplication
--Input: Two polynomials f,g and a degree bound of N=(N_0,N_1)
--Output: A dictionary with keys corresponding to degrees and values being the 
--    	  correct degree part in the multiplication

fasterMult = (f,g) ->(
    deg1 = (degree(f))_0;
    deg2 = (degree(g))_0;
    A = delete(0,apply(toList(0..deg1), i->{i,degPart(i,f)}));
    B = delete(0,apply(toList(0..deg2), i->{i,degPart(i,g)}));
    T = new MutableHashTable from apply(toList(0..(deg1+deg2)),i->(i=>0));
    m0 := ideal(t_0,t_1);
    m1 := ideal(t_2,t_3);
    scan(A, i-> scan(B, j-> (T#(i_0+j_0) = T#(i_0+j_0) + (i_1*j_1))));
    T
    )
    
--    scan(toList(0..deg1), i -> (scan(toList(0..deg2), j-> (
--		    if keys(T)#?(i*j) then(
--		    	T#(i*j) = T#(i*j)+(A_i*B_j);
--		    	)
--		    else T#(i*j) = A_i*B_j);
--	    	);	  
--	    );
--	);


--- Input: tuple is a 2-tuple, D is also a 2-tuple
--- Output: true if tuple is a multiple of D

isTupleMultiple = (tuple, D) -> (
    alpha = tuple_0//D_0;
    alpha*D == tuple)


--- Input: a 4-tuple multidegree as a list and the multidegree
--- Output: List of all the multidegrees that divide the given multidegree as a list of lists
--- HARDCODED FOR B= {0,0}

div = (L,D) -> (
    output := delete(,flatten apply(toList(0..L_0), i ->(
	    flatten apply(toList(0..L_1), j->(
		    flatten apply(toList(0..L_2), k->(
			    apply(toList(0..L_3),l->(
				    if isTupleMultiple({L_0+L_1-i-j,L_2+L_3-k-l},D) then {i,j,k,l})))))))));
    output)

newDiv = (L,D) ->(
    output := {};
    for i0 from 0 to L_0 do(
	for i2 from 0 to L_2 do(
	    j0 := (L_0+L_1-i0)%D_0;
	    for j1 from 0 to floor((L_1-j0)/D_0) do(
	    	if (L_0+L_1-i0-j0-D_0*j1)%D_0 != 0 then print "Oh no";
	    	kp := (L_0+L_1-i0-j0-D_0*j1)//D_0;
		lastEntry := L_2+L_3-i2-kp*D_1;
		if lastEntry >= 0 then output = output|{{i0, j0+D_0*j1,i2,lastEntry}});
	    );
	);
    output)
	
			    



--- Input: l is a 4-tuple representing a multidegree as a list, B is a polynomial
--- Output: Coefficient of AA

coeffAA = (l,Bp,D,N) -> (
    if ((l_0+l_1) > N_0 or (l_2+l_3)>N_1) then 0
    else sum apply(newDiv(l,D), i -> (
	(coefficient(t_0^(i_0)*t_1^(i_1)*t_2^(i_2)*t_3^(i_3), Bp)*(t_0^(l_0)*t_1^(l_1)*t_2^(l_2)*t_3^(l_3))))))
      
    
--- Input: Total degree
--- Output: All 4-tuple multidegrees that sum to that total degree

totalToMulti = (num) -> (
    output := flatten apply(toList(0..num), i->(
	    flatten apply(toList(0..num-i), j->(
		    apply(toList(0..num-i-j), k->(
			    {i,j,k,num-i-j-k}))))));
    output)
    

newTotalToMulti = (k,D,B) -> (
    J:={k*D_0+B_0,k*D_1+B_1};
    apply(toList({0,0}..J), E -> {E_0,J_0-E_0,E_1,J_1-E_1})
    )
    	

--  The analogue of Algorithm 3.3 for P1xP1.
--Input: B,D as a lists {B},{D}
--Output:  The multigraded Betti table if there 
--        were no ``overlaps''
naiveMultiBetti = (B,D) ->(
    S := QQ[x_0,x_1,y_0,y_1,Degrees => {{1,0},{1,0},{0,1},{0,1}}];
    --  N is supposed to be the biggest ZZ^2 degree
    --  we would need to worry about.
    --  CAVEAT: This comes from a regularity bound, but requires various unspecified assumptions on (B,D).  Needs to be examined in detail later.
    N := (D_0+1)*(D_1+1)*{D_0,D_1}+{B_0,B_1};
    Lpre0 := flatten apply(toList(0..(N_0)),i->(
       	    apply(toList(0..(N_0)-i),j->(
		    {i,j}))));
    Lpre1 := flatten apply(toList(0..(N_1)),i->(
       	    apply(toList(0..(N_1)-i),j->(
		    {i,j}))));
    L := {};	   
    --- This scan is a bottleneck
    scan(Lpre0, a->(
	    scan(Lpre1, a'->(
		    K = {(a_0+a_1-B_0),(a'_0+a'_1-B_1)};
		    alpha = K_0//D_0;
		    if alpha*D == K then L = L|{flatten {a,a'}};
		    ))));
    -- C is also bottleneck
    C := sum(L,l-> 1*t_0^(l_0)*t_1^(l_1)*t_2^(l_2)*t_3^(l_3));
    Bexps0 := flatten apply(toList(0..D_0),i->({{i,D_0-i}}));
    Bexps1 := flatten apply(toList(0..D_1),i->({{i,D_1-i}}));
    Bexps := flatten apply(Bexps0,a->apply(Bexps1,b-> flatten{a,b}));
    Bpoly := product apply(Bexps,l-> 1 - 1*t_0^(l_0)*t_1^(l_1)*t_2^(l_2)*t_3^(l_3));
    m0 := ideal(t_0,t_1);
    m1 := ideal(t_2,t_3);
    QR := A/(m0^(N_0+1)+m1^(N_1+1));
    -- A multiplication is big bottleneck
    AA := substitute(Bpoly,QR)* substitute(C,QR);
    topGuy := (N_0+N_1)//(D_0+D_1);
    --  This will have issues if not Cohen-Macaulay.
    --  Fix later ???
    numCols :=  (D_0+1)*(D_1+1)-3;
    BBkeys0 := apply(toList(0..numCols),i-> (i,0)=>0);
    BBkeys1 := apply(toList(1..numCols),i-> (i,1)=>0);
    BBkeys2 := apply(toList(2..numCols),i-> (i,2)=>0);
--    BBkeys3 := apply(toList(3..topGuy),i-> (i,3)=>0);
    emptyTable := BBkeys0|BBkeys1|BBkeys2;
    T := new MutableHashTable from emptyTable;
    --  Not exactly sure how high k should go, we put a guess in for now
    col = 0;
    row = 0;
    for k from 0 to (((degree AA)_0 - (B_0+B_1))//(D_0+D_1)) do(
	coef = (-1)^(col)*degPart((D_0+D_1)*(k)+(B_0+B_1),AA);
	--print coef;
	if coef == 0 then sg = 0_ZZ;
	if coef != 0 then sg = sub(coefficient(leadMonomial coef,coef),ZZ);
	if sg >= 0 then ((T#(col,row) = coef); 
	    (col= col+1));
	if sg < 0 then (
	    row = row +1;
	    T#(col-1,row) = -coef);	
	);	    
    new HashTable from apply(toList(keys(T)), i-> i=>T#i)
    );

--Faster Niave
-- Same inputs and outputs

naiveMultiBettiFaster = (B,D) ->(
    -- S := QQ[x_0,x_1,y_0,y_1,Degrees => {{1,0},{1,0},{0,1},{0,1}}];
    --  N is supposed to be the biggest ZZ^2 degree
    --  we would need to worry about.
    --  CAVEAT: This comes from a regularity bound, but requires various unspecified assumptions on (B,D).  Needs to be examined in detail later.
    N := (D_0+1)*(D_1+1)*{D_0,D_1}+{B_0,B_1};
---    Lpre0 := flatten apply(toList(0..(N_0)),i->(
---       	    apply(toList(0..(N_0)-i),j->(
---		    {i,j}))));
---    Lpre1 := flatten apply(toList(0..(N_1)),i->(
---       	    apply(toList(0..(N_1)-i),j->(
---		    {i,j}))));
---    L := {};	   
    --- This scan is a bottleneck
---    scan(Lpre0, a->(
---	    scan(Lpre1, a'->(
---		    K = {(a_0+a_1-B_0),(a'_0+a'_1-B_1)};
---		    alpha = K_0//D_0;
---		    if alpha*D == K then L = L|{flatten {a,a'}};
---		    ))));
    -- C is also bottleneck
---    C := sum(L,l-> 1*t_0^(l_0)*t_1^(l_1)*t_2^(l_2)*t_3^(l_3));
    Bexps0 := flatten apply(toList(0..D_0),i->({{i,D_0-i}}));
    Bexps1 := flatten apply(toList(0..D_1),i->({{i,D_1-i}}));
    Bexps := flatten apply(Bexps0,a->apply(Bexps1,b-> flatten{a,b}));
    Bpoly := product apply(Bexps,l-> 1 - 1*t_0^(l_0)*t_1^(l_1)*t_2^(l_2)*t_3^(l_3));
---    m0 := ideal(t_0,t_1);
---    m1 := ideal(t_2,t_3);
    -- A multiplication is big bottleneck
---    AA := Bpoly*C % (m0^(N_0+1)+m1^(N_1+1));
    topGuy := (N_0+N_1)//(D_0+D_1);
    --  This will have issues if not Cohen-Macaulay.
    --  Fix later ???
    numCols :=  (D_0+1)*(D_1+1)-3;
    BBkeys0 := apply(toList(0..numCols),i-> (i,0)=>0);
    BBkeys1 := apply(toList(1..numCols),i-> (i,1)=>0);
    BBkeys2 := apply(toList(2..numCols),i-> (i,2)=>0);
--    BBkeys3 := apply(toList(3..topGuy),i-> (i,3)=>0);
    emptyTable := BBkeys0|BBkeys1|BBkeys2;
    T := new MutableHashTable from emptyTable;
    --  Not exactly sure how high k should go, we put a guess in for now
    col = 0;
    row = 0;
    for k from 0 to (topGuy) do(
---	coef = (-1)^(col)*sum(totalToMulti((D_0+D_1)*(k)+(B_0+B_1)), l-> coeffAA(l,Bpoly,D,N));
    	coef = (-1)^(col)*sum(newTotalToMulti(k,D,B), l->(
		 coeffAA(l,Bpoly,D,N)));
	--print coef;
	if coef == 0 then sg = 0_ZZ;
	if coef != 0 then sg = sub(coefficient(leadMonomial coef,coef),ZZ);
	if sg >= 0 then ((T#(col,row) = coef); 
	    (col= col+1));
	if sg < 0 then (
	    row = row +1;
	    T#(col-1,row) = -coef);	
	);	    
    new HashTable from apply(toList(keys(T)), i-> i=>T#i)
    );


--  The analogue of Algorithm 3.3 for P1xP1.
--Input: B,D
--Output:  The Betti table if there were no ``overlaps'',
--    indexed with K_{p,q} conventions
naiveBetti = (B,D)->(
    K := naiveMultiBetti(B,D);
    new HashTable from apply(keys K,k-> (k => totalRank(K#k)))
    )

naiveBettiTally = (B,D)->(
    K := naiveBetti(B,D);
    new BettiTally from apply(keys K,k-> (k#0,{k#0+k#1},k#0+k#1)=> K#k)
    )


naiveSchur = (B,D)->(
    K := naiveMultiBetti(B,D);
    new HashTable from apply(keys K,k-> (k => (decomposeHilb(K#k))_0))
    )

--Part 3: Relevant Range

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
    if isCM(a,B,D) == false then return "error not CM";
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
    if isCM(a,B,D) == false then return "error not CM";
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

--- This is a function so that fixMultiBetti can take a path to the BettiSeries files as an input
--- Input: B,D as usual, startPath is the beginning of the path where you are storing the bettiSeries files
---  i = {p,q} corresponding to the K_{p,q}-group.
---    	   we always assume the files are stored in a file called step1_(D_0)_(D_1)/b_(B_0)_(B_1)/betti/bettiSeries etc.
---    	   so the startPath is just the path to where this folder is stored
--- Output: The correct path to the file containing the Betti Series for B, D fixed and i=(i_0,i_1) the betti number.

pathToBettiSeries = (startPath,B,D,i)->(
    startPath|"/step1_"|toString(D_0)|"_"|toString(D_1)|"/b_"|toString(B_0)|"_"|toString(B_1)|"/betti/bettiSeries_"|toString(i#0)|"_"|toString(i#1)|".txt"
    )


--- WARNING THIS IS HARDCODED FOR B={0,0},AND 
--- ASSUMES WE READ IN THE FIRST ROW.
--- THE PATH IS DEFINITELY WRONG
--- Input: B,D as usual, pathIn is the path to where you are storing
---    	   the bettiSeries files. See input description for pathToBettiSeries for further formatting details.
fixMultiBetti = (B,D,pathIn)->(
    --H = new MutableHashTable from naiveMultiBetti(B,D);
    H = new MutableHashTable from buildAHash(B,D);
    L1 = delete(,apply(relevantRange(0,B,D),i->(if i#1 == 1 then i)));
    scan(L1,i->(
	    f := value get (pathToBettiSeries(pathIn,B,D,i));
	    fOLD := H#i;
	    H#i = f;
	    H#(i#0-1,i#1+1) = H#(i#0-1,i#1+1) + (f-fOLD);
	    ));
    H = new HashTable from H
    )

fixedTotalBetti = (H)->(
    applyValues(H,v->if v != 0 then sub(v,{t_0=>1,t_1=>1,t_2=>1,t_3=>1}) else 0)
    )

fixedSchurBetti = (H)->(
    applyValues(H,v->if v != 0 then (decomposeHilb(v))_0 else 0)
    )

fileName = (B,D)->(
    toString(B#0)|"_"|toString(B#1)|"_"|toString(D#0)|"_"|toString(D#1)
    )


A = QQ[t_0,t_1,t_2,t_3,MonomialOrder => Lex];

--Input: k an integer (between 0 and topGuy), D and B as usual
--Output: the possible ZZ^4-degrees of total degree k*D+B.
ZZ4degs = (k,D,B)->(
    E = k*D+B;
    flatten apply(E_0+1,i->(
	    apply(E_1+1,j->(
		    {i,E_0-i,j,E_1-j}
		    ))))
    );

--Input: D,B as usual
--Output:  a hash table BpH where BpH#k is the multigraded part of Bpoly of total degree k*D+B.
buildBPolyHash = (D,B)->( 
    Bexps0 := flatten apply(toList(0..D_0),i->({{i,D_0-i}}));
    Bexps1 := flatten apply(toList(0..D_1),i->({{i,D_1-i}}));
    Bexps := flatten apply(Bexps0,a->apply(Bexps1,b-> flatten{a,b}));
    Bpoly := product apply(Bexps,l-> 1 - 1*t_0^(l_0)*t_1^(l_1)*t_2^(l_2)*t_3^(l_3));
    N := (D_0+1)*(D_1+1)*{D_0,D_1}+{B_0,B_1};
    topGuy := (N_0+N_1)//(D_0+D_1);
    BpHM = new MutableHashTable;
    scan(topGuy+1,k-> BpHM#k=0_A);
    --tally apply(terms Bpoly, m-> (degree(m)//(D_0+D_1)) )
    scan(terms Bpoly, m->(
	    BpHM#((degree m)_0//(D_0+D_1)) = BpHM#((degree m)_0//(D_0+D_1)) + m));
    new HashTable from BpHM
    )

--Input: D,B as usual
--Output:  a hashtable AH
--WATCHOUT!!!  Dan and Bobby will check soon.
buildAHash = (B,D)->(
    BpH := buildBPolyHash(D,B);
    N := (D_0+1)*(D_1+1)*{D_0,D_1}+{B_0,B_1};
    topGuy := (N_0+N_1)//(D_0+D_1);
    CH := hashTable apply(topGuy+1,k->(
	    k=>sum(ZZ4degs(k,D,B),l->(
		    t_0^(l_0)*t_1^(l_1)*t_2^(l_2)*t_3^(l_3)))));
    mm := (ideal(t_0,t_1))^(N_0+1)+ (ideal(t_2,t_3))^(N_1+1); 
    AH = hashTable apply(topGuy+1,k->(
	k => ((sum apply(k,i->(
	(BpH#i)*(CH#(k-i))))%mm))
	));
    numCols :=  (D_0+1)*(D_1+1)-3;
    BBkeys0 := apply(toList(0..numCols),i-> (i,0)=>0);
    BBkeys1 := apply(toList(1..numCols),i-> (i,1)=>0);
    BBkeys2 := apply(toList(2..numCols),i-> (i,2)=>0);
--    BBkeys3 := apply(toList(3..topGuy),i-> (i,3)=>0);
    emptyTable := BBkeys0|BBkeys1|BBkeys2;
    T := new MutableHashTable from emptyTable;
    --  Not exactly sure how high k should go, we put a guess in for now
    col = 0;
    row = 0;
    for k from 0 to (max keys AH) do(
	coef = (-1)^(col)*AH#k;
	--print coef;
	if coef == 0 then sg = 0_ZZ;
	if coef != 0 then sg = sub(coefficient(leadMonomial coef,coef),ZZ);
	if sg >= 0 then ((T#(col,row) = coef); 
	    (col= col+1));
	if sg < 0 then (
	    row = row +1;
	    T#(col-1,row) = -coef);	
	);	    
    new HashTable from apply(toList(keys(T)), i-> i=>T#i)
    );


---- H SHOULD BE THE OUTPUT OF fixMultiBetti
makeOutputFiles =  (B,D,H)->(
    --  Writing the files
    g = openOut ("M2OutputFiles/bettiF0_"|fileName(B,D)|".m2");
    --g<< "--This file computes Betti tables for P^1P^1 for d = "|toString d|" and b = "|toString b;
    --g<< endl;
    g<< "A := QQ[t_0,t_1,t_2,t_3];";
    -- need to add bit to initialize A = QQ[t_0,t_1,t_2];
    g<< endl;
    g<< "--tb stands for Total Betti numbers";
    g<< endl;
    g<< "tb"|fileName(B,D)|" = ";
    g<< toExternalString fixedTotalBetti(H);
    g<< ";";
    g<< endl; 
    g<< "--mb stands for Multigraded Betti numbers";
    g<< endl ;
    g<< "mb"|fileName(B,D)|" = ";
    g<< toExternalString H;
    g<< ";";
    g<< endl;
    g<< "--sb represents the betti numbers as sums of Schur functors";
    g<< endl ;
    g<< "sb"|fileName(B,D)|" = ";
    g<< toExternalString fixedSchurBetti(H);
    g<< ";";
    g<< endl;
    g<< "end;";    
    close g;    
    )

end;




restart


load "Making_Faster.m2"
(B,D) = ({1,0},{2,2})
time mB = naiveMultiBetti({0,0},{2,3});
mB#(2,1)
nS = naiveSchur({0,0},{2,2});
nS#(4,1)

new BettiTally from apply(keys nS,k->(
	(k#0,{k#0+k#1},k#0+k#1)=> #(nS#k)
	))

new BettiTally from apply(keys nS,k->(
	(k#0,{k#0+k#1},k#0+k#1)=>(
	    if nS#k == {} then 0 else max apply(nS#k, i-> i_1)
	    )
	))

naiveBettiTally({1,0},{2,2})
o21
nS = oo;
nS#(3,1)
-- indexed for a Betti tally
naiveBettiTally(B,D)


K = naiveMultiBetti(B,D)
sub(K#(1,1),{t_0=>1,t_1=>1,t_2=>1,t_3=>1})
nB = o173

L#(2,1)
naiveBetti(B,D)
    Betti = new BettiTally from L;
    Betti
    );


degPart(2,f)


--  The analogue of Algorithm 3.3 for P1xP1.
--Input: B,D
--Output:  The Betti table if there were no ``overlaps''
naiveBetti = (B,D) ->(
    S := QQ[x_0,x_1,y_0,y_1,Degrees => {{1,0},{1,0},{0,1},{0,1}}];
    --  (B,D) = ({1,0},{1,2})
    --  N is supposed to be the biggest ZZ^2 degree
    --  we would need to worry about.
    --  CAVEAT: This comes from a regularity bound, but requires various unspecified assumptions on (B,D).  Needs to be examined in detail later.
    N := (D_0+1)*(D_1+1)*{D_0,D_1}+{B_0,B_1};
    Lpre0 := flatten apply(toList(0..(N_0)),i->(
       	    apply(toList(0..(N_0)-i),j->(
		    {i,j}))));
    Lpre1 := flatten apply(toList(0..(N_1)),i->(
       	    apply(toList(0..(N_1)-i),j->(
		    {i,j}))));
    L := {};	   
    scan(Lpre0, a->(
	    scan(Lpre1, a'->(
		    K = {(a_0+a_1-B_0),(a'_0+a'_1-B_1)};
		    alpha = K_0//D_0;
		    if alpha*D == K then L = L|{flatten {a,a'}};
		    ))));
    R := QQ[t_0,t_1,t_2,t_3];
    C := sum(L,l-> 1*t_0^(l_0)*t_1^(l_1)*t_2^(l_2)*t_3^(l_3));
    Bexps0 := flatten apply(toList(0..D_0),i->({{i,D_0-i}}));
    Bexps1 := flatten apply(toList(0..D_1),i->({{i,D_1-i}}));
    Bexps := flatten apply(Bexps0,a->apply(Bexps1,b-> flatten{a,b}));
    Bpoly := product apply(Bexps,l-> 1 - 1*t_0^(l_0)*t_1^(l_1)*t_2^(l_2)*t_3^(l_3));
    m0 := ideal(t_0,t_1);
    m1 := ideal(t_2,t_3);
    A := Bpoly*C % (m0^(N_0+1)*m1^(N_1+1));
    --print A;
--    sub(A,{t_1=>t_0,t_3=>t_2});
    Asg := sub(A,{t_1=>t_0,t_3=>t_0,t_2=>t_0});
--    print Asg;
    topGuy := (N_0+N_1)//(D_0+D_1);
    BBkeys0 := apply(toList(0..topGuy),i-> (i,{i},i)=>0);
    BBkeys1 := apply(toList(1..topGuy),i-> (i,{i+1},i+1)=>0);
    BBkeys2 := apply(toList(2..topGuy),i-> (i,{i+2},i+2)=>0);
    BBkeys3 := apply(toList(3..topGuy),i-> (i,{i+3},i+3)=>0);
    emptyTable := BBkeys0|BBkeys1|BBkeys2|BBkeys3;
    T := new MutableHashTable from emptyTable;
    --  Not exactly sure how high k should go, we put a guess in for now
    col = 0;
    row = 0;
    for k from 0 to topGuy+3 do(
	coef = sub(((-1)^(col)*coefficient(t_0^((D_0+D_1)*(k)+(B_0+B_1)),Asg)),ZZ);
	--print coef;
	if coef >= 0 then ((T#(col,{col+row},col+row) = coef); 
	    (col= col+1));
	if coef < 0 then (
	    row = row +1;
	    --print row;
	    --print col;
	    T#(col-1,{col+row-1},col+row-1) = -coef);
	);
    L := apply(toList(keys(T)), i-> i=>T#i);
    Betti := new BettiTally from L;
    Betti
    );



end

f= t_0^4*t_2^2*t_3^2+t_0^7*t_1^3*t_2^6*t_3^4+2*t_0^7*t_1^3*t_2^5*t_3^5

g = t_0^3*t_1^2*t_3^4 + 2*t_1^2*t_2^5*t_3^2

"/Users/BobbyLaudone/Documents/University of Wisconsin/Research/P1P1/Making_Faster/step1_"|toString(D_0)|"_"|toString(D_1)|"/b_"|toString(B_0)|"_"|toString(B_1)|"/betti/bettiSeries_"|toString(i#0)|"_"|toString(i#1)|".txt"
