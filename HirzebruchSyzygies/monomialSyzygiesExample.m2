installPackage "HirzebruchSyzygies";
needsPackage "Posets";


--------------------------------------------------------------------
--------------------------------------------------------------------
----- INPUT: (a,B,D,P)
-----
----- OUPUT: A list of the monomial syzygies given by a list of 
----- multidegrees of the form {{m1.... mp},f} to represent the 
----- syzygy m1 ^.... ^mp \otimes f.  
-----
----- DESCRIPTION: This function constructs all of the non-zero monomial
----- syzygies constructed in an alagous fashion to EEL18 for the entry 
----- K_(p,q)(F_a,B;D) where the input sequence P=(p,q).
-----
----- CAVEAT: This function is currently only made for a=0, 
----- B={0,0}, and q != 0. If any other case is tested this function will
----- just return true. 
-----
--------------------------------------------------------------------
--------------------------------------------------------------------
monomialSyzygies = method();
monomialSyzygies (ZZ,List,List,Sequence) := (a,B,D,P) ->(
    p := P#0;
    q := P#1;
    if a != 0 or B != {0,0} then error;
    ---
    S = QQ[x_0,x_1,y_0,y_1,Degrees=>{{1,0},{1,0},{0,1},{0,1}}];
    I = ideal(x_0^(D#0)*y_0^(D#1), x_0^(D#0)*y_1^(D#1)+x_1^(D#0)*y_0^(D#1), x_1^(D#0)*y_1^(D#1));
    ---
    if q == 0 then return {};
    --if q == 1 then (F = {x_0^(D#0-1)*x_1^(B#0+1)*y_0^((D#1+B#1))});
    if q == 1 then (F = {x_0^(D#0+B#0)*y_0^((D#1)-1)*y_1^((B#1)+1)});
    --if q == 1 then (F = {x_0^(D#0-1)*x_1^(B#0+1)*y_0^((D#1+B#1)),x_0^(D#0+B#0)*y_0^((D#1)-1)*y_1^((B#1)+1)}) ;
    if q == 2 then (F = {x_0^(2*(D#0)-1)*x_1^((B#0)+1)*y_0^((D#1)-1)*y_1^((D#1)+(B#1)+1)});
    ---	
    flatten apply(F, f->(
	    Q := quotient(I,sub(f,S));
    	    ---
    	    B1 := unique flatten entries basis(D,S);
    	    B2 := unique flatten entries basis(q*D,S);
    	    B3 := delete(,apply(B1, m->(if m%Q == 0 then m)));
	    L0 := apply({x_0^(D#0)*y_1^(D#1), x_1^(D#0)*y_1^(D#1), x_0^(D#0)*y_0^(D#1)},i->sub(i,S));
    	    L1 := toList (set B3 - set L0);
    	    ----
    	    R := QQ[x_0,x_1,y_0,y_1,Degrees=>{{1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1}}];
    	    J := ideal(x_0^(D#0)*y_0^(D#1), x_0^(D#0)*y_1^(D#1)+x_1^(D#0)*y_0^(D#1), x_1^(D#0)*y_1^(D#1));
    	    --T = R/J;
    	    L2 := apply(L1,i->sub(i,R));
    	    B4 := apply(B1,i->sub(i,R));
    	    g := sub(f,R);
    	    B5 := apply(B2,i->(sub(i,R)));
    	    H := unique delete(,apply(B5,h->(if (g-h)%J == 0 or (g+h)%J == 0 then h)));	  
    	    ---
    	    D1 := unique delete(,flatten apply(H,g->(apply(B4,m->(if g%m == 0 then m)))));	
    	    L3 := apply({x_0^(D#0)*y_1^(D#1), x_1^(D#0)*y_1^(D#1), x_0^(D#0)*y_0^(D#1)},i->sub(i,R));
    	    D2 := toList (set D1 - set L3);
    	    L4 := toList (set L2 - set D2);
    	    apply(subsets(L4,p-(#D2)),L->(
		    {apply(L,i->degree i)|apply(D2,j->degree j),degree g}
		    ))
	    ))
	)

monomialSyzygies1 = method();
monomialSyzygies1 (ZZ,List,List,Sequence) := (a,B,D,P) ->(
    p := P#0;
    q := P#1;
    if a != 0 or B != {0,0} then error;
    ---
    S := QQ[x_0,x_1,y_0,y_1,Degrees=>{{1,0},{1,0},{0,1},{0,1}}];
    I := ideal(x_0^(D#0)*y_0^(D#1), x_0^(D#0)*y_1^(D#1)+x_1^(D#0)*y_0^(D#1), x_1^(D#0)*y_1^(D#1));
    ---
    if q == 0 then return {};
    if q == 1 then (F = {x_0^(D#0-1)*x_1^(B#0+1)*y_0^((D#1+B#1)),x_0^(D#0+B#0)*y_0^((D#1)-1)*y_1^((B#1)+1)}) ;
    if q == 2 then (F = {x_0^(2*(D#0)-1)*x_1^((B#0)+1)*y_0^((D#1)-1)*y_1^((D#1)+(B#1)+1)});
    ---	
    flatten apply(F,f->(
	    Q := quotient(I,f);
    	    ---
    	    B1 := unique flatten entries basis(D,S);
    	    B2 := delete(,apply(B1, m->(if m%Q == 0 then m)));
    	    L1 := toList (set B2 - set {x_0^(D#0)*y_1^(D#1), x_1^(D#0)*y_1^(D#1), x_0^(D#0)*y_0^(D#1)});
    	    ---
    	    R := QQ[x_0,x_1,y_0,y_1,Degrees=>{{1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1}}];
    	    L2 := apply(L1,i->sub(i,R));
    	    B3 := apply(B1,i->sub(i,R));
    	    g := sub(f,R);
    	    ---
    	    D1 := delete(,apply(B3,m->(if g%m == 0 then m)));
    	    D2 := toList (set D1 - set {x_0^(D#0)*y_1^(D#1), x_1^(D#0)*y_1^(D#1), x_0^(D#0)*y_0^(D#1)});
    	    L3 := toList (set L2 - set D2);
    	    apply(subsets(L3,p-(#D2)),L->(
		{apply(L,i->degree i)|apply(D2,j->degree j),degree g}
		))
	    ))
	)
    
--------------------------------------------------------------------
--------------------------------------------------------------------
----- INPUT: (a,B,D,P)
-----
----- OUPUT: A list of the total weights monomial syzygies 
----- 
----- DESCRIPTION: This function uses the monomialSyzygies function to
----- constructs all of the non-zero monomial syzygies for the entry 
----- K_(p,q)(F_a,B;D) where the input sequence P=(p,q), and then computes
----- the total weight of each syzygy. 
-----
----- CAVEAT: This function is currently only made for a=0, 
----- B={0,0}, and q != 0. If any other case is tested this function will
----- just return true. 
-----
--------------------------------------------------------------------
--------------------------------------------------------------------
monomialSyzygiesWeights = method();
monomialSyzygiesWeights (ZZ,List,List,Sequence) := (a,B,D,P) ->(
    L := monomialSyzygies(a,B,D,P);
    apply(L,i->(
	    sum (i#0) + (i#1)
	    ))
    )

--------------------------------------------------------------------
--------------------------------------------------------------------
----- INPUT: (L,P)
-----
----- OUPUT: A Boolean value 
-----
----- DESCRIPTION: Given lists L and P representing weights on P^1xP^1
----- (e.g. L={a0,a1,b0,b1} and P likewise) this function checks whether
----- L dominants P with respect to the product dominance order. If 
----- L dominats P true is returned, otherwise false is returned. 
-----
--------------------------------------------------------------------
--------------------------------------------------------------------
dominanceComparison = method();
dominanceComparison (List,List) := (L,P) ->(
    C1 := (L#0)>=(P#0) and ((L#0)+(L#1))>=((P#0)+(P#1));
    C2 := (L#2)>=(P#2) and ((L#2)+(L#3))>=((P#2)+(P#3));
    C1 and C2
    ) 

dominanceComparison (Sequence,Sequence) := (L,P) ->(
    dominanceComparison(L#0,P#0)
    )
 
--------------------------------------------------------------------
--------------------------------------------------------------------
----- INPUT: (L,P)
-----
----- OUPUT: A Boolean value 
----- 
----- DESCRIPTION: Given a list L representing a weight on P^1xP^1 and
----- a list of weights P this function checks whether L is dominanted
----- by any of the weights in the list P. If L is dominated by some
----- element of P false is returned, otherwise true is returned. 
-----
--------------------------------------------------------------------
--------------------------------------------------------------------
dominanceComparisonList = method();
dominanceComparisonList (List,List) := (L,P) ->(
    for p in P do if dominanceComparison(p,L) == true then return false;
    true
    ) 

dominanceComparisonList (Sequence,List) := (L,P) ->(
    for p in P do if dominanceComparison((p#0),(L#0)) == true then return false;
    true
    ) 
--------------------------------------------------------------------
--------------------------------------------------------------------
----- INPUT: (W)
-----
----- OUPUT: The domeinance poset of W.
----- 
----- DESCRIPTION: Given a list W representing weights on P1xP1 this
----- constriucts the dominance poset on W where weights are 
----- compared using the product dominance order on P1xP1. 
-----
----- CAVEAT: Note that the internal comparision functions do the 
----- opposite of what you expect. This is because of how the 
----- Poset package is implmented.
--------------------------------------------------------------------
--------------------------------------------------------------------
dominantWeightsPoset = method(Options => {Schur => false})
dominantWeightsPoset (List) := opts -> (W)->(
    if W == {} then return G;
    if opts.Schur == false then (
	dominanceComparison1 := (P,L) -> (
            C1 := (L#0)>=(P#0) and ((L#0)+(L#1))>=((P#0)+(P#1));
    	    C2 := (L#2)>=(P#2) and ((L#2)+(L#3))>=((P#2)+(P#3));
    	    C1 and C2
        );
    	return poset(W,dominanceComparison1,AntisymmetryStrategy => "none")
    );
    if opts.Schur == true then (
	dominanceComparison2 := (V,W) -> (
	    L := W#0;
	    P := V#0;    
            C1 := (L#0)>=(P#0) and ((L#0)+(L#1))>=((P#0)+(P#1));
    	    C2 := (L#2)>=(P#2) and ((L#2)+(L#3))>=((P#2)+(P#3));
    	    C1 and C2
        );
    	return poset(W,dominanceComparison2,AntisymmetryStrategy => "none")
    );
    )

--------------------------------------------------------------------
--------------------------------------------------------------------
----- INPUT: (W)
-----
----- OUPUT: A list of maximal elements in the dominance post on W.
----- 
----- DESCRIPTION: Given a list W representing weights on P1xP1 this
----- returns a list of the maximal elements in dominance poset on W.
-----
--------------------------------------------------------------------
--------------------------------------------------------------------
maximalDominantWeights = method(Options => {Schur => false})
maximalDominantWeights (List) := opts -> (L)->(
    maximalElements dominantWeightsPoset(L,Schur=>opts.Schur)
    );


--------------------------------------------------------------------
--------------------------------------------------------------------
----- INPUT: (W)
-----
----- OUPUT: A list of the dominant weights in W.
----- 
----- DESCRIPTION: Given a list W representing weights on P1xP1 this
----- returns a list of the dominant weights in W. This is similar
----- to maximalDominantWeights, but is somewhat quicker as it does not
----- require computing the entire dominance poset on W. Instead
----- a quick comparision is used to rule out some elements of W 
----- before constructing the dominance poset on this subset. 
-----
----- CAVEAT: Note we need to take unique of D1 becuase of the Poset
----- package does not consider the post {a,a} to have a max. 
--------------------------------------------------------------------
--------------------------------------------------------------------
dominantWeights = method(Options => {Schur => false})
dominantWeights (List) := opts -> (L)->(
    if L == {} then return {};
    D1 = {first L};
    scan(L,l->(
	    if dominanceComparisonList(l,D1) == true then D1 = D1|{l};
	    ));
	maximalDominantWeights((unique D1),Schur=>opts.Schur) 
    );


--------------------------------------------------------------------
--------------------------------------------------------------------

entryToHashWrapperDomMonomial = method();
entryToHashWrapperDomMonomial (ZZ,List,List) := (a,B,D) ->(
    H1 := schurBetti(a,B,D);
    applyPairs(H1,(k,v)->(k,dominantWeights monomialSyzygiesWeights(a,B,D,k)))
    )

H1 = entryToHashWrapperDomMonomial(0,{0,0},{2,4})
H2 = applyValues(dominantWeightsBetti(0,{0,0},{2,4}),v->(apply(v,i->i#0)))
M1 = new MutableHashTable
apply(keys H2,k->M1#k = #toList( set (H1#k) - set (H2#k)))
peek M1
M2 = new MutableHashTable
delete(,apply(keys H2,k->(if  ((sort (H1#k)) == (sort (H2#k))) == false then k
	)))

testCases = {{2,2},{2,3},{2,4},{2,5},{2,7},{3,3},{3,4},{3,5}};
M = new MutableHashTable;
apply(testCases, D->(
	H1 = entryToHashWrapperDomMonomial(0,{0,0},D);
	H2 = applyValues(dominantWeightsBetti(0,{0,0},D),v->(apply(v,i->i#0)));
	M#D = delete(,apply(keys H2,k->(if  ((sort (H1#k)) == (sort (H2#k))) == false then k
	)))
	))

 A = H1#(10,2)  
 B = H2#(10,2)
