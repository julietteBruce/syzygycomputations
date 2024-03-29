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
-*monomialSyzygies = method();
monomialSyzygies (ZZ,List,List,Sequence) := (a,B,D,P) ->(
    p = P#0;
    q = P#1;
    if a != 0 or B != {0,0} then error;
    ---
    S = QQ[x_0,x_1,y_0,y_1,Degrees=>{{1,0},{1,0},{0,1},{0,1}}];
<<<<<<< HEAD
    I = ideal(x_0^(D#0)*y_0^(D#1), x_0^(D#0)*y_1^(D#1)+x_1^(D#0)*y_0^(D#1), x_1^(D#0)*y_1^(D#1));
=======
    I := ideal(x_0^(D#0)*y_0^(D#1), x_0^(D#0)*y_1^(D#1)+x_1^(D#0)*y_0^(D#1), x_1^(D#0)*y_1^(D#1));
>>>>>>> parent of 67e86b7... updating data
    ---
    if q == 0 then return {};
    if q == 1 then (f = x_0^(D#0+B#0)*y_0^((D#1)-1)*y_1^((B#1)+1)) ;
    --if q == 2 then (f = x_0^(2*(D#0)-1)*x_1^((B#0)+1)*y_0^((D#1)-1)*y_1^((D#1)+(B#1)+1));
    if q == 2 then (f = x_0*x_1^(2*(D#0)-1)*y_0^((D#1)+1)*y_1^((D#1)+(B#1)-1));
    --if q == 2 then (f = x_0^((D#0)+1)*x_1^((D#0)-1)*y_0*y_1^(2*(D#1)+(B#1)-1));
    ---	
    Q = quotient(I,f);
    ---
    B1 = unique flatten entries basis(D,S);
<<<<<<< HEAD
    B2 := delete(,apply(B1, m->(if m%Q == 0 then m)));
    L1 := toList (set B2 - set {x_0^(D#0)*y_1^(D#1), x_1^(D#0)*y_1^(D#1), x_0^(D#0)*y_0^(D#1)});
=======
    B2 = delete(,apply(B1, m->(if m%Q == 0 then m)));
    L1  = toList (set B2 - set {x_1^(D#0)*y_0^(D#1), x_1^(D#0)*y_1^(D#1), x_0^(D#0)*y_0^(D#1)});
>>>>>>> parent of 67e86b7... updating data
    ---
    R = QQ[x_0,x_1,y_0,y_1,Degrees=>{{1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1}}];
    J = ideal(x_0^(D#0)*y_0^(D#1), x_0^(D#0)*y_1^(D#1)+x_1^(D#0)*y_0^(D#1), x_1^(D#0)*y_1^(D#1));
    T = R/J;
    L2 := apply(L1,i->sub(i,T));
    B3 := apply(B1,i->sub(i,T));
    g := sub(f,T);
    ---
    D1 = delete(,apply(B3,m->(if g%m == 0 then m)));
<<<<<<< HEAD
    L4 = apply({x_0^(D#0)*y_1^(D#1), x_1^(D#0)*y_1^(D#1), x_0^(D#0)*y_0^(D#1)},i->sub(i,T));
    D2 = toList (set D1 - set L4);
    L3 = toList (set L2 - set D2);
    apply(subsets(L3,p-(#D2)),L->(
=======
    D2 = toList (set D1 - set {x_1^(D#0)*y_0^(D#1), x_1^(D#0)*y_1^(D#1), x_0^(D#0)*y_0^(D#1)});
    L3 = toList (set L2 - set D2);
    apply(subsets(L3,p-(#D1)),L->(
>>>>>>> parent of 67e86b7... updating data
		{apply(L,i->degree i)|apply(D2,j->degree j),degree g}
		))
	) *-

<<<<<<< HEAD
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
    if q == 1 then (f = x_0^(D#0+B#0)*y_0^((D#1)-1)*y_1^((B#1)+1)) ;
    if q == 2 then (f = x_0^(2*(D#0)-1)*x_1^((B#0)+1)*y_0^((D#1)-1)*y_1^((D#1)+(B#1)+1));
    --if q == 2 then (f = x_0*x_1^(2*(D#0)-1)*y_0^((D#1)+1)*y_1^((D#1)+(B#1)-1));
    --if q == 2 then (f = x_0^((D#0)+1)*x_1^((D#0)-1)*y_0*y_1^(2*(D#1)+(B#1)-1));
    ---	
    Q = quotient(I,f);
    ---
    B1 = unique flatten entries basis(D,S);
    B2 = unique flatten entries basis(q*D,S);
    B3 = delete(,apply(B1, m->(if m%Q == 0 then m)));
    L1 = toList (set B3 - set {x_0^(D#0)*y_1^(D#1), x_1^(D#0)*y_1^(D#1), x_0^(D#0)*y_0^(D#1)});
    ---
    R = QQ[x_0,x_1,y_0,y_1,Degrees=>{{1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1}}];
    J = ideal(x_0^(D#0)*y_0^(D#1), x_0^(D#0)*y_1^(D#1)+x_1^(D#0)*y_0^(D#1), x_1^(D#0)*y_1^(D#1));
    --T = R/J;
    L2 = apply(L1,i->sub(i,R));
    B4 = apply(B1,i->sub(i,R));
    g = sub(f,R);
    B5 = apply(B2,i->(sub(i,R)));
    F = unique delete(,apply(B5,h->(if (g-h)%J == 0 or (g+h)%J == 0 then h)));
    ---
    D1 = unique delete(,flatten apply(F,g->(apply(B4,m->(if g%m == 0 then m)))));	
    L3 = apply({x_0^(D#0)*y_1^(D#1), x_1^(D#0)*y_1^(D#1), x_0^(D#0)*y_0^(D#1)},i->sub(i,R));
    D2 = toList (set D1 - set L3);
    L4 = toList (set L2 - set D2);
    apply(subsets(L4,p-(#D2)),L->(
		{apply(L,i->degree i)|apply(D2,j->degree j),degree g}
		))
	)

----- true for Q1=>1 but not Q1=0,2
H1 = entryToHashWrapperDomMonomial(0,{0,0},{2,3},Q1=>1);
H2 = applyValues(dominantWeightsBetti(0,{0,0},{2,3}),v->(apply(v,i->i#0)));
M = new MutableHashTable
apply(keys H2,k->M#k = isSubset((H1#k),(H2#k)))
apply(keys H2,k->M#k = (H1#k)==(H2#k))


----- true for Q1=>1 but not Q1=0,2
H1 = entryToHashWrapperDomMonomial(0,{0,0},{2,6},Q1=>1);
H2 = applyValues(dominantWeightsBetti(0,{0,0},{2,6}),v->(apply(v,i->i#0)));
M = new MutableHashTable
apply(keys H2,k->M#k = isSubset((H1#k),(H2#k)))  
apply(keys H2,k->M#k = (H1#k)==(H2#k))

----- true for Q1=>1 but not Q1=0,2
H1 = entryToHashWrapperDomMonomial(0,{0,0},{2,3},Q1=>1);
H2 = applyValues(dominantWeightsBetti(0,{0,0},{2,3}),v->(apply(v,i->i#0)));
M = new MutableHashTable
apply(keys H2,k->M#k = isSubset((H1#k),(H2#k)))

----- false for everything but best for Q1=>1
H1 = entryToHashWrapperDomMonomial(0,{0,0},{3,3},Q1=>1);
H2 = applyValues(dominantWeightsBetti(0,{0,0},{3,3}),v->(apply(v,i->i#0)));
M = new MutableHashTable
apply(keys H2,k->M#k = isSubset((H1#k),(H2#k)))

----- false for everything but best for Q1=>1
H1 = entryToHashWrapperDomMonomial(0,{0,0},{3,4},Q1=>1);
H2 = applyValues(dominantWeightsBetti(0,{0,0},{3,4}),v->(apply(v,i->i#0)));
M = new MutableHashTable
apply(keys H2,k->M#k = isSubset((H1#k),(H2#k)))
  
  
----- false for everything but best for Q1=>1
H1 = entryToHashWrapperDomMonomial(0,{0,0},{3,5},Q1=>1);
H2 = applyValues(dominantWeightsBetti(0,{0,0},{3,5}),v->(apply(v,i->i#0)));
M = new MutableHashTable
apply(keys H2,k->M#k = isSubset((H1#k),(H2#k)))  	
		
		
		
		
		
monomialSyzygies = method(Options => {Q1 => 0});
monomialSyzygies (ZZ,List,List,Sequence) := opts -> (a,B,D,P) ->(
    p := P#0;
    q := P#1;
    if a != 0 or B != {0,0} then error;
    ---
    S = QQ[x_0,x_1,y_0,y_1,Degrees=>{{1,0},{1,0},{0,1},{0,1}}];
    I = ideal(x_0^(D#0)*y_0^(D#1), x_0^(D#0)*y_1^(D#1)+x_1^(D#0)*y_0^(D#1), x_1^(D#0)*y_1^(D#1));
    ---
    if q == 0 then return {};
    if q == 1 then (
	if opts.Q1 == 0 then (F = {x_0^(D#0-1)*x_1^(B#0+1)*y_0^((D#1+B#1))});
	if opts.Q1 == 1 then (F = {x_0^(D#0+B#0)*y_0^((D#1)-1)*y_1^((B#1)+1)});
	if opts.Q1 == 2 then (F = {x_0^(D#0-1)*x_1^(B#0+1)*y_0^((D#1+B#1)),x_0^(D#0+B#0)*y_0^((D#1)-1)*y_1^((B#1)+1)});
	);
    --if q == 1 then (F = {x_0^(D#0-1)*x_1^(B#0+1)*y_0^((D#1+B#1))});
     --if q == 1 then (F = {x_0^(D#0+B#0)*y_0^((D#1)-1)*y_1^((B#1)+1)});
    --if q == 1 then (F = {x_0^(D#0-1)*x_1^(B#0+1)*y_0^((D#1+B#1)),x_0^(D#0+B#0)*y_0^((D#1)-1)*y_1^((B#1)+1)}) ;
    if q == 2 then (F = {x_0^(2*(D#0)-1)*x_1^((B#0)+1)*y_0^((D#1)-1)*y_1^((D#1)+(B#1)+1)});
    --if q == 2 then (f = x_0*x_1^(2*(D#0)-1)*y_0^((D#1)+1)*y_1^((D#1)+(B#1)-1));
    --if q == 2 then (f = x_0^((D#0)+1)*x_1^((D#0)-1)*y_0*y_1^(2*(D#1)+(B#1)-1));
    ---	
    flatten apply(F, f->(
	    Q := quotient(I,sub(f,S));
	    print Q;
    	    ---
    	    B1 := unique flatten entries basis(D,S);
    	    B2 := unique flatten entries basis(q*D,S);
    	    B3 := delete(,apply(B1, m->(if m%Q == 0 then m)));
	    print (#B3);
	    L0 := apply({x_0^(D#0)*y_1^(D#1), x_1^(D#0)*y_1^(D#1), x_0^(D#0)*y_0^(D#1)},i->sub(i,S));
    	    L1 := toList (set B3 - set L0);
	    print (#L1);
    	    ----
    	    R := QQ[x_0,x_1,y_0,y_1,Degrees=>{{1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1}}];
    	    J := ideal(x_0^(D#0)*y_0^(D#1), x_0^(D#0)*y_1^(D#1)+x_1^(D#0)*y_0^(D#1), x_1^(D#0)*y_1^(D#1));
    	    --T = R/J;
    	    L2 := apply(L1,i->sub(i,R));
    	    B4 := apply(B1,i->sub(i,R));
    	    g := sub(f,R);s
    	    B5 := apply(B2,i->(sub(i,R)));
    	    H := unique delete(,apply(B5,h->(if (g-h)%J == 0 or (g+h)%J == 0 then h)));	  
    	    ---
    	    D1 := unique delete(,flatten apply(H,g->(apply(B4,m->(if g%m == 0 then m)))));	
    	    L3 := apply({x_0^(D#0)*y_1^(D#1), x_1^(D#0)*y_1^(D#1), x_0^(D#0)*y_0^(D#1)},i->sub(i,R));
    	    D2 := toList (set D1 - set L3);
    	    L4 := toList (set L2 - set D2);
	    print (#L4);
    	    apply(subsets(L4,p-(#D2)),L->(
		    {apply(L,i->degree i)|apply(D2,j->degree j),degree g}
		    ))
	    ))
	)
{8, 7, 19, 1}
M1 = unique monomialSyzygiesWeights(0,{0,0},{2,3},(4,1),Q1=>0);
M2 = unique monomialSyzygiesWeights(0,{0,0},{2,3},(4,1),Q1=>1);
M3 = unique M1|M2;
M4 = unique monomialSyzygiesWeights(0,{0,0},{2,3},(4,1),Q1=>2);

{{2, 0, 2, 1}, {1, 1, 1, 2}, {2, 0, 1, 2}, {1, 1, 3, 0}}|{1, 1, 3, 0}|


(unique M3) == (unique M4)
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
    
=======
>>>>>>> parent of 67e86b7... updating data
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
monomialSyzygiesWeights = method(Options => {Q1 => 0});
monomialSyzygiesWeights (ZZ,List,List,Sequence) := opts -> (a,B,D,P) ->(
    L := monomialSyzygies(a,B,D,P,Q1=>opts.Q1);
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

entryToHashWrapperDomSchur = method();
entryToHashWrapperDomSchur (ZZ,List,List) := (a,B,D) ->(
    H1 := schurBetti(a,B,D);
    applyPairs(H1,(k,v)->(k,dominantWeights(v,Schur=>true)))
    )

entryToHashWrapperMonomial = method();
entryToHashWrapperMonomial (ZZ,List,List) := (a,B,D) ->(
    H1 := schurBetti(a,B,D);
    applyPairs(H1,(k,v)->(k,monomialSyzygiesWeights(a,B,D,k)))
    )

entryToHashWrapperDomMonomial = method(Options => {Q1 => 0});
entryToHashWrapperDomMonomial (ZZ,List,List) := opts -> (a,B,D) ->(
    H1 := schurBetti(a,B,D);
    applyPairs(H1,(k,v)->(k,dominantWeights monomialSyzygiesWeights(a,B,D,k,Q1=>opts.Q1)))
    )


fileName1 = (B,D)->(
    toString(B#0)|toString(B#1)|toString(D#0)|toString(D#1)
    )
fileName = (B,D)->(
    toString(B#0)|"_"|toString(B#1)|"_"|toString(D#0)|"_"|toString(D#1)
    )


updateOutputFilesDW =  (a,B,D)->(
    --  Writing the files
    g = openOut ("HirzebruchSyzygiesNew/bettiF0_"|fileName(B,D)|".m2");
    g<< "A := QQ[t_0,t_1,t_2,t_3];";
    g<< endl;
    g<< "--tb stands for Total Betti numbers";
    g<< endl;
    g<< "tb"|fileName1(B,D)|" = ";
    g<< toExternalString totalBetti(a,B,D);
    g<< ";";
    g<< endl; 
    g<< "--sb represents the betti numbers as sums of Schur functors";
    g<< endl ;
    g<< "sb"|fileName1(B,D)|" = ";
    g<< toExternalString schurBetti(a,B,D);
    g<< ";";
    g<< endl;
    g<< "--dw stands for dominant weights";
    g<< endl;
    g<< "dw"|fileName1(B,D)|" = ";
    g<< toExternalString entryToHashWrapperDomSchur(a,B,D);
    g<< ";";
    g<< endl;
    --g<< "--mw stands for monomial weights";
    --g<< endl;
    --g<< "mw"|fileName1(B,D)|" = ";
    --g<< toExternalString entryToHashWrapperMonomial(a,B,D);
    --g<< ";";
    --g<< endl;
    g<< "--dmw stands for dominant monomial weights";
    g<< endl;
    g<< "dmw"|fileName1(B,D)|" = ";
    g<< toExternalString entryToHashWrapperDomMonomial(a,B,D);
    g<< ";";
    g<< endl;
    g<< "end;";
    close g;    
    )


dataRange = value get "dataRange.m2"
dataRange2 = delete({{2,0},{3,11}},delete({{2,1},{3,11}},delete({{1,2},{2,15}},delete({{1,2},{2,14}},delete({{1,0},{2,13}},dataRange)))))
apply(dataRange2,i->(
	print i;
	updateOutputFilesDW(0,(i#0),(i#1));
	    ))