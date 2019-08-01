restart
needsPackage "SchurRings"
-- A = QQ[t_0,t_1,t_2,MonomialOrder => Lex];
-- --  we globally define a hilbert series ring to avoid issues with "sub"

--  Conjecture:  is that the Betti number
--               in row 1 and column binomial(d+1,2)
--               equals the Schur functor
--              associated to (a,b,c) where
--      a =  binomial(d+2,3)-1
--      b = d*(d^2+5)/6
--  and c = binomial(d+1,3)+1
--   Note that these add up to d*binomial(d+1,2).

-- --  This is a subroutine for taking a representation
-- --  and removing the Schur functor corresponding
-- --  to its highest weight vector.
-- --  Input:  A multigraded Hilbert series P
-- --  Output:  The hilbert series of P - h, where
-- --  h is the Hilbert series of the Schur functor
-- --  of the highest weight vector of P.
-- decomposeHilb1 = P ->(
--     if P == 0 then return "P is zero";
--     m := expM(leadTerm P);
--     h := hilbSeries(m);
--     (P - h,m)
--     )


-- --  By iterating decomposeHilb1, this should decompose
-- --  a representation into Schur functors.
-- --  It spits out an error if we ever get a negative
-- --  leading coefficient in the Hilbert series.
-- --  This allows us to "catch" errors in our Hilbert series computation.
-- --  Input:  A multigraded Hilbert series P
-- --  Output:  A pair (L,P) where L is the corresponding
-- --        list of Schur functors and P is the "error".
-- --       P should be 0 unless there was an error
-- --       our computation of multigraded numbers.
-- decomposeHilb = (P) -> (
--     L := {};
--     if P == 0 then return {{},0};
--     while leadCoefficient P > 0 do(
--       	  K := decomposeHilb1(P);
-- 	  L = L|{K#1};
-- 	  P = K#0;
-- 	  );
--     (L,P)
--     )


-- --  Input:  a partition L with 3 entries (hard coded!!!)
-- --  Output:  the hilbert series of the corresponding Schur functor
-- hilbSeries = L ->(
--     D := schurRing(s,3);
--     a := L#0;
--     b := L#1;
--     c := L#2;
--     pole := toE D_{a-c,b-c};
--     phi := map(A,ring pole, {t_0+t_1+t_2,t_0*t_1 + t_1*t_2 + t_2*t_0,t_1*t_2*t_0}|apply(6,i->0));
--     (t_0*t_1*t_2)^c * phi(pole)
--     )

--  A minor improvement where you can input a ring A = QQ[t_0,t_1,t_2]
--  so the output always lands in A and thus you can do arithmetic
--  with the Hilbert series more easily.  Very hacky.
--betterHilb = (L,A)->(
--    P = hilbSeries(L);
--    phi = map(A,ring P, vars A);
--    phi(P)
--    );


-- --  Just a subroutine for builiding the functions.
-- --  Input is K an entry of schurBetti
-- --  Output the multigraded Hilbert series.
-- totalHilb = (K)->(
--     if #K == 0 then return 0_A;
--     sum apply(K, i-> i_1*hilbSeries(i_0))  
--     );

-- --  Input:  A hilbert series P
-- --  Output:  P(1,1,...,1), i.e. the corresponding total Betti number
-- --  Caveat: doesn't seem to work w quotient rings
-- totalRank = P ->(
--     phi = map(ZZ,ring P, apply(numgens ring P, i-> 1));
--     phi(P)
--     )
--
--
--  The next set of codes compares our computations with
--  the results from Ein-Erman-Lazarsfeld
--  This had issues!!!
--
OldEELWeights = (d,p,q,b,n)->(
    S = QQ[x_0..x_n, MonomialOrder => Lex];
    R = S/ideal apply(n+1,i-> x_i^d);
    f = product apply(q,i-> x_i^(d-1))* x_(q)^(q+b);
    I = ann(f);
    L = flatten entries gens image basis(d,I);
    L' = apply(L, l-> sub(l,S));
    K = subsets(L',p);
    apply(K,k->  expM(sub(f,S)*product(k)))
    )    

--  Input:  A polynomial
--  Output:  the exponent vector of the leading term (??)
-- expM = m->(
--     if m == 0 then return {};
--     (exponents m)#0
--     );

topWeight = L ->(
    first reverse sort L
    );
--
--  Multigraded Betti Numbers for non-cancellative entries.
--  We've embedded P^n into P^N by O(d).
--
--  This codes produces the denominator of the hilbert Series 
--  of P^N itself. 
denom = (n,d,b) ->(
    N = binomial(n+d,d);
    B = A/(ideal gens A)^(d*(N-1)+b+1);
    weights = flatten entries basis(d,B);
    product apply(weights, w-> (1-w))
    )

--  This code produces the Hilbert series (as a truncated series)
--  of P^n with respect to O(d), so only including every d'th entry
--  the ring A is hardcoded which is bad.
veroHilb = (d,b,n) ->(
    N = binomial(n+d,d);
    --weights = sum flatten entries basis(d,A);
    sum apply(N,i-> sum flatten entries basis(d*i+b,B))
    )

--  Mutiplying denom and veroHilb should give the 
--  numerator of the rational function representation
--  of the Hilbert series of P^n with resepct to O(d)
--  which is the alternating sum of the Betti numbers.
--  The truncation removes the "noise" introdcued
--  by the fact that veroHilb was not actually an infinite series.
Knum = (d,b,n) -> (
    D1 = denom(n,d,b);
    V1 = veroHilb(d,b,n);
    D1*V1
    )  


--  This code list the various
--  multigraded Betti numbers,
--  avoiding the need to recompute at each step.
--  TODO:  we want the output to be in the ambient ring of P
eulerBettiList = (d,L,b,n) ->(
    K2 = Knum(d,b,n);
    OUT = apply(L,l->(
    	    P = sum apply(terms K2, i-> if degree(i) == {d*(l)+b} then i else 0 );
	    if P != 0 then P = sub(P, ambient(ring P));
    	    if P == 0 then return P;
    	    if leadCoefficient P < 0 then P = (-1)*P;
    	    P
	    ));
	OUT
    )

--  Finally, this code picks off the actual multigraded
--  Betti number, at least up to alternating sums.
eulerBetti = (d,p,q,b,n) ->(
    K2 = Knum(d,b,n);
    P = sum apply(terms K2, i-> if degree(i) == {d*(p+q)} then i else 0 );
    if P == 0 then return P;
    if leadCoefficient P < 0 then P = (-1)*P;
    P
    )

-- --  Input:  A mutable hash table
-- --  Output:  A hash table
-- --  Used only for writing out results to file.
-- unMute = H ->(
--     hashTable apply(keys H,h-> h => H#h)
--     )

-- --  This is annoying internal code to translate between vero.m2
-- --  and the form we were using
-- standardize = P -> (
--     if P == 0 then return {};
--     K := listForm P;
--     apply(K, k-> if #(k#0) == 3 then k else (k#0|apply(3-#(k#0),i-> 0), k#1))
--     )


-- needsPackage "Posets";
-- maximalDominantWeights = method()
-- maximalDominantWeights (List) := (G)->(
--     if G == {} then return G;
--     cmp = (a, b) -> (
--         if #b > #a then return false;
--         sa := 0;
--         sb := 0;
--         for k from 0 to #b - 1 do (
--             sa = sa + a_k;
--             sb = sb + b_k;
--             if sa > sb then return false;
--             );
--         true
--         );
--     P = poset(G,cmp,AntisymmetryStrategy => "none");
--     M = maximalElements(P)
--     );


--  Turns a Hash Table into a Betti Tally
--  for easy viewing.
makeBetti = H ->(
    new BettiTally from apply(keys H, h-> (h_0,{h_0+h_1},h_0+h_1)=> H#h)
    )




betterDominantWeights = K->(
    if K == {} then return {};
    L = reverse sort K;
    dWts = {first L};
    scan(L,l->(
	    if cmpL(l,dWts) == false then dWts = dWts|{l};
	    ));
    maximalDominantWeights(dWts)    
    );

--  Input: two weight
--  Output:  true if b dominates a
cmp = (a, b) -> (
        if #b > #a then return false;
        sa := 0;
        sb := 0;
        for k from 0 to #b - 1 do (
            sa = sa + a_k;
            sb = sb + b_k;
            if sa > sb then return false;
            );
        true
        );
    
--  Input: a weight a and a list of weights K
--  Output:  true if any element in K dominates a
cmpL = (a,K) -> (
    for k in K do if cmp(a,k) == true then return true;
    false
    );




end;







