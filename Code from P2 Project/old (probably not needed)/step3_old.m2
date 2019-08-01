restart
needsPackage "Posets";
needsPackage "SchurRings"
load "veroP2.m2";

-- This post-processing file takes the necesary multigraded Hilbert Series
-- and returns the correct Betti tables. The heart of the code is 
-- makeOutputFiles the other functions are built into that one.
-- e.g. to create the betti table output for n=2,d=5,b=0 just run
-- makeOutputFiles(5,2,0). 
-- There may be issues if n!= 2.

--  We globally define a hilbert series ring to avoid issues with "sub"
A = QQ[t_0,t_1,t_2,MonomialOrder => Lex];

--  This is annoying internal code to translate between vero.m2
--  and the form we were using
standardize = P -> (
    if P == 0 then return {};
    K := listForm P;
    apply(K, k-> if #(k#0) == 3 then k else (k#0|apply(3-#(k#0),i-> 0), k#1))
    )

--  Input:  A polynomial
--  Output:  the exponent vector of the leading term (??)
expM = m->(
    if m == 0 then return {};
    (exponents m)#0
    );

--  Input:  a partition L with 3 entries (hard coded!!!)
--  Output:  the hilbert series of the corresponding Schur functor
hilbSeries = L ->(
    D := schurRing(s,3);
    a := L#0;
    b := L#1;
    c := L#2;
    pole := toE D_{a-c,b-c};
    phi := map(A,ring pole, {t_0+t_1+t_2,t_0*t_1 + t_1*t_2 + t_2*t_0,t_1*t_2*t_0}|apply(6,i->0));
    (t_0*t_1*t_2)^c * phi(pole)
    )

--  This is a subroutine for taking a representation
--  and removing the Schur functor corresponding
--  to its highest weight vector.
--  Input:  A multigraded Hilbert series P
--  Output:  The hilbert series of P - h, where
--  h is the Hilbert series of the Schur functor
--  of the highest weight vector of P.
decomposeHilb1 = P ->(
    if P == 0 then return "P is zero";
    m := expM(leadTerm P);
    h := hilbSeries(m);
    (P - h,m)
    )

--  By iterating decomposeHilb1, this should decompose
--  a representation into Schur functors.
--  It spits out an error if we ever get a negative
--  leading coefficient in the Hilbert series.
--  This allows us to "catch" errors in our Hilbert series computation.
--  Input:  A multigraded Hilbert series P
--  Output:  A pair (L,P) where L is the corresponding
--        list of Schur functors and P is the "error".
--       P should be 0 unless there was an error
--       our computation of multigraded numbers.
decomposeHilb = (P) -> (
    L := {};
    if P == 0 then return {{},0};
    while leadCoefficient P > 0 do(
      	  K := decomposeHilb1(P);
	  L = L|{K#1};
	  P = K#0;
	  );
    (L,P)
    )

--  Just a subroutine for builiding the functions.
--  Input: a list of pairs (partition, multiplicity)
--  Output: the multigraded Hilbert series of the corresponding representation
totalHilb = (K)->(
    if #K == 0 then return 0_A;
    sum apply(K, i-> i_1*hilbSeries(i_0))  
    );

--  Input:  A multigraded hilbert series P
--  Output: P(1,1,...,1), i.e. the corresponding total Betti number
--  Caveat: doesn't seem to work w quotient rings
totalRank = P ->(
    phi := map(ZZ,ring P, apply(numgens ring P, i-> 1));
    phi(P)
    )

--  Input:  A mutable hash table
--  Output:  A hash table
--  Used only for writing out results to file.
unMute = H ->(
    hashTable apply(keys H,h-> h => H#h)
    )

maximalDominantWeights = method()
maximalDominantWeights (List) := (G)->(
    if G == {} then return G;
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
    P := poset(G,cmp,AntisymmetryStrategy => "none");
    M := maximalElements(P)
    );

--  Insanely long code for building all Betti tables
--  Input is e = degree of embedding; 
--           n = dimension of projective space (needs to be 2 for now)
--           b = degree of line bundle you're resolving
--  Output a pair (multigraded Betti table, schur Betti table),
--     where we have assumed no cancelation amoungst Schur functors.
startBetti = (d,n,b)->(
    k = n+1;
    r = b;
    if k != 3 then error "This code only works for P^2";
    koz := veroP2(k,d,r);
    -- running vero clears value of e apparently.
    schurBetti = new MutableHashTable;
    keyList := flatten apply(#koz, p-> apply(#(koz#0),j-> (p,j)));
    scan(keyList, i ->(  schurBetti#i = standardize(koz#(i_0))#(i_1)));
    --  Some weird step to avoid conflicts of variable names.
    --  Building the Betti tables
    multiBetti = new MutableHashTable;
    scan(keyList, i->( multiBetti#i = totalHilb(schurBetti#i)));
    (multiBetti,schurBetti)
    )


--  Input: a multigraded hilbert series, up to perhaps a minor error
--     we try to "fix" the error by rewriting the hilbert series as
--     a sum of hilbert series of schur functors
--  Output: a triple (K,P',E) where K is the guess at the schur functor decomposition,
--     P' = hilb(K) is the "fixed" hilbert series, and E = P - P' is "error"
fixMultiHilb = P ->(
    L := decomposeHilb(P);
    K := pairs tally (L_0);
    (K,totalHilb(K),L_1)
    );

lowerBound = (d,n,b,q)->(
    m := (q*d+b)//(d-1);
    r := (q*d+b)%(d-1);
    binomial(m+d,m)-binomial(m+d-r-1,m)-m
    )

upperBound = (d,n,b,q)->(
    m := (q*d+b)//(d-1);
    r := (q*d+b)%(d-1);
    binomial(n+d,n)+binomial(n-m+r,n-m)-binomial(n-m+d,n-m)-m-1
    )

--- Finds the relevant range based upon EEL non-vanishing.
--- For n>0 this does not quite handle the edge cases correctly. 
--- Daniel: Do you mean n>2 not n>0???
--- e.g. for (3,5,0) this includes (1,1) or (3,5,3) includes (-1,0)
relevantRange = (d,n,b) ->(
    if n == 2 then (
    	flatten apply(ceiling(n+1-(n+b)/d),q->(
		toList apply((lowerBound(d,n,b,q+1)+1..upperBound(d,n,b,q)),i->(i,q))
		))
    )
    else (
	overlap := flatten apply(ceiling(n+1-(n+b)/d),q->(
		toList apply((lowerBound(d,n,b,q+1)+1..upperBound(d,n,b,q)),i->(i,q))
		));
	ends := flatten apply(ceiling(n+1-(n+b)/d),q->({(lowerBound(d,n,b,q)-1,q),(upperBound(d,n,b,q)+1,q)}
		));
	overlap|ends	
    )
)



--adjusts the multigraded Betti numbers in the entries in the relevant range
adjustMultiBetti = (d,n,b,mb)->(
    --L :=  value get concatenate{"correctionList/P",toString(n),"/", toString(n), "_", toString(d), "_", toString(b),".m2"};
    L := relevantRange(d,n,b);
    cL := apply(L, l->((l#0-1,l#1+1)));
    apply(#L,l->(
	    --newEntry := fixMultiHilb  value get concatenate{"mgHilbertSeries/P",toString(n),"/", toString(n), "_", toString(d), "_", toString(b), "-", toString((L#l)#0), "_", toString((L#l)#1), ".txt"};
	    newPart := mb#(L#l) -  value get concatenate{"mgHilbertSeries/P",toString(n),"/", toString(n), "_", toString(d), "_", toString(b), "-", toString((L#l)#0), "_", toString((L#l)#1), ".txt"};
	    oldEntry := mb#(L#l);
	    newPartFixed := fixMultiHilb newPart;
	    mb#(L#l)= oldEntry + newPartFixed_1;
	    mb#(cL#l) = mb#(cL#l) + newPartFixed_1;
    	    ));
    mb 
    )

--adjusts the multigraded Betti numbers in the entries in the relevant range
adjustSchurBetti = (d,n,b,sb)->(
    --L :=  value get concatenate{"correctionList/P",toString(n),"/", toString(n), "_", toString(d), "_", toString(b),".m2"};
    L := relevantRange(d,n,b);
    cL := apply(L, l->((l#0-1,l#1+1)));
    apply(#L,l->(
	    newEntry := fixMultiHilb  value get concatenate{"mgHilbertSeries/P",toString(n),"/", toString(n), "_", toString(d), "_", toString(b), "-", toString((L#l)#0), "_", toString((L#l)#1), ".txt"};
	    oldEntry := mb#(L#l);
	    mb#(L#l)= newEntry_0;
	    mb#(cL#l) = mb#(cL#l) + newEntry_0 - oldEntry;
    	    ));
    mb 
    )


--  adjust multiBetti and schurBetti and creates error file,
--  all without duplicate computation.
--  Input:  parameters d,n,b and startB which is the output of startBetti
--     and so startB_0 is a multigraded mutable hashtable and 
--       startB_1 is a schur mutable hash table.
--  Output:  a triplet of mutable hash tables (multi, schur, error)
--          which incorporate ranks from our computations and correct for error
--          using fixMultiHilb
adjustBetti = (d,n,b,startB) -> (
    mb := startB_0;
    sb := startB_1;
    er := new MutableHashTable;    
    L := relevantRange(d,n,b);
    cL := apply(L, l->((l#0-1,l#1+1)));
    apply(#L,l->(
--	     l = 0
	    newPart := - mb#(L#l) +  value get concatenate{"mgHilbertSeries/P",toString(n),"/", toString(n), "_", toString(d), "_", toString(b), "-", toString((L#l)#0), "_", toString((L#l)#1), ".txt"};
	    oldEntry := mb#(L#l);
	    newPartFixed := fixMultiHilb newPart;
	    --  we adjust multigraded entries
	    mb#(L#l) = oldEntry + newPartFixed_1;
	    mb#(cL#l) = mb#(cL#l) + newPartFixed_1;
            --  we adjust schur entries
    	    sb#(L#l) = sb#(L#l) | newPartFixed_0;
	    sb#(cL#l) = sb#(cL#l) | newPartFixed_0;
	    --  we encode errors
	    er#(L#l) = newPartFixed_2;
    	    ));
    (mb,sb,er)     
    )

--computes the error terms in our numerical methods
erBetti = (d,n,b)->(
    --L :=  value get concatenate{"correctionList/P",toString(n),"/", toString(n), "_", toString(d), "_", toString(b),".m2"};
    L := relevantRange(d,n,b);
    cL := apply(L, l->((l#0-1,l#1+1)));
    er := new MutableHashTable;
    apply(#L,l->(
	    newEntry := fixMultiHilb  value get concatenate{"mgHilbertSeries/P",toString(n),"/", toString(n), "_", toString(d), "_", toString(b), "-", toString((L#l)#0), "_", toString((L#l)#1), ".txt"};
	    er#(L#l) = newEntry_1;
	    ));
    er 
    )

makeOutputFiles =  (d,n,b)->(
    (mb,sb,er) := adjustBetti(d,n,b,startBetti(d,n,b));
    multiBetti := mb;
    schurBetti := sb;
    keyList := keys multiBetti;
    errorBetti := new MutableHashTable;
    scan(keyList, i->  errorBetti#i = 0_A);
    scan(keys er, i-> errorBetti#i = er#i);
    totalBetti := new MutableHashTable;
    scan(keyList, i->  totalBetti#i = totalRank(multiBetti#i));
    domWeights  := new MutableHashTable;
    scan(keyList, i->  domWeights#i = maximalDominantWeights(apply(schurBetti#i, l -> l_0)));
    numReps  = new MutableHashTable;
    scan(keyList, i->  numReps#i = #(schurBetti#i));
    numRepsMulti := new MutableHashTable;
    scan(keyList, i->  numRepsMulti#i = sum apply(schurBetti#i, j-> j_1));
    lexLeadingWeights := new MutableHashTable;
    scan(keyList, i->  lexLeadingWeights#i = if schurBetti#i == {} then {} else last sort apply(schurBetti#i,j-> j_0));
    --errorBetti := erBetti(d,n,b);
    --  Writing the files
    g = openOut ("bettiP2_"|toString d|"_"|toString b|".m2");
    g<< "--This file computes Betti tables for P^2 for d = "|toString d|" and b = "|toString b;
    -- need to add bit to initialize A = QQ[t_0,t_1,t_2];
    g<< endl;
    g<< "--tb stands for Total Betti numbers";
    g<< endl;
    g<< "tb"|toString d|toString b|" = ";
    g<< toExternalString unMute(totalBetti);
    g<< ";";
    g<< endl; 
    g<< "--mb stands for Multigraded Betti numbers";
    g<< endl ;
    g<< "mb"|toString d|toString b|" = ";
    g<< toExternalString unMute(multiBetti);
    g<< ";";
    g<< endl;
    g<< "--sb represents the betti numbers as sums of Schur functors";
    g<< endl ;
    g<< "sb"|toString d|toString b|" = ";
    g<< toExternalString unMute(schurBetti);
    g<< ";";
    g<< endl;
    g<< "--dw encodes the dominant weights in each entry";
    g<< endl ;
    g<< "dw"|toString d|toString b|" = ";
    g<< toExternalString unMute(domWeights);
    g<< ";";
    g<< endl;
    g<< "--lw encodes the lex leading weight in each entry";
    g<< endl ;
    g<< "lw"|toString d|toString b|" = ";
    g<< toExternalString unMute(lexLeadingWeights);
    g<< ";";
    g<< endl;
    g<< "--nr encodes the number of disctinct reprsentations in each entry";
    g<< endl ;
    g<< "nr"|toString d|toString b|" = ";
    g<< toExternalString unMute(numReps);
    g<< ";";
    g<< endl;
    g<< "--nrm encodes the number of representations with multiplicity in each entry";
    g<< endl ;
    g<< "nrm"|toString d|toString b|" = ";
    g<< toExternalString unMute(numRepsMulti);
    g<< ";";
    g<< endl;
    g<< "--er encodes the errors in the computed multigraded Hilbert series via our Schur method in each entry";
    g<< endl ;
    g<< "er"|toString d|toString b|" = ";
    g<< toExternalString unMute(errorBetti);
    g<< ";";
    g<< endl;
    g<< "end;";
    close g;    
    )

time makeOutputFiles(5,2,0)
end;



---
restart
load "step3.m2"
--- Test
makeOutputFiles(4,2,0)
makeBetti = H ->(
    new BettiTally from apply(keys H, h-> (h_0,{h_0+h_1},h_0+h_1)=> H#h)
    )
makeBetti tb50
peek er50
