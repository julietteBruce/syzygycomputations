restart
needsPackage "Posets";
needsPackage "SchurRings"

-- This file combines the pre-processing and post-processing steps
-- of our code.  In particular, it takes the necesary multigraded 
-- Hilbert Series, which is computed via Matlab, and it incorporates
-- that to produce a correct Betti table.
--
-- The heart of the code is makeOutputFiles which hardcodes the finalized
-- data on the Betti tables.  All other functions feed into that one.
--
-- To create the betti table output for n=2,d=5,b=0 just run
-- makeOutputFiles(5,2,0). It assumes that we have placed corrected
-- multigraded Hilbert series for entries in the relevant as .txt files
-- in a certain directory, and thus this code cannot be run without the
-- correct directories.
--
-- There are many issues, even for subroutines if n != 2.  Most of these
-- issues stem from the veroP2 code, though there may be others.


--
--  Table of Contents
--
--  Section 0: Some global definitions
--  Section 1: Basic Functions for Manipulating Hilbert Series 
--             and Schur Functors
--  Section 2: Functions for Building schur and multigraded Betti tables
--             outside of relevant range
--  Section 3: Functions for correcting multigraded/schur entries based on
--             our large scale rank computations.
--  Section 4: Functions for other statistics, like dominant weights,
--             Boij-Soederberg coefficients, etc.
--  Section 5: Functions to create the actual output files
--
--  Appendix 1: Claudiu's code for estimating schur Betti table assuming no cancellations
--  Appendix 2: Helper functions for Boij-Soederberg coefficients.  Some
--              were written by Courtney Gibbons
--
--




--
--
--
--  Section 0: Some global definitions
--
--
--
--  We globally define a hilbert series ring to avoid issues with "sub"
A = QQ[t_0,t_1,t_2,MonomialOrder => Lex];

--  Input:  A mutable hash table
--  Output:  A hash table
--  Used only for writing out results to file.
unMute = H ->(
    hashTable apply(keys H,h-> h => H#h)
    )

--  Input: the hash table for total betti numbers
--  Output:  the corresponding BettiTally
makeBetti = H ->(
    new BettiTally from apply(keys H, h-> (h_0,{h_0+h_1},h_0+h_1)=> H#h)
    )



--
--
--
--  Section 1: Basic Functions for Manipulating Hilbert Series 
--             and Schur Functors
--
--
--



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


--  Input: a list of pairs (partition, multiplicity)
--  Output: the multigraded Hilbert series of the corresponding representation
--  Note:  this code was taking too long on large examples, so we created
--         some hacky multistep stuff below.
totalHilb = (K)->(
    if #K == 0 then return 0_A;
    sum apply(K, i-> i_1*hilbSeries(i_0))  
    );

--  Input: a list of pairs (partition, multiplicity)
--  Output:  A mutable hash table M where M#k is the sublist of the 
--    input grouped by last coordinate of the partition.
--    This was done to streamline the Hilbert Series code, which requires 
--    multiplication by t_0*t_1*t_2^(last coordinate)
groupByLastCoord = K ->(
    a := min apply(K,k-> min (k_0));
    b := max apply(K,k-> min (k_0));
    MM := new MutableHashTable;
    scan(toList(a..b), i -> MM#i = {});
    scan(K, k-> MM#(min (k_0)) = MM#(min (k_0))|{k});
    MM
    )


--  Input:  a list of pairs (parition,multiplicity) where all partitions
--          have the same last coordinate
--  Output:  the total Hilbert series of those partitions
hilbLastCoord = NN ->(
    if NN == {} then return 0_A;
    D := schurRing(s,3);
    c := min(NN#0#0);
    pole := sum apply(NN,k-> k_1*(toE (D_{k#0#0-k#0#2,k#0#1-k#0#2})));
--    pole := toE(sum apply(NN,k-> k_1*( (D_{k#0#0-k#0#2,k#0#1-k#0#2}))));
    phi := map(A,ring pole, {t_0+t_1+t_2,t_0*t_1 + t_1*t_2 + t_2*t_0,t_1*t_2*t_0}|apply(6,i->0));
    (t_0*t_1*t_2)^c * phi(pole)    
    )


--  Input:  a mutable hash table where all outputs are lists of pairs
--        (parition,multiplicity) where all partitions have the same 
--        last coordinate
--  Output:  the total Hilbert series of all partitions
--  Note:  those was done to speed up totalHilb, which is a bottleneck step
--         due to the large number of multiplications required
newTotalHilb = K ->(
    if K == {} then return 0_A;
    MM = groupByLastCoord(K);
    sum apply(keys MM, i-> hilbLastCoord(MM#i))
    )


--  Input:  A multigraded hilbert series P
--  Output: P(1,1,...,1), i.e. the corresponding total Betti number
--  Caveat: doesn't seem to work w quotient rings
totalRank = P ->(
    phi := map(ZZ,ring P, apply(numgens ring P, i-> 1));
    phi(P)
    )


--  Input:  An integer n and a list of pairs (partition, multiplicity)
--  Output: The dimension of the corresponding representation
schurRank = (n,L) ->(
    S := schurRing(QQ,s,n+1);
    sum apply(L,l-> l#1*(dim s_(l#0)))
    )


--  Input:  A polynomial
--  Output:  the exponent vector of the leading term
expM = m->(
    if m == 0 then return {};
    (exponents m)#0
    );




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
	  L = L|{K#1};
	  P = K#0;
	  );
    (L,P)
    )




--
--
--
--
--  Section 2:  Functions for Building schur and multigraded Betti tables
--              outside of relevant range
--              See Appendix 1 for complicated code underlying startBetti
--
--
--





--  Input is d = degree of embedding; 
--           n = dimension of projective space (needs to be 2 for now)
--           b = degree of line bundle you're resolving
--  Output a pair (multigraded Betti table, schur Betti table),
--     where we have assumed no cancelation amoungst Schur functors.
startBetti = (d,n,b)->(
    k := n+1;
    r := b;
    if k != 3 then error "This code only works for P^2";
    koz := veroP2(k,d,r);
    -- running vero clears value of e apparently.
    schBetti := new MutableHashTable;
    keyList := flatten apply(#koz, p->(
	    apply(#(koz#0),j-> (p,j))));
    scan(keyList, i ->(  schBetti#i = standardize(koz#(i_0))#(i_1)));
    --  Some weird step to avoid conflicts of variable names.
    --  Building the Betti tables
    multiBetti := new MutableHashTable;
    --scan(keyList, i->( multiBetti#i = totalHilb(schBetti#i)));
    scan(keyList, i->( multiBetti#i = newTotalHilb(schBetti#i)));
    (multiBetti,schBetti)
    )


--  Input: d = degree of embedding; 
--         n = dimension of projective space (needs to be 2 for now)
--         b = degree of line bundle you're resolving
--  Output: Schur Betti table, where we have assumed no cancelation 
--        amoungst Schur functors.
startBettiSchur = (d,n,b)->(
    k := n+1;
    r := b;
    if k != 3 then error "This code only works for P^2";
    koz := veroP2(k,d,r);
    -- running vero clears value of e apparently.
    schBetti := new MutableHashTable;
    keyList := flatten apply(#koz, p->(
	    apply(#(koz#0),j-> (p,j))));
    scan(keyList, i ->(  schBetti#i = standardize(koz#(i_0))#(i_1)));
    schBetti
    )



--  Input: (d,n,b,q)
--  Output:  lower bound of nonvanishing range for K_(p,q)(PP^n, d;b)
lowerBound = (d,n,b,q)->(
    m := (q*d+b)//(d-1);
    r := (q*d+b)%(d-1);
    binomial(m+d,m)-binomial(m+d-r-1,m)-m
    )

--  Input: (d,n,b,q)
--  Output:  upper bound of nonvanishing range for K_(p,q)(PP^n, d;b)
upperBound = (d,n,b,q)->(
    m := (q*d+b)//(d-1);
    r := (q*d+b)%(d-1);
    binomial(n+d,n)+binomial(n-m+r,n-m)-binomial(n-m+d,n-m)-m-1
    )



--- Finds the relevant range based upon EEL non-vanishing.
--- For n>0 this does not quite handle the edge cases correctly. 
--
--- Daniel: Do you mean n>2 not n>0???
--
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







--
--
--
--  Section 3:  Functions for correcting multigraded/schur entries based on
--              our large scale rank computations.
--
--
--



--  Input: a multigraded hilbert series, up to perhaps a minor error
--     we try to fix the error by rewriting the hilbert series as
--     a sum of hilbert series of schur functors
--  Output: a triple (K,P',E) where K is the guess at the schur functor decomposition,
--     P' = hilb(K) is the fixed hilbert series, and E = P - P' is error
fixMultiHilb = P ->(
    L := decomposeHilb(P);
    K := pairs tally (L_0);
    (K,totalHilb(K),L_1)
    );


--  Input:  parameters d,n,b and startB which is the output of startBetti
--     and so startB_0 is a multigraded mutable hashtable and 
--       startB_1 is a schur mutable hash table.
--  Output:  a triplet of mutable hash tables (multi, schur, error)
--          which incorporate ranks from our computations and correct for error
--          using fixMultiHilb
--
--  NB: We correct mb/sb and create er simultaneously to avoid duplicate computation
--
--  something is wrong in adjustBetti.
--
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
	    -- this was just concatenation; we need to add them
	    --sb#(L#l) addPartitions(sb#(L#l),newPartFixed_0);
	    --sb#(cL#l) = sb#(cL#l) | newPartFixed_0;
	    sb#(L#l) = sb#(L#l) | newPartFixed_0;
	    sb#(cL#l) = sb#(cL#l) | newPartFixed_0;
	    --  we encode errors
	    er#(L#l) = newPartFixed_2;
    	    ));
    (mb,sb,er)     
    )



-- Input: parameters, d,n,b
-- Output: Computes the dual of O(b) on P^n embedded by O(d)
--  Note:  may not depend on b.
dualMonomial = (d,n,b) ->(
   b' := ((-1)*b - 3) % d;
   q := regTriplet(d,n,b);
   bigR := ((binomial(d+n,d)-3+q)*d+b'+b);
   assert(bigR % 3 == 0);
   r := (bigR)//(n+1);    
   product apply(n+1, i-> t_i^r)
   )

-- Input: parameters, d,n,b
--  Output: the regularity
--  Note:  you can use this to compute the r parameter in dualHilb below.
--         but then you'll have to change some of your code.
regTriplet = (d,n,b)->(
    b' := ((-1)*b - 3) % d;
--    q = (3 - (b+b'+3) // d)
     (3 - (b+b'+3) // d)
    )
--

-- Input: parameters, d,n,b,r where r is regularity
-- Output: Computes the multi-graaded Hilbert series of O(b) on P^n embedded by O(d)
-- from the the multi-graaded Hilbert series of the dual O(b').
-- This is hacky.
-- Note: We must have an M2OutputFile for the dual table already saved.
dualHilb = (d,n,b)->(
   r := regTriplet(d,n,b);
   L := relevantRange(d,n,b);
   dD := d;
   dN := n;
   dB := (-b-n-1)%d;
   apply(L, l->(
	   dP := binomial(n+d,d)-n-l#0-1;
	   dQ := r-l#1;
   	   --P := value get concatenate{"mgHilbertSeries/P",toString(dN),"/", toString(dN), "_", toString(dD), "_", toString(dB), "-", toString(dP), "_", toString(dQ), ".txt"};
   	   load concatenate{"M2OutputFiles/bettiP", toString(dN),"_",toString(dD),"_",toString(dB),".m2"};
	   P := value concatenate{"mb",toString(dD),toString(dB),"#(",toString(dP),",",toString(dQ),")"};
	   dG := dualMonomial(d,n,b);
   	   (M,C) := coefficients P;
   	   ML := matrix {apply(flatten entries M, m-> sub(dG/m, ring m))};
   	   Q := (ML*C)_(0,0);
   	   --  write Q to a file.
   	   g = openOut ("mgHilbertSeries/P"|toString(n)|"/"|toString(n)|"_"|toString d|"_"|toString b|"-"|toString((l#0))|"_"|toString(l#1)|".txt");
   	   g << toExternalString Q;
   	   close g
	   ))
   )



--
--
--
--  Section 4:  Functions for other statistics, like dominant weights,
--              Boij-Soederberg coefficients, etc.
--
--
--


--  Input: A list of weights
--  Output: The dominant weights from that list
--  NB:  this is multi-step code, done to speed it up.
betterDominantWeights = K->(
    if K == {} then return {};
    L = reverse sort K;
    dWts = {first L};
    scan(L,l->(
	    if cmpL(l,dWts) == false then dWts = dWts|{l};
	    ));
    maximalDominantWeights(dWts)    
    );


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



fixDuplicateSchurs = L->(
    oldGuys = {};
    newGuys = {};
    scan(L, i-> if ring(i#1) === QQ then oldGuys = oldGuys|{i} else newGuys = newGuys|{i});
    oldKeys = apply(oldGuys, i-> i#0);
    OK = new MutableList from oldGuys;
    scan(newGuys, i-> scan(#oldKeys, j->(
	    	if i#0 == oldKeys#j then OK#j = (OK#j#0, OK#j#1 + i#1)
		)));
    nonDups = {};
    scan(newGuys, i->(
	    if isSubset({i#0},oldKeys) == false then nonDups = nonDups|{i}));
    (toList OK)|nonDups
    );

fixAllDuplicateSchur = SB->(
    KEYS = keys SB;
    apply(KEYS, k-> SB#k = fixDuplicateSchurs(SB#k));
    SB
    );



--
--
--
--  Section 5:  Functions that create the actual Output files
--
--
--


--  Overview:  The main function is "makeOutputFiles" which builds
--             the output fields, assuming we've created the multigraded
--             Hilbert series as a doc in the correct folder, for all entries
--             in the relevant range.
--             In cases where we were not able to compute the multigraded Hilb
--             series, like d>6, we've created partial output files.
--             There are two routines "makePartialOutput" and 
--             "makePartialNoMultiOutput" for those cases.  The second skips
--             the multigraded Betti table entirely, because that is the most
--             expensive computation.
--             There is also a routine for P1 called "makeP1OutputFiles".



--  Input:  parameters, d,n,b
--  Output: M2 file of syzygy data for O(b) on P^n embedded by O(d)
makeOutputFiles =  (d,n,b)->(
    (mb,sb,er) := adjustBetti(d,n,b,startBetti(d,n,b));
    multiBetti := mb;
    schBetti := fixAllDuplicateSchur(sb);
    --schBetti := sb;    
    keyList := keys multiBetti;
    errorBetti := new MutableHashTable;
    scan(keyList, i->  errorBetti#i = 0_A);
    scan(keys er, i-> errorBetti#i = er#i);
    totBetti := new MutableHashTable;
    scan(keyList, i->  totBetti#i = totalRank(multiBetti#i));
    domWeights  := new MutableHashTable;
    scan(keyList, i->  domWeights#i = maximalDominantWeights(apply(schBetti#i, l -> l_0)));
    numReps  = new MutableHashTable;
    scan(keyList, i->  numReps#i = #(schBetti#i));
    numRepsMulti := new MutableHashTable;
    scan(keyList, i->  numRepsMulti#i = sum apply(schBetti#i, j-> j_1));
    lexLeadingWeights := new MutableHashTable;
    scan(keyList, i->  lexLeadingWeights#i = if schBetti#i == {} then {} else last sort apply(schBetti#i,j-> j_0));
    BS := BSCoeffs makeBetti totBetti;
    --  Writing the files
    g = openOut ("M2OutputFiles/bettiP2_"|toString d|"_"|toString b|".m2");
    g<< "--This file computes Betti tables for P^2 for d = "|toString d|" and b = "|toString b;
    g<< endl;
    g<< "A := QQ[t_0,t_1,t_2];";
    -- need to add bit to initialize A = QQ[t_0,t_1,t_2];
    g<< endl;
    g<< "--tb stands for Total Betti numbers";
    g<< endl;
    g<< "tb"|toString d|toString b|" = ";
    g<< toExternalString unMute(totBetti);
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
    g<< toExternalString unMute(schBetti);
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
    g<< "--bs encodes the Boij-Soederberg coefficients each entry";
    g<< endl ;
    g<< "bs"|toString d|toString b|" = ";
    g<< toExternalString BS;
    g<< ";";
    g<< endl;
    g<< "end;";    
    close g;    
    )




--  This just reformats previously computed output files for P1
--  to create new cases you'd need other code.
makeP1Output =  (d,n,b)->(
    needsPackage "SchurVeronese";
    if n != 1 then error "n needs to be 1 for this command";
    mb := multiBetti(d,n,b);
    sb := schurBetti(d,n,b);
    multBetti := mb;
    schBetti := sb;
    keyList := keys multBetti;
    errBetti := new MutableHashTable;
    scan(keyList, i->  errBetti#i = 0_A);
    --scan(keys er, i-> errorBetti#i = er#i);
    --no error in P1 cases because no relevant range
    totBetti := new MutableHashTable;
    scan(keyList, i->  totBetti#i = totalRank(multBetti#i));
    domWeights  := new MutableHashTable;
    scan(keyList, i->  domWeights#i = maximalDominantWeights(apply(schBetti#i, l -> l_0)));
    numReps  = new MutableHashTable;
    scan(keyList, i->  numReps#i = #(schBetti#i));
    numRepsMulti := new MutableHashTable;
    scan(keyList, i->  numRepsMulti#i = sum apply(schBetti#i, j-> j_1));
    lexLeadingWeights := new MutableHashTable;
    scan(keyList, i->  lexLeadingWeights#i = if schBetti#i == {} then {} else last sort apply(schBetti#i,j-> j_0));
    BS := BSCoeffs makeBetti totBetti;
    --  Writing the files
    g = openOut ("M2OutputFiles/bettiP"|toString n|"_"|toString d|"_"|toString b|".m2");
    g<< "--This file computes Betti tables for P^2 for d = "|toString d|" and b = "|toString b;
    g<< endl;
    g<< "A := QQ[t_0,t_1,t_2];";
    -- need to add bit to initialize A = QQ[t_0,t_1,t_2];
    g<< endl;
    g<< "--tb stands for Total Betti numbers";
    g<< endl;
    g<< "tb"|toString d|toString b|" = ";
    g<< toExternalString unMute(totBetti);
    g<< ";";
    g<< endl; 
    g<< "--mb stands for Multigraded Betti numbers";
    g<< endl ;
    g<< "mb"|toString d|toString b|" = ";
    g<< toExternalString unMute(multBetti);
    g<< ";";
    g<< endl;
    g<< "--sb represents the betti numbers as sums of Schur functors";
    g<< endl ;
    g<< "sb"|toString d|toString b|" = ";
    g<< toExternalString unMute(schBetti);
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
    g<< toExternalString unMute(errBetti);
    g<< ";";
    g<< endl;
    g<< "--bs encodes the Boij-Soederberg coefficients each entry";
    g<< endl ;
    g<< "bs"|toString d|toString b|" = ";
    g<< toExternalString BS;
    g<< ";";
    g<< endl;
    g<< "end;";    
    close g;    
    )


--  Input:  parameters, d,n,b
--  Output: M2 file of syzygy data for O(b) on P^n embedded by O(d)
--          but only for entries outside of relevant range.
makePartialOutput =  (d,n,b)->(
    L := relevantRange(d,n,b);
    cL := apply(L,l-> (l#0-1,l#1+1));
    (multBetti,schBetti) := startBetti(d,n,b);
    keyList := keys multBetti;
    altList := toList(set(keys multBetti) - set(L) - set(cL));
    scan(L|cL, i-> multBetti#i = infinity);
    scan(L|cL, i-> schBetti#i = ({0,0,0},infinity));
    errorBetti := new MutableHashTable;
    scan(altList, i->  errorBetti#i = 0_A);
    scan(L|cL, i-> errorBetti#i = infinity);
    --errorBetti := new MutableHashTable;
    --scan(keyList, i->  errorBetti#i = 0_A);
    --scan(keys er, i-> errorBetti#i = er#i);
    totBetti := new MutableHashTable;
    scan(altList, i->  totBetti#i = totalRank(multBetti#i));
    scan(L|cL, i->  totBetti#i = infinity);
    domWeights  := new MutableHashTable;
    scan(altList, i->  domWeights#i = maximalDominantWeights(apply(schBetti#i, l -> l_0)));
    numReps  = new MutableHashTable;
    scan(altList, i->  numReps#i = #(schBetti#i));
    numRepsMulti := new MutableHashTable;
    scan(altList, i->  numRepsMulti#i = sum apply(schBetti#i, j-> j_1));
    lexLeadingWeights := new MutableHashTable;
    scan(altList, i->  lexLeadingWeights#i = if schBetti#i == {} then {} else last sort apply(schBetti#i,j-> j_0));
    --  Can't do BS coeffs for a partial table.
    --
    --  Writing the files
    g = openOut ("M2OutputFiles/bettiP2_"|toString d|"_"|toString b|".m2");
    g<< "--This file computes Betti tables for P^2 for d = "|toString d|" and b = "|toString b;
    g<< endl;
    g<< "A := QQ[t_0,t_1,t_2];";
    -- need to add bit to initialize A = QQ[t_0,t_1,t_2];
    g<< endl;
    g<< "--tb stands for Total Betti numbers";
    g<< endl;
    g<< "tb"|toString d|toString b|" = ";
    g<< toExternalString unMute(totBetti);
    g<< ";";
    g<< endl; 
    g<< "--mb stands for Multigraded Betti numbers";
    g<< endl ;
    g<< "mb"|toString d|toString b|" = ";
    g<< toExternalString unMute(multBetti);
    g<< ";";
    g<< endl;
    g<< "--sb represents the betti numbers as sums of Schur functors";
    g<< endl ;
    g<< "sb"|toString d|toString b|" = ";
    g<< toExternalString unMute(schBetti);
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


--  Input:  parameters, d,n,b
--  Output: M2 file of syzygy data for O(b) on P^n embedded by O(d)
--          but only for entries outside of relevant range.
--          And no multigraded Hilbert series (since that's a computational
--          bottleneck)
makePartialNoMultiOutput =  (d,n,b)->(
    L := relevantRange(d,n,b);
    cL := apply(L,l-> (l#0-1,l#1+1));
    schBetti := startBettiSchur(d,n,b);
    keyList := keys schBetti;
    altList := toList(set(keys schBetti) - set(L) - set(cL));
    scan(L|cL, i-> schBetti#i = ({0,0,0},infinity));
    errorBetti := new MutableHashTable;
    scan(altList, i->  errorBetti#i = 0_A);
    scan(L|cL, i-> errorBetti#i = infinity);
    --errorBetti := new MutableHashTable;
    --scan(keyList, i->  errorBetti#i = 0_A);
    --scan(keys er, i-> errorBetti#i = er#i);
    totBetti := new MutableHashTable;
    scan(altList, i->  totBetti#i = schurRank(n,schBetti#i));
    scan(L|cL, i->  totBetti#i = infinity);
    domWeights  := new MutableHashTable;
    scan(altList, i->  domWeights#i = betterDominantWeights(apply(schBetti#i, l -> l_0)));
    numReps  = new MutableHashTable;
    scan(altList, i->  numReps#i = #(schBetti#i));
    numRepsMulti := new MutableHashTable;
    scan(altList, i->  numRepsMulti#i = sum apply(schBetti#i, j-> j_1));
    lexLeadingWeights := new MutableHashTable;
    scan(altList, i->  lexLeadingWeights#i = if schBetti#i == {} then {} else last sort apply(schBetti#i,j-> j_0));
    --  Can't do BS coeffs for a partial table.
    --
    --  Writing the files
    g = openOut ("M2OutputFiles/bettiP2_"|toString d|"_"|toString b|".m2");
    g<< "--This file computes Betti tables for P^2 for d = "|toString d|" and b = "|toString b;
    g<< endl;
    g<< "A := QQ[t_0,t_1,t_2];";
    -- need to add bit to initialize A = QQ[t_0,t_1,t_2];
    g<< endl;
    g<< "--tb stands for Total Betti numbers";
    g<< endl;
    g<< "tb"|toString d|toString b|" = ";
    g<< toExternalString unMute(totBetti);
    g<< ";";
    g<< endl; 
    g<< "--sb represents the betti numbers as sums of Schur functors";
    g<< endl ;
    g<< "sb"|toString d|toString b|" = ";
    g<< toExternalString unMute(schBetti);
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




--
--
--  Appendix 1:  Claudiu's code for estimating schur Betti table assuming no cancellationsw
--
--
--

--  This is internal code to translate between the vero function
--  and the format we were using
standardize = P -> (
    if P == 0 then return {};
    K := listForm P;
    apply(K, k-> if #(k#0) == 3 then k else (k#0|apply(3-#(k#0),i-> 0), k#1))
    )




recsyz := method();
recsyz (Thing) := (el) -> min(el,0);
recsyz (RingElement) := (el) ->
(
    T := ring el;
    listForm el/((u,v)->T_u*recsyz(v))//sum
    );


--  This is Claudiu Raicu's code for computing Schur Betti tables
--  assuming no consecutive cancellations
--  It computes it for for O(r) on PP^{k-1} embedded by O(d)
--  So coordinating with our notation, (k,d,r) = (n+1,d,b)
veroP2 = (k,d,r)->(
    S := schurRing(QQ,s,k);
    T := schurRing(S,u,1);
    aS := schurRing(QQ,s,10);
    aR := symmetricRing aS;
    parn := (partitions d)/(i->toList i);
    m := matrix{apply(parn,i->i/(j->h_j)//product)};
    Mat := transpose matrix{apply(parn,i->toH s_i)};
    ikost := lift(contract(Mat,m),QQ);
    ind = select(for i from 0 to #parn-1 list if #(parn#i)<=k then i,j->j=!=null);
    pars = select(parn,i->#i<=k);
    aQ = schurRing(QQ,q,k);
    subkost = ikost^ind_ind;
    qmat = matrix{pars/(i->q_i)};
    mpols = flatten entries (qmat * subkost);
    mpols = drop(mpols,1);
    aT := schurRing(aQ,u,1);
    factors = for pol in mpols list(
     	di = lift(dim(pol),ZZ);
     	sq = splice{k:di*d//k};
     	rez = 1_aT + (-1)^di*u_{di}*q_(sq);
     	rez = rez + sum for i from 1 to di//2 list 
     	(
	    ePow = exteriorPower(i,pol);
	    if i!=di-i then 
	    (
	       	transp = listForm(ePow)/(i->i#1*q_(reverse select(sq-(i#0|splice{k-#(i#0):0}),j->j!=0)))//sum;
	       	(-1)^i*u_{i}*ePow + (-1)^(di-i)*u_{di-i}*transp
	       	)
	    else (-1)^i*u_{i}*ePow
	    );
     	rez
     	);
    genfun = product factors;
    listForm genfun;
    wedpowers = (listForm genfun)/(i->lift(last i,aQ))//reverse;
    wedpowers = wedpowers/(i->toS(i,S));
    p1 = 1 + sum for i from 1 to #wedpowers-1 list wedpowers#i * T_{i};
    A := QQ[x];
    cyclo = x^d-1;
    for i from 1 to d-1 do if d%i == 0 then cyclo = cyclo // gcd(cyclo,x^i-1);
    B := A/ideal(cyclo);
    SS := schurRing(B,ss,k);
    TS := schurRing(SS,ts,1);
    pol = sum for i from 0 to d-1 list(
     	x^(i*(d-r)) * (product for j from 0 to d-1 list
     	    if j == i then 1 else
	    1 + sum for p from 1 to k list (-1)^p*x^(j*p)*SS_(splice{p:1})*TS_p)
     	);
    pol = pol/d;
    pol = (listForm pol)/(i->(try(TS_(i#0#0//d)) else 1_TS)*(last i))//sum;
    p2 = toS(pol,T);
    reso = p1*p2;
    a := binomial(d+2,2)-3;
    b := 2;
    koz = new MutableList from for i from 0 to a list new MutableList from {0,0,0};
    rev = reverse listForm reso;
    maxq0 = binomial(r+2,2)-1;
    minq2 = binomial(d+2,2)-binomial(d-r-1,2)-2;
    for syzy in rev do(
     	deg := sum (first syzy);
     	rep := (-1)^deg * (last syzy);
     	reppos := rep - recsyz(rep);
     	repneg := reppos - rep;
     	if repneg != 0 then koz#(deg-1)#1 = repneg;
     	if deg <= maxq0 then koz#deg#0 = reppos
     	else if deg >= minq2 then koz#(deg-2)#2 = reppos
     	else if reppos != 0 then error"Something's wrong...";
     	);
    koz
    );








--
--
--
--    Appendix 2:  helper functions for Boij-Soederberg coefficients
--                 some were written by Courtney Gibbons
--
--

needsPackage "BoijSoederberg";
--The "HK" suffix indicates a different rescaling
--of the choice of pureBetti diagram for each degree
--sequence.  We've chosen the good "uniform" one
--where the ith entry is 1/(\prod_i\ne j |d_i-d_j|).
isStrictlyIncreasing = method();
isStrictlyIncreasing List := L -> (
     t:=true;
     for i from 0 to #L-2 do t=(t and (L_i<L_(i+1)));
     t)

--Input: Degs must be a strictly increasing list of positive integers
--Output: List of ranks of the minimal integral betti sequence that satisfy the
--"Peskine-Szpiro" equations
pureBettiHK = method()
pureBettiHK List := (Degs) -> (
     c := # Degs;
     p := 1;
     for i from 1 to c-1 do (
	  if Degs#i <= Degs#(i-1) then error "--pureBetti: expected an increasing list of integers";
	  for j from 0 to i-1 do p=p*(Degs_i-Degs_j)
	  );
     D := for i from 0 to c-1 list (-1)^i * product(i, j->Degs_j-Degs_i) * product(i+1..c-1, j->Degs_j-Degs_i);
--     Bettis := for i from 0 to c-1 list (p/D_i);
     Bettis := for i from 0 to c-1 list (1/D_i);
--     Bettis = Bettis / gcd Bettis;
     apply(Bettis, x -> lift(x,QQ)))

pureBettiDiagramHK = method()
pureBettiDiagramHK List := (degs) -> (
     B := pureBettiHK degs;
     new BettiTally from apply(#degs, i -> (i, {degs#i}, degs#i) => B#i)
     )

decompose1 = method();
decompose1 BettiTally := B -> (
     L:=lowestDegrees B;
     if not isStrictlyIncreasing L then print "NOT IN THIS SIMPLEX OF PURE BETTI DIAGRAMS";
     C:=pureBettiDiagramHK L;
     ratio:=min apply(#L, i->(B#(i,{L_i}, L_i))/(C#(i,{L_i},L_i)));
     (C,ratio,merge(B,C, (i,j)->i-ratio*j))
     )



--
--  
--  
decomposeHK = method();
decomposeHK BettiTally := B-> (
     Components:={};
     B1:= new MutableHashTable from B;
     while min values B1 >= 0 and max values B1 > 0 do (
	  X:=decompose1(new BettiTally from B1);
	  B1=new MutableHashTable from X_2;
	  --change the type of the values in X_0 to ZZ
	  Y:=new BettiTally from apply(pairs X_0, i->{first i, lift(last i, QQ)});
	  Components = append(Components, hold(X_1) * Y));
     sum Components
     )

---Input: Betti table
---Output: List of the Boij-Soederberg coefficients
---Caveat: for simplicity we omit the degree sequences.
--         so this is only meaningful if you know the corresponding
--         degree sequences.
BSCoeffs = method();
BSCoeffs BettiTally := B-> (
     Components := {};
     Coeffs :={};
     B1:= new MutableHashTable from B;
     while min values B1 >= 0 and max values B1 > 0 do (
	  X:=decompose1(new BettiTally from B1);
	  B1=new MutableHashTable from X_2;
	  --change the type of the values in X_0 to ZZ
	  Y:=new BettiTally from apply(pairs X_0, i->{first i, lift(last i, QQ)});
	  Components = append(Components, hold(X_1) * Y);
      	  Coeffs = append(Coeffs,X_1);
	  );
     Coeffs
     )

---Input: Betti table
---Output: List of the Boij-Soederberg coefficients, rescaled by a
--         power of 10 so that the largest coefficient is between 1 and 10.
rescaledBSCoeffs = method();
rescaledBSCoeffs BettiTally := B->(
    L := BSCoeffs B;
    m := max apply(L,l-> floor log(10,l));
    apply(L,l-> (l/10^m)_RR)
    )





end;



---
restart
load "step3.m2"
--- Test
--dualHilb(5,2,0,2)
--dualHilb(5,2,4,1)
--dualHilb(6,2,0,2)
--dualHilb(6,2,1,2)
dualHilb(5,2,0)
dualHilb(5,2,4)
dualHilb(6,2,0)
dualHilb(6,2,1)
makeOutputFiles(5,2,3)
makeOutputFiles(5,2,4)
makeOutputFiles(6,2,0)
makeOutputFiles(6,2,1)
makeBetti = H ->(
    new BettiTally from apply(keys H, h-> (h_0,{h_0+h_1},h_0+h_1)=> H#h)
    )
BSCoeffs makeBetti tb41
break
peek er50
