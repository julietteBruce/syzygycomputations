
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

--  The analogue of Algorithm 3.3 for P1xP1.
--Input: B,D
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
    scan(Lpre0, a->(
	    scan(Lpre1, a'->(
		    K = {(a_0+a_1-B_0),(a'_0+a'_1-B_1)};
		    alpha = K_0//D_0;
		    if alpha*D == K then L = L|{flatten {a,a'}};
		    ))));
    C := sum(L,l-> 1*t_0^(l_0)*t_1^(l_1)*t_2^(l_2)*t_3^(l_3));
    Bexps0 := flatten apply(toList(0..D_0),i->({{i,D_0-i}}));
    Bexps1 := flatten apply(toList(0..D_1),i->({{i,D_1-i}}));
    Bexps := flatten apply(Bexps0,a->apply(Bexps1,b-> flatten{a,b}));
    Bpoly := product apply(Bexps,l-> 1 - 1*t_0^(l_0)*t_1^(l_1)*t_2^(l_2)*t_3^(l_3));
    m0 := ideal(t_0,t_1);
    m1 := ideal(t_2,t_3);
    A := Bpoly*C % (m0^(N_0+1)*m1^(N_1+1));
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
    for k from 0 to (((degree A)_0 - (B_0+B_1))//(D_0+D_1)) do(
	coef = (-1)^(col)*degPart((D_0+D_1)*(k)+(B_0+B_1),A);
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
end;


--- WARNING THIS IS HARDCODED FOR B={0,0},AND 
--- ASSUMES WE READ IN THE FIRST ROW.
--- THE PATH IS DEFINITELY WRONG
fixMultiBetti = (B,D)->(
    H = new MutableHashTable from naiveMultiBett(B,D);
    L1 = delete(,apply(relevantRange(0,B,D),i->(if i#1 == 1 then i)));
    scan(L1,i->(
	    f = value get "ranksToBetti/mgBettiData/bettiSeries_"|toString(i#0)|"_"|toString(i#1);
	    fOLD = H#i;
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
    toString(B#0)|"_"|toString(B#1)|"_"|toString(D#0)|"_"|toString(D#1)"
    )

---- H SHOULD BE THE OUTPUT OF fixMultiBetti
makeOutputFiles =  (B,D,H)->(
    --  Writing the files
    g = openOut ("M2OutputFiles/bettiF0_"|toString d|"_"|toString b|".m2");
    --g<< "--This file computes Betti tables for P^2 for d = "|toString d|" and b = "|toString b;
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




restart


load "hilbP1P1_V2.m2"
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

--Setup
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
--Output:  a hash table BpH where BpH#k is the multigraded part of Bpoly of total degree k*D.
buildBPolyHash = (D,B)->( 
    Bexps0 := flatten apply(toList(0..D_0),i->({{i,D_0-i}}));
    Bexps1 := flatten apply(toList(0..D_1),i->({{i,D_1-i}}));
    Bexps := flatten apply(Bexps0,a->apply(Bexps1,b-> flatten{a,b}));
    Bpoly := product apply(Bexps,l-> 1 - 1*t_0^(l_0)*t_1^(l_1)*t_2^(l_2)*t_3^(l_3));
    BpHM = new MutableHashTable;
    scan(topGuy+1,k-> BpHM#k=0_A);
    --tally apply(terms Bpoly, m-> (degree(m)//(D_0+D_1)) )
    scan(terms Bpoly, m->(
	    BpHM#((degree m)_0//(D_0+D_1)) = BpHM#((degree m)_0//(D_0+D_1)) + m));
    new HashTable from BpHM
    )

--Input: D,B as usual
--Output:  a hashtable AH
buildAHash = (D,B)->(
    BpH := buildBPolyHash(D,B);
    N := (D_0+1)*(D_1+1)*{D_0,D_1}+{B_0,B_1};
    topGuy := (N_0+N_1)//(D_0+D_1);
    CH := hashTable apply(topGuy+1,k->(
	    k=>sum(ZZ4degs(k,D,B),l->(
		    t_0^(l_0)*t_1^(l_1)*t_2^(l_2)*t_3^(l_3)))));
    mm := (ideal(t_0,t_1))^(N_0+1)+ (ideal(t_2,t_3))^(N_1+1); 
    hashTable apply(topGuy+1,k->(
	k => ((sum apply(k,i->(
	(BpH#i)*(CH#(k-i))))%mm))
	))  
    );

end
restart
load "hilbP1P1_V2.m2"
(D,B) = ({2,4},{0,0})
time AH = buildAHash(D,B);

