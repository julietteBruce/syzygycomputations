
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
    h := hilbSeries22(m);
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


------------------------
------------------------
--Part 2: naiveBetti code
------------------------
------------------------


--  The analogue of Algorithm 3.3 for P1xP1.
--Input B,D
--Output:

naiveBetti = (B,D) ->(
    S = QQ[x_0,x_1,y_0,y_1,Degrees => {{1,0},{1,0},{0,1},{0,1}}];
    --  (B,D) = ({1,0},{1,2})
    --  N is supposed to be the biggest ZZ^2 degree
    --  we would need to worry about.
    --  CAVEAT: This comes from a regularity bound, but requires various unspecified assumptions on (B,D).  Needs to be examined in detail later.
    N =(D_0+1)*(D_1+1)*{D_0,D_1}+{B_0,B_1};
    Lpre0 = flatten apply(toList(0..(N_0)),i->(
       	    apply(toList(0..(N_0)-i),j->(
		    {i,j}))));
    Lpre1 = flatten apply(toList(0..(N_1)),i->(
       	    apply(toList(0..(N_1)-i),j->(
		    {i,j}))));
    L = {};	   
    scan(Lpre0, a->(
	    scan(Lpre1, a'->(
		    K = {(a_0+a_1-B_0),(a'_0+a'_1-B_1)};
		    alpha = K_0//D_0;
		    if alpha*D == K then L = L|{flatten {a,a'}};
		    ))));
    R = QQ[t_0,t_1,t_2,t_3];
    C = sum(L,l-> 1*t_0^(l_0)*t_1^(l_1)*t_2^(l_2)*t_3^(l_3));
    Bexps0 = flatten apply(toList(0..D_0),i->({{i,D_0-i}}));
    Bexps1 = flatten apply(toList(0..D_1),i->({{i,D_1-i}}));
    Bexps = flatten apply(Bexps0,a->apply(Bexps1,b-> flatten{a,b}));
    Bpoly = product apply(Bexps,l-> 1 - 1*t_0^(l_0)*t_1^(l_1)*t_2^(l_2)*t_3^(l_3));
    --  DONT REUSE B
    m0 = ideal(t_0,t_1);
    m1 = ideal(t_2,t_3);
    A = Bpoly*C % (m0^(N_0+1)*m1^(N_1+1));
    --print A;
    sub(A,{t_1=>t_0,t_3=>t_2});
    Asg = sub(A,{t_1=>t_0,t_3=>t_0,t_2=>t_0});
    print Asg;
    topGuy = (N_0+N_1)//(D_0+D_1);
    BBkeys0 = apply(toList(0..topGuy),i-> (i,{i},i)=>0);
    BBkeys1 = apply(toList(1..topGuy),i-> (i,{i+1},i+1)=>0);
    BBkeys2 = apply(toList(2..topGuy),i-> (i,{i+2},i+2)=>0);
    BBkeys3 = apply(toList(3..topGuy),i-> (i,{i+3},i+3)=>0);
    emptyTable = BBkeys0|BBkeys1|BBkeys2|BBkeys3;
    T = new MutableHashTable from emptyTable;
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
    L = apply(toList(keys(T)), i-> i=>T#i);
    Betti = new BettiTally from L;
    Betti
    );

end
