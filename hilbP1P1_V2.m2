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

--------All of the below was us checking code:
topGuy=10
BBkeys0 = apply(toList(0..topGuy),i-> (i,{i},i)=>0);
    BBkeys1 = apply(toList(1..topGuy),i-> (i,{i+1},i+1)=>0);
    BBkeys2 = apply(toList(2..topGuy),i-> (i,{i+2},i+2)=>0);
    BBkeys3 = apply(toList(3..topGuy),i-> (i,{i+3},i+3)=>0);
    emptyTable = BBkeys0|BBkeys1|BBkeys2|BBkeys3;
    Betti = new BettiTally from emptyTable;
    T = new MutableHashTable from emptyTable;
    Betti
    emptyTable
    Asg = -30*t_0^9 - 2*t_0^3 + 2*t_0
    col = 0
    row = 0
    D_0 = 1
    D_1 = 1
    B_0 = 1
    B_1 = 0
    for k from 0 to 4 do(
	coef = sub(((-1)^(col)*coefficient(t_0^((D_0+D_1)*(k)+(B_0+B_1)),Asg)),ZZ);
	print coef;
	if coef >= 0 then ((T#(col,{col+row},col+row) = coef); 
	    (col= col+1);
	    print col);
	if coef < 0 then (
	    row = row +1;
	    print row;
	    print col;
	    T#(col-1,{col+row-1},col+row-1) = -coef);
	)
    L = apply(toList(keys(T)), i-> i=>T#i)
    new BettiTally from L
T#(0,{0},0)	
T#(1,{1},1)  
T#(2,{2},2)
T#(3,{3},3)
T#(3,{4},4)
T#(4,{4},4)    

--Extra code   
    che = true;
    row0 = apply(BBkeys0, K-> 
	if row = 0 do(
	    if ((-1)^(K_0)*coefficient(t_0^((D_0+D_1)*(K_2)+(B_0+B_1)),Asg)) >= 0 then K=>((-1)^(K_0)*coefficient(t_0^((D_0+D_1)*(K_2)+(B_0+B_1)),Asg))
	    if 
	else K=>0
	);
    row1 = apply(BBkeys1, K-> 
	if ((-1)^(K_0)*coefficient(t_0^((D_0+D_1)*(K_2)+(B_0+B_1)),Asg)) >= 0 then K=>((-1)^(K_0)*coefficient(t_0^((D_0+D_1)*(K_2)+(B_0+B_1)),Asg))
    	else K=>0
	);
    row2 = apply(BBkeys2, K-> 
	if ((-1)^(K_0)*coefficient(t_0^((D_0+D_1)*(K_2)+(B_0+B_1)),Asg)) >= 0 then K=>((-1)^(K_0)*coefficient(t_0^((D_0+D_1)*(K_2)+(B_0+B_1)),Asg))
    	else K=>0
	);
    row3 = apply(BBkeys3, K-> 
	if ((-1)^(K_0)*coefficient(t_0^((D_0+D_1)*(K_2)+(B_0+B_1)),Asg)) >= 0 then K=>((-1)^(K_0)*coefficient(t_0^((D_0+D_1)*(K_2)+(B_0+B_1)),Asg))
    	else K=>0
	);--Can create function for the if statement, this is super inefficient
    
    --print row0;
    --print row1;
    --print row2;
    --row0 = {(0,{0},0)=>1}; 


---------TESTING:
restart
load "hilbP1P1_V2.m2"
--Segre Embedding, this case is alright
naiveBetti((({0,0},{1,1}))

--Errors when we make B_0 != {0,0}
naiveBetti(({0,0},{2,3}))

--
naiveBetti(({1,1},{4,1}))

--
naiveBetti(({1,1},{3,2}))

