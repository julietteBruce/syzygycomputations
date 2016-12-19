EELWeights = (d,p,q,b,n)->(
    x = symbol x;
    S := QQ[x_0..x_n, MonomialOrder => Lex];
    R := S/ideal apply(n+1,i-> x_i^d);
    --f := product apply(q,i-> x_i^(d-1))* x_(q)^(q+b);
    if basis(q*d+b,R) == 0 then return {};
    f = (flatten entries basis(q*d+b,R))_0;    
    I = ann(sub(f,R));
    L := flatten entries gens image basis(d,I);
    if #L < p then return {};
    --This next line is very hacky
    if L == {} then (if p == 0 then return {{0,0,0}} else return {});
    Div = {};
    NotDiv = {};
    scan(L, l -> if f % l == 0 then Div = Div|{l} else NotDiv = NotDiv|{l});
    DivE = apply(Div, l->  (exponents sub(l,S))_0);
    if #DivE > p then return {};
    if #DivE == p then(
	if p == 0 then return {(exponents f)_0} else return  {((exponents f)_0+ sum DivE)};
	);
    NotDivE = apply(NotDiv, l->  (exponents sub(l,S))_0);
    K := subsets(NotDivE,p-#DivE);
    --if NotDivE == {} then NotDivE = {{0,0,0}};
    if DivE == {} then DivE = {{0,0,0}};
    E = (exponents f)_0 + sum DivE;
--  some hacky stuff to try and deal w boundary cases but causes problems.
    if K == {} then K =  {{{0,0,0}}};
    if K == {{}} then K = {{{0,0,0}}};
    reverse sort unique apply(K,k->  sum k + E)
    --we screwed up the boundary cases
    )  

end;

load "Rep_Theory.m2"
load "EEL.m2"
L = betterDominantWeights EELWeights(6,5,0,3,2);
K = unique apply(flatten apply(L,l-> {l+{1,-1,0},l+{1,0,-1},l+{0,1,-1}}), k-> reverse sort k);
--K is the list for the (5,0) entry where d=6 and b=3
--We expect all multigraded Betti numbers from K to be 0.
#K
K
--Confirmed by Jay in several seconds!!

L = betterDominantWeights EELWeights(7,6,0,4,2);
K = unique apply(flatten apply(L,l-> {l+{1,-1,0},l+{1,0,-1},l+{0,1,-1}}), k-> reverse sort k);
#K
K
--K is the list for the (6,0) entry where d=7 and b=4
--We expect all multigraded Betti numbers from K to be 0.

L = betterDominantWeights EELWeights(7,5,0,4,2);
K = unique apply(flatten apply(L,l-> {l+{1,-1,0},l+{1,0,-1},l+{0,1,-1}}), k-> reverse sort k);
#K
K

weightsToCheck = (d,p,q,b)->(
    L := betterDominantWeights EELWeights(d,p,q,b,2);
    K := unique apply(flatten apply(L,l-> {l+{1,-1,0},l+{1,0,-1},l+{0,1,-1}}), k-> reverse sort k);
    K
    )

weightsToCheck(7,6,0,4)
apply(weightsToCheck(7,7,0,4),reverse)
apply(weightsToCheck(7,8,0,4),reverse)



printWeightsToCheck(7,6,0,4,"test.txt")
