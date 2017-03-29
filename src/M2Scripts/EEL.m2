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

load "Rep_Theory.m2"

--returns false if the lengths disagree!
--return true if f is true for all elements in the list
compareList = (f,a,b) -> (
    if #a!=#b then return false;
    for i from 0 to (#a-1) do (
        if not f(a#i,b#i) then return false
        );
    true
    )

--find all wedges of m things with everything less than d not divisible by f. For this to work properly,
--currently works only for n=2
dominantEELWedges = (f,m,d) -> (
    if m<0 then return {};
    if m==0 then return {{}};
    isGen := lst -> all(lst, x-> x<d);
    dividesF := lst -> compareList((x,y) -> x<=y,lst,f);
    inAnnF := lst -> any(lst+f,x-> x>=d);
    --successors in the dominance order that are in the annihilator of f in R/(x_1^d,...)
    --currently we assume that if a child is not in the annihlator then that entire branch can be ignored
    succ := lst -> select(apply({{ -1,1,0 }, {0,-1,1}},offset -> lst + offset),x->all(x,y-> y>=0) and inAnnF x);
    --successors ignoring things divisible by f and ignoring those that are redundent
    --do I need guarentees about ordering?
    filteredSucc := (lst) -> (
        p := partition(x -> dividesF x or (not isGen x),succ(lst));
        good := if p#?false then p#false else {};
        bad := if p#?true then p#true else {};
        good|select(unique flatten apply(bad,filteredSucc),lst -> not cmpL(lst,good)));

    --do depth first search on the poset to find all of the wedges.
    --current is the current set of wedges,
    --visible is the current layer of the poset under consideration
    --usableIndex tells us where the first usuable element
    --  i.e. the first element not prohibited by uniqueness
    dfs := (current,visible,usableIndex) -> (
        if #current == m then return {current};
        flatten (for idx from usableIndex to #visible - 1 list (
            newElement := visible#idx;
            rest := drop(visible,{idx,idx});
            newVisible := select(filteredSucc(newElement),lst -> not cmpL(lst,rest) );
            dfs(current | {newElement},unique (rest|newVisible),idx)
            ))
        );
    start = {d-1,1,0};
    if dividesF start or (not isGen start)
    then dfs({},filteredSucc(start),0)
    else dfs({},{start},0)
    )

dominantEELWeights = (d,p,q,b,n) -> (
    x = symbol x;
    S := QQ[x_0..x_n, MonomialOrder => Lex];
    R := S/ideal apply(n+1,i-> x_i^d);
    --f := product apply(q,i-> x_i^(d-1))* x_(q)^(q+b);
    if basis(q*d+b,R) == 0 then return {};
    f = (flatten entries basis(q*d+b,R))_0;
    I = ann(sub(f,R));
    L := flatten entries gens image basis(d,I);
    DivE = apply(select(L,l -> f % l==0),g -> (exponents g)#0);
    E = (exponents f)#0 + if #DivE==0 then {0,0,0} else sum DivE; --this won't work for n!=2
    betterDominantWeights reverse sort unique apply(dominantEELWedges((exponents f)#0,p-#DivE,d),k->  sum (k|{{0,0,0}}) +E))
end;

restart
load "Rep_Theory.m2"
load "EEL.m2"
for i from 0 to 30 do (
    for j from 0 to 2 do (
        if i==0 and j==0 then continue;
        L := dominantEELWeights(7,i,j,0,2);
        M := betterDominantWeights dominantEELWeights(7,i,j,0,2);
        if L!=M
        then (print("weights aren't dominant " | toString(i) | "," | toString (j));
              print(toString(L) | " " | toString(M)));
        )
    )
--for 6,2 we need to compute 2,3,4,5
--for 6,1 we need to compute 1,2
EELWeights(6,6,1,1,2)
dominantEELWeights(6,6,1,1,2)
betterDominantWeights EELWeights(6,6,1,1,2)

L = betterDominantWeights EELWeights(6,5,0,3,2);
L
K = unique apply(flatten apply(L,l-> {l+{1,-1,0},l+{1,0,-1},l+{0,1,-1}}), k-> reverse sort k);
--K is the list for the (5,0) entry where d=6 and b=3
--We expect all multigraded Betti numbers from K to be 0.
#K
K
--Confirmed by Jay in several seconds!!

dominantEELWeights(6,5,0,3,2)

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


weightsToCheck(9,8,1,6)

betterDominantWeights EELWeights(9,8,1,6,2)
dominantEELWeights(9,8,1,6,2)

printWeightsToCheck(7,6,0,4,"test.txt")

--This case is too big with EELWeights, but works with dominantEELWeights
EELWeights(7,19,1,4,2)
dominantEELWeights(7,19,1,4,2)
