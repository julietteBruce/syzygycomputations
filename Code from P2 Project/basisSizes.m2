
--basisSizes returns a tally mapping multidegrees to the basis sizes
basisSizes := (n,d,b,p,q) -> (
    S := QQ[x_0..x_n];
    x = symbol x;
    L1 := flatten entries basis(d,S);
    L2 := flatten entries basis(q+b,S);
    L3 :=
    flatten for l1 in subsets(L1,p) list (
        for l2 in L2 list product(l1)*l2);
    print tally(L3);
    )
basisSizes(2,5,0,3,0)
