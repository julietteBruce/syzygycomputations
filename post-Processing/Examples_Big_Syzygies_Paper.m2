needsPackage "SchurVeronese"
load "step3.m2"
viewHelp SchurVeronese
tex totalBetti2(5,2,3)
apply(bsCoeffs(5,2,3),n-> (10^(-16)*n)_RR)
(schurBetti(5,2,3))#(9,0)


(d,n,b)= (5,2,0);
nr = numRepsBetti(d,n,b);
L ={};
scan(keys nr,k-> if k#1==1 then L = L|{k});
L = sort L
apply(L,k-> (nr#k))
dw = dominantWeightsBetti(d,n,b);
apply(L,k-> #(dw#k))

findRedundants = (d,n,b)->(
    sb = schurBetti(d,n,b);
    R = relevantRange(d,n,b);
    if #R == 0 then return {};
    r = R#0;
    r' = (r#0-1,r#1+1);
    justWeights := apply(sb#(r'),i-> i#0);
    OV = {};
    scan(sb#r,i->(
	    if isSubset({i#0},justWeights) then OV = OV|{i#0}
	    ));
    OV
    );
findRedundants(5,2,0)
break
