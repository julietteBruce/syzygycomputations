needsPackage"SchurRings"
recsyz := method();
recsyz (Thing) := (el) -> min(el,0);
recsyz (RingElement) := (el) ->
(
    T := ring el;
    listForm el/((u,v)->T_u*recsyz(v))//sum
    );

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
end

----------------------------------
----------------------------------
----------------------------------
----------------------------------
----------------------------------
----------------------------------

restart
d = 5 --d-th Veronese embedding
k = 3 --dim(V) = k (should be 3), called n in Ein-Lazarsfeld
r = 0 --twist of structure sheaf, called b in Ein-Lazarsfeld
time load"veroP2.m2"
time mat = matrix(toList koz / (i -> (toList i / (j-> if j != 0 then dim j else 0))))
koz#3#1
value toString last s_{23,8,4}
new List from koz#6
