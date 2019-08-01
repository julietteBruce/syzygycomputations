restart
needsPackage "SchurVeronese"
(schurBetti(5,2,3))#(6,0)
#oo

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



SB = new MutableHashTable from (schurBetti(5,2,3));
KEYS = keys (sb)
apply(KEYS, k-> sb#k = fixDuplicateSchurs(sb#k));

apply({5,6,7,8,9},k->(
L1 = rsort apply((schurBetti(5,2,3))#(k,0), i->i#0);
{(k,0),unique delete(,apply(L1, i->(
            c = 0;
            c = #select(L1,j->(j==i));
            if c>1 then i
            )))}
))
 
 
 
 
tally (schurBetti(5,2,3))#(6,0)
needsPackage "SchurVeronese"



#fixDuplicateSchurs(L)



OK | (new MutableList from newGuys)
OK
viewHelp MutableList

(schurBetti(5,2,3))#(6,0)
SB#{15,14,4}
l = (keys SB)#0
SB#
scan(keys SB, i-> scan(keys SB, j-> if 



apply({4,5,6,7,8},k->(
L1 = rsort apply((schurBetti(5,2,3))#(k,1), i->i#0);
{(k,1),unique delete(,apply(L1, i->(
            c = 0;
            c = #select(L1,j->(j==i));
            if c>1 then i
            )))}
))
 