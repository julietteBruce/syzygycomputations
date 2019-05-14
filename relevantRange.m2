--------------------------------------------------------------------
--------------------------------------------------------------------
----- INPUT: (q,B) 
-----
----- OUPUT: an integer
-----
----- DESCRIPTION: This function computes h^0(F_a, O(B)).
-----
----- WARNINGS: Everything is hardcoded to only handle a = 0.
--------------------------------------------------------------------
--------------------------------------------------------------------
h0 = (a,B)->(
    if a != 0 then return "error: a != 0";
    product(B+{1,1})
    )

-------------------------------------------------------------------
--------------------------------------------------------------------
----- INPUT: (q,B,D) 
-----
----- OUPUT: a boolean value
-----
----- DESCRIPTION: This function computes h^0(F_a, O(B)).
-----
----- WARNINGS: Everything is hardcoded to only handle a = 0.
--------------------------------------------------------------------
--------------------------------------------------------------------
isCM = (a,B)->(
    if a != 0 then return "error: a != 0";
    product(B+{1,1})
    )


--------------------------------------------------------------------
--------------------------------------------------------------------
----- INPUT: (q,a,D,B) 
-----
----- OUPUT: an integer
-----
----- DESCRIPTION: This function gives the largest N such that we
----- are certain that the Koszul cohomology groups 
----- K_{p,q}(F_a, O(B), O(D))=0 for all p<N. Here F_a is the a-th
----- Hirzebruch surface. Note this lower bound is not necesarily sharp. 
-----
----- CAVEATS: This function will at time assume certain conditions
----- on the relationship between B and D in order to apply results 
----- from the literature. 
-----
----- WARNINGS: Everything is hardcoded to only handle a = 0.
--------------------------------------------------------------------
--------------------------------------------------------------------
lowerBound = (q,a,D,B)->(
    if a != 0 then return "error: a != 0";
    -- for q = 0 see Prop 5.1 in EL
    if q == 0 then (0) 
    else if q == 1 then (1)
    -- for q = 2 apply duality (see Prop 3.5 in EL)
    else if q == 2 then (
	if ((B#0+D#0)< -1 and (B#0+D#0)>0) or ((B#0+D#0)>0 and (B#1+D#1)< -1) then return "error: intermediate cohomology"
	else h0(a,{-2,-2}-B+D)
	)
    )

a = 0
D = {4,5}
B = {2,3}
lowerBound(0,a,D,B)
lowerBound(1,a,D,B)
lowerBound(2,a,D,B)


--------------------------------------------------------------------
--------------------------------------------------------------------
----- INPUT: (q,a,D,B) 
-----
----- OUPUT: an integer
-----
----- DESCRIPTION: This function gives the smallest N such that we
----- are certain that the Koszul cohomology groups 
----- K_{p,q}(F_a, O(B), O(D))=0 for all p>N. Here F_a is the a-th
s----- Hirzebruch surface. Note this lower bound is not necesarily sharp. 
-----
----- CAVEATS: This function will at time assume certain conditions
----- on the relationship between B and D in order to apply results 
----- from the literature. 
-----
----- WARNINGS: Everything is hardcoded to only handle a = 0.
--------------------------------------------------------------------
--------------------------------------------------------------------
upperBound = (q,a,D,B)->(
    if a != 0 then return "error: a != 0";
    -- for q = 0 see Prop 5.1 in EL
    if q == 0 then ( 
	if (B#0 > D#0) or (B#1 > D#1) then return "error: B>D"
	else h0(a,B)
	)  
    else if q == 1 then (
	h0(a,D)
	)
    -- for q = 2 apply duality (see Prop 3.5 in EL)
    else if q == 2 then (
	h0(a,D)
	)
    )

a = 0
D = {4,5}
B = {2,3}
upperBound(0,a,D,B)
upperBound(1,a,D,B)
upperBound(2,a,D,B)


--- Finds the relevant range based upon EEL non-vanishing.
--- For n>0 this does not quite handle the edge cases correctly. 
--
--- Daniel: Do you mean n>2 not n>0???
--
--- e.g. for (3,5,0) this includes (1,1) or (3,5,3) includes (-1,0)
relevantRange = (a,B,D) ->(
    if a != 0 then return "error: a != 0";
    apply(3,q->(
	    toList apply(lowerBound(
	    ))
    
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