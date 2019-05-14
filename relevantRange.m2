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
----- DESCRIPTION: This function checks whether the pushforward of
----- O(B) by O(D) on F_{a} is Cohen-Macaulay. Here F_a is the a-th
----- Hirzebruch surface.
-----
----- WARNINGS: Everything is hardcoded to only handle a = 0.
--------------------------------------------------------------------
--------------------------------------------------------------------
isCM = (a,B,D)->(
    if a != 0 then return "error: a != 0";
    if ((D#0)/(D#1)*(B#1)-B#0 < 2) and ((D#1)/(D#0)*(B#0)-B#1 < 2) then true
    else false
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
lowerBound = (q,a,B,D)->(
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
lowerBound(0,a,B,D)
lowerBound(1,a,B,D)
lowerBound(2,a,B,D)


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
upperBound = (q,a,B,D)->(
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
upperBound(0,a,B,D)
upperBound(1,a,B,D)
upperBound(2,a,B,D)


apply(3,q->( (lowerBound(q,0,B,D)..upperBound(q,0,B,D))))

relevantRange = (a,B,D) ->(
    if a != 0 then return "error: a != 0";
    apply(3,q->(
	    toList apply((lowerBound(q,a,B,D)..upperBound(q,0,B,D)),i->(i,q))
	    ))
    

    else (
	overlap := flatten apply(ceiling(n+1-(n+b)/d),q->(
		toList apply((lowerBound(d,n,b,q+1)+1..upperBound(d,n,b,q)),i->(i,q))
		));
	ends := flatten apply(ceiling(n+1-(n+b)/d),q->({(lowerBound(d,n,b,q)-1,q),(upperBound(d,n,b,q)+1,q)}
		));
	overlap|ends	
    )
)