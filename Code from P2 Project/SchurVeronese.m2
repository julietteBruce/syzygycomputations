-- -*- coding: utf-8 -*-
--------------------------------------------------------------------------------
-- Copyright 2017  Juliette Bruce, Daniel Erman, Jay yang
--
-- This program is free software: you can redistribute it and/or modify it under
-- the terms of the GNU General Public License as published by the Free Software
-- Foundation, either version 3 of the License, or (at your option) any later
-- version.
--
-- This program is distributed in the hope that it will be useful, but WITHOUT
-- ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
-- FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
-- details.
--
-- You should have received a copy of the GNU General Public License along with
-- this program.  If not, see <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------------------
--------------------------------------------------------------------------------
-- PURPOSE : Data for Veronese embeddings of projective space
--
--
-- PROGRAMMERS : Juliette Bruce, Daniel Erman, Steve Goldstein, Jay yang
--
--
-- UPDATE HISTORY #0 - 
--
--
-- UPDATE HISTORY #1 - April 2019 - Juliette Bruce: Began trying to prepare
-- package for eventual publication. Adding tests, comments, documentation,
-- cleaning up code, etc.
--
--
-- TO DO LIST : 
-- #0 - Is Steve going to be an author on the package or JSAG paper?
-- #2 - Is the globally defined ring A really necesary, and will it be allowed
-- in a published M2 package?
-- #3 - Delete sandbox.
-- #4 - Add documentation for all of the exported functions.
-- #5 - Add documentation for the package itself.
--------------------------------------------------------------------------------



newPackage("SchurVeronese",
    Version => "1.1",
    Date => "13 June 2017",
    Headline => "Data for Veronese embeddings of projective space",
    Authors => {
        {Name => "Juliette Bruce",           Email => "jebruce2@wisc.edu",       HomePage => "https://juliettebruce.github.io"},
        {Name => "Daniel Erman",             Email => "derman@math.wisc.edu",    HomePage => "http://www.math.wisc.edu/~derman/"},	     
        {Name => "Steve Goldstein",          Email => "",   HomePage => ""},
	{Name => "Jay Yang",                 Email => "yangjay@math.wisc.edu",   HomePage => "http://www.math.wisc.edu/~yangjay/"}
	},
  DebuggingMode => true,
  AuxiliaryFiles => true
  )

export {
  "makeBettiTally", -- docs
  "multiBetti", -- docs
  "schurBetti", -- docs
  "totalBetti", -- docs
  "totalBetti2", -- docs
  "errorBetti", -- docs
  "dominantWeightsBetti", -- docs
  "lexWeightsBetti", --  docs
  "numDistinctRepsBetti", -- docs
  "numRepsBetti", -- docs
  "bsCoeffs" -- docs
  }

--------------------------------------------------------------------
--------------------------------------------------------------------
----- CODE
--------------------------------------------------------------------
--------------------------------------------------------------------



--  We've had issues with multigraded Hilbert series
--  when A is not a globally defined ring.
if A === symbol A then A = QQ[t_0,t_1,t_2,MonomialOrder => Lex];


--------------------------------------------------------------------
--------------------------------------------------------------------
----- INPUT: A hash table representing a Betti table 
-----
----- OUPUT: The BettiTally for the inputed Betti table
-----
----- DESCRIPTION: This function is used to convert a hash table
----- representing a Betti Table into the more familar BettiTally.
----- The input hash table is assumed to have its keys being pairs
----- (p,q) such that H#(p,q)=>{K_{p,q}(M)}.
--------------------------------------------------------------------
--------------------------------------------------------------------
makeBettiTally = method();
makeBettiTally HashTable := H ->(
    new BettiTally from apply(keys H, h-> (h_0,{h_0+h_1},h_0+h_1)=> H#h)
    )


--------------------------------------------------------------------
--------------------------------------------------------------------
----- INPUT: (d,n,b) 
-----
----- OUPUT: (true,null) or (false, error message)
-----
----- DESCRIPTION: This unexported function is used to check whether 
----- we have  data for a particular case. If we have data for O(b)
----- on P^n embedded by O(d) then it returns (true,null) otherwise it 
----- returns (false, error message). This is primarily used as a
----- way forother functions to quickly produce error messages
--------------------------------------------------------------------
--------------------------------------------------------------------
rangeCheck = method();
rangeCheck (ZZ,ZZ,ZZ) :=(d,n,b) ->(
    message := (true,null);
    if n > 2 or n < 1 then message = (false,"Need n = 1 or 2");
    if b >= d or b < 0 then message = (false,"Need 0 <= b < d");    
    if n == 1 and d > 10 or d < 2 then message = (
	false, "If n = 1 then we need 2 <= d <= 10");
    if n == 2 and d > 8 or d < 2 then message = (
	false,"If n = 2 then we need 2 <= d <= 8");
    message
    )


--------------------------------------------------------------------
--------------------------------------------------------------------
----- INPUT: (d,n,b) 
-----
----- OUPUT: string
-----
----- DESCRIPTION: This non-exported function outputs the path for 
----- the auxilary file containing the data for O(b) on P^n embedded
----- by O(d) as a string. It is only used internally to ensure that
----- other functions look for the data files in the correct path.
----- It seems curDir needs to be defined outside of the function.
--------------------------------------------------------------------
--------------------------------------------------------------------
curDir := currentFileDirectory;
getFileName := (d,n,b) ->(curDir|"SchurVeronese/bettiP"|toString(n)|"_" | toString(d) | "_"|toString b|".m2")


--------------------------------------------------------------------
--------------------------------------------------------------------
----- INPUT: (d,n,b) 
-----
----- OUPUT: A hash table containing the multigraded Betti data for 
----- O(b) on P^n embedded by O(d).
-----
----- DESCRIPTION: This function returns a hashtable containing the 
----- multigraded Betti data for O(b) on P^n embedded by O(d). The
----- keys for this hash table are pairs (p,q) corresponding to the
----- Betti number K_{p,q}(n,b;d). Notice that the multigraded Betti
----- data is stored via a mulitgraded Hilbert series. See
----- [Sec 1.1, BEGY] for definitions.
--------------------------------------------------------------------
--------------------------------------------------------------------
multiBetti = method();
multiBetti (ZZ,ZZ,ZZ) := (d,n,b) ->(
    message := rangeCheck(d,n,b);
    if message_0 == false then error message_1;    
    load getFileName(d,n,b);
    value("mb"|toString d|toString b)
    )


--------------------------------------------------------------------
--------------------------------------------------------------------
----- INPUT: (d,n,b) 
-----
----- OUPUT: A hash table containing the Schur Betti data for 
----- O(b) on P^n embedded by O(d).
-----
----- DESCRIPTION: This function returns a hashtable containing the 
----- Schur Betti data for O(b) on P^n embedded by O(d). The
----- keys for this hash table are pairs (p,q) corresponding to the
----- Betti number K_{p,q}(n,b;d). Notice that the Schur Betti data
----- is stored as a list of pairs ({a,b,c},m) where {a,b,c} is the
----- partition describing the Schur functor and m is the multiplicity
----- with which this Shur functor appears in K_{p,q}(n,b;d). See
----- [Sec 1.2, BEGY] for definitions.
---------------------------------------------------------------------
---------------------------------------------------------------------
schurBetti = method();
schurBetti (ZZ,ZZ,ZZ) := (d,n,b) ->(
    message := rangeCheck(d,n,b);
    if message_0 == false then error message_1; 
    load getFileName(d,n,b);
    value("sb"|toString d|toString b)
    )


--------------------------------------------------------------------
--------------------------------------------------------------------
----- INPUT: (d,n,b) 
-----
----- OUPUT: A hash table containing the total graded Betti data for 
----- O(b) on P^n embedded by O(d).
-----
----- DESCRIPTION: This function returns a hash table containing the 
----- total graded Betti data for O(b) on P^n embedded by O(d). The
----- keys for this hash table are pairs (p,q) with the corresponding
----- value being the Betti number dim K_{p,q}(n,b;d). This function
----- is the same as totalBetti2 expect that the ouptut is a hash table
----- instead of a BettiTally. See [Sec 1.1, BEGY] for definitions.
---------------------------------------------------------------------
---------------------------------------------------------------------
totalBetti = method();
totalBetti (ZZ,ZZ,ZZ) := (d,n,b) ->(
    message := rangeCheck(d,n,b);
    if message_0 == false then error message_1;    
    load getFileName(d,n,b);
    value("tb"|toString d|toString b)
    )


--------------------------------------------------------------------
--------------------------------------------------------------------
----- INPUT: (d,n,b) 
-----
----- OUPUT: A BettiTally containing the total graded Betti data for 
----- O(b) on P^n embedded by O(d).
-----
----- DESCRIPTION: This function returns a BettiTally containing the 
----- total graded Betti data for O(b) on P^n embedded by O(d). This 
----- function is the same as totalBetti expect that the ouptut is a
----- BettiTally instead of a hash table. See [Sec 1.1, BEGY] for 
----- definitions.
---------------------------------------------------------------------
---------------------------------------------------------------------
totalBetti2 = method();
totalBetti2 (ZZ,ZZ,ZZ) := (d,n,b) ->(
    message := rangeCheck(d,n,b);
    if message_0 == false then error message_1;    
    load getFileName(d,n,b);
    makeBettiTally value("tb"|toString d|toString b)
    )

--------------------------------------------------------------------
--------------------------------------------------------------------
----- INPUT: (d,n,b) 
-----
----- OUPUT: A hash table containing lists of the doinant Schur 
----- functors appearing in the decomposition of the graded Betti 
----- numbers for O(b) on P^n embedded by O(d).
-----
----- DESCRIPTION: This function returns a hash table whose keys 
----- are pairs (p,q) with the corresponding value being a list
----- of the Schur functors of dominant weight appearing in the 
----- decomposition of K_{p,q}(n,b;d). Note the Schur functors 
----- are recorded via their weights {a,b,c}. See [Sec 1.3, BEGY] 
----- for definitions.
---------------------------------------------------------------------
---------------------------------------------------------------------
dominantWeightsBetti = method();
dominantWeightsBetti (ZZ,ZZ,ZZ) := (d,n,b) ->(
    message := rangeCheck(d,n,b);
    if message_0 == false then error message_1;    
    load getFileName(d,n,b);
    value("dw"|toString d|toString b)
    )


--------------------------------------------------------------------
--------------------------------------------------------------------
----- INPUT: (d,n,b) 
-----
----- OUPUT: A hash table containing lists of the lex-leading
----- weights of the Schur functors appearing in the decomposition 
----- of the graded Betti numbers for O(b) on P^n embedded by O(d).
-----
----- DESCRIPTION: This function returns a hash table whose keys 
----- are pairs (p,q) with the corresponding value being a list
----- of the Schur functors of lex-leading weight appearing in the 
----- decomposition of K_{p,q}(n,b;d). Note the Schur functors 
----- are recorded via their weights {a,b,c}. See [Sec 1.3, BEGY] 
----- for definitions.
---------------------------------------------------------------------
---------------------------------------------------------------------
lexWeightsBetti = method();
lexWeightsBetti (ZZ,ZZ,ZZ) := (d,n,b) ->(
    message := rangeCheck(d,n,b);
    if message_0 == false then error message_1;    
    A := QQ[t_0,t_1,t_2, MonomialOrder => Lex];
    load getFileName(d,n,b);
    value("lw"|toString d|toString b)
    )


--------------------------------------------------------------------
--------------------------------------------------------------------
----- INPUT: (d,n,b) 
-----
----- OUPUT: A hash table containing the number of distinct Schur 
----- functors appearing in the decomposition of the graded Betti 
----- numbers for O(b) on P^n embedded by O(d).
-----
----- DESCRIPTION: This function returns a hash table whose keys 
----- are pairs (p,q) with the corresponding value being the number
----- of distinct Schur functors appearing in the decomposition of
----- K_{p,q}(n,b;d). See [Sec 1.2, BEGY] for definitions.
---------------------------------------------------------------------
---------------------------------------------------------------------
numDistinctRepsBetti = method();
numDistinctRepsBetti  (ZZ,ZZ,ZZ) := (d,n,b) ->(
    message := rangeCheck(d,n,b);
    if message_0 == false then error message_1;    
    load getFileName(d,n,b);
    value("nr"|toString d|toString b)
    )


--------------------------------------------------------------------
--------------------------------------------------------------------
----- INPUT: (d,n,b) 
-----
----- OUPUT: A hash table containing the number of Schur functors
----- appearing in the decomposition of the graded Betti numbers 
----- for O(b) on P^n embedded by O(d).
-----
----- DESCRIPTION: This function returns a hash table whose keys 
----- are pairs (p,q) with the corresponding value being the number
----- of Schur functors counted with multiplicity appearing in the 
----- decomposition of K_{p,q}(n,b;d). See [Sec 1.2, BEGY] for
----- definitions.
---------------------------------------------------------------------
---------------------------------------------------------------------
numRepsBetti = method();
numRepsBetti  (ZZ,ZZ,ZZ) := (d,n,b) ->(
    message := rangeCheck(d,n,b);
    if message_0 == false then error message_1;    
    load getFileName(d,n,b);
    value("nrm"|toString d|toString b)
    )

--------------------------------------------------------------------
--------------------------------------------------------------------
----- INPUT: (d,n,b) 
-----
----- OUPUT: A hash table containing the errors found when computing
----- the multigrded Betti numbers for O(b) on P^n embedded by O(d).
-----
----- DESCRIPTION: As the methods we use to compute multigraded 
----- Betti numbers are numerical in nature there is room for error,
----- and so, we have implemented post processing to catch errors. 
----- This function returns a hash table whose keys are pairs (p,q)
----- with the corresponding value being a multigraded Hilber 
----- series recording the errors encounted when computing the
----- multigraded Betti numbers for K_{p,q}(n,b;d). See [Sec 5.2, BEGY]
----- for a discussion on error processing.
---------------------------------------------------------------------
---------------------------------------------------------------------
errorBetti = method();
errorBetti  (ZZ,ZZ,ZZ) := (d,n,b) ->(
    message := rangeCheck(d,n,b);
    if message_0 == false then error message_1;    
    load getFileName(d,n,b);
    value("er"|toString d|toString b)
    )

--------------------------------------------------------------------
--------------------------------------------------------------------
----- INPUT: (d,n,b) 
-----
----- OUPUT: A list of the Boij-Soederberg coefficents for the 
----- decomposition of the Betti table of O(b) on P^n embedded by O(d).
-----
----- DESCRIPTION: This funnction returns a ist of the Boij-Soederberg 
----- coefficents for the  decomposition of the Betti table of O(b) on
----- P^n embedded by O(d). See [Sec 6.3, BEGY] for definitions. 
---------------------------------------------------------------------
---------------------------------------------------------------------
bsCoeffs = method();
bsCoeffs  (ZZ,ZZ,ZZ) := (d,n,b) ->(
    message := rangeCheck(d,n,b);
    if message_0 == false then error message_1;    
    load getFileName(d,n,b);
    value("bs"|toString d|toString b)
    )


--------------------------------------------------------------------
--------------------------------------------------------------------
----- Begining of the tests and the documentation
--------------------------------------------------------------------
--------------------------------------------------------------------

load ("./tests.m2")
beginDocumentation()
load ("./doc.m2")


end--;



--------------------------------------------------------------------
--------------------------------------------------------------------
----- Begining of sandbox
--------------------------------------------------------------------
--------------------------------------------------------------------


---
---
---
restart
uninstallPackage "SchurVeronese"
restart
installPackage "SchurVeronese"
check "SchurVeronese"
viewHelp SchurVeronese
viewHelp totalBetti2










needsPackage "SchurVeronese"
redundantEntryTest = (B,p,q) -> abs(1 - (B#(p,q) / B#(p-1,q+1)))_RR
redundantEntryTest(tb,8,0)
redundantTest = B ->(
    P := max apply(keys B, l-> l_0);
    rL := {};
    scan(toList(1..P), p->(
	    if (B#(p,0))*(B#(p-1,1)) != 0 then rL = rL|{redundantEntryTest(B,p,0)}
	));
    scan(toList(1..P), p->(
	    if (B#(p,1))*(B#(p-1,2)) != 0 then rL = rL|{redundantEntryTest(B,p,1)}
	));
    min rL
    )
redundantTest(B)



--d=4
L = apply(toList(0..3),b-> totalBetti(4,2,b));
apply(L,l-> redundantTest(l))
1 - (252/450)
--d=5
L = apply(toList(0..4),b-> totalBetti(5,2,b));
apply(L,l-> redundantTest(l))
--d=6
L = apply(toList(0..3),b-> totalBetti(6,2,b));
apply(L,l-> redundantTest(l))



--BS coefficients
L = bsCoeffs(6,2,0);
m = max apply(L,l-> floor log(10,l));
bL = apply(L,l-> (l/10^m)_RR)




--Error analysis
--d=5
L = apply(toList(0..4),b-> makeBettiTally errorBetti(5,2,b));
netList L
--d=6
L = apply(toList(0..3),b-> makeBettiTally errorBetti(6,2,b));
netList L


-- we had serious issues with (6,2,3) and (7,0), (8,0), and (9,0) entries
break
restart
load "step3.m2"
startB = startBetti(6,2,3);
(d,n,b) = (6,2,3)
(p,q) = (8,0)
adjustBetti = (d,n,b,startB) -> (
    mb := startB_0;
    sb := startB_1;
    er := new MutableHashTable;    
    L := relevantRange(d,n,b);
    cL := apply(L, l->((l#0-1,l#1+1)));
    apply(#L,l->(
--	     l = 4
	    newPart := - mb#(L#l) +  value get concatenate{"mgHilbertSeries/P",toString(n),"/", toString(n), "_", toString(d), "_", toString(b), "-", toString((L#l)#0), "_", toString((L#l)#1), ".txt"};
sb#(8,0)
	    oldEntry := mb#(L#l);
	    newPartFixed := fixMultiHilb newPart;
	    --  we adjust multigraded entries
	    mb#(L#l) = oldEntry + newPartFixed_1;
	    mb#(cL#l) = mb#(cL#l) + newPartFixed_1;
            --  we adjust schur entries
    	    sb#(L#l) = sb#(L#l) | newPartFixed_0;
	    sb#(cL#l) = sb#(cL#l) | newPartFixed_0;
	    --  we encode errors
	    er#(L#l) = newPartFixed_2;
    	    ));
    (mb,sb,er)     
    )

needsPackage "SchurVeronese"
TB = totalBetti(6,2,3);
TB#(8,0) == rankP((multiBetti(6,2,3))#(8,0))
altP = (multiBetti(6,2,3))#(8,0);
netList sort newFixed_0
newFixed_2

anotherTry = newFixed_2 + (1/2)*5*symmP({23,17,17})
anotherFixed = fixMultiHilb(anotherTry)


symmP = L -> sum(permutations(length L),K-> t_(K_0)^(L_0)*t_(K_1)^(L_1)*t_(K_2)^(L_2))
guessAtNewPart = newPart - 42*(symmP({24,18,15}));
newFixed = fixMultiHilb guessAtNewPart;
rankP(newFixed_1)
rankP(newPart)


rankP = P -> sub(P,{t_0=>1,t_1=>1,t_2=>1})
rankP(newPart)
netList newPartFixed_0
rankP(newPartFixed_1)


decomposeHilb = (P) -> (
    L := {};
    if P == 0 then return {{},0};
    while leadCoefficient P > 0 do(
      	  K := decomposeHilb1(P);
	  print(K_1);
	  L = L|{K#1};
	  P = K#0;
	  );
    (L,P)
    )


H = hilbSeries({25,21,11})

decomposeHilb(newPart)
newPart = newPart - 38*t_0^24*t_1^18*t_2^15



rankP(oldEntry)
P = value get concatenate{"mgHilbertSeries/P",toString(n),"/", toString(n), "_", toString(d), "_", toString(b), "-", toString((L#l)#0), "_", toString((L#l)#1), ".txt"};
rankP(P)
P
(mb,sb,er) := adjustBetti(d,n,b,startBetti(d,n,b));

E1 = 
(errorBetti(6,2,3))#(8,0) == 0

mb = (multiBetti(6,2,3))#(8,0);
sb = (schurBetti(6,2,3))#(8,0);
sb

B = tb
apply(keys B, l-> if (B#l)*(B#(l_0-1,l_1+1)) != 0 then redundantEntryTest(B,l_0,l_1))
max

dw = dominantWeightsBetti(4,2,3);
nr = numRepsBetti(4,2,3);
ndr = numDistinctRepsBetti(4,2,3);



--  This is all just junk
viewHelp SchurVeronese

needsPackage "TensorComplexes"
viewHelp tensorComplexes
totalBetti2(4,2,2)
totalBetti2(8,2,0)
(schurBetti(3,2,0))#(1,1)
(schurBetti(4,2,0))#(1,1)
(schurBetti(5,2,0))#(1,1)
(schurBetti(6,2,0))#(1,1)
(schurBetti(7,2,0))#(1,1)
(schurBetti(8,2,0))#(24,1)

tb = totalBetti(5,2,0);
mean = 8.6
-- want to write b = mean + (a/2)* sqrt(mean)





--  I entered tb60 into Macaulay2.
H = tb60;
--  I wanted a reasonable range of input values, avoiding the tail.
--  So I chose {6,7,...,16}, since that avoid the relevant range.
--  For the mean, I chose exactly 12 because it looked good.
--  The "mean" choice is hacky.
R = toList(6..16);
mean = 12;
dist = (b,mean) ->(
    (b - mean)*sqrt(2)/sqrt(mean)
    )

transformToNormal = (H,R,q)->(
    mean := (1/2)*(min R + max R);
    total = sum apply(R, i->tb#(i,q));
    apply(R, i->(dist(i,mean),(tb#(i,q)/total)_RR))
    )


nd = a ->(
    (1/(2*pi))*exp(-(a^2)/2)
    )
chiSquared = L ->(
    sum apply(L, l-> (l#1 - nd(l#0))^2)
    )
 L = transformToNormal(H,R,1)
altL = apply(L,l-> (l#0,nd(l#0)))
chiSquared(L)

netList L
netList altL

-- nonzero from 1 to 21 but "infinite" at 16.
R = toList(6..16)

L = transformToNormal(H,R,1)

tikz = toExternalString(L#0)
scan(#L,i-> tikz = tikz|"--"
apply(keys H,i-> if i#1 == 1 and H#i != 0 then print i#0)
dist(8,mean)
total = sum apply(toList(4..12), i->tb#(i,1))
L = apply(toList(4..12), i->(dist(i,mean),(tb#(i,1)/total)_RR))


--  figuring out duality
needsPackage "SchurVeronese"
viewHelp SchurVeronese
mb = multiBetti(5,2,2);
mb#(18,2)




(d,n,b) = (5,2,2)
--  write in terms of d/b
whatsTheDualGuy = (d,n,b) ->(
    r := ((binomial(d+n,d)-1)*d+b)//(n+1);    
    product apply(n+1, i-> t_i^r)
    )

produceTheDualHilbSeries = (d,n,b)->(
    P := --some crazy code drawing from a file
    dG := whatsTheDualGuy(d,n,b);
    (M,C) := coefficients P;
    ML := matrix {apply(flatten entries M, m-> sub(dualGuy/m, ring m))};
    Q := ML*C;
    --  write Q to a file.
    )
