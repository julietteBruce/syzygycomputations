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
-- #6 - fix makeBettiTally to workquit with funkier things.  e.g makeBettiTally 
--      Daniel:  I worked on this for a while and I can't get it to work.
-- #7 - What's wrong with the \mathcal{O} in the doc?
--      We can go ahead without fixing this one.
--------------------------------------------------------------------------------



newPackage("SchurVeronese",
    Version => "1.1",
    Date => "13 May 2019",
    Headline => "Data for Veronese embeddings of projective space",
    Authors => {
        {Name => "Juliette Bruce",           Email => "jebruce2@wisc.edu",       HomePage => "https://juliettebruce.github.io"},
        {Name => "Daniel Erman",             Email => "derman@math.wisc.edu",    HomePage => "http://www.math.wisc.edu/~derman/"},	     
        {Name => "Steve Goldstein",          Email => "sgoldstein@wisc.edu",   HomePage => ""},
	{Name => "Jay Yang",                 Email => "jkyang@umn.edu",   HomePage => "http://www-users.math.umn.edu/~jkyang/"}
	},
  DebuggingMode => true,
  AuxiliaryFiles => true
  )

export {
  "makeBettiTally", -- docs
  "multiBetti", -- docs
  "schurBetti", -- docs
  "totalBetti", -- docs
  "totalBettiTally", -- docs
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
--if A === symbol A then A = QQ[t_0,t_1,t_2,MonomialOrder => Lex];


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
    --if n > 2 or n < 1 then message = (false,"Need n = 1 or 2");
    --if b >= d or b < 0 then message = (false,"Need 0 <= b < d");    
    --if n == 1 and d > 10 or d < 2 then message = (
	--false, "If n = 1 then we need 2 <= d <= 10");
    --if n == 2 and d > 8 or d < 2 then message = (
	--false,"If n = 2 then we need 2 <= d <= 8");
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
getFileName := (a,B,D) ->(curDir|"M2OutputFiles/bettiF"|toString(a)|"_" | toString(B#0) | "_"| toString(B#1) | "_"| toString(D#0) | "_"| toString(D#1) |".m2")

bettiF0_0_0_2_4
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
    if message_0 == false then return message_1;
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
    if message_0 == false then return message_1; 
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
----- is the same as totalBettiTally expect that the ouptut is a hash table
----- instead of a BettiTally. See [Sec 1.1, BEGY] for definitions.
---------------------------------------------------------------------
---------------------------------------------------------------------
totalBetti = method();
totalBetti (ZZ,ZZ,ZZ) := (d,n,b) ->(
    message := rangeCheck(d,n,b);
    if message_0 == false then return message_1;    
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
totalBettiTally = method();
totalBettiTally (ZZ,ZZ,ZZ) := (d,n,b) ->(
    message := rangeCheck(d,n,b);
    if message_0 == false then return message_1;
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
    if message_0 == false then return message_1;    
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
    if message_0 == false then return message_1;    
    --A := QQ[t_0,t_1,t_2, MonomialOrder => Lex];
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
    if message_0 == false then return message_1;    
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
    if message_0 == false then return message_1;    
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
    if message_0 == false then return message_1;    
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
    if message_0 == false then return message_1;    
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
load ("./doc2.m2")
--Daniel: I screwed something up so I temporarily renamed the doc file.

end--;



--------------------------------------------------------------------
--------------------------------------------------------------------
----- Begining of sandbox
--------------------------------------------------------------------
--------------------------------------------------------------------

---
---
restart
uninstallPackage "SchurVeronese"
restart
installPackage "SchurVeronese"
check "SchurVeronese"
viewHelp SchurVeronese
