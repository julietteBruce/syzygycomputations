--  This provides the code for building all of the Betti tables and 
--  storing them as output files.


--
--  Building for P2 for d = 2,3,4 and all b values
--  Time: 45 seconds
--
restart
load "step3.m2"
L = flatten apply(toList(2..4),d-> apply(toList(0..d-1),b-> (d,2,b)))
time scan(L,l-> makeOutputFiles(l))


--
--  Building for P2 for d = 5 and b = 0,1,2
--  Time: ??
--  Still awaiting missing data for b = 3/4 and d = 5
--
restart
load "step3.m2"
L = apply({0,1,2}, b-> (5,2,b));
time scan(L,l-> makeOutputFiles(l))

--
--  Building d = 5 and b = 3/4 but w/o corrections.
--  Time: 38 seconds
--
restart
load "step3.m2"
time scan({3,4}, b-> makePartialOutput(5,2,b));



--  Building partial outputs for d = 6.
--  "no Multi" takes about 3 seconds for each one
--  with multi takes about 270 seconds for each one (not fully run yet)
restart
load "step3.m2"
time scan({0,1,2,3,4,5}, b-> makePartialOutput(6,2,b));
--time scan({0,1,2,3,4,5}, b-> makePartialNoMultiOutput(6,2,b));



--
--  Building outputs w/ no multigraded Hilbert series for d = 7
--  Total time: 67 seconds
--
restart
load "step3.m2"
time scan({0,1,2,3,4,5,6},b-> makePartialNoMultiOutput(7,2,b))


--
--  Building outputs w/ no multigraded Hilbert series for d = 8
--  Total time: 562 seconds
restart
load "step3.m2"
time scan({0,1,2,3,4,5,6,7},b-> makePartialNoMultiOutput(8,2,b))





--  
--  Building stuff for P1
--  We'd already built them once so we just filled in missing data.
--  Time: 1.5 seconds
--  
restart
needsPackage "SchurVeronese"
load "step3.m2"
L = flatten apply(toList(2..10),d-> apply(toList(0..d-1),b-> (d,1,b)))
time scan(L,l-> makeP1Output(l))


























