uninstallPackage "HirzebruchSyzygies"
restart
installPackage "HirzebruchSyzygies"
needsPackage "BoijSoederberg"
needsPackage "Posets";



dataRange = value get "dataRange.m2"

tally apply(dataRange,D->(if (D#1)#0 != 1 then D#0))

L = delete(,apply(dataRange,D->(if (D#1)#0 != 1 and D#0 == {1,1} then D)))

g = openOut("tablesForPaper11.txt")
apply(L,D->(
	g << "\\[";
	g<< endl;
	g << "\\beta((1,1);("|toString((D#1)#0)|","|toString((D#1)#1)|"))=";
	g<< endl;
	g<< tex totalBettiTally(0,D#0,D#1);
	g<< endl;
	g << "\\]";
	g<< endl;
	g<< endl;
	))
close g
