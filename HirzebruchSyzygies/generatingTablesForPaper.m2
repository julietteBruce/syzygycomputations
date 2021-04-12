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

H = schurBetti(0,{1,2},{2,3})
g = openOut("schurExampleForPaper.txt")
apply(sort keys H,k->(
	if H#k != {} then (
	    g << "K_{"|toString(k#0)|","|toString(k#1)|"}(\\PP^{1}\\times\\PP^{1},(1,2);(2,3))\\;=\\;&";
	    apply(sort H#k,i->(
		    if i#1 == 1 then (
			j = i#0;
			g<< "\\bS_{("|toString(j#0)|","|toString(j#1)|","|toString(j#2)|","|toString(j#3)|")}\\oplus";
			);
		    if i#1 != 1 then (
			j = i#0;
			g<< "\\bS_{("|toString(j#0)|","|toString(j#1)|","|toString(j#2)|","|toString(j#3)|")}^{\\oplus"|toString(i#1)|"}\\oplus";
			);
		    ));
	    g<< endl;
	    g<< endl;
	    );
	))
close g
