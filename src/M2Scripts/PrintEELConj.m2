load "Rep_Theory.m2"
load "EEL.m2"


weightsToCheck = (d,p,q,b)->(
    L := betterDominantWeights EELWeights(d,p,q,b,2);
    K := unique apply(flatten apply(L,l-> {l+{1,-1,0},l+{1,0,-1},l+{0,1,-1}}), k-> reverse sort k);
    K
    )


printWeightsToCheck = (d,p,q,b,outfile) -> (
    outfile = outfile << "";--some magic that opens the file
    K := apply(weightsToCheck(d,p,q,b),reverse);
    for md in K do (
        for i in md do (outfile << i << " ");
        outfile << endl);
    outfile << close
    )

--Note, this is a bit evil, since value evaluates arbitrary strines
printWeightsToCheck(value(scriptCommandLine#1),value(scriptCommandLine#2),value(scriptCommandLine#3),value(scriptCommandLine#4),scriptCommandLine#5);
