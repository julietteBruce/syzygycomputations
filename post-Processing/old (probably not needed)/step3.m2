restart
needsPackage "SchurRings"
load "Rep_Theory.m2"

--  Insanely long code for building all Betti tables
--  Input is e = degree of embedding; 
--           n = dimension of projective space (needs to be 2 for now)
--           b = degree of line bundle you're resolving
--  Output a mutligraded Betti table assuming no cancelation amoungst Schur functors.
startBetti = (e,n,b)->(
    k = n+1;
    d = e;
    r = b;
    if k != 3 then error "This code only works for P^2";
    load "veroP2.m2";
    schurBetti = new MutableHashTable;
    keyList = flatten apply(#koz, p-> apply(#(koz#0),j-> (p,j)));
    scan(keyList, i ->(  schurBetti#i = standardize(koz#(i_0))#(i_1)));
    --  Some weird step to avoid conflicts of variable names.
    t = symbol t;
    A = QQ[t_0,t_1,t_2];
    --  Building the Betti tables
    multiBetti = new MutableHashTable;
    scan(keyList, i->( multiBetti#i = totalHilb(schurBetti#i,A)));
    multiBetti
    )

mb = startBetti(5,2,0);
MB = mb;
L = {(14,1),(15,1)}
cancelBetti = (e,n,b,L)->(
    apply(L,l->(
	mb#l= value get concatenate{"../mg-generating-functions/mgHilbertSeries-", toString(n), "_", toString(e), "_", toString(b), "-", toString(l#0), "_", toString(l#1), ".txt"};
    ));
)

f = value get  "../mg-generating-functions/poly-2_5_0-multigrad_15_1.txt"
A = QQ[t_0,t_1,t_2]


    totalBetti = new MutableHashTable;
    scan(keyList, i->  totalBetti#i = totalRank(multiBetti#i));
    domWeights  = new MutableHashTable;
    scan(keyList, i->  domWeights#i = maximalDominantWeights(apply(schurBetti#i, l -> l_0)));
    numReps  = new MutableHashTable;
    scan(keyList, i->  numReps#i = #(schurBetti#i));
    numRepsMulti  = new MutableHashTable;
    scan(keyList, i->  numRepsMulti#i = sum apply(schurBetti#i, j-> j_1));
    lexLeadingWeights = new MutableHashTable;
    scan(keyList, i->  lexLeadingWeights#i = if schurBetti#i == {} then {} else last sort apply(schurBetti#i,j-> j_0));
    --  Writing the files
    b = r;
    g = openOut ("bettiP2_"|toString d|"_"|toString b|".m2");
    g<< "--This file computes Betti tables for P^2 for d = "|toString d|" and b = "|toString b;
    -- need to add bit to initialize A = QQ[t_0,t_1,t_2];
    g<< endl;
    g<< "--tb stands for Total Betti numbers";
    g<< endl;
    g<< "tb"|toString d|toString b|" = ";
    g<< toExternalString unMute(totalBetti);
    g<< ";";
    g<< endl; 
    g<< "--mb stands for Multigraded Betti numbers";
    g<< endl ;
    g<< "mb"|toString d|toString b|" = ";
    g<< toExternalString unMute(multiBetti);
    g<< ";";
    g<< endl;
    g<< "--sb represents the betti numbers as sums of Schur functors";
    g<< endl ;
    g<< "sb"|toString d|toString b|" = ";
    g<< toExternalString unMute(schurBetti);
    g<< ";";
    g<< endl;
    g<< "--dw encodes the dominant weights in each entry";
    g<< endl ;
    g<< "dw"|toString d|toString b|" = ";
    g<< toExternalString unMute(domWeights);
    g<< ";";
    g<< endl;
    g<< "--lw encodes the lex leading weight in each entry";
    g<< endl ;
    g<< "lw"|toString d|toString b|" = ";
    g<< toExternalString unMute(lexLeadingWeights);
    g<< ";";
    g<< endl;
    g<< "--nr encodes the number of disctinct reprsentations in each entry";
    g<< endl ;
    g<< "nr"|toString d|toString b|" = ";
    g<< toExternalString unMute(numReps);
    g<< ";";
    g<< endl;
    g<< "--nrm encodes the number of representations with multiplicity in each entry";
    g<< endl ;
    g<< "nrm"|toString d|toString b|" = ";
    g<< toExternalString unMute(numRepsMulti);
    g<< ";";
    g<< endl;
    g<< "end;";
    close g;    
    )