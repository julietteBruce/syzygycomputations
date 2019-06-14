## Compute the rank of a sparse matrix over the rationals.

## Run like this to avoid hard-coding the matrix file.
##  (export matrixFile=multidegree_4_7_28.dat;
##             cat ranks.rational.gap | gap -q 1> rank.out 2> rank.err ) &

LoadPackage("Gauss");;
matrixFile := GAPInfo.SystemEnvironment.matrixFile;;
## Or hard-code value
#matrixFile := "multidegree_4_7_28.dat";;

d := SplitString(StringFile(matrixFile),"\n");;

SMC := NewDictionary(1,true);;
SMV := NewDictionary(1,true);;

maxC := 1;;
maxR := 1;;

for i in [1..Length(d)] do
    e := SplitString(d[i], " ", "\t");
    
    r := Int(e[1]);
    c := Int(e[2]);
    ## Coerce the sparse matrix into Q;
    v := Rat(Rat(e[3])/2);
    
    if KnowsDictionary(SMC, r) then
        Add(LookupDictionary(SMC, r), c);
        Add(LookupDictionary(SMV, r), v);
    else
	AddDictionary(SMC, r, [c]);
	AddDictionary(SMV, r, [v]);
    fi;
    
    if r>maxR then
        maxR := r;
    fi;
    if c>maxC then
        maxC := c;
    fi;
od;

cols := [];;
vals := [];;

for r in [1..maxR] do
    Add(cols, LookupDictionary(SMC,r));
    Add(vals, LookupDictionary(SMV,r));
    SortParallel(cols[Length(cols)], vals[Length(cols)]);
od;

SM := SparseMatrix(maxR, maxC, cols, vals);;

Rank(SM);
maxR;
quit;
