## Compute the rank of a sparse matrix over a finite field.

## Run like this to avoid hard-coding the matrix file and the prime.
##  (export p=32003;  export matrixFile=multidegree_4_7_28.dat;
##                 cat ranks.gap | gap -q 1> rank.out 2> rank.err ) &

LoadPackage("Gauss");;

matrixFile := GAPInfo.SystemEnvironment.matrixFile;;
p := GAPInfo.SystemEnvironment.p;;
p := Int(p);;

## Or hard-code values 
#matixFile := "multidegree_4_7_28.dat";;
#p := 32003;

d := SplitString(StringFile(matrixFile),"\n");;

SMC := NewDictionary(1,true);;
SMV := NewDictionary(1,true);;

maxC := 1;;
maxR := 1;;

for i in [1..Length(d)] do
    e := SplitString(d[i], " ", "\t");
    
    r := Int(e[1]);
    c := Int(e[2]);
    v := Int(e[3]);
    
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

SM := SparseMatrix(maxR, maxC, cols, vals*One(GF(p)));;

Rank(SM);
maxR;
quit;
