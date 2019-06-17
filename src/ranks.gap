## Compute the rank of a sparse matrix over a finite field.

## Run like this to avoid hard-coding the matrix file and the prime.
##  (export p=32003;  export matrixFile=multidegree_4_7_28.dat;
##                 cat ranks.gap | gap -q 1> rank.out 2> rank.err ) &

LoadPackage("Gauss");;

matrixFile := GAPInfo.SystemEnvironment.matrixFile;;

if (IsBound(GAPInfo.SystemEnvironment.p)) then
    p := GAPInfo.SystemEnvironment.p;;
    p := Int(p);;
else
    p := 0;
fi;

## Or hard-code values, e.g.
#matrixFile := "multidegree_4_7_28.dat";;
#p := 32003;

d := SplitString(StringFile(matrixFile),"\n");;

SMR := NewDictionary(1,true);;
SMC := NewDictionary(1,true);;
SMV := NewDictionary(1,true);;

maxC := 1;;
maxR := 0;;

for i in [1..Length(d)] do
    e := SplitString(d[i], " ", "\t");
    
    r := Int(e[1]);
    c := Int(e[2]);
    v := Int(e[3]);
    if p = 0 then
        v := Rat(e[3])/2;
    else
        v := Int(e[3]);
    fi;
    
    if KnowsDictionary(SMR,r) then
        rowIndex := LookupDictionary(SMR,r);
    else
        maxR := maxR +1;
        rowIndex := maxR;
        AddDictionary(SMR,r,rowIndex);
    fi;
        
    if KnowsDictionary(SMC, rowIndex) then
        Add(LookupDictionary(SMC, rowIndex), c);
        Add(LookupDictionary(SMV, rowIndex), v);
    else
	AddDictionary(SMC, rowIndex, [c]);
	AddDictionary(SMV, rowIndex, [v]);
    fi;
    
    if c>maxC then
        maxC := c;
    fi;
od;

cols := [];;
vals := [];;

for rowIndex in [1..maxR] do
    Add(cols, LookupDictionary(SMC,rowIndex));
    Add(vals, LookupDictionary(SMV,rowIndex));
    SortParallel(cols[Length(cols)], vals[Length(cols)]);
od;

if p = 0 then
    SM := SparseMatrix(Length(vals), maxC, cols, vals);
else
    SM := SparseMatrix(Length(cols), maxC, cols, vals*One(GF(p)));
fi;

Rank(SM);
maxR;

if (IsBound(GAPInfo.SystemEnvironment.RankStats)) then 
    time := StringTime(Runtime());
    memBytes := TotalMemoryAllocated();;
    memMB := Int(memBytes/(1024*1024));
    Print(time, " seconds\t", memMB, "MB RAM\n");
fi;

quit;
