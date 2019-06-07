LoadPackage("Gauss");

d := SplitString(StringFile("out-23-test-new/matrices/map_8_2/multidegree_4_16_12_18.dat"),"\n");

SMC := NewDictionary(1,true);
SMV := NewDictionary(1,true);

maxC := 1;
maxR := 1;

for i in [1..Length(d)] do
    e := SplitString(d[i], " ");
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

cols := [];
vals := [];

for r in [1..maxR] do
    Add(cols, LookupDictionary(SMC,r));
    Add(vals, LookupDictionary(SMV,r));
od;

SM := SparseMatrix(maxR, maxC, cols, vals*One(GF(32003))) ;
Rank(SM);
maxR;