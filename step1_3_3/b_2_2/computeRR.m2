load "src/relevantRange.m2"
d={3,3}
b={2,2}
rR=relevantRange(0,b,d)
g=openOut "step1_3_3/b_2_2/relevantpq.txt"
g<< "1 1" << endl;
g<< concatenate(toString(d#0), " ", toString(d#1)) << endl;
g<< concatenate(toString(b#0), " ", toString(b#1)) << endl;
for i from 0 to #rR-1 do (
    g<<concatenate(toString(rR#i#0), " ", toString(rR#i#1)) << endl;
    );
g<<close;