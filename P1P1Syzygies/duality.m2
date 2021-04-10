needsPackage "P1P1Syzygies"

dualTotalBetti = (a,B,D) -> (
    TB := totalBettiP1P1(a,B,D);
    return  hashTable apply(keys TB, PQ -> {((D#0+1)*(D#1+1)-PQ#0-3, 2-PQ#1), TB#PQ} )
    )


dualSchurBetti = (a,B,D) -> (
    TB := totalBettiP1P1(a,B,D);
    dTB := dualTotalBetti(a,B,D);
    SB := schurBettiP1P1(a,B,D);
    D0D1 := (D#0+1)*(D#1+1);
    alpha := {(D0D1*D#0-2)//2, (D0D1*D#0-2)//2, (D0D1*D#1-2)//2, (D0D1*D#1-2)//2 };
    dualPQ := PQ -> ((D#0+1)*(D#1+1)-PQ#0-3, 2-PQ#1);
    wOpp := w -> ( v := alpha - w; return {v#1,v#0,v#3,v#2} );
    --Lt = (SB#(12,2))#0;
    --Lt#
    dualSw := PQ -> apply(SB#PQ, L -> (wOpp(L#0), L#1) );
    dualSB := hashTable apply(keys SB, PQ -> {dualPQ(PQ), dualSw(PQ)} );   
    return dualSB
    )


-- dualDominantWeight = (a,B,D) -> (
--     TB := totalBetti(a,B,D);
--     dTB := dualTotalBetti(a,B,D);
--     SB := schurBetti(a,B,D);
--     D0D1 := (D#0+1)*(D#1+1);
--     alpha := {(D0D1*D#0-2)//2, (D0D1*D#0-2)//2, (D0D1*D#1-2)//2, (D0D1*D#1-2)//2 };
--     dualPQ := PQ -> ((D#0+1)*(D#1+1)-PQ#0-3, 2-PQ#1);
--     wOpp := w -> ( v := alpha - w; return {v#1,v#0,v#3,v#2} );
--     --Lt = (SB#(12,2))#0;
--     --Lt#
--     dualSw := PQ -> apply(SB#PQ, L -> (wOpp(L#0), L#1) );
--     dualSB := hashTable apply(keys SB, PQ -> {dualPQ(PQ), dualSw(PQ)} );   
--     return dualSB
--     )


makeDualBettiFile = (a,B,D) -> (
    dB := D-{2,2}-B;
    dD := D;
    fName := "P1P1Syzygies/bettiF0_"|dB#0|"_"|dB#1|"_"|dD#0|"_"|dD#1|".m2";
    f := openOut(fName);
    f << "A := QQ[t_0,t_1,t_2,t_3];" << endl;
    f << "--tb stands for Total Betti numbers" << endl;
    f << "tb"|dB#0|dB#1|dD#0|dD#1|" = new HashTable from "|toString(dualTotalBetti(a,B,D)) << endl
    f << "--sb represents the betti numbers as sums of Schur functors" << endl
    f << "sb"|dB#0|dB#1|dD#0|dD#1|" = new HashTable from "|toString(dualSchurBetti(a,B,D)) << endl 
--    f << "--dw stands for dominant weights" << endl
--    f << "dw"|dB#0|dB#1|dD#0|dD#1|" = new HashTable from "|toString(dualDominantWeight(a,B,D)) << endl 
    f << end;
    close f;
    return fName;
    )



--to add
--0 0 2 10
--0 0 2 9
--0 0 3 6
--0 1 3 6
--0 1 4 5
--0 2 3 6
--0 2 3 7
--0 3 3 6
--0 3 3 7
--0 3 4 6
--0 4 3 6
--1 0 3 8
--1 0 3 9
--1 0 4 5
--1 1 3 8
--1 1 4 5
--1 1 4 6
--2 0 4 6
--2 0 4 7

















