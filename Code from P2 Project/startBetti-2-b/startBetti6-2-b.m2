load "step3.m2"
(d,n,b) = (6,2,0)
startB = timing startBetti(6,2,0);
ans = {startB#0, (unMute (startB#1)#0, unMute (startB#1)#1)};
g = openOut "schurBetti6-2-0.txt";
g << toExternalString ans
close g
quit
