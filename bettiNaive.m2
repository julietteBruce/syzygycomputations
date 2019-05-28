S=QQ[x0,x1]**QQ[y0,y1]
segreRing = (d1,d2) -> QQ[apply((0..d1), i-> apply((0..d2), j-> value("z"|d1-i|i|d2-j|j)  )) ]
segreVeronese = (d1, d2) -> map(S,segreRing(d1,d2), 
    flatten toList apply((0..d1), i-> toList apply((0..d2), j-> x0^(d1-i)*x1^i*y0^(d2-j)*y1^j)  ))

bettiP1P1 = (d1,d2) -> (
    phi := segreVeronese(d1,d2);
    C := resolution kernel phi;
    betti C 
    )

