needsPackage "SchurVeronese"
viewHelp SchurVeronese
needsPackage "NormalToricVarieties"
viewHelp NormalToricVarieties
code smoothFanoToricVariety


g = openOut "schurExamplesPapert.tex"
g << "\\begin{align*}"
g << endl;
g<< "K_{14,1}(5;0)\\cong";
L1 = rsort (schurBetti(5,2,0))#(14,1)
apply(L1,i->(
	g<< "\\bS_{(";
	g<< concatenate{toString((i#0)#0), ",", toString((i#0)#1), ",", toString((i#0)#2)};
	g<< ")}^{";
	g<< toString(i#1);
	g<< "}\\oplus";
	))
g << endl;
g << "\\end{align*}"
g << endl;
g << endl;

g << "\\begin{align*}"
g << endl;
g<< "K_{15,1}(5;0)\\cong";
L1 = rsort (schurBetti(5,2,0))#(15,1)
apply(L1,i->(
	g<< "\\bS_{(";
	g<< concatenate{toString((i#0)#0), ",", toString((i#0)#1), ",", toString((i#0)#2)};
	g<< ")}^{";
	g<< toString(i#1);
	g<< "}\\oplus";
	))
g << endl;
g << "\\end{align*}"
g << endl;
g << endl;

g << "\\begin{align*}"
g << endl;
g<< "K_{13,2}(5;0)\\cong";
L1 = rsort (schurBetti(5,2,0))#(13,2)
apply(L1,i->(
	g<< "\\bS_{(";
	g<< concatenate{toString((i#0)#0), ",", toString((i#0)#1), ",", toString((i#0)#2)};
	g<< ")}^{";
	g<< toString(i#1);
	g<< "}\\oplus";
	))
g << endl;
g << "\\end{align*}"
g << endl;
g << endl;

g << "\\begin{align*}"
g << endl;
g<< "K_{14,2}(5;0)\\cong";
L1 = rsort (schurBetti(5,2,0))#(14,2)
apply(L1,i->(
	g<< "\\bS_{(";
	g<< concatenate{toString((i#0)#0), ",", toString((i#0)#1), ",", toString((i#0)#2)};
	g<< ")}^{";
	g<< toString(i#1);
	g<< "}\\oplus";
	))
g << endl;
g << "\\end{align*}"
g << endl;
g << endl;

g << "\\begin{align*}"
g << endl;
g<< "K_{5,0}(5;3)\\cong";
L1 = rsort (schurBetti(5,2,3))#(5,0)
apply(L1,i->(
	g<< "\\bS_{(";
	g<< concatenate{toString((i#0)#0), ",", toString((i#0)#1), ",", toString((i#0)#2)};
	g<< ")}^{";
	g<< toString(i#1);
	g<< "}\\oplus";
	))
g << endl;
g << "\\end{align*}"
g << endl;
g << endl;

g << "\\begin{align*}"
g << endl;
g<< "K_{6,0}(5;3)\\cong";
L1 = rsort (schurBetti(5,2,3))#(6,0)
apply(L1,i->(
	g<< "\\bS_{(";
	g<< concatenate{toString((i#0)#0), ",", toString((i#0)#1), ",", toString((i#0)#2)};
	g<< ")}^{";
	g<< toString(i#1);
	g<< "}\\oplus";
	))
g << endl;
g << "\\end{align*}"
g << endl;
g << endl;

g << "\\begin{align*}"
g << endl;
g<< "K_{7,0}(5;3)\\cong";
L1 = rsort (schurBetti(5,2,3))#(7,0)
apply(L1,i->(
	g<< "\\bS_{(";
	g<< concatenate{toString((i#0)#0), ",", toString((i#0)#1), ",", toString((i#0)#2)};
	g<< ")}^{";
	g<< toString(i#1);
	g<< "}\\oplus";
	))
g << endl;
g << "\\end{align*}"
g << endl;
g << endl;

g << "\\begin{align*}"
g << endl;
g<< "K_{8,0}(5;3)\\cong";
L1 = rsort (schurBetti(5,2,3))#(8,0)
apply(L1,i->(
	g<< "\\bS_{(";
	g<< concatenate{toString((i#0)#0), ",", toString((i#0)#1), ",", toString((i#0)#2)};
	g<< ")}^{";
	g<< toString(i#1);
	g<< "}\\oplus";
	))
g << endl;
g << "\\end{align*}"
g << endl;
g << endl;

g << "\\begin{align*}"
g << endl;
g<< "K_{9,0}(5;3)\\cong";
L1 = rsort (schurBetti(5,2,3))#(9,0)
apply(L1,i->(
	g<< "\\bS_{(";
	g<< concatenate{toString((i#0)#0), ",", toString((i#0)#1), ",", toString((i#0)#2)};
	g<< ")}^{";
	g<< toString(i#1);
	g<< "}\\oplus";
	))
g << endl;
g << "\\end{align*}"
g << endl;
g << endl;

g << "\\begin{align*}"
g << endl;
g<< "K_{4,1}(5;3)\\cong";
L1 = rsort (schurBetti(5,2,3))#(4,1)
apply(L1,i->(
	g<< "\\bS_{(";
	g<< concatenate{toString((i#0)#0), ",", toString((i#0)#1), ",", toString((i#0)#2)};
	g<< ")}^{";
	g<< toString(i#1);
	g<< "}\\oplus";
	))
g << endl;
g << "\\end{align*}"
g << endl;
g << endl;

g << "\\begin{align*}"
g << endl;
g<< "K_{5,1}(5;3)\\cong";
L1 = rsort (schurBetti(5,2,3))#(5,1)
apply(L1,i->(
	g<< "\\bS_{(";
	g<< concatenate{toString((i#0)#0), ",", toString((i#0)#1), ",", toString((i#0)#2)};
	g<< ")}^{";
	g<< toString(i#1);
	g<< "}\\oplus";
	))
g << endl;
g << "\\end{align*}"
g << endl;
g << endl;

g << "\\begin{align*}"
g << endl;
g<< "K_{6,1}(5;3)\\cong";
L1 = rsort (schurBetti(5,2,3))#(6,1)
apply(L1,i->(
	g<< "\\bS_{(";
	g<< concatenate{toString((i#0)#0), ",", toString((i#0)#1), ",", toString((i#0)#2)};
	g<< ")}^{";
	g<< toString(i#1);
	g<< "}\\oplus";
	))
g << endl;
g << "\\end{align*}"
g << endl;
g << endl;

g << "\\begin{align*}"
g << endl;
g<< "K_{7,1}(5;3)\\cong";
L1 = rsort (schurBetti(5,2,3))#(7,1)
apply(L1,i->(
	g<< "\\bS_{(";
	g<< concatenate{toString((i#0)#0), ",", toString((i#0)#1), ",", toString((i#0)#2)};
	g<< ")}^{";
	g<< toString(i#1);
	g<< "}\\oplus";
	))
g << endl;
g <<"\\end{align*}"
g << endl;
g << endl;

g << "\\begin{align*}"
g << endl;
g<< "K_{8,1}(5;3)\\cong";
L1 = rsort (schurBetti(5,2,3))#(8,1)
apply(L1,i->(
	g<< "\\bS_{(";
	g<< concatenate{toString((i#0)#0), ",", toString((i#0)#1), ",", toString((i#0)#2)};
	g<< ")}^{";
	g<< toString(i#1);
	g<< "}\\oplus";
	))
g << endl;
g << "\\end{align*}"
g << endl;
g << endl;

L1 = rsort apply((schurBetti(5,2,0))#(14,1),i->i#0)
unique delete(,apply(L1, i->(
	c = 0;
	c = #select(L1,j->(j==i));
	if c>1 then i
	)))

L1 = rsort apply((schurBetti(5,2,0))#(15,1),i->i#0)
unique delete(,apply(L1, i->(
	c = 0;
	c = #select(L1,j->(j==i));
	if c>1 then i
	)))

L1 = rsort apply((schurBetti(5,2,0))#(13,2),i->i#0)
unique delete(,apply(L1, i->(
	c = 0;
	c = #select(L1,j->(j==i));
	if c>1 then i
	)))

L1 = rsort apply((schurBetti(5,2,0))#(14,2),i->i#0)
unique delete(,apply(L1, i->(
	c = 0;
	c = #select(L1,j->(j==i));
	if c>1 then i
	)))

needsPackage "SchurVeronese"

apply({5,6,7,8,9},k->(
L1 = rsort apply((schurBetti(5,2,3))#(k,0), i->i#0);
{(k,0),unique delete(,apply(L1, i->(
	c = 0;
	c = #select(L1,j->(j==i));
	if c>1 then i
	)))}
))


apply({4,5,6,7,8},k->(
L1 = rsort apply((schurBetti(5,2,3))#(k,1), i->i#0);
{(k,1),unique delete(,apply(L1, i->(
	c = 0;
	c = #select(L1,j->(j==i));
	if c>1 then i
	)))}
))

