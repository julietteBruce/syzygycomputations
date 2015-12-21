function [] = qr_single(input,output)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
S = load(input);
T = spconvert(S);
tic
[Q,R,P] = qr(T); 
[m,n] = size(T);
tol = max(size(T))*eps*abs(R(1,1));
r = 0;
while (r<min(m,n) & abs(R(r+1,r+1)) >= tol)
    r = r+1;
end
t=toc;
outf = fopen(output,'w');
fprintf(outf,'%d %d\n',r,t);
end
