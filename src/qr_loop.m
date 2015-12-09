function [] = qr_loop(x)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
cd(char(x));
files = dir('*.dat');
L = [0 0];
for file = files'
    S = load(file.name);
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
    r;
    L = [L; r t];  
end
save('ranks.mat', 'L')

end

