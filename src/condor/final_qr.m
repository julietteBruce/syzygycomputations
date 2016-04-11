


%%% Loads Matrix From Given File %%%
sparse_data load(inFile)
T = spconvert(sparse_data);

%%% Computes QR Factorization %%%
[m,n] = size(T);
q=colamd(T, [n m]);
[Q,R,P] = qr(T(:,q));

%%% Computes Rank from QR Factorization %%%
tol = max(size(T))*eps*abs(R(1,1));
r = 0;
while ( abs(R(r+1,r+1)) >= tol & r < min(n,m)+1 )
r = r+1;
end

%%% Saves Rank to MatLab File %%%
save('test_2_5_7_2_multidegree_20_18_7.mat','r')
r
