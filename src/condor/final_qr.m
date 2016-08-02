
function final_qr(inFile,outFile)


function y = qr_rank(inFile)
%%% Loads Matrix From Given File %%%
sparse_data = load(inFile);
T = spconvert(sparse_data);

%%% Computes QR Factorization %%%
#[m,n] = size(T);
#q=colamd(T, [n m]);
#T=T(:,q);
T(:,~any(T,1)) = [];
[m,n] = size(T);
[Q,R,P] = qr(T);


%%% Computes Rank from QR Factorization %%%
tol = max(size(T))*eps*abs(R(1,1));
r = 0;
while (r < min(n,m) & abs(R(r+1,r+1)) >= tol)
    r = r+1;
end

%%% Saves Rank to MatLab File %%%
%save('test_2_5_7_2_multidegree_20_18_7.mat','r')
y=r;
end

outFd = fopen(outFile,'w');

res = qr_rank(inFile);
fprintf(outFd,'%d',res);

fclose(outFd);



end