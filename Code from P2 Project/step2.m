%% Input: The path where the multigraded rank data is stored AND
%% a file name for output. These need to be manually changed.
%% Output: A .txt file that contains the multigraded Hilbert Series.
%% WARNING: This only works for P2 and the (b,d) need to be manually adj.
cd('./mgBettiData/P2/6_3')
files = dir('*.txt')
for file = files'
    A = importdata(file.name)
    [pathstr,name,ext] = fileparts(file.name) 
     B = [];
     for i=1:size(A,1)
         P = unique(perms(A(i,1:3)),'rows');
         [n,m] = size(P);
         B = [B;[P,A(i,4)*ones(n,1)]];
     end
     B=A;
     B = B';
     [n,m] = size(B);
     C = [B(1,1:(m-1)); B(2,1:(m-1)); B(3,1:(m-1)); B(4,1:(m-1))];
     formatSpec = '(t_0^%d*t_1^%d*t_2^%d)*%d+'
     fileID = fopen(strcat('../../../mgHilbertSeries/P2/2_6_3-',sprintf(name),'.txt'),'w');
     fprintf(fileID,formatSpec,C)
     formatSpec = '(t_0^%d*t_1^%d*t_2^%d)*%d'
     fprintf(fileID,formatSpec,B(:,m))
     fclose(fileID);
 end

% A = importdata('./mgBettiData/P2/4_3/12_1.txt')
% 
% B = [];
% for i=1:size(A,1)
%     P = unique(perms(A(i,1:3)),'rows');
%     [n,m] = size(P);
%     B = [B;[P,A(i,4)*ones(n,1)]];
% end
% B=A;
% B = B';
% [n,m] = size(B);
% C = [B(1,1:(m-1)); B(2,1:(m-1)); B(3,1:(m-1)); B(4,1:(m-1))];
% formatSpec = '(t_0^%d*t_1^%d*t_2^%d)*%d+'
% fileID = fopen('mgHilbertSeries/P2/2_4_3-12_1.txt','w');
% fprintf(fileID,formatSpec,C)
% formatSpec = '(t_0^%d*t_1^%d*t_2^%d)*%d'
% fprintf(fileID,formatSpec,B(:,m))
% fclose(fileID);