function [score]=IMC(A,X,Y,r)
%[M,N]=size(A);
[~,m]=size(X);
[~,n]=size(Y);
W=rand(m,r);
H=rand(n,r);
k=1;
while k<1000
    H=H.*(Y'*A'*(X*W))./(Y'*(Y*H)*(W'*(X')*(X*W))+H);
    W=W.*(X'*A*(Y*H))./(X'*(X*W)*(H'*(Y')*(Y*H))+W);
    %error(k)=0.5*norm(A-X*W*H'*Y','fro')^0.5;%+0.5*norm(W,'fro')^2+0.5*norm(H,'fro')^2;
    k=k+1;
end
score=(X*W)*(Y*H)';
end

