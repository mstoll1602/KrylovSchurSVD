function [H,Q,P]=bidiagredrowup(A);

[n,nn]=size(A);
Q = eye(n);
P = eye(n);
oldA = A;

for i=1:n-1
    
    v = housegen(A(n-i+1,1:n-i+1)',n-i+1);
    W1 = speye(n-i+1)-v*v';
    W  = speye(n);
    W(1:n-i+1,1:n-i+1) = W1;
    A = A*W;
    Q = Q*W;

    
    v = housegen(A(1:n-i,n-i+1),n-i);
    W  = speye(n);
    W1 = eye(n-i)-v*v';
    W(1:n-i,1:n-i) = W1;
    P = P*W;
    A = W'*A;
    
end;

H=A;
index = find(abs(H)<1e-10);
H(index) =0;
% % disp('hurrrrrrrrrrrrrrrrrrrra')
% % P'*oldA*Q

