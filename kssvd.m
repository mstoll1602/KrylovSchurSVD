function [s,UU,VV,residuals]=kssvdfinal4(A,v0)

% Lanczos bidiagonalization with upper triangular matrix

%% Initialize all the variables
lockindex = 0;
[n,m] = size(A); % dimension of A
sigma = 'PS';
ldes =10; % number of wanted eigenvalues
l = ldes+3; % adjusted desired number of singular values
ldim = 20;% dimension of the search space
ldim_org = ldim;
k = 1;
maxit = min(1000,n); % maximal iteration number
tol1  = 1e-10; % tolerance for the singular values
iter  = 1;
SVTol = min(sqrt(eps),tol1);  % Tolerance to determine when a singular vector has converged.
reorth = 'one';
mprod = 0;
tau = 1.5


s = zeros(ldes,1);
VV = zeros(m,ldim);
UU = zeros(n,ldim);

V = zeros(m,ldim+1);  
U = zeros(n,ldim+1);
B = zeros(ldim,ldim);

% forward sequence
% tmp=rand(n,1);
V(:,k) = v0;
betakplus1 = norm(V(:,k));
V(:,k) = V(:,k)/betakplus1;

% adjoint sequence
U(:,k) = A*V(:,k);
alphak = norm(U(:,k));
mprod = mprod+1;

% check for linear dependency
if alphak <= SVTol
            disp('hier')
   U(:,k) = randn(size(U,1),1);
   U(:,k) = orthog(U(:,k),U(:,1:k-1));
   U(:,k) = U(:,k)/norm(U(:,k));
   alphak = 0;
else
   U(:,k) = U(:,k)/alphak;
end
B(1,1) = alphak;




%% initial Lanczos   bidiagonalisation of dimension ldim
for k = 1:l
    % one step of the Lanczos bidiagonalisation procedure
    V(:,k+1) = A'*U(:,k)-alphak*V(:,k);
    V(:,k+1) = orthog(V(:,k+1),V(:,1:k));
    betakplus1 = norm(V(:,k+1));
    if betakplus1 <= SVTol
        disp('beta=0')
        V(:,k+1) = randn(size(V,1),1);
        V(:,k+1) = orthog(V(:,k+1),V(:,1:k));
        V(:,k+1) = V(:,k+1)/norm(V(:,k+1));
        betakplus1 = 0;
    else
        V(:,k+1) = V(:,k+1)/betakplus1;
    end
    mprod = mprod+1;
    
    U(:,k+1) = (V(:,k+1)'*A')'-betakplus1*U(:,k);
    U(:,k+1) = orthog(U(:,k+1),U(:,1:k));
    alphakplus1 = norm(U(:,k+1));
    % check for linear dependency
    if alphakplus1 <= SVTol
        disp('alpha=0')
        U(:,k+1) = randn(size(U,1),1);
        U(:,k+1) = orthog(U(:,k+1),U(:,1:k));
        U(:,k+1) = U(:,k+1)/norm(U(:,k+1));
        alphakplus1 = 0;
    else
        U(:,k+1) = U(:,k+1)/alphakplus1;
    end
    mprod = mprod+1;
    
    B(k+1,k+1) = alphakplus1;
    B(k,k+1) = betakplus1;
    alphak = alphakplus1;
end



% % % Adjust K to include more vectors as the number of vectors converge.
% % Len_res = length(find(residuals(1:K_org)' < SVTol*Smax));
% % K = max(K, K_org + Len_res); if K > Bsz - 3, K = Bsz - 3; end

% Change the dimension to have less matrix multiplications!
%% MAIN ITERATION
tt = l+1;
while  iter <= maxit
    %%  enlarge decomposition from l to ldim
    for k = tt:ldim
        % one step of the Lanczos bidiagonalisation procedure
        V(:,k+1) = A'*U(:,k)-alphak*V(:,k);
        V(:,k+1) = orthog(V(:,k+1),V(:,1:k));
%         V(:,k+1) = orthog(V(:,k+1),V(:,lockindex+1:k));
        betakplus1 = norm(V(:,k+1));
        if betakplus1 <= SVTol
            disp('beta=0')
            V(:,k+1) = randn(size(V,1),1);
            V(:,k+1) = orthog(V(:,k+1),V(:,1:k));
            V(:,k+1) = V(:,k+1)/norm(V(:,k+1));
            betakplus1 = 0;
        else
            V(:,k+1) = V(:,k+1)/betakplus1;
        end
        mprod = mprod+1;
        U(:,k+1) = (V(:,k+1)'*A')'-betakplus1*U(:,k);
%         U(:,k+1) = A*V(:,k+1)-betakplus1*U(:,k);
        if  strcmp(reorth,'two')
            U(:,k+1) = orthog(U(:,k+1),U(:,1:k));
%             U(:,k+1) = orthog(U(:,k+1),U(:,lockindex+1:k));
        end
        alphakplus1 = norm(U(:,k+1));
        % check for linear dependency
        if alphakplus1 <= SVTol
            disp('alpha')
            U(:,k+1) = randn(size(U,1),1);
            U(:,k+1) = orthog(U(:,k+1),U(:,1:k));
            U(:,k+1) = U(:,k+1)/norm(U(:,k+1));
            alphakplus1 = 0;
        else
            U(:,k+1) = U(:,k+1)/alphakplus1;
        end
        mprod = mprod+1;

        B(k+1-lockindex,k+1-lockindex) = alphakplus1;
        B(k-lockindex,k+1-lockindex) = betakplus1;
        alphak = alphakplus1;
    end

    % Computing the SVD of B
    [P,S,Q] = svd(B(1:k-lockindex,1:k-lockindex));
    ss = diag(S);
%     pause
   
    % Updating the bidiagonal factorisation
    U(:,lockindex+1:k) = U(:,lockindex+1:k)*P;    
    V(:,lockindex+1:k) = V(:,lockindex+1:k)*Q;    
   
    % residual vector    pause
    res_vec = P(k-lockindex,1:k-lockindex);


    if strcmp(sigma, 'LS')
        [sorted, index] = sort(ss,'descend');
    elseif strcmp(sigma, 'SS')
        [sorted, index] = sort(ss,'ascend');
    else 
%         abs(ss-tau)
        [sorted, index] = sort(abs(ss-tau),'ascend');
%         ss(index)
%         pause
%         pause
    end

%% Computing the stopping criteria and locking
    stopping = abs(betakplus1*res_vec(index(1:l-lockindex)));   
    indexstop = find(stopping<tol1);
    lockplus = length(indexstop);
    if lockplus > 0
        s(lockindex+1:lockindex+lockplus) = ss(index(indexstop));
        UU(:,lockindex+1:lockindex+lockplus)  = U(:,lockindex+index(indexstop));
        VV(:,lockindex+1:lockindex+lockplus)  = V(:,lockindex+index(indexstop));
%         norm(A*VV(:,lockindex+1:lockindex+lockplus)-UU(:,lockindex+1:lockindex+lockplus)*diag(s(lockindex+1:lockindex+lockplus)))
        
    end
    [rubbish,index2] = sort(stopping,'ascend');
    
    S2 = S(index(index2(lockplus+1:l-lockindex)),index(index2(lockplus+1:l-lockindex)));

    V(:,lockindex+1:l) = V(:,lockindex+index(index2(1:l-lockindex)));
    U(:,lockindex+1:l) = U(:,lockindex+index(index2(1:l-lockindex)));
    U(:,l+1) = U(:,k+1);
    V(:,l+1) = V(:,k+1);
    res_vec = res_vec(index(index2(lockplus+1:l-lockindex)));
    lockindex = lockindex + lockplus; % increase lockindex
    residuals(iter,:) = abs([tol1*ones(1,lockindex) res_vec ]);

%     tt = max(tt,ldes+lockplus);   
    tt = max(tt,ldes+lockindex);
%     ldim = ldim_org+lockindex;
    
%     if lockplus >= ldes
    if lockindex >= ldes
        break;
    end;
 
%% transforming the residual vector to a multiple of e_{k}^{T}
    v = housegen(res_vec',l-lockindex);
    W = speye(l-lockindex)-v*v';
    temp1 = res_vec*W;
    betakplus1 = temp1(end)*betakplus1;
    C_k = W*S2*W';

%% transforming the dense matrix C_k to bidiagonal form
    [B,Q_l,P_l]=bidiagredrowup(C_k);
    P_l = W*P_l;
    Q_l = W*Q_l;
    U(:,lockindex+1:l) = U(:,lockindex+1:l)*P_l;
    V(:,lockindex+1:l) = V(:,lockindex+1:l)*Q_l;   

    U(:,l+1) = A*V(:,l+1)-betakplus1*U(:,l);
    if  strcmp(reorth,'two')
        U(:,l+1) = orthog(U(:,k+1),U(:,1:k));
%         U(:,l+1) = orthog(U(:,k+1),U(:,lockindex+1:k));
    end
    alphakplus1 = norm(U(:,l+1));
    % check for linear dependency
    if alphakplus1 <= SVTol
        disp('hier')
        U(:,l+1) = randn(size(U,1),1);
        U(:,l+1) = orthog(U(:,l+1),U(:,1:l));
        U(:,l+1) = U(:,l+1)/norm(U(:,l+1));
        alphakplus1 = 0;
    else
        U(:,l+1) = U(:,l+1)/alphakplus1;
    end
    mprod = mprod+1;
    B(l+1-lockindex,l+1-lockindex) = alphakplus1;
    B(l-lockindex,l+1-lockindex) = betakplus1;
    alphak = alphakplus1;
    clear lockplus;
%     lockindex
    iter = iter +1 ;
    
end
iter
mprod

%% Orthogonalisation function
function Y = orthog(Y,X)
% Orthogonalize vectors Y against vectors X.

if size(X,2) < size(Y,2), dotY = X'*Y; else dotY = (Y'*X)'; end
Y = Y - X*dotY; 