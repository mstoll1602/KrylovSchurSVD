clear all
clc
cd examples
A=mmread('ash958.mtx');
% A=mmread('well1850.mtx');
% A=mmread('well1033.mtx');
% load lp_d6cube
% load aircraft 
% load neos1
% load lp_nug20
% load nug08-3rd
% A = Problem.A';
[n,m] = size(A);
if m > n, A=A';end
% load G
% A = G;
% load testsvd1
cd ..


v0 = randn(size(A,2),1);

% parameters for irlba
opts.k=10;
opts.tol=1e-10;
opts.M_B=20;
opts.V0 =v0;
opts.sigma='LS';
opts.AUG = 'RITZ';
% opts.AUG = 'HARM';
% parameters for svds
options.tol = 1e-10;
options.p = opts.M_B;
sigma =1.5;
% sigma ='L';

tic
[s1,UU,VV,residuals]=kssvd(A,v0);
% check the quality of the truncated SVD
for i=1:size(s1)
res1(i)=norm(A*VV(:,i)-s1(i)*UU(:,i));
end
norm(A*VV(:,1:size(s1))-UU(:,1:size(s1))*diag(s1))
toc

disp('##########')
[UU,s2,VV]=irlba(A,opts);
for i=1:size(s2)
norm(A*VV(:,i)-s2(i,i)*UU(:,i))
end
s2 = diag(s2);

disp('##########')
[UU,s3,VV,flag]=svds(A,opts.k,sigma,options);
% check the quality of the truncated SVD
for i=1:size(s3)
res2(i)=norm(A*VV(:,i)-s3(i,i)*UU(:,i));
end
norm(A*VV(:,1:size(s3))-UU(:,1:size(s3))*s3)
disp('##########')
s3=diag(s3);
res1(1:opts.k) 
format long e
% [sort(s1(1:opts.k),'descend') s2 s3]
[sort(s1(1:opts.k),'descend') res1(1:opts.k)' s3 res2(1:opts.k)']
sort(s3,'descend')
sort(s1,'descend')

