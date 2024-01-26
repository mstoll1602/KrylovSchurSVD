function v = housegen(a,n)

% This is a complex version of the house function
% H=I-uu^T is computed such that Ha=nu e_1.
%
% Martin Stoll 2004

v = a;
nu = norm(a);
if nu == 0, 
    v(n)=sqrt(2);
    return;
end;
if (v(n)~=0),
    rho = conj(v(n))/abs(v(n));
else 
    rho = 1;
end;
v = (rho/nu)*v;
v(n) = 1+v(n);
v = v/sqrt(v(n));
nu = -conj(rho)*nu;
return    