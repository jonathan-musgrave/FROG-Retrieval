function [Et1, Z] = MinZerr_sd(Esig, Et, dZ)

p1 = MinZerrKern_sd(Esig,Et,dZ);

p = polyder(p1);

r = roots(p);

X = r(find(imag(r) == 0))';


Z1 = polyval(p1,X);

[Z, minZIndx] = min(Z1);

% Assign eps value to calculated Z values that are lower
Zmin = eps * p1(end);

if Z < Zmin
    Z = Zmin;
end

Z = sqrt(Z);
X = X(minZIndx);

Et1 = Et + X * dZ;