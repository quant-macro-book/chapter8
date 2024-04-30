function [T dT d2T] = polybas(xnodes,m,a,b)

n = length(xnodes);
z = 2.*(xnodes-a)./(b-a) - ones(n,1);
dzdx = 2/(b-a);
T = zeros(n,m);
dT = zeros(n,m);
d2T = zeros(n,m);

T1       = z;
T2       = 2.*z.*z - ones(n,1);
T(:,1)   = T1;
dT(:,1)  = ones(n,1);
d2T(:,1) = zeros(n,1);
T(:,2)   = T2;
dT(:,2)  = 4.*z;
d2T(:,2) = 4.*ones(n,1);

for i=3:m
    T(:,i) = 2.*z.*T(:,i-1) - T(:,i-2);
    dT(:,i) = 2.*T(:,i-1) + 2.*z.*dT(:,i-1) - dT(:,i-2);
    d2T(:,i) = 2.*dT(:,i-1) + 2.*dT(:,i-1) + 2.*z.*d2T(:,i-1) - d2T(:,i-2);
end

T   = [ones(n,1) T(:,1:(m-1))];
dT  = dzdx.*[zeros(n,1) dT(:,1:(m-1))];
d2T = dzdx.*dzdx.*[zeros(n,1) d2T(:,1:(m-1))];