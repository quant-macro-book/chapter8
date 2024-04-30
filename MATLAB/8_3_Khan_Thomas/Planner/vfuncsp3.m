function f0 = vfuncsp3(n0,cfz,knotsk,k0,z0)

global GAMY BETA DELTA THETA NU ETA;

rk = length(knotsk)-2;

yterm = z0.*k0.^THETA.*n0.^NU + (1-DELTA).*k0;

% first-order condition of n
c0 = (NU/ETA).*z0.*k0.^THETA.*n0.^(NU-1);
kp = (yterm-c0)/GAMY;
vp = speva(cfz,kp,rk,knotsk);

util = log(c0) + ETA*(1-n0);
f0 = util + BETA*vp;
f0 = -f0; % for golden.m