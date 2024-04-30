function f0 = nfunc(n0,k0,kp,z0);

global GAMY BETA DELTA THETA NU ETA B;

yterm = z0.*k0.^THETA.*n0.^NU + (1-DELTA).*k0;
c0 = yterm - GAMY.*kp;

% first-order condition of n
f0 = -NU.*z0.*k0.^THETA.*n0.^(NU-1) + ETA.*c0;