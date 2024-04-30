function f = vfuncsp2(kp,mp,p,cmat,knotsk,knotsm)

global GAMY BETA;

rk = length(knotsk)-2;
rm = length(knotsm)-2;

ev = speva2(cmat,kp,mp,rk,rm,knotsk,knotsm);
% ev = diag(ev);
f = -GAMY.*p.*kp + BETA.*ev;
f = -f; % for golden.m