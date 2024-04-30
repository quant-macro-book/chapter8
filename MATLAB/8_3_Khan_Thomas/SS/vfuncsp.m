function f0 = vfuncsp(k0,cfz,knotsk)

global GAMY BETA;

r = length(knotsk)-2;
v0 = speva(cfz,k0,r,knotsk);
f0 = -GAMY*k0 + BETA*v0;
f0 = -f0; % for golden.m