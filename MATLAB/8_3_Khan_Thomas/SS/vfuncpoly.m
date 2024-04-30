function [f0 df0 d2f0] = vfuncsp(kp,cfz,knotsk,kbounds)

global GAMY BETA;

nk = size(knotsk,1);

[Phi dPhi d2Phi] = polybas(kp,nk,kbounds(1),kbounds(end));
EV   = Phi*cfz;
EDV  = dPhi*cfz;
ED2V = d2Phi*cfz;

f0   = -GAMY*kp + BETA*EV;
df0  = -GAMY + BETA*EDV;
d2f0 = BETA*ED2V;

f0 = -f0; % for golden.m