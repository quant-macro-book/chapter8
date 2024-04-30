% LumpySSWagemy.m
% solve for the steady state of Khan and Thomas (2003)
% based on the section 3.3 of Khan and Thomas (2009)
% w is bisected.
% January 5 2010

clear all;

t = cputime;

% convergences
global critv critg;
critbw  = 1e-5;
critv   = 1e-6;
critg   = 1e-10;

% parameters, taken from p.341 of Khan and Thomas (2003)
global GAMY BETA DELTA THETA NU ETA B;

GAMY  = 1.016;
BETA  = 0.954;
DELTA = 0.060;
THETA = 0.325;
NU    = 0.580;
ETA   = 3.614;
RHO   = 0.9225;
SIGMA = 0.0134;
B     = 0.002;

kbounds = [0.1 3.0];
nk = 50;
knotsk = logspace(log(kbounds(1) + -1.0*kbounds(1) + 1.0)/log(10.0), log(kbounds(2) + -1.0*kbounds(1) + 1.0)/log(10.0), nk)';
knotsk = knotsk + ones(size(knotsk))*(kbounds(1) - 1.0);

rk = nk-2;
T = spbas(rk,knotsk);
invT = inv(T);

% ss value of z
znow = 1.0;

% bisection method
wL = 1.0;
wLnew = LumpySSmy(wL,znow,knotsk,nk);
BL = wL-wLnew;
disp('  ');

wH = 1.2;
wHnew = LumpySSmy(wH,znow,knotsk,nk);
BH = wH-wHnew;
disp('  ');

diff = 1e+4;
iter = 0;

while (diff>critbw)

    w0 = (wL+wH)/2;
    wnew = LumpySSmy(w0,znow,knotsk,nk);
    B0 = w0-wnew; % g(w) = w-f(w)

    if B0*BL>0
        wL = w0;
        BL = B0;
    else
        wH = w0;
    end

    diff = wH-wL;
    iter = iter + 1;
    s = sprintf('  bisection %2d,  wH-wL = %6.8f',iter,diff);
    disp(s);
    disp('  ');

end

[w v kw e0 xi theta kdist] = LumpySSmy(w0,znow,knotsk,nk);
werr = w-w0;
p = ETA/w;
disp('  ');
disp(sprintf('  Total computation time = %6.8f',cputime-t));
disp('  Equilibrium prices');
disp('    Wage      Price     Wage error');
disp([w p werr]);

figure;
subplot(211);
plot(knotsk,xi/B);
xlim([0.7 1.6]);
subplot(212);
plot(kdist,theta);
xlim([min(kdist) max(kdist)]);
ylim([0 0.3]);

save LumpySS v kw e0 xi theta kdist;