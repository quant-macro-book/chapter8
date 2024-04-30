% LumpySSWagemy.m
% solve for the steady state of Khan and Thomas (2003)
% based on the section 3.3 of Khan and Thomas (2009)
% w is bisected.
% January 5 2010

clear all;

t = cputime;

% convergences
global critv critg;
critbw  = 1e-4;
critv   = 1e-5;
critg   = 1e-10;

% parameters, taken from p.341 of Khan and Thomas (2003)
global GAMY BETA DELTA THETA NU ETA B kSS;

GAMY  = 1.016;
BETA  = 0.954;
DELTA = 0.060;
THETA = 0.325;
NU    = 0.580;
ETA   = 3.614;
RHO   = 0.9225;
SIGMA = 0.0134;
B     = 0.002;

ykSS = (GAMY-BETA*(1-DELTA))/BETA/THETA;
ckSS = ykSS + (1-GAMY-DELTA);
ycSS = ykSS/ckSS;
nSS = NU/ETA*ycSS;
kSS = (ykSS*nSS^(-NU))^(1/(THETA-1));

kbounds = [0.1 3.0];
nk = 201;
knotsk = logspace(log(kbounds(1) + -1.0*kbounds(1) + 1.0)/log(10.0), log(kbounds(2) + -1.0*kbounds(1) + 1.0)/log(10.0), nk)';
knotsk = knotsk + ones(size(knotsk))*(kbounds(1) - 1.0);
%knotsk = polygrid(nk,kbounds(1),kbounds(2));
%Phi = polybas(knotsk,nk,kbounds(1),kbounds(2));

% ss value of z
znow = 1.0;

% bisection method
wL = 1.0;
wLnew = LumpySSmy(wL,znow,knotsk,nk); %,kbounds,Phi);
BL = wL-wLnew;
disp('  ');

wH = 1.2;
wHnew = LumpySSmy(wH,znow,knotsk,nk); %,kbounds,Phi);
BH = wH-wHnew;
disp('  ');

diff = 1e+4;
iter = 0;

while (diff>critbw)

    w0 = (wL+wH)/2;
    wnew = LumpySSmy(w0,znow,knotsk,nk); %,kbounds,Phi);
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

[w v kw e0 xi theta kdist] = LumpySSmy(w0,znow,knotsk,nk); %,kbounds,Phi);
werr = w-w0;
p = ETA/w;
disp('  ');
disp(sprintf('  Total computation time = %6.8f',cputime-t));
disp('  Equilibrium prices');
disp('    Wage      Price     Wage error');
disp([w p werr]);

% plot policy function
gvec = zeros(nk,1);
for ik = 1:nk
    know = knotsk(ik);
    alpha = xi(ik)/B;
    kp = alpha*kw + (1-alpha)*(1-DELTA)*know/GAMY;
    gvec(ik) = kp;
end

figure;
subplot(211);
plot(knotsk,gvec);
xlim([0.5 2.0]);

figure;
subplot(211);
plot(knotsk,xi/B,'b-','LineWidth',2.0);
title('資本調整確率');
xlim([0.7 1.6]);
box on;
grid on;
set(gca,'Fontsize',12)
subplot(212);
bar(kdist,theta,'blue');
title('企業分布');
xlabel('k');
ylabel('\mu');
%xlim([0.7 1.15]);
xlim([0.7 1.6]);
%xlim([min(kdist) max(kdist)+0.1]);
ylim([0 0.3]);
box on;
grid on;
set(gca,'Fontsize',12)

print -depsc2 KT.eps

save LumpySS v kw e0 xi theta kdist;