% ceinner1my.m

clear all;

global critin critbp critg2;
critin = 1e-4;
critbp = 1e-10;
critg2 = 1e-10;

% parameters, taken from p.341 of Khan and Thomas (2003)
global GAMY BETA DELTA THETA NU ETA B kSS;

GAMY  = 1.0160;
BETA  = 0.9540;
DELTA = 0.060;
THETA = 0.3250;
NU    = 0.580;
ETA   = 3.6142;
% RHO   = 0.9225;
% SIGMA = 0.0134;
B     = 0.002;

ykSS = (GAMY-BETA*(1-DELTA))/BETA/THETA;
ckSS = ykSS + (1-GAMY-DELTA);
ycSS = ykSS/ckSS;
nSS = NU/ETA*ycSS;
kSS = (ykSS*nSS^(-NU))^(1/(THETA-1));

% aggregate shock process
% Z = [0.99 1.01]';
% nz = 2;
% Pi = [0.9 0.1;
%     0.1 0.9];
Z = [0.9328 0.9658 1.0000 1.0354 1.0720]';
nz = 5;
Pi = [0.8537    0.1377    0.0083    0.0002    0.0000;
    0.0344    0.8579    0.1035    0.0042    0.0001;
    0.0014    0.0690    0.8593    0.0690    0.0014;
    0.0001    0.0042    0.1035    0.8579    0.0344;
    0.0000    0.0002    0.0083    0.1377    0.8537];

kbounds = [0.1 3.0];
nk = 25;
% knotsk = logspace(log(kbounds(1))/log(10), log(kbounds(2))/log(10), nk)';
knotsk = logspace(log(kbounds(1) + -1.0*kbounds(1) + 1.0)/log(10.0), log(kbounds(2) + -1.0*kbounds(1) + 1.0)/log(10.0), nk)';
knotsk = knotsk + ones(size(knotsk))*(kbounds(1) - 1.0);
rk = nk-2;

% mbounds = [0.9*kSS 1.1*kSS];
mbounds = [0.8 1.2];
nm = 5;
knotsm = linspace(mbounds(1),mbounds(2),nm)';
rm = nm-2;

% realization of shocks are the same among outer loops
simT = 2500;
cumPi= cumsum(Pi')';
izvec = zeros(simT,1);
izvec(1) = 1;
for tt = 1:simT-1
    cumPivec = cumPi(izvec(tt),:);
    izvec(tt+1) = sum(rand-cumPivec>=0);
    izvec(tt+1) = min(izvec(tt+1)+1,nz);
end

load ../Xpa/lumpyeqxparesult v mpmat pmat;
[BetaK,Betap] = mat2coef(mpmat,pmat,knotsm);
% [Kvec1 Kpvec1 Cvec1] = ceouter1my(v,BetaK,Betap,knotsk,knotsm,Z,Pi,izvec);
[Yvec1 Ivec1 Cvec1 Nvec1 Wvec1 Zvec1 Kvec1 Kpvec1] = lumpyeqoutersim(v,mpmat,knotsk,knotsm,Z,Pi,izvec);
Rvec1 = GAMY/BETA*Cvec1(2:end)./Cvec1(1:end-1);
data1 = [Yvec1 Ivec1 Cvec1 Nvec1 Wvec1];
data1 = [data1(2:end,:) Rvec1];
sd1 = 100*std(hpfilter(log(data1),100));
corr1 = corr(hpfilter(log(data1),100));
corr1 = corr1(1,:);
% sd1(6) = 100*std(hpfilter(log(Rvec1),100));
disp([sd1; corr1]);

load ../KS/lumpyeqksresult v BetaK Betap;
mpmat = coef2mat(BetaK,Betap,knotsm);
[Yvec2 Ivec2 Cvec2 Nvec2 Wvec2 Zvec2 Kvec2 Kpvec2] = lumpyeqouter(v,BetaK,Betap,knotsk,knotsm,Z,Pi,izvec);
% [Kvec2 Kpvec2 Cvec2] = ceouter1mysim(v,mpmat,knotsk,knotsm,Z,Pi,izvec);
Rvec2 = GAMY/BETA*Cvec2(2:end)./Cvec2(1:end-1);
data2 = [Yvec2 Ivec2 Cvec2 Nvec2 Wvec2];
data2 = [data2(2:end,:) Rvec2];
sd2 = 100*std(hpfilter(log(data2),100));
corr2 = corr(hpfilter(log(data2),100));
corr2 = corr2(1,:);
% sd1(6) = 100*std(hpfilter(log(Rvec1),100));
disp([sd2; corr2]);

% time series
tt = [1001:1:2000]';
figure;
% subplot(311);
plot(tt,Yvec1(tt),'b-',tt,Yvec2(tt),'g-');
legend('Xpa','KS');
print -depsc2 lumpyyvec.eps;
figure;
% subplot(312);
plot(tt,Ivec1(tt),'b-',tt,Ivec2(tt),'g-');
print -depsc2 lumpyivec.eps;
figure;
% subplot(313);
plot(tt,Cvec1(tt),'b-',tt,Cvec2(tt),'g-');
legend('Xpa','KS');
print -depsc2 lumpycvec.eps;
% print -depsc2 lumpysimdata1.eps;
figure;
% subplot(311);
plot(tt,Nvec1(tt),'b-',tt,Nvec2(tt),'g-');
print -depsc2 lumpynvec.eps;
figure;
% subplot(312);
plot(tt,Wvec1(tt),'b-',tt,Wvec2(tt),'g-');
legend('Xpa','KS');
print -depsc2 lumpywvec.eps;
figure;
% subplot(313);
plot(tt,Rvec1(tt),'b-',tt,Rvec2(tt),'g-');
print -depsc2 lumpyrvec.eps;
% print -depsc2 lumpysimdata2.eps;
figure;
plot(tt,Kvec1(tt),'b-',tt,Kvec2(tt),'g-');
print -depsc2 lumpykvec.eps;

% forecasting rules
load ../Xpa/lumpyeqxparesult mpmat pmat;
mpmat1 = mpmat;
pmat1 = pmat;
clear mpmat pmat;
load ../KS/lumpyeqksresult v BetaK Betap;
[mpmat2 pmat2] = coef2mat(BetaK,Betap,knotsm);

figure;
for iz = 1:5
    subplot(3,2,iz);
    plot(knotsm,mpmat1(:,iz),'bo-',knotsm,mpmat2(:,iz),'g-');
    xlim([knotsm(1) knotsm(end)]);
    ylim([0.7 1.2]);
end
print -depsc2 lumpympmat.eps

figure;
for iz = 1:5
    subplot(3,2,iz);
    plot(knotsm,pmat1(:,iz),'bo-',knotsm,pmat2(:,iz),'g-');
    xlim([knotsm(1) knotsm(end)]);
    ylim([2.5 4.0]);
end
print -depsc2 lumpypmat.eps