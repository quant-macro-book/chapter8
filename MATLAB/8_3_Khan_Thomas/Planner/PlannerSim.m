% PlannerSim.m
% solve the control model and generate simulated data
% January 5 2010

clear all;

% parameters, taken from p.341 of Khan and Thomas (2003)
global GAMY BETA DELTA THETA NU ETA B;

GAMY  = 1.016;
BETA  = 0.954;
DELTA = 0.060;
THETA = 0.325;
NU    = 0.580;
ETA   = 3.614;
% RHO   = 0.9225;
% SIGMA = 0.0134;
B     = 0.002;

% agg. shock process
Z = [0.9328 0.9658 1.0000 1.0354 1.0720]';
nz = 5;
Pi = [0.8537    0.1377    0.0083    0.0002    0.0000;
    0.0344    0.8579    0.1035    0.0042    0.0001;
    0.0014    0.0690    0.8593    0.0690    0.0014;
    0.0001    0.0042    0.1035    0.8579    0.0344;
    0.0000    0.0002    0.0083    0.1377    0.8537];

kbounds = [0.1 3.0];
nk = 25;
knotsk = logspace(log(kbounds(1) + -1.0*kbounds(1) + 1.0)/log(10.0), log(kbounds(2) + -1.0*kbounds(1) + 1.0)/log(10.0), nk)';
knotsk = knotsk + ones(size(knotsk))*(kbounds(1) - 1.0);

simT = 2500;
burn = 0;
Kvec  = zeros(simT,1);
Zvec  = zeros(simT,1);
Kpvec = zeros(simT,1);
Yvec  = zeros(simT,1);
Ivec  = zeros(simT,1);
Cvec  = zeros(simT,1);
Nvec  = zeros(simT,1);
Wvec  = zeros(simT,1);

cumPi= cumsum(Pi')';
izvec = zeros(simT,1);
izvec(1) = 1;
for tt = 1:simT-1
    cumPivec = cumPi(izvec(tt),:);
    izvec(tt+1) = sum(rand-cumPivec>=0);
    izvec(tt+1) = min(izvec(tt+1)+1,nz);
end

Kvec(1) = 1;

load Planner n kp;

for time = 1:simT

    know = Kvec(time);
    iznow = izvec(time);
    znow = Z(iznow);
    
    nnow = interp1(knotsk,n(:,iznow),know);

    yterm = znow*know^THETA*nnow^NU + (1-DELTA)*know;

    % first-order condition of n
    cnow = (NU/ETA)*znow*know^THETA*nnow^(NU-1);
    kp = (yterm-cnow)/GAMY;
    
    % record aggregate variables  
    inow = GAMY*kp - (1-DELTA)*know;
    ynow = cnow + GAMY*kp - (1-DELTA)*know;
    pnow = 1/cnow;
    wnow = ETA*cnow;
    
    Kpvec(time,1)  = kp;
    Yvec(time,1)   = ynow;
    Ivec(time,1)   = inow;
    Cvec(time,1)   = cnow;
    Nvec(time,1)   = nnow;
    Wvec(time,1)   = wnow;
    Kvec(time+1,1) = kp;
    
    disp(sprintf('  time = %4d: pnow = %3.4f',time,pnow));
    disp([ynow inow cnow nnow wnow know znow]);
       
end

sd = 100*std(hpfilter(log([Yvec Ivec Cvec Nvec Wvec]),100));
Rvec = GAMY/BETA*Cvec(2:end)./Cvec(1:end-1);
sd(6) = 100*std(hpfilter(Rvec,100));
disp('    Standard deviations');
disp('    Y         I         C         N         W         R');
disp(sd);

save PlannerSim Kvec Kpvec Cvec izvec nz simT;