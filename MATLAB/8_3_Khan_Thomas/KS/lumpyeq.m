% Lumpyeq.m
% solve the lumpy investment model with Krusell-Smith algorithm
% January 5 2010

clear all;

tic;

% criteria for convergence
global critin critbp critg;
critout = 1e-4;  % for outer loop
critin  = 1e-5;  % for inner loop
critbp  = 1e-10; % for bisection to solve p = F(p)
critg   = 1e-10; % for golden section search to solve for kw

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
rk = nk-2;

mbounds = [0.8 1.2];
nm = 5;
knotsm = linspace(mbounds(1),mbounds(2),nm)';
rm = nm-2;

% % realization of shocks are the same among outer loops
% simT = 500;
% cumPi= cumsum(Pi')';
% izvec = zeros(simT,1);
% izvec(1) = 1;
% for tt = 1:simT-1
%     cumPivec = cumPi(izvec(tt),:);
%     izvec(tt+1) = sum(rand-cumPivec>=0);
%     izvec(tt+1) = min(izvec(tt+1)+1,nz);
% end

load ../Planner/PlannerSim Kvec Kpvec Cvec izvec nz simT;
[BetaK Betap] = CEpolicylin(Z(izvec),izvec,Kvec',Kpvec',1./Cvec',nz,simT);

diff = 1e+4;
iter = 0;

while(diff>critout)

    t = cputime;
    
    v = lumpyeqinner(BetaK,Betap,knotsk,knotsm,Z,Pi);

    [Yvec Ivec Cvec Nvec Wvec Zvec Kvec Kpvec] = lumpyeqouter(v,BetaK,knotsk,knotsm,Z,Pi,izvec);

    for iz = 1:nz

        x = Kvec(izvec==iz);
        y = Kpvec(izvec==iz);
        p = 1./Cvec(izvec==iz);

        X = [ones(length(x),1) log(x)];
        y = [log(y) log(p)];
        Beta = (X'*X)\(X'*y);
        BetaKnew(iz,:) = Beta(:,1)';
        Betapnew(iz,:) = Beta(:,2)';

    end

    diffmp = max(max(abs(BetaKnew-BetaK)));
    diffp = max(max(abs(Betapnew-Betap)));
    diff = max(diffmp,diffp);
    iter = iter + 1;
    s = sprintf( '  iteration  %4d:  ||Tmp-mp|| = %3.4f  ||Tp-p|| = %3.4f  Elapsed time = %3.4f'  , ...
        iter,diffmp,diffp,cputime-t);
    disp(s);
    
    BetaK = BetaKnew;
    Betap = Betapnew;

    for iz = 1:nz
    
        figure;
        subplot(3,2,iz);
        plot(log(Kvec(izvec==iz)),log(Kpvec(izvec==iz)),'b*');
        figure;
        subplot(3,2,iz);
        plot(log(Kvec(izvec==iz)),log(1./Cvec(izvec==iz)),'r*');
        
    end
    
    pause;    
    
end

CEpolicylin(Z(izvec)', izvec, Kvec', Kpvec', 1./Cvec', nz, simT);

toc;

save lumpyeqksresult v BetaK Betap Kvec Kpvec Cvec izvec;
save lumpyeqksresult2 Yvec Ivec Cvec Nvec Wvec Zvec Kvec Kpvec;