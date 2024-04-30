% Planner.m
% solve the control model and generate simulated data
% January 5 2010

clear all;

% convergences
critv = 1e-4;
critn = 1e-5;
critg = 1e-5;
% parameters, taken from p.341 of Khan and Thomas (2003)
global GAMY BETA DELTA THETA NU ETA B;

GAMY  = 1.016;
BETA  = 0.954;
DELTA = 0.060;
THETA = 0.325;
NU    = 0.580;
ETA   = 3.614;
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

% stage 1. solve for the last-period value function as initial value.
v0 = zeros(nk,1);
n0 = zeros(nk,1);
nL = zeros(nk,1);
nH = zeros(nk,1);

for iz = 1:nz
    
    z0 = Z(iz);
    
    for ik = 1:nk
        
        k0 = knotsk(ik);
        % the Newton-Rhapson, to be replaced by my own code?
        % 20 percent of available time is spent in market work
        [n0temp rc] = csolve(@nfunc,0.2,[],critn,100,k0,0,z0);
        if rc ~= 0; disp(sprintf('something''s wrong at (k,z)=(%3.4f,%3.4f)',k0,z0)); end;
        % bounds on n
        kL = knotsk(1);
        [nLtemp rc] = csolve(@nfunc,0.2,[],critn,100,k0,kL,z0);
        if rc ~= 0; disp(sprintf('something''s wrong at (k,z)=(%3.4f,%3.4f)',know,znow)); end;
        kH = knotsk(end);
        [nHtemp rc] = csolve(@nfunc,0.2,[],critn,100,k0,kH,z0);
        if rc ~= 0; disp(sprintf('something''s wrong at (k,z)=(%3.4f,%3.4f)',know,znow)); end;

        yterm = z0*k0^THETA*n0temp^NU + (1-DELTA)*k0;
        
        c0 = yterm;
        util = log(c0) + ETA*(1-n0temp);

        v0(ik,iz) = util;
        n0(ik,iz) = n0temp;
        nL(ik,iz) = nLtemp;
        nH(ik,iz) = nHtemp;
    
    end
end

n = n0;
v = v0;
kp = zeros(nk,nz);

% stage 2. iteration to obtain v
nnew = zeros(nk,nz);
vnew = zeros(nk,nz);
kpnew = zeros(nk,nz);
diff = 1e+4;
iter = 0;
s1 = 0;

rk = nk-2;
T = spbas(rk,knotsk);
invT = inv(T);

while(diff>critv)
    
    for iz = 1:nz
        
        znow = Z(iz);
        
        vcond = v*Pi(iz,:)';
        cfz = spfit(invT,vcond,rk,knotsk);
        
        for ik = 1:nk
            
            know = knotsk(ik);
            
            if s1 == 0

                % solve for n
                nnow = golden('vfuncsp3',nL(ik,iz),nH(ik,iz),critg,cfz,knotsk,know,znow);
                nnew(ik,iz) = nnow;
                % solve for kp
                yterm = znow*know^THETA*nnow^NU + (1-DELTA)*know;
                % first-order condition of n
                cnow = (NU/ETA)*znow*know^THETA*nnow^(NU-1);
                kptemp = (yterm-cnow)/GAMY;
                kpnew(ik,iz) = kptemp;
                
            end

            vnew(ik,iz) = -vfuncsp3(nnew(ik,iz),cfz,knotsk,know,znow);
            
        end
        
    end
    
    diffv = max(max(abs(vnew-v)));
    diffkp = max(max(abs(kpnew-kp)));
    diff = diffv;
    iter = iter+1;
    s = sprintf( '  iteration  %4d:  ||Tkp-kp|| = %6.8f  ||Tv-v|| = %6.8f', ...
        iter,diffkp,diffv);
    disp(s);
    if (diffkp<1e-4 && s1==0);
        s1 = 1;
    elseif (s1>0 && s1<20);
        s1 = s1+1;
    elseif s1>=20
        s1 = 0;
    end    

    n = nnew;
    v = vnew;
    kp = kpnew;
    
end

save Planner v n kp;