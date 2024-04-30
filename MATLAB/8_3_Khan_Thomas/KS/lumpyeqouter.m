function [Yvec Ivec Cvec Nvec Wvec Zvec Kvec Kpvec] = lumpyeqouter(v,BetaK,knotsk,knotsm,Z,Pi,izvec)

global GAMY BETA kSS;

nk = length(knotsk);
nm = length(knotsm);
nz = length(Z);
rk = nk-2;
rm = nm-2;

disp('  OUTER LOOP');
t = cputime;

simT = size(izvec,1);
Zvec  = zeros(simT,1);
Kvec  = zeros(simT,1);
Kpvec = zeros(simT,1);
Yvec  = zeros(simT,1);
Ivec  = zeros(simT,1);
Cvec  = zeros(simT,1);
Nvec  = zeros(simT,1);
Wvec  = zeros(simT,1);

Thetanow = 1;
Kvecnow = kSS;

cmat = zeros(16,rk+1,rm+1,nz);

for iz = 1:nz
    
    vcond = zeros(nk,nm);   
    for jz = 1:nz
        vcond = vcond + Pi(iz,jz)*v(:,:,jz); % E[V(k,K,z')|z]
    end
    cmat(:,:,:,iz) = spfit2(vcond,rk,rm,knotsk,knotsm);

end

    
for time = 1:simT;

    mnow = Thetanow'*Kvecnow;
    iz = izvec(time);

    X = [1 log(mnow)];
    mp = exp(BetaK(iz,:)*X');

    znow = Z(iz);
    cmatz = cmat(:,:,:,iz);

    klow = 0.5*mp;
    khigh = 1.5*mp;
    [ev1 edv1] = speva2(cmatz,klow,mp,rk,rm,knotsk,knotsm);
    phigh = BETA*edv1/GAMY;
    [ev2 edv2] = speva2(cmatz,khigh,mp,rk,rm,knotsk,knotsm);
    plow = BETA*edv2/GAMY;

    [pnew Thetanew Kvecnew Yagg Iagg Cagg Nagg wnew] = bisectp(plow,phigh,cmatz,knotsk,knotsm,mp,znow,Thetanow,Kvecnow);

    % update distribution
    Thetanow = Thetanew;
    Kvecnow  = Kvecnew;

    % record aggregate variables  
    Kagg = mnow;
    Zagg = znow;
    Kpagg = Thetanew'*Kvecnew;
%     Iagg = GAMY*Kpagg - (1-DELTA)*mnow;
%     Yagg = Cagg + Iagg;    

    Zvec(time)  = znow;
    Kvec(time)  = mnow;
    Kpvec(time) = Kpagg;
    Yvec(time)  = Yagg;
    Ivec(time)  = Iagg;
    Cvec(time)  = Cagg;
    Nvec(time)  = Nagg;
    Wvec(time)  = wnew;
    
    if (mod(time,100)==0); disp(sprintf('  time = %4d: pnow = %3.4f, pl = %3.4f',time,pnew,plow)); end;
    %disp([Yagg Iagg Cagg Nagg wnew Zagg Kagg]);
    
end

disp(sprintf('  Elapsed time = %6.8f',cputime-t));