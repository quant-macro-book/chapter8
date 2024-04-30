function [pnew Thetanew Kvecnew Yagg Iagg Cagg Nagg wnew] = pricemapmy(pvec,cmat,knotsk,knotsm,mp,znow,Thetanow,Kvecnow)

global GAMY BETA DELTA THETA NU ETA B;

rk = length(knotsk)-2;
rm = length(knotsm)-2;

np = length(pvec);
pnew = zeros(np,1);

for ip = 1:np

    p = pvec(ip);
    w = ETA/p;

    kpnew = golden('vfuncsp2',knotsk(1),knotsk(end),mp,p,cmat,knotsk,knotsm);
    e0 = -vfuncsp2(kpnew,mp,p,cmat,knotsk,knotsm);
    
    nk = length(Kvecnow);
    alpha = zeros(nk,1);
    Ivec = zeros(nk,1);
    Yvec = zeros(nk,1);
    Nvec = zeros(nk,1);

    for ik = 1:nk

        know = Kvecnow(ik);
        % solve for n
        yterm = znow*know^THETA;
        nnow = (NU*yterm/w)^(1/(1-NU));
        % solve for xi
        v1 = speva2(cmat,(1-DELTA)*know/GAMY,mp,rk,rm,knotsk,knotsm);
        e1 = -p*(1-DELTA)*know + BETA*v1;
        xitemp = (e0-e1)/p/w;
        xi = min(B,max(0,xitemp));       
        % adjusting probability
        alpha(ik) = xi/B;

        inow = alpha(ik)*(GAMY*kpnew - (1-DELTA)*know);
        ynow = znow*know^THETA*nnow^NU;
        
        % distribution
        Ivec(ik) = inow;
        Yvec(ik) = ynow;
        Nvec(ik) = nnow + xi^2/2/B;

    end
    
    % next theta and kvec
    J = max(find(alpha<1));
    
    Thetanew = zeros(J+1,1);
    Kvecnew = zeros(J+1,1);
    Thetanew(1,1) = alpha'*Thetanow;
    Thetanew(2:J+1,1) = (1-alpha(1:J)).*Thetanow(1:J);
    Kvecnew(1,1) = kpnew;
    Kvecnew(2:J+1,1) = (1-DELTA)/GAMY*Kvecnow(1:J);   
    
    Kpagg = Thetanew'*Kvecnew;
    Iagg = Thetanow'*Ivec;
    Yagg = Thetanow'*Yvec;   
    Cagg = Yagg-Iagg;
    Nagg = Thetanow'*Nvec;
    pnew(ip) = 1/Cagg;
    wnew = ETA*Cagg;
    
end