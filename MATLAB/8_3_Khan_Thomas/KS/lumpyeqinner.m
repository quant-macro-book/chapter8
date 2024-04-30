function v = lumpyeqinner(BetaK,Betap,knotsk,knotsm,Z,Pi)

global critin;
global GAMY BETA DELTA THETA NU ETA B;

nk = length(knotsk);
nm = length(knotsm);
nz = length(Z);
rk = nk-2;
rm = nm-2;

disp('  INNER LOOP');
t = cputime;

% Part I. The Initial Value Function

v0 = zeros(nk,nm,nz);
nmat  = zeros(nk,nm,nz);
mpmat = zeros(nm,nz);
pmat  = zeros(nm,nz);
wmat  = zeros(nm,nz);

for im = 1:nm
    for iz = 1:nz
        
        mnow = knotsm(im);
        znow = Z(iz);
        
        X = [1 log(mnow)];
        mp = exp(BetaK(iz,:)*X');
        p0  = exp(Betap(iz,:)*X');
        
        mpmat(im,iz) = mp;
        pmat(im,iz) = p0;
        w0 = ETA/p0;
        wmat(im,iz) = w0;
        
        for ik=1:nk
            
            know = knotsk(ik);
            yterm = znow*know^THETA;
            n = (NU*yterm/w0)^(1/(1-NU));
            nmat(ik,im,iz) = n;
            y = yterm*n^NU;
            v0temp = y - w0*n + (1-DELTA)*know;
            v0(ik,im,iz) = v0temp*p0;

        end
        
    end
end

v = v0;
kp = zeros(nm,nz);

% Part II. Iteration on Contraction Mapping

vnew   = zeros(nk,nm,nz);
kpnew  = zeros(nm,nz);
e0     = zeros(nm,nz);
e1     = zeros(nk,nm,nz);
xi     = zeros(nk,nm,nz);
% critin = 1e-4;
diff   = 1e+4;
iter   = 0;
% f      = figure;
s1 = 0;

while(diff>critin)

    for iz = 1:nz

%         znow = Z(iz);

        vcond = zeros(nk,nm);
        for jz = 1:nz
            vcond = vcond + Pi(iz,jz)*v(:,:,jz); % E[V(k,K,z')|z]
        end
      
        cmat = spfit2(vcond,rk,rm,knotsk,knotsm);

        % target k
        if s1==0

            for im = 1:nm

                % current prices and next period's aggregate capital
                mp = mpmat(im,iz);
                p = pmat(im,iz);
                w  = wmat(im,iz);

                kpnew(im,iz) = golden('vfuncsp2',knotsk(1),knotsk(end),mp,p,cmat,knotsk,knotsm);
                
            end % loop for m'

        end

        % solve for xi(k,K,z)
        for im = 1:nm
            
            % current prices and next period's aggregate capital
            mp = mpmat(im,iz);
            p = pmat(im,iz);
            w  = wmat(im,iz);
            e0(im,iz) = -vfuncsp2(kpnew(im,iz),mp,p,cmat,knotsk,knotsm);
            
            for ik = 1:nk
                
                know = knotsk(ik);
                v1 = speva2(cmat,(1-DELTA)/GAMY*know,mp,rk,rm,knotsk,knotsm);
                e1(ik,im,iz) = -p*(1-DELTA)*know + BETA*v1;
                xitemp = (e0(im,iz)-e1(ik,im,iz))/p/w;
                xi(ik,im,iz) = min(B,max(0,xitemp));
                % solve for vnew
                vnew(ik,im,iz) = v0(ik,im,iz) ...
                - p*w*xi(ik,im,iz)^2/2/B ...
                + xi(ik,im,iz)/B*e0(im,iz) + (1-xi(ik,im,iz)/B)*e1(ik,im,iz);
            
            end % loop for k
            
        end % loop for m

    end % loop for z
    
    diffkp = max(max(abs(kpnew-kp)));
    diffv = max(max(max(abs(vnew-v))));
    diff = diffv;
    iter = iter + 1;
    s = sprintf( '  iteration  %4d:  ||Tkp-kp|| = %6.8f  ||Tv-v|| = %6.8f', ...
        iter,diffkp,diffv);
    if (mod(iter,100)==0); disp(s); end;
    if (diffkp<1e-4 && s1==0);
        s1 = 1;
    elseif (s1>0 && s1<20);
        s1 = s1+1;
    elseif s1>=20
        s1 = 0;
    end
    kp = kpnew;
    v = vnew;
    
end

disp(sprintf('  Elapsed time = %6.8f',cputime-t));
% save lumpyeqinner1mydat v kp;