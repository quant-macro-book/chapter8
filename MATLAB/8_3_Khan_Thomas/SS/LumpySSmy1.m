% function [wnew v kw e0 xi theta kdist] = LumpySSmy1(w,znow,knotsk,nk,kbounds,Phi)
% LumpySSmy.m
% solve for new steady state w and distribution given old w
% January 5 2010

% convergences
global critv;

global GAMY BETA DELTA THETA NU ETA B kSS;

disp(sprintf('  INNER LOOP at w = %6.5f',w));

t = cputime;

% Part I. The Initial Value Function

n0 = zeros(nk,1);
v0 = zeros(nk,1);
for ik = 1:nk
            
    know = knotsk(ik);
    yterm = znow*know^THETA;
    n = (NU*yterm/w)^(1/(1-NU));
    n0(ik,1) = n;
    y = yterm*n^NU;
    v0temp = y - w*n + (1-DELTA)*know;
    v0(ik,1) = v0temp; % normalized by p0

end

v = v0;
kw = kSS;

% Part II. Iteration on Contraction Mapping

vnew = zeros(nk,1);

diff = 1e+4;
iter = 0;
s1 = 0;

% rk = nk-2;
% T = spbas(rk,knotsk);
% invT = inv(T);

c = Phi\v;

while(diff>critv)
    
    % stage 2a. calcurate the threshold adjustment costs for each k
%     cf = spfit(invT,v,rk,knotsk);
    cold = c;

    % solve for kp using the golden section search
    % note that the value of adjustment is independent of current k
    % target k
    if s1 == 0
        kwnew = golden('vfuncpoly',knotsk(1),knotsk(end),c,knotsk,kbounds);
    end
    
    e0 = -vfuncpoly(kwnew,c,knotsk,kbounds);
 
    for ik = 1:nk
        
        know = knotsk(ik,1);
        p1 = polybas((1-DELTA)*know/GAMY,nk,kbounds(1),kbounds(end));
        v1 = p1*c;
%         v1 = speva(cf,(1-DELTA)*know/GAMY,rk,knotsk);
        e1 = -(1-DELTA)*know + BETA*v1;
        xitemp = (e0-e1)/w;
        xi(ik,1) = min(B,max(0,xitemp));
        alpha(ik,1) = xi(ik,1)/B;
        % solve for vnew
        vnew(ik,1) = v0(ik,1) ...
        - w*xi(ik,1)^2/2/B ...
        + alpha(ik,1)*e0 + (1-alpha(ik,1))*e1;

        vjac(ik,:) = alpha(ik,1)*BETA*polybas(kwnew,nk,kbounds(1),kbounds(end)) ...
            + (1-alpha(ik,1))*BETA*polybas((1-DELTA)*know/GAMY,nk,kbounds(1),kbounds(end));   
    
    end
    
%     c = Phi\vnew;
%     vjac = BETA*polybas(gnew,nk,knotsk(1),knotsk(end));
    c = c - (Phi-vjac)\(Phi*c-vnew);
    
    diffkw = abs(kwnew-kw);  
    diffv = max(abs(vnew-v));
    diff = diffv;
    iter = iter+1;
%     s = sprintf( '  iteration  %4d:  ||Tkp-kp|| = %6.8f  ||Tv-v|| = %6.8f', ...
%         iter,diffkw,diffv);
%     disp(s);
    if (diffkw<1e-4 && s1==0);
        s1 = 1;
    elseif (s1>0 && s1<20);
        s1 = s1+1;
    elseif s1>=20
        s1 = 0;
    end
    kw = kwnew;
    v = vnew;
    
end

disp(sprintf('  Elapsed time = %6.8f',cputime-t));

% stage 3. calcurate the distribution of capital
% v and kw implies ss distribution

% pick up J, to be rewritten in a better way?
kidx = min(find(xi<B))-1; % upward hazard
kidx = min(nk,max(1,kidx));
kl = knotsk(kidx); % all firms adjust below this level of k
kdist(1) = kw;

for j=1:100
    kdist(j+1,1) = (1-DELTA)*kdist(j);
    if kdist(j+1,1) <= kl; break; end;
end

% adjustment probability alpha(jj) jj-1 to jj, jj = 1,...,j
xidist = interp1(knotsk,xi,kdist);
alpha1 = xidist./B; % G(xi), uniform distribution

% calcurate theta
a = 1;
A = 1;
for jj=1:j
    a = (1-alpha1(jj))*a;
    A = A + a;
end

theta = zeros(j+1,1);
theta(1) = 1/A;
for jj=1:j
    theta(jj+1) = (1-alpha1(jj))*theta(jj);
end

% work through distribution
for ik = 1:length(kdist)

    know = kdist(ik);
    % solve for n
    yterm = znow*know^THETA;
    nnow = (NU*yterm/w)^(1/(1-NU));

    inow = alpha1(ik)*(GAMY*kwnew - (1-DELTA)*know);
    ynow = znow*know^THETA*nnow^NU;

    % distribution
    idist(ik,1) = inow;
    ydist(ik,1) = ynow;
    ndist(ik,1) = nnow + xidist(ik)^2/2/B;

end

Kagg = theta'*kdist;
Iagg = theta'*idist;
Yagg = theta'*ydist;
Nagg = theta'*ndist;
Cagg = Yagg-Iagg;
pnew = 1/Cagg;

% disp('  ');
% disp('  aggregate quantities');
% disp([Yagg Iagg Cagg Nagg Kagg]);

% new w
pnew = 1/Cagg;
wnew = ETA*Cagg;