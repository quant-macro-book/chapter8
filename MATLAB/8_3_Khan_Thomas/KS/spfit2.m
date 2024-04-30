function cmat2 = spfit2(fmat,rx,ry,xknots,yknots)
% fit 2-dimensional spline function on the function values fmat
% Nov 18 2010 Takeki Sunakawa

Tx = spbas(rx,xknots);
invTx = inv(Tx);

dmat = zeros(4*(rx+1),ry+2);
for iy=1:ry+2

    fval = fmat(:,iy);
    % fit univariate x-spline
    d = spfit(invTx,fval,rx,xknots);
    dmat(:,iy) = reshape(d,4*(rx+1),1);

end

Ty = spbas(ry,yknots);
invTy = inv(Ty);

cmat = zeros(4,(ry+1),4*(rx+1));
for ix=1:4*(rx+1)

    fval = dmat(ix,:)';
    % fit univariate y-spline,
    c = spfit(invTy,fval,ry,yknots); % 4x(ry+1)
    cmat(:,:,ix) = c; %reshape(c,4*(rx+1),1);

end

% reshaping cmat, 16x(rx+1)x(ry+1)
cmat2 = zeros(16,rx+1,ry+1);
for iy=1:ry+1

    for ix=1:rx+1

        index = (4*(ix-1)+1:4*ix);
        tempcvec = reshape(cmat(:,iy,index),16,1);
        cmat2(:,ix,iy) = tempcvec; 
        % c_{1,tx,1,ty},c_{1,tx,2,ty},c_{1,tx,3,ty},c_{1,tx,4,ty},c_{2,tx,1
        % ,ty},c_{2,tx,2,ty},c_{2,tx,3,ty},c_{2,tx,4,ty},...

    end
    
end

function T = spbas(r,knots)
% evaluate basis function of cubic spline on r+2 knots
% Feb 4 2010 Takeki Sunakawa

T = zeros(r,r);

for i = 1:r
    
    k = i+1; % index on knots
    
    dt0 = knots(k)-knots(k-1);
    dt1 = knots(k+1)-knots(k);
    
    if i==1
        
        T(i,i) = 2*(dt0+dt1)-(dt1-dt0^2/dt1);
        T(i,i+1) = dt1+dt0^2/dt1;    

    elseif i==r
        
        T(i,i-1) = dt1+dt1^2/dt0;
        T(i,i) = 2*(dt0+dt1)-(dt0-dt1^2/dt0);
        
    else
        
        T(i,i-1) = dt0;
        T(i,i) = 2*(dt0+dt1);
        T(i,i+1) = dt1; 
    
    end
        
end

function c = spfit(invT,fval,r,knots)
% fit spline function on the function values given inversed basis matrix
% invT
% Feb 4 2010 Takeki Sunakawa

F = zeros(r,1);

for i = 1:r
    
    k = i+1; % index on knots
    
    f0 = fval(k-1);
    f1 = fval(k);
    f2 = fval(k+1);
    dt0 = knots(k)-knots(k-1);
    dt1 = knots(k+1)-knots(k);
    df0 = (f1-f0)/dt0;
    df1 = (f2-f1)/dt1;
    
    if i==1
        
        F(i) = 3*(dt1*df0+dt0*df1)-2*(dt1*df0-dt0^2/dt1*df1);
        
        a1 = 1-(dt0/dt1)^2;
        a2 = (dt0/dt1)^2;
        a0 = 2*(df0-(dt0/dt1)^2*df1);

    elseif i==r
        
        F(i) = 3*(dt1*df0+dt0*df1)-2*(dt0*df1-dt1^2/dt0*df0);

        b1 = 1-(dt1/dt0)^2;
        b2 = (dt1/dt0)^2;
        b0 = 2*(df1-(dt1/dt0)^2*df0);
        
    else
        
        F(i) = 3*(dt1*df0+dt0*df1);
    
    end
        
end

s = invT*F;

% not-a-knot
SVEC(1) = a0 - a1*s(1) + a2*s(2);
SVEC(2) = b0 - b1*s(r) + b2*s(r-1);
% natural spline
% SVEC = [0 0];
s = [SVEC(1);s;SVEC(2)];

c = zeros(4,r+1);

for i = 1:r+1
    
    k = i+1; % index on knots
    
    f0 = fval(k-1);
    f1 = fval(k);
    dt0 = knots(k)-knots(k-1);
    df0 = (f1-f0)/dt0;
    
    c(1,i) = f0;
    c(2,i) = s(i);
    c(3,i) = 3*df0/dt0 - 2*s(i)/dt0 - s(i+1)/dt0;
    c(4,i) = -2*df0/(dt0^2) + s(i)/(dt0^2) + s(i+1)/(dt0^2);
    
end