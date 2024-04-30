function [f0,dxf0,dyf0] = speva2(cmat,xnodes,ynodes,rx,ry,xknots,yknots)
% evaluate 2-dim spline function for input vector x,y and fitted coefficient c
% Nov 18 2010 Takeki Sunakawa

ns   = length(xnodes);
% ny   = length(ynodes);
f0   = zeros(ns,1);
% df0  = zeros(n,1);
% d2f0 = zeros(n,1);

for is = 1:ns

    x = xnodes(is);
    kx = 0;
    for jx = 1:rx+1
        if(xknots(jx)>x); break; end;
        kx = kx+1;
    end
%     k = sum(knots<=x);
    kx = max(kx,1);
    kx = min(kx,rx+1);
    tx = xknots(kx); % leftknot

    y = ynodes(is);
    ky = 0;
    for jy = 1:ry+1
        if(yknots(jy)>y); break; end;
        ky = ky+1;
    end
%     k = sum(knots<=x);
    ky = max(ky,1);
    ky = min(ky,ry+1);
    ty = yknots(ky); % leftknot

    c = cmat(:,kx,ky); % 16x1
%     kx
%     ky
%     c
%     pause;

    f0(is,1) = c(1) + c(2)*(y-ty) + c(3)*(y-ty)^2 + c(4)*(y-ty)^3 ...
        + (x-tx)*(c(5) + c(6)*(y-ty) + c(7)*(y-ty)^2 + c(8)*(y-ty).^3) ...
        + (x-tx)^2*(c(9) + c(10)*(y-ty) + c(11)*(y-ty)^2 + c(12)*(y-ty)^3) ...
        + (x-tx)^3*(c(13) + c(14)*(y-ty) + c(15)*(y-ty)^2 + c(16)*(y-ty)^3);

    dxf0(is,1) = c(5) + c(6)*(y-ty) + c(7)*(y-ty)^2 + c(8)*(y-ty).^3 ...
        + 2*(x-tx)*(c(9) + c(10)*(y-ty) + c(11)*(y-ty)^2 + c(12)*(y-ty)^3) ...
        + 3*(x-tx)^2*(c(13) + c(14)*(y-ty) + c(15)*(y-ty)^2 + c(16)*(y-ty)^3);

    dyf0(is,1) = c(2) + c(6)*(x-tx) + c(10)*(x-tx)^2 + c(14)*(x-tx).^3 ...
        + 2*(y-ty)*(c(3) + c(7)*(x-tx) + c(11)*(x-tx)^2 + c(15)*(x-tx)^3) ...
        + 3*(y-ty)^2*(c(4) + c(8)*(x-tx) + c(12)*(x-tx)^2 + c(16)*(x-tx)^3);
    
    end
    
end