% SPFIT Fit spline function on the function values given inversed basis matrix
% invT
% USAGE
%   [c] = spbas(invT,fval,r,knots);
% INPUTS
%   invT  :
%   fval  :
%   r     : the number of interior knot points
%   knots : the knot points
%
% OUTPUTS
%   c :  4x(r+1) matrix which store coefficient of cubic splines
%
% USES: none
%
% See also: SPBAS, SPEVA

% February 4 2010, Takeki Sunakawa
% takeki.sunakawa@gmail.com

function c = spfit(invT,fval,r,knots)

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