% SPEVA Evaluate spline function for input vector x and fitted coefficient c
% USAGE
%   [f0 df0 d2f0] = speva(c,xnodes,r,knots);
% INPUTS
%   c      : the coefficients of cubic splines
%   xnodes : evaluation points
%   r      : the number of interior knot points
%   knots  : the knot points
%
% OUTPUTS
%   f0   :  
%   df0  :  
%   d2f0 :  
%
% USES: none
%
% See also: SPBAS, SPFIT

% February 4 2010, Takeki Sunakawa
% takeki.sunakawa@gmail.com

function [f0 df0 d2f0] = speva(c,xnodes,r,knots)

n    = length(xnodes);
f0   = zeros(n,1);
df0  = zeros(n,1);
d2f0 = zeros(n,1);

for i = 1:n

    x = xnodes(i);
    k = 0;
    for j = 1:r+2
        if(knots(j)>x); break; end;
        k = k+1;
    end
    
    k = min(k,r+1);
    k = max(k,1);

% extrapolation: linear
%     if k<=1
%         x1 = knots(1);
%         x2 = knots(2);
%         f1 = c(1,1);
%         f2 = c(1,2);
%         s = (f2-f1)/(x2-x1);
%         f0(i) = s*(x-x1)+f1;
%         df0(i) = s;
%         d2f0(i) = 0;
%     elseif k>=r+2
%         x1 = knots(r+1); % leftknot
%         x2 = knots(r+2);
%         f1 = c(1,r+1);
%         f2 = c(1,r+1) + c(2,r+1).*(x2-x1) + c(3,r+1).*(x2-x1).^2 + c(4,r+1).*(x2-x1).^3;
%         s = (f2-f1)/(x2-x1);
%         f0(i) = s*(x-x1)+f1;
%         df0(i) = s;
%         d2f0(i) = 0;
%     else
        t = knots(k); % leftknot
        f0(i)   = c(1,k) + c(2,k).*(x-t) + c(3,k).*(x-t).^2 + c(4,k).*(x-t).^3;
        df0(i)  = c(2,k) + 2*c(3,k).*(x-t) + 3*c(4,k).*(x-t).^2;
        d2f0(i) = 2*c(3,k) + 6*c(4,k).*(x-t);
%     end
    
end