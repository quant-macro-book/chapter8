% SPBAS Evaluate basis function of cubic spline T where T*s = F on r+2 knots
% USAGE
%   [T] = spbas(r,knots);
% INPUTS
%   r     : the name of the function to be maximized
%   knots : the left point of x
%
% OUTPUTS
%   T :  basis function T(k)*s = F where s is slope and F is ...
%
% USES: none
%
% See also:  SPFIT, SPEVA

% February 4 2010, Takeki Sunakawa
% takeki.sunakawa@gmail.com

function T = spbas(r,knots)

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