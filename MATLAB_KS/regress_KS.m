function [icept, slope, R2] = regress_KS(k_path, z_path)
% Function regress
%  [icept, slope, R2] = regress_KS( )
%
% Purpose:
%  Given k and z path, compute regression coefficients.
%
%  log(K(t+1)) = a(z) + b(z)*log(K(t))
%
% Input:
%  No input
%
% Output:
%  icept: intercepts of regression equation
%  slope: slopes of regression equation
%  R2:    R square of regression equation
%
%  Record of revisions:
%    Date      Programmer  Description of change
%  ==========  ==========  =====================
%  11/04/2019  T. Yamada   Original

global Params

% define variable
icept = zeros(Params.nz, 1);
slope = zeros(Params.nz, 1);
R2 = zeros(Params.nz, 1);

% count realized states(discard first 1000 periods)
nobs = zeros(Params.nz, 1);

for it = 1001:Params.nums-1
    for z = 1:Params.nz
        if z_path(it) == z
            nobs(z) = nobs(z) + 1;
        end
    end
end

% construct data
kk_path = zeros(nobs(1), 1);

% find regression parameters

counter = ones(Params.nz, 1);

for z = 1:Params.nz

    if z ~= 1
        kk_path = zeros(nobs(z), 2);
    end

    for it = 1001:Params.nums-1
        if z_path(it) == z
            kk_path(counter(z), 1) = log(k_path(it));
            kk_path(counter(z), 2) = log(k_path(it+1));
            counter(z)             = counter(z) + 1;
        end
    end
    
    % regression: use Matlab's library
    coef = fitlm(kk_path(:,1), kk_path(:,2));
    icept(z) = table2array(coef.Coefficients(1, 1));
    slope(z) = table2array(coef.Coefficients(2, 1));
    R2(z)    = coef.Rsquared.Ordinary(1);

end

return;
