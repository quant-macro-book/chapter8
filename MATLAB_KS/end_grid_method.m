 function policy = end_grid_method(icept, slope)
% Function end_grid_method
%  [policy] = end_grid_method( )
%
% Purpose:
%  Compute a policy function of Krusell & Smith model
%  by the Endogenous Gridpoint Method.
%
% Record of revisions:
%    Date     Programmer  Description of change
% ==========  ==========  =====================
% 09/11/2018  T. Yamada   Original code

global Params

%% ****** DEFINE VARIABLES ******
con0 = zeros(Params.na+1, Params.ne, Params.nk, Params.nz);
con1 = zeros(Params.na+1, Params.ne, Params.nk, Params.nz);
coh0 = zeros(Params.na+1, Params.ne, Params.nk, Params.nz);
coh1 = zeros(Params.na+1, Params.ne, Params.nk, Params.nz);
asset = zeros(Params.na, Params.ne, Params.nk, Params.nz);
%================================

%% INITIAL METRIC OF ERROR

it = 1;
maxit = 1000;
metric = 1.0;
toler = 1.0e-004;

%% INITIAL GUESS (HAND-TO-MOUTH)

con0(1, :, :, :) = 0.0;
coh0(1, :, :, :) = 0.0;

for z = 1:Params.nz
    for k = 1:Params.nk
        [rent, wage] = factor_price(Params.kgrid(k), Params.L_agg(z), Params.tfp(z), Params.alpha, Params.delta);
        for e = 1:Params.ne
            coh0(2:Params.na+1, e, k, z) = wage*Params.endow(e) + (1.0+rent)*Params.aprime(:);
            con0(2:Params.na+1, e, k, z) = coh0(2:Params.na+1, e, k, z);
        end
    end
end

% given current capital and approximate aggregation, compute the next period's capital

agg_cap_next = zeros(Params.nk, Params.nz);

for z = 1:Params.nz
    for k = 1:Params.nk
        agg_cap_next(k, z) = approx_agg(Params.kgrid(k), z, icept, slope);
    end
end

%% CONSUMPTION FUNCTION OVER CASH ON HAND

disp('--- now computing policy function ---');

while (it<=maxit && metric>toler)
    
%     disp(it)
%     disp(metric)

    con1(1, :, :, :) = 0.0;
    coh1(1, :, :, :) = 0.0;

    for z = 1:Params.nz
        for k = 1:Params.nk
            for e = 1:Params.ne
                for i = 1:Params.na
                    rhs = RHS_Euler(Params.aprime(i), e, agg_cap_next(k, z), z, con0, coh0);
                    con1(i+1, e, k, z) = rhs^(-1/Params.gamma);
                    coh1(i+1, e, k, z) = con1(i+1, e, k, z) + Params.aprime(i);
                end
            end
        end
    end

    % check convergence
    metric = max(abs((con1(2:Params.na+1, :, :, :)-con0(2:Params.na+1, :, :, :))./con0(2:Params.na+1, :, :, :)), [], 'all');

    % update policy funciton
    con0 = con1;
    coh0 = coh1;
    it   = it+1;

end

% for debug use only
%save tmp_egm.mat;
%load('tmp_egm.mat', 'coh0')

%% RETRIEVE CURRENT FINANCIAL ASSET

for z = 1:Params.nz
    for k = 1:Params.nk
        [rent, wage] = factor_price(Params.kgrid(k), Params.L_agg(z), Params.tfp(z), Params.alpha, Params.delta);
        for e = 1:Params.ne
            asset(:, e, k, z) = (coh0(2:Params.na+1, e, k, z) - wage.*Params.endow(e))./(1.0+rent);
        end
    end
end

policy = zeros(Params.nd, Params.ne, Params.nk, Params.nz);

for z = 1:Params.nz
    for k = 1:Params.nk
        for e = 1:Params.ne
            for i = 1:Params.nd
                % out of range exception
                if asset(1, e, k, z) > Params.grid(i)
                    policy(i, e, k, z) = Params.amin;
                elseif Params.grid(i) > asset(Params.na, e, k, z)
                    policy(i, e, k, z) = interp1(asset(:, e, k, z), Params.aprime, Params.grid(i), 'linear', 'extrap');
                else
                    policy(i, e, k, z) = interp1(asset(:, e, k, z), Params.aprime, Params.grid(i), 'linear', 'extrap');
                end
                if (policy(i, e, k, z) < Params.amin)
                    policy(i, e, k, z) = Params.amin;
                end
            end
        end
    end
end

return;
