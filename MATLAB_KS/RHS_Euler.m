function RHS = RHS_Euler(asset, e_state, K_agg, z_state, consf, xgrid)
% Function RHS_Euler
%  [RHS] = RHS_Euler(conf, xgrid, asset, z, rent, wage)
%
% Purpose:
%  Compute the RHS of the Euler equation.
%  See Carroll (2006,EL).
%
% Record of revisions:
%    Date     Programmer  Description of change
% ==========  ==========  =====================
% 09/11/2019  T. Yamada   Original code

global Params

%% ***** LOCAL VARIABLES *****
cons = zeros(Params.ne, Params.nz);
mu = zeros(Params.ne, Params.nz);
%=============================

for z = 1:Params.nz
    for e = 1:Params.ne

        % facor prices in the next period
        [rent, wage] = factor_price(K_agg, Params.L_agg(z), Params.tfp(z), Params.alpha, Params.delta);
        coh = wage*Params.endow(e) + (1.0+rent)*asset;

        % use linear interpolation over aggregate capital grid
        kloc = locate(Params.kgrid, K_agg);
        if kloc >= Params.nk
            kloc = Params.nk - 1;
        elseif kloc < 1
            kloc = 1;
        end
        weight = (K_agg - Params.kgrid(kloc))/(Params.kgrid(kloc+1) - Params.kgrid(kloc));
        if weight > 1.0
            weight = 1.0;
        elseif weight < 0.0
            weight = 0.0;
        end

        % lower capital grid
        if coh < xgrid(1, e, kloc, z)
            cons0 = coh;
        elseif coh > xgrid(Params.na, e, kloc, z)
            cons0 = interp1(xgrid(2:Params.na+1, e, kloc, z), consf(2:Params.na+1, e, kloc, z), coh, 'linear', 'extrap');
        else
            cons0 = interp1(xgrid(2:Params.na+1, e, kloc, z), consf(2:Params.na+1, e, kloc, z), coh, 'linear', 'extrap');
        end

        % upper capital grid
        if coh < xgrid(1, e, kloc+1, z)
            cons1 = coh;
        elseif coh > xgrid(Params.na, e, kloc+1, z)
            cons1 = interp1(xgrid(2:Params.na+1, e, kloc+1, z), consf(2:Params.na+1, e, kloc+1, z), coh, 'linear', 'extrap');
        else
            cons1 = interp1(xgrid(2:Params.na+1, e, kloc+1, z), consf(2:Params.na+1, e, kloc+1, z), coh, 'linear', 'extrap');
        end

        cons(e, z) = (1.0-weight)*cons0 + weight*cons1;

        % marginal utility x gross interest rate
        mu(e, z)   = mu_CRRA(cons(e, z), Params.gamma)*(1.0+rent);

    end
end

% expected value

exp_value = 0.0;

for z = 1:Params.nz
    for e = 1:Params.ne
        exp_value = Params.prob(e_state, z_state, e, z)*mu(e, z) + exp_value;
    end
end

RHS = Params.beta*exp_value;

return;
