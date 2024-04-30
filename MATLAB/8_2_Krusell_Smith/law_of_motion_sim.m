function [k_path, z_path] = law_of_motion_sim(policy)
% Function law_of_motion
%  [k_path, z_path] = law_of_motion( )
%
% Purpose:
%  given policy function, compute wealth distribution.
%
% Record of revisions:
%    Date     Programmer  Description of change
% ==========  ==========  =====================
% 11/01/2019  T. Yamada   Original code

global Params

%% ***** DEFINE VARIABLES *****
emp0 = zeros(Params.numi, 1);
emp1 = zeros(Params.numi, 1);
cap0 = zeros(Params.numi, 1);
cap1 = zeros(Params.numi, 1);
u_path = zeros(Params.nums, 1);
%==============================

%% SET RANDOM SEED

rng(1024);

%% GENERATE A SEQUENCE OF AGGREGATE SHOCK

k_path = zeros(Params.nums, 1);
z_path = markov_sim(Params.tran_zz, Params.nz, Params.nums);

%% INITIAL DISTRIBUTION: 10% OF INDIVIDIALS START AS UNEMPLOYED

cap0(:) = 35.0;
emp0(1:round(0.9*Params.numi)) = 1; % employed
emp0(round(0.9*Params.numi)+1:Params.numi) = 2; % unemployed
k_path(1) = sum(cap0)/Params.numi;
u_path(1) = 0.1;

%% ITERATE DISTRIBUTION FORWARDLY

disp('--- now computing law of motion by simulation ---');

for it = 1:Params.nums-1

%     disp(it)
%     disp(k_path(it))
%     disp(u_path(it))

    % generate random variable for each period
    rndw = rand(1, Params.numi);

    for i = 1:Params.numi

        % worker         
        if emp0(i) == 1

            % savings in the next period
            kloc = locate(Params.kgrid, k_path(it));
            if kloc >= Params.nk
                kloc = Params.nk - 1;
            elseif kloc <= 1
                kloc = 1;
            end
            savings0 = interp1(Params.grid, policy(:, 1, kloc, z_path(it)), cap0(i), 'linear', 'extrap');
            savings1 = interp1(Params.grid, policy(:, 1, kloc+1, z_path(it)), cap0(i), 'linear', 'extrap');
            wght = (k_path(it) - Params.kgrid(kloc)) / (Params.kgrid(kloc+1) - Params.kgrid(kloc));
            if wght > 1.0
                wght = 1.0;
            elseif Params.kgrid(kloc+1) - Params.kgrid(kloc) < 0.0
                wght = 0.0;
            end
            cap1(i) = (1.0-wght)*savings0 + wght*savings1;

            % transition of employment status
            pr_next = Params.prob(1, z_path(it), 1, z_path(it+1)) / Params.tran_zz(z_path(it), z_path(it+1));
            for e = 1:Params.ne
                if rndw(i) <= pr_next
                    emp1(i) = e;
                    if emp1(i) ~= 0
                        break
                    end
                else
                    pr_next = pr_next + Params.prob(1, z_path(it), e+1, z_path(it+1)) / Params.tran_zz(z_path(it), z_path(it+1));
                end
            end

        end

        % unemployed
        if emp0(i) == 2

            % savings in the next period
            kloc = locate(Params.kgrid, k_path(it));
            if kloc >= Params.nk
                kloc = Params.nk - 1;
            elseif kloc <= 1
                kloc = 1;
            end
            savings0 = interp1(Params.grid, policy(:, 2, kloc, z_path(it)), cap0(i), 'linear', 'extrap');
            savings1 = interp1(Params.grid, policy(:, 2, kloc+1, z_path(it)), cap0(i), 'linear', 'extrap');
            wght = (k_path(it) - Params.kgrid(kloc)) / (Params.kgrid(kloc+1) - Params.kgrid(kloc));
            if wght > 1.0
                wght = 1.0;
            elseif Params.kgrid(kloc+1) - Params.kgrid(kloc) < 0.0
                wght = 0.0;
            end
            cap1(i) = (1.0-wght)*savings0 + wght*savings1;

            % transition of employment status
            pr_next = Params.prob(2, z_path(it), 1, z_path(it+1)) / Params.tran_zz(z_path(it), z_path(it+1));
            for e = 1:Params.ne
                if rndw(i) <= pr_next
                    emp1(i) = e;
                    if emp1(i) ~= 0
                        break
                    end
                else
                    pr_next = pr_next + Params.prob(2, z_path(it), e+1, z_path(it+1)) / Params.tran_zz(z_path(it), z_path(it+1));
                end
            end

        end

    end

    % aggregate capital and unemployment rate
    count = 0;
    for i = 1:Params.numi
        if emp1(i) == 2
            count = count + 1;
        end
    end
    u_path(it+1) = count/Params.numi;
    k_path(it+1) = sum(cap1)/Params.numi;

    % update
    emp0 = emp1;
    cap0 = cap1;
    emp1 = 0;
    cap1 = 0.0;

end

return;
