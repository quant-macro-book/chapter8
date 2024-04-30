function [pf_c pf_n pf_sav] = HH_opt_EGM(beta,p,grid_a,grid_z,prob_z,wToday,RToday,tauToday,dToday,pf_c_init,pf_n_init,pf_sav_init)

    %
    pf_c = pf_c_init;
    pf_n = pf_n_init;
    pf_sav = pf_sav_init;

    % WHILE LOOP: iterate until the policy functions converge.
    for EGM_iter = 1:p.MAX_EGM_ITER

        pf_c_int = zeros(p.N_A+1,p.N_Z);
        endog_grid_a = zeros(p.N_A+1,p.N_Z);

        % Assuming that the borrowing constraint is not binding, solve for the policy functions.
        for i_z = 1:p.N_Z % Fix z

            zToday      = grid_z(i_z);
            cond_prob_z = prob_z(i_z,:);

            for i_a = 1:p.N_A

                % EGM: construct a grid for a', rather than a grid for a.
                %     solve for the value of a that would have led to the choice a'.
                aTomorrow = grid_a(i_a);

                % Compute the expected marginal utility of consumption using a', and derive consumption today from the Euler equation.
                % Note that here we skip a root-finding algorithm, thereby making the algorithm much more efficient than the usual policy function iteration!
                EMUc = 0.0;
                for j_z = 1:p.N_Z
                    EMUc = EMUc + beta*RToday*(cond_prob_z(j_z)*pf_c(i_a,j_z)^(-p.GAMMA));
                end

                cToday = EMUc^(-1.0/p.GAMMA);

                % Derive labor supply today from the labor supply equation.
                nToday = (wToday*zToday*cToday^(-p.GAMMA))^(1.0/p.PSI);

                % Derive today's asset from the budget constraint. This current asset level is called the endogenous gridpoints.
                aToday = cToday + aTomorrow/RToday - wToday*zToday*nToday + tauToday*p.TAU_BAR(i_z) - dToday;

                % Store the values
                pf_c_int(i_a,i_z)     = cToday;
                endog_grid_a(i_a,i_z) = aToday;

            end % a

        end % z

        % The following is used for evaluating the new policy functions at the exogenous gridpoints using the pchip interpolation if the exogenous gridpoints are beyond the maximum endogenous gridpoints.
        % Here, pf_c_int(N_A+1,:) is linearly extrapolated.
        for i_z = 1:p.N_Z
            endog_grid_a(p.N_A+1,i_z) = 1e+5; % Note(20190825): Nakata-san set this as 1e8, rather than 1d5.
            pf_c_int(p.N_A+1,i_z) = pf_c_int(p.N_A,i_z) + ...
                ((pf_c_int(p.N_A,i_z) - pf_c_int(p.N_A-1,i_z)) / (endog_grid_a(p.N_A,i_z) - endog_grid_a(p.N_A-1,i_z))) ...
                * (endog_grid_a(p.N_A+1,i_z)-endog_grid_a(p.N_A,i_z));
            % The coefficient of the second term is just a slope.
        end

        % Evaluate the new policy functions at the exogenous gridpoints (or the original gridpoints).
        for i_z = 1:p.N_Z % Fix z

            zToday      = grid_z(i_z);
            cond_prob_z = prob_z(i_z,:);

            for i_a = 1:p.N_A

                aToday = grid_a(i_a);

                % In this case, endog_grid_a(1,:) is the value of bond holdings that induces the borrowing constraint to bind next period.
                % This is because the far left gridpoint in this program is set to the borrowing limit.
                if (aToday > endog_grid_a(1,i_z))

                    % The borrowing constraint does not bind.
                    % shape-preserving spline!!!
                    cToday = pchip(endog_grid_a(:,i_z),pf_c_int(:,i_z),aToday);

                    nToday = (wToday*zToday*cToday^(-p.GAMMA))^(1.0/p.PSI);

                    aTomorrow = (aToday + wToday*zToday*nToday - tauToday*p.TAU_BAR(i_z) + dToday - cToday)*RToday;

                else

                    % The borrowing costraint binds. Use the subroutine 'EGMConstrained' to compute cToday and nToday when the borrowing constraint is binding.
                    aTomorrow = p.BL;

                    % Call 'EGMConstrained' to obtain the values of cToday and nToday when the borrowing constraint is binding.
                    [cToday nToday] = EGMconstrained(p,aToday,zToday,RToday,wToday,tauToday,dToday,i_z);

                end

                % Obtain the new policy functions for consumption, labor supply, and savings as well as the associated value function.
                % Prepation for the value function
    %             Ev_int = 0d0
    %             do j_z = 1, p%N_Z
    %                 call interp1_pchip(vTomorrow, p%N_A, grid_a, aTomorrow, values(:,j_z))
    %                 Ev_int = Ev_int + cond_prob_z(j_z)*vTomorrow
    %             end do
    %             vToday = cToday**(1d0 - p%GAMMA)/(1d0 - p%GAMMA) - nToday**(1d0 + p%PSI)/(1d0 + p%PSI) + beta*Ev_int
    %             disp([i_a i_z grid_a(i_a) cToday nToday aTomorrow]);
    %             pause

                pf_c_new(i_a,i_z)   = cToday;
                pf_n_new(i_a,i_z)   = nToday;
                pf_sav_new(i_a,i_z) = aTomorrow;
                %values_new(i_a,i_z) = vToday

            end % a

        end % z

        % Evaluate convergence.
        EGM_err_c   = max(max(abs(pf_c_new - pf_c)));
        EGM_err_n   = max(max(abs(pf_n_new - pf_n)));
        EGM_err_sav = max(max(abs(pf_sav_new - pf_sav)));
        EGM_err     = max([EGM_err_c, EGM_err_n, EGM_err_sav]);

        if (mod(EGM_iter,50)==0); disp([EGM_iter EGM_err]); end;
    % 
    %     ! write(*,*) "------------------------------------------------------"
    %     ! write(*,*) "AT ITERATION   = ", EGM_iter
    %     ! write(*,*) "MAX DIFFERENCE = ", EGM_err
    %     ! write(*,*) "------------------------------------------------------"
    % 
        if ( EGM_err < p.EGM_ERR_TOL ); break; end;
    % 
    %     ! Update the policy functions for consumption, labor supply, and savings as well as the associated value function
        pf_c   = pf_c_new;
        pf_n   = pf_n_new;
        pf_sav = pf_sav_new;
    %     values = values_new
    % 
    end % End of WHILE LOOP over the policy functions

end


function [cToday nToday] = EGMconstrained(p,aToday,zToday,RToday,wToday,tauToday,dToday,i_z)

    % Initial guess of the value of labor supply and the associated consumption value.
    nToday = 0.6;

    cToday = aToday + wToday*zToday*nToday - tauToday*p.TAU_BAR(i_z) + dToday - p.BL/RToday;

    labor_eq_diff = (cToday^(-p.GAMMA))*wToday*zToday - nToday^p.PSI; % labor_eq_diff denotes the difference of the labor supply equation.

    % Initialization
    EGM_const_err = 100.0;

    % WHILE LOOP: iterate until we find a pair of labor supply and consumption that satisfies the labor supply equation.
    %             Here we use the Newton-Raphson method.
    for EGM_const_iter = 1:p.MAX_EGM_CONST_ITER

        labor_eq_adj  = -p.GAMMA*(cToday^(-p.GAMMA-1.0))*(wToday*zToday)^2.0 - p.PSI*nToday^(p.PSI-1.0);

        nToday        = nToday - labor_eq_diff/labor_eq_adj;

        cToday        = aToday + zToday*wToday*nToday - tauToday*p.TAU_BAR(i_z) + dToday - p.BL/RToday;

        labor_eq_diff = (cToday^(-p.GAMMA))*wToday*zToday - nToday^p.PSI;

        EGM_const_err = abs(labor_eq_diff);

    %     ! write(*,*) "------------------------------------------------------"
    %     ! write(*,*) "AT ITERATION   = ", EGM_const_iter
    %     ! write(*,*) "DIFFERENCE = ", EGM_const_err
    %     ! write(*,*) "------------------------------------------------------"
    %     disp([EGM_const_iter EGM_const_err labor_eq_adj nToday cToday]);
        if (EGM_const_err < p.EGM_CONST_ERR_TOL); break; end;

    end

end

% if (EGM_const_err >= p%EGM_CONST_ERR_TOL) then
%     write(*,*) "EGMConstrained did not converge"
% end if