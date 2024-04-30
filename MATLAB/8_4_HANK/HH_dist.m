function [dist meanC meanN meanA] = HH_dist(beta,p,grid_a,grid_a_NS,grid_z,prob_z,wToday,RToday,tauToday,dToday,pf_c,pf_n,pf_sav);

AA = sparse(p.N_A_NS*p.N_Z,p.N_A_NS*p.N_Z);

% Non-stochastic simulations for period t.
for i_z = 1:p.N_Z % Fix i_z

    zToday = grid_z(i_z);

    for i_a_NS = 1:p.N_A_NS % Fix i_a_NS

        aToday = grid_a_NS(i_a_NS);

        % Compute savings of households whose state variables are a_{i_a_NS} and z_{i_z}, using the policy function for savings.
        aTomorrow = pchip(grid_a,pf_sav(:,i_z),aToday);

        indexToday = p.N_A_NS*(i_z-1) + i_a_NS;

        % Find the index j such that grid_a_NS(j) < aTomorrow <= grid_a_NS(j+1).
        for k_a_NS = 1:p.N_A_NS-1

            if (aTomorrow <= p.MIN_A_NS) % Make sure that MIN_A_NS = grid_a_NS(1) (= BL)
                j_a_NS = 0;
                break;
            elseif (grid_a_NS(k_a_NS) < aTomorrow && aTomorrow <= grid_a_NS(k_a_NS+1))
                j_a_NS = k_a_NS;
                break;
            elseif (aTomorrow > p.MAX_A_NS)
                j_a_NS = p.N_A_NS;
                break;
            end

        end

        % Redistribute the current mass, Hist(i_a_NS,i_z), to the points (grid_a_NS(j), z(1)), (grid_a_NS(j), z(2)), (grid_a_NS(j), z(3)),
        % (grid_a_NS(j+1), z(1)), (grid_a_NS(j+1), z(2)), and (grid_a_NS(j+1), z(3)) according to the weights, weight_NS*prob(i_z, 1),  weight_NS*prob(i_z, 2),  weight_NS*prob(i_z, 3),
        % (1d0 - weight_NS)*prob(i_z, 1),  (1d0 - weight_NS)*prob(i_z, 2), and (1d0 - weight_NS)*prob(i_z, 3).
        % Note that 'Hist_up' denotes the end-of-period t distribution, while 'Hist' denotes the beginning-of-period t distribution.
        if ( j_a_NS == 0 )

            for j_z = 1:p.N_Z %

                indexTomorrow = p.N_A_NS*(j_z-1) + 1;
                AA(indexToday,indexTomorrow) = prob_z(i_z,j_z);

            end

        elseif ( j_a_NS == p.N_A_NS )

            for j_z = 1:p.N_Z %

                indexTomorrow = p.N_A_NS*(j_z-1) + p.N_A_NS;
                AA(indexToday,indexTomorrow) = prob_z(i_z,j_z);

            end

        else

            weight_NS = 1.0 - ((aTomorrow - grid_a_NS(j_a_NS))/(grid_a_NS(j_a_NS+1) - grid_a_NS(j_a_NS)));

            for j_z = 1:p.N_Z

                indexTomorrow = p.N_A_NS*(j_z-1) + j_a_NS;
                AA(indexToday,indexTomorrow) = prob_z(i_z,j_z)*weight_NS;
                AA(indexToday,indexTomorrow+1) = prob_z(i_z,j_z)*(1.0-weight_NS);

            end

        end

    end % a

end % z

opts.tol = p.NS_ERR_TOL;
[v d] = eigs(AA',1,'LM',opts);
dist = v/sum(v);
dist = max(0,reshape(dist,[p.N_A_NS p.N_Z]));

meanC = 0.0;
meanN = 0.0;
meanA = 0.0;
%	meanV = 0.0;

for i_z = 1:p.N_Z

    zToday = grid_z(i_z);

    for i_a_NS = 1:p.N_A_NS

        aToday = grid_a_NS(i_a_NS);

        cToday = pchip(grid_a,pf_c(:,i_z),aToday);
        nToday = pchip(grid_a,pf_n(:,i_z),aToday);
        aTomorrow = pchip(grid_a,pf_sav(:,i_z),aToday);
        if (aTomorrow < p.BL)
            aTomorrow = p.BL;
        end

%         if (i_a_NS == 1)
%             disp([cToday nToday aTomorrow]); 
%             pause;
%         end

        meanC = meanC + dist(i_a_NS,i_z)*cToday;
        meanN = meanN + dist(i_a_NS,i_z)*zToday*nToday; % Note: meanN is calculated as the mean value of the 'effective' labor supply, labor supply multiplied by idiosyncratic labor productivity.
        meanA = meanA + dist(i_a_NS,i_z)*aTomorrow; % Note: The original code uses aToday instead of aTomorrow, which was wrong.
%			meanV = meanV + (dist(i_a_NS,i_z)/sum(dist))*vToday

    end

end